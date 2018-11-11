//  BSD 3-Clause License
//
//  Copyright (c) 2018, Bailin Deng, Yue Peng
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
//  * Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "GeometrySolver.h"
#include "MeshTypes.h"
#include "FileParser.h"
#include <iostream>
#include "Parameters.h"

bool read_handle_files(const char *index_filename, const char *xy_filename,
                       std::vector<int> &handle_indices,
                       Matrix2X &coordinates) {
  if (!read_sequence_file(index_filename, handle_indices)) {
    return false;
  }

  if (handle_indices.empty()) {
    std::cerr << "Error: no indices read" << std::endl;
    return false;
  }

  std::vector<Scalar> x_coords, y_coords;

  if (read_sequence_file(xy_filename, x_coords, y_coords)) {
    if (x_coords.size() != handle_indices.size()
        || y_coords.size() != handle_indices.size()) {
      std::cerr << "Error: inconsistent number of handles and coordinates"
                << std::endl;
      return false;
    }

    coordinates.resize(2, x_coords.size());
    for (int i = 0; i < static_cast<int>(x_coords.size()); ++i) {
      coordinates.col(i) = Vector2(x_coords[i], y_coords[i]);
    }

    std::cout << "Successfully read handle files" << std::endl;
    return true;
  }

  std::cerr << "Error reading coordinate files" << std::endl;
  return false;
}

void ARAP2D_with_handles(const TriMesh &mesh,
                         const Matrix2X &problem_def_coords,
                         const Matrix2X &init_coords,
                         const std::vector<int> &handle_indices,
                         const Matrix2X &handle_coords, int iter,
                         Parameters::GlobalUpdateType global_update_type,
                         int Anderson_m, SPDSolverType spd_solver,
                         const char *output_filename = NULL) {
  GeometrySolver<2> solver;

  for (TriMesh::ConstFaceIter cf_it = mesh.faces_begin();
      cf_it != mesh.faces_end(); ++cf_it) {
    std::vector<int> vtx_idx;
    for (TriMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*cf_it);
        cfv_it.is_valid(); ++cfv_it) {
      vtx_idx.push_back(cfv_it->idx());
    }

    solver.add_constraint(
        new ARAP2DTriangleConstraint(vtx_idx, 1.0, problem_def_coords));
  }

  if (solver.initialize(mesh.n_vertices(), global_update_type, spd_solver,
                        handle_indices)) {
    Matrix2X input_coords = init_coords;

    for (int i = 0; i < static_cast<int>(handle_indices.size()); ++i) {
      input_coords.col(handle_indices[i]) = handle_coords.col(i);
    }

    solver.solve(input_coords, iter, Anderson_m);
    solver.output_iteration_history(AA_SOLVER);

    if (output_filename) {
      TriMesh new_mesh = mesh;
      Matrix3X opt_coords(3, mesh.n_vertices());
      opt_coords.block(0, 0, 2, opt_coords.cols()) = solver.get_points();
      opt_coords.row(2).setZero();
      set_vertex_points(new_mesh, opt_coords);
      if (!OpenMesh::IO::write_mesh(new_mesh, output_filename,
                                    OpenMesh::IO::Options::Default, 16)) {
        std::cerr << "Error: unable to save mesh to file " << output_filename
                  << std::endl;
      }
    }
  } else {
    std::cerr << "Error: unable to initialize solver" << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 7) {
    std::cerr
        << "Usage:	 <ARAP2D> <PROBLEM_DEFINITION_MESH> <INIT_MESH> <HANDLE_INDEX_FILE> <HANDLE_XY_COORDS_FILE> <OPTIONS_FILE> <OUTPUT_MESH>"
        << std::endl;
    return 1;
  }

  TriMesh problem_def_mesh;
  if (!OpenMesh::IO::read_mesh(problem_def_mesh, argv[1])) {
    std::cerr << "Error: unable to read problem definition mesh from the file "
              << argv[1] << std::endl;
    return 1;
  }

  TriMesh init_mesh;
  if (!OpenMesh::IO::read_mesh(init_mesh, argv[2])) {
    std::cerr << "Error: unable to read initial mesh from the file " << argv[1]
              << std::endl;
    return 1;
  }

  std::vector<int> handle_indices;
  Matrix2X handle_coords;
  if (!read_handle_files(argv[3], argv[4], handle_indices, handle_coords)) {
    return 1;
  }

  int n_vtx = problem_def_mesh.n_vertices();
  for (int i = 0; i < static_cast<int>(handle_indices.size()); ++i) {
    if (handle_indices[i] < 0 || handle_indices[i] >= n_vtx) {
      std::cerr << "Error: invalid handle index: " << handle_indices[i]
                << std::endl;
      return 1;
    }
  }

  Parameters param;
  if (!param.load(argv[5])) {
    std::cerr << "Error: unable to load option file " << argv[5] << std::endl;
    return 1;
  }
  if (!param.valid_parameters()) {
    std::cerr << "Invalid filter options. Aborting..." << std::endl;
    return 1;
  }
  param.output();

  Matrix3X coords_3d;
  get_vertex_points(init_mesh, coords_3d);
  Matrix2X init_coords_2d = coords_3d.block(0, 0, 2, coords_3d.cols());

  get_vertex_points(problem_def_mesh, coords_3d);
  Matrix2X problem_def_coords_2d = coords_3d.block(0, 0, 2, coords_3d.cols());

  if (!param.acceleration)
    param.anderson_m = 0;

  std::cout << "========== Accelerate solver begins... ============="
            << std::endl;
  ARAP2D_with_handles(problem_def_mesh, problem_def_coords_2d, init_coords_2d,
                      handle_indices, handle_coords, param.iter,
                      param.global_update_type, param.anderson_m, LDLT_SOLVER,
                      argv[6]);
  std::cout << std::endl;

  return 0;
}

