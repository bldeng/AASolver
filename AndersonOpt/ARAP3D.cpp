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

#include "MeshTypes.h"
#include "GeometrySolver.h"
#include "FileParser.h"
#include "TetMeshIO.h"
#include <iostream>
#include "Parameters.h"

bool read_handle_files(const char *index_filename, const char *xyz_filename,
                       std::vector<int> &handle_indices,
                       Matrix3X &coordinates) {
  if (!read_sequence_file(index_filename, handle_indices)) {
    return false;
  }

  if (handle_indices.empty()) {
    std::cerr << "Error: no indices read" << std::endl;
    return false;
  }

  std::vector<Scalar> x_coords, y_coords, z_coords;

  if (read_sequence_file(xyz_filename, x_coords, y_coords, z_coords)) {
    if (x_coords.size() != handle_indices.size()
        || y_coords.size() != handle_indices.size()
        || z_coords.size() != handle_indices.size()) {
      std::cerr << "Error: inconsistent number of handles and coordinates"
                << std::endl;
      return false;
    }

    coordinates.resize(3, x_coords.size());
    for (int i = 0; i < static_cast<int>(x_coords.size()); ++i) {
      coordinates.col(i) = Vector3(x_coords[i], y_coords[i], z_coords[i]);
    }

    std::cout << "Successfully read handle files" << std::endl;
    return true;
  }

  std::cerr << "Error reading coordinate files" << std::endl;
  return false;
}

bool ARAP_tetrahedrons_with_handles(
    const Matrix3X &problem_def_poisitions, const Matrix3X &init_positions,
    const Eigen::Matrix4Xi &tet_connectivity,
    const Eigen::Matrix3Xi &tet_boundary,
    const std::vector<int> &handle_indices, const Matrix3X &handle_coords,
    int iter, Parameters::GlobalUpdateType global_update_type, int Anderson_m,
    bool use_fast_svd, SPDSolverType spd_solver, Matrix3X &opt_positions,
    const std::string &output_filename) {
  GeometrySolver<3> solver;

  int n_tets = tet_connectivity.cols();
  for (int i = 0; i < n_tets; ++i) {
    std::vector<int> vtx_idx;
    for (int j = 0; j < 4; ++j) {
      vtx_idx.push_back(tet_connectivity(j, i));
    }

    solver.add_constraint(
        new ARAP3DTetConstraint(vtx_idx, 1.0, problem_def_poisitions,
                                use_fast_svd));
  }

  int n_vtx = init_positions.cols();
  if (solver.initialize(n_vtx, global_update_type, spd_solver,
                        handle_indices)) {
    Matrix3X input_coords = init_positions;
    for (int i = 0; i < static_cast<int>(handle_indices.size()); ++i) {
      input_coords.col(handle_indices[i]) = handle_coords.col(i);
    }

    solver.solve(input_coords, iter, Anderson_m);
    solver.output_iteration_history(AA_SOLVER);
    opt_positions = solver.get_points();

    TetMeshIO::save_boundary_mesh(output_filename, opt_positions, tet_boundary);

    return true;
  }

  std::cerr << "Error: unable to initialize solver" << std::endl;
  return false;
}

int main(int argc, char **argv) {
  if (argc != 7) {
    std::cerr
        << "Usage:	 <ARAP3D> <PROBLEM_DEFINITION_MODEL_NAME> <INIT_COORD_MODEL_NAME> <HANDLE_INDEX_FILE> <HANDLE_XYZ_COORDS_FILE> <OPTIONS_FILE> <OUTPUT_MODEL_NAME>"
        << std::endl;
    return 1;
  }

  Matrix3X problem_definitin_coords_3d, init_coords_3d;
  Eigen::Matrix4Xi tet_connectivity, tet_connectivity_2;
  Eigen::Matrix3Xi tet_boundary, tet_boundary_2;

  if (!TetMeshIO::load(argv[1], problem_definitin_coords_3d, tet_connectivity,
                       tet_boundary)) {
    std::cerr << "Error: unable to read problem definition mesh from the file "
              << argv[1] << std::endl;
    return 1;
  }

  if (!TetMeshIO::load(argv[2], init_coords_3d, tet_connectivity_2,
                       tet_boundary_2)) {
    std::cerr << "Error: unable to read initial mesh from the file " << argv[1]
              << std::endl;
    return 1;
  }

  std::vector<int> handle_indices;
  Matrix3X handle_coords, opt_positions;
  if (!read_handle_files(argv[3], argv[4], handle_indices, handle_coords)) {
    return 1;
  }

  int n_vtx = problem_definitin_coords_3d.cols();
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

  if (!param.acceleration)
    param.anderson_m = 0;

  bool use_fast_svd = true;  // Change this to false in order to use more accurate Jacobi SVD
  if (use_fast_svd) {
    std::cout << "Using fast SVD solver..." << std::endl;
  } else {
    std::cout << "Using Jacobi SVD solver..." << std::endl;
  }

  ARAP_tetrahedrons_with_handles(problem_definitin_coords_3d, init_coords_3d,
                                 tet_connectivity, tet_boundary, handle_indices,
                                 handle_coords, param.iter,
                                 param.global_update_type, param.anderson_m,
                                 use_fast_svd, LDLT_SOLVER, opt_positions,
                                 argv[6]);

  std::cout << std::endl;

  return 0;
}
