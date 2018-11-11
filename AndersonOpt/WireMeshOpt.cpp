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
#include "Constraint.h"
#include "MeshTypes.h"
#include <iostream>
#include <random>
#include "Parameters.h"

void check_wiremesh_error(const PolyMesh &mesh, double target_edge_length,
                          double min_angle_radian, double max_angle_radian) {
  VectorX edge_err(mesh.n_edges());
  VectorX angle_err(mesh.n_faces() * 4);

  for (PolyMesh::ConstEdgeIter ce_it = mesh.edges_begin();
      ce_it != mesh.edges_end(); ++ce_it) {
    edge_err(ce_it->idx()) = std::fabs(
        mesh.calc_edge_length(*ce_it) - target_edge_length);
  }
  edge_err /= target_edge_length;

  for (PolyMesh::ConstFaceIter cf_it = mesh.faces_begin();
      cf_it != mesh.faces_end(); ++cf_it) {
    std::vector<Vector3> vtx;
    for (PolyMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*cf_it);
        cfv_it.is_valid(); ++cfv_it) {
      vtx.push_back(to_eigen_vec3(mesh.point(*cfv_it)));
    }

    assert(static_cast<int>(vtx.size()) == 4);
    for (int i = 0; i < 4; ++i) {
      Vector3 e1 = (vtx[(i + 1) % 4] - vtx[i]).normalized(), e2 = (vtx[(i + 3)
          % 4] - vtx[i]).normalized();
      Scalar angle = std::acos(e1.dot(e2));

      if (angle < min_angle_radian) {
        angle_err[4 * cf_it->idx() + i] = min_angle_radian - angle;
      } else if (angle >= max_angle_radian) {
        angle_err[4 * cf_it->idx() + i] = angle - max_angle_radian;
      } else {
        angle_err[4 * cf_it->idx() + i] = 0;
      }
    }
  }

  angle_err *= (180.0 / M_PI);

  std::cout << "Edge length error (relative): max " << edge_err.maxCoeff()
      << ",  average " << edge_err.mean() << std::endl;
  std::cout << "Angle error: max " << angle_err.maxCoeff() << ",  average "
      << angle_err.mean() << std::endl;
}

void check_ref_surface_distance(const PolyMesh &mesh, const TriMesh &ref_mesh) {
  Matrix3X ref_mesh_points, mesh_pts;
  Eigen::Matrix3Xi ref_mesh_faces;
  get_vertex_points(ref_mesh, ref_mesh_points);
  get_face_vertex_index(ref_mesh, ref_mesh_faces);
  get_vertex_points(mesh, mesh_pts);

  MatrixX3 V_ref = ref_mesh_points.transpose(), P = mesh_pts.transpose();
  Eigen::MatrixX3i F_ref = ref_mesh_faces.transpose();

  igl::AABB<MatrixX3, 3> tree;
  tree.init(V_ref, F_ref);
  VectorX sqrD;
  Eigen::VectorXi I;
  MatrixX3 C;
  tree.squared_distance(V_ref, F_ref, P, sqrD, I, C);
  Scalar edge_length = average_edge_length(mesh);
  VectorX dist_err = sqrD.array().sqrt().matrix() / edge_length;
  std::cout << "Reference surface distance (normalized by edge length): ";
  std::cout << "Max " << dist_err.maxCoeff() << ", ";
  std::cout << "Average " << dist_err.mean() << std::endl;
}

// Anderson_M = 0 means no acceleration
void optimize_mesh(const PolyMesh &mesh, const TriMesh &ref_mesh, int iter,
                   Parameters::GlobalUpdateType global_update_type,
                   int Anderson_M, SPDSolverType spd_solver,
                   Scalar angle_weight, Scalar min_angle_radian,
                   Scalar max_angle_radian, Scalar edge_length_weight,
                   Scalar edge_length, Scalar closeness_weight,
                   const char *output_filename = NULL) {
  Matrix3X p;
  get_vertex_points(mesh, p);

  Matrix3X ref_mesh_points;
  Eigen::Matrix3Xi ref_mesh_faces;
  get_vertex_points(ref_mesh, ref_mesh_points);
  get_face_vertex_index(ref_mesh, ref_mesh_faces);

  GeometrySolver<3> solver;

  if (closeness_weight > 0) {
    solver.add_constraint(
        new ReferenceSurfceConstraint(p.cols(), closeness_weight,
                                      ref_mesh_points, ref_mesh_faces));
  }

  for (PolyMesh::ConstFaceIter f_it = mesh.faces_begin();
      f_it != mesh.faces_end(); ++f_it) {
    std::vector<int> id_vector;
    for (PolyMesh::ConstFaceVertexIter fv_it = mesh.cfv_iter(*f_it);
        fv_it.is_valid(); ++fv_it) {
      id_vector.push_back(fv_it->idx());
    }

    assert(static_cast<int>(id_vector.size()) == 4);
    for (int i = 0; i < 4; ++i) {
      solver.add_constraint(
          new AngleConstraint<3>(id_vector[i], id_vector[(i + 1) % 4],
                                 id_vector[(i + 3) % 4], angle_weight,
                                 min_angle_radian, max_angle_radian));
    }
  }

  for (PolyMesh::ConstEdgeIter ce_it = mesh.edges_begin();
      ce_it != mesh.edges_end(); ++ce_it) {
    PolyMesh::HalfedgeHandle heh = mesh.halfedge_handle(*ce_it, 0);
    int v1 = mesh.from_vertex_handle(heh).idx(), v2 = mesh.to_vertex_handle(heh)
        .idx();
    solver.add_constraint(
        new EdgeLengthConstraint<3>(v1, v2, edge_length_weight, edge_length));
  }

  if (solver.initialize(p.cols(), global_update_type, spd_solver)) {
    solver.solve(p, iter, Anderson_M);
    solver.output_iteration_history(AA_SOLVER);
    PolyMesh new_mesh = mesh;
    set_vertex_points(new_mesh, solver.get_points());
    check_wiremesh_error(new_mesh, edge_length, min_angle_radian,
                         max_angle_radian);
    check_ref_surface_distance(new_mesh, ref_mesh);

    if (output_filename) {
      if (!OpenMesh::IO::write_mesh(new_mesh, output_filename,
                                    OpenMesh::IO::Options::Default, 16)) {
        std::cerr << "Error: unable to save perturbed mesh to file "
                  << output_filename << std::endl;
      }
    }
  } else {
    std::cerr << "Error: unable to initialize solver" << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cout
        << "Usage:   <WireMeshOpt>  <INPUT_POLY_MESH>  <REF_TRI_MESH>  <OPTIONS_FILE>  <OUTPUT_MESH>"
        << std::endl;
    return 1;
  }

  PolyMesh mesh;
  if (!OpenMesh::IO::read_mesh(mesh, argv[1])) {
    std::cerr << "Error: unable to read input mesh from the file " << argv[1]
              << std::endl;
    return 1;
  }

  TriMesh ref_mesh;
  if (!OpenMesh::IO::read_mesh(ref_mesh, argv[2])) {
    std::cerr << "Error: unable to read referece mesh from file " << argv[2]
              << std::endl;
    return 1;
  }

  Scalar edge_length = average_edge_length(mesh);
  Scalar min_angle_radian = M_PI * 0.25, max_angle_radian = M_PI * 0.75;

  PolyMesh sub_mesh = subdivide_and_smooth_mesh(mesh);
  edge_length *= 0.5;
  check_wiremesh_error(sub_mesh, edge_length, min_angle_radian,
                       max_angle_radian);
  check_ref_surface_distance(sub_mesh, ref_mesh);

  if (argc == 5) {
    if (!OpenMesh::IO::write_mesh(sub_mesh, argv[3],
                                  OpenMesh::IO::Options::Default, 16)) {
      std::cerr << "Error: unable to save subdvided mesh to file " << argv[4]
                << std::endl;
    }
  }

  Parameters param;
  if (!param.load(argv[3])) {
    std::cerr << "Error: unable to load option file " << argv[3] << std::endl;
    return 1;
  }
  if (!param.valid_parameters()) {
    std::cerr << "Invalid filter options. Aborting..." << std::endl;
    return 1;
  }
  param.output();

  if (!param.acceleration)
    param.anderson_m = 0;

  Scalar edge_length_weight = 100, angle_weight = 100, closeness_weight = 1;

  std::cout << "========== Accelerate solver begin... ============="
            << std::endl;
  optimize_mesh(sub_mesh, ref_mesh, param.iter, param.global_update_type,
                param.anderson_m, LDLT_SOLVER, angle_weight, min_angle_radian,
                max_angle_radian, edge_length_weight, edge_length,
                closeness_weight, argv[4]);

  return 0;
}

