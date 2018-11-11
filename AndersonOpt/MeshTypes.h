//  BSD 3-Clause License
//
//  Copyright (c) 2018, Bailin Deng
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

#ifndef MESHTYPES_H
#define MESHTYPES_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <string>
#include "Types.h"
#include "LinearRegularization.h"

/// Default traits for polygonal mesh
struct PolyTraits : public OpenMesh::DefaultTraits {
  /// Use double precision points
  typedef OpenMesh::Vec3d Point;
  /// Use double precision Normals
  typedef OpenMesh::Vec3d Normal;

  /// Use RGBA Color
  typedef OpenMesh::Vec4f Color;
};

/// Simple Name for Mesh
typedef OpenMesh::PolyMesh_ArrayKernelT<PolyTraits> PolyMesh;

/// Default traits for triangle mesh
struct TriTraits : public OpenMesh::DefaultTraits {
  /// Use double precision points
  typedef OpenMesh::Vec3d Point;
  /// Use double precision Normals
  typedef OpenMesh::Vec3d Normal;

  /// Use RGBA Color
  typedef OpenMesh::Vec4f Color;

};

/// Simple name for triangle mesh
typedef OpenMesh::TriMesh_ArrayKernelT<TriTraits> TriMesh;

// Collect all point positions into a matrix.
// Argument points must be have been initialized with the correct number of columns.
template<typename MeshT>
void get_vertex_points(const MeshT &mesh, Matrix3X &points) {
  points.resize(3, mesh.n_vertices());
  for (typename MeshT::ConstVertexIter cv_it = mesh.vertices_begin();
      cv_it != mesh.vertices_end(); ++cv_it) {
    points.col(cv_it->idx()) = to_eigen_vec3(mesh.point(*cv_it));
  }
}

template<typename MeshT>
void set_vertex_points(MeshT &mesh, const Matrix3X &pos) {
  assert(
      static_cast<int>(mesh.n_vertices()) * 3 == static_cast<int>(pos.size()));

  for (typename MeshT::ConstVertexIter cv_it = mesh.vertices_begin();
      cv_it != mesh.vertices_end(); ++cv_it) {
    Vector3 pt = pos.col(cv_it->idx());
    mesh.set_point(*cv_it, from_eigen_vec3<typename MeshT::Point>(pt));
  }
}

template<typename MeshT>
void get_vertex_points(const MeshT &mesh, std::vector<double> &points) {
  points.assign(3 * mesh.n_vertices(), 0.0);
  for (typename MeshT::ConstVertexIter cv_it = mesh.vertices_begin();
      cv_it != mesh.vertices_end(); ++cv_it) {
    typename MeshT::Point pt = mesh.point(*cv_it);
    int v_idx = cv_it->idx();

    assert(v_idx >= 0 && v_idx < static_cast<int>(mesh.n_vertices()));

    for (int i = 0; i < 3; ++i) {
      points[3 * v_idx + i] = pt[i];
    }
  }
}

template<typename MeshT>
void set_vertex_points(MeshT &mesh, const std::vector<double> &pos) {
  assert(mesh.n_vertices() * 3 == pos.size());

  for (typename MeshT::ConstVertexIter cv_it = mesh.vertices_begin();
      cv_it != mesh.vertices_end(); ++cv_it) {
    int addr = cv_it->idx() * 3;
    mesh.set_point(
        *cv_it, typename MeshT::Point(pos[addr], pos[addr + 1], pos[addr + 2]));
  }
}

// Write a mesh to an ASCII file with high accuracy
template<typename MeshT>
inline bool write_mesh_high_accuracy(const MeshT &mesh,
                                     const std::string &filename) {
  return OpenMesh::IO::write_mesh(mesh, filename,
                                  OpenMesh::IO::Options::Default, 16);
}

template<typename MeshT>
inline Vector3 bbox_dimension(const MeshT &mesh) {
  Matrix3X vtx_pos;
  get_vertex_points(mesh, vtx_pos);

  return vtx_pos.rowwise().maxCoeff() - vtx_pos.rowwise().minCoeff();
}

template<typename MeshT>
inline double bbox_diag_length(const MeshT &mesh) {
  return bbox_dimension(mesh).norm();
}

template<typename MeshT>
inline double average_edge_length(const MeshT &mesh) {
  if (static_cast<int>(mesh.n_edges()) == 0) {
    return 0.0;
  }

  double length = 0;

  for (typename MeshT::ConstEdgeIter ce_it = mesh.edges_begin();
      ce_it != mesh.edges_end(); ++ce_it) {
    length += mesh.calc_edge_length(*ce_it);
  }

  return length / mesh.n_edges();
}

inline PolyMesh quad_subdivision(const PolyMesh &input_mesh) {
  PolyMesh out_mesh;

  std::vector<PolyMesh::VertexHandle> vertex_vhs(input_mesh.n_vertices()),
      edge_vhs(input_mesh.n_edges());

  for (PolyMesh::ConstVertexIter cv_it = input_mesh.vertices_begin();
      cv_it != input_mesh.vertices_end(); ++cv_it) {
    vertex_vhs[cv_it->idx()] = out_mesh.add_vertex(input_mesh.point(*cv_it));
  }

  for (PolyMesh::ConstEdgeIter ce_it = input_mesh.edges_begin();
      ce_it != input_mesh.edges_end(); ++ce_it) {
    PolyMesh::HalfedgeHandle heh = input_mesh.halfedge_handle(*ce_it, 0);
    PolyMesh::Point mid_pt = (input_mesh.point(
        input_mesh.from_vertex_handle(heh))
        + input_mesh.point(input_mesh.to_vertex_handle(heh))) * 0.5;
    edge_vhs[ce_it->idx()] = out_mesh.add_vertex(mid_pt);
  }

  for (PolyMesh::ConstFaceIter cf_it = input_mesh.faces_begin();
      cf_it != input_mesh.faces_end(); ++cf_it) {
    PolyMesh::Point pt(0, 0, 0);
    int n_pts = 0;
    for (PolyMesh::ConstFaceVertexIter cfv_it = input_mesh.cfv_iter(*cf_it);
        cfv_it.is_valid(); ++cfv_it) {
      pt += input_mesh.point(*cfv_it);
      n_pts++;
    }

    PolyMesh::VertexHandle face_vh = out_mesh.add_vertex(
        pt / PolyMesh::Scalar(n_pts));

    std::vector<PolyMesh::VertexHandle> current_edge_vhs, current_vtx_vhs;
    for (PolyMesh::ConstFaceHalfedgeIter cfh_it = input_mesh.cfh_iter(*cf_it);
        cfh_it.is_valid(); ++cfh_it) {
      current_edge_vhs.push_back(
          edge_vhs[input_mesh.edge_handle(*cfh_it).idx()]);
      current_vtx_vhs.push_back(
          vertex_vhs[input_mesh.to_vertex_handle(*cfh_it).idx()]);
    }

    for (int i = 0; i < n_pts; ++i) {
      std::vector<PolyMesh::VertexHandle> new_face;
      new_face.push_back(current_edge_vhs[i]);
      new_face.push_back(current_vtx_vhs[i]);
      new_face.push_back(current_edge_vhs[(i + 1) % n_pts]);
      new_face.push_back(face_vh);

      out_mesh.add_face(new_face);
    }
  }

  return out_mesh;
}

inline PolyMesh subdivide_and_smooth_mesh(const PolyMesh &input_mesh) {
  PolyMesh out_mesh;

  std::vector<PolyMesh::VertexHandle> vertex_vhs(input_mesh.n_vertices()),
      edge_vhs(input_mesh.n_edges());

  for (PolyMesh::ConstVertexIter cv_it = input_mesh.vertices_begin();
      cv_it != input_mesh.vertices_end(); ++cv_it) {
    vertex_vhs[cv_it->idx()] = out_mesh.add_vertex(input_mesh.point(*cv_it));
  }

  for (PolyMesh::ConstEdgeIter ce_it = input_mesh.edges_begin();
      ce_it != input_mesh.edges_end(); ++ce_it) {
    PolyMesh::HalfedgeHandle heh = input_mesh.halfedge_handle(*ce_it, 0);
    PolyMesh::Point mid_pt = (input_mesh.point(
        input_mesh.from_vertex_handle(heh))
        + input_mesh.point(input_mesh.to_vertex_handle(heh))) * 0.5;
    edge_vhs[ce_it->idx()] = out_mesh.add_vertex(mid_pt);
  }

  for (PolyMesh::ConstFaceIter cf_it = input_mesh.faces_begin();
      cf_it != input_mesh.faces_end(); ++cf_it) {
    PolyMesh::Point pt(0, 0, 0);
    int n_pts = 0;
    for (PolyMesh::ConstFaceVertexIter cfv_it = input_mesh.cfv_iter(*cf_it);
        cfv_it.is_valid(); ++cfv_it) {
      pt += input_mesh.point(*cfv_it);
      n_pts++;
    }

    PolyMesh::VertexHandle face_vh = out_mesh.add_vertex(
        pt / PolyMesh::Scalar(n_pts));

    std::vector<PolyMesh::VertexHandle> current_edge_vhs, current_vtx_vhs;
    for (PolyMesh::ConstFaceHalfedgeIter cfh_it = input_mesh.cfh_iter(*cf_it);
        cfh_it.is_valid(); ++cfh_it) {
      current_edge_vhs.push_back(
          edge_vhs[input_mesh.edge_handle(*cfh_it).idx()]);
      current_vtx_vhs.push_back(
          vertex_vhs[input_mesh.to_vertex_handle(*cfh_it).idx()]);
    }

    for (int i = 0; i < n_pts; ++i) {
      std::vector<PolyMesh::VertexHandle> new_face;
      new_face.push_back(current_edge_vhs[i]);
      new_face.push_back(current_vtx_vhs[i]);
      new_face.push_back(current_edge_vhs[(i + 1) % n_pts]);
      new_face.push_back(face_vh);

      out_mesh.add_face(new_face);
    }
  }

  // Construct Laplacian system
  LinearRegularization<3> regularizer;
  for (PolyMesh::ConstVertexIter cv_it = out_mesh.vertices_begin();
      cv_it != out_mesh.vertices_end(); ++cv_it) {
    std::vector<int> lap_vhs;
    lap_vhs.push_back(cv_it->idx());

    if (out_mesh.is_boundary(*cv_it)) {
      std::vector<PolyMesh::FaceHandle> fhs;
      for (PolyMesh::ConstVertexOHalfedgeIter cvoh_it = out_mesh.cvoh_iter(
          *cv_it); cvoh_it.is_valid(); ++cvoh_it) {
        PolyMesh::EdgeHandle eh = out_mesh.edge_handle(*cvoh_it);
        if (out_mesh.is_boundary(eh)) {
          lap_vhs.push_back(out_mesh.to_vertex_handle(*cvoh_it).idx());
          for (int k = 0; k < 2; ++k) {
            PolyMesh::HalfedgeHandle cur_heh = out_mesh.halfedge_handle(eh, k);
            if (!out_mesh.is_boundary(cur_heh)) {
              fhs.push_back(out_mesh.face_handle(cur_heh));
            }
          }
        }
      }

      if (static_cast<int>(fhs.size()) == 2 && fhs[0] != fhs[1]) {
        regularizer.add_uniform_laplacian(lap_vhs, 1.0);
      }
    } else {
      for (PolyMesh::ConstVertexVertexIter cvv_it = out_mesh.cvv_iter(*cv_it);
          cvv_it.is_valid(); ++cvv_it) {
        lap_vhs.push_back(cvv_it->idx());
      }

      regularizer.add_uniform_laplacian(lap_vhs, 1.0);
    }
  }

  int n_vtx = out_mesh.n_vertices();
  ColMajorSparseMatrix L;
  MatrixX3 rhs;
  if (!regularizer.get_regularization_system(n_vtx, L, rhs)) {
    std::cerr << "Error: unable to construct regularization system"
        << std::endl;
    return out_mesh;
  }

  MatrixX3 fixed_pt_coords(n_vtx, 3);
  fixed_pt_coords.setZero();
  std::vector<bool> is_fixed(n_vtx, false);
  for (int i = 0; i < static_cast<int>(vertex_vhs.size()); ++i) {
    int v_idx = vertex_vhs[i].idx();
    is_fixed[v_idx] = true;
    fixed_pt_coords.row(v_idx) = to_eigen_vec3(out_mesh.point(vertex_vhs[i]))
        .transpose();
  }

  std::vector<Triplet> var_selection_triplets;
  int var_row = 0;
  for (int i = 0; i < n_vtx; ++i) {
    if (!is_fixed[i]) {
      var_selection_triplets.push_back(Triplet(var_row++, i, Scalar(1)));
    }
  }
  ColMajorSparseMatrix var_selection(var_row, n_vtx);
  var_selection.setFromTriplets(var_selection_triplets.begin(),
                                var_selection_triplets.end());
  ColMajorSparseMatrix A = L * var_selection.transpose();
  ColMajorSparseMatrix M = A.transpose() * A;
  MatrixX3 reduced_rhs = -(A.transpose() * L) * fixed_pt_coords;
  Eigen::SimplicialLDLT<ColMajorSparseMatrix> ldlt_solver(M);
  MatrixX3 sol = ldlt_solver.solve(reduced_rhs);
  Matrix3X opt_pos = (var_selection.transpose() * sol + fixed_pt_coords)
      .transpose();
  set_vertex_points(out_mesh, opt_pos);

  return out_mesh;
}

void get_face_vertex_index(const TriMesh &mesh,
                           Eigen::Matrix3Xi &face_vertex_idx) {
  face_vertex_idx.resize(3, mesh.n_faces());
  for (TriMesh::ConstFaceIter cf_it = mesh.faces_begin();
      cf_it != mesh.faces_end(); ++cf_it) {
    std::vector<int> vtx;
    for (TriMesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(*cf_it);
        cfv_it.is_valid(); ++cfv_it) {
      vtx.push_back(cfv_it->idx());
    }

    face_vertex_idx.col(cf_it->idx()) = Eigen::Map < Eigen::Vector3i
        > (vtx.data(), 3);
  }
}

#endif // MESHTYPES_H
