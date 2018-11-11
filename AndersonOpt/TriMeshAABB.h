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

#ifndef TRIMESHAABB_H_
#define TRIMESHAABB_H_

#include "Types.h"
#include "MeshTypes.h"
#include <igl/AABB.h>

class TriMeshAABB {
 public:
  TriMeshAABB(const TriMesh &ref_mesh) {
    Matrix3X ref_mesh_points;
    Eigen::Matrix3Xi ref_mesh_faces;
    get_vertex_points(ref_mesh, ref_mesh_points);
    get_face_vertex_index(ref_mesh, ref_mesh_faces);
    V_ref = ref_mesh_points.transpose();
    F_ref = ref_mesh_faces.transpose();

    F_normals_.resize(3, ref_mesh.n_faces());
    for (TriMesh::ConstFaceIter cf_it = ref_mesh.faces_begin();
        cf_it != ref_mesh.faces_end(); ++cf_it) {
      F_normals_.col(cf_it->idx()) = to_eigen_vec3(
          ref_mesh.calc_face_normal(*cf_it)).normalized();
    }

    AABB_tree_.init(V_ref, F_ref);
  }

  Scalar get_closest_point(const Vector3 &point, Vector3 &closest_point,
                           Vector3* closest_point_normal = NULL) const {
    int i;
    igl::AABB<MatrixX3, 3>::RowVectorDIMS p = point.transpose(), c;
    Scalar sqr_dist = AABB_tree_.squared_distance(V_ref, F_ref, p, i, c);
    closest_point = c.transpose();

    if (closest_point_normal) {
      *closest_point_normal = F_normals_.col(i);
    }

    return sqr_dist;
  }

 private:
  MatrixX3 V_ref;
  Eigen::MatrixX3i F_ref;
  Matrix3X F_normals_;
  igl::AABB<MatrixX3, 3> AABB_tree_;
};

#endif /* TRIMESHAABB_H_ */
