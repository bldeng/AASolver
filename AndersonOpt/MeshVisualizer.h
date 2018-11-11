//  BSD 3-Clause License
//
//  Copyright (c) 2018, Yue Peng
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

#ifndef MESHVISUALIZER_H_
#define MESHVISUALIZER_H_

#include "Types.h"

class MeshVisualizer {
 public:
  virtual ~MeshVisualizer() {
  }

  void render(bool smooth_shading) {
    update_normals(smooth_shading);

  }

 protected:
  virtual int n_surface_vertices() = 0;
  virtual int n_surface_faces() = 0;
  virtual const Matrix3X& get_surface_vertices() = 0;
  virtual const Eigen::Matrix3Xi& get_surface_connectivity() = 0;

  Matrix3X face_normals_, vertex_normals_;
  bool normal_arrays_initialized_;

  // triangle normal without normalization
  Vector3 triangle_normal(const Vector3 &p0, const Vector3 &p1,
                          const Vector3 &p2) {
    return (p2 - p1).cross(p0 - p1);
  }

  void update_normals(bool smooth_shading) {
    int n_vertices = n_surface_vertices(), n_faces = n_surface_faces();

    if (!normal_arrays_initialized_) {
      face_normals_.setZero(3, n_faces);
      vertex_normals_.setZero(3, n_vertices);
      normal_arrays_initialized_ = true;
    } else {
      face_normals_.setZero();
      vertex_normals_.setZero();
    }

    const Matrix3X& points = get_surface_vertices();
    const Eigen::Matrix3Xi& face_vtx_idx = get_surface_connectivity();

    for (int i = 0; i < n_faces; ++i) {
      Vector3 N = triangle_normal(points.col(face_vtx_idx(0, i)),
                                  points.col(face_vtx_idx(1, i)),
                                  points.col(face_vtx_idx(2, i)));
      face_normals_.col(i) = N;

      if (smooth_shading) {
        for (int k = 0; k < 3; ++k) {
          vertex_normals_.col(face_vtx_idx(k, i)) += N;
        }
      }
    }

    face_normals_.colwise().normalize();

    if (smooth_shading) {
      vertex_normals_.colwise().normalize();
    }
  }
};

#endif /* MESHVISUALIZER_H_ */
