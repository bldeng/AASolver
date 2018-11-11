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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Types.h"
#include "TriMeshAABB.h"
#include <igl/svd3x3.h>
#include <igl/AABB.h>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

#ifndef M_PI
/** Defining pi.*/
#define M_PI 3.14159265358979323846
#endif

template<unsigned int N>
class Constraint {
 protected:
  typedef MatrixT<N, 1> VectorN;
  typedef MatrixT<N, Eigen::Dynamic> MatrixNX;

  // The type of transform before projection
  enum InvarianceTransformType {
    MEAN_CENTERING,
    SUBTRACT_FIRST,
    IDENTITY,
  };

 public:
  Constraint(const std::vector<int> &idI, Scalar weight,
             InvarianceTransformType transform_type)
      : idI_(Eigen::Map<const Eigen::VectorXi>(idI.data(), idI.size())),
        weight_(std::sqrt(weight)),
        transform_type_(transform_type),
        idO_(-1) {
    transformed_points_.resize(N, num_transformed_points());
  }

  virtual ~Constraint() {
  }

  // Perform projection, and return the weighted squared distance from the transformed points to the projection points
  virtual Scalar project(const MatrixNX &positions, MatrixNX &projections) {
    if (transform_type_ == SUBTRACT_FIRST) {
      VectorN first_pt = positions.col(idI_[0]);

      for (int i = 1; i < static_cast<int>(idI_.size()); ++i) {
        transformed_points_.col(i - 1) = positions.col(idI_[i]) - first_pt;
      }
    } else {
      for (int i = 0; i < static_cast<int>(idI_.size()); ++i) {
        transformed_points_.col(i) = positions.col(idI_[i]);
      }

      if (transform_type_ == MEAN_CENTERING) {
        VectorN mean_vector = transformed_points_.rowwise().mean();
        transformed_points_.colwise() -= mean_vector;
      }
    }

    Eigen::Map<MatrixNX> proj_block(&(projections(0, idO_)), N,
                                    transformed_points_.cols());

    project_impl(transformed_points_, proj_block);

    Scalar squared_dist = (transformed_points_ - proj_block).squaredNorm();

    proj_block *= weight_;

    return squared_dist * (weight_ * weight_) * 0.5;
  }

  virtual void add_constraint(std::vector<Triplet> &triplets, int &idO) {
    idO_ = idO;
    int n_idx = static_cast<int>(idI_.size());

    if (transform_type_ == MEAN_CENTERING) {
      Scalar coef1 = (1.0 - 1.0 / n_idx) * weight_;
      Scalar coef2 = -weight_ / n_idx;
      for (int i = 0; i < n_idx; ++i) {
        for (int j = 0; j < n_idx; ++j) {
          triplets.push_back(Triplet(idO, idI_[j], (i == j ? coef1 : coef2)));
        }

        idO++;
      }
    } else if (transform_type_ == SUBTRACT_FIRST) {
      for (int i = 1; i < n_idx; ++i) {
        triplets.push_back(Triplet(idO, idI_[0], -weight_));
        triplets.push_back(Triplet(idO, idI_[i], weight_));
        idO++;
      }
    } else {
      for (int i = 0; i < n_idx; ++i) {
        triplets.push_back(Triplet(idO++, idI_[i], weight_));
      }
    }
  }

  int num_indices() const {
    return idI_.size();
  }

 protected:

  int num_transformed_points() const {
    if (transform_type_ == SUBTRACT_FIRST) {
      return num_indices() - 1;
    } else {
      return num_indices();
    }
  }

  Eigen::VectorXi idI_;

  Scalar weight_;

  InvarianceTransformType transform_type_;

  int idO_;

  MatrixNX transformed_points_;

  // Implementation of the projection, each subclass should override this method
  virtual void project_impl(const MatrixNX &transformed_pts,
                            Eigen::Map<MatrixNX> &projection) {
    projection = transformed_pts;
  }

  Scalar clamp(Scalar val, Scalar min_value, Scalar max_value) {
    return std::min(std::max(min_value, val), max_value);
  }
};

template<unsigned int N>
class EdgeLengthConstraint : public Constraint<N> {
 protected:
  using typename Constraint<N>::VectorN;
  using typename Constraint<N>::MatrixNX;

 public:
  EdgeLengthConstraint(int idx1, int idx2, Scalar weight, Scalar target_length)
      : Constraint<N>(std::vector<int>( { idx1, idx2 }), weight,
                      Constraint<N>::SUBTRACT_FIRST),
        target_length_(target_length) {
  }

  virtual ~EdgeLengthConstraint() {
  }

 protected:
  virtual void project_impl(const MatrixNX &transformed_pts,
                            Eigen::Map<MatrixNX> &projection) {
    projection.col(0) = transformed_pts.col(0).normalized() * target_length_;
  }

 private:
  Scalar target_length_;
};

template<unsigned int N>
class AngleConstraint : public Constraint<N> {
 protected:
  using typename Constraint<N>::VectorN;
  using typename Constraint<N>::MatrixNX;
  using Constraint<N>::clamp;

 public:
  AngleConstraint(int tip_idx, int side_idx1, int side_idx2, Scalar weight,
                  Scalar min_radian, Scalar max_radian)
      : Constraint<N>(std::vector<int>( { tip_idx, side_idx1, side_idx2 }),
                      weight, Constraint<N>::SUBTRACT_FIRST),
        min_angle_(std::max(Scalar(0), min_radian)),
        max_angle_(std::min(Scalar(M_PI), max_radian)) {
    min_angle_cos_ = clamp(std::cos(min_angle_), Scalar(-1), Scalar(1));
    max_angle_cos_ = clamp(std::cos(max_angle_), Scalar(-1), Scalar(1));
    assert(max_angle_cos_ <= min_angle_cos_);
  }

  virtual ~AngleConstraint() {
  }

 protected:
  virtual void project_impl(const MatrixNX &transformed_pts,
                            Eigen::Map<MatrixNX> &projection) {
    projection = transformed_pts;

    VectorN v1 = transformed_pts.col(0), v2 = transformed_pts.col(1);
    Scalar epsilon = 1e-14;
    Scalar v1_sqrnorm = v1.squaredNorm(), v2_sqrnorm = v2.squaredNorm();
    Scalar v1_norm = v1.norm(), v2_norm = v2.norm();
    VectorN unit_v1 = v1.normalized(), unit_v2 = v2.normalized();

    // cosine value of the angle gamma between v1 and v2
    Scalar cos_gamma = clamp(unit_v1.dot(unit_v2), Scalar(-1), Scalar(1));

    // Proceed only when the current angle lies outside the target range, and v1, v2 are not colinear
    if ((Scalar(1) - std::abs(cos_gamma) > epsilon)
        && (cos_gamma > min_angle_cos_ || cos_gamma < max_angle_cos_)) {
      Scalar gamma = std::acos(cos_gamma);

      // Angle eta: the sum of displacement angles from v1 and v2
      Scalar eta =
          cos_gamma > min_angle_cos_ ?
              (min_angle_ - gamma) : (gamma - max_angle_);
      eta = (std::max)(eta, 0.0);

      // Compute angle theta between v1 and its projection, and angle phi between v2 and its projection
      Scalar theta = Scalar(0.5)
          * std::atan2(v2_sqrnorm * std::sin(2 * eta),
                       v1_sqrnorm + v2_sqrnorm * std::cos(2 * eta));
      theta = (std::max)(0.0, (std::min)(eta, theta));
      Scalar phi = eta - theta;

      // Compute unit vectors that are coplanar with v1, v2, and orthogonal to one of them.
      // They form orthogonal frames with v1 and v2 respectively, within which we compute the projection using the above angles
      VectorN unit_v3 = (unit_v2 - unit_v1 * cos_gamma).normalized(), unit_v4 =
          (unit_v1 - unit_v2 * cos_gamma).normalized();

      // Determine if v1, v2 should move away from each other or towards each other
      if (cos_gamma > min_angle_cos_) {
        unit_v3 *= -1.0;
        unit_v4 *= -1.0;
      }

      projection.col(0) =
          (unit_v1 * std::cos(theta) + unit_v3 * std::sin(theta))
              * (v1_norm * std::cos(theta));
      projection.col(1) = (unit_v2 * std::cos(phi) + unit_v4 * std::sin(phi))
          * (v2_norm * std::cos(phi));
    }
  }

 private:
  Scalar min_angle_, max_angle_;
  Scalar min_angle_cos_, max_angle_cos_;
};

// Use this class for dynamic handle constraint: handle position is updated using set_target_pos
template<unsigned int N>
class ClosenessConstraint : public Constraint<N> {
 protected:
  using typename Constraint<N>::VectorN;

 public:
  ClosenessConstraint(int idx, Scalar weight, const VectorN &target_pos)
      : Constraint<N>(std::vector<int>( { idx }), weight,
                      Constraint<N>::IDENTITY),
        target_pos_(target_pos) {
  }

  virtual ~ClosenessConstraint() {
  }

  void set_target_pos(const VectorN &target_pos) {
    target_pos_ = target_pos;
  }

 protected:
  virtual void proj_impl(const Matrix3X&, Eigen::Map<Matrix3X> &projection) {
    projection.col(0) = target_pos_;
  }

 private:
  VectorN target_pos_;
};

class PointToRefSurfaceConstraint : public Constraint<3> {
 public:
  PointToRefSurfaceConstraint(int pt_idx, Scalar weight,
                              const std::shared_ptr<TriMeshAABB> &aabb)
      : Constraint(std::vector<int>( { pt_idx }), weight, IDENTITY),
        AABB_tree_(aabb) {
  }

  virtual ~PointToRefSurfaceConstraint() {
  }

 protected:
  virtual void project_impl(const Matrix3X &transformed_pts,
                            Eigen::Map<Matrix3X> &projection) {
    Vector3 closest_point;
    AABB_tree_->get_closest_point(transformed_pts.col(0), closest_point);
    projection.col(0) = closest_point;
  }

 private:
  std::shared_ptr<TriMeshAABB> AABB_tree_;
};

class ReferenceSurfceConstraint : public Constraint<3> {
 public:
  ReferenceSurfceConstraint(int n_points, Scalar weight,
                            const Matrix3X &ref_surface_vtx,
                            const Eigen::Matrix3Xi &ref_surface_faces)
      : Constraint(std::vector<int>(), weight, IDENTITY) {
    n_points_ = n_points;

    idI_.resize(n_points);
    for (int i = 0; i < n_points; ++i) {
      idI_(i) = i;
    }

    sqrD_.setZero(n_points);
    I_.setZero(n_points);
    C_.setZero(n_points, 3);

    ref_surface_vtx_ = ref_surface_vtx.transpose();
    ref_surface_faces_ = ref_surface_faces.transpose();
    AABB_tree_.init(ref_surface_vtx_, ref_surface_faces_);
  }

  virtual ~ReferenceSurfceConstraint() {
  }

  // Perform projection, and return the weighted squared distance from the transformed points to the projection points
  virtual Scalar project(const Matrix3X &positions, Matrix3X &projections) {
    pos_ = positions.transpose();
    AABB_tree_.squared_distance(ref_surface_vtx_, ref_surface_faces_, pos_,
                                sqrD_, I_, C_);
    projections.block(0, idO_, 3, n_points_) = C_.transpose() * weight_;
    return sqrD_.sum() * (weight_ * weight_) * 0.5;
  }

 private:
  int n_points_;
  MatrixX3 pos_;
  MatrixX3 ref_surface_vtx_;
  Eigen::MatrixX3i ref_surface_faces_;
  VectorX sqrD_;
  Eigen::VectorXi I_;
  MatrixX3 C_;
  igl::AABB<MatrixX3, 3> AABB_tree_;
};

class PlaneConstraint : public Constraint<3> {
 public:
  PlaneConstraint(const std::vector<int> &idI, Scalar weight)
      : Constraint(idI, weight, MEAN_CENTERING) {
  }

  virtual ~PlaneConstraint() {
  }

 protected:
  virtual void project_impl(const Matrix3X &transformed_pts,
                            Eigen::Map<Matrix3X> &projection) {
    Eigen::JacobiSVD<Matrix3X> jSVD;
    jSVD.compute(transformed_pts, Eigen::ComputeFullU);
    Vector3 best_fit_normal = jSVD.matrixU().col(2).normalized();
    projection = transformed_pts
        - best_fit_normal * (best_fit_normal.transpose() * transformed_pts);
  }
};

class ARAP2DTriangleConstraint : public Constraint<2> {
 public:
  ARAP2DTriangleConstraint(const std::vector<int> &idI, Scalar weight,
                           const Matrix2X& positions)
      : Constraint(idI, weight, MEAN_CENTERING) {
    assert(idI.size() == 3);
    Matrix23 points;
    for (int i = 0; i < 3; ++i) {
      points.col(i) = positions.col(idI[i]);
    }
    points.colwise() -= points.rowwise().mean().eval();

    Matrix22 edge_vecs;
    edge_vecs.col(0) = points.col(1) - points.col(0);
    edge_vecs.col(1) = points.col(2) - points.col(0);
    Scalar area = 0.5 * std::fabs(edge_vecs.determinant());
    weight_ *= std::sqrt(area);

    // Transformation matrix
    transform_ = points.transpose().jacobiSvd(
        Eigen::ComputeFullU | Eigen::ComputeFullV).solve(
        Matrix33::Identity() - Matrix33::Constant(Scalar(1) / 3.0)).transpose();
  }

  virtual void add_constraint(std::vector<Triplet> &triplets, int &idO) {
    idO_ = idO;

    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 3; ++j) {
        triplets.push_back(Triplet(idO, idI_[j], weight_ * transform_(j, i)));
      }

      idO++;
    }
  }

  virtual Scalar project(const Matrix2X &positions, Matrix2X &projections) {
    for (int i = 0; i < 3; ++i) {
      transformed_points_.col(i) = positions.col(idI_[i]);
    }

    Matrix22 M = transformed_points_ * transform_;
    Vector2 b(M(0, 0) + M(1, 1), M(0, 1) - M(1, 0));
    b.normalize();
    Matrix22 P;
    P(0, 0) = P(1, 1) = b(0);
    P(0, 1) = b(1);
    P(1, 0) = -b(1);
    Scalar weighted_dist = (M - P).squaredNorm() * (weight_ * weight_)
        * Scalar(0.5);
    projections.block(0, idO_, 2, 2) = weight_ * P;

    return weighted_dist;
  }

 private:
  // Use unaligned Eigen matrix to avoid alignment issues
  typedef Eigen::Matrix<Scalar, 3, 2, Eigen::ColMajor | Eigen::DontAlign> UnAlignedMatrix32;
  UnAlignedMatrix32 transform_;  // Transform of the mean-centered vertex positions to compute deformation gradient
};

class ARAP3DTetConstraint : public Constraint<3> {
 public:
  ARAP3DTetConstraint(const std::vector<int> &idI, Scalar weight,
                      const Matrix3X& positions, bool fast_svd)
      : Constraint(idI, weight, MEAN_CENTERING),
        use_fast_svd_(fast_svd) {
    assert(idI.size() == 4);
    Matrix34 points;
    for (int i = 0; i < 4; ++i) {
      points.col(i) = positions.col(idI[i]);
    }
    points.colwise() -= points.rowwise().mean().eval();

    Matrix33 edge_vecs;
    for (int i = 0; i < 3; ++i) {
      edge_vecs.col(i) = points.col(i + 1) - points.col(0);
    }

    volume_ = std::fabs(edge_vecs.determinant()) / 6.0;
    weight_ *= std::sqrt(volume_);

    // Transformation matrix
    transform_ = points.transpose().jacobiSvd(
        Eigen::ComputeFullU | Eigen::ComputeFullV).solve(
        Matrix44::Identity() - Matrix44::Constant(Scalar(1) / 4.0)).transpose();
  }

  virtual void add_constraint(std::vector<Triplet> &triplets, int &idO) {
    idO_ = idO;

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 4; ++j) {
        triplets.push_back(Triplet(idO, idI_[j], weight_ * transform_(j, i)));
      }

      idO++;
    }
  }

  virtual Scalar project(const Matrix3X &positions, Matrix3X &projections) {
    for (int i = 0; i < 4; ++i) {
      transformed_points_.col(i) = positions.col(idI_[i]);
    }

    Matrix33 M = transformed_points_ * transform_;
    Matrix33 P;

    if (use_fast_svd_) {
      Eigen::Matrix<float, 3, 3> U, V, A, R;

      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          A(i, j) = M(i, j);
        }
      }

      Eigen::Matrix<float, 3, 1> S;
      igl::svd3x3(A, U, S, V);
      R = U * V.transpose();

      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          P(i, j) = R(i, j);
        }
      }
    } else {
      Eigen::JacobiSVD<Matrix33> jSVD(
          M, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Matrix33 I = Matrix33::Identity();
      if ((jSVD.matrixU() * jSVD.matrixV().transpose()).determinant() < 0) {
        I(2, 2) = -1;
      }
      P = jSVD.matrixU() * I * jSVD.matrixV().transpose();
    }

    Scalar weighted_dist = (M - P).squaredNorm() * (weight_ * weight_)
        * Scalar(0.5);
    projections.block(0, idO_, 3, 3) = weight_ * P;
    return weighted_dist;
  }

  Scalar get_volume() const {
    return volume_;
  }

 private:
  // Use unaligned Eigen matrix to avoid alignment issues
  typedef Eigen::Matrix<Scalar, 4, 3, Eigen::ColMajor | Eigen::DontAlign> UnAlignedMatrix43;
  UnAlignedMatrix43 transform_;  // Transform of the mean-centered vertex positions to compute deformation gradient
  Scalar volume_;		// Volume of the original tetrahedron
  bool use_fast_svd_;
};

class TetStrainConstraint : public Constraint<3> {
 public:
  TetStrainConstraint(const std::vector<int> &idI, Scalar weight,
                      const Matrix3X& positions, bool fast_svd)
      : Constraint(idI, weight, IDENTITY),
        use_fast_svd_(fast_svd) {
    assert(idI.size() == 4);

    Matrix33 edge_vecs;
    for (int i = 0; i < 3; ++i) {
      edge_vecs.col(i) = positions.col(idI[i + 1]) - positions.col(idI[0]);
    }

    Matrix34 T;
    T.setZero();
    for (int i = 0; i < 3; ++i) {
      T(i, 0) = -1;
      T(i, i + 1) = 1;
    }

    volume_ = std::fabs(edge_vecs.determinant()) / 6.0;
    weight_ *= std::sqrt(volume_);

    // Transformation matrix
    transform_ = edge_vecs.transpose().jacobiSvd(
        Eigen::ComputeFullU | Eigen::ComputeFullV).solve(T).transpose();
  }

  virtual void add_constraint(std::vector<Triplet> &triplets, int &idO) {
    idO_ = idO;

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 4; ++j) {
        triplets.push_back(Triplet(idO, idI_[j], weight_ * transform_(j, i)));
      }

      idO++;
    }
  }

  virtual Scalar project(const Matrix3X &positions, Matrix3X &projections) {
    for (int i = 0; i < 4; ++i) {
      transformed_points_.col(i) = positions.col(idI_[i]);
    }

    Matrix33 M = transformed_points_ * transform_;
    Matrix33 P;

    if (use_fast_svd_) {
      Eigen::Matrix<float, 3, 3> U, V, A, R;

      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          A(i, j) = M(i, j);
        }
      }

      Eigen::Matrix<float, 3, 1> S;
      igl::svd3x3(A, U, S, V);
      R = U * V.transpose();

      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          P(i, j) = R(i, j);
        }
      }
    } else {
      Eigen::JacobiSVD<Matrix33> jSVD(
          M, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Matrix33 I = Matrix33::Identity();
      if ((jSVD.matrixU() * jSVD.matrixV().transpose()).determinant() < 0) {
        I(2, 2) = -1;
      }
      P = jSVD.matrixU() * I * jSVD.matrixV().transpose();
    }

    Scalar weighted_dist = (M - P).squaredNorm() * (weight_ * weight_)
        * Scalar(0.5);
    projections.block(0, idO_, 3, 3) = weight_ * P;
    //std::cout << "proj_block = " << P << std::endl;
    return weighted_dist;
  }

  Scalar get_volume() const {
    return volume_;
  }
 private:
  // Use unaligned Eigen matrix to avoid alignment issues
  typedef Eigen::Matrix<Scalar, 4, 3, Eigen::ColMajor | Eigen::DontAlign> UnAlignedMatrix43;
  UnAlignedMatrix43 transform_;  // Transform of the mean-centered vertex positions to compute deformation gradient
  Scalar volume_;		// Volume of the original tetrahedron
  bool use_fast_svd_;
};
#endif // CONSTRAINT_H
