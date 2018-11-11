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

#ifndef LINEARREGULARIZATION_H_
#define LINEARREGULARIZATION_H_

#include "Types.h"
#include <vector>

template<unsigned int N>
class LinearRegularization {
 protected:
  typedef MatrixT<N, 1> VectorN;
  typedef MatrixT<N, Eigen::Dynamic> MatrixNX;
  typedef MatrixT<Eigen::Dynamic, N> MatrixXN;

 public:
  // Add Laplacian regularization.
  // The first index in the index array corresponds to the central point
  void add_uniform_laplacian(const std::vector<int> &indices, Scalar weight) {
    int n_pts = indices.size();
    std::vector<Scalar> coefs;
    coefs.push_back(Scalar(1));
    coefs.insert(coefs.end(), n_pts - 1, Scalar(-1.0 / double(n_pts - 1)));
    add_laplacian_helper(indices, coefs, weight);
  }

  void add_laplacian(const std::vector<int> &indices,
                     const std::vector<Scalar> coefs, Scalar weight) {
    add_laplacian_helper(indices, coefs, weight);
  }

  void add_relative_uniform_laplacian(const std::vector<int> &indices,
                                      Scalar weight,
                                      const MatrixNX &ref_points) {
    int n_pts = indices.size();
    std::vector<Scalar> coefs;
    coefs.push_back(Scalar(1));
    coefs.insert(coefs.end(), n_pts - 1, Scalar(-1.0 / double(n_pts - 1)));
    add_laplacian_helper(indices, coefs, weight, &ref_points);
  }

  void add_relative_laplacian(const std::vector<int> &indices,
                              const std::vector<Scalar> coefs, Scalar weight,
                              const MatrixNX &ref_points) {
    add_laplacian_helper(indices, coefs, weight, &ref_points);
  }

  void add_closeness(int idx, Scalar weight, const VectorN &target_pt) {
    Scalar sqrt_weight(std::sqrt(double(weight)));

    Eigen::VectorXi index_vec(1);
    index_vec(0) = idx;
    point_indices_.push_back(index_vec);

    VectorX coef_vec(1);
    coef_vec(0) = sqrt_weight;
    coefficients_.push_back(coef_vec);

    add_to_target_values(target_pt * sqrt_weight);
  }

  // Set up the linear system AX = b for the regularization
  bool get_regularization_system(int n_pts, ColMajorSparseMatrix &L,
                                 MatrixXN &rhs) {
    int n_rows = point_indices_.size();
    rhs.setZero(n_rows, N);
    L.resize(n_rows, n_pts);

    if (n_rows == 0) {
      return false;
    }

    std::vector<Triplet> triplets;
    for (int i = 0; i < n_rows; ++i) {
      int n_coefs = point_indices_[i].size();
      for (int j = 0; j < n_coefs; ++j) {
        triplets.push_back(
            Triplet(i, point_indices_[i](j), coefficients_[i](j)));
      }

      rhs.row(i) = Eigen::Map < VectorX
          > (&(target_values_[i * N]), N).transpose();
    }

    L.setFromTriplets(triplets.begin(), triplets.end());
    L.makeCompressed();

    return true;
  }

 private:
  std::vector<Eigen::VectorXi> point_indices_;	// Indices of points subject to regularization constraints
  std::vector<VectorX> coefficients_;  // Linear combination coefficients
  std::vector<Scalar> target_values_;	// Target value of the linear combination, flatten into a single array to avoid alignment issue with std::vector<Vector2>

  void add_laplacian_helper(const std::vector<int> &indices,
                            const std::vector<Scalar> coefs, Scalar weight,
                            const MatrixNX *ref_points = NULL) {
    assert(indices.size() == coefs.size());

    Scalar sqrt_weight(std::sqrt(double(weight)));

    Eigen::VectorXi index_vec = Eigen::Map<const Eigen::VectorXi>(
        indices.data(), indices.size());
    point_indices_.push_back(index_vec);

    VectorX coef_vec = Eigen::Map<const VectorX>(coefs.data(), coefs.size())
        * sqrt_weight;
    coefficients_.push_back(coef_vec);

    VectorN target_value = VectorN::Zero();
    if (ref_points) {
      for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
        target_value += ref_points->col(indices[i]) * coefs[i];
      }
    }

    add_to_target_values(target_value * sqrt_weight);
  }

  void add_to_target_values(const VectorN &val) {
    for (unsigned int i = 0; i < N; ++i) {
      target_values_.push_back(val(i));
    }
  }
};

#endif /* LINEARREGULARIZATION_H_ */
