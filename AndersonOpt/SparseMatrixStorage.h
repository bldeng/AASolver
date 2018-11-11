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


#ifndef SPARSEMATRIXSTORAGE_H_
#define SPARSEMATRIXSTORAGE_H_

#include "Types.h"
#include <cassert>

// Contiguous row-wise storage of the element in a sparse matrix
class RowSparseStorage {
 public:
  RowSparseStorage()
      : n_rows_(0) {
  }

  void init_from_transposed_matrix(const ColMajorSparseMatrix &mat_transposed) {
    elements_.clear();
    row_addr_.clear();
    n_rows_ = mat_transposed.cols();

    for (int k = 0; k < mat_transposed.outerSize(); ++k) {
      row_addr_.push_back(elements_.size());

      for (ColMajorSparseMatrix::InnerIterator it(mat_transposed, k); it;
          ++it) {
        elements_.push_back(std::make_pair(int(it.row()), Scalar(it.value())));
      }
    }

    row_addr_.push_back(elements_.size());
  }

  void clear() {
    elements_.clear();
    row_addr_.clear();
    n_rows_ = 0;
  }

  template<typename MatrixType, typename VectorType>
  void evaluate(int row_idx, const MatrixType &dense_mat, VectorType &res) {
    int addr_begin = row_addr_[row_idx], addr_end = row_addr_[row_idx + 1];
    for (int i = addr_begin; i < addr_end; ++i) {
      std::pair<int, Scalar> &current_entry = elements_[i];
      res += dense_mat.col(current_entry.first) * current_entry.second;
    }
  }

  int n_rows() {
    return n_rows_;
  }

 private:
  std::vector<std::pair<int, Scalar> > elements_;
  std::vector<int> row_addr_;
  int n_rows_;
};

// Row-wise storage of the diagonally dominant matrix for Jacobi iteration
class JacobiRowSparseStorage {
 public:
  JacobiRowSparseStorage()
      : n_rows_(0) {
  }

  bool init_from_transposed_matrix(const ColMajorSparseMatrix &mat_transposed) {
    elements_.clear();
    row_addr_.clear();
    n_rows_ = mat_transposed.cols();

    for (int k = 0; k < mat_transposed.outerSize(); ++k) {
      row_addr_.push_back(elements_.size());

      Scalar diag_value = 1.0;
      Scalar off_diag_values = 0;
      bool has_diag = false;

      for (ColMajorSparseMatrix::InnerIterator it(mat_transposed, k); it;
          ++it) {
        if (it.row() == it.col()) {
          diag_value = it.value();
          has_diag = true;
        } else {
          elements_.push_back(
              std::make_pair(int(it.row()), Scalar(it.value())));
          off_diag_values += std::fabs(it.value());
        }
      }

      if (!has_diag) {
        std::cerr << "Error: zero diagonal element for Jacobi system matrix"
            << std::endl;
        return false;
      }

      /*
       if(diag_value <= off_diag_values){
       std::cout << "Warning: the matrix is not diagonally dominant on row " <<
       k <<  ": diag " << diag_value << ", off-diag values " << off_diag_values << std::endl;
       }
       */

      elements_.push_back(std::pair<int, Scalar>(-1, diag_value));
    }

    row_addr_.push_back(elements_.size());

    return true;
  }

  void clear() {
    elements_.clear();
    row_addr_.clear();
    n_rows_ = 0;
  }

  template<typename MatrixType, typename VectorType>
  void evaluate(int row_idx, const MatrixType &X, VectorType &res) {
    int addr_begin = row_addr_[row_idx], addr_end = row_addr_[row_idx + 1] - 1;
    for (int i = addr_begin; i < addr_end; ++i) {
      std::pair<int, Scalar> &current_entry = elements_[i];
      res -= X.col(current_entry.first) * current_entry.second;
    }

    res /= elements_[addr_end].second;
  }

  int n_rows() {
    return n_rows_;
  }

  Scalar get_diagonal(int row_idx) {
    return elements_[row_addr_[row_idx + 1] - 1].second;
  }

  void set_diagonal(int row_idx, Scalar value) {
    elements_[row_addr_[row_idx + 1] - 1].second = value;
  }

  void get_diagonals(VectorX &diag) {
    diag.resize(n_rows_);
    for (int i = 0; i < n_rows_; ++i) {
      diag(i) = get_diagonal(i);
    }
  }

  void set_diagonals(const VectorX &diag) {
    for (int i = 0; i < n_rows_; ++i) {
      set_diagonal(i, diag(i));
    }
  }

 private:
  std::vector<std::pair<int, Scalar> > elements_;
  std::vector<int> row_addr_;
  int n_rows_;
};

#endif /* SPARSEMATRIXSTORAGE_H_ */
