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

#ifndef GEOMETRYSOLVER_H
#define GEOMETRYSOLVER_H

#include "SolverCommon.h"
#include "Constraint.h"
#include "SPDSolver.h"
#include "LinearRegularization.h"
#include "SparseMatrixStorage.h"
#include "OMPHelper.h"
#include "AndersonAcceleration.h"
#include "Parameters.h"
#include <vector>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <limits>
#include <deque>
#include <iomanip>
#include <set>

template<unsigned int N>
class GeometrySolver {
 protected:
  typedef MatrixT<N, 1> VectorN;
  typedef MatrixT<N, Eigen::Dynamic> MatrixNX;
  typedef MatrixT<Eigen::Dynamic, N> MatrixXN;

 public:
  GeometrySolver()
      : current_points_var_(NULL),
        alt_points_var_(NULL),
        prev_points_var_(NULL),
        full_coords_ptr_(NULL),
        projections_ptr_(NULL),
        SPD_solver_(NULL),
        global_update_type_(Parameters::FULL_SOLVE),
        solver_initialized_(false),
        handle_solver_initialized_(false) {
  }

  ~GeometrySolver() {
    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i) {
      if (constraints_[i]) {
        delete constraints_[i];
      }
    }
  }

  bool initialize(int n_points, Parameters::GlobalUpdateType global_update_type,
                  SPDSolverType spd_solver_type = LDLT_SOLVER,
                  const std::vector<int> &handle_indices = std::vector<int>()) {
    Timer timer;
    Timer::EventID t_begin = timer.get_time();

    if (std::set<int>(handle_indices.begin(), handle_indices.end()).size()
        != handle_indices.size()) {
      std::cerr << "Error: duplicate indices in handle array" << std::endl;
      return false;
    }

    int n_handles = handle_indices.size();
    handle_indices_.resize(n_handles);
    if (n_handles > 0) {
      handle_indices_ = Eigen::Map<const Eigen::VectorXi>(
          handle_indices.data(), handle_indices.size());
    }

    int n_constraints = static_cast<int>(constraints_.size());
    assert(n_points != 0);
    assert(n_constraints != 0);
    int n_var_points = n_points - n_handles;
    std::cout << "n_var_points =  " << n_var_points << std::endl;

    points_var_1_.setZero(N, n_var_points);
    points_var_2_.setZero(N, n_var_points);
    points_var_3_.setZero(N, n_var_points);

    if (n_handles > 0) {
      full_coords_1_.setZero(N, n_points);
    }

    // Set up full constraint matrix
    std::vector<Triplet> triplets;
    int idO = 0;
    for (int i = 0; i < n_constraints; ++i) {
      constraints_[i]->add_constraint(triplets, idO);
    }
    int projection_dim = idO;
    std::cout << "projection_dim =  " << projection_dim << std::endl;

    // Set up full global update matrix
    ColMajorSparseMatrix A = ColMajorSparseMatrix(idO, n_points);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    ColMajorSparseMatrix At = A.transpose();
    ColMajorSparseMatrix global_mat = At * A;

    // Set up regularization terms
    rhs_fixed_.setZero(N, n_var_points);
    MatrixXN full_rhs_fixed(n_points, N);
    full_rhs_fixed.setZero();
    ColMajorSparseMatrix L;
    if (regularization_.get_regularization_system(n_points, L,
                                                  regularization_rhs_)) {
      ColMajorSparseMatrix Lt = L.transpose();
      global_mat += Lt * L;
      full_rhs_fixed = Lt * regularization_rhs_;
      regularization_matrix_storage_.init_from_transposed_matrix(Lt);
    } else {
      regularization_matrix_storage_.clear();
    }

    projections_1_.setZero(N, projection_dim);
    constraint_err_.setZero(constraints_.size());
    reg_err_.setZero(L.rows());
    lbfgs_err_.setZero(n_var_points);

    // Reduce the full system to a system about variable points
    ColMajorSparseMatrix At_reduced;
    if (n_handles == 0) {
      At_reduced = At;
      rhs_fixed_ = full_rhs_fixed.transpose();

      var_point_indices_.resize(n_points);
      for (int i = 0; i < n_points; ++i) {
        var_point_indices_(i) = i;
      }
    } else {
      std::vector<bool> is_handle(n_points, false);
      for (int i = 0; i < n_handles; ++i) {
        is_handle[handle_indices[i]] = true;
      }

      // Set up selection matrix for handles and variable points
      std::vector<int> var_point_idx;
      std::vector<Triplet> handle_triplets, var_point_triplets;
      int handle_row = 0, var_row = 0;
      for (int i = 0; i < n_points; ++i) {
        if (is_handle[i]) {
          handle_triplets.push_back(Triplet(handle_row++, i, Scalar(1)));
        } else {
          var_point_idx.push_back(i);
          var_point_triplets.push_back(Triplet(var_row++, i, Scalar(1)));
        }
      }

      var_point_indices_ = Eigen::Map < Eigen::VectorXi
          > (var_point_idx.data(), var_point_idx.size());
      assert(n_var_points == var_row);
      assert(n_handles == handle_row);
      variables_selection_.resize(n_var_points, n_points);
      variables_selection_.setFromTriplets(var_point_triplets.begin(),
                                           var_point_triplets.end());
      variables_selection_.makeCompressed();
      handles_selection_.resize(n_handles, n_points);
      handles_selection_.setFromTriplets(handle_triplets.begin(),
                                         handle_triplets.end());
      handles_selection_.makeCompressed();

      ColMajorSparseMatrix var_to_full = variables_selection_.transpose();
      ColMajorSparseMatrix handle_to_full = handles_selection_.transpose();

      rhs_handle_contribution_ = -variables_selection_ * global_mat
          * handle_to_full;
      global_mat = variables_selection_ * (global_mat * var_to_full);
      rhs_fixed_ = (variables_selection_ * full_rhs_fixed).transpose();
      At_reduced = variables_selection_ * At;
    }

    if (global_update_type == Parameters::FULL_SOLVE) {
      At_ = At_reduced;
      if (spd_solver_type == LDLT_SOLVER) {
        SPD_solver_ = std::make_shared<SimplicialLDLTSolver>();
      } else if (spd_solver_type == LLT_SOLVER) {
        SPD_solver_ = std::make_shared<SimplicialLLTSolver>();
      } else {
        SPD_solver_ = std::make_shared<CGSolver>();
      }

      SPD_solver_->initialize(global_mat);
      if (SPD_solver_->info() != Eigen::Success) {
        std::cerr << "Error: SPD solver initialization failed" << std::endl;
        return false;
      }
    } else {
      At_storage_.init_from_transposed_matrix(At_reduced.transpose());
      if (!global_system_Jacobi_storage_.init_from_transposed_matrix(
          global_mat)) {
        std::cerr << "Error: Jacobi system initialization failed" << std::endl;
        return false;
      }
    }

    handle_solver_initialized_ = n_handles > 0;
    solver_initialized_ = n_handles == 0;
    global_update_type_ = global_update_type;

    std::cout << "predecomposition time = "
              << timer.elapsed_time(t_begin, timer.get_time()) << std::endl;

    return true;
  }

  void solve(const MatrixNX &init_points, int iter, int Anderson_m) {
    if (!(handle_solver_initialized_ || solver_initialized_)) {
      std::cerr << "Error: solver not initialized yet" << std::endl;
    }
    clear_iteration_history();

    bool has_handles = static_cast<int>(handle_indices_.size()) > 0;
    MatrixNX rhs_fixed_part = rhs_fixed_;
    if (has_handles) {
      rhs_fixed_part += (rhs_handle_contribution_
          * (handles_selection_ * init_points.transpose())).transpose();
      points_var_1_ =
          (variables_selection_ * init_points.transpose()).transpose();
      points_var_2_ = points_var_1_;
      points_var_3_ = points_var_1_;
    } else {
      points_var_1_ = init_points;
      points_var_2_ = init_points;
      points_var_3_ = init_points;
    }

    bool accelerate = Anderson_m > 0;
    //bool alt_storage_required = false;

    current_points_var_ = &points_var_1_;
    alt_points_var_ = &points_var_2_;
    prev_points_var_ = &points_var_3_;
    projections_ptr_ = &projections_1_;
    if (has_handles) {
      full_coords_1_ = init_points;
      full_coords_ptr_ = &full_coords_1_;
    } else {
      full_coords_ptr_ = current_points_var_;
    }

    int n_jacobi_rows =
        (global_update_type_ == Parameters::JACOBI_STEP) ?
            At_storage_.n_rows() : 0;
    int n_direct_columns =
        (global_update_type_ == Parameters::FULL_SOLVE) ? N : 0;
    int n_points_var = var_point_indices_.size();

    Scalar new_err = 0;
    Scalar current_err = std::numeric_limits<Scalar>::max();

    AndersonAcceleration accelerator;
    if (accelerate) {
      accelerator.init(Anderson_m, current_points_var_->size(),
                       current_points_var_->data());
    }

    Timer timer;
    Timer::EventID t_begin = timer.get_time();

    for (int it = 0; it <= iter; ++it) {
      int n_constraints = constraints_.size();
      int n_reg_terms = regularization_matrix_storage_.n_rows();
      int n_writeback_var_coords = has_handles ? n_points_var : 0;

      int n_post_constraints = 0, n_post_reg_terms = 0,
          n_post_writeback_var_coords = 0;
      bool new_err_needs_update = false;

      OMP_PARALLEL
      {
        // Write the variable point coordinates into the full array
        OMP_FOR
        for (int i = 0; i < n_writeback_var_coords; ++i) {
          full_coords_ptr_->col(var_point_indices_(i)) = current_points_var_
              ->col(i);
        }

        OMP_FOR
        for (int i = 0; i < n_constraints; ++i) {
          constraint_err_(i) = constraints_[i]->project(*full_coords_ptr_,
                                                        *projections_ptr_);
        }

        OMP_FOR
        for (int i = 0; i < n_reg_terms; ++i) {
          reg_err_(i) = evaluate_regularization_term(i, *full_coords_ptr_);
        }

        OMP_SECTIONS
        {
          OMP_SECTION
          {
            new_err = constraint_err_.sum() + reg_err_.sum();
          }
        }

        OMP_SINGLE
        {
          bool require_swap = (it > 0) && accelerate && (new_err > current_err);

          if (require_swap) {
            std::swap(current_points_var_, alt_points_var_);
            accelerator.replace(current_points_var_->data());

            if (!has_handles) {
              full_coords_ptr_ = current_points_var_;
            }

            // Request re-computation of projections and error values
            n_post_constraints = n_constraints;
            n_post_reg_terms = n_reg_terms;
            n_post_writeback_var_coords = n_writeback_var_coords;
            new_err_needs_update = true;

          } else {
            current_err = new_err;
          }

          Anderson_reset_.push_back(require_swap);

        }

        // Re-projection in case the new error is not smaller than the previous one
        OMP_FOR
        for (int i = 0; i < n_post_writeback_var_coords; ++i) {
          full_coords_ptr_->col(var_point_indices_(i)) = current_points_var_
              ->col(i);
        }

        OMP_FOR
        for (int i = 0; i < n_post_constraints; ++i) {
          constraint_err_(i) = constraints_[i]->project(*full_coords_ptr_,
                                                        *projections_ptr_);
        }

        OMP_FOR
        for (int i = 0; i < n_post_reg_terms; ++i) {
          reg_err_(i) = evaluate_regularization_term(i, *full_coords_ptr_);
        }

        OMP_SINGLE
        {
          if (new_err_needs_update) {
            current_err = constraint_err_.sum() + reg_err_.sum();
          }

          // If we reach the last iteration, set n_rows to zero to skip the following for loops of update
          if (it == iter) {
            n_jacobi_rows = n_direct_columns = 0;
          } else {
            // Only swap the pointer if we have more iterations,
            // This ensure current_points_var_ points to the final solution
            std::swap(prev_points_var_, current_points_var_);
          }

          Timer::EventID t_iter = timer.get_time();
          elapsed_time_.push_back(timer.elapsed_time(t_begin, t_iter));
          function_values_.push_back(current_err);
        }

        // Compute new points using Jacobi iteration
        OMP_FOR
        for (int i = 0; i < n_jacobi_rows; ++i) {
          // Construct right hand sides
          VectorN rhs = rhs_fixed_part.col(i);
          At_storage_.evaluate(i, *projections_ptr_, rhs);

          // Perform Jacobi iteration
          global_system_Jacobi_storage_.evaluate(i, *prev_points_var_, rhs);

          // use the alternative point array to store the result
          alt_points_var_->col(i) = rhs;
        }

        // Compute new points using direct solve
        OMP_FOR
        for (int i = 0; i < n_direct_columns; ++i) {
          alt_points_var_->row(i) = (SPD_solver_->solve(
              rhs_fixed_part.row(i).transpose()
                  + At_ * (projections_ptr_->row(i).transpose()),
              prev_points_var_->row(i).transpose())).transpose();
        }
      }

      // Anderson acceleration
      if (it != iter) {
        if (accelerate) {
          (*current_points_var_) = Eigen::Map<const MatrixNX>(
              accelerator.compute(alt_points_var_->data()).data(), N,
              current_points_var_->cols());
        } else {
          // If not using acceleration, copy alt_points to current_points by swapping pointers
          std::swap(current_points_var_, alt_points_var_);
        }
      }

      if (!has_handles) {
        full_coords_ptr_ = current_points_var_;
      }
    }
  }

  const MatrixNX& get_points() {
    return *full_coords_ptr_;
  }

  void add_constraint(Constraint<N>* constraint) {
    constraints_.push_back(constraint);
  }

  void add_closeness(int idx, Scalar weight, const VectorN &target_pt) {
    regularization_.add_closeness(idx, weight, target_pt);
  }

  void add_uniform_laplacian(const std::vector<int> &indices, Scalar weight) {
    regularization_.add_uniform_laplacian(indices, weight);
  }

  void add_laplacian(const std::vector<int> &indices,
                     const std::vector<Scalar> coefs, Scalar weight) {
    regularization_.add_laplacian(indices, coefs, weight);
  }

  void add_relative_uniform_laplacian(const std::vector<int> &indices,
                                      Scalar weight,
                                      const MatrixNX &ref_points) {
    regularization_.add_relative_uniform_laplacian(indices, weight, ref_points);
  }

  void add_relative_laplacian(const std::vector<int> &indices,
                              const std::vector<Scalar> coefs, Scalar weight,
                              const MatrixNX &ref_points) {
    regularization_.add_relative_laplacian(indices, coefs, weight, ref_points);
  }

  void output_iteration_history(SolverType solver_type) {
    if (solver_type == AA_SOLVER)
      assert(function_values_.size() == Anderson_reset_.size());

    assert(function_values_.size() == elapsed_time_.size());
    int n_iter = function_values_.size();
    for (int i = 0; i < n_iter; ++i) {
      std::cout << "Iteration " << i << ": ";
      std::cout << std::setprecision(6) << elapsed_time_[i] << " secs, ";
      std::cout << " target value " << std::setprecision(16)
                << function_values_[i];
      if ((solver_type == AA_SOLVER) && (Anderson_reset_[i])) {
        std::cout << " (reject accelerator)";
      }

      std::cout << std::endl;
    }

    std::cout << std::endl;
  }

 protected:
  MatrixNX points_var_1_, points_var_2_, points_var_3_;	// Vertex coordinate matrices

  // Pointers to the coordinate matrices
  // - prev_points_var: points in previous iteration
  // - alt_pionts_var: alternative point coordinates, for storing unaccelerated point positions
  // - current_points_var: final point coordinates for the current iteration
  MatrixNX *current_points_var_, *alt_points_var_, *prev_points_var_;

  MatrixNX full_coords_1_, full_coords_2_;	// Coordinates of all points
  MatrixNX *full_coords_ptr_;

  // Constraints and regularization terms
  std::vector<Constraint<N>*> constraints_;
  MatrixNX projections_1_, projections_2_;
  MatrixNX *projections_ptr_;
  LinearRegularization<N> regularization_;
  RowSparseStorage regularization_matrix_storage_;
  MatrixXN regularization_rhs_;
  VectorX constraint_err_, reg_err_, alt_constraint_err_, alt_reg_err_,
      lbfgs_err_;	// Storage for error term values from constraints and regularization

  // Matrices mapping the full points matrix to the variable points and the handle points
  ColMajorSparseMatrix variables_selection_, handles_selection_;
  Eigen::VectorXi handle_indices_, var_point_indices_;	// The indices of handles and variable points within the full point array
  ColMajorSparseMatrix rhs_handle_contribution_;	// Matrix that encodes the handle positions contribution to the right hand side of global update system

  // Data structures used for direct solve update
  std::shared_ptr<SPDSolver> SPD_solver_;
  ColMajorSparseMatrix At_;
  MatrixNX rhs_fixed_;

  // Data structures used for Jacobi update
  RowSparseStorage At_storage_;
  JacobiRowSparseStorage global_system_Jacobi_storage_;

  Parameters::GlobalUpdateType global_update_type_;

 public:
  // History of iterations
  std::deque<bool> Anderson_reset_;
  std::vector<double> function_values_;
  std::vector<double> elapsed_time_;

 private:
  bool solver_initialized_;
  bool handle_solver_initialized_;

  void clear_iteration_history() {
    Anderson_reset_.clear();
    function_values_.clear();
    elapsed_time_.clear();
  }

  Scalar evaluate_regularization_term(int term_idx,
                                      const MatrixNX &full_coords) {
    VectorN reg_term = -regularization_rhs_.row(term_idx).transpose();
    regularization_matrix_storage_.evaluate(term_idx, full_coords, reg_term);
    return reg_term.squaredNorm() * 0.5;
  }
};

#endif	// GEOMETRYSOLVER_H
