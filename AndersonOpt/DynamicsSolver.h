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

#ifndef DYNAMICSSOLVER_H_
#define DYNAMICSSOLVER_H_

#include "SolverCommon.h"
#include "Types.h"
#include "Constraint.h"
#include "SparseMatrixStorage.h"
#include "OMPHelper.h"
#include "AndersonAcceleration.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include <cassert>

using namespace std;

#define LARGER_EPSILON 1e-6

class DynamicsSolver {
 public:

  DynamicsSolver()
      : current_points_(NULL),
        alt_points_(NULL),
        old_points_(NULL),
        prev_points_(NULL),
        n_points_(0),
        time_step_(0),
        air_damping_(1),
        pos_constraint_weight_(1),
        handle_idx_(-1),
        handle_active_(false),
        SPD_solver_ptr_(NULL),
        accelerate_(false),
        Anderson_m_(0),
        initialized_(false),
        initialized_with_handle_(false),
        elasticity_(0) {
  }

  ~DynamicsSolver() {
    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i) {
      if (constraints_[i]) {
        delete constraints_[i];
      }
    }
  }

  void add_constraint(Constraint<3>* constraint) {
    constraints_.push_back(constraint);
  }

  Scalar get_fixed_point_weight() {
    return pos_constraint_weight_;
  }

  bool initialize(const Matrix3X &init_pos, const Matrix3X &init_velocities,
                  int Anderson_m, const VectorX &mass, Scalar time_step,
                  const std::vector<int> &fixed_points, Scalar air_damping = 1,
                  Scalar elasticity = 0, bool with_handle = false,
                  Scalar fixed_points_relative_weight = 5) {
    assert(init_pos.cols() == init_velocities.cols());
    assert(init_pos.cols() == mass.size());

    n_points_ = init_pos.cols();
    assert(n_points_ > 0);
    points_1_ = init_pos;
    points_2_ = init_pos;
    points_3_ = init_pos;
    points_4_ = init_pos;

    current_points_ = &points_1_;
    alt_points_ = &points_2_;
    old_points_ = &points_3_;
    prev_points_ = &points_4_;

    mass_ = mass;
    velocities_ = init_velocities;
    time_step_ = time_step;
    elasticity_ = elasticity;

    momentum_.setZero(3, n_points_);
    momentum_and_pos_rhs_ = momentum_;
    air_damping_ = air_damping;
    unit_gravity_ = Vector3(Scalar(0), Scalar(-9.8), Scalar(0));	// Define gravity force along -y direction, to be compatible with OpenGL
    gravity_contribution_to_momentum_.resize(3, n_points_);
    gravity_contribution_to_momentum_.colwise() = unit_gravity_
        * (time_step_ * time_step_);
    fixed_point_idx_ = fixed_points;
    fixed_point_pos_.clear();
    for (int i = 0; i < static_cast<int>(fixed_point_idx_.size()); ++i) {
      int vtx_idx = fixed_point_idx_[i];
      assert(vtx_idx >= 0 && vtx_idx < n_points_);
      fixed_point_pos_.push_back(init_pos.col(vtx_idx));
    }
    handle_idx_ = -1;

    // Set up full constraint matrix
    int n_constraints = static_cast<int>(constraints_.size());
    assert(n_constraints > 0);
    std::vector<Triplet> constraint_mat_triplets;
    int idO = 0;
    for (int i = 0; i < n_constraints; ++i) {
      constraints_[i]->add_constraint(constraint_mat_triplets, idO);
    }
    projections_.setZero(3, idO);
    constraint_err_.setZero(n_constraints);
    momentum_and_pos_err_.setZero(n_points_);

    // Construct constraint matrix
    ColMajorSparseMatrix A = ColMajorSparseMatrix(idO, n_points_);
    A.setFromTriplets(constraint_mat_triplets.begin(),
                      constraint_mat_triplets.end());
    A.makeCompressed();
    ColMajorSparseMatrix At = A.transpose();
    ColMajorSparseMatrix global_mat = (At * A) * (time_step_ * time_step_);
    VectorX global_mat_diag = global_mat.diagonal();
    global_mat_diag += mass_;

    pos_constraint_weight_ = std::max(
        Scalar(1), fixed_points_relative_weight * global_mat_diag.maxCoeff());

    // Construct the diagonal matrix for mass and fixed point weights
    ColMajorSparseMatrix M(n_points_, n_points_);
    std::vector<Triplet> M_triplets;
    for (int i = 0; i < n_points_; ++i) {
      M_triplets.push_back(Triplet(i, i, mass_(i)));
    }
    for (int i = 0; i < static_cast<int>(fixed_point_idx_.size()); ++i) {
      int vtx_idx = fixed_point_idx_[i];
      M_triplets.push_back(Triplet(vtx_idx, vtx_idx, pos_constraint_weight_));
    }
    M.setFromTriplets(M_triplets.begin(), M_triplets.end());
    global_mat += M;
    if (with_handle) {
      system_mat_diag_values_ = global_mat.diagonal();
    }

    At_ = At * (time_step_ * time_step_);
    SPD_solver_.compute(global_mat);
    if (SPD_solver_.info() != Eigen::Success) {
      std::cerr << "Error: System SPD solver initialization failed"
                << std::endl;
      return false;
    }

    if (with_handle) {
      handle_system_mat_ = global_mat;
      handle_system_mat_.makeCompressed();

      handle_system_mat_diag_entry_addr_.assign(n_points_, NULL);
      for (int k = 0; k < handle_system_mat_.outerSize(); ++k) {
        for (ColMajorSparseMatrix::InnerIterator it(handle_system_mat_, k); it;
            ++it) {
          if (it.row() == it.col()) {
            handle_system_mat_diag_entry_addr_[it.row()] = &(it.valueRef());
          }
        }
      }

      SPD_solver_with_handles_.analyzePattern(handle_system_mat_);
      if (SPD_solver_with_handles_.info() != Eigen::Success) {
        std::cerr << "Error: Handle system SPD solver initialization failed"
                  << std::endl;
        return false;
      }
    }

    Anderson_m_ = Anderson_m;
    accelerate_ = Anderson_m_ > 0;

    initialized_ = true;
    initialized_with_handle_ = with_handle;

    return true;
  }

  bool set_handle(int handle_vtx_idx) {
    if (!initialized_with_handle_) {
      std::cerr << "Error: solver not initialized for interactive simulation"
                << std::endl;
      return false;
    }

    if (handle_vtx_idx < 0 || handle_vtx_idx >= n_points_) {
      std::cerr << "Error: invalid handle index" << std::endl;
      return false;
    }

    if (handle_vtx_idx == handle_idx_) {
      // No need to update
      handle_active_ = true;
      return true;
    }

    // Restore the diagonal entry for the previous handle vertex
    if (handle_idx_ >= 0
        && std::find(fixed_point_idx_.begin(), fixed_point_idx_.end(),
                     handle_idx_) == fixed_point_idx_.end()) {
      *(handle_system_mat_diag_entry_addr_[handle_idx_]) =
          system_mat_diag_values_[handle_idx_];
    }

    handle_idx_ = handle_vtx_idx;
    Scalar new_diag_value =
        (std::find(fixed_point_idx_.begin(), fixed_point_idx_.end(),
                   handle_idx_) == fixed_point_idx_.end()) ?
            (system_mat_diag_values_[handle_idx_] + pos_constraint_weight_) :
            system_mat_diag_values_[handle_idx_];

    // Update system matrix diagonal entry for the handle
    *(handle_system_mat_diag_entry_addr_[handle_idx_]) = new_diag_value;
    SPD_solver_with_handles_.factorize(handle_system_mat_);

    handle_active_ = true;

    return true;
  }

  void disable_handle() {
    handle_active_ = false;
  }

  void solve(int iter, const Vector3 *handle_pos = NULL) {

    if (handle_pos) {
      if (!initialized_with_handle_) {
        std::cerr << "Error: solver not initialized yet" << std::endl;
        return;
      }

      if (handle_idx_ < 0) {
        std::cerr << "Error: handle position required" << std::endl;
        return;
      }

      if (!handle_active_) {
        std::cerr << "Error: handle is not active" << std::endl;
      }

      // If we are using initial fixed points as handles, update their target positions
      std::vector<int>::iterator i = std::find(fixed_point_idx_.begin(),
                                               fixed_point_idx_.end(),
                                               handle_idx_);
      if (i != fixed_point_idx_.end()) {
        fixed_point_pos_[i - fixed_point_idx_.begin()] = *handle_pos;
      }
      pos_constrained_vtx_idx_ = fixed_point_idx_;
      constrianed_pos_ = fixed_point_pos_;

      // If the current handle is not an initial fixed point, add it to the position-constrained vertex list
      if (i == fixed_point_idx_.end()) {
        pos_constrained_vtx_idx_.push_back(handle_idx_);
        constrianed_pos_.push_back(*handle_pos);
      }

      SPD_solver_ptr_ = &SPD_solver_with_handles_;
    } else {
      if (!initialized_) {
        std::cerr << "Error: solver not initialized yet" << std::endl;
        return;
      }

      pos_constrained_vtx_idx_ = fixed_point_idx_;
      constrianed_pos_ = fixed_point_pos_;

      SPD_solver_ptr_ = &SPD_solver_;
    }

    // Compute momentum term
    int n_pos_constraints = pos_constrained_vtx_idx_.size();
    momentum_ = (*current_points_) + velocities_ * time_step_
        + gravity_contribution_to_momentum_;
    std::swap(current_points_, old_points_);
    (*current_points_) = momentum_;
    momentum_and_pos_rhs_ = momentum_ * mass_.asDiagonal();
    for (int i = 0; i < n_pos_constraints; ++i) {
      momentum_and_pos_rhs_.col(pos_constrained_vtx_idx_[i]) +=
          pos_constraint_weight_ * constrianed_pos_[i];
    }

    if (accelerate_) {
      accelerator_.init(Anderson_m_, old_points_->size(), old_points_->data());
    }

    int n_direct_columns = 3;

    Scalar new_err = 0;
    Scalar current_err = std::numeric_limits<Scalar>::max();

    for (int it = 0; it <= iter; ++it) {
      int n_constraints = constraints_.size();
      int n_momentum_terms = n_points_;
      int n_pos_terms = pos_constrained_vtx_idx_.size();
      ;
      bool new_err_needs_update = false;

      OMP_PARALLEL
      {
        OMP_FOR
        for (int i = 0; i < n_constraints; ++i) {
          constraint_err_(i) = constraints_[i]->project(*current_points_,
                                                        projections_);
        }

        OMP_FOR
        for (int i = 0; i < n_momentum_terms; ++i) {
          momentum_and_pos_err_(i) = evaluate_momentum_term(i,
                                                            *current_points_);
        }

        OMP_FOR
        for (int i = 0; i < n_pos_terms; ++i) {
          momentum_and_pos_err_(pos_constrained_vtx_idx_[i]) +=
              evaluate_position_term(i, *current_points_);
        }

        OMP_SINGLE
        {
          new_err = constraint_err_.sum() * (time_step_ * time_step_)
              + momentum_and_pos_err_.sum();
          bool require_swap = (it > 0) && accelerate_
              && (new_err > current_err);

          if (require_swap) {
            std::swap(current_points_, alt_points_);
            accelerator_.replace(current_points_->data());
            new_err_needs_update = true;
          } else {
            // Disable further function check for alternative solution
            n_constraints = 0;
            n_momentum_terms = 0;
            n_pos_terms = 0;
            current_err = new_err;
          }

        }

        // Re-projection in case the new error is not smaller than the previous one
        OMP_FOR
        for (int i = 0; i < n_constraints; ++i) {
          constraint_err_(i) = constraints_[i]->project(*current_points_,
                                                        projections_);
        }

        OMP_FOR
        for (int i = 0; i < n_momentum_terms; ++i) {
          momentum_and_pos_err_(i) = evaluate_momentum_term(i,
                                                            *current_points_);
        }

        OMP_FOR
        for (int i = 0; i < n_pos_terms; ++i) {
          momentum_and_pos_err_(pos_constrained_vtx_idx_[i]) +=
              evaluate_position_term(i, *current_points_);
        }

        OMP_SINGLE
        {
          if (new_err_needs_update) {
            current_err = constraint_err_.sum() * (time_step_ * time_step_)
                + momentum_and_pos_err_.sum();
          }

          // If we reach the last iteration, skip the following for-loops of update
          if (it == iter)
            n_direct_columns = 0;
        }

        // Compute new points using direct solve
        OMP_FOR
        for (int i = 0; i < n_direct_columns; ++i) {
          alt_points_->row(i) = (SPD_solver_ptr_->solve(
              momentum_and_pos_rhs_.row(i).transpose()
                  + At_ * (projections_.row(i).transpose()))).transpose();
        }
      }

      // Anderson acceleration
      if (it != iter) {
        if (accelerate_) {
          (*current_points_) = Eigen::Map<const Matrix3X>(
              accelerator_.compute(alt_points_->data()).data(), 3,
              current_points_->cols());
        } else {
          // If not using acceleration, copy alt_points to current_points by swapping pointers
          std::swap(current_points_, alt_points_);
        }
      }
    }

    update_velocities();
  }

  const Matrix3X& get_points() const {
    return *current_points_;
  }

 protected:
  Matrix3X points_1_, points_2_, points_3_, points_4_;	// Point coordinate matrices

  // Pointers to the coordinate matrices
  // - prev_points_var: points in previous solver iteration
  // - alt_pionts_var: alternative point coordinates, for storing unaccelerated point positions
  // - current_points_var: final point coordinates for the current iteration
  // - old_points_: old positions before solving
  Matrix3X *current_points_, *alt_points_, *old_points_, *prev_points_;

  Matrix3X momentum_, momentum_and_pos_rhs_;

  int n_points_;
  VectorX mass_;
  Matrix3X velocities_;
  Scalar time_step_;
  Scalar air_damping_;
  Vector3 unit_gravity_;
  Matrix3X gravity_contribution_to_momentum_;

  // Constraints terms
  std::vector<Constraint<3>*> constraints_;
  Matrix3X projections_;
  VectorX constraint_err_;

  // Data for evaluating momentum and position error
  Scalar pos_constraint_weight_;
  std::vector<int> pos_constrained_vtx_idx_;
  std::vector<Vector3> constrianed_pos_;
  VectorX momentum_and_pos_err_;

  // Data structures used for direct solve update
  Eigen::SimplicialLDLT<ColMajorSparseMatrix, Eigen::Lower> SPD_solver_;
  ColMajorSparseMatrix system_mat_;  // Original system matrix without handles
  ColMajorSparseMatrix At_;

  // Fixed points: points that are fixed during the whole simulation process. Their target positions might change if used as handles
  std::vector<int> fixed_point_idx_;
  std::vector<Vector3> fixed_point_pos_;

  // Data structures for handling interactive handle
  int handle_idx_;
  bool handle_active_;
  Eigen::SimplicialLDLT<ColMajorSparseMatrix, Eigen::Lower> SPD_solver_with_handles_;
  VectorX system_mat_diag_values_;	// Diagonal entries for the original system matrix
  ColMajorSparseMatrix handle_system_mat_;	// System matrix with added diagonal entry for handle
  std::vector<Scalar*> handle_system_mat_diag_entry_addr_;	// Address of diagonal entries within the handle system matrix array

  Eigen::SimplicialLDLT<ColMajorSparseMatrix, Eigen::Lower>* SPD_solver_ptr_;

  bool accelerate_;
  AndersonAcceleration accelerator_;
  int Anderson_m_;

  Scalar elasticity_;

 private:
  bool initialized_;
  bool initialized_with_handle_;

  Scalar evaluate_momentum_term(int vtx_idx, const Matrix3X &coords) {
    return 0.5 * mass_(vtx_idx)
        * (coords.col(vtx_idx) - momentum_.col(vtx_idx)).squaredNorm();
  }

  Scalar evaluate_position_term(int term_idx, const Matrix3X &vtx_coords) {
    return 0.5 * pos_constraint_weight_
        * (vtx_coords.col(pos_constrained_vtx_idx_[term_idx])
            - constrianed_pos_[term_idx]).squaredNorm();
  }

  Vector3 evaluate_position_term_Gradient(int term_idx,
                                          const Matrix3X &vtx_coords) {
    return pos_constraint_weight_
        * (vtx_coords.col(pos_constrained_vtx_idx_[term_idx])
            - constrianed_pos_[term_idx]);
  }

  void update_velocities() {
    velocities_ = ((*current_points_) - (*old_points_))
        * (air_damping_ / time_step_);
  }

};

#endif
