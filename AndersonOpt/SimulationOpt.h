//  BSD 3-Clause License
//
//  Copyright (c) 2018, Yue Peng, Bailin Deng
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


#ifndef SIMULATIONOPT_H_
#define SIMULATIONOPT_H_

#include "MeshTypes.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <Glut/glut.h>
#elif defined(_WIN32) || defined(_WIN64)
#include <GLUT/glcmolorut.h>
#else
#include <GL/glut.h>
#endif
#include "DynamicsSolver.h"
#include "Parameters.h"

class SimulationOpt {
 public:
  SimulationOpt()
      : time_step_(0.033),
        n_frames_(2),
        n_nodes_(0),
        countFrame(0),
        elasticity(sqrt(5000000.0)),
        is_solve_begin(false) {
  }

  ~SimulationOpt() {
  }
  void simulate(int solver_iter) {
    countFrame++;

    solver.solve(solver_iter);

    if (countFrame > 1)
      node_positions_.col(0) = node_positions_.col(1);
    node_positions_.col(1) = Eigen::Map<const Eigen::VectorXd>(
        solver.get_points().data(), 3 * n_nodes_);

    is_solve_begin = true;
  }

  void setTetMessage(const Eigen::Matrix4Xi &tet_connectivity_,
                     const Eigen::Matrix3Xi &tet_boundary_,
                     const std::vector<int> &handle_indices) {
    tet_boundary = tet_boundary_;
    tet_connectivity = tet_connectivity_;

    fixed_node_idx_.setZero(handle_indices.size());
    for (size_t i = 0; i < handle_indices.size(); i++) {
      fixed_node_idx_(i) = handle_indices[i];
    }
  }

  void initialize_tet_data(const Matrix3X &problem_def_poisitions) {
    double g = 9.8;

    n_nodes_ = problem_def_poisitions.cols();
    node_positions_.resize(3 * n_nodes_, n_frames_);
    node_positions_.setZero();
    external_forces_.resize(3, n_nodes_);

    std::cout << "Initializing masses and forces..." << std::endl;

    double density = 10.0;
    mass_.resize(n_nodes_);
    mass_.setZero();
    for (int id = 0; id < tet_connectivity.cols(); ++id) {
      Matrix33 edge_vecs;
      for (int i = 0; i < 3; ++i) {
        edge_vecs.col(i) = problem_def_poisitions.col(
            tet_connectivity(i + 1, id))
            - problem_def_poisitions.col(tet_connectivity(0, id));
      }

      double Vi = std::fabs(edge_vecs.determinant()) / 6.0;
      for (int i = 0; i < 4; ++i) {
        mass_(tet_connectivity(i, id)) += Vi * density;
      }

    }

    const Eigen::Vector3d G(0.0, -g, 0.0);
    for (int j = 0; j < n_nodes_; ++j) {
      external_forces_.col(j) = G * mass_(j);
    }
    init_node_positions_ = problem_def_poisitions;
  }

  void initSolver(const Matrix3X &problem_def_poisitions, bool use_fast_svd,
                  Parameters param) {
    initialize_tet_data(problem_def_poisitions);

    std::cout << "Setting up solver..." << std::endl;

    //set elasticity
    elasticity = sqrt(param.elasticity);
    time_step_ = param.time_step;

    node_positions_.col(0) = Eigen::Map<const Eigen::VectorXd>(
        problem_def_poisitions.data(), 3 * n_nodes_);
    Matrix3X init_velocities(3, problem_def_poisitions.cols());
    init_velocities.setZero();

    // corotated strain for each tetrahedren
    int n_tets = tet_connectivity.cols();
    std::cout << "n_tets = " << n_tets << std::endl;
    for (int i = 0; i < n_tets; ++i) {
      std::vector<int> vtx_idx;
      for (int j = 0; j < 4; ++j) {
        vtx_idx.push_back(tet_connectivity(j, i));
      }

      solver.add_constraint(
          new TetStrainConstraint(vtx_idx, elasticity, problem_def_poisitions,
                                  use_fast_svd));
    }

    double damping = 0.99;

    std::vector<int> fixed_points;
    for (int i = 0; i < static_cast<int>(fixed_node_idx_.size()); ++i) {
      fixed_points.push_back(fixed_node_idx_(i));
    }

    if (!param.acceleration)
      param.anderson_m = 0;

    if (!solver.initialize(problem_def_poisitions, init_velocities,
                           param.anderson_m, mass_, time_step_, fixed_points,
                           damping, elasticity, false)) {
      std::cerr << "Initialization error" << std::endl;
      return;
    }
  }

  void Render(int visual_mode = 0, int render_mode = 0) {
    //Build_VN();

    Matrix3X m(3, init_node_positions_.cols());
    if (is_solve_begin) {
      m = Eigen::Map < Matrix3X > (node_positions_.col(1).data(), 3, n_nodes_);
    } else {
      m = Eigen::Map < Matrix3X > (node_positions_.col(0).data(), 3, n_nodes_);
    }

    if (visual_mode == 0) {
      glEnable(GL_LIGHTING);
      float diffuse_color[3] = { 0.52, 0.52, 0.52 };
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
      for (int i = 0; i < tet_boundary.cols(); i++) {
        Vector3 p0 = m.col(tet_boundary(0, i));
        Vector3 p1 = m.col(tet_boundary(1, i));
        Vector3 p2 = m.col(tet_boundary(2, i));
        Vector3 e1 = p1 - p0;
        Vector3 e2 = p2 - p0;
        Vector3 n = e1.cross(e2);
        n.normalize();

        if (render_mode == 0)
          glNormal3f(n(0), n(1), n(2));
        glBegin(GL_TRIANGLES);
        glVertex3d(p0(0), p0(1), p0(2));
        glVertex3d(p1(0), p1(1), p1(2));
        glVertex3d(p2(0), p2(1), p2(2));
        glEnd();
      }
    }
  }

 private:
  double time_step_;
  int n_frames_;
  int n_nodes_;
  int countFrame;
  double elasticity;

  bool is_solve_begin;

  Matrix3X init_node_positions_;
  Eigen::MatrixXd node_positions_;
  Eigen::VectorXd mass_;
  Eigen::Matrix3Xd external_forces_;
  Eigen::VectorXi fixed_node_idx_;

  Eigen::Matrix4Xi tet_connectivity;
  Eigen::Matrix3Xi tet_boundary;

  DynamicsSolver solver;

};

#endif /* SIMULATIONOPT_H_ */
