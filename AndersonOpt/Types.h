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

#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef USE_FLOAT_SCALAR
typedef float Scalar
#else
typedef double Scalar;
#endif

typedef Eigen::SparseMatrix<Scalar, Eigen::ColMajor> ColMajorSparseMatrix;
typedef Eigen::SparseMatrix<Scalar, Eigen::RowMajor> RowMajorSparseMatrix;
typedef Eigen::Triplet<Scalar> Triplet;

#ifdef EIGEN_DONT_ALIGN
#define EIGEN_ALIGNMENT Eigen::DontAlign
#else
#define EIGEN_ALIGNMENT Eigen::AutoAlign
#endif

template<int Rows, int Cols, int Options = (Eigen::ColMajor | EIGEN_ALIGNMENT) >
using MatrixT = Eigen::Matrix<Scalar, Rows, Cols, Options>;  ///< A typedef of the dense matrix of Eigen.
typedef MatrixT<2, 1> Vector2;								///< A 2d column vector.
typedef MatrixT<2, 2> Matrix22;								///< A 2 by 2 matrix.
typedef MatrixT<2, 3> Matrix23;								///< A 2 by 3 matrix.
typedef MatrixT<3, 1> Vector3;								///< A 3d column vector.
typedef MatrixT<3, 2> Matrix32;								///< A 3 by 2 matrix.
typedef MatrixT<3, 3> Matrix33;								///< A 3 by 3 matrix.
typedef MatrixT<3, 4> Matrix34;								///< A 3 by 4 matrix.
typedef MatrixT<4, 1> Vector4;								///< A 4d column vector.
typedef MatrixT<4, 4> Matrix44;								///< A 4 by 4 matrix.
typedef MatrixT<4, Eigen::Dynamic> Matrix4X;				///< A 4 by n matrix.
typedef MatrixT<3, Eigen::Dynamic> Matrix3X;				///< A 3 by n matrix.
typedef MatrixT<Eigen::Dynamic, 3> MatrixX3;				///< A n by 3 matrix.
typedef MatrixT<2, Eigen::Dynamic> Matrix2X;				///< A 2 by n matrix.
typedef MatrixT<Eigen::Dynamic, 2> MatrixX2;				///< A n by 2 matrix.
typedef MatrixT<Eigen::Dynamic, 1> VectorX;					///< A nd column vector.
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> MatrixXX;  ///< A n by m matrix.

// Conversion between a 3d vector type to Eigen::Vector3d
template<typename Vec_T>
inline Vector3 to_eigen_vec3(const Vec_T &vec) {
  return Vector3(vec[0], vec[1], vec[2]);
}

template<typename Vec_T>
inline Vec_T from_eigen_vec3(const Vector3 &vec) {
  Vec_T v;
  v[0] = vec(0);
  v[1] = vec(1);
  v[2] = vec(2);

  return v;
}

#endif // TYPES_H
///////////////////////////////////////////////////////////////////////////////
