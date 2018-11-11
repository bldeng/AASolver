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

#ifndef TETMESHIO_H_
#define TETMESHIO_H_

#include "Types.h"
#include "MeshTypes.h"
#include <utility>
#include <sstream>
#include <iostream>
#include <array>
#include <string>

class TetMeshIO {
 public:

  static bool load(const std::string base_name, Matrix3X &node_positions,
                   Eigen::Matrix4Xi &tet_vtx_indices,
                   Eigen::Matrix3Xi &boundary_faces) {
    std::string node_file_name = base_name + ".node";

    bool zero_based_index;
    if (!load_node_file(node_file_name, node_positions, zero_based_index)) {
      return false;
    }

    std::string ele_file_name = base_name + ".ele";
    if (!load_ele_file(ele_file_name, tet_vtx_indices)) {
      return false;
    }
    if (!zero_based_index) {
      tet_vtx_indices.array() -= 1;
    }

    build_boundary_triangles(tet_vtx_indices, boundary_faces);

    std::cout << "Successfully read tetrahedral mesh info: "
        << node_positions.cols() << " vertices, " << tet_vtx_indices.cols()
        << " tetrahedrons, " << boundary_faces.cols() << " boundary faces"
        << std::endl;

    return true;
  }

  static bool save_boundary_mesh(const std::string& filename,
                                 const Matrix3X &node_positions,
                                 const Eigen::Matrix3Xi &boundary_faces) {
    std::ofstream ofile(filename);
    if (!ofile.is_open()) {
      std::cerr << "Error: unable to open file " << filename << std::endl;
      return false;
    }

    for (int i = 0; i < static_cast<int>(node_positions.cols()); ++i) {
      ofile << "v " << node_positions(0, i) << " " << node_positions(1, i)
          << " " << node_positions(2, i) << std::endl;
    }

    for (int i = 0; i < static_cast<int>(boundary_faces.cols()); ++i) {
      int v1 = boundary_faces(0, i) + 1, v2 = boundary_faces(1, i) + 1, v3 =
          boundary_faces(2, i) + 1;
      ofile << "f " << v1 << " " << v2 << " " << v3 << std::endl;
    }

    return true;
  }

 private:

  static bool load_node_file(const std::string &filename,
                             Matrix3X &vertex_positions,
                             bool &zero_based_index) {
    std::ifstream ifile(filename);
    if (!ifile.is_open()) {
      std::cerr << "Error while opening file " << filename << std::endl;
      return false;
    }

    bool has_meta_info = false;
    Matrix3X all_nodes;

    zero_based_index = false;
    int n_nodes = -1;

    std::string line;
    while (std::getline(ifile, line)) {
      std::string::size_type pos = line.find_first_not_of(' ');
      if (pos == std::string::npos) {
        continue;
      }

      // Check for comment line
      else if (line.at(pos) == '#') {
        continue;
      }

      std::stringstream line_stream(line);

      if (!has_meta_info) {
        int dim = -1, attrib, boundary;
        if (!(line_stream >> n_nodes >> dim >> attrib >> boundary)) {
          std::cerr << "Unable to parse line: " << line << std::endl;
          return false;
        }

        assert(dim == 3);
        assert(n_nodes >= 0);
        all_nodes.resize(3, n_nodes + 1);  // Initialize with an additional column, to accommodate both 0- and 1-based indexing
        has_meta_info = true;
      } else {
        int node_idx = -1;
        Scalar x, y, z;
        if (!(line_stream >> node_idx >> x >> y >> z)) {
          std::cerr << "Unable to parse line: " << line << std::endl;
          return false;
        }

        if (node_idx == 0) {
          zero_based_index = true;
        }

        all_nodes.col(node_idx) = Vector3(x, y, z);
      }
    }

    if (zero_based_index) {
      vertex_positions = all_nodes.block(0, 0, 3, n_nodes);
    } else {
      vertex_positions = all_nodes.block(0, 1, 3, n_nodes);
    }

    return true;
  }

  static bool load_ele_file(const std::string &filename,
                            Eigen::Matrix4Xi &tet_vtx_indices) {
    std::ifstream ifile(filename);
    if (!ifile.is_open()) {
      std::cerr << "Error while opening file " << filename << std::endl;
      return false;
    }

    bool has_meta_info = false;

    Eigen::Matrix4Xi all_tets;
    bool zero_based_index = false;
    int n_tets = -1;

    std::string line;
    int n_read_tets = 0;
    while (std::getline(ifile, line)) {
      std::string::size_type pos = line.find_first_not_of(' ');
      if (pos == std::string::npos) {
        continue;
      }

      // Check for comment line
      else if (line.at(pos) == '#') {
        continue;
      }

      std::stringstream line_stream(line);

      if (!has_meta_info) {
        int nodes_per_tet, num_attrib;
        if (!(line_stream >> n_tets >> nodes_per_tet >> num_attrib)) {
          std::cerr << "Unable to parse line: " << line << std::endl;
          return false;
        }

        assert(nodes_per_tet == 4);
        assert(n_tets > 0);
        all_tets.resize(4, n_tets + 1);
        has_meta_info = true;
      } else {
        int tet_idx = -1;
        int v1, v2, v3, v4;
        if (!(line_stream >> tet_idx >> v1 >> v2 >> v3 >> v4)) {
          std::cerr << "Unable to parse line: " << line << std::endl;
          return false;
        }

        if (tet_idx == 0) {
          zero_based_index = true;
        }

        all_tets.col(tet_idx) = Eigen::Vector4i(v1, v2, v3, v4);
        n_read_tets++;
      }
    }

    assert(n_read_tets == n_tets);

    if (zero_based_index) {
      tet_vtx_indices = all_tets.block(0, 0, 4, n_tets);
    } else {
      tet_vtx_indices = all_tets.block(0, 1, 4, n_tets);
    }

    return true;
  }

  static void build_boundary_triangles(const Eigen::Matrix4Xi &tet_vtx_indices,
                                       Eigen::Matrix3Xi &boundary_faces) {
    typedef std::array<int, 3> HalfFaceVertexIndex;
    typedef std::pair<HalfFaceVertexIndex, int> HalfFaceInfo;

    std::vector<HalfFaceVertexIndex> half_face_vertex_indces;  // List of all half-faces with correct orientation
    std::vector<HalfFaceInfo> half_faces_info;

    int n_tets = tet_vtx_indices.cols();
    for (int i = 0; i < n_tets; ++i) {
      int v1 = tet_vtx_indices(0, i), v2 = tet_vtx_indices(1, i), v3 =
          tet_vtx_indices(2, i), v4 = tet_vtx_indices(3, i);

      // Counter-clockwise orientation for each face
      HalfFaceVertexIndex f1 { { v3, v2, v1 } };
      HalfFaceVertexIndex f2 { { v2, v3, v4 } };
      HalfFaceVertexIndex f3 { { v1, v4, v3 } };
      HalfFaceVertexIndex f4 { { v4, v1, v2 } };

      half_faces_info.push_back(
          HalfFaceInfo(f1, half_face_vertex_indces.size()));
      half_face_vertex_indces.push_back(f1);

      half_faces_info.push_back(
          HalfFaceInfo(f2, half_face_vertex_indces.size()));
      half_face_vertex_indces.push_back(f2);

      half_faces_info.push_back(
          HalfFaceInfo(f3, half_face_vertex_indces.size()));
      half_face_vertex_indces.push_back(f3);

      half_faces_info.push_back(
          HalfFaceInfo(f4, half_face_vertex_indces.size()));
      half_face_vertex_indces.push_back(f4);
    }

    // Sort the face vertex indices for each half-face, to prepare for half-face sorting
    int n_half_faces = half_faces_info.size();
    for (int i = 0; i < n_half_faces; ++i) {
      std::sort(half_faces_info[i].first.begin(),
                half_faces_info[i].first.end());
    }

    // Sort all half-faces
    std::sort(half_faces_info.begin(), half_faces_info.end());

    // Find out non-duplicate face vertex indices
    std::vector<int> boundary_faces_addr;
    for (int i = 0; i < n_half_faces; ++i) {
      if ((i - 1 < 0 || half_faces_info[i - 1].first != half_faces_info[i].first)
          && (i + 1 >= n_half_faces
              || half_faces_info[i + 1].first != half_faces_info[i].first)) {
        boundary_faces_addr.push_back(half_faces_info[i].second);
      }
    }

    int n_boundary_faces = boundary_faces_addr.size();
    boundary_faces.resize(3, n_boundary_faces);
    for (int i = 0; i < n_boundary_faces; ++i) {
      HalfFaceVertexIndex &current_face =
          half_face_vertex_indces[boundary_faces_addr[i]];
      boundary_faces.col(i) = Eigen::Vector3i(current_face[0], current_face[1],
                                              current_face[2]);
    }
  }
};

#endif /* TETMESHIO_H_ */
