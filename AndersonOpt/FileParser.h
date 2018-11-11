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

#ifndef FILEPARSER_H_
#define FILEPARSER_H_

#include <string>
#include <iostream>
#include <fstream>

inline bool load_value(const std::string &str, double &value) {
  try {
    value = std::stod(str);
  } catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return false;
  } catch (const std::out_of_range &oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    return false;
  }

  return true;
}

inline bool load_value(const std::string &str, float &value) {
  try {
    value = std::stof(str);
  } catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return false;
  } catch (const std::out_of_range &oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    return false;
  }

  return true;
}

inline bool load_value(const std::string &str, int &value) {
  try {
    value = std::stoi(str);
  } catch (const std::invalid_argument& ia) {
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return false;
  } catch (const std::out_of_range &oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
    return false;
  }

  return true;
}

template<typename T>
bool read_sequence_file(const char *filename, std::vector<T> &values) {
  std::ifstream ifile(filename);
  if (!ifile.is_open()) {
    std::cerr << "Error while opening file " << filename << std::endl;
    return false;
  }

  values.clear();
  std::string line;
  while (std::getline(ifile, line)) {
    T val;
    if (load_value(line, val)) {
      values.push_back(val);
    } else {
      break;
    }
  }

  std::cout << "Read " << values.size() << " numbers from " << filename
      << std::endl;

  return true;
}

template<typename T>
bool read_sequence_file(const char *filename, std::vector<T> &values_x,
                        std::vector<T> &values_y) {
  std::ifstream ifile(filename);
  if (!ifile.is_open()) {
    std::cerr << "Error while opening file " << filename << std::endl;
    return false;
  }

  values_x.clear();
  values_y.clear();
  while (ifile) {
    T val, val2;
    if (ifile >> val >> val2) {
      values_x.push_back(val);
      values_y.push_back(val2);
    } else {
      break;
    }
  }

  std::cout << "Read " << values_x.size() << " numbers from " << filename
      << std::endl;

  return true;
}

template<typename T>
bool read_sequence_file(const char *filename, std::vector<T> &values_x,
                        std::vector<T> &values_y, std::vector<T> &values_z) {
  std::ifstream ifile(filename);
  if (!ifile.is_open()) {
    std::cerr << "Error while opening file " << filename << std::endl;
    return false;
  }

  values_x.clear();
  values_y.clear();
  values_z.clear();
  while (ifile) {
    T val, val2, val3;
    if (ifile >> val >> val2 >> val3) {
      values_x.push_back(val);
      values_y.push_back(val2);
      values_z.push_back(val3);
    } else {
      break;
    }
  }

  std::cout << "Read " << values_x.size() << " numbers from " << filename
      << std::endl;

  return true;
}

#endif /* FILEPARSER_H_ */
