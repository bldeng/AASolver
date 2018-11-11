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

#ifndef OMPHELPER_H_
#define OMPHELPER_H_

#ifdef USE_OPENMP
#include <omp.h>
#ifdef USE_MSVC
#define OMP_PARALLEL __pragma(omp parallel)
#define OMP_FOR __pragma(omp for)
#define OMP_SINGLE __pragma(omp single)
#define OMP_SECTIONS __pragma(omp sections)
#define OMP_SECTION __pragma(omp section)
#else
#define OMP_PARALLEL _Pragma("omp parallel")
#define OMP_FOR _Pragma("omp for")
#define OMP_SINGLE _Pragma("omp single")
#define OMP_SECTIONS _Pragma("omp sections")
#define OMP_SECTION _Pragma("omp section")
#endif
#else
#include <ctime>
#define OMP_PARALLEL
#define OMP_FOR
#define OMP_SINGLE
#define OMP_SECTIONS
#define OMP_SECTION
#endif

#include <cassert>
#include <vector>

class Timer {
 public:

  typedef int EventID;

  EventID get_time() {
    EventID id = time_values_.size();

#ifdef USE_OPENMP
    time_values_.push_back(omp_get_wtime());
#else
    time_values_.push_back(clock());
#endif

    return id;
  }

  double elapsed_time(EventID event1, EventID event2) {
    assert(event1 >= 0 && event1 < static_cast<EventID>(time_values_.size()));
    assert(event2 >= 0 && event2 < static_cast<EventID>(time_values_.size()));

#ifdef USE_OPENMP
    return time_values_[event2] - time_values_[event1];
#else
    return double(time_values_[event2] - time_values_[event1]) / CLOCKS_PER_SEC;
#endif
  }

  void reset() {
    time_values_.clear();
  }

 private:
#ifdef USE_OPENMP
  std::vector<double> time_values_;
#else
  std::vector<clock_t> time_values_;
#endif
};

#endif /* OMPHELPER_H_ */
