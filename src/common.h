#pragma once
/**
 * @file primitives.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <math.h>
#include <vector>
#include <deque>
#include <limits.h>
#include <flann/flann.hpp>
#include "point-types.h"
#include "vector-types.h"
#include "constants.h"

#define PROBLEM_DIMENSION   static_cast<size_t>(this->problem.dimension)

#define MIN(X, Y) ((X < Y) ? (X) : (Y))
#define MAX(X, Y) ((X > Y) ? (X) : (Y))

struct FileStruct;

int ParseString(std::string &inp, std::string &outp1, std::string &outp2, std::string &delimiter);
FileStruct PrefixFileName(const FileStruct &path, const std::string &insert);
std::string Ltrim(const std::string &s);
std::string Rtrim(const std::string &s);
std::string Trim(const std::string &s);

enum Dimensions {
  D2 = 2,
  D2Dubins = 3,
  D3 = 6,
  D3Dubins = 7
};

enum FileType {
  Map,
  Obj
};

enum SaveOptions {
  None = 0,
  SaveGoals = 1,
  SaveTree = 2,
  SaveRaw = 4,
  SaveSmooth = 8,
  SaveParams = 16,
  SaveTSP = 32,
  SaveFrontiers = 64,
  Invalid = 128
};

enum SolverType{
  SFF,
  RRT,
  Lazy
};

struct FileStruct {
  std::string fileName;
  FileType type;
};

template <class T>
struct Range {
  T minX;
  T maxX;
  T minY;
  T maxY;
  T minZ;
  T maxZ;
};

// FLANN FUNCTOR
template<class T>
struct D6Distance {
  typedef bool is_vector_space_distance;

  typedef T ElementType;
  typedef typename flann::Accumulator<T>::Type ResultType;

  template <typename Iterator1, typename Iterator2>
  ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType worst_dist = -1) const {
    ResultType result = ResultType();
    ResultType diff;

    for (int i{0}; i < 3; ++i) {
      diff = (ResultType)(*a++ - *b++);
      result = diff * diff;
    }

    Quaternion<T> q1, q2;
    for (int i{3}; i < 7; ++i) {
      q1.Set(*a++);
      q2.Set(*b++);
    }

    diff = q1.Distance(q2);
    result += diff * diff;
    return result;
  }
};

template<class T>
class SymmetricMatrix {
  public:
    SymmetricMatrix(const int size) {
      holder.resize(size * (size+1) / 2);
      this->size = size;
    }

    T& operator()(int i, int j) {
      int index;
      if (i <= j) {
        index = i * size - (i - 1) * i / 2 + j - i;
      } else {
        index = j * size - (j - 1) * j / 2 + i - j;
      }
      return holder[index];
    }

    const bool Exists(int i, int j) {
      return this->operator()(i, j).Exists();
    }
  private:
    std::deque<T> holder;
    int size;
};

template<class T, class R>
struct DistanceHolder {
  R *node1;
  R *node2;
  T distance;
  std::deque<R *> plan;

  DistanceHolder() : node1{NULL}, node2{NULL}, distance{std::numeric_limits<T>::max()} {
  }

  DistanceHolder(R *first, R *second) : node1{first}, node2{second} {
    if (*first < *second) {
      node1 = first;
      node2 = second;
    } else {
      node1 = second;
      node2 = first;
    }
    distance = first->DistanceToRoot + second->DistanceToRoot + first->Position.distance(second->Position);
  }

  DistanceHolder(R *first, R *second, T dist) : distance{dist} {
    if (*first < *second) {
      node1 = first;
      node2 = second;
    } else {
      node1 = second;
      node2 = first;
    }
  }

  DistanceHolder(R *first, R *second, T dist, std::deque<R *> &plan) : distance{dist}, plan{plan} {
    if (*first < *second) {
      node1 = first;
      node2 = second;
    } else {
      node1 = second;
      node2 = first;
      std::reverse(this->plan.begin(), this->plan.end());
    }
  }

  friend bool operator<(const DistanceHolder<T, R> &l, const DistanceHolder<T, R> &r) {
    return l.distance < r.distance;
  }  

  friend bool operator==(const DistanceHolder<T, R> &l, const DistanceHolder<T, R> &r) {
    return l.node1 == r.node1 && l.node2 == r.node2;
  }

  const bool Exists() const {
    return node1 != nullptr;
  }

  void UpdateDistance() {
    distance = node1->DistanceToRoot + node2->DistanceToRoot + node1->Position.distance(node2->Position);
  }
};

template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <class T>
T NormalizeAngle(T angle) {
  if (angle < -M_PI) {
    return angle + 2 * M_PI;
  } else if (angle >= M_PI) {
    return angle - 2 * M_PI;
  } else {
    return angle;
  }
}

template <class T, class R>
T AngleDifference(T a1, R a2) {
  T diff{(T)a2 - a1};
  return NormalizeAngle(diff);
}

#endif /* __COMMON_H__ */
