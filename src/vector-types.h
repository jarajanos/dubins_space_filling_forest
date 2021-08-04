#pragma once
/**
 * @file vector-types.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __VECTOR_TYPES_H__
#define __VECTOR_TYPES_H__

#include <math.h>

template <class T>
class Vector {
 public:
  Vector(int capacity) {
    this->capacity = capacity;
    coords = new T[capacity];
  }

  ~Vector() {
    delete[] coords;
  }

  inline T* operator()() {
    return coords;
  }

  inline T operator[](int i) const {
    return coords[i];
  }

  void Set(int pos, T val) {
    this->coords[pos] = val;
  }

  friend Vector<T> operator*(Vector<T> v1, const T scale) {
    Vector<T> retVal{v1.capacity};
    for (int i{0}; i < v1.capacity; ++i) {
      retVal.coords[i] = v1.coords[i] * scale;
    }

    return retVal;
  }

  friend T operator|(Vector<T> v1, Vector<T> v2) {
    assert(v1.capacity == v2.capacity);
    double accum{0};

    for (int i{0}; i < v1.capacity; ++i) {
      accum += v1[i] * v2[i];
    }

    return accum;
  }

  friend Vector<T> operator*(Vector<T> v1, Vector<T> v2) {
    assert((v1.capacity == 3 && v2.capacity == 3));
    T  cache[3];
    cache[0] = v1[1] * v2[2] - v1[2] * v2[1];
    cache[1] = v1[2] * v2[0] - v1[0] * v2[2];
    cache[2] = v1[0] * v2[1] - v1[1] * v2[0];

    Vector<T> retVal{cache[0], cache[1], cache[2]};
    return retVal;
  }

  friend Vector<T> operator-(const Vector<T> v1, const Vector<T> v2) {
    assert(v1.capacity == v2.capacity);
    Vector<T> retVal{v1.capacity};
    for (int i{0}; i < v1.capacity; ++i) {
      retVal.coords[i] = v1[i] - v2[i];
    }

    return retVal;
  }

  friend Vector<T> operator+(const Vector<T> v1, const Vector<T> v2) {
    assert(v1.capacity == v2.capacity);
    Vector<T> retVal{v1.capacity};
    for (int i{0}; i < v1.capacity; ++i) {
      retVal.coords[i] = v1[i] + v2[i];
    }

    return retVal;
  }

  void operator+=(const Vector<T> &v) {
    assert(this->capacity == v.capacity);
    for (int i{0}; i < this->capacity; ++i) {
      this->coords[i] += v[i];
    }
  }

  void operator-=(const Vector<T> &v) {
    assert(this->capacity == v.capacity);
    for (int i{0}; i < this->capacity; ++i) {
      this->coords[i] -= v[i];
    }
  }

  friend bool operator==(const Vector<T> &a, const Vector<T> &b) {
    if (a.capacity != b.capacity) {
      return false;
    }

    bool correct{true};
    for (int i{0}; i < a.capacity && correct; ++i) {
      correct &= (a[i] == b[i]);
    }

    return correct;
  }

  friend bool operator<(const Vector<T> &a, const Vector<T> &b) {
    assert(a.capacity == b.capacity);
    bool equal{true};
    int i{0};

    for ( ; i < a.capacity && equal; ++i) {
      equal = (a[i] == b[i]);
    }

    return a[i - 1] < b[i - 1];
  }

  inline Vector<T> Inverse(const Vector<T> &vec) {
    Vector<T> retVal(this->capacity);
    for (int i{0}; i < this->capacity; ++i) {
      retVal.coords[i] = -vec[i];
    }

    return retVal;
  } 

  T Size() const {
    T retVal{0};

    for (int i{0}; i < this->capacity; ++i) {
      retVal += coords[i] * coords[i];
    }

    return sqrt(retVal);
  }

  void Normalize() {
    T s{this->Size()};
    for (int i{0}; i < this->capacity; ++i) {
      coords[i] /= s;
    }
  }

 protected:
  int capacity;
  T *coords;
};

template<class T>
class PointVector3D : public Vector<T> {
  public:
    PointVector3D(T x, T y, T z) : Vector<T>(3) {
      this->coords[0] = x;
      this->coords[1] = y;
      this->coords[2] = z;
    }

    PointVector3D(Point3D<T> p) : PointVector3D<T>(p[0], p[1], p[2]) {
    }

    PointVector3D(Point3D<T> p1, Point3D<T> p2) : Vector<T>(3) {
      for (int i{0}; i < 3; ++i) {
        this->coords[i] = p2[i] - p1[i];
      }
    }
};

template<class T>
class PointVector2D : public Vector<T> {
  public:
    PointVector2D(T x, T y) : PointVector2D<T>(2) {
      this->coords[0] = x;
      this->coords[1] = y;
    }

    PointVector2D(Point2D<T> p) : PointVector2D<T>(p[0], p[1]) {
    }

    PointVector2D(Point2DDubins<T> p) : PointVector2D<T>(p[0], p[1]) {
    }

    PointVector2D(Point2D<T> p1, Point2D<T> p2) : PointVector2D<T>(2) {
      for (int i{0}; i < 2; ++i) {
        this->coords[i] = p2[i] - p1[i];
      }
    }

    PointVector2D(Point2DDubins<T> p1, Point2DDubins<T> p2) : PointVector2D<T>(2) {
      for (int i{0}; i < 2; ++i) {
        this->coords[i] = p2[i] - p1[i];
      }
    }
};

template<class T>
class Quaternion : public Vector<T> {
  public:
    Quaternion() : Vector<T>(4) {
    }

    Quaternion(T w, T x, T y, T z) : Vector<T>(4) {
      this->coords[0] = w;
      this->coords[1] = x;
      this->coords[2] = y;
      this->coords[3] = z;
    }

    Quaternion(T yaw, T pitch, T roll) : Vector<T>(4) {
      T c_yaw = cos(yaw * 0.5);
      T s_yaw = sin(yaw * 0.5);
      T c_pitch = cos(pitch * 0.5);
      T s_pitch = sin(pitch * 0.5);
      T c_roll = cos(roll * 0.5);
      T s_roll = sin(roll * 0.5);

      this->coords[0] = c_roll * c_pitch * c_yaw + s_roll * s_pitch * s_yaw;
      this->coords[1] = s_roll * c_pitch * c_yaw - c_roll * s_pitch * s_yaw;
      this->coords[2] = c_roll * s_pitch * c_yaw + s_roll * c_pitch * s_yaw;
      this->coords[3] = c_roll * c_pitch * s_yaw - s_roll * s_pitch * c_yaw;
    }

    friend Quaternion<T> operator*(Quaternion<T> a, Quaternion<T> b) {
      Quaternion<T> retVal;
      retVal.coords[0] = a.coords[0] * b.coords[0] - a.coords[1] * b.coords[1] - a.coords[2] * b.coords[2] - a.coords[3] * b.coords[3];
      retVal.coords[1] = a.coords[0] * b.coords[1] + a.coords[1] * b.coords[0] + a.coords[2] * b.coords[3] - a.coords[3] * b.coords[2];
      retVal.coords[2] = a.coords[0] * b.coords[2] - a.coords[1] * b.coords[3] + a.coords[2] * b.coords[0] + a.coords[3] * b.coords[1];
      retVal.coords[3] = a.coords[0] * b.coords[3] + a.coords[1] * b.coords[2] - a.coords[2] * b.coords[1] + a.coords[3] * b.coords[0];

      return retVal;
    }

    Quaternion<T> Inverse() const {
      Quaternion<T> retVal;

      retVal.coords[0] = this->coords[0];
      for (int i{1}; i < 4; ++i) {
        retVal.coords[i] = -this->coords[i];
      }

      return retVal;
    }

    void ToRotationMatrix(T (&matrix)[3][3]) {
      matrix[0][0] = 2 * this->coords[0] * this->coords[0] + 2 * this->coords[1] * this->coords[1] - 1;
      matrix[0][1] = 2 * this->coords[1] * this->coords[2] - 2 * this->coords[0] * this->coords[3];
      matrix[0][2] = 2 * this->coords[0] * this->coords[2] + 2 * this->coords[1] * this->coords[3];
      matrix[1][0] = 2 * this->coords[0] * this->coords[3] + 2 * this->coords[1] * this->coords[2];
      matrix[1][1] = 2 * this->coords[0] * this->coords[0] + 2 * this->coords[2] * this->coords[2] - 1;
      matrix[1][2] = 2 * this->coords[2] * this->coords[3] - 2 * this->coords[0] * this->coords[1];
      matrix[2][0] = 2 * this->coords[1] * this->coords[3] - 2 * this->coords[0] * this->coords[2];
      matrix[2][1] = 2 * this->coords[2] * this->coords[3] + 2 * this->coords[0] * this->coords[1];
      matrix[2][2] = 2 * this->coords[0] * this->coords[0] + 2 * this->coords[3] * this->coords[3] - 1;
    }

    void Combine(const Quaternion<T> &q) {
      T tempCoords[4];
    
      tempCoords[0] = q.coords[0] * this->coords[0] - q.coords[1] * this->coords[1] - q.coords[2] * this->coords[2] - q.coords[3] * this->coords[3];
      tempCoords[1] = q.coords[0] * this->coords[1] + q.coords[1] * this->coords[0] + q.coords[2] * this->coords[3] - q.coords[3] * this->coords[2];
      tempCoords[2] = q.coords[0] * this->coords[2] - q.coords[1] * this->coords[3] + q.coords[2] * this->coords[0] + q.coords[3] * this->coords[1];
      tempCoords[3] = q.coords[0] * this->coords[3] + q.coords[1] * this->coords[2] - q.coords[2] * this->coords[1] + q.coords[3] * this->coords[0];
    
      for (int i{0}; i < 4; ++i) {
        this->coords[i] = tempCoords[i];
      }
    }

    void RotatePoint(Point3D<T> &p) {
      Quaternion<T> q{0, p[0], p[1], p[2]};
      Quaternion<T> temp = *this * q;
      q = temp * Inverse();

      p.set(q[1], q[2], q[3]);
    }

    T Distance(const Quaternion<T> &b) const {
      Quaternion<T> diff{Inverse() * b};
      Vector<T> reduced(3);
      reduced.Set(0, diff[1]);
      reduced.Set(1, diff[2]);
      reduced.Set(2, diff[3]);

      return 2 * atan2(reduced.Size(), diff[0]);
    }
};

#endif
