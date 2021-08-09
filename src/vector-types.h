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

class Point2D;
class Point2DDubins;
class Point3D;

#include <math.h>
#include <cassert>

class Vector {
 public:
  Vector(int capacity);
  virtual ~Vector();

  virtual double* GetRawCoords();
  double operator[](int i) const;
  virtual void Set(int pos, double val);
  friend Vector operator*(Vector v1, const double scale);
  friend double operator|(Vector v1, Vector v2);
  friend Vector operator*(Vector v1, Vector v2);
  friend Vector operator-(const Vector v1, const Vector v2);
  friend Vector operator+(const Vector v1, const Vector v2);
  virtual void operator+=(const Vector &v);
  virtual void operator-=(const Vector &v);
  friend bool operator==(const Vector &a, const Vector &b);
  friend bool operator<(const Vector &a, const Vector &b);
  virtual Vector Inverse(const Vector &vec);
  virtual double Size() const;
  virtual void Normalize();
 protected:
  int capacity;
  double *coords;
};

class Quaternion : public Vector {
  public:
    Quaternion();
    Quaternion(double w, double x, double y, double z);
    Quaternion(double yaw, double pitch, double roll);
    friend Quaternion operator*(Quaternion a, Quaternion b);
    Quaternion Inverse() const;
    void ToRotationMatrix(double (&matrix)[3][3]);
    void Combine(const Quaternion &q);
    double Distance(const Quaternion &b) const;
};

#endif
