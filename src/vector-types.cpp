/**
 * @file vector-types.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "vector-types.h"

Vector::Vector(int capacity) {
  this->capacity = capacity;
  coords = new double[capacity];
}

Vector::~Vector() {
  delete[] coords;
}

double* Vector::GetRawCoords() {
  return coords;
}

double Vector::operator[](int i) const {
  return coords[i];
}

void Vector::Set(int pos, double val) {
  this->coords[pos] = val;
}

Vector operator*(Vector v1, const double scale) {
  Vector retVal{v1.capacity};
  for (int i{0}; i < v1.capacity; ++i) {
    retVal.coords[i] = v1.coords[i] * scale;
  }

  return retVal;
}

double operator|(Vector v1, Vector v2) {
  assert(v1.capacity == v2.capacity);
  double accum{0};

  for (int i{0}; i < v1.capacity; ++i) {
    accum += v1[i] * v2[i];
  }

  return accum;
}

Vector operator*(Vector v1, Vector v2) {
  assert((v1.capacity == 3 && v2.capacity == 3));
  Vector retVal(3);
  retVal.Set(0, v1[1] * v2[2] - v1[2] * v2[1]);
  retVal.Set(1, v1[2] * v2[0] - v1[0] * v2[2]);
  retVal.Set(2, v1[0] * v2[1] - v1[1] * v2[0]);

  return retVal;
}

Vector operator-(const Vector v1, const Vector v2) {
  assert(v1.capacity == v2.capacity);
  Vector retVal{v1.capacity};
  for (int i{0}; i < v1.capacity; ++i) {
    retVal.coords[i] = v1[i] - v2[i];
  }

  return retVal;
}

Vector operator+(const Vector v1, const Vector v2) {
  assert(v1.capacity == v2.capacity);
  Vector retVal{v1.capacity};
  for (int i{0}; i < v1.capacity; ++i) {
    retVal.coords[i] = v1[i] + v2[i];
  }

  return retVal;
}

void Vector::operator+=(const Vector &v) {
  assert(this->capacity == v.capacity);
  for (int i{0}; i < this->capacity; ++i) {
    this->coords[i] += v[i];
  }
}

void Vector::operator-=(const Vector &v) {
  assert(this->capacity == v.capacity);
  for (int i{0}; i < this->capacity; ++i) {
    this->coords[i] -= v[i];
  }
}

bool operator==(const Vector &a, const Vector &b) {
  if (a.capacity != b.capacity) {
    return false;
  }

  bool correct{true};
  for (int i{0}; i < a.capacity && correct; ++i) {
    correct &= (a[i] == b[i]);
  }

  return correct;
}

bool operator<(const Vector &a, const Vector &b) {
  assert(a.capacity == b.capacity);
  bool equal{true};
  int i{0};

  for ( ; i < a.capacity && equal; ++i) {
    equal = (a[i] == b[i]);
  }

  return a[i - 1] < b[i - 1];
}

Vector Vector::Inverse(const Vector &vec) {
  Vector retVal(this->capacity);
  for (int i{0}; i < this->capacity; ++i) {
    retVal.coords[i] = -vec[i];
  }

  return retVal;
} 

double Vector::Size() const {
  double retVal{0};

  for (int i{0}; i < this->capacity; ++i) {
    retVal += coords[i] * coords[i];
  }

  return sqrt(retVal);
}

void Vector::Normalize() {
  double s{this->Size()};
  for (int i{0}; i < this->capacity; ++i) {
    coords[i] /= s;
  }
}

Quaternion::Quaternion() : Vector(4) {
}

Quaternion::Quaternion(double w, double x, double y, double z) : Vector(4) {
  this->coords[0] = w;
  this->coords[1] = x;
  this->coords[2] = y;
  this->coords[3] = z;
}

Quaternion::Quaternion(double yaw, double pitch, double roll) : Vector(4) {
  double c_yaw = cos(yaw * 0.5);
  double s_yaw = sin(yaw * 0.5);
  double c_pitch = cos(pitch * 0.5);
  double s_pitch = sin(pitch * 0.5);
  double c_roll = cos(roll * 0.5);
  double s_roll = sin(roll * 0.5);

  this->coords[0] = c_roll * c_pitch * c_yaw + s_roll * s_pitch * s_yaw;
  this->coords[1] = s_roll * c_pitch * c_yaw - c_roll * s_pitch * s_yaw;
  this->coords[2] = c_roll * s_pitch * c_yaw + s_roll * c_pitch * s_yaw;
  this->coords[3] = c_roll * c_pitch * s_yaw - s_roll * s_pitch * c_yaw;
}

Quaternion operator*(Quaternion a, Quaternion b) {
  Quaternion retVal;
  retVal.coords[0] = a.coords[0] * b.coords[0] - a.coords[1] * b.coords[1] - a.coords[2] * b.coords[2] - a.coords[3] * b.coords[3];
  retVal.coords[1] = a.coords[0] * b.coords[1] + a.coords[1] * b.coords[0] + a.coords[2] * b.coords[3] - a.coords[3] * b.coords[2];
  retVal.coords[2] = a.coords[0] * b.coords[2] - a.coords[1] * b.coords[3] + a.coords[2] * b.coords[0] + a.coords[3] * b.coords[1];
  retVal.coords[3] = a.coords[0] * b.coords[3] + a.coords[1] * b.coords[2] - a.coords[2] * b.coords[1] + a.coords[3] * b.coords[0];

  return retVal;
}

Quaternion Quaternion::Inverse() const {
  Quaternion retVal;

  retVal.coords[0] = this->coords[0];
  for (int i{1}; i < 4; ++i) {
    retVal.coords[i] = -this->coords[i];
  }

  return retVal;
}

void Quaternion::ToRotationMatrix(double (&matrix)[3][3]) {
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

void Quaternion::Combine(const Quaternion &q) {
  double tempCoords[4];

  tempCoords[0] = q.coords[0] * this->coords[0] - q.coords[1] * this->coords[1] - q.coords[2] * this->coords[2] - q.coords[3] * this->coords[3];
  tempCoords[1] = q.coords[0] * this->coords[1] + q.coords[1] * this->coords[0] + q.coords[2] * this->coords[3] - q.coords[3] * this->coords[2];
  tempCoords[2] = q.coords[0] * this->coords[2] - q.coords[1] * this->coords[3] + q.coords[2] * this->coords[0] + q.coords[3] * this->coords[1];
  tempCoords[3] = q.coords[0] * this->coords[3] + q.coords[1] * this->coords[2] - q.coords[2] * this->coords[1] + q.coords[3] * this->coords[0];

  for (int i{0}; i < 4; ++i) {
    this->coords[i] = tempCoords[i];
  }
}

double Quaternion::Distance(const Quaternion &b) const {
  Quaternion diff{Inverse() * b};
  Vector reduced(3);
  reduced.Set(0, diff[1]);
  reduced.Set(1, diff[2]);
  reduced.Set(2, diff[3]);

  return 2 * atan2(reduced.Size(), diff[0]);
}
