/**
 * @file point_types.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "point-types.h"

Point2D::Point2D() : coords{0, 0} {
}

Point2D::Point2D(double x, double y) : coords{x, y} {
}

Point2D::Point2D(const std::string &s, double scale) {
  std::regex r("\\[(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*)\\]");
  std::smatch m;
  std::regex_search(s, m, r);
  if (m.size() != 3) {
    throw std::invalid_argument("Unknown format of point");
  }

  //position
  for (int i{0}; i < 2; ++i) {
    coords[i] = std::stod(m[i + 1]) * scale;
  }
}

Point2D& Point2D::operator=(const Point2D &p) {
  if(*this != p) {
    for (int i{0}; i < 2; ++i) {
      coords[i] = p[i];
    }
  }

  return *this;
}

void Point2D::SetPosition(double x, double y) {
  coords[0] = x;
  coords[1] = y;
}

const double* Point2D::GetPosition() const {
  return coords;
}

const double* Point2D::GetRawCoords() const {
  return coords;
}

const double Point2D::operator[](int i) const {
  if (i < 2) {
    return coords[i];
  } else {
    return 1;
  }
}

void Point2D::operator+=(const Vector &translate) {
  for (int i{0}; i < 2; ++i) {
    coords[i] += translate[i];
  }
}

bool operator==(const Point2D &p1, const Point2D &p2) {
  return p1[0] == p2[0] && p1[1] == p2[1];
}

bool operator!=(Point2D &p1, const Point2D &p2) {
  return !(p1 == p2);
}

bool operator<(const Point2D &p1, const Point2D &p2) {
  return p1[0] < p2[0] ||
        (p1[0] == p2[0] && p1[1] < p2[1]);
}

// scale position
Point2D operator/(const Point2D &p1, const double scale) {
  Point2D newPoint{p1};
  for (int i{0}; i < 2; ++i) {
    newPoint.coords[i] /= scale;
  }
  
  return newPoint;
}

double Point2D::Distance(const Point2D &other) const {
  double sum{0};
  for (int i{0}; i < 2; ++i) {
    double diff{(*this)[i] - other[i]};
    sum += diff * diff;
  }

  return sqrt(sum);
}

Point2D Point2D::GetStateInDistance(Point2D &other, double dist) const {
  double realDist{this->Distance(other)};
  PointVector2D direction(*this, other);
  Point2D retVal;
  double ratio{dist / realDist};

  retVal.SetPosition(coords[0] + direction[0] * ratio, coords[1] + direction[1] * ratio);

  return retVal;
}

void Point2D::FillRotationMatrix(double (&matrix)[3][3]) const {
  matrix[0][0] = 1;
  matrix[0][1] = 0;
  matrix[0][2] = 0;
  matrix[1][0] = 0;
  matrix[1][1] = 1;
  matrix[1][2] = 0;
  matrix[2][0] = 0;
  matrix[2][1] = 0;
  matrix[2][2] = 1;
}

void Point2D::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1];
}

Point2DDubins::Point2DDubins() : coords{0, 0}, phi{0} {
}

Point2DDubins::Point2DDubins(double x, double y) : coords{x, y} {
}

Point2DDubins::Point2DDubins(const std::string &s, double scale) {
  std::regex r("\\[(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*)\\]");
  std::smatch m;
  std::regex_search(s, m, r);
  if (m.size() != 4) {
    throw std::invalid_argument("Unknown format of point");
  }

  //position
  for (int i{0}; i < 2; ++i) {
    coords[i] = std::stod(m[i + 1]) * scale;
  }

  //angle
  phi = std::stod(m[3]);
}

Point2DDubins::Point2DDubins(opendubins::State dubinsState) : coords{dubinsState.getPoint().getX(), dubinsState.getPoint().getY()}, phi{dubinsState.getAngle()} {
}

Point2DDubins& Point2DDubins::operator=(const Point2DDubins &p) {
  if(*this != p) {
    for (int i{0}; i < 2; ++i) {
      coords[i] = p[i];
    }
  }

  phi = p.GetAngle();

  return *this;
}

void Point2DDubins::SetAngle(double angle) {
  this->phi = angle;
}

double Point2DDubins::GetAngle() const {
  return phi;
}

void Point2DDubins::SetPosition(double x, double y) {
  coords[0] = x;
  coords[1] = y;
}

const double* Point2DDubins::GetPosition() const {
  return coords;
}

const double* Point2DDubins::GetRawCoords() const {
  return coords;
}

const double Point2DDubins::operator[](int i) const {
  if (i < 2) {
    return coords[i];
  } else if (i == 3) {
    return phi;
  } else {
    return 1;
  }
}

void Point2DDubins::operator+=(const Vector &translate) {
  for (int i{0}; i < 2; ++i) {
    coords[i] += translate[i];
  }
}

bool operator==(const Point2DDubins &p1, const Point2DDubins &p2) {
  bool equal{true};
  for (int i{0}; i < 3; ++i) {
    equal &= (p1[i] == p2[i]);
  }

  return equal;
}

bool operator!=(Point2DDubins &p1, const Point2DDubins &p2) {
  return !(p1 == p2);
}

bool operator<(const Point2DDubins &p1, const Point2DDubins &p2) {
  bool equal{true};
  int i{0};
  for (; i < 3 && equal; ++i) {
    equal &= (p1[i] == p2[i]);
  }

  return p1[i - 1] < p2[i - 1];
}

// scale position
Point2DDubins operator/(const Point2DDubins &p1, const double scale) {
  Point2DDubins newPoint{p1};
  for (int i{0}; i < 2; ++i) {
    newPoint.coords[i] /= scale;
  }
  
  return newPoint;
}

double Point2DDubins::Distance(const Point2DDubins &other) const {
  opendubins::State a{coords[0], coords[1], GetAngle()};
  opendubins::State b{other[0], other[1], other.GetAngle()};
  opendubins::Dubins dubPath{a, b, DubinsRadius};

  return dubPath.length;
}

Point2DDubins Point2DDubins::GetStateInDistance(Point2DDubins &other, double dist) const {
  double realDist{Distance(other)};
  PointVector2D direction(*this, other);
  Point2DDubins retVal;
  double ratio{dist / realDist};

  retVal.SetPosition(coords[0] + direction[0] * ratio, coords[1] + direction[1] * ratio);
  retVal.SetAngle(GetAngle() + AngleDifference(this->GetAngle(), other.GetAngle()) * ratio);

  return retVal;
}

void Point2DDubins::FillRotationMatrix(double (&matrix)[3][3]) const {
  matrix[0][0] = cos(phi);
  matrix[0][1] = -sin(phi);
  matrix[0][2] = 0;
  matrix[1][0] = sin(phi);
  matrix[1][1] = cos(phi);
  matrix[1][2] = 0;
  matrix[2][0] = 0;
  matrix[2][1] = 0;
  matrix[2][2] = 1;
}

void Point2DDubins::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1];
}

Point3D::Point3D() : coords{0, 0, 0} {
}

Point3D::Point3D(double x, double y, double z, double yaw, double pitch, double roll) : coords{x, y, z}, rotation{yaw, pitch, roll} {
}

Point3D::Point3D(const std::string &s, double scale) {
  std::regex r("\\[(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*)\\]");
  std::smatch m;
  std::regex_search(s, m, r);
  if (m.size() != 7) {
    throw std::invalid_argument("Unknown format of point");
  }

  //position
  for (int i{0}; i < 3; ++i) {
    coords[i] = std::stod(m[i + 1]) * scale;
  }

  // rotation
  rotation = Quaternion(std::stod(m[4]), std::stod(m[5]), std::stod(m[6]));
}

Quaternion Point3D::GetRotation() const {
  return rotation;
}

void Point3D::SetRotation(Quaternion q) {
  rotation = q;
}

void Point3D::Set(double x, double y, double z, double yaw, double pitch, double roll) {
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
  rotation = Quaternion(yaw, pitch, roll);
}

void Point3D::SetPosition(double x, double y, double z) {
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}

const double* Point3D::GetPosition() const {
  return coords;
}

const double* Point3D::GetRawCoords() const {
  return coords;
}

const double Point3D::operator[](int i) const {
  if (i < 3) {
    return coords[i];
  } else if (i < 6) {
    return rotation[i - 3];
  } else {
    return 1;
  }
}

void Point3D::operator+=(const Vector &translate) {
  for (int i{0}; i < 3; ++i) {
    coords[i] += translate[i];
  }
}

bool operator==(const Point3D &p1, const Point3D &p2) {
  return p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1.rotation == p2.rotation;
}

bool operator!=(const Point3D &p1, const Point3D &p2) {
  return !(p1 == p2);
}

bool operator<(const Point3D &p1, const Point3D &p2) {
  return p1[0] < p2[0] ||
        (p1[0] == p2[0] && p1[1] < p2[1]) ||
        (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] < p2[2]) ||
        (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1.rotation < p2.rotation);
}

// scale position (NOT the rotation)
Point3D operator/(const Point3D &p1, const double scale) {
  Point3D newPoint{p1};
  for (int i{0}; i < 3; ++i) {
    newPoint.coords[i] /= scale;
  }
  
  return newPoint;
}

double Point3D::Distance(const Point3D &other) const {
  double sum{0};
  for (int i{0}; i < 3; ++i) {
    double diff{(*this)[i] - other[i]};
    sum += diff * diff;
  }
  
  double diff{rotation.Distance(other.GetRotation())};
  sum += diff * diff;

  return sqrt(sum);
}

// Interpolation of rotation done according to https://ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf
Point3D Point3D::GetStateInDistance(Point3D &other, double dist) const {
  double realDist{Distance(other)};
  PointVector3D direction(*this, other);
  Point3D retVal;
  double ratio{dist / realDist};

  retVal.SetPosition(coords[0] + direction[0] * ratio, coords[1] + direction[1] * ratio, coords[2] + direction[2] * ratio);

  Quaternion first{this->GetRotation()};
  Quaternion second{other.GetRotation()};
  double product{first | second}; // inner product of quaternions
  if (product < 0) {
    second = second.Inverse();
    second.Set(0, -second[0]);
  }

  Quaternion resultRot{first[0] + ratio * (second[0] - first[0]), first[1] + ratio * (second[1] - first[1]), first[2] + ratio * (second[2] - first[2]), first[3] + ratio * (second[3] - first[3])};
  resultRot.Normalize();
  retVal.SetRotation(resultRot);
  return retVal;
}

void Point3D::FillRotationMatrix(double (&matrix)[3][3]) const {
  this->GetRotation().ToRotationMatrix(matrix);
}

Point3D Point3D::RotatePoint(Quaternion rotation) {
  Point3D retVal;
  Quaternion q{0, (*this)[0], (*this)[1], (*this)[2]};
  
  Quaternion temp = rotation * q;
  q = temp * rotation.Inverse();

  retVal.SetPosition(q[1], q[2], q[3]);
  return retVal;
}

void Point3D::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1] << DELIMITER_OUT << (*this)[2];
}

PointVector3D::PointVector3D() : Vector(3) {
  this->coords[0] = 0;
  this->coords[1] = 0;
  this->coords[2] = 0;     
}

PointVector3D::PointVector3D(double x, double y, double z) : Vector(3) {
  this->coords[0] = x;
  this->coords[1] = y;
  this->coords[2] = z;
}

PointVector3D::PointVector3D(Point3D p) : PointVector3D(p[0], p[1], p[2]) {
}

PointVector3D::PointVector3D(Point3D p1, Point3D p2) : Vector(3) {
  for (int i{0}; i < 3; ++i) {
    this->coords[i] = p2[i] - p1[i];
  }
}

PointVector2D::PointVector2D() : Vector(2) {
  this->coords[0] = 0;
  this->coords[1] = 0;
}

PointVector2D::PointVector2D(double x, double y) : Vector(2) {
  this->coords[0] = x;
  this->coords[1] = y;
}

PointVector2D::PointVector2D(Point2D p) : PointVector2D(p[0], p[1]) {
}

PointVector2D::PointVector2D(Point2DDubins p) : PointVector2D(p[0], p[1]) {
}

PointVector2D::PointVector2D(Point2D p1, Point2D p2) : Vector(2) {
  for (int i{0}; i < 2; ++i) {
    this->coords[i] = p2[i] - p1[i];
  }
}

PointVector2D::PointVector2D(Point2DDubins p1, Point2DDubins p2) : Vector(2) {
  for (int i{0}; i < 2; ++i) {
    this->coords[i] = p2[i] - p1[i];
  }
}

PointVector3D PointVector2D::To3DVector() const {
  return PointVector3D(this->coords[0], this->coords[1], 0);
}

// STREAM OUTPUTS
std::ostream& operator<<(std::ostream &out, const Point2D &p) {
  return out << p[0] << DELIMITER_OUT << p[1];
}

std::ostream& operator<<(std::ostream &out, const Point2DDubins &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << p.GetAngle();
}

std::ostream& operator<<(std::ostream &out, const Point3D &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << p[2] << DELIMITER_OUT << p.GetRotation().GetYaw() << DELIMITER_OUT << p.GetRotation().GetPitch() << DELIMITER_OUT << p.GetRotation().GetRoll();
}
