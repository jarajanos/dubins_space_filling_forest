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
#include "common.h"

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

  phi = p.GetHeading();

  return *this;
}

void Point2DDubins::SetHeading(double angle) {
  this->phi = angle;
}

void Point2DDubins::SetHeading(int angleId, int angleResolution) {
  this->phi = NormalizeAngle(this->phi + (angleId * 2 * M_PI) / angleResolution);
}

double Point2DDubins::GetHeading() const {
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
  } else if (i == 2) {
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
    equal &= inBounds(p1[i], p2[i], EQ_TOLERANCE);
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
  opendubins::State a{coords[0], coords[1], GetHeading()};
  opendubins::State b{other[0], other[1], other.GetHeading()};
  opendubins::Dubins dubPath{a, b, DubinsRadius};

  return dubPath.length;
}

Point2DDubins Point2DDubins::GetStateInDistance(Point2DDubins &other, double dist) const {
  opendubins::State a{coords[0], coords[1], GetHeading()};
  opendubins::State b{other[0], other[1], other.GetHeading()};
  opendubins::Dubins dubPath{a, b, DubinsRadius};

  return Point2DDubins(dubPath.getState(dist));
}

std::deque<Point2DDubins> Point2DDubins::SampleDubinsPathTo(const Point2DDubins &other, double dist) {
  std::deque<Point2DDubins> retVal;
  opendubins::State startDub{coords[0], coords[1], this->GetHeading()};
  opendubins::State finishDub{other[0], other[1], other.GetHeading()};
  opendubins::Dubins pathDub{startDub, finishDub, DubinsRadius};
  
  double pathDist{pathDub.length};
  double parts{pathDist / dist};

  retVal.emplace_back(*this);
  for (int index{1}; index < parts; ++index) {
    opendubins::State temp{pathDub.getState(index * pathDist / parts)};
    retVal.emplace_back(temp);
  }

  return retVal;
}

Point2DDubins Point2DDubins::GetInvertedPoint() {
  Point2DDubins retVal{*this};
  retVal.SetHeading(NormalizeAngle(this->GetHeading() + M_PI));

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

Point2DPolynom::Point2DPolynom() : Point2DPolynom(0,0) {
}

Point2DPolynom::Point2DPolynom(double x, double y) : coords{x, y}, velocity{0, 0}, acceleration{0, 0} {
}

Point2DPolynom::Point2DPolynom(const std::string &s, double scale) {
  std::regex r("\\[(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*)\\]");
  std::smatch m;
  std::regex_search(s, m, r);
  if (m.size() != 7) {
    throw std::invalid_argument("Unknown format of point");
  }

  //position
  for (int i{0}; i < 2; ++i) {
    coords[i] = std::stod(m[i + 1]) * scale;
    velocity[i] = std::stod(m[i + 3]) * scale;
    acceleration[i] = std::stod(m[i + 5]) * scale;
  }
}

Point2DPolynom::Point2DPolynom(const Vec2 pos) : Point2DPolynom(pos.x, pos.y) {

}

void Point2DPolynom::SetPosition(double x, double y) {
  coords[0] = x;
  coords[1] = y;
}

void Point2DPolynom::SetPosition(Vec2 vec) {
  coords[0] = vec.x;
  coords[1] = vec.y;
}

void Point2DPolynom::SetVelocity(double x, double y) {
  velocity[0] = x;
  velocity[1] = y;
}

void Point2DPolynom::SetVelocity(Vec2 vec) {
  velocity[0] = vec.x;
  velocity[1] = vec.y;
}

void Point2DPolynom::SetAcceleration(double x, double y) {
  acceleration[0] = x;
  acceleration[1] = y;
}

void Point2DPolynom::SetAcceleration(Vec2 vec) {
  acceleration[0] = vec.x;
  acceleration[1] = vec.y;
}

const double* Point2DPolynom::GetPosition() const {
  return coords;
}

const double* Point2DPolynom::GetRawCoords() const {
  return coords;
}

const double Point2DPolynom::operator[](int i) const {
  if (i < 2) {
    return coords[i];
  } else if (i < 4) {
    return velocity[i - 2];
  } else if (i < 6) {
    return acceleration[i - 4];
  }

  return 1;
}

void Point2DPolynom::operator+=(const Vector &translate) {
  for (int i{0}; i < 2; ++i) {
    coords[i] += translate[i];
  }
}

bool operator==(const Point2DPolynom &p1, const Point2DPolynom &p2) {
  bool equal{true};
  for (int i{0}; i < 6 && equal; ++i) {
    equal &= inBounds(p1[i], p2[i], EQ_TOLERANCE);
  }
  return equal;
}

bool operator!=(const Point2DPolynom &p1, const Point2DPolynom &p2) {
  return !(p1 == p2);
}

bool operator<(const Point2DPolynom &p1, const Point2DPolynom &p2) {
  for (int i{0}; i < 6; ++i) {
    if (p1[i] == p2[i]) {
      continue;
    }

    return p1[i] < p2[i];
  }
  return false;
}

Point2DPolynom operator/(const Point2DPolynom &p1, const double scale) {
  Point2DPolynom newPoint{p1};
  for (int i{0}; i < 2; ++i) {
    newPoint.coords[i] /= scale;
  }
  
  return newPoint;
}

double Point2DPolynom::Distance(const Point2DPolynom &other) const {
  double totalDist{EuclideanDistance(other)};
  double totalTime{totalDist / AverageVelocity};

  Vec2 pos0{(*this)[0], (*this)[1]};
  Vec2 vel0{(*this)[2], (*this)[3]};
  Vec2 acc0{(*this)[4], (*this)[5]};

  Vec2 posf{other[0], other[1]};
  Vec2 velf{other[2], other[3]};
  Vec2 accf{other[4], other[5]};

  Vec2 gravity{0,-Gravity};
  RapidTrajectoryGenerator2D traj(pos0, vel0, acc0, gravity);
  traj.SetGoalPosition(posf);
  traj.SetGoalVelocity(velf);
  traj.SetGoalAcceleration(accf);

  traj.Generate(totalTime);

  return traj.GetCost();
}

// Euclidean 6D distance
double Point2DPolynom::EuclideanDistance(const Point2DPolynom &other) const {
  double distance{0};
  for (int i{0}; i < 6; ++i) {
    distance += POW((*this)[i] - other[i]);
  }
  
  distance = sqrt(distance);
  return distance;
}

Point2DPolynom Point2DPolynom::GetStateInDistance(Point2DPolynom &other, double dist) const {
  Point2DPolynom retVal;

  Point2D start{(*this)[0], (*this)[1]};
  Point2D finish{other[0], other[1]};
  double totalDist{start.Distance(finish)};
  double totalTime{totalDist / AverageVelocity};

  Vec2 pos0{(*this)[0], (*this)[1]};
  Vec2 vel0{(*this)[2], (*this)[3]};
  Vec2 acc0{(*this)[4], (*this)[5]};

  Vec2 posf{other[0], other[1]};
  Vec2 velf{other[2], other[3]};
  Vec2 accf{other[4], other[5]};

  Vec2 gravity{0,-this->Gravity};
  RapidTrajectoryGenerator2D traj(pos0, vel0, acc0, gravity);
  traj.SetGoalPosition(posf);
  traj.SetGoalVelocity(velf);
  traj.SetGoalAcceleration(accf);

  traj.Generate(totalTime);

  retVal.SetPosition(traj.GetPosition(dist / AverageVelocity));
  retVal.SetVelocity(traj.GetVelocity(dist / AverageVelocity));
  retVal.SetAcceleration(traj.GetAcceleration(dist / AverageVelocity));

  return retVal;
}

std::deque<Point2DPolynom> Point2DPolynom::SampleTrajectory(Point2DPolynom &other, double interval) {
  std::deque<Point2DPolynom> retVal;

  Point2D start{GetPositionOnly()};
  Point2D finish{other.GetPositionOnly()};
  double totalDist{start.Distance(finish)};
  double totalTime{totalDist / AverageVelocity};

  Vec2 pos0{(*this)[0], (*this)[1]};
  Vec2 vel0{(*this)[2], (*this)[3]};
  Vec2 acc0{(*this)[4], (*this)[5]};

  Vec2 posf{other[0], other[1]};
  Vec2 velf{other[2], other[3]};
  Vec2 accf{other[4], other[5]};

  Vec2 gravity{0,-this->Gravity};
  RapidTrajectoryGenerator2D traj(pos0, vel0, acc0, gravity);
  traj.SetGoalPosition(posf);
  traj.SetGoalVelocity(velf);
  traj.SetGoalAcceleration(accf);

  traj.Generate(totalTime);

  double parts{totalTime / interval};
  for (int i{0}; i < parts; ++i) {
    Vec2 pos{traj.GetPosition(i * interval)};
    retVal.emplace_back(pos);
  }

  return retVal;
}

void Point2DPolynom::FillRotationMatrix(double (&matrix)[3][3]) const {
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

void Point2DPolynom::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1] << DELIMITER_OUT << "0";
}

Point2D Point2DPolynom::GetPositionOnly() {
  return Point2D((*this)[0], (*this)[1]);
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

// for conversion of 3D Dubins problem to 3D for LazySFF 
Point3D::Point3D(const Point3DDubins &point) : coords{point[0], point[1], point[2]}, rotation{point.GetHeading(), point.GetPitch(), 0} {
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
  retVal.SetRotation(this->rotation * rotation);
  return retVal;
}

void Point3D::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1] << DELIMITER_OUT << (*this)[2];
}

Point3DDubins::Point3DDubins(opendubins::State3D dubinsState) 
  : coords{dubinsState.getPoint().getX(), dubinsState.getPoint().getY(), dubinsState.getPoint().getZ()},
    yaw{dubinsState.getHeading()}, pitch{dubinsState.getPitch()} {
}

Point3DDubins::Point3DDubins() : coords{0, 0, 0} {
}

Point3DDubins::Point3DDubins(double x, double y, double z, double yaw, double pitch) : coords{x, y, z}, yaw{yaw}, pitch{pitch} {
}

Point3DDubins::Point3DDubins(const std::string &s, double scale) {
  std::regex r("\\[(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*)\\]");
  std::smatch m;
  std::regex_search(s, m, r);
  if (m.size() != 6) {
    throw std::invalid_argument("Unknown format of point");
  }

  //position
  for (int i{0}; i < 3; ++i) {
    coords[i] = std::stod(m[i + 1]) * scale;
  }

  // rotation
  this->yaw = std::stod(m[4]);
  this->pitch = std::stod(m[5]);
}

Point3DDubins::Point3DDubins(const Point3D &point) : coords{point[0], point[1], point[2]} {
  this->yaw = point.GetRotation().GetYaw();
  this->pitch = point.GetRotation().GetPitch();
}

void Point3DDubins::SetHeading(double yaw) {
  this->yaw = yaw;
}

void Point3DDubins::SetHeading(int angleId, int angleResolution) {
  SetHeading(yaw + (angleId * 2 * M_PI) / angleResolution);
}

double Point3DDubins::GetHeading() const {
  return yaw;
}

void Point3DDubins::SetPitch(double pitch) {
  this->pitch = pitch;
}

double Point3DDubins::GetPitch() const {
  return pitch;
}

void Point3DDubins::Set(double x, double y, double z, double yaw, double pitch) {
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
  this->yaw = yaw;
  this->pitch = pitch;
}

void Point3DDubins::SetPosition(double x, double y, double z) {
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}

const double* Point3DDubins::GetPosition() const {
  return coords;
}

const double* Point3DDubins::GetRawCoords() const {
  return coords;
}

const double Point3DDubins::operator[](int i) const {
  if (i < 3) {
    return coords[i];
  } else if (i == 3) {
    return yaw;
  } else if (i == 4) {
    return pitch;
  } else {
    return 1;
  }
}

void Point3DDubins::operator+=(const Vector &translate) {
  for (int i{0}; i < 3; ++i) {
    coords[i] += translate[i];
  }
}

bool operator==(const Point3DDubins &p1, const Point3DDubins &p2) {
  bool equal{true};
  for (int i{0}; i < 5 && equal; ++i) {
    equal &= inBounds(p1[i], p2[i], EQ_TOLERANCE);
  }
  return equal;
}

bool operator!=(const Point3DDubins &p1, const Point3DDubins &p2) {
  return !(p1 == p2);
}

bool operator<(const Point3DDubins &p1, const Point3DDubins &p2) {
  return p1[0] < p2[0] ||
        (p1[0] == p2[0] && p1[1] < p2[1]) ||
        (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] < p2[2]) ||
        (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1.yaw < p2.yaw) ||
        (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1.yaw == p2.yaw && p1.pitch < p2.pitch);
}

// scale position (NOT the rotation)
Point3DDubins operator/(const Point3DDubins &p1, const double scale) {
  Point3DDubins newPoint{p1};
  for (int i{0}; i < 3; ++i) {
    newPoint.coords[i] /= scale;
  }
  
  return newPoint;
}

double Point3DDubins::Distance(const Point3DDubins &other) const {
  bool equal{true};
  for (int i{0}; i < 3; ++i) {
    equal &= inBounds((*this)[i], other[i], EQ_TOLERANCE);
  }

  if (equal) {
    return 0;
  }

  opendubins::State3D startDub{coords[0], coords[1], coords[2], GetHeading(), GetPitch()};
  opendubins::State3D finishDub{other[0], other[1], other[2], other.GetHeading(), other.GetPitch()};
  opendubins::Dubins3D pathDub{startDub, finishDub, DubinsRadius, -MaxPitch, MaxPitch};
  
  return pathDub.getLength();
}

Point3DDubins Point3DDubins::GetStateInDistance(Point3DDubins &other, double dist) const {
  opendubins::State3D startDub{coords[0], coords[1], coords[2], GetHeading(), GetPitch()};
  opendubins::State3D finishDub{other[0], other[1], other[2], other.GetHeading(), other.GetPitch()};
  opendubins::Dubins3D pathDub{startDub, finishDub, DubinsRadius, -MaxPitch, MaxPitch};
  
  if (pathDub.getLength() == std::numeric_limits<double>::max()) {
    return Point3DDubins();
  }

  opendubins::State3D temp{pathDub.getState(dist)};
  return Point3DDubins(temp);
}

std::deque<Point3DDubins> Point3DDubins::SampleDubinsPathTo(const Point3DDubins &other, double dist) {
  std::deque<Point3DDubins> retVal;
  opendubins::State3D startDub{coords[0], coords[1], coords[2], GetHeading(), GetPitch()};
  opendubins::State3D finishDub{other[0], other[1], other[2], other.GetHeading(), other.GetPitch()};
  opendubins::Dubins3D pathDub{startDub, finishDub, DubinsRadius, -MaxPitch, MaxPitch};

  double pathDist{pathDub.length};
  if (pathDist == std::numeric_limits<double>::max()) {
    return retVal;
  }

  double parts{pathDist / dist};
  retVal.emplace_back(*this);
  for (int index{1}; index < parts; ++index) {
    opendubins::State3D temp{pathDub.getState(index * pathDist / parts)};
    retVal.emplace_back(temp);
  }

  return retVal;
}

void Point3DDubins::FillRotationMatrix(double (&matrix)[3][3]) const {
  matrix[0][0] = cos(yaw) * cos(pitch);
  matrix[0][1] = -sin(yaw);
  matrix[0][2] = cos(yaw) * sin(pitch);
  matrix[1][0] = sin(yaw) * cos(pitch);
  matrix[1][1] = cos(yaw);
  matrix[1][2] = sin(yaw) * sin(pitch);
  matrix[2][0] = -sin(pitch);
  matrix[2][1] = 0;
  matrix[2][2] = cos(pitch);
}

Point3DDubins Point3DDubins::RotatePoint(double deltaYaw, double deltaPitch) {
  Point3DDubins pos{*this};
  pos.SetHeading(NormalizeAngle(this->yaw + deltaYaw));
  pos.SetPitch(NormalizeAngle(this->pitch + deltaPitch));

  return pos;
}

Point3DDubins Point3DDubins::GetInvertedPoint() {
  Point3DDubins pos{RotatePoint(M_PI, 0)};
  pos.SetPitch(-this->pitch);
  return pos;
}

void Point3DDubins::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1] << DELIMITER_OUT << (*this)[2];
}

Point3DPolynom::Point3DPolynom() : Point3DPolynom(0,0,0) {
}

Point3DPolynom::Point3DPolynom(double x, double y, double z) : coords{x, y, z}, acceleration{0,0,0}, velocity{0,0,0} {
}

Point3DPolynom::Point3DPolynom(const std::string &s, double scale) {
  std::regex r("\\[(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*);\\s*(\\-?[\\d]+[\\.]?[\\d]*)\\]");
  std::smatch m;
  std::regex_search(s, m, r);
  if (m.size() != 10) {
    throw std::invalid_argument("Unknown format of point");
  }

  //position
  for (int i{0}; i < 3; ++i) {
    coords[i] = std::stod(m[i + 1]) * scale;
    velocity[i] = std::stod(m[i + 4]) * scale;
    acceleration[i] = std::stod(m[i + 7]) * scale;
  }
}

Point3DPolynom::Point3DPolynom(const Vec3 pos) : Point3DPolynom(pos.x, pos.y, pos.z) {
}

void Point3DPolynom::SetPosition(double x, double y, double z) {
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}

void Point3DPolynom::SetPosition(Vec3 vec) {
  coords[0] = vec.x;
  coords[1] = vec.y;
  coords[2] = vec.z;
}

void Point3DPolynom::SetVelocity(double x, double y, double z) {
  velocity[0] = x;
  velocity[1] = y;
  velocity[2] = z;
}

void Point3DPolynom::SetVelocity(Vec3 vec) {
  velocity[0] = vec.x;
  velocity[1] = vec.y;
  velocity[2] = vec.z;
}
    
void Point3DPolynom::SetAcceleration(double x, double y, double z) {
  acceleration[0] = x;
  acceleration[1] = y;
  acceleration[2] = z;
}

void Point3DPolynom::SetAcceleration(Vec3 vec) {
  acceleration[0] = vec.x;
  acceleration[1] = vec.y;
  acceleration[2] = vec.z;
}

const double* Point3DPolynom::GetPosition() const {
  return coords;
}

const double* Point3DPolynom::GetRawCoords() const {
  return coords;
}

const double Point3DPolynom::operator[](int i) const {
  if (i < 3) {
    return coords[i];
  } else if (i < 6) {
    return velocity[i - 3];
  } else if (i < 9) {
    return acceleration[i - 6];
  }

  return 1;
}

void Point3DPolynom::operator+=(const Vector &translate) {
  for (int i{0}; i < 3; ++i) {
    coords[i] += translate[i];
  }
}

bool operator==(const Point3DPolynom &p1, const Point3DPolynom &p2) {
  bool equal{true};
  for (int i{0}; i < 9 && equal; ++i) {
    equal &= inBounds(p1[i], p2[i], EQ_TOLERANCE);
  }
  return equal;
}

bool operator!=(const Point3DPolynom &p1, const Point3DPolynom &p2) {
  return !(p1 == p2);
}

bool operator<(const Point3DPolynom &p1, const Point3DPolynom &p2) {
  for (int i{0}; i < 9; ++i) {
    if (p1[i] == p2[i]) {
      continue;
    }

    return p1[i] < p2[i];
  }
  return false;
}

Point3DPolynom operator/(const Point3DPolynom &p1, const double scale) {
  Point3DPolynom newPoint{p1};
  for (int i{0}; i < 3; ++i) {
    newPoint.coords[i] /= scale;
  }
  
  return newPoint;
}

// trajectory distance
double Point3DPolynom::Distance(const Point3DPolynom &other) const {
  double totalDist{EuclideanDistance(other)};
  double totalTime{totalDist / AverageVelocity};

  Vec3 pos0{(*this)[0], (*this)[1], (*this)[2]};
  Vec3 vel0{(*this)[3], (*this)[4], (*this)[5]};
  Vec3 acc0{(*this)[6], (*this)[7], (*this)[8]};

  Vec3 posf{other[0], other[1], other[2]};
  Vec3 velf{other[3], other[4], other[5]};
  Vec3 accf{other[6], other[7], other[8]};

  Vec3 gravity{0,0,-Gravity};
  RapidTrajectoryGenerator traj(pos0, vel0, acc0, gravity);
  traj.SetGoalPosition(posf);
  traj.SetGoalVelocity(velf);
  traj.SetGoalAcceleration(accf);

  traj.Generate(totalTime);

  return traj.GetCost();
}

// Euclidean 9D distance
double Point3DPolynom::EuclideanDistance(const Point3DPolynom &other) const {
  double distance{0};
  for (int i{0}; i < 9; ++i) {
    distance += POW((*this)[i] - other[i]);
  }
  
  distance = sqrt(distance);
  return distance;
}

Point3DPolynom Point3DPolynom::GetStateInDistance(Point3DPolynom &other, double dist) const {
  Point3DPolynom retVal;

  Point3D start{(*this)[0], (*this)[1], (*this)[2], 0, 0, 0};
  Point3D finish{other[0], other[1], other[2], 0, 0, 0};
  double totalDist{start.Distance(finish)};
  double totalTime{totalDist / AverageVelocity};

  Vec3 pos0{(*this)[0], (*this)[1], (*this)[2]};
  Vec3 vel0{(*this)[3], (*this)[4], (*this)[5]};
  Vec3 acc0{(*this)[6], (*this)[7], (*this)[8]};

  Vec3 posf{other[0], other[1], other[2]};
  Vec3 velf{other[3], other[4], other[5]};
  Vec3 accf{other[6], other[7], other[8]};

  Vec3 gravity{0,0,-this->Gravity};
  RapidTrajectoryGenerator traj(pos0, vel0, acc0, gravity);
  traj.SetGoalPosition(posf);
  traj.SetGoalVelocity(velf);
  traj.SetGoalAcceleration(accf);

  traj.Generate(totalTime);

  retVal.SetPosition(traj.GetPosition(dist / AverageVelocity));
  retVal.SetVelocity(traj.GetVelocity(dist / AverageVelocity));
  retVal.SetAcceleration(traj.GetAcceleration(dist / AverageVelocity));

  return retVal;
}

std::deque<Point3DPolynom> Point3DPolynom::SampleTrajectory(Point3DPolynom &other, double interval) {
  std::deque<Point3DPolynom> retVal;

  Point3D start{GetPositionOnly()};
  Point3D finish{other.GetPositionOnly()};
  double totalDist{start.Distance(finish)};
  double totalTime{totalDist / AverageVelocity};

  Vec3 pos0{(*this)[0], (*this)[1], (*this)[2]};
  Vec3 vel0{(*this)[3], (*this)[4], (*this)[5]};
  Vec3 acc0{(*this)[6], (*this)[7], (*this)[8]};

  Vec3 posf{other[0], other[1], other[2]};
  Vec3 velf{other[3], other[4], other[5]};
  Vec3 accf{other[6], other[7], other[8]};

  Vec3 gravity{0,0,-this->Gravity};
  RapidTrajectoryGenerator traj(pos0, vel0, acc0, gravity);
  traj.SetGoalPosition(posf);
  traj.SetGoalVelocity(velf);
  traj.SetGoalAcceleration(accf);

  traj.Generate(totalTime);

  double parts{totalTime / interval};
  for (int i{0}; i < parts; ++i) {
    Vec3 pos{traj.GetPosition(i * interval)};
    retVal.emplace_back(pos);
  }

  return retVal;
}

void Point3DPolynom::FillRotationMatrix(double (&matrix)[3][3]) const {
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

void Point3DPolynom::PrintPosition(std::ostream &out) {
  out << (*this)[0] << DELIMITER_OUT << (*this)[1] << DELIMITER_OUT << (*this)[2];
}

Point3D Point3DPolynom::GetPositionOnly() {
  return Point3D((*this)[0], (*this)[1], (*this)[2], 0, 0, 0);
}

PointVector3D::PointVector3D() : Vector(3) {
  this->coords.get()[0] = 0;
  this->coords.get()[1] = 0;
  this->coords.get()[2] = 0;     
}

PointVector3D::PointVector3D(double x, double y, double z) : Vector(3) {
  this->coords.get()[0] = x;
  this->coords.get()[1] = y;
  this->coords.get()[2] = z;
}

PointVector3D::PointVector3D(Point3D p) : PointVector3D(p[0], p[1], p[2]) {
}

PointVector3D::PointVector3D(Point3DDubins p) : PointVector3D(p[0], p[1], p[2]) {
}

PointVector3D::PointVector3D(Point3DPolynom p) : PointVector3D(p[0], p[1], p[2]) {
}

PointVector3D::PointVector3D(Point3D p1, Point3D p2) : Vector(3) {
  for (int i{0}; i < 3; ++i) {
    this->coords.get()[i] = p2[i] - p1[i];
  }
}

PointVector3D::PointVector3D(Point3DDubins p1, Point3DDubins p2) : Vector(3) {
  for (int i{0}; i < 3; ++i) {
    this->coords.get()[i] = p2[i] - p1[i];
  }
}

PointVector3D::PointVector3D(Point3DPolynom p1, Point3DPolynom p2) : Vector(3) {
  for (int i{0}; i < 3; ++i) {
    this->coords.get()[i] = p2[i] - p1[i];
  }
}

PointVector2D::PointVector2D() : Vector(2) {
  this->coords.get()[0] = 0;
  this->coords.get()[1] = 0;
}

PointVector2D::PointVector2D(double x, double y) : Vector(2) {
  this->coords.get()[0] = x;
  this->coords.get()[1] = y;
}

PointVector2D::PointVector2D(Point2D p) : PointVector2D(p[0], p[1]) {
}

PointVector2D::PointVector2D(Point2DDubins p) : PointVector2D(p[0], p[1]) {
}

PointVector2D::PointVector2D(Point2DPolynom p) : PointVector2D(p[0], p[1]) {
}

PointVector2D::PointVector2D(Point2D p1, Point2D p2) : Vector(2) {
  for (int i{0}; i < 2; ++i) {
    this->coords.get()[i] = p2[i] - p1[i];
  }
}

PointVector2D::PointVector2D(Point2DDubins p1, Point2DDubins p2) : Vector(2) {
  for (int i{0}; i < 2; ++i) {
    this->coords.get()[i] = p2[i] - p1[i];
  }
}

PointVector2D::PointVector2D(Point2DPolynom p1, Point2DPolynom p2) : Vector(2) {
  for (int i{0}; i < 2; ++i) {
    this->coords.get()[i] = p2[i] - p1[i];
  }
}

PointVector3D PointVector2D::To3DVector() const {
  return PointVector3D(this->coords.get()[0], this->coords.get()[1], 0);
}

// STREAM OUTPUTS
std::ostream& operator<<(std::ostream &out, const Point2D &p) {
  return out << p[0] << DELIMITER_OUT << p[1];
}

std::ostream& operator<<(std::ostream &out, const Point2DDubins &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << p.GetHeading();
}

std::ostream& operator<<(std::ostream &out, const Point2DPolynom &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << "0";
}

std::ostream& operator<<(std::ostream &out, const Point3D &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << p[2] << DELIMITER_OUT << p.GetRotation().GetYaw() << DELIMITER_OUT << p.GetRotation().GetPitch() << DELIMITER_OUT << p.GetRotation().GetRoll();
}

std::ostream& operator<<(std::ostream &out, const Point3DDubins &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << p[2] << DELIMITER_OUT << p.GetHeading() << DELIMITER_OUT << p.GetPitch() << DELIMITER_OUT << "0";
}

std::ostream& operator<<(std::ostream &out, const Point3DPolynom &p) {
  return out << p[0] << DELIMITER_OUT << p[1] << DELIMITER_OUT << p[2] << DELIMITER_OUT << "0" << DELIMITER_OUT << "0" << DELIMITER_OUT << "0";
}
