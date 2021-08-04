/**
 * @file point-types.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __POINT_TYPES_H__
#define __POINT_TYPES_H__

template <class T> class Point2D;
template <class T> class Point2DDubins;
template <class T> class Point3D;

#include <regex>
#include <math.h>
#include "constants.h"
#include "vector-types.h"

template<class T>
class Point2D {
  public:
    Point2D<T>() : coords{0, 0} {
    }

    Point2D<T>(T x, T y) : coords{x, y} {
    }

    Point2D<T>(const std::string &s, T scale=1) {
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

    Point2D<T>& operator=(const Point2D<T> &p) {
      if(*this != p) {
        for (int i{0}; i < 2; ++i) {
          coords[i] = p[i];
        }
      }

      return *this;
    }

    void SetPosition(T x, T y) {
      coords[0] = x;
      coords[1] = y;
    }

    const T* GetPosition() const {
      return coords;
    }

    const T operator[](int i) const {
      if (i < 2) {
        return coords[i];
      } else {
        return 1;
      }
    }

    inline void operator+=(const Vector<T> &translate) {
      for (int i{0}; i < 2; ++i) {
        coords[i] += translate[i];
      }
    }

    friend bool operator==(const Point2D<T> &p1, const Point2D<T> &p2) {
      return p1.x() == p2.x() && p1.y() == p2.y();
    }

    friend bool operator<(const Point2D<T> &p1, const Point2D<T> &p2) {
      return p1.x() < p2.x() ||
            (p1.x() == p2.x() && p1.y() < p2.y());
    }

    // scale position
    friend Point2D<T> operator/(const Point2D<T> &p1, const T scale) {
      Point2D<T> newPoint{p1};
      for (int i{0}; i < 2; ++i) {
        newPoint.coords[i] /= scale;
      }
      
      return newPoint;
    }

    T Distance(const Point2D<T> &other) const {
      T sum{0};
      for (int i{0}; i < 2; ++i) {
        T diff{(*this)[i] - other[i]};
        sum += diff * diff;
      }

      return sqrt(sum);
    }

    Point2D<T> GetStateInDistance(Point2D<T> &other, T dist) const {
      T realDist{Distance(other)};
      PointVector2D<T> direction(*this, other);
      Point2D<T> retVal;
      T ratio{dist / realDist};

      retVal.SetPosition(coords[0] + direction[0] * ratio, coords[1] + direction[1] * ratio);

      return retVal;
    }

    void PrintPosition(std::ostream &out) {
      out << (*this)[0] << DELIMITER_OUT << (*this)[1];
    }

  private:
    T coords[2];
};

template<class T>
class Point2DDubins {
  public:
    Point2DDubins<T>() : coords{0, 0}, phi{0} {
    }

    Point2DDubins<T>(T x, T y) : coords{x, y} {
    }

    Point2DDubins<T>(const std::string &s, T scale=1) {
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

    Point2DDubins<T>& operator=(const Point2DDubins<T> &p) {
      if(*this != p) {
        for (int i{0}; i < 2; ++i) {
          coords[i] = p[i];
        }
      }

      phi = p.GetAngle();

      return *this;
    }

    void SetAngle(T angle) {
      this->phi = angle;
    }

    T GetAngle() const {
      return phi;
    }

    void SetPosition(T x, T y) {
      coords[0] = x;
      coords[1] = y;
    }

    const T* GetPosition() const {
      return coords;
    }

    const T operator[](int i) const {
      if (i < 2) {
        return coords[i];
      } else if (i == 3) {
        return phi;
      } else {
        return 1;
      }
    }

    inline void operator+=(const Vector<T> &translate) {
      for (int i{0}; i < 2; ++i) {
        coords[i] += translate[i];
      }
    }

    friend bool operator==(const Point2DDubins<T> &p1, const Point2DDubins<T> &p2) {
      bool equal{true};
      for (int i{0}; i < 3; ++i) {
        equal &= (p1[i] == p2[i]);
      }

      return equal;
    }

    friend bool operator<(const Point2DDubins<T> &p1, const Point2DDubins<T> &p2) {
      bool equal{true};
      int i{0};
      for (; i < 3 && equal; ++i) {
        equal &= (p1[i] == p2[i]);
      }

      return p1[i - 1] < p2[i - 1];
    }

    // scale position
    friend Point2DDubins<T> operator/(const Point2DDubins<T> &p1, const T scale) {
      Point2DDubins<T> newPoint{p1};
      for (int i{0}; i < 2; ++i) {
        newPoint.coords[i] /= scale;
      }
      
      return newPoint;
    }

    T Distance(const Point2DDubins<T> &other) const {
      T sum{0};
      for (int i{0}; i < 2; ++i) {
        T diff{(*this)[i] - other[i]};
        sum += diff * diff;
      }

      T diff{AngleDifference(GetAngle(), other.GetAngle())};
      sum += diff * diff;

      return sqrt(sum);
    }

    Point2D<T> GetStateInDistance(Point2DDubins<T> &other, T dist) const {
      T realDist{Distance(other)};
      PointVector2D<T> direction(*this, other);
      Point2DDubins<T> retVal;
      T ratio{dist / realDist};

      retVal.SetPosition(coords[0] + direction[0] * ratio, coords[1] + direction[1] * ratio);
      
      T angleDirection{GetAngle() - other.GetAngle()};
      if (GetAngle() >= other.GetAngle()) {
        if (angleDirection < M_PI) {
          angleDirection *= -1;
        } else {
          angleDirection = 2 * M_PI - angleDirection;
        }
      } else {
        if (angleDirection < -M_PI) {
          angleDirection = -2 * M_PI - angleDirection;
        } else {
          angleDirection *= -1; 
        }
      }
      retVal.SetAngle(GetAngle() + angleDirection * ratio);

      return retVal;
    }

    void PrintPosition(std::ostream &out) {
      out << (*this)[0] << DELIMITER_OUT << (*this)[1];
    }

  private:
    T coords[2];
    T phi;
};

template<class T>
class Point3D {
  public:
    Point3D<T>() : coords{0, 0, 0} {
    }

    Point3D<T>(T x, T y, T z, T yaw, T pitch, T roll) : coords{x, y, z}, rotation{yaw, pitch, roll} {
    }

    Point3D<T>(const std::string &s, T scale=1) {
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
      rotation = Quaternion<T>(m[4], m[5], m[6]);
    }

    Quaternion<T> GetRotation() const {
      return rotation;
    }

    void SetRotation(Quaternion<T> q) {
      rotation = q;
    }

    void Set(T x, T y, T z, T yaw, T pitch, T roll) {
      coords[0] = x;
      coords[1] = y;
      coords[2] = z;
      rotation = Quaternion<T>(yaw, pitch, roll);
    }

    void SetPosition(T x, T y, T z) {
      coords[0] = x;
      coords[1] = y;
      coords[2] = z;
    }

    const T* GetPosition() const {
      return coords;
    }

    const T operator[](int i) const {
      if (i < 3) {
        return coords[i];
      } else if (i < 6) {
        return rotation[i - 3];
      } else {
        return 1;
      }
    }

    inline void operator+=(const Vector<T> &translate) {
      for (int i{0}; i < 3; ++i) {
        coords[i] += translate[i];
      }
    }

    friend bool operator==(const Point3D<T> &p1, const Point3D<T> &p2) {
      return p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1.rotation == p2.rotation;
    }

    friend bool operator<(const Point3D<T> &p1, const Point3D<T> &p2) {
      return p1[0] < p2[0] ||
            (p1[0] == p2[0] && p1[1] < p2[1]) ||
            (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] < p2[2]) ||
            (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1.rotation < p2.rotation);
    }

    // scale position (NOT the rotation)
    friend Point3D<T> operator/(const Point3D<T> &p1, const T scale) {
      Point3D<T> newPoint{p1};
      for (int i{0}; i < 3; ++i) {
        newPoint.coords[i] /= scale;
      }
      
      return newPoint;
    }

    T Distance(const Point3D<T> &other) const {
      T sum{0};
      for (int i{0}; i < 3; ++i) {
        T diff{(*this)[i] - other[i]};
        sum += diff * diff;
      }
      
      T diff{rotation.Distance(other.GetRotation())};
      sum += diff * diff;

      return sqrt(sum);
    }

    // Interpolation of rotation done according to https://ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf
    Point3D<T> GetStateInDistance(Point3D<T> &other, T dist) const {
      T realDist{Distance(other)};
      PointVector3D<T> direction(*this, other);
      Point3D<T> retVal;
      T ratio{dist / realDist};

      retVal.SetPosition(coords[0] + direction[0] * ratio, coords[1] + direction[1] * ratio, coords[2] + direction[2] * ratio);

      Quaternion<T> first{this->GetRotation()};
      Quaternion<T> second{other.GetRotation()};
      double product{first | second}; // inner product of quaternions
      if (product < 0) {
        second = second.Inverse();
        second.Set(0, -second[0]);
      }

      Quaternion<T> resultRot{first[0] + ratio * (second[0] - first[0]), first[1] + ratio * (second[1] - first[1]), first[2] + ratio * (second[2] - first[2]), first[3] + ratio * (second[3] - first[3])};
      resultRot.Normalize();
      retVal.SetRotation(resultRot);
      return retVal;
    }

    void PrintPosition(std::ostream &out) {
      out << (*this)[0] << DELIMITER_OUT << (*this)[1] << DELIMITER_OUT << (*this)[2];
    }

  private:
    T coords[3];
    Quaternion<T> rotation;
};

#endif
