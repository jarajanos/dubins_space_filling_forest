#pragma once
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

class Point2D;
class Point2DDubins;
class Point3D;
class Point3DDubins;

#include <regex>
#include <math.h>
#include "constants.h"
#include "vector-types.h"

#include "opendubins/dubins.h"
#include "opendubins/dubins3D.h"

#include "RapidTrajectoryGenerator.h"
#include "RapidTrajectoryGenerator2D.h"
#include "Vec3.h"
#include "Vec2.h"

using namespace RapidQuadrocopterTrajectoryGenerator;

class Point2D {
  public:
    Point2D();
    Point2D(double x, double y);
    Point2D(const std::string &s, double scale=1);

    Point2D& operator=(const Point2D &p);
    void SetPosition(double x, double y);
    const double* GetPosition() const;
    const double* GetRawCoords() const;
    const double operator[](int i) const;
    void operator+=(const Vector &translate);
    friend bool operator==(const Point2D &p1, const Point2D &p2);
    friend bool operator!=(Point2D &p1, const Point2D &p2);
    friend bool operator<(const Point2D &p1, const Point2D &p2);
    friend Point2D operator/(const Point2D &p1, const double scale);    // scale position
    double Distance(const Point2D &other) const;
    Point2D GetStateInDistance(Point2D &other, double dist) const;
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    void PrintPosition(std::ostream &out);
  protected:
    double coords[2];
};

class Point2DDubins {
  public:
    inline static double DubinsRadius;

    Point2DDubins();
    Point2DDubins(double x, double y);
    Point2DDubins(const std::string &s, double scale=1);
    Point2DDubins(opendubins::State dubinsState);

    Point2DDubins& operator=(const Point2DDubins &p);
    void SetHeading(double angle);
    double GetHeading() const;
    void SetHeading(int angleId, int angleResolution);
    void SetPosition(double x, double y);
    const double* GetPosition() const;
    const double* GetRawCoords() const;
    const double operator[](int i) const;
    void operator+=(const Vector &translate);
    friend bool operator==(const Point2DDubins &p1, const Point2DDubins &p2);
    friend bool operator!=(Point2DDubins &p1, const Point2DDubins &p2);
    friend bool operator<(const Point2DDubins &p1, const Point2DDubins &p2);
    friend Point2DDubins operator/(const Point2DDubins &p1, const double scale);    // scale position
    double Distance(const Point2DDubins &other) const;
    Point2DDubins GetStateInDistance(Point2DDubins &other, double dist) const;
    std::deque<Point2DDubins> SampleDubinsPathTo(const Point2DDubins &other, double dist);
    Point2DDubins GetInvertedPoint();
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    void PrintPosition(std::ostream &out);
  protected:
    double coords[2];
    double phi;
};

class Point2DPolynom {
  public:
    inline static double AverageVelocity;
    inline static double Gravity;

    Point2DPolynom();
    Point2DPolynom(double x, double y);
    Point2DPolynom(const std::string &s, double scale=1);
    Point2DPolynom(const Vec2 pos);

    void SetPosition(double x, double y);
    void SetPosition(Vec2 vec);
    void SetVelocity(double x, double y);
    void SetVelocity(Vec2 vec);
    void SetAcceleration(double x, double y);
    void SetAcceleration(Vec2 vec);
    const double* GetPosition() const;
    const double* GetRawCoords() const;
    const double operator[](int i) const;
    void operator+=(const Vector &translate);
    friend bool operator==(const Point2DPolynom &p1, const Point2DPolynom &p2);
    friend bool operator!=(const Point2DPolynom &p1, const Point2DPolynom &p2);
    friend bool operator<(const Point2DPolynom &p1, const Point2DPolynom &p2);
    friend Point2DPolynom operator/(const Point2DPolynom &p1, const double scale);    // scale position (NOT the rotation)
    double Distance(const Point2DPolynom &other) const;
    double EuclideanDistance(const Point2DPolynom &other) const;
    Point2DPolynom GetStateInDistance(Point2DPolynom &other, double dist) const;
    std::deque<Point2DPolynom> SampleTrajectory(Point2DPolynom &other, double interval);
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    void PrintPosition(std::ostream &out);

    Point2D GetPositionOnly();
  protected:
    double coords[2];
    double velocity[2];
    double acceleration[2];
};

class Point3D {
  public:
    Point3D();
    Point3D(double x, double y, double z, double yaw, double pitch, double roll);
    Point3D(const std::string &s, double scale=1);
    Point3D(const Point3DDubins &point);
    
    Quaternion GetRotation() const;
    void SetRotation(Quaternion q);
    void Set(double x, double y, double z, double yaw, double pitch, double roll);
    void SetPosition(double x, double y, double z);
    const double* GetPosition() const;
    const double* GetRawCoords() const;
    const double operator[](int i) const;
    void operator+=(const Vector &translate);
    friend bool operator==(const Point3D &p1, const Point3D &p2);
    friend bool operator!=(const Point3D &p1, const Point3D &p2);
    friend bool operator<(const Point3D &p1, const Point3D &p2);
    friend Point3D operator/(const Point3D &p1, const double scale);    // scale position (NOT the rotation)
    double Distance(const Point3D &other) const;
    Point3D GetStateInDistance(Point3D &other, double dist) const;      // Interpolation of rotation done according to https://ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    Point3D RotatePoint(Quaternion rotation);
    void PrintPosition(std::ostream &out);
  protected:
    double coords[3];
    Quaternion rotation;
};

class Point3DDubins {
  public:
    inline static double DubinsRadius;
    inline static double MaxPitch;

    Point3DDubins();
    Point3DDubins(double x, double y, double z, double yaw, double pitch);
    Point3DDubins(const std::string &s, double scale=1);
    Point3DDubins(opendubins::State3D dubinsState);
    Point3DDubins(const Point3D &point);

    void SetRotation(double yaw, double pitch);
    void SetHeading(double yaw);
    double GetHeading() const;
    void SetPitch(double pitch);
    double GetPitch() const;
    void SetHeading(int angleId, int angleResolution);
    void Set(double x, double y, double z, double yaw, double pitch);
    void SetPosition(double x, double y, double z);
    const double* GetPosition() const;
    const double* GetRawCoords() const;
    const double operator[](int i) const;
    void operator+=(const Vector &translate);
    friend bool operator==(const Point3DDubins &p1, const Point3DDubins &p2);
    friend bool operator!=(const Point3DDubins &p1, const Point3DDubins &p2);
    friend bool operator<(const Point3DDubins &p1, const Point3DDubins &p2);
    friend Point3DDubins operator/(const Point3DDubins &p1, const double scale);    // scale position (NOT the rotation)
    double Distance(const Point3DDubins &other) const;
    Point3DDubins GetStateInDistance(Point3DDubins &other, double dist) const;
    std::deque<Point3DDubins> SampleDubinsPathTo(const Point3DDubins &other, double dist);
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    Point3DDubins RotatePoint(double deltaYaw, double deltaPitch);
    Point3DDubins GetInvertedPoint();
    void PrintPosition(std::ostream &out);
  protected:
    double coords[3];
    double yaw;
    double pitch;
};

class Point3DPolynom {
  public:
    inline static double AverageVelocity;
    inline static double Gravity;

    Point3DPolynom();
    Point3DPolynom(double x, double y, double z);
    Point3DPolynom(const std::string &s, double scale=1);
    Point3DPolynom(const Vec3 pos);

    void SetPosition(double x, double y, double z);
    void SetPosition(Vec3 vec);
    void SetVelocity(double x, double y, double z);
    void SetVelocity(Vec3 vec);
    void SetAcceleration(double x, double y, double z);
    void SetAcceleration(Vec3 vec);
    const double* GetPosition() const;
    const double* GetRawCoords() const;
    const double operator[](int i) const;
    void operator+=(const Vector &translate);
    friend bool operator==(const Point3DPolynom &p1, const Point3DPolynom &p2);
    friend bool operator!=(const Point3DPolynom &p1, const Point3DPolynom &p2);
    friend bool operator<(const Point3DPolynom &p1, const Point3DPolynom &p2);
    friend Point3DPolynom operator/(const Point3DPolynom &p1, const double scale);    // scale position (NOT the rotation)
    double Distance(const Point3DPolynom &other) const;
    double EuclideanDistance(const Point3DPolynom &other) const;
    Point3DPolynom GetStateInDistance(Point3DPolynom &other, double dist) const;
    std::deque<Point3DPolynom> SampleTrajectory(Point3DPolynom &other, double interval);
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    void PrintPosition(std::ostream &out);

    Point3D GetPositionOnly();
  protected:
    double coords[3];
    double velocity[3];
    double acceleration[3];
};

class PointVector3D : public Vector {
  public:
    PointVector3D();
    PointVector3D(double x, double y, double z);
    PointVector3D(Point3D p);
    PointVector3D(Point3DDubins p);
    PointVector3D(Point3DPolynom p);
    PointVector3D(Point3D p1, Point3D p2);
    PointVector3D(Point3DDubins p1, Point3DDubins p2);
    PointVector3D(Point3DPolynom p1, Point3DPolynom p2);
};

class PointVector2D : public Vector {
  public:
    PointVector2D();
    PointVector2D(double x, double y);
    PointVector2D(Point2D p);
    PointVector2D(Point2DDubins p);
    PointVector2D(Point2DPolynom p);
    PointVector2D(Point2D p1, Point2D p2);
    PointVector2D(Point2DDubins p1, Point2DDubins p2);
    PointVector2D(Point2DPolynom p1, Point2DPolynom p2);
    PointVector3D To3DVector() const;
};

// STREAM OUTPUTS
std::ostream& operator<<(std::ostream &out, const Point2D &p);
std::ostream& operator<<(std::ostream &out, const Point2DDubins &p);
std::ostream& operator<<(std::ostream &out, const Point2DPolynom &p);
std::ostream& operator<<(std::ostream &out, const Point3D &p);
std::ostream& operator<<(std::ostream &out, const Point3DDubins &p);
std::ostream& operator<<(std::ostream &out, const Point3DPolynom &p);

#endif
