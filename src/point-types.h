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

#include <regex>
#include <math.h>
#include "constants.h"
#include "vector-types.h"
#include "common.h"

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
  private:
    double coords[2];
};

class Point2DDubins {
  public:
    Point2DDubins();
    Point2DDubins(double x, double y);
    Point2DDubins(const std::string &s, double scale=1);

    Point2DDubins& operator=(const Point2DDubins &p);
    void SetAngle(double angle);
    double GetAngle() const;
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
    void FillRotationMatrix(double (&matrix)[3][3]) const;
    void PrintPosition(std::ostream &out);
  private:
    double coords[2];
    double phi;
};

class Point3D {
  public:
    Point3D();
    Point3D(double x, double y, double z, double yaw, double pitch, double roll);
    Point3D(const std::string &s, double scale=1);
    
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
  private:
    double coords[3];
    Quaternion rotation;
};

class PointVector3D : public Vector {
  public:
    PointVector3D();
    PointVector3D(double x, double y, double z);
    PointVector3D(Point3D p);
    PointVector3D(Point3D p1, Point3D p2);
};

class PointVector2D : public Vector {
  public:
    PointVector2D();
    PointVector2D(double x, double y);
    PointVector2D(Point2D p);
    PointVector2D(Point2DDubins p);
    PointVector2D(Point2D p1, Point2D p2);
    PointVector2D(Point2DDubins p1, Point2DDubins p2);
    PointVector3D To3DVector() const;
};

#endif
