/**
 * @file environment.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __ENVIRONMENT_H__
#define __ENVIRONMENT_H__

#include <deque>
#include <vector>

#include "common.h"
#include "RAPID.H"

template<class R> class Environment;
template<class R> class Obstacle;

template <class R>
class Triangle {
 public:
  Triangle(R x, R y, R z) : vertices{x, y, z} {
  }

  inline R operator[](int i) {
    return vertices[i];
  }

 protected:
  R vertices[3];
};

template <class R>
class Environment {
  public:
    std::deque<Obstacle<R>> Obstacles;
    Obstacle<R> *Robot;
    Range<double> limits{__DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__};
    bool HasMap{true};
    T ScaleFactor;

    Environment() : Robot{nullptr} {
    }

    ~Environment() {
      delete Robot;
    }

    bool Collide(R position); // whether robot in such position collides with any known obstacle
    
    void processLimits(Range<double> &limits) {
      this->limits.minX = MIN(this->limits.minX, limits.minX);
      this->limits.maxX = MAX(this->limits.maxX, limits.maxX);
      this->limits.minY = MIN(this->limits.minY, limits.minY);
      this->limits.maxY = MAX(this->limits.maxY, limits.maxY);
      this->limits.minZ = MIN(this->limits.minZ, limits.minZ);
      this->limits.maxZ = MAX(this->limits.maxZ, limits.maxZ);
    }

};

template<class R>
class Obstacle {
  public:
    inline static std::string Delimiter = " ";
    inline static std::string NameDelimiter = "_";

    R Position;
    Obstacle() {}

    Obstacle(const std::string fileName, const bool isObj, const R position, const double scaleFactor);
    Obstacle(const std::string fileName, const bool isObj, const double scaleFactor);

    virtual ~Obstacle();
    void ParseOBJFile(const std::string fileName);
    void ParseMapFile(const std::string fileName);
  
    R getPosition();
    RAPID_model *getRapidModel();
    R &getRange();

    static bool Collide(Obstacle<R> &object1, R pos1, Obstacle<R> &object2, R pos2);
    static bool Collide(Obstacle<R> &object1, Obstacle<R> &object2);
    static bool Collide(Obstacle<R> &object, Obstacle<R> &robot, R robPos);

  protected:
    std::vector<R> facePoints;
    std::vector<Triangle<T>> faces;
    RAPID_model *rapidModel = NULL;
    Range<double> localRange{__DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__};
    double scale{1};
    inline static int rapidId = 0;
    inline static double eyeRotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; 

    virtual void addPoint(int objId, double coords[3]);
    virtual void addFacet(int objId, int offset, int faceInts[3]);
};

#endif
