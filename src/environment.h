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
#include <iostream>
#include <fstream>

#include "common.h"
#include "RAPID.H"

template<class R> class Environment;
template<class R> class Obstacle;

template <class R>
class Triangle {
 public:
  Triangle(R x, R y, R z) : vertices{x, y, z} {
  }

  Triangle(const Triangle<Point3DDubins> &t) {
    ERROR("Unsupported conversion");
    exit(1);
  }

  inline R operator[](int i) const {
    return vertices[i];
  }

 protected:
  R vertices[3];
};

template<> Triangle<Point3D>::Triangle(const Triangle<Point3DDubins> &t);
template<> Triangle<Point3DDubins>::Triangle(const Triangle<Point3DDubins> &t);

template <class R>
class Environment {
  public:
    std::deque<Obstacle<R>> Obstacles;
    std::shared_ptr<Obstacle<R>> Robot;
    Range Limits{__DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__, __DBL_MAX__};
    bool HasMap{true};
    double ScaleFactor;

    Environment() {
    }

    Environment(const Environment<Point3DDubins> &env) {
      ERROR("Unsupported conversion of environment");
      exit(1);
    }

    bool Collide(R position); // whether robot in such position collides with any known obstacle
    
    void ProcessLimits(Range &limits) {
      for (int i{0}; i < 3; ++i) {
        this->Limits.mins[i] = MIN(this->Limits.mins[i], limits.mins[i]);
        this->Limits.maxs[i] = MAX(this->Limits.maxs[i], limits.maxs[i]);
      }
    }

};

template<> Environment<Point3D>::Environment(const Environment<Point3DDubins> &env);

template<class R>
class Obstacle {
  public:
    inline static std::string Delimiter = " ";
    inline static std::string NameDelimiter = "_";

    R Position;
    Obstacle() {}
    Obstacle(const Obstacle<Point3DDubins> &obst) {
      ERROR("Unsupported conversion of obstacle");
      exit(1);
    }

    Obstacle(const std::string fileName, const FileType type, const R position, const double scaleFactor);
    Obstacle(const std::string fileName, const FileType type, const double scaleFactor);

    //virtual ~Obstacle();
    void ParseOBJFile(const std::string fileName);
    void ParseMapFile(const std::string fileName);
  
    R GetPosition();
    std::shared_ptr<RAPID_model> GetRapidModel() const;
    Range GetRange() const;
    double GetScale() const;

    std::vector<R> GetFacePoints() const;
    std::vector<Triangle<R>> GetFaces() const;

    static bool Collide(Obstacle<R> &object1, R pos1, Obstacle<R> &object2, R pos2);
    static bool Collide(Obstacle<R> &object1, Obstacle<R> &object2);
    static bool Collide(Obstacle<R> &object, Obstacle<R> &robot, R robPos);

  protected:
    std::vector<R> facePoints;
    std::vector<Triangle<R>> faces;
    std::shared_ptr<RAPID_model> rapidModel;
    Range localRange{__DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__, __DBL_MAX__, -__DBL_MAX__, __DBL_MAX__};
    double scale{1};
    inline static int rapidId = 0;
    inline static double eyeRotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; 

    virtual void addPoint(int objId, double coords[3]);
    virtual void addFacet(int objId, int offset, int faceInts[3]);
    virtual void updateLocalRange(double coords[]);
};

template<> Obstacle<Point3D>::Obstacle(const Obstacle<Point3DDubins> &obst);

template <class R>
Obstacle<R>::Obstacle(const std::string fileName, const FileType type, const double scaleFactor) : Obstacle<R>(fileName, type, R(), scaleFactor) {
}

template <class R>
Obstacle<R>::Obstacle(const std::string fileName, const FileType type, const R position, const double scaleFactor) : Position{position}, scale{scaleFactor} {
  this->rapidModel = std::make_shared<RAPID_model>();
  this->rapidModel->BeginModel();

  if (type == Obj) {
    // OBJ file
    ParseOBJFile(fileName);
  } else {
    // other map file
    ParseMapFile(fileName);
  }  

  this->rapidModel->EndModel();
}

// template <class R>
// Obstacle<R>::~Obstacle() {
//   if (rapidModel != NULL) {
//     delete rapidModel;
//   }
// }

template <class R>
void Obstacle<R>::ParseOBJFile(const std::string fileName) {
  std::ifstream modelFile(fileName);
  unsigned int elem{0}, offset{0};
  int faceCache[3], objId{0};
  double pointCache[3];
  std::string line, value;

  while (getline(modelFile, line)) {
    ParseString(line, value, line, this->Delimiter);
    switch (value[0]) {
      case 'v':
        for (int i{0}; i < 3; ++i) {
          double d;
          ParseString(line, value, line, this->Delimiter);
          d = std::stod(value);
          pointCache[i] = d + Position[i];
        }
        addPoint(objId, pointCache);
        ++elem;
        break;
      case 'f':
        for (int i{0}; i < 3; ++i) {
          int k;
          ParseString(line, value, line, this->Delimiter);
          k = std::stoi(value);
          faceCache[i] = k;
        }

        addFacet(objId, offset, faceCache);
        break;
      case 'o':
        if (objId != 0) {
          offset += elem;
          elem = 0;
        }

        ParseString(line, value, line, this->NameDelimiter);
        break;
    }
  }
  modelFile.close();
}

template<class R>
void Obstacle<R>::updateLocalRange(double coords[]) {
  for (int i{0}; i < 3; ++i) {
    localRange.mins[i] = MIN(localRange.mins[i], coords[i]);
    localRange.maxs[i] = MAX(localRange.maxs[i], coords[i]);
  }
}

template <class R>
R Obstacle<R>::GetPosition() {
  return this->position;
}

template <class R>
std::shared_ptr<RAPID_model> Obstacle<R>::GetRapidModel() const {
  return this->rapidModel;
}

template <class R>
Range Obstacle<R>::GetRange() const {
  return this->localRange;
}

template <class R>
double Obstacle<R>::GetScale() const {
  return this->scale;
}

template <class R>
std::vector<R> Obstacle<R>::GetFacePoints() const {
  return this->facePoints;
}

template <class R>
std::vector<Triangle<R>> Obstacle<R>::GetFaces() const {
  return this->faces;
}

template <class R>
bool Environment<R>::Collide(R position) {
  if (!this->HasMap) {
    return false;
  } 
  
  bool retVal{false};
  for (Obstacle<R> &obs : this->Obstacles) {
    retVal |= Obstacle<R>::Collide(obs, *(this->Robot.get()), position);
  }
  return retVal;
}

/**
 * @brief Handler of the main collide function, assuming null positions of the objects
 * 
 * BEWARE THAT THE DEFAULT POSITION OF OBJECT IS ALREADY CONSIDERED AND SHOULD NOT BE
 * ADDITIONALLY PASSED TO RAPID
 */
template <class R>
bool Obstacle<R>::Collide(Obstacle<R> &object1, Obstacle<R> &object2) {
  R nullPnt;
  return Collide(object1, nullPnt, object2, nullPnt);
}

#endif
