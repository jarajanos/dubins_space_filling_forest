/**
 * @file environment.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "environment.h"

template <>
void Obstacle<Point2D>::addPoint(int objId, double coords[3]) {
  // scale
  for (int i{0}; i < 2; ++i) {
    coords[i] *= scale;
  }

  coords[2] = 0;
  this->facePoints.emplace_back(coords[0], coords[1]);
  updateLocalRange(coords);
}

template <>
void Obstacle<Point2DDubins>::addPoint(int objId, double coords[3]) {
  // scale
  for (int i{0}; i < 2; ++i) {
    coords[i] *= scale;
  }

  coords[2] = 0;
  this->facePoints.emplace_back(coords[0], coords[1]);
  updateLocalRange(coords);
}

template <>
void Obstacle<Point3D>::addPoint(int objId, double coords[3]) {
  // scale
  for (int i{0}; i < 3; ++i) {
    coords[i] *= scale;
  }

  this->facePoints.emplace_back(coords[0], coords[1], coords[2], 0, 0, 0);
  updateLocalRange(coords);
}

template <>
void Obstacle<Point2D>::ParseMapFile(const std::string fileName) {
  std::ifstream modelFile(fileName);
  int index{1};
  int faceCache[3];
  double pointCache[3] = {0, 0, 0};
  std::string line, value;

  while (getline(modelFile, line)) {
    line = Trim(line);
    if (!strcmp(line.c_str(), "")) {
      continue;
    }

    for (int i{0}; i < 3; ++i) {
      for (int j{0}; j < 2; ++j) {
        ParseString(line, value, line, this->Delimiter);
        pointCache[j] = std::stod(value) + Position[j];
      }
      addPoint(0, pointCache);
      faceCache[i] = index + i;
    }
    index += 3;

    addFacet(0, 0, faceCache);
  }
  modelFile.close();
}

template <>
void Obstacle<Point2DDubins>::ParseMapFile(const std::string fileName) {
  std::ifstream modelFile(fileName);
  int index{1};
  int faceCache[3];
  double pointCache[3] = {0, 0, 0};
  std::string line, value;

  while (getline(modelFile, line)) {
    line = Trim(line);
    if (!strcmp(line.c_str(), "")) {
      continue;
    }

    for (int i{0}; i < 3; ++i) {
      for (int j{0}; j < 2; ++j) {
        ParseString(line, value, line, this->Delimiter);
        pointCache[j] = std::stod(value) + Position[j];
      }
      addPoint(0, pointCache);
      faceCache[i] = index + i;
    }
    index += 3;

    addFacet(0, 0, faceCache);
  }
  modelFile.close();
}

template <>
void Obstacle<Point3D>::ParseMapFile(const std::string fileName) {
  std::ifstream modelFile(fileName);
  int index{1};
  int faceCache[3];
  double pointCache[3] = {0, 0, 0};
  std::string line, value;

  while (getline(modelFile, line)) {
    line = Trim(line);
    if (!strcmp(line.c_str(), "")) {
      continue;
    }

    for (int i{0}; i < 3; ++i) {
      for (int j{0}; j < 3; ++j) {
        ParseString(line, value, line, this->Delimiter);
        pointCache[j] = std::stod(value) + Position[j];
      }
      addPoint(0, pointCache);
      faceCache[i] = index + i;
    }
    index += 3;

    addFacet(0, 0, faceCache);
  }
  modelFile.close();
}

/**
 * @brief Handler for the RAPID functions for the collision checking
 * 
 * @param object1 First object to be tested for collisions
 * @param pos1 Position of the first object
 * @param object2 Second object to be tested
 * @param pos2 Position of the second object
 * @return true If objects collide
 * @return false Otherwise
 * 
 * BEWARE THAT THE DEFAULT POSITION OF OBJECT IS ALREADY CONSIDERED AND SHOULD NOT BE
 * ADDITIONALLY PASSED TO RAPID
 */

template <>
bool Obstacle<Point2D>::Collide(Obstacle<Point2D> &object1, Point2D pos1, Obstacle<Point2D> &object2, Point2D pos2) {  
  PointVector2D vecPos1{pos1}, vecPos2{pos2};
  double rotMat1[3][3], rotMat2[3][3];
  pos1.FillRotationMatrix(rotMat1);
  pos2.FillRotationMatrix(rotMat2);

  RAPID_Collide(rotMat1, (vecPos1.To3DVector()).GetRawCoords(), object1.GetRapidModel(), rotMat2, (vecPos2.To3DVector()).GetRawCoords(), object2.GetRapidModel());
  return RAPID_num_contacts != 0;
}

template <>
bool Obstacle<Point2DDubins>::Collide(Obstacle<Point2DDubins> &object1, Point2DDubins pos1, Obstacle<Point2DDubins> &object2, Point2DDubins pos2) {  
  PointVector2D vecPos1{pos1}, vecPos2{pos2};
  double rotMat1[3][3], rotMat2[3][3];
  pos1.FillRotationMatrix(rotMat1);
  pos2.FillRotationMatrix(rotMat2);

  RAPID_Collide(rotMat1, (vecPos1.To3DVector()).GetRawCoords(), object1.GetRapidModel(), rotMat2, (vecPos2.To3DVector()).GetRawCoords(), object2.GetRapidModel());
  return RAPID_num_contacts != 0;
}

template <>
bool Obstacle<Point3D>::Collide(Obstacle<Point3D> &object1, Point3D pos1, Obstacle<Point3D> &object2, Point3D pos2) {  
  PointVector3D vecPos1{pos1}, vecPos2{pos2};
  double rotMat1[3][3], rotMat2[3][3];
  pos1.FillRotationMatrix(rotMat1);
  pos2.FillRotationMatrix(rotMat2);

  RAPID_Collide(rotMat1, vecPos1.GetRawCoords(), object1.GetRapidModel(), rotMat2, vecPos2.GetRawCoords(), object2.GetRapidModel());
  return RAPID_num_contacts != 0;
}

/**
 * @brief Handler of the base Collide function, tailored for the robot model and its position
 * 
 * BEWARE THAT THE DEFAULT POSITION OF OBJECT IS ALREADY CONSIDERED AND SHOULD NOT BE
 * ADDITIONALLY PASSED TO RAPID
 */
template <>
bool Obstacle<Point2D>::Collide(Obstacle<Point2D> &object, Obstacle<Point2D> &robot, Point2D robPos) {
  PointVector2D vecPos1, vecPos2{robPos};
  double rotMat2[3][3];
  robPos.FillRotationMatrix(rotMat2);

  RAPID_Collide(Obstacle<Point2D>::eyeRotation, (vecPos1.To3DVector()).GetRawCoords(), object.GetRapidModel(), rotMat2, (vecPos2.To3DVector()).GetRawCoords(), robot.GetRapidModel());
  return RAPID_num_contacts != 0;
}

template <>
bool Obstacle<Point2DDubins>::Collide(Obstacle<Point2DDubins> &object, Obstacle<Point2DDubins> &robot, Point2DDubins robPos) {
  PointVector2D vecPos1, vecPos2{robPos};
  double rotMat2[3][3];
  robPos.FillRotationMatrix(rotMat2);

  RAPID_Collide(Obstacle<Point2DDubins>::eyeRotation, (vecPos1.To3DVector()).GetRawCoords(), object.GetRapidModel(), rotMat2, (vecPos2.To3DVector()).GetRawCoords(), robot.GetRapidModel());
  return RAPID_num_contacts != 0;
}

template <>
bool Obstacle<Point3D>::Collide(Obstacle<Point3D> &object, Obstacle<Point3D> &robot, Point3D robPos) {
  PointVector3D vecPos1, vecPos2{robPos};
  double rotMat2[3][3];
  robPos.FillRotationMatrix(rotMat2);

  RAPID_Collide(Obstacle<Point3D>::eyeRotation, vecPos1.GetRawCoords(), object.GetRapidModel(), rotMat2, vecPos2.GetRawCoords(), robot.GetRapidModel());
  return RAPID_num_contacts != 0;
}