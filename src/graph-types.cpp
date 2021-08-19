/**
 * @file graph-types.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "graph-types.h"

void FlannHolder<Node<Point2D>>::CreateIndex(flann::Matrix<float> &matrix) {
  Index = new flann::Index<flann::L2<float>>(matrix, flann::KDTreeIndexParams(FLANN_NUM_KD_TREES));
  Index->buildIndex();
  this->PtrsToDel.push_back(matrix.ptr());
}

FlannHolder<Node<Point2D>>::~FlannHolder() {
  if (Index != nullptr) {
    delete Index;
  }

  for (auto ptr : PtrsToDel) {
    if (ptr != nullptr) {
      delete[] ptr;
    }
  }
}

void FlannHolder<Node<Point2DDubins>>::CreateIndex(flann::Matrix<float> &matrix) {
  Index = new flann::Index<flann::L2Dubins<float>>(matrix, flann::KDTreeIndexParams(FLANN_NUM_KD_TREES));
  Index->buildIndex();
  this->PtrsToDel.push_back(matrix.ptr());
}

FlannHolder<Node<Point2DDubins>>::~FlannHolder() {
  if (Index != nullptr) {
    delete Index;
  }
  
  for (auto ptr : PtrsToDel) {
    if (ptr != nullptr) {
      delete[] ptr;
    }
  }
}

void FlannHolder<Node<Point3D>>::CreateIndex(flann::Matrix<float> &matrix) {
  Index = new flann::Index<D6Distance<float>>(matrix, flann::KDTreeIndexParams(FLANN_NUM_KD_TREES));
  Index->buildIndex();
  this->PtrsToDel.push_back(matrix.ptr());
}

FlannHolder<Node<Point3D>>::~FlannHolder() {
  if (Index != nullptr) {
    delete Index;
  }
  
  for (auto ptr : PtrsToDel) {
    if (ptr != nullptr) {
      delete[] ptr;
    }
  }
}

// does not make sense = only purpose is correct ordering in DistanceHolder - position of Node in argument MAKES sense here, the problem is not symmetric
template<>
bool NodeBase<Point2DDubins>::operator<(const NodeBase<Point2DDubins> &l) {
  return true;
}
