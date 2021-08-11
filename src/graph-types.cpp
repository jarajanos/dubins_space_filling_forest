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

FlannHolder<Node<Point2D>>::~FlannHolder() {
  for (auto ptr : PtrsToDel) {
    if (ptr != nullptr) {
      delete[] ptr;
    }
  }
}

FlannHolder<Node<Point2DDubins>>::~FlannHolder() {
  for (auto ptr : PtrsToDel) {
    if (ptr != nullptr) {
      delete[] ptr;
    }
  }
}

FlannHolder<Node<Point3D>>::~FlannHolder() {
  for (auto ptr : PtrsToDel) {
    if (ptr != nullptr) {
      delete[] ptr;
    }
  }
}

// does not make sense = only purpose is correct ordering in DistanceHolder - position of Node in argument MAKES sense here, the problem is not symmetric
template<>
bool Node<Point2DDubins>::operator<(const Node<Point2DDubins> &l) {
  return true;
}
