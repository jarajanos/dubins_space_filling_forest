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
