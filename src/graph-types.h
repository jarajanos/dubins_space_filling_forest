/**
 * @file graph-types.h
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __GRAPH_TYPES_H__
#define __GRAPH_TYPES_H__

#include "common.h"
#include "heap.h"
#include <flann/flann.hpp>
#include <deque>

template<class R> class Node;
template<class R> class FlannHolder;
template<class R> class Tree;

template<class R>
class Node {
  public:
    int ID;
    R Position;

    Node<R> *Root;
    Node<R> *Closest;
    std::deque<Node<R> *> Children;
    double DistanceToRoot;
    double DistanceToClosest;

    bool operator<(const Node<R> &l);
    
    unsigned int GetAge() const {
      return generation;
    }
  private:
    inline static int globID = 0;
    unsigned int generation;
};

template<>
class FlannHolder<Node<Point2D>> {
  public:
    flann::Index<flann::L2<float>> *Index;
    std::deque<float *> PtrsToDel;

  ~FlannHolder();
};

template<>
class FlannHolder<Node<Point2DDubins>> {
  public:
    flann::Index<flann::L2Dubins<float>> *Index;
    std::deque<float *> PtrsToDel;

  ~FlannHolder();
};

template<>
class FlannHolder<Node<Point3D>> {
  public:
    flann::Index<D6Distance<float>> *Index;
    std::deque<float *> PtrsToDel;

  ~FlannHolder();
};

template<class R>
class Tree {
  public:
    std::deque<R> Leaves;
    FlannHolder<R> Flann;
    R *Root;

    std::vector<Heap<R>> Frontiers;
    std::deque<DistanceHolder<R>> Links;
};

template<class R>
bool Node<R>::operator<(const Node<R> &l) {
  return this->ID < l.ID;
}

#endif
