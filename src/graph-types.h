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
#include "opendubins/dubins.h"
#include <flann/flann.hpp>
#include <deque>

template<class R> class NodeBase;
template<class R> class Node;
template<class R> class FlannHolder;
template<class R> class Tree;

template<class R>
class NodeBase {
  public:
    int ID;
    R Position;

    Tree<Node<R>> *Root;
    Node<R> *Closest;
    std::deque<Node<R> *> Children;
    double DistanceToRoot;
    double DistanceToClosest;

    NodeBase(R position, Tree<Node<R>> *root, Node<R> *closest, double distanceToClosest, double distanceToRoot, unsigned int iteration) : Position{position},
      Root{root}, Closest{closest}, DistanceToClosest{distanceToClosest}, DistanceToRoot{distanceToRoot}, generation{iteration} {
        ID = globID++;
      }

    friend bool operator==(const NodeBase<R> &r, const NodeBase<R> &l) {
      return r.ID == l.ID;
    }

    bool operator<(const NodeBase<R> &l);

    unsigned int GetAge() const {
      return generation;
    }

    bool IsRoot() const {
      return this->Root->Root->ID == this->ID;
    }
  private:
    inline static int globID = 0;
    unsigned int generation;
};


template<class R>
class Node : public NodeBase<R> {
  using NodeBase<R>::NodeBase;
};

template<>
class Node<Point2DDubins> : public NodeBase<Point2DDubins> {
  using NodeBase<Point2DDubins>::NodeBase;
  public:
    // opendubins::Dubins PathToClosest;
    // opendubins::Dubins PathFromClosest;
};

template<>
class FlannHolder<Node<Point2D>> {
  public:
    flann::Index<flann::L2<float>> *Index;
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<>
class FlannHolder<Node<Point2DDubins>> {
  public:
    flann::Index<flann::L2Dubins<float>> *Index;
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<>
class FlannHolder<Node<Point3D>> {
  public:
    flann::Index<D6Distance<float>> *Index;
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<class R>
class Tree {
  public:
    std::deque<R> Leaves;
    FlannHolder<R> Flann;
    R *Root;

    std::vector<Heap<R>> Frontiers;
    std::deque<DistanceHolder<R>> Links;

    void AddFrontier(R *reference) {
      Frontiers.emplace_back(this->Leaves, reference);
    }

    const bool EmptyFrontiers() {
      bool retVal{true};
      for (auto &prior : Frontiers) {
        retVal &= prior.Empty();
      }
      
      return retVal;
    }

    friend bool operator<(const Tree<R> &a, const Tree<R> &b) {
      return a.Root->ID < b.Root->ID;
    }
};

template<class R>
bool NodeBase<R>::operator<(const NodeBase<R> &l) {
  return this->ID < l.ID;
}

#endif
