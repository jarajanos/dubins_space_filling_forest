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
#include "opendubins/dubins3D.h"
#include <flann/flann.hpp>
#include <deque>

template<class R> class NodeBase;
template<class R> class FlannHolder;
template<class R> class Tree;

template<class R>
class NodeBase {
  public:
    int ID;
    R Position;
    Node<R> *Closest{nullptr};
    Tree<R> *SourceTree{nullptr};
    std::deque<Node<R> *> Children;
    virtual double DistanceToRoot() = 0;

    NodeBase(R position, Tree<R> *root, Node<R> *closest, unsigned int iteration) : Position{position}, SourceTree{root}, 
      Closest{closest}, generation{iteration} {
        ID = globID++;
        isRoot = !iteration;
      }

    virtual double Distance(const NodeBase<R> &b) {
      return this->Position.Distance(b.Position);
    }

    virtual double Distance(const R &point) {
      return this->Position.Distance(point);
    }

    friend bool operator==(const NodeBase<R> &r, const NodeBase<R> &l) {
      return r.ID == l.ID;
    }

    bool operator<(const NodeBase<R> &l);

    unsigned int GetAge() const {
      return generation;
    }

    bool IsRoot() const {
      return isRoot;
    }

  private:
    inline static int globID = 0;
    unsigned int generation;
    bool isRoot;
};

template<> bool NodeBase<Point2DDubins>::operator<(const NodeBase<Point2DDubins> &l);
template<> bool NodeBase<Point3DDubins>::operator<(const NodeBase<Point3DDubins> &l);

template<class R>
class Node<R, false> : public NodeBase<R> {
  using NodeBase<R>::NodeBase;
  public:
    double DistanceToClosest;

    Node(R position, Tree<R> *root, Node<R> *closest, double distanceToClosest, unsigned int iteration) : NodeBase<R>(position, root, closest, iteration),
      DistanceToClosest{distanceToClosest} {
      }

    double DistanceToRoot() override {
      Node<R> *previous{this->Closest};
      double distance{this->DistanceToClosest};
      while (previous != nullptr) {
        distance += previous->DistanceToClosest;
        previous = previous->Closest;
      }

      return distance;
    }
};

template<class R>
class PrmNode : public NodeBase<R> {
  using NodeBase<R>::NodeBase;
  public:
    std::map<PrmNode<R>*, double> VisibleNodes;

    PrmNode() {
    }

    double DistanceToRoot() {
      return distanceToRoot;
    }

    void SetDistanceToRoot(double distance) {
      distanceToRoot = distance;
    }
    
    void Reset() {
      this->SetDistanceToRoot(std::numeric_limits<double>::max());
      this->Closest = nullptr;
    }
  private:
    double distanceToRoot{std::numeric_limits<double>::max()};
};

template<class R>
class Node<R, true> : public NodeBase<R> {
  using NodeBase<R>::NodeBase;
  public:
    std::deque<double> DistanceToClosest;
    std::deque<int> ExpandedAngles;

    Node(R position, Tree<R> *root, Node<R> *closest, double distanceToClosest, unsigned int iteration) : NodeBase<R>(position, root, closest, iteration) {
      this->DistanceToClosest = std::deque<double>(1, distanceToClosest);
    }

    // special constructor for root's ancestors
    Node(R position, Tree<R> *root, Node<R> *closest, std::deque<double> &distanceToClosest, std::deque<int> &expandedAngles, unsigned int iteration) : NodeBase<R>(position, root, closest, iteration) {
      this->DistanceToClosest = std::move(distanceToClosest);
      isRootChild = true;
      this->ExpandedAngles = std::move(expandedAngles);
    }

    double DistanceToRoot() override {
      Node<R> *previous{this};
      double distance{0};
      while (previous != nullptr) {
        if (!previous->isRootChild) {
          distance += previous->DistanceToClosest[0];
        } else {
          double min{std::numeric_limits<double>::max()};
          for (auto &item : previous->DistanceToClosest) {
            if (item < min) {
              min = item;
            }
          }
          distance += min;
        }
        previous = previous->Closest;
      }

      return distance;
    }

    double DistanceToRoot(const int angleID) {
      Node<R> *previous{this};
      double distance{0};
      while (previous != nullptr) {
        if (!previous->isRootChild) {
          distance += previous->DistanceToClosest[0];
        } else {
          distance += previous->DistanceToClosest[angleID];
        }
        previous = previous->Closest;
      }

      return distance;
    }

    std::deque<int> &GetExpandedAngles() {
      if (!this->isRootChild && !this->IsRoot()) {
        return this->Closest->GetExpandedAngles();
      } else {
        return ExpandedAngles;
      }
    }

    R DubinsPosition(int angleId, int angleResolution, bool reverse) {
      R pos{this->Position};
      if (this->IsRoot()) {
        pos.SetHeading(angleId, angleResolution);
      } else if (reverse) {
        pos = pos.GetInvertedPoint();
      }

      return pos;
    }

  private:
    bool isRootChild{false};
};

template<>
class FlannHolder<Node<Point2D>> {
  public:
    flann::Index<flann::L2<float>> *Index{nullptr};
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<>
class FlannHolder<Node<Point2DDubins>> {
  public:
    flann::Index<flann::L2Dubins<float>> *Index{nullptr};
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<>
class FlannHolder<Node<Point3D>> {
  public:
    //flann::Index<D6Distance<float>> *Index;
    flann::Index<flann::L2_3D<float>> *Index{nullptr};
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<>
class FlannHolder<Node<Point3DDubins>> {
  public:
    //flann::Index<D6Distance<float>> *Index;
    flann::Index<flann::L2_3D<float>> *Index{nullptr};
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<>
class FlannHolder<Node<Point3DPolynom>> {
  public:
    //flann::Index<D9Distance<float>> *Index{nullptr};
    flann::Index<flann::L2_3D<float>> *Index{nullptr};
    std::deque<float *> PtrsToDel;

    ~FlannHolder();

    void CreateIndex(flann::Matrix<float> &matrix);
};

template<class R>
class Tree {
  public:
    std::deque<Node<R>> Leaves;
    FlannHolder<Node<R>> Flann;
    Node<R> *Root;

    std::vector<Heap<Node<R>>> Frontiers;
    std::deque<DistanceHolder<R>> Links;

    void AddFrontier(Node<R> *reference) {
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
