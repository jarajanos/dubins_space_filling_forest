/**
 * @file heap.h
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __HEAP_H__
#define __HEAP_H__

#include <deque>
#include "common.h"

#define LCHILD(x) 2 * (x) + 1
#define RCHILD(x) 2 * (x) + 2
#define PARENT(x) ((x) - 1) / 2

template<class R>
struct HeapNode {
  R *Node;
  double distance;

  HeapNode(R* node, R *refPoint) : Node{node} {
    distance = Node->Distance(*refPoint);
  }

  double Distance() {
    return distance;
  }

  friend bool operator==(const HeapNode<R> &a, const HeapNode<R> &b) {
    return a.Node == b.Node;
  }
};

template<class R>
struct HeapNode<PrmNode<R>> {
  PrmNode<R> *Node;

  HeapNode(PrmNode<R>* node, PrmNode<R> *refPoint) : Node{node} {
  }

  double Distance() {
    return Node->DistanceToRoot();
  }

  friend bool operator==(const HeapNode<PrmNode<R>> &a, const HeapNode<PrmNode<R>> &b) {
    return a.Node == b.Node;
  }
};

template<class R>
class Heap {
  public:
    Heap(std::deque<R> &leaves, R *refPoint, bool sort=false);
    
    R *Get();
    R *Get(int index);
    void Pop();
    void Pop(int index);
    void Push(R *item);
    void BubbleUp(R *item);

    std::deque<HeapNode<R>> &GetHeapVector();

    int Size();
    bool Empty();
  private:
    std::deque<HeapNode<R>> heapVector;
    R *goalNode;

    void bubbleDown(int index);
    void bubbleUp(int index);
    void Sort();
    double getCost(int Index);
};

template<class R>
Heap<R>::Heap(std::deque<R> &leaves, R *refPoint, bool sort) : goalNode{refPoint} {
  for (int i{0}; i < leaves.size(); ++i) {
    heapVector.emplace_back(&(leaves[i]), refPoint);
  }

  if (sort) {
    this->Sort();
  }
}

template<class R>
R* Heap<R>::Get() {
  if (Empty()) {
    return nullptr;
  }

  return heapVector[0].Node;
}

template<class R>  
R* Heap<R>::Get(int index) {
  if (index >= Size()) {
    return nullptr;
  }

  return heapVector[index].Node;
}

template<class R>
void Heap<R>::Pop() {
  if (heapVector.empty()) {
    return;
  }

  heapVector[0] = std::move(heapVector[Size() - 1]);
  heapVector.pop_back();
  bubbleDown(0);
}

template<class R>
void Heap<R>::Pop(int index) {
  int maxSize{Size()};

  if (index < maxSize - 1) {
    double oldCost{heapVector[index].Distance()};
    double newCost{heapVector[maxSize - 1].Distance()};

    std::swap(heapVector[index], heapVector[maxSize - 1]);
    heapVector.pop_back();
    if (oldCost > newCost) {
      bubbleUp(index);
    } else if (index <= Size() / 2) {
      bubbleDown(index);
    }
  } else if (index == maxSize - 1) {
    heapVector.pop_back();
  }
}

template<class R>
void Heap<R>::Push(R *item) {
  int index{Size()};
  heapVector.emplace_back(item, this->goalNode);
  bubbleUp(index);
}

template<class R>
void Heap<R>::BubbleUp(R *item) {
  int index;
  HeapNode<R> dummy{item, this->goalNode};
  auto iter{std::find(heapVector.begin(), heapVector.end(), dummy)};
  if (iter != heapVector.end()) {
    index = std::distance(heapVector.begin(), iter);
    this->bubbleUp(index);
  }
}

template<class R>
std::deque<HeapNode<R>> &Heap<R>::GetHeapVector() {
  return heapVector;
}

template<class R>
int Heap<R>::Size() {
  return static_cast<int>(heapVector.size());
}

template<class R>
bool Heap<R>::Empty() {
  return Size() == 0;
}

template<class R>
void Heap<R>::bubbleDown(int index) {
  int leftChildIndex{LCHILD(index)};
  int rightChildIndex{RCHILD(index)};
  int minIndex{index};
  int maxSize{Size()};

  if (leftChildIndex < maxSize && getCost(leftChildIndex) < getCost(minIndex)) {
    minIndex = leftChildIndex;
  }

  if (rightChildIndex < maxSize && getCost(rightChildIndex) < getCost(minIndex)) {
    minIndex = rightChildIndex;
  }

  if (minIndex != index) {
    // need to swap
    std::swap(heapVector[minIndex], heapVector[index]);
    bubbleDown(minIndex);
  }
}

template<class R>
void Heap<R>::bubbleUp(int index) {
  if (index == 0) {
    return;
  }

  int parentIndex = PARENT(index);
  if (getCost(parentIndex) > getCost(index)) {
    // swap
    std::swap(heapVector[parentIndex], heapVector[index]);
    bubbleUp(parentIndex);
  }
}

template<class R>
void Heap<R>::Sort() {
  for (int i{Size() / 2}; i >= 0; --i) {
    bubbleDown(i);
  }
}

template<class R>
double Heap<R>::getCost(int index) {
  return heapVector[index].Distance();
}

#endif
