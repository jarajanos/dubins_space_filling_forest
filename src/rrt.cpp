/**
 * @file rrt.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "rrt.h"

template<>
void RapidExpTree<Point2DDubins>::rewireNodes(Node<Point2DDubins> *newNode, Node<Point2DDubins> &neighbor, double newDistance) {
  std::deque<Node<Point2DDubins> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring RRT: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest[0] = newDistance;
  newNode->Children.push_back(&neighbor);
}

template<>
void RapidExpTree<Point2DDubins>::getPaths() {
  for (DistanceHolder<Point2DDubins> &link : centralRoot->Links) {
    if (link.Distance > 1e100) {
      ERROR("Fatal error: max distance reached");
      exit(1);
    } 
    std::deque<Point2DDubins> &plan{link.Plan};

    // one tree
    Node<Point2DDubins> *nodeToPush{link.Node1};
    plan.push_front(nodeToPush->Position);
    while (!nodeToPush->IsRoot()) {
      nodeToPush = nodeToPush->Closest;
      plan.push_front(nodeToPush->Position);
    }

    // second tree
    nodeToPush = link.Node2;
    plan.push_back(nodeToPush->Position);
    while(!nodeToPush->IsRoot()) {
      nodeToPush = nodeToPush->Closest;
      plan.push_back(nodeToPush->Position);
    }

    link.UpdateDistance();

    this->neighboringMatrix.AddLink(link, link.Node1->SourceTree->Root->ID, link.Node2->SourceTree->Root->ID, 0, 0);
  }
}
