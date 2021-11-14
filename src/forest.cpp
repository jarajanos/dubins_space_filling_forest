/**
 * @file forest.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "forest.h"

/************************************************** 2D Dubins specialization *******************************************************/

template<>
void SpaceForest<Point2DDubins>::initNodeTypeSpecific(Node<Point2DDubins> &node) {
  // important because of checkOtherTrees - need to have correctly set up expanded angles
  node.ExpandedAngles = std::deque<int>(this->problem.DubinsResolution);
  for (int i{0}; i < this->problem.DubinsResolution; ++i) {
    node.ExpandedAngles.push_back(i);
  } 
}

template<>
bool SpaceForest<Point2DDubins>::pathToGoalExists(Point2DDubins &newPoint) {
  bool solved{false};
  for (int j{0}; j < this->problem.DubinsResolution && !solved; ++j) {
    Point2DDubins temp{this->problem.Goal};
    temp.SetAngle(j, this->problem.DubinsResolution);
    solved |= this->isPathFree(newPoint, temp);
  }

  return solved;
}

template<>
bool SpaceForest<Point2DDubins>::expandNode(Node<Point2DDubins> *expanded, bool &solved, const unsigned int iteration) {
  Point2DDubins newPoint;
  Node<Point2DDubins> *newNode; 
  std::deque<int> expandedAngles;
  std::deque<double> distances;

  Tree<Point2DDubins> *expandedTree{expanded->SourceTree};

  // sample new point and check basic collisions with obstacles
  if (getAndCheckNewPoint(expanded, &newPoint)) {
    return true;
  }
  
  flann::Matrix<float> pointToAdd{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION}; // just one point to add
  for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
    pointToAdd[0][i] = newPoint[i];
  }

  // ensure that the node does not "grow" into the tree
  if (checkExpandedTree(expanded, &newPoint, pointToAdd)) {
    delete[] pointToAdd.ptr();
    return true;
  }

  // finally check clear path to expanded node
  if (!expanded->IsRoot()) {
    if (!this->isPathFree(expanded->Position, newPoint)) {
      delete[] pointToAdd.ptr();
      return true;
    }
  } else {
    // check all angles - root only
    bool oneSuitable{false};
    for (int i{0}; i < this->problem.DubinsResolution; ++i) {
      Point2DDubins temp{expanded->Position};
      temp.SetAngle(i, this->problem.DubinsResolution);
      if (this->isPathFree(temp, newPoint)) {
        oneSuitable = true;
        expandedAngles.push_back(i);
        distances.push_back(temp.Distance(newPoint));
      } else {
        distances.push_back(std::numeric_limits<double>::max());
      }
    }

    if (!oneSuitable) {
      delete[] pointToAdd.ptr();
      return true;
    }
  }

  // check other trees, link them together or check availability of the goal node
  if (checkOtherTrees(expanded, &newPoint, pointToAdd, solved)) {
    delete[] pointToAdd.ptr();
    return true;
  }

  double parentDistance{newPoint.Distance(expanded->Position)};
  if (this->problem.Optimize && !expanded->IsRoot()) {
    // sff star -- do not optimize root children (optimization not needed, would mess up the expandedAngles list)
    optimizeConnections(expanded, &newPoint, newNode, pointToAdd, iteration);
  } else {
    if (!expanded->IsRoot()) {
      newNode = &(expandedTree->Leaves.emplace_back(newPoint, expandedTree, expanded, parentDistance, iteration));
    } else {
      newNode = &(expandedTree->Leaves.emplace_back(newPoint, expandedTree, expanded, distances, expandedAngles, iteration));
    }
    expanded->Children.push_back(newNode);
  }
  this->allNodes.push_back(newNode);
  
  // finally add the new node to flann and frontiers
  if (this->problem.PriorityBias != 0) {
    for (Heap<Node<Point2DDubins>> &prior : expandedTree->Frontiers) {
      prior.Push(newNode); 
    }
  } else {
    frontier.push_back(newNode);
  }
  expandedTree->Flann.Index->addPoints(pointToAdd);
  expandedTree->Flann.PtrsToDel.push_back(pointToAdd.ptr());

  if (solved) {
    // // this means the problem has goal and a path has been found
    // int goalRootId{this->problem.GetNumRoots() - 1};
    // int expandedRootId{expandedTree->Root->ID};

    // // temporal solution - select root angle with shortest path
    // double bestDist{std::numeric_limits<double>::max()};
    // int selectedRootAngleId;
    // for (auto &angleId : expanded->GetExpandedAngles()) {
    //   double dist{expanded->DistanceToRoot(angleId)};
    //   if (dist < bestDist) {
    //     bestDist = dist;
    //     selectedRootAngleId = angleId;
    //   }
    // }
    
    // // again - select the shortest possible path
    // bestDist = std::numeric_limits<double>::max();
    // int selectedGoalAngle, selectedGoalAngleId;
    // for (int i{0}; i < this->problem.DubinsResolution; ++i) {
    //   Point2DDubins tempGoal{this->problem.Goal};
    //   tempGoal.SetAngle(i, this->problem.DubinsResolution);

    //   if (isPathFree(newNode->Position, tempGoal)) {
    //     double dist{newNode->Position.Distance(tempGoal)};
    //     if (dist < bestDist) {
    //       bestDist = dist;
    //       selectedGoalAngle = tempGoal.GetAngle();
    //       selectedGoalAngleId = i;
    //     }
    //   }
    // }
    
    // DistanceHolder<Point2DDubins> holder{newNode, goalNode};
    // this->borders.AddLink(holder, expandedRootId, goalRootId, selectedRootAngleId, selectedGoalAngleId, true);  // true = second angle is "inlet" and shall not be inverted
    this->borders(expandedTree->Root->ID, this->problem.GetNumRoots() - 1).emplace_back(newNode, goalNode);
  }

  return false;
}

//optimize connections - same as for 2D and 3D, overriden only rewire and node emplacement
template<> 
bool SpaceForest<Point2DDubins>::filterNode(Node<Point2DDubins> &neighbor) {
  return neighbor.IsRoot();
}

template<>
void SpaceForest<Point2DDubins>::emplaceNewNode(Node<Point2DDubins> *expanded, Point2DDubins* newPoint, Node<Point2DDubins>* &newNode, const int iteration) {
  newNode = &(expanded->SourceTree->Leaves.emplace_back(*newPoint, expanded->SourceTree, expanded, expanded->Distance(*newPoint), iteration));
}

template<>
void SpaceForest<Point2DDubins>::rewireNodes(Node<Point2DDubins> *newNode, Node<Point2DDubins> &neighbor, double newDistance) {
  std::deque<Node<Point2DDubins> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring SFF: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest[0] = newDistance;
  newNode->Children.push_back(&neighbor);
}

// check borders for checkConnectivity function = check, whether the problem was solved
// solved = there exists ANY path connecting all roots, BUT IT MIGHT NOT BE CONTINUOUS
template<>
bool SpaceForest<Point2DDubins>::checkBorders(int id1, int id2) {
  auto link{this->borders(id1, id2)}; 
  if (link.empty()) {
    return false;
  } else if (link[0].IsValid) {
    return true;
  } else {
    while (!link.empty()) {
      DistanceHolder<Point2DDubins> holder{link.front()};

      if (!holder.Node1->IsRoot() && !holder.Node2->IsRoot()) {
        // shortcut solution
        if (this->isPathFree(holder.Node1->Position, holder.Node2->Position.GetInvertedPoint())) {
          link[0].IsValid = true;
          return true;
        } 
      } else {
        for (auto angle1 : holder.Node1->GetExpandedAngles()) {
          for (auto angle2 : holder.Node2->GetExpandedAngles()) {
            Point2DDubins point1{holder.Node1->DubinsPosition(angle1, this->problem.DubinsResolution, false)};
            Point2DDubins point2{holder.Node2->DubinsPosition(angle2, this->problem.DubinsResolution, true)};
            if (this->isPathFree(point1, point2)) {
              link[0].IsValid = true;
              return true;
            } 
          }
        }
      }

      link.pop_front();
    }
  }
  return false;
}

template<>
void SpaceForest<Point2DDubins>::getPaths() {
  // due to optimizations, the distance might have changed -> update it
  int numRoots{this->problem.GetNumRoots()};
  int numAngles{this->problem.DubinsResolution};
  for (int i{0}; i < numRoots; ++i) {
    for (int j{0}; j < numRoots; ++j) {
      std::deque<DistanceHolder<Point2DDubins>> &borderPoints{this->borders(i, j)};
      for (DistanceHolder<Point2DDubins> holder : borderPoints) { // intentional copy
        bool checkAngles{true};
        if (!holder.Node1->IsRoot() && !holder.Node2->IsRoot()) {
          if (holder.IsValid || this->isPathFree(holder.Node1->Position, holder.Node2->Position.GetInvertedPoint())) {
            checkAngles = false;
          } else {
            continue;
          }
        } 

        for (int angle1 : holder.Node1->GetExpandedAngles()) {
          for (int angle2 : holder.Node2->GetExpandedAngles()) {
            Point2DDubins point1{holder.Node1->DubinsPosition(angle1, numAngles, false)};
            Point2DDubins point2{holder.Node2->DubinsPosition(angle2, numAngles, true)};
            if (!checkAngles || this->isPathFree(point1, point2)) {
              // update distance for angle pair angle1 & angle2, but the "angle2" part is flown in opposite direction -> 
              // save it for oposite expanded angle (which should be the same trajectory) (ensured by AddLink, with secondIsInlet==false)
              holder.UpdateDistance(angle1, angle2);
              this->neighboringMatrix.AddLink(holder, i, j, angle1, angle2, false);
            }
          }
        }
      }

     for (int k{0}; k < numAngles; ++k) {
        for (int l{0}; l < numAngles; ++l) {
          if (!this->neighboringMatrix.Exists(i, j, k, l)) {
            continue;
          }

          DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(i, j, k, l)};
          std::deque<Point2DDubins> &plan{holder.Plan};
          // one tree
          Node<Point2DDubins> *nodeToPush{holder.Node1};
          plan.push_front(nodeToPush->DubinsPosition(k, numAngles, false));
          while (!nodeToPush->IsRoot()) {
            nodeToPush = nodeToPush->Closest;
            plan.push_front(nodeToPush->DubinsPosition(k, numAngles, false));
          }

          // second tree
          nodeToPush = holder.Node2;
          plan.push_back(nodeToPush->DubinsPosition(l, numAngles, true));
          while(!nodeToPush->IsRoot()) {
            nodeToPush = nodeToPush->Closest;
            plan.push_back(nodeToPush->DubinsPosition(l, numAngles, true));
          }
        }
      }
    }
  }
}
