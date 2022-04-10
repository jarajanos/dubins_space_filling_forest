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

template<>
bool SpaceForestBase<Point3DPolynom>::checkExpandedTree(Node<Point3DPolynom> *expanded, Point3DPolynom* newPoint, flann::Matrix<float> &matrix) {
  double parentDistance{expanded->Position.EuclideanDistance(*newPoint)};
  Tree<Point3DPolynom> *expandedTree{expanded->SourceTree};
  int expandedRootID{expandedTree->Root->ID};
  Point3DPolynom newPosition{newPoint->GetPosition()[0], newPoint->GetPosition()[1], newPoint->GetPosition()[2]};

  //find nearest neighbors - perform radius search
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;
  double checkDist{FLANN_PREC_MULTIPLIER * this->problem.SamplingDist};
  int neighbours{expandedTree->Flann.Index->radiusSearch(matrix, indices, dists, SQR(checkDist), 
    flann::SearchParams(FLANN_NUM_SEARCHES))}; // better and faster than knnSearch - tested by perf
  
  std::vector<int> &indRow{indices[0]}; // just one point
  for (int i{0}; i < neighbours; ++i) {
    int neighID{indRow[i]};
    Node<Point3DPolynom> *neighbor{&(expandedTree->Leaves[neighID])};
    
    // distance betweeen positions only -- checking EUCLIDEAN DISTANCE not the trajectory cost
    double realDist{neighbor->Position.EuclideanDistance(*newPoint)};
    Point3DPolynom neighPoint{neighbor->Position};
    Point3DPolynom neighPosition{neighPoint.GetPosition()[0], neighPoint.GetPosition()[1], neighPoint.GetPosition()[2]};
    
    if (realDist < (parentDistance - SFF_TOLERANCE)) {
      // local planner checking ONLY for position (not trajectory between points)
      double parts{realDist / problem.CollisionDist};
      bool isFree{true};
      PointVector3D direction{newPosition, neighPosition};
      Point3DPolynom position;
      for (int index{1}; index < parts && isFree; ++index) {
        position.SetPosition(newPosition[0] + direction[0] * index / parts, newPosition[1] + direction[1] * index / parts, newPosition[2] + direction[2] * index / parts);

        isFree &= !problem.Env.Collide(position);
      }
      
      if (isFree) {
        return true;
      }
    }


  }

  return false;
}

template<> 
bool SpaceForestBase<Point3DPolynom>::checkOtherTrees(Node<Point3DPolynom> *expanded, Point3DPolynom* newPoint, flann::Matrix<float> &matrix, bool &solved) {
  int expandedRootID{expanded->SourceTree->Root->ID};
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;

  double checkDist{FLANN_PREC_MULTIPLIER * this->problem.DistTree};
  for (int j{0}; j < this->trees.size(); ++j) {
    Tree<Point3DPolynom> &tree{this->trees[j]};
    int neighbourRootID{tree.Root->ID};
    if (neighbourRootID == expandedRootID) {
      continue;
    } 

    int neighbours{tree.Flann.Index->radiusSearch(matrix, indices, dists, SQR(checkDist), 
      flann::SearchParams(FLANN_NUM_SEARCHES))};
    
    std::vector<int> &indRow{indices[0]}; // just one point
    for (int i{0}; i < neighbours; ++i) {
      int neighID{indRow[i]};
      Node<Point3DPolynom> *neighbour{&(tree.Leaves[neighID])};
      
      // checking EUCLIDEAN distance to the tree (more like a tolerance)
      double realDist{neighbour->Position.EuclideanDistance(*newPoint)};
  
      if (realDist < (this->problem.DistTree - SFF_TOLERANCE)) {   // if realDist is bigger than tree distance, it might get expanded later
        // neighbouring trees, add to neighboring matrix and below dTree distance -> invalid point, but save expanded
        // check whether the goal was achieved
        if (this->problem.HasGoal && neighbour->Position == this->problem.Goal) {
          solved = pathToGoalExists(*newPoint);  // when true, break the solve loop and after the creation of the new node, add it to the neighbouring matrix
        } else if (!this->problem.HasGoal) {   
          std::deque<DistanceHolder<Point3DPolynom>> &borderPoints{this->borders(neighbourRootID, expandedRootID)};
          DistanceHolder<Point3DPolynom> holder{neighbour, expanded};
          if (std::find(borderPoints.begin(), borderPoints.end(), holder) == borderPoints.end()) {  // is it necessary? maybe use sorted set?
            this->borders(neighbourRootID, expandedRootID).push_back(holder);
          }
        }

        if (!solved) {  // if solved, the new node lies near the goal and should be created
          return true;
        }
      }
    }
  }

  return false;
}

