/**
 * @file rrt.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-15
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __RRT_H__
#define __RRT_H__

#include "solver.h"

template<class R>
class RapidExpTree : public Solver<R> {
  public:
    RapidExpTree(Problem<R> &problem); 

    void Solve() override;
  private:
    int actNumTrees;
    Tree<R> *centralRoot{nullptr};

    Node<R> *goalNode;
    std::deque<Tree<R> *> treeFrontier;
    UnionFind<Tree<R> *> expandedTrees;

    void expandNode(Tree<R> *tree, bool &solved, const unsigned int iteration);
    bool getAndCheckNewPoint(Tree<R> *treeToExpand, R *newPoint, Node<R>* &parent);
    bool optimizeConnections(Tree<R> *treeToExpand, R *newPoint, Node<R>* &parent, Node<R>* &newNode, const unsigned iteration);
    bool checkOtherRewire(Tree<R> *treeToExpand, R *newPoint, Node<R>* &newNode);
    void rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance);
    
    void getPaths() override;
    void getConnectedTrees();
};

template<> void RapidExpTree<Point2DDubins>::rewireNodes(Node<Point2DDubins> *newNode, Node<Point2DDubins> &neighbor, double newDistance);
template<> void RapidExpTree<Point2DDubins>::getPaths();

template<class R>
RapidExpTree<R>::RapidExpTree(Problem<R> &problem) : Solver<R>(problem) {
  for (int j{0}; j < this->problem.Roots.size(); ++j) {
    Tree<R> &tree{this->trees.emplace_back()};
    Node<R> &node{tree.Leaves.emplace_back(this->problem.Roots[j], &tree, nullptr, 0, 0)};
    tree.Root = &node;
    this->allNodes.push_back(&node);
    flann::Matrix<float> rootMat{new float[1 * PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      rootMat[0][i] = node.Position[i];
    }
    tree.Flann.CreateIndex(rootMat);

    treeFrontier.push_back(&tree);
  }    
  actNumTrees = this->problem.Roots.size() - 1;

  // add goal, which is not expanded
  if (this->problem.HasGoal) {
    Tree<R> &tree{this->trees.emplace_back()};
    flann::Matrix<float> rootMat{new float[1 * PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
    goalNode = &(tree.Leaves.emplace_back(this->problem.Goal, &tree, nullptr, 0, 0));
    this->allNodes.push_back(goalNode);
    tree.Root = goalNode;
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      rootMat[0][i] = goalNode->Position[i];
    }

    tree.Flann.CreateIndex(rootMat);

    treeFrontier.push_back(&tree);
  }   

  expandedTrees = UnionFind(treeFrontier);
}

template<class R>
void RapidExpTree<R>::Solve() {
  StopWatch watch;
  if (SaveGoals <= this->problem.SaveOpt) {
    this->saveCities(this->problem.FileNames[SaveGoals]);
  }

  watch.Start();

  int iter{0};
  bool solved{false};
  while (!solved && iter != this->problem.MaxIterations) {
    ++iter;
    Tree<R> *treeToExpand{this->treeFrontier[this->rnd.RandomIntMinMax(0, actNumTrees)]};
    expandNode(treeToExpand, solved, iter);

    this->checkIterationSaves(iter);
  }
  
  watch.Stop();

  getConnectedTrees();
  if (SaveTree <= this->problem.SaveOpt) {
    this->saveTrees(this->problem.FileNames[SaveTree]);
  }

  if (this->centralRoot != nullptr) {
    this->getPaths();
    this->getAllPaths();
    if (SaveRoadmap <= this->problem.SaveOpt) {
      this->savePaths(this->problem.FileNames[SaveRoadmap]);
    }
  }
  
  this->getConnected();
  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.FileNames[SaveParams], iter, solved, watch.GetElapsed());
  }

  if (SaveTSPFile <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSPFile]);
  }
}

template<class R>
bool RapidExpTree<R>::getAndCheckNewPoint(Tree<R> *treeToExpand, R *newPoint, Node<R>* &parent) {
  R rndPoint;
  if (this->problem.PriorityBias > 0 && this->rnd.RandomProbability() <= this->problem.PriorityBias) {
    rndPoint = goalNode->Position; 
  } else {
    this->rnd.RandomPointInSpace(rndPoint);  
  }
  
  // find nearest neighbour
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;
  flann::Matrix<float> rndPointMat{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION}; // just one point to add
  for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
    rndPointMat[0][i] = rndPoint[i];
  }
  treeToExpand->Flann.Index->knnSearch(rndPointMat, indices, dists, 1, flann::SearchParams(FLANN_NUM_SEARCHES));
  delete[] rndPointMat.ptr();

  parent = &(treeToExpand->Leaves[indices[0][0]]);
  
  // get point in this direction, check for collisions
  *newPoint = parent->Position.GetStateInDistance(rndPoint, this->problem.SamplingDist);
  if (this->problem.Env.Collide(*newPoint) || !this->isPathFree(parent->Position, *newPoint)) {
    return true;    
  }

  return false;
}

template<class R>
bool RapidExpTree<R>::optimizeConnections(Tree<R> *treeToExpand, R *newPoint, Node<R>* &parent, Node<R>* &newNode, const unsigned iteration) {
  double bestDist{newPoint->Distance(parent->Position) + parent->DistanceToRoot()};
  double krrt{2 * M_E * log10(this->allNodes.size())};
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;

  flann::Matrix<float> newPointMat{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
  for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
    newPointMat[0][i] = (*newPoint)[i];
  }
  treeToExpand->Flann.Index->knnSearch(newPointMat, indices, dists, krrt, flann::SearchParams(FLANN_NUM_SEARCHES));

  delete newPointMat.ptr();

  std::vector<int> &indRow{indices[0]};
  for (int &ind : indRow) {
    Node<R> &neighbor{treeToExpand->Leaves[ind]};
    double neighDist{neighbor.Distance(*newPoint) + neighbor.DistanceToRoot()};
    if (neighDist < bestDist - SFF_TOLERANCE && this->isPathFree(neighbor.Position, *newPoint)) {
      bestDist = neighDist;
      parent = &neighbor;
    }
  }

  newNode = &(treeToExpand->Leaves.emplace_back(*newPoint, parent->SourceTree, parent, parent->Distance(*newPoint), iteration));
  parent->Children.push_back(newNode);

  for (int &ind : indRow) {
    Node<R> &neighbor{treeToExpand->Leaves[ind]}; // offset goal node
    double newPointDist{newPoint->Distance(neighbor.Position)};
    double proposedDist{bestDist + newPointDist};
    if (proposedDist < neighbor.DistanceToRoot() - SFF_TOLERANCE && this->isPathFree(*newPoint, neighbor.Position)) {
      // rewire
      rewireNodes(newNode, neighbor, proposedDist);
    }
  }

  return false;
}

template<class R>
void RapidExpTree<R>::rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) {
  std::deque<Node<R> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring RRT: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest = newDistance;
  newNode->Children.push_back(&neighbor);
}

template<class R>
bool RapidExpTree<R>::checkOtherRewire(Tree<R> *treeToExpand, R *newPoint, Node<R>* &newNode) {
  flann::Matrix<float> refPoint{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
  for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
    refPoint[0][i] = (*newPoint)[i];
  }
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;

  // check distance to other trees and rewire
  for (int i{0}; i < this->treeFrontier.size(); ++i) {
    Tree<R> *tree{this->treeFrontier[i]};
    if (tree->Root->ID == treeToExpand->Root->ID) {
      continue;
    }

    tree->Flann.Index->knnSearch(refPoint, indices, dists, 1, flann::SearchParams(FLANN_NUM_SEARCHES));
    Node<R> &neighbor{tree->Leaves[indices[0][0]]}; // just one-to-one connection
    double neighDist{neighbor.Distance(*newPoint)};

    if (neighDist < this->problem.DistTree && this->isPathFree(*newPoint, neighbor.Position)) {
      // create link
      treeToExpand->Links.emplace_back(newNode, &neighbor);

      // transfer nodes to tree with lower index
      Tree<R> *neighborExpanded{expandedTrees.Find(neighbor.SourceTree)};
      Tree<R> *to{*treeToExpand < *neighborExpanded ? treeToExpand : neighborExpanded};
      Tree<R> *from{*treeToExpand < *neighborExpanded ? neighborExpanded : treeToExpand};
      expandedTrees.Union(from, to);

      flann::Matrix<float> fromNodes{new float[from->Leaves.size() * PROBLEM_DIMENSION], from->Leaves.size(), PROBLEM_DIMENSION};
      int rewirePos{static_cast<int>(to->Leaves.size())};
      for (int i{0}; i < from->Leaves.size(); ++i) {
        for (int j{0}; j < 2; ++j) {
          fromNodes[i][j] = from->Leaves[i].Position[j];
        }
        to->Leaves.push_back(from->Leaves[i]);
      }

      // update new points
      for (int j{rewirePos}; j < to->Leaves.size(); ++j) {
        Node<R> &node{to->Leaves[j]};
        
        if (node.Closest != nullptr) {
          auto subIter{find(to->Leaves.begin(), to->Leaves.end(), *(node.Closest))};
          node.Closest = &(*subIter);
        }

        std::deque<Node<R> *> temp{node.Children};
        node.Children.clear();
        for(Node<R> *&child : temp) {
          auto subIter{find(to->Leaves.begin(), to->Leaves.end(), *child)};
          node.Children.push_back(&(*subIter));
        }
      }

      // updating "to"'s links -- pointer to the other node might have changed
      std::deque<DistanceHolder<R>> temp{to->Links};
      to->Links.clear();
      for (DistanceHolder<R> &link : temp) {
        auto iter{find(to->Leaves.begin(), to->Leaves.end(), *link.Node1)};
        auto iter2{find(to->Leaves.begin(), to->Leaves.end(), *link.Node2)};
        to->Links.emplace_back(&(*iter), &(*iter2));
      }

      // updating "from"'s links
      for(DistanceHolder<R> &link : from->Links) {
        auto iter{find(to->Leaves.begin(), to->Leaves.end(), *link.Node1)};
        auto iter2{find(to->Leaves.begin(), to->Leaves.end(), *link.Node2)};
        to->Links.emplace_back(&(*iter), &(*iter2));
      }

      to->Flann.Index->addPoints(fromNodes);
      to->Flann.PtrsToDel.push_back(fromNodes.ptr());

      // delete tree 
      auto iter{find(this->treeFrontier.begin(), this->treeFrontier.end(), from)};
      if (iter == this->treeFrontier.end()) {
        ERROR("Fatal error during tree merging (RRT)");
        exit(1);
      }
      this->treeFrontier.erase(iter);
      treeToExpand = to;
      
      --actNumTrees;
      --i;
    }
  }
  delete[] refPoint.ptr();
  return false;
}

template<class R>
void RapidExpTree<R>::expandNode(Tree<R> *treeToExpand, bool &solved, const unsigned int iteration) {
  R newPoint;
  Node<R> *parent, *newNode;

  if (getAndCheckNewPoint(treeToExpand, &newPoint, parent)) {
    return;
  }

  // rrt star
  if (this->problem.Optimize) {
    if (optimizeConnections(treeToExpand, &newPoint, parent, newNode, iteration)) {
      return;
    }
  } else {
    newNode = &(treeToExpand->Leaves.emplace_back(newPoint, parent->SourceTree, parent, this->problem.SamplingDist, iteration));
    parent->Children.push_back(newNode);
  }
  this->allNodes.push_back(newNode);

  // add the point to flann
  flann::Matrix<float> pointToAdd{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
  for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
    pointToAdd[0][i] = newPoint[i];
  }
  treeToExpand->Flann.Index->addPoints(pointToAdd);
  treeToExpand->Flann.PtrsToDel.push_back(pointToAdd.ptr());

  checkOtherRewire(treeToExpand, &newPoint, newNode);

  solved = (this->treeFrontier.size() == 1);
}

template<class R>
void RapidExpTree<R>::getPaths() {
  for (DistanceHolder<R> &link : centralRoot->Links) {
    if (link.Distance > 1e100) {
      ERROR("Fatal error: max distance reached");
      exit(1);
    } 
    std::deque<R> &plan{link.Plan};

    // one tree
    Node<R> *nodeToPush{link.Node1};
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

    this->neighboringMatrix(link.Node1->SourceTree->Root->ID, link.Node2->SourceTree->Root->ID) = link;
  }
}

template <class R>
void RapidExpTree<R>::getConnectedTrees() {
  int maxConn{0};
  int numRoots{this->problem.HasGoal ? actNumTrees + 1 : actNumTrees};
  for (int i{0}; i < numRoots + 1; ++i) {   // +1 because it is rather the maximum index
    int count{this->expandedTrees.GetCountOf(this->treeFrontier[i])};
    if (count > maxConn) {
      maxConn = count;
      this->centralRoot = this->treeFrontier[i];
      this->connectedTrees = this->expandedTrees.GetAllWithParent(this->treeFrontier[i]);
    }
  }
}

#endif 
