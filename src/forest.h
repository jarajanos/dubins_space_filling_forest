/**
 * @file forest.h
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-11
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __FOREST_H__
#define __FOREST_H__

#include "solver.h"
#include "opendubins/dubins.h"

template<class R>
class SpaceForest : public Solver<R> {
  public:
    SpaceForest(Problem<R> &problem);
    void Solve() override;
  
  protected:
    std::deque<Node<R>*> frontier;
    std::deque<Node<R>*> closed_list;
    Node<R> *goalNode{nullptr};

    DistanceMatrix<std::deque<DistanceHolder<R>>> borders;

    void initBorders();
    void initNodeTypeSpecific(Node<R> &node);

    bool expandNode(Node<R> *expanded, bool &solved, const unsigned int iteration);
    bool getAndCheckNewPoint(Node<R> *expanded, R* newPoint);
    bool checkExpandedTree(Node<R> *expanded, R* newPoint, flann::Matrix<float> &matrix);
    bool checkOtherTrees(Node<R> *expanded, R* newPoint, flann::Matrix<float> &matrix, bool &solved);

    void optimizeConnections(Node<R> *expanded, R* newPoint, Node<R>* &newNode, flann::Matrix<float> &matrix, int iteration);
    bool filterNode(Node<R> &neighbor);
    void emplaceNewNode(Node<R> *expanded, R* newPoint, Node<R>* &newNode, int iteration);
    void rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance);

    int checkConnected();
    bool checkBorders(int id1, int id2);

    void getPaths() override;
    void saveFrontiers(const FileStruct file);
    void checkIterationSaves(const int iter) override;
};

template<> void SpaceForest<Point2DDubins>::initBorders();
template<> void SpaceForest<Point2DDubins>::initNodeTypeSpecific(Node<Point2DDubins> &node);
template<> bool SpaceForest<Point2DDubins>::expandNode(Node<Point2DDubins> *expanded, bool &solved, const unsigned int iteration);
template<> bool SpaceForest<Point2DDubins>::checkOtherTrees(Node<Point2DDubins> *expanded, Point2DDubins* newPoint, flann::Matrix<float> &matrix, bool &solved);
template<> bool SpaceForest<Point2DDubins>::filterNode(Node<Point2DDubins> &neighbor);
template<> void SpaceForest<Point2DDubins>::emplaceNewNode(Node<Point2DDubins> *expanded, Point2DDubins* newPoint, Node<Point2DDubins>* &newNode, int iteration);
template<> void SpaceForest<Point2DDubins>::rewireNodes(Node<Point2DDubins> *newNode, Node<Point2DDubins> &neighbor, double newDistance);
template<> bool SpaceForest<Point2DDubins>::checkBorders(int id1, int id2);
template<> void SpaceForest<Point2DDubins>::getPaths();

template<class R>
SpaceForest<R>::SpaceForest(Problem<R> &problem) : Solver<R>(problem) {
  initBorders();
  for (int j{0}; j < this->problem.Roots.size(); ++j) {
    Tree<R> &tree{this->trees.emplace_back()};
    Node<R> &node{tree.Leaves.emplace_back(this->problem.Roots[j], &tree, nullptr, 0, 0)};
    this->allNodes.push_back(&node);
    tree.Root = &node;
    flann::Matrix<float> rootMat{new float[1 * PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
    frontier.push_back(&node);
    
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      rootMat[0][i] = node.Position[i];
    }
    tree.Flann.CreateIndex(rootMat);
  }

  if (this->problem.PriorityBias != 0 && !this->problem.HasGoal) {
    for (int i{0}; i < this->trees.size(); ++i) {
      Tree<R> &tree{this->trees[i]};
      for (int j{0}; j < this->trees.size(); ++j) {
        if (i == j) {
          continue;
        }
        tree.AddFrontier(this->trees[j].Root);
      }
    }
  }

  // add goal - it is not expanded
  if (this->problem.HasGoal) {
    Tree<R> &tree{this->trees.emplace_back()};
    Node<R> &node{tree.Leaves.emplace_back(this->problem.Goal, &tree, nullptr, 0, 0)};
    this->allNodes.push_back(&node);
    flann::Matrix<float> rootMat{new float[1 * PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      rootMat[0][i] = this->problem.Goal[i];
    }
    tree.Flann.CreateIndex(rootMat);
    tree.Root = &node;
    goalNode = &node;

    if (this->problem.PriorityBias != 0) {  // must be done after goalNode assignment
      for (int i{0}; i < this->trees.size() - 1; ++i) {
        this->trees[i].AddFrontier(goalNode);
      }
    }
  }
}

template<class R>
void SpaceForest<R>::Solve() {
  Node<R> *nodeToExpand;
  StopWatch timeMeasure;

  if (SaveGoals <= this->problem.SaveOpt) {
    this->saveCities(this->problem.FileNames[SaveGoals]);
  }

  timeMeasure.Start();

  int iter{0};
  bool solved{false};
  bool emptyFrontier{false};
  // iterate
  while (!solved && iter < this->problem.MaxIterations) {
    Tree<R> *rndTree{nullptr};
    Heap<Node<R>> *prior{nullptr};
    int priorPos{-1};
    // get a random priority queue if applicable (priority bias is set and the open list is not empty)
    if (this->problem.PriorityBias != 0 && !emptyFrontier) {
      while (rndTree == nullptr || rndTree->EmptyFrontiers()) {
        rndTree = &(this->trees[this->rnd.RandomIntMinMax(0, this->trees.size() - 1)]);
      }
      while(priorPos == -1 || rndTree->Frontiers[priorPos].Empty()) {
        priorPos = this->rnd.RandomIntMinMax(0, rndTree->Frontiers.size() - 1);
      }
      prior = &(rndTree->Frontiers[priorPos]);
    }
    
    int pos;
    bool fromClosed{false};
    // get a random node to expand
    if (!closed_list.empty() && emptyFrontier) {
      // open list is empty, but the closed list is NOT empty
      pos = this->rnd.RandomIntMinMax(0, closed_list.size() - 1);
      nodeToExpand = closed_list[pos];
      fromClosed = true;
    } else {
      // get random node from open list, i. e. selected priority queue (when bias is given) or central list
      if (this->problem.PriorityBias != 0 && this->rnd.RandomProbability() <= this->problem.PriorityBias) {
        nodeToExpand = prior->Get();
      } else if (this->problem.PriorityBias != 0) {
        pos = this->rnd.RandomIntMinMax(0, prior->Size() - 1);
        nodeToExpand = prior->Get(pos);
      } else {
        pos = this->rnd.RandomIntMinMax(0, frontier.size() - 1);
        nodeToExpand = frontier[pos];
      }
    }

    // try to expand the selected node, given maximum number of tries and maximum number of all expansion-tries
    bool expandResult{true};
    for (int i{0}; i < this->problem.MaxMisses && expandResult && iter < this->problem.MaxIterations; ++i) {
      ++iter;
      expandResult &= expandNode(nodeToExpand, solved, iter);
      checkIterationSaves(iter);
    }
    if (expandResult && !fromClosed) {
      // expansion was NOT successfull -> delete the expanded node from the open list and insert it into closed list
      if (this->problem.PriorityBias == 0) {
        // case without priority bias
        auto iter{frontier.begin() + pos};
        frontier.erase(iter);
      } else {
        // priority bias -> delete from all priority queues
        for (int i{(int)(rndTree->Frontiers.size()) - 1}; i > -1; --i) {
          Heap<Node<R>> &p{rndTree->Frontiers[i]};
          for (int j{p.Size() - 1}; j > -1; --j) {
            if (p.Get(j) == nodeToExpand) {
              p.Pop(j);
              break;
            }
          }
        }
      }
      closed_list.push_back(nodeToExpand);
    }

    if (!solved) {
      // check, whether all frontiers are empty
      if (this->problem.PriorityBias != 0) {  // uses priority bias
        emptyFrontier = true;
        for (int i{0}; i < this->trees.size() && emptyFrontier; ++i) {
          emptyFrontier &= this->trees[i].EmptyFrontiers();
        }
      } else {
        emptyFrontier = frontier.empty();
      }

      bool connected{false};
      if (emptyFrontier) {
        // check whether all trees are connected
        connected = checkConnected() == this->problem.GetNumRoots();
      }

      solved = (!this->problem.HasGoal && emptyFrontier && connected);
    } else if (solved) {
      // update the connected trees to correct params output
      checkConnected(); 
    }
  }

  // consider as solved also the case, where frontiers are not empty, but all trees are connected
  if (!solved && !this->problem.HasGoal) {
    solved = checkConnected() == this->problem.GetNumRoots();
  }

  timeMeasure.Stop();

  if (SaveTree <= this->problem.SaveOpt) {
    this->saveTrees(this->problem.FileNames[SaveTree]);
  }
  
  getPaths();
  this->getAllPaths();

  if (this->problem.ComputeTSP) {
    this->computeTsp();
  }

  if (SaveRoadmap <= this->problem.SaveOpt) {
    this->savePaths(this->problem.FileNames[SaveRoadmap]);
  }

  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.FileNames[SaveParams], iter, solved, timeMeasure.GetElapsed());
  }

  if (SaveTSPFile <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSPFile]);
  }

  if (this->problem.ComputeTSP && SaveTSPPaths <= this->problem.SaveOpt) {
    this->saveTspPaths(this->problem.FileNames[SaveTSPPaths]);
  }

  if (SaveFrontiers <= this->problem.SaveOpt) {
    this->saveFrontiers(this->problem.FileNames[SaveFrontiers]);
  }
}

template<class R>
bool SpaceForest<R>::expandNode(Node<R> *expanded, bool &solved, const unsigned int iteration) {
  R newPoint;
  Node<R> *newNode; 

  Tree<R> *expandedTree{expanded->SourceTree};

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
  if (!this->isPathFree(expanded->Position, newPoint)) {
    delete[] pointToAdd.ptr();
    return true;
  }

  // check other trees, link them together or check availability of the goal node
  if (checkOtherTrees(expanded, &newPoint, pointToAdd, solved)) {
    delete[] pointToAdd.ptr();
    return true;
  }

  double parentDistance{expanded->Position.Distance(newPoint)};
  if (this->problem.Optimize) {
    // sff star
    optimizeConnections(expanded, &newPoint, newNode, pointToAdd, iteration);
  } else {
    newNode = &(expandedTree->Leaves.emplace_back(newPoint, expandedTree, expanded, parentDistance, iteration));
    expanded->Children.push_back(newNode);
  }
  this->allNodes.push_back(newNode);
  
  // finally add the new node to flann and frontiers
  if (this->problem.PriorityBias != 0) {
    for (Heap<Node<R>> &prior : expandedTree->Frontiers) {
      prior.Push(newNode); 
    }
  } else {
    frontier.push_back(newNode);
  }
  expandedTree->Flann.Index->addPoints(pointToAdd);
  expandedTree->Flann.PtrsToDel.push_back(pointToAdd.ptr());

  if (solved) {
    this->borders(this->problem.GetNumRoots() - 1, expandedTree->Root->ID).emplace_back(newNode, goalNode);
  }

  return false;
}

template<class R>
bool SpaceForest<R>::getAndCheckNewPoint(Node<R> *expanded, R* newPoint) {
  // sample new point
  if (!this->rnd.RandomPointInDistance(expanded->Position, *newPoint, this->problem.SamplingDist)) {
    return true;
  }

  // check limits and collisions
  if (this->problem.Env.Collide(*newPoint)) {     
    return true;
  }

  return false;
}

template<class R>
bool SpaceForest<R>::checkExpandedTree(Node<R> *expanded, R* newPoint, flann::Matrix<float> &matrix) {
  double parentDistance{expanded->Position.Distance(*newPoint)};
  Tree<R> *expandedTree{expanded->SourceTree};
  int expandedRootID{expandedTree->Root->ID};

  //find nearest neighbors - perform radius search
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;
  double checkDist{FLANN_PREC_MULTIPLIER * this->problem.SamplingDist};
  int neighbours{expandedTree->Flann.Index->radiusSearch(matrix, indices, dists, SQR(checkDist), 
    flann::SearchParams(FLANN_NUM_SEARCHES))}; // better and faster than knnSearch - tested by perf
  
  std::vector<int> &indRow{indices[0]}; // just one point
  for (int i{0}; i < neighbours; ++i) {
    int neighID{indRow[i]};
    Node<R> *neighbour{&(expandedTree->Leaves[neighID])};
    
    double realDist{neighbour->Position.Distance(*newPoint)};
    
    if (realDist < (parentDistance - SFF_TOLERANCE) && this->isPathFree(neighbour->Position, *newPoint)) {
      // closer node available in the same tree
      return true;
    }
  }

  return false;
}

template<class R>
bool SpaceForest<R>::checkOtherTrees(Node<R> *expanded, R* newPoint, flann::Matrix<float> &matrix, bool &solved) {
  int expandedRootID{expanded->SourceTree->Root->ID};
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;

  double checkDist{FLANN_PREC_MULTIPLIER * this->problem.DistTree};
  for (int j{0}; j < this->trees.size(); ++j) {
    Tree<R> &tree{this->trees[j]};
    int neighbourRootID{tree.Root->ID};
    if (neighbourRootID == expandedRootID) {
      continue;
    } 

    int neighbours{tree.Flann.Index->radiusSearch(matrix, indices, dists, SQR(checkDist), 
      flann::SearchParams(FLANN_NUM_SEARCHES))};
    
    std::vector<int> &indRow{indices[0]}; // just one point
    for (int i{0}; i < neighbours; ++i) {
      int neighID{indRow[i]};
      Node<R> *neighbour{&(tree.Leaves[neighID])};
      
      double realDist{neighbour->Position.Distance(*newPoint)};
  
      if (realDist < (this->problem.DistTree - SFF_TOLERANCE)) {   // if realDist is bigger than tree distance, it might get expanded later
        // neighbouring trees, add to neighboring matrix and below dTree distance -> invalid point, but save expanded
        // check whether the goal was achieved
        if (this->problem.HasGoal && neighbour->Position == this->problem.Goal) {
          solved = this->isPathFree(*newPoint, this->problem.Goal);  // when true, break the solve loop and after the creation of the new node, add it to the neighbouring matrix
        } else if (!this->problem.HasGoal) {   
          std::deque<DistanceHolder<R>> &borderPoints{this->borders(neighbourRootID, expandedRootID)};
          DistanceHolder<R> holder{neighbour, expanded};
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

template<class R>
void SpaceForest<R>::optimizeConnections(Node<R> *expanded, R* newPoint, Node<R>* &newNode, flann::Matrix<float> &matrix, int iteration) {
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;
  double bestDist{newPoint->Distance(expanded->Position) + expanded->DistanceToRoot()};
  double ksff{2 * M_E * log10(this->allNodes.size())};

  expanded->SourceTree->Flann.Index->knnSearch(matrix, indices, dists, ksff, flann::SearchParams(FLANN_NUM_SEARCHES));

  std::vector<int> &indRow{indices[0]};
  for (int &ind : indRow) {
    Node<R> &neighbor{expanded->SourceTree->Leaves[ind]};
    if (filterNode(neighbor)) {   // additional check of possible parent -- mainly becaused of Dubins (forbid connection with root)
      continue;
    }
    double neighborDist{neighbor.Position.Distance(*newPoint) + neighbor.DistanceToRoot()};
    if (neighborDist < bestDist - SFF_TOLERANCE && this->isPathFree(neighbor.Position, *newPoint)) {
      bestDist = neighborDist;
      expanded = &neighbor;
    }
  }

  emplaceNewNode(expanded, newPoint, newNode, iteration);
  expanded->Children.push_back(newNode);

  for (int &ind : indRow) {
    Node<R> &neighbor{expanded->SourceTree->Leaves[ind]};
    double newPointDist{newPoint->Distance(neighbor.Position)};
    double proposedDist{bestDist + newPointDist};
    if (proposedDist < neighbor.DistanceToRoot() - SFF_TOLERANCE && this->isPathFree(*newPoint, neighbor.Position)) {
      // rewire
      rewireNodes(newNode, neighbor, proposedDist);
    }
  }
}

template<class R> 
bool SpaceForest<R>::filterNode(Node<R> &neighbor) {
  return false;
}

template<class R>
void SpaceForest<R>::emplaceNewNode(Node<R> *expanded, R* newPoint, Node<R>* &newNode, const int iteration) {
  newNode = &(expanded->SourceTree->Leaves.emplace_back(*newPoint, expanded->SourceTree, expanded, newPoint->Distance(expanded->Position), iteration));
}

template<class R>
void SpaceForest<R>::rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) {
  std::deque<Node<R> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring SFF: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest = newDistance;
  newNode->Children.push_back(&neighbor);
}

template<class R>
int SpaceForest<R>::checkConnected() {
  int maxConn{0};
  this->connectedTrees.clear();
  int remaining{this->problem.GetNumRoots()};
  bool connMat[this->problem.GetNumRoots()];
  for (int i{1}; i < remaining; ++i) {
    connMat[i] = false;
  }
  int unconnected{0};
  
  while (maxConn < remaining) {
    this->connectedTrees.clear();
    std::deque<int> stack{unconnected};
    connMat[unconnected] = true;

    while(!stack.empty()) {
      int root{stack.front()};
      stack.pop_front();
      this->connectedTrees.push_back(&(this->trees[root]));
      for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
        if (root == i) {
          continue;
        }

        if (checkBorders(root, i) && !connMat[i]) {
          connMat[i] = true;
          stack.push_front(i);
        }
      }
    }
    maxConn = this->connectedTrees.size();
    for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
      if (!connMat[i]) {
        unconnected = i;
        break;
      }
    }
    remaining -= maxConn;
  }

  return maxConn;
}

template<class R>
bool SpaceForest<R>::checkBorders(int id1, int id2) {
  auto link{this->borders(id1, id2)}; 
  if (link.empty()) {
    return false;
  } else if (link[0].IsValid) {
    return true;
  } else {
    while (!link.empty()) {
      DistanceHolder<R> holder{link[0]};
      if (this->isPathFree(holder.Node1->Position, holder.Node2->Position)) {
        link[0].IsValid = true;
        return true;
      } 
      link.pop_front();
    }

    return false;
  }
}

template<class R>
void SpaceForest<R>::getPaths() {
  // due to optimizations, the distance might have changed -> update it
  int numRoots{this->problem.GetNumRoots()};
  for (int i{0}; i < numRoots; ++i) {
    for (int j{i + 1}; j < numRoots; ++j) {
      if (this->borders(i, j).empty()) {
        continue;
      }

      std::deque<DistanceHolder<R>> &borderPoints{this->borders(i, j)};
      for (DistanceHolder<R> &dist : borderPoints) {
        dist.UpdateDistance();
      }
      std::sort(borderPoints.begin(), borderPoints.end());
    }
  }

  // go through borders, select one with the shortest distance, then update neighboring matrix and create paths
  for (int i{0}; i < numRoots; ++i) {
    for (int j{i + 1}; j < numRoots; ++j) {
      if (this->borders(i, j).empty()) {
        continue;
      }

      std::deque<DistanceHolder<R>> &borderPoints{this->borders(i, j)};
      for (DistanceHolder<R> &dist : borderPoints) {
        if (dist.IsValid || this->isPathFree(dist.Node1->Position, dist.Node2->Position)) {
          this->neighboringMatrix(i, j) = dist;
          break;
        }
      }

      if (!this->neighboringMatrix.Exists(i,j)) {
        continue;
      }

      DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
      std::deque<R> &plan{holder.Plan};
      // one tree
      Node<R> *nodeToPush{holder.Node1};
      plan.push_front(nodeToPush->Position);
      while (!nodeToPush->IsRoot()) {
        nodeToPush = nodeToPush->Closest;
        plan.push_front(nodeToPush->Position);
      }

      // second tree
      nodeToPush = holder.Node2;
      plan.push_back(nodeToPush->Position);
      while(!nodeToPush->IsRoot()) {
        nodeToPush = nodeToPush->Closest;
        plan.push_back(nodeToPush->Position);
      }
    }
  }
}

template<class R>
void SpaceForest<R>::saveFrontiers(const FileStruct file) {
  INFO("Saving frontiers");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
  }

  if (fileStream.is_open()) {
    int numRoots{(int)this->connectedTrees.size()};
    if (file.type == Obj) {
      fileStream << "o Open nodes\n";
      if (this->problem.PriorityBias != 0) {
        for (int i{0}; i < this->trees.size(); ++i) {
          Heap<Node<R>> &heap{this->trees[i].Frontiers[0]}; // all points are in all heaps, printing only one of them is sufficient
          for (HeapNode<Node<R>> &heapNode : heap.GetHeapVector()) {
            R temp{heapNode.Node->Position / this->problem.Env.ScaleFactor};
            fileStream << "v" << DELIMITER_OUT;
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }
        }
      } else {
        for (Node<R> *node : this->frontier) {
          fileStream << "v" << DELIMITER_OUT << node->Position / this->problem.Env.ScaleFactor << "\n";
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Frontiers" << DELIMITER_OUT << this->problem.Dimension << "\n";
      if (this->problem.PriorityBias != 0) {
        for (int i{0}; i < this->trees.size(); ++i) {
          Heap<Node<R>> &heap{this->trees[i].Frontiers[0]};   // all nodes are in all heaps, so printing out the first heap is sufficient
          for (HeapNode<Node<R>> &heapNode : heap.GetHeapVector()) {
            fileStream << heapNode.Node->Position / this->problem.Env.ScaleFactor << DELIMITER_OUT << "1\n"; // the one in the end is just for plotting purposes
          }
        }
      } else {
        for (Node<R> *node : this->frontier) {
          fileStream << node->Position / this->problem.Env.ScaleFactor << DELIMITER_OUT << "1\n";
        }
      }
    } else {
      throw std::string("Not implemented");
    }

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }
}

template<class R>
void SpaceForest<R>::checkIterationSaves(const int iter) {
  Solver<R>::checkIterationSaves(iter);

  if (this->problem.SaveFreq[SaveFrontiers] != 0 && iter % this->problem.SaveFreq[SaveFrontiers] == 0) {
    std::string prefix{"iter_" + std::to_string(iter) + "_"};
    this->saveFrontiers(PrefixFileName(this->problem.FileNames[SaveFrontiers], prefix));
  }
}

template<class R>
void SpaceForest<R>::initBorders() {
  this->borders = DistanceMatrix<std::deque<DistanceHolder<R>>>(this->problem.GetNumRoots());
}

template<class R>
void SpaceForest<R>::initNodeTypeSpecific(Node<R> &node) {
}

#endif