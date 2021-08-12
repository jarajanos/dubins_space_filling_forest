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

    DistanceMatrix<std::deque<DistanceHolder<Node<R>>>> borders;

    bool expandNode(Node<R> *expanded, bool &solved, const unsigned int iteration);
    int checkConnected();

    void getPaths() override;
    void saveFrontiers(const FileStruct file);
    void checkIterationSaves(const int iter) override;
};

template<class R>
SpaceForest<R>::SpaceForest(Problem<R> &problem) : Solver<R>(problem), borders{this->problem.GetNumRoots()} {
  for (int j{0}; j < this->problem.Roots.size(); ++j) {
    Tree<Node<R>> &tree{this->trees.emplace_back()};
    Node<R> &node{tree.Leaves.emplace_back(this->problem.Roots[j], &tree, nullptr, 0, 0, 0)};
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
      Tree<Node<R>> &tree{this->trees[i]};
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
    Tree<Node<R>> &tree{this->trees.emplace_back()};
    Node<R> &node{tree.Leaves.emplace_back(this->problem.Goal, &tree, nullptr, 0, 0, 0)};
    this->allNodes.push_back(&node);
    flann::Matrix<float> rootMat{new float[1 * PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      rootMat[0][i] = this->problem.Goal[i];
    }
    tree.Flann.CreateIndex(rootMat);
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
    Tree<Node<R>> *rndTree{nullptr};
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
      // check whether all trees are connected
      bool connected{checkConnected() == this->problem.GetNumRoots()};
      
      // check, whether all frontiers are empty
      if (this->problem.PriorityBias != 0) {  // uses priority bias
        emptyFrontier = true;
        for (int i{0}; i < this->trees.size() && emptyFrontier; ++i) {
          emptyFrontier &= this->trees[i].EmptyFrontiers();
        }
        solved = (!this->problem.HasGoal && emptyFrontier && connected);
      } else {
        emptyFrontier = frontier.empty();
      }
    } else {
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

  if (SaveRoadmap <= this->problem.SaveOpt) {
    this->savePaths(this->problem.FileNames[SaveRoadmap]);
  }

  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.FileNames[SaveParams], iter, solved, timeMeasure.GetElapsed());
  }

  if (SaveTSP <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSP]);
  }

  if (SaveFrontiers <= this->problem.SaveOpt) {
    this->saveFrontiers(this->problem.FileNames[SaveFrontiers]);
  }
}

template<class R>
bool SpaceForest<R>::expandNode(Node<R> *expanded, bool &solved, const unsigned int iteration) {

  return true;
}

template<class R>
int SpaceForest<R>::checkConnected() {
  int maxConn{0};
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

        // TODO: non-emptinness does NOT assure correct connectivity due to lazy behaviour
        if (!this->borders(root, i).empty() && !connMat[i]) {
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
void SpaceForest<R>::getPaths() {

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

#endif