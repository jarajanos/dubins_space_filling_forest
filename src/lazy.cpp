/**
 * @file lazy.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "lazy.h"

template<>
LazyTSP<Point2DDubins>::LazyTSP(Problem<Point2DDubins> &problem) : Solver<Point2DDubins>(problem) {
  // create adjacency matrix
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    Node<Point2DDubins> &node{rootNodes.emplace_back(this->problem.Roots[i], nullptr, nullptr, 0, 0)};
    this->allNodes.push_back(&node);
  }
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{0}; j < this->problem.GetNumRoots(); ++j) {
      if (this->neighboringMatrix(i, j, 0, 0).Exists()) {
        continue;
      }
      DistanceHolder<Point2DDubins> dist{&(rootNodes[i]), &(rootNodes[j]), rootNodes[i].Distance(rootNodes[j])};
      this->neighboringMatrix.AddLink(dist, i, j, 0, 0);
    }
  }
}

template<>
double LazyTSP<Point2DDubins>::updateEdge(std::tuple<int, int> selectedEdge, int iter) {
  int first, second;
  std::tie(first, second) = selectedEdge;
  DistanceHolder<Point2DDubins> &edge{this->neighboringMatrix(first, second, 0, 0)};
  if (edge.Plan.empty()) {
    runRRT(&edge, iter);
  } 
  return edge.Distance;
}

template<> 
void LazyTSP<Point2DDubins>::rewireNodes(Node<Point2DDubins> *newNode, Node<Point2DDubins> &neighbor, double newDistance) {
  std::deque<Node<Point2DDubins> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring LazyTSP: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest[0] = newDistance;
  newNode->Children.push_back(&neighbor);
}

template<>
void LazyTSP<Point2DDubins>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) {
  INFO("Saving paths");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    int numRoots{this->problem.GetNumRoots()};
    if (file.type == Obj) {
      ERROR("Not implemented yet");
      // fileStream << "o Paths\n";
      // for (int i{0}; i < this->allNodes.size(); ++i) {
      //   R temp{this->allNodes[i]->Position / problem.Env.ScaleFactor};
      //   fileStream << "v" << DELIMITER_OUT;
      //   temp.PrintPosition(fileStream);
      //   fileStream << "\n";
      // }
      
      // for (int i{0}; i < numRoots; ++i) {
      //   for (int j{i + 1}; j < numRoots; ++j) {
      //     DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
      //     if (holder.Node1 == NULL) {
      //       continue;
      //     }

      //     std::deque<Node<R> *> &plan{holder.Plan};
      //     for (int k{0}; k < plan.size() - 1; ++k) {
      //       fileStream << "l" << DELIMITER_OUT << plan[k]->ID + 1 << DELIMITER_OUT << plan[k+1]->ID + 1 << "\n";
      //     }
      //   }
      // }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
      
      for (auto &pair : selectedPaths) {
        int first, second;
        std::tie(first, second) = pair;
        DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(first, second, 0, 0)};

        std::deque<Point2DDubins> &plan{holder.Plan};
        for (int k{0}; k < plan.size() - 1; ++k) {
          opendubins::State finishDub{plan[k][0], plan[k][1], plan[k].GetHeading()};
          opendubins::State startDub{plan[k+1][0], plan[k+1][1], plan[k+1].GetHeading()};
          opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
          
          Point2DDubins lastPoint{startDub};
          Point2DDubins actPoint;
          double length{pathFromClosest.length};
          double parts{length / this->problem.CollisionDist};
          for (int index{1}; index < parts; ++index) {
            actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
            fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
            lastPoint = actPoint;
          }
          actPoint = Point2DDubins(finishDub);
          fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
        }
        fileStream << "\n";
      }
    } else {
      throw std::string("Unimplemented file type");
    }

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }  
}

template<>
void LazyTSP<Point2DDubins>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime, const std::deque<std::tuple<int,int>> &selectedEdges) {
  INFO("Saving parameters");
  std::ofstream fileStream{file.fileName.c_str(), std::ios_base::openmode::_S_app};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;
    ERROR(message.str());
    return;
  }

  if (fileStream.is_open()) {
    fileStream << this->problem.ID << CSV_DELIMITER;
    fileStream << this->problem.Repetition << CSV_DELIMITER;
    fileStream << iterations << CSV_DELIMITER;
    fileStream << (solved ? "solved" : "unsolved") << CSV_DELIMITER;
    fileStream << "[";
    int iter{0};
    for (auto &pair : selectedEdges) {
      ++iter;
      int first, second;
      std::tie(first, second) = pair;
      fileStream << first;
      if (iter != this->problem.GetNumRoots()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
    iter = 0;
    for (auto &pair : selectedEdges) {
      ++iter;
      int first, second;
      std::tie(first, second) = pair;
      fileStream << this->neighboringMatrix(first,second, 0, 0).Distance / this->problem.Env.ScaleFactor;
      if (iter != this->problem.GetNumRoots()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER;
    fileStream << elapsedTime.count() << "\n";
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    ERROR(message.str());
  }
}
