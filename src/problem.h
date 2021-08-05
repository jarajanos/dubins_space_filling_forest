/**
 * @file problem.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PROBLEM_H__
#define __PROBLEM_H__

#include <deque>

#include "common.h"

template<class R>
struct Problem {
  int iteration{0};
  Dimensions dimension{D3};

  SolverType solver;
  bool optimal;
  bool smoothing;

  Environment<R> environment;
  std::deque<R> roots;
  R goal;
  bool hasGoal{false};
  
  bool autoRange{false};

  double distTree;
  double collisionDist;
  
  int maxIterations;
  double priorityBias{0};
  int saveTreeIter{0};
  int saveFrontiersIter{0};
  
  SaveOptions saveOptions{None};
  std::map<SaveOptions, FileStruct> fileNames;
  std::string id{ "Solver" };

  std::string tspSolver;
  std::string tspType;

  int GetNumRoots() {
    if (hasGoal) {
      return roots.size() + 1;
    } else {
      return roots.size();
    }
  }
};

#endif /* __PROBLEM_H__ */
