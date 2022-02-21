/**
 * @file problem.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 1.0
 * @date 18. 02. 2022
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "problem.h"

template<> 
Problem<Point3D>::Problem(const Problem<Point3DDubins> &problem) {
  this->Repetition = problem.Repetition;
  this->Dimension = D3;

  this->Solver = problem.Solver;
  this->Optimize = problem.Optimize;

  this->Env = Environment<Point3D>(problem.Env);

  for (auto &root : problem.Roots) {
    this->Roots.emplace_back(root);
  }
  
  this->Goal = Point3D(problem.Goal);
  this->HasGoal = problem.HasGoal;
  this->ComputeTSP = problem.ComputeTSP;

  this->AutoRange = problem.AutoRange;

  this->DistTree = problem.DistTree;
  this->SamplingDist = problem.SamplingDist;
  this->CollisionDist = problem.CollisionDist;
  this->DubinsRadius = problem.DubinsRadius;
  
  this->MaxIterations = problem.MaxIterations;
  this->PriorityBias = problem.PriorityBias;
  this->MaxMisses = problem.MaxMisses;
  this->DubinsResolution = problem.DubinsResolution;
  this->PrmConnections = problem.PrmConnections;
  
  this->SaveOpt = problem.SaveOpt;
  this->FileNames = problem.FileNames;
  this->SaveFreq = problem.SaveFreq;
  this->ID = problem.ID;

  this->TspSolver = problem.TspSolver;
  this->TspType = problem.TspType;
}
