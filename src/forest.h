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

    SymmetricMatrix<std::deque<DistanceHolder<Node<R>>>> borders;

    void getPaths() override;
};

template<class R>
SpaceForest<R>::SpaceForest(Problem<R> &problem) : Solver<R>(problem), borders{this->problem.GetNumRoots()} {
   std::cout << "Hello, I am OpenDubins library!" << std::endl;
    opendubins::Dubins d = opendubins::Dubins(opendubins::State(0, 0, 0), opendubins::State(1, 1, 1), 1);
    std::cout << d << std::endl;

}

template<class R>
void SpaceForest<R>::getPaths() {

}


#endif