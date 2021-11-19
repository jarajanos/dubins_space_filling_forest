/**
 * @file dijkstra.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief Implementation of the Dijkstra algorithm, as proposed by Robert Penicka
 * @version 1.0
 * @date 23/04/2020
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#ifndef __DIJKSTRA_H__
#define __DIJKSTRA_H__

#include <set>
#include "common.h"
#include "heap.h"

#define DIJKSTRA_MAX std::numeric_limits<double>::max()

template <class R>
class Dijkstra {
 public:
  Dijkstra(){
	}

	std::deque<DistanceHolder<R>> findPath(int start, std::vector<int> &goals, 
			std::deque<PrmNode<R>> &connectedPoints);

 private:
	PrmNode<R> *expandBest(std::deque<PrmNode<R>> &connectedPoints);
	std::unique_ptr<Heap<PrmNode<R>>> heap;
};

template <class R>
std::deque<DistanceHolder<R>> Dijkstra<R>::findPath(int start, std::vector<int> &goals, 
			std::deque<PrmNode<R>> &connectedPoints) {
	// 1. initialize nodes and heap
	std::deque<DistanceHolder<R>> retVal;

  for (auto &node : connectedPoints) {
    node.Reset();
  }

  connectedPoints[start].SetDistanceToRoot(0);
  connectedPoints[start].Closest = (Node<R> *)&(connectedPoints[start]);
	heap = std::make_unique<Heap<PrmNode<R>>>(connectedPoints, &(connectedPoints[start]), true);

	std::set<int> unvisitedGoals(goals.begin(), goals.end());
	
	// 2. dijkstra
	while(true) {
		PrmNode<R> *best{expandBest(connectedPoints)};

		if (best == NULL || best->DistanceToRoot() == DIJKSTRA_MAX) {
			break;
		}

		// verify goals (goals are only reward points, so...)
		if (best->IsRoot()) {
			auto goalIter{unvisitedGoals.find(best->ID)};
			if (goalIter != unvisitedGoals.end()) {
				unvisitedGoals.erase(best->ID);
				if (unvisitedGoals.empty()) {
					break;
				}
			}
		}
	}

	// 3. create list of plans
	for (int goal : goals) {
		DistanceHolder<R> &plan{retVal.emplace_back((Node<R> *)&(connectedPoints[start]), (Node<R> *)&(connectedPoints[goal]))};

		PrmNode<R> *actNode{&(connectedPoints[goal])};
		if (actNode->Closest != nullptr) {
			plan.Distance = actNode->DistanceToRoot();
			plan.Plan.push_front(actNode->Position);

			PrmNode<R> *startNode{&(connectedPoints[start])};
			while (actNode != startNode) {
				actNode = (PrmNode<R> *)actNode->Closest;
				plan.Plan.push_front(actNode->Position);
			}
		} else {
			//invalidate
			plan.Node1 = nullptr;
			plan.Node2 = nullptr;
			plan.Distance = DIJKSTRA_MAX;
		}
	}

	return retVal;
}

template <class R>
PrmNode<R>* Dijkstra<R>::expandBest(std::deque<PrmNode<R>> &connectedPoints) {
	PrmNode<R> *actual{heap->Get()};
	heap->Pop();

	if (actual != NULL) {
		for (auto &nodeDistPair : actual->VisibleNodes) {
      auto [ node, dist ] = nodeDistPair;
			double newDist{actual->DistanceToRoot() + dist};

			if (newDist < node->DistanceToRoot()) {
				node->Closest = (Node<R> *)actual;
				node->SetDistanceToRoot(newDist);
        heap->BubbleUp(node);
			}
		}
	}

	return actual;
}

#endif
