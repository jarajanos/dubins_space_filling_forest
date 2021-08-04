/**
 * @file random-generator.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __RANDOM_GENERATOR_H__
#define __RANDOM_GENERATOR_H__

#include <random>
#include <math.h>
#include <vector>
#include <chrono>

#include "common.h"

typedef std::mt19937_64 randomEngine;   // set used random engine

template <class R>
class RandomGenerator {
  public:
    RandomGenerator(const Range<double> samplingRange);

    bool RandomPointInDistance(const R& center, R& point, const double distance);
    void RandomPointInSpace(R& point);
    int RandomIntMinMax(const int min, const int max);
    double RandomProbability();

  private:
    Range<double> limits;
    double Distance;
    randomEngine rndEng;

    std::uniform_real_distribution<double> uniDistAngle;
    std::uniform_real_distribution<double> uniSpaceX;
    std::uniform_real_distribution<double> uniSpaceY;
    std::uniform_real_distribution<double> uniSpaceZ;
    std::uniform_real_distribution<double> uniProb;

    bool isInLimits(R& p);
};

template<class R>
RandomGenerator<R>::RandomGenerator(const Range<double> samplingRange) : limits{samplingRange} {
  // seed with actual time
  std::chrono::time_point<std::chrono::high_resolution_clock, std::chrono::nanoseconds> tSeed{std::chrono::high_resolution_clock::now()};
  std::uint_fast64_t uSeed{static_cast<uint_fast64_t>(tSeed.time_since_epoch().count())};

  rndEng = randomEngine(uSeed);

  uniDistAngle = std::uniform_real_distribution<double>(-M_PI, M_PI);
  uniSpaceX = std::uniform_real_distribution<double>(limits.minX, limits.maxX);
  uniSpaceY = std::uniform_real_distribution<double>(limits.minY, limits.maxY);
  uniSpaceZ = std::uniform_real_distribution<double>(limits.minZ, limits.maxZ);
  uniProb = std::uniform_real_distribution<double>(0, 1);
}

template<class R>
int RandomGenerator<R>::RandomIntMinMax(const int min, const int max) {
  std::uniform_int_distribution<int> uniInt(min, max);
  return uniInt(rndEng);
}

template <class R>
double RandomGenerator<R>::RandomProbability() {
  return uniProb(rndEng);
}

#endif /*__RANDOM_GENERATOR_H__*/