#pragma once
/**
 * @file constants.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <string>

#ifndef DELIMITER_OUT
#define DELIMITER_OUT (" ")
#endif

#ifndef CSV_DELIMITER
#define CSV_DELIMITER (",")
#endif

#ifndef CSV_DELIMITER_2
#define CSV_DELIMITER_2 (";")
#endif

#ifndef TSP_DELIMITER
#define TSP_DELIMITER (" ")
#endif

// used for LazyTSP
#ifndef TEMP_TSP
#define TEMP_TSP      "tempTsp.tsp"
#endif

#ifndef TEMP_RESULT
#define TEMP_RESULT   "tempTsp.result"
#endif


#define DEFAULT_THRES_MISS 3
#define SFF_TOLERANCE 1e-9

#define DEFAULT_SAMP_DIST 1

#define FLANN_NUM_KD_TREES    4
#define FLANN_NUM_SEARCHES    128
#define FLANN_PREC_MULTIPLIER 1.1

const std::string WHITESPACE = " \n\r\t\f\v";

#endif 
