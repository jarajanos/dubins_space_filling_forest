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

#define DEFAULT_THRES_MISS 3
#define DEFAULT_OBSTAC_MISS 10
#define TOLERANCE 1e-9

#define SAMPLING_ANGLE_DISP 1 // something about 60 degrees
#define DEFAULT_SAMP_DIST 1
#define EXPANSION_MULTIPLIER 2
#define DEFAULT_DIST_DIV 1
#define ANGLE_MOVE 0.8

const std::string WHITESPACE = " \n\r\t\f\v";

#endif 
