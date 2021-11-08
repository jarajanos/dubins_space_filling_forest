/**
 * @file common.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "common.h"


void DistanceHolder<Point2DDubins>::UpdateDistance(int angleId1, int angleId2) {
  Distance = Node1->Position.Distance(Node2->Position);
  if (Distance >= std::numeric_limits<double>::max()) {
    ERROR("Infinity distance between points");
    return;
  }

  if (angleId1 != -1) {
    Distance += Node1->DistanceToRoot(angleId1);
  } else {
    Distance += Node1->DistanceToRoot();
  }

  if (Distance >= std::numeric_limits<double>::max()) {
    ERROR("Infinity distance to root 1");
    return;
  }

  if (angleId2 != -1) {
    Distance += Node2->DistanceToRoot(angleId2);
  } else {
    Distance += Node2->DistanceToRoot();
  }

  if (Distance >= std::numeric_limits<double>::max()) {
    ERROR("Infinity distance to root 2");
    return;
  }
}

void StopWatch::Start() {
  this->startTime = std::chrono::high_resolution_clock::now();
}

void StopWatch::Stop() {
  this->stopTime = std::chrono::high_resolution_clock::now();
}

std::chrono::duration<double> StopWatch::GetElapsed() {
  return this->stopTime - this->startTime;
}

FileStruct PrefixFileName(const FileStruct &path, const std::string &insert) {
	FileStruct retVal{path};

	auto pos{retVal.fileName.find_last_of("//")};
  if (pos != std::string::npos) {
    retVal.fileName.insert(pos + 1, insert);
  } else {
    retVal.fileName.insert(0, insert);
  }

	return retVal;
}

std::string Ltrim(const std::string &s) {
  size_t start{s.find_first_not_of(WHITESPACE)};
  return (start == std::string::npos) ? "" : s.substr(start);
}

std::string Rtrim(const std::string &s) {
  size_t end{s.find_last_not_of(WHITESPACE)};
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string Trim(const std::string &s) {
  return Ltrim(Rtrim(s));
}

std::string ToLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), 
                   [](unsigned char c){ return std::tolower(c); } 
                  );
    return s;
}

int ParseString(std::string &inp, std::string &outp1, std::string &outp2, std::string &delimiter) {
  size_t pos = inp.find(delimiter);
  int delimSize{static_cast<int>(delimiter.size())};
  int miss{1};
  if (pos != std::string::npos) {
    while (inp[pos + miss] == delimiter[miss]) {
      ++miss;
    }
    outp1 = inp.substr(0, pos);
    outp2 = inp.substr(pos + miss, inp.length());
    return pos;
  } else {
    outp1 = inp.substr(0, inp.length());
    outp2 = "";
    return -1;
  }
}
