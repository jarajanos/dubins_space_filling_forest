/**
 * @file main.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "main.h"

int main(int argc, char *argv[]) {
  YAML::Node config;
  int repetition{0};

  if (argc < 2) {
    return 2;
  } else if (argc == 3) {
    repetition = std::stoi(argv[2]);
  }

  try {
    config = YAML::LoadFile(".//test_building.yaml");
  } catch (...) {
    ERROR("Error loading configuration file!");
    exit(1);
  }

  std::string dim{config["problem"]["dimension"].as<std::string>()};
  if (dim == "2D") {
    Problem<Point2D> problem;
    problem.repetition = repetition;
    problem.dimension = D2;
    SolveProblem(config, problem);
  } else if (dim == "2DDubins") {
    Problem<Point2DDubins> problem;
    problem.repetition = repetition;
    problem.dimension = D2Dubins;
    SolveProblem(config, problem);
  } else if (dim == "3D") {
    Problem<Point3D> problem;
    problem.repetition = repetition;
    problem.dimension = D3;
    SolveProblem(config, problem);
  } else if (dim == "3DDubins") {

  } else {
    ERROR("Invalid value of \"dimension\" node in \"problem\" root node.");
    exit(1);
  }

  return 0;
}

template<class R>
void SolveProblem(YAML::Node &config, Problem<R> &problem) {
  ParseFile(config, problem);
  // TODO
}

template<class R>
void ParseFile(YAML::Node &config, Problem<R> &problem) {
  double scale;
  R pos;
  FileStruct tempFile;
  try {
    // problem
    YAML::Node node{config["problem"]};
    if (node.IsNull()) {
      throw std::invalid_argument("invalid problem node");
    }

    YAML::Node subNode{node["solver"]};
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid solver node in \"problem\" root node");
    }
    std::string value{subNode.as<std::string>()};
    if (!strcmp(value.c_str(),"sff")) {
      problem.solver = SFF;
    } else if (!strcmp(value.c_str(), "rrt")) {
      problem.solver = RRT;
    } else if (!strcmp(value.c_str(), "lazy")) {
      problem.solver = Lazy;
    } else {
      throw std::invalid_argument("unknown solver type in \"problem\" root node, use either sff, rrt or lazy");
    }

    subNode = node["optimize"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid optimize node in \"problem\" root node");
    }
    problem.optimize = subNode.as<bool>();

    subNode = node["iterations"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid iterations node in \"problem\" root node");
    }
    problem.maxIterations = subNode.as<int>();

    // delimiters
    node = config["delimiters"];
    if (!node.IsNull()) {
      subNode = node["standard"];
      if (!subNode.IsNull()) {
        Obstacle<R>::Delimiter = subNode.as<std::string>();
      }

      subNode = node["name"];
      if (!subNode.IsNull()) {
        Obstacle<R>::NameDelimiter = subNode.as<std::string>();
      }
    }

    // TSP solver for Lazy-RRT
    node = config["TSP-solver"];
    if (node.IsNull() && problem.solver == Lazy) {
      throw std::invalid_argument("missing TSP solver parameters for Lazy solver!");
    } else if (!node.IsNull()) {
      if (problem.solver != Lazy) {
        WARN("TSP solver is called only in Lazy solver algorithm -- defined TSP parameters are redundant and will not be used.");
      }

      subNode = node["path"];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid path node in \"TSP-solver\" root node");
      }
      problem.tspSolver = subNode.as<std::string>();

      subNode = node["type"];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid type node in \"TSP-solver\" root node");
      }
      problem.tspType = subNode.as<std::string>();
    }

    // algorithm node
    node = config["algorithm"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalid \"algorithm\" root node");
    }
    subNode = node["m2r-ratio"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid m2r-ratio node in \"algorithm\" root node");
    }
    scale = subNode.as<double>();
    subNode = node["misses"];
    if (!subNode.IsNull()) {
      problem.maxMisses = subNode.as<int>();
    }
    subNode = node["bias"];
    if (!subNode.IsNull()) {
      problem.priorityBias = subNode.as<double>();
    }

    if (problem.solver == Lazy && problem.priorityBias != 0) {
      throw std::invalid_argument("priority bias for Lazy solver is not implemented!");
    }

    // robot node
    node = config["robot"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalid \"robot\" root node");
    }
    if (GetFile(node, tempFile)) {
      throw std::invalid_argument("invalid file node in \"robot\" root node");
    }
    problem.environment.Robot = new Obstacle<R>(tempFile.fileName, tempFile.type, scale);

    // parse range
    node = config["range"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalide \"range\" root node");
    }
    subNode = node["autoDetect"];
    problem.autoRange = (!subNode.IsNull() && subNode.as<bool>());

    if (!problem.autoRange) {
      std::string tempText;
      subNode = node['x'];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid x node in \"range\" root node");
      }
      tempText = subNode.as<std::string>();
      problem.environment.limits.Parse(tempText, scale, 0);

      subNode = node['y'];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid y node in \"range\" root node");
      }
      tempText = subNode.as<std::string>();
      problem.environment.limits.Parse(tempText, scale, 1);

      subNode = node['z'];
      if (subNode.IsNull() && (problem.dimension == D3 || problem.dimension == D3Dubins)) {
        throw std::invalid_argument("invalid z node in \"range\" root node");
      } else if (!subNode.IsNull()) {
        tempText = subNode.as<std::string>();
        problem.environment.limits.Parse(tempText, scale, 2);
      }
    }

    // parse obstacles
    node = config["obstacles"];
    if (node.IsNull()) {
      problem.environment.HasMap = false;
    } else {
      problem.environment.ScaleFactor = scale;

      for (YAML::const_iterator it=node.begin(); it != node.end(); ++it) {
        YAML::Node tempNode{*it};
        if (GetFile(tempNode, tempFile)) {
          throw std::invalid_argument("invalid file attribute in \"obstacles\" root sequence");
        }

        subNode = (*it)["position"];
        if (subNode.IsNull()) {
          pos = R();
        } else {
          pos = R(subNode.as<std::string>());
        }

        Obstacle<R> &obst{problem.environment.Obstacles.emplace_back(tempFile.fileName, tempFile.type, pos, scale)};

        if (problem.autoRange) {
          problem.environment.processLimits(obst.GetRange());
        }
      }
    }

    // parse cities
    node = config["cities"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalid \"cities\" root node");
    } else {
      int numCities{0};
      for (YAML::const_iterator it=node.begin(); it != node.end(); ++it) {
        problem.roots.emplace_back(it->as<std::string>(), scale);
        ++numCities;
      }

      if (problem.solver == RRT && problem.optimize && numCities > 1) {
        throw std::invalid_argument("Multi-T-RRT* is undefined");
      } 
    }

    // parse goal
    node = config["goal"];
    if (!node.IsNull()) {
      if (problem.solver == Lazy) {
        throw std::invalid_argument("single point path planning not defined for Lazy solver (use RRT/RRT* solver instead)");
      } else if (problem.roots.size() > 1) {
        WARN("Multi-source planning with one goal has not been tested");
      }
      problem.hasGoal = true;
      problem.goal = R(node.as<std::string>(), scale);
    }
    
    if (!problem.hasGoal && problem.priorityBias != 0 && problem.solver == RRT) {
      throw std::invalid_argument("Multi-T-RRT with bias is undefined");
    }

    // parse distances
    node = config["distances"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalid \"distances\" root node");
    }
    subNode = node["collision"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid collision node in \"distances\" root node");
    }
    problem.collisionDist = subNode.as<double>() * scale;

    subNode = node["sampling"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid sampling node in \"distances\" root node");
    }
    problem.samplingDist = subNode.as<double>() * scale;

    subNode = node["tree"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid tree node in \"distances\" root node");
    }
    problem.distTree = subNode.as<double>() * scale;

    // parse save
    node = config["save"];
    if (node.IsNull()) {
      subNode = node["goals"];
      if (!GetFile(subNode, tempFile, problem.repetition)) {
        problem.saveOptions = problem.saveOptions | SaveGoals;
        problem.fileNames[SaveGoals] = tempFile;
      }

      subNode = node["tree"];
      if (!GetFile(subNode, tempFile, problem.repetition)) {
        problem.saveOptions = problem.saveOptions | SaveTree;
        problem.fileNames[SaveTree] = tempFile;

        subNode = subNode["frequency"];
        if (!subNode.IsNull() && subNode.as<int>() != 0) {
          problem.saveFreq[SaveTree] = subNode.as<int>();
        }
      }

      subNode = node["roadmap"];
      if (!GetFile(subNode, tempFile, problem.repetition)) {
        problem.saveOptions = problem.saveOptions | SaveRoadmap;
        problem.fileNames[SaveRoadmap] = tempFile;
      }

      subNode = node["params"];
      if (!GetFile(subNode, tempFile, problem.repetition, false)) {
        problem.saveOptions = problem.saveOptions | SaveParams;
        problem.fileNames[SaveParams] = tempFile;

        subNode = subNode["id"];
        if (!subNode.IsNull()) {
          problem.id = subNode.as<std::string>();
        }
      }

      subNode = node["TSP"];
      if (!GetFile(subNode, tempFile, problem.repetition)) {
        problem.saveOptions = problem.saveOptions | SaveTSP;
        problem.fileNames[SaveTSP] = tempFile;
      }

      subNode = node["frontiers"];
      if (!GetFile(subNode, tempFile, problem.repetition)) {
        if (problem.solver != SFF) {
          throw std::invalid_argument("frontier output defined only for SFF-based solvers");
        }
        problem.saveOptions = problem.saveOptions | SaveFrontiers;
        problem.fileNames[SaveFrontiers] = tempFile;

        subNode = subNode["frequency"];
        if (!subNode.IsNull() && subNode.as<int>() != 0) {
          problem.saveFreq[SaveFrontiers] = subNode.as<int>();
        }
      } 
    }


  } catch (const std::invalid_argument e) {
    std::stringstream message;
    message << "Problem loading configuaration file: " << e.what();
    ERROR(message.str());
    exit(1);
  }
}

bool GetFile(YAML::Node &node, FileStruct &file, int repetition, bool includeIter) {
  YAML::Node subNode;
  if (node.IsNull()) {
    return true;
  }
  
  subNode = node["file"];
  if (subNode.IsNull()) { 
    return true;
  }  
  file.fileName = subNode.as<std::string>();
  if (repetition != 0 && includeIter) {
    size_t iter{file.fileName.find_last_of('.')};
    file.fileName.insert(iter, '_' + std::to_string(repetition));
  }

  subNode = node["type"];
  if (subNode.IsNull()) {
    file.type = Map;
  } 
  std::string value{subNode.as<std::string>()};
  if (!strcmp(value.c_str(), "map")) {
    file.type = Map;
  } else if (!strcmp(value.c_str(), "obj")) {
    file.type = Obj;
  } else {
    throw std::invalid_argument("invalid type node in file node");
  }

  return false;
}
