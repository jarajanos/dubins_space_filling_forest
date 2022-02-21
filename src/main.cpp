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
    std::stringstream ss;
    ss << "Solving problem file " << argv[1];
    INFO(ss.str());
    
    config = YAML::LoadFile(argv[1]);
  } catch (...) {
    ERROR("Error loading configuration file!");
    exit(1);
  }

  std::string dim{config["problem"]["dimension"].as<std::string>()};
  if (dim == "2D") {
    Problem<Point2D> problem;
    problem.Repetition = repetition;
    problem.Dimension = D2;
    SolveProblem(config, problem);
  } else if (dim == "2DDubins") {
    Problem<Point2DDubins> problem;
    problem.Repetition = repetition;
    problem.Dimension = D2Dubins;
    SolveProblem(config, problem);
  } else if (dim == "3D") {
    Problem<Point3D> problem;
    problem.Repetition = repetition;
    problem.Dimension = D3;
    SolveProblem(config, problem);
  } else if (dim == "3DDubins") {
    Problem<Point3DDubins> problem;
    problem.Repetition = repetition;
    problem.Dimension = D3Dubins;
    SolveProblem(config, problem);
  } else {
    ERROR("Invalid value of \"dimension\" node in \"problem\" root node.");
    exit(1);
  }

  return 0;
}

template<class R>
void SolveProblem(YAML::Node &config, Problem<R> &problem) {
  ParseFile(config, problem);
  std::unique_ptr<Solver<R>> solver;
  if (problem.Solver == SFF) {
    solver = std::make_unique<SpaceForest<R>>(problem);
  } else if (problem.Solver == RRT) {
    solver = std::make_unique<RapidExpTree<R>>(problem);
  } else if (problem.Solver == LazyRRT) {
    solver = std::make_unique<LazyTSP<R>>(problem);
  } else if (problem.Solver == PRM) {
    solver = std::make_unique<ProbRoadMaps<R>>(problem);
  } else if (problem.Solver == LazySFF) {
    solver = std::make_unique<LazySpaceForest<R>>(problem);
  } else {
    ERROR("Unimplemented problem solver");
  }

  solver->Solve();
}

template<class R>
void ParseFile(YAML::Node &config, Problem<R> &problem) {
  double scale;
  R pos;
  FileStruct tempFile;
  try {
    // problem
    YAML::Node node{config["problem"]};
    if (!node.IsDefined()) {
      throw std::invalid_argument("invalid problem node");
    }

    YAML::Node subNode{node["solver"]};
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid solver node in \"problem\" root node");
    }
    std::string value{subNode.as<std::string>()};
    if (!strcmp(value.c_str(),"sff")) {
      problem.Solver = SFF;
    } else if (!strcmp(value.c_str(), "rrt")) {
      problem.Solver = RRT;
    } else if (!strcmp(value.c_str(), "lazy-rrt")) {
      problem.Solver = LazyRRT;
    } else if (!strcmp(value.c_str(), "prm")) {
      problem.Solver = PRM;
    } else if (!strcmp(value.c_str(), "lazy-sff")) {
      problem.Solver = LazySFF;
    } else {
      throw std::invalid_argument("unknown solver type in \"problem\" root node, use either sff, rrt or lazy");
    }

    subNode = node["optimize"];
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid optimize node in \"problem\" root node");
    }
    problem.Optimize = subNode.as<bool>();

    subNode = node["iterations"];
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid iterations node in \"problem\" root node");
    }
    problem.MaxIterations = subNode.as<int>();

    // delimiters
    node = config["delimiters"];
    if (node.IsDefined()) {
      subNode = node["standard"];
      if (subNode.IsDefined()) {
        Obstacle<R>::Delimiter = subNode.as<std::string>();
      }

      subNode = node["name"];
      if (subNode.IsDefined()) {
        Obstacle<R>::NameDelimiter = subNode.as<std::string>();
      }
    }

    // TSP solver
    node = config["TSP-solver"];
    if (!node.IsDefined() && problem.Solver == LazyRRT) {
      throw std::invalid_argument("missing TSP solver parameters for Lazy solver!");
    } else if (node.IsDefined()) {
      problem.ComputeTSP = true;

      subNode = node["type"];
      if (!subNode.IsDefined()) {
        throw std::invalid_argument("invalid type node in \"TSP-solver\" root node");
      }

      std::string inpTsp{ToLower(subNode.as<std::string>())};
      if (inpTsp == "concorde") {
        problem.TspType = Concorde;
      } else if (inpTsp == "lkh") {
        problem.TspType = LKH;
      } else {
        throw std::invalid_argument("unknown type node in \"TSP-solver\" root node");
      }
      
      subNode = node["path"];
      if (subNode.IsDefined()) {
        problem.TspSolver = subNode.as<std::string>();
      } else {
        switch (problem.TspType) {
          case Concorde:
            problem.TspSolver = DEFAULT_CONCORDE;
            break;
          case LKH:
            problem.TspSolver = DEFUALT_LKH;
            break;
        }
      }
    }

    // algorithm node
    node = config["algorithm"];
    if (!node.IsDefined()) {
      throw std::invalid_argument("invalid \"algorithm\" root node");
    }
    subNode = node["m2r-ratio"];
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid m2r-ratio node in \"algorithm\" root node");
    }
    scale = subNode.as<double>();
    subNode = node["misses"];
    if (subNode.IsDefined()) {
      problem.MaxMisses = subNode.as<int>();
    }
    subNode = node["dubins-radius"];
    if ((problem.Dimension == D2Dubins || problem.Dimension == D3Dubins) && !subNode.IsDefined()) {
      throw std::invalid_argument("dubins radius missing");
    } else if (subNode.IsDefined()) {
      problem.DubinsRadius = subNode.as<double>() * scale;
      Point2DDubins::DubinsRadius = subNode.as<double>() * scale;
      Point3DDubins::DubinsRadius = subNode.as<double>() * scale;
    }
    subNode = node["dubins-resolution"];
    if ((problem.Dimension == D2Dubins || problem.Dimension == D3Dubins) && !subNode.IsDefined()) {
      throw std::invalid_argument("dubins resolution missing");
    } else if (subNode.IsDefined()) {
      problem.DubinsResolution = subNode.as<int>();
    }
    subNode = node["bias"];
    if (subNode.IsDefined()) {
      problem.PriorityBias = subNode.as<double>();
    }

    if (problem.Solver == LazyRRT && problem.PriorityBias != 0) {
      throw std::invalid_argument("priority bias for Lazy solver is not implemented");
    }

    subNode = node["prm-connections"];
    if (problem.Solver == PRM && !problem.Optimize && !subNode.IsDefined()) {
      throw std::invalid_argument("number of connections for classic sPRM must be defined!");
    } else if (subNode.IsDefined()) {
      problem.PrmConnections = subNode.as<int>();
    }

    // robot node
    node = config["robot"];
    if (!node.IsDefined()) {
      throw std::invalid_argument("invalid \"robot\" root node");
    }
    if (GetFile(node, tempFile)) {
      throw std::invalid_argument("invalid file node in \"robot\" root node");
    }
    problem.Env.Robot = std::make_shared<Obstacle<R>>(tempFile.fileName, tempFile.type, scale);

    // parse range
    node = config["range"];
    if (!node.IsDefined()) {
      throw std::invalid_argument("invalid \"range\" root node");
    }
    subNode = node["autodetect"];
    problem.AutoRange = (subNode.IsDefined() && subNode.as<bool>());

    if (!problem.AutoRange) {
      std::string tempText;
      subNode = node['x'];
      if (!subNode.IsDefined()) {
        throw std::invalid_argument("invalid x node in \"range\" root node");
      }
      tempText = subNode.as<std::string>();
      problem.Env.Limits.Parse(tempText, scale, 0);

      subNode = node['y'];
      if (!subNode.IsDefined()) {
        throw std::invalid_argument("invalid y node in \"range\" root node");
      }
      tempText = subNode.as<std::string>();
      problem.Env.Limits.Parse(tempText, scale, 1);

      subNode = node['z'];
      if (!subNode.IsDefined() && (problem.Dimension == D3 || problem.Dimension == D3Dubins)) {
        throw std::invalid_argument("invalid z node in \"range\" root node");
      } else if (subNode.IsDefined()) {
        tempText = subNode.as<std::string>();
        problem.Env.Limits.Parse(tempText, scale, 2);
      }

      subNode = node["pitch"];
      if (!subNode.IsDefined()) {
        if (problem.Dimension == D3Dubins) {
          WARN("missing pitch node in \"range\" root node, settings to +/- (3.14 / 2)");
          problem.Env.Limits.maxs[3] = M_PI_2;
          problem.Env.Limits.mins[3] = -M_PI_2;
        } else if (problem.Dimension == D3) {
          problem.Env.Limits.maxs[3] = __DBL_MAX__;
          problem.Env.Limits.mins[3] = -__DBL_MAX__;
        }
      } else {
        tempText = subNode.as<std::string>();
        problem.Env.Limits.Parse(tempText, scale, 3);
        if (problem.Env.Limits.mins[3] != -problem.Env.Limits.maxs[3]) {
          std::invalid_argument("invalid pitch node in \"range\" root node: angles must be symmetric");
        }
      }
      Point3DDubins::MaxPitch = problem.Env.Limits.maxs[3];
    }

    // parse obstacles
    node = config["obstacles"];
    if (!node.IsDefined()) {
      problem.Env.HasMap = false;
      problem.Env.ScaleFactor = scale;
    } else {
      problem.Env.ScaleFactor = scale;

      for (YAML::const_iterator it=node.begin(); it != node.end(); ++it) {
        YAML::Node tempNode{*it};
        if (GetFile(tempNode, tempFile)) {
          throw std::invalid_argument("invalid file attribute in \"obstacles\" root sequence");
        }

        subNode = (*it)["position"];
        if (!subNode.IsDefined()) {
          pos = R();
        } else {
          pos = R(subNode.as<std::string>());
        }

        Obstacle<R> &obst{problem.Env.Obstacles.emplace_back(tempFile.fileName, tempFile.type, pos, scale)};

        if (problem.AutoRange) {
          Range obstRange{obst.GetRange()};
          problem.Env.ProcessLimits(obstRange);
        }
      }
    }

    // parse cities
    node = config["cities"];
    if (!node.IsDefined()) {
      throw std::invalid_argument("invalid \"cities\" root node");
    } else {
      int numCities{0};
      for (YAML::const_iterator it=node.begin(); it != node.end(); ++it) {
        problem.Roots.emplace_back(it->as<std::string>(), scale);
        ++numCities;
      }

      if (problem.Solver == RRT && problem.Optimize && numCities > 1) {
        throw std::invalid_argument("Multi-T-RRT* is undefined");
      } 

      if (numCities > 1 && problem.DubinsResolution % 2 == 1) {
        throw std::invalid_argument("The dubins resolution for multi-goal planning has to be an even number greater than 1, to ensure continuity of the path");
      }
    }

    // parse goal
    node = config["goal"];
    if (node.IsDefined()) {
      if (problem.Solver == LazyRRT) {
        throw std::invalid_argument("single point path planning not defined for Lazy solver (use RRT/RRT* solver instead)");
      } else if (problem.Roots.size() > 1) {
        WARN("Multi-source planning with one goal has not been tested");
      }

      if (problem.ComputeTSP) {
        WARN("TSP solver specified for single goal -- TSP computation switched off!");
        problem.ComputeTSP = false;
      }
      problem.HasGoal = true;
      problem.Goal = R(node.as<std::string>(), scale);
    }
    
    if (!problem.HasGoal && problem.PriorityBias != 0 && problem.Solver == RRT) {
      throw std::invalid_argument("Multi-T-RRT with bias is undefined");
    }

    // parse distances
    node = config["distances"];
    if (!node.IsDefined()) {
      throw std::invalid_argument("invalid \"distances\" root node");
    }
    subNode = node["collision"];
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid collision node in \"distances\" root node");
    }
    problem.CollisionDist = subNode.as<double>() * scale;

    subNode = node["sampling"];
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid sampling node in \"distances\" root node");
    }
    problem.SamplingDist = subNode.as<double>() * scale;

    subNode = node["tree"];
    if (!subNode.IsDefined()) {
      throw std::invalid_argument("invalid tree node in \"distances\" root node");
    }
    problem.DistTree = subNode.as<double>() * scale;

    // parse save
    node = config["save"];
    if (node.IsDefined()) {
      subNode = node["goals"];
      if (!GetFile(subNode, tempFile, problem.Repetition)) {
        problem.SaveOpt = problem.SaveOpt + SaveGoals;
        problem.FileNames[SaveGoals] = tempFile;
      }

      subNode = node["tree"];
      if (!GetFile(subNode, tempFile, problem.Repetition)) {
        problem.SaveOpt = problem.SaveOpt + SaveTree;
        problem.FileNames[SaveTree] = tempFile;

        subNode = subNode["frequency"];
        if (subNode.IsDefined() && subNode.as<int>() != 0) {
          problem.SaveFreq[SaveTree] = subNode.as<int>();
        }
      }

      subNode = node["roadmap"];
      if (!GetFile(subNode, tempFile, problem.Repetition)) {
        problem.SaveOpt = problem.SaveOpt + SaveRoadmap;
        problem.FileNames[SaveRoadmap] = tempFile;
      }

      subNode = node["params"];
      if (!GetFile(subNode, tempFile, problem.Repetition, false)) {
        problem.SaveOpt = problem.SaveOpt + SaveParams;
        problem.FileNames[SaveParams] = tempFile;

        subNode = subNode["id"];
        if (subNode.IsDefined()) {
          problem.ID = subNode.as<std::string>();
        }
      }

      subNode = node["TSP-file"];
      if (!GetFile(subNode, tempFile, problem.Repetition)) {
        problem.SaveOpt = problem.SaveOpt + SaveTSPFile;
        problem.FileNames[SaveTSPFile] = tempFile;
      }

      subNode = node["TSP-paths"];
      if (!GetFile(subNode, tempFile, problem.Repetition)) {
        if (!problem.ComputeTSP) {
          WARN("TSP paths specified although the TSP is not solved -- maybe you forgot to specify TSP solver?");
        }
        problem.SaveOpt = problem.SaveOpt + SaveTSPPaths;
        problem.FileNames[SaveTSPPaths] = tempFile;
      }

      subNode = node["frontiers"];
      if (!GetFile(subNode, tempFile, problem.Repetition)) {
        if (problem.Solver != SFF) {
          throw std::invalid_argument("frontier output defined only for SFF-based solvers");
        }
        problem.SaveOpt = problem.SaveOpt + SaveFrontiers;
        problem.FileNames[SaveFrontiers] = tempFile;

        subNode = subNode["frequency"];
        if (subNode.IsDefined() && subNode.as<int>() != 0) {
          problem.SaveFreq[SaveFrontiers] = subNode.as<int>();
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
  if (!node.IsDefined()) {
    return true;
  }
  
  subNode = node["path"];
  if (!subNode.IsDefined()) { 
    return true;
  }  
  file.fileName = subNode.as<std::string>();
  if (repetition != 0 && includeIter) {
    size_t iter{file.fileName.find_last_of('.')};
    file.fileName.insert(iter, '_' + std::to_string(repetition));
  }

  subNode = node["type"];
  if (!subNode.IsDefined()) {
    file.type = Map;
  } else {
    std::string value{subNode.as<std::string>()};
    if (!strcmp(value.c_str(), "map")) {
      file.type = Map;
    } else if (!strcmp(value.c_str(), "obj")) {
      file.type = Obj;
    } else {
      throw std::invalid_argument("invalid type node in file node");
    }
  }

  return false;
}
