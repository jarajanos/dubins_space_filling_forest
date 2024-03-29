TARGET := main
OBJDIR := obj
DBGDIR := debug
RELDIR := release

CPLEX_ROOT_DIR := /opt/ibm/ILOG/CPLEX_Studio_Community201/cplex

SOURCES := $(wildcard src/*.cpp)
#SOURCES := src/main.cpp src/random-generator.cpp
SOURCES := $(SOURCES:src/%=%)
OBJECTS := $(patsubst %.c,%.o, $(patsubst %.cpp,%.o,$(SOURCES)))

INCLUDE := -I. -I./lib/rapid-2.01 -I./lib/flann/src/cpp -I./lib/ann/ann/include -I./lib/yaml-cpp/include -I./lib/gdip/gdip/include -I./lib/RapidQuadrocopterTrajectories/C++
LIBPATH := -L./lib/ -L./lib/ann/ -L./lib/flann/ -L./lib/rapid-2.01 -L./lib/yaml-cpp/build -L./lib/gdip/gdip/lib -L./lib/RapidQuadrocopterTrajectories/lib
LIBS := -lgmp -lRAPID -llz4 -lyaml-cpp -lopendubins_core -lstdc++fs -lquad-trajectories

CXXFLAGS := -std=c++17
CXX := g++-8

DBGEXE := $(DBGDIR)/$(TARGET)
DBGOBJS := $(addprefix $(OBJDIR)/, $(addprefix $(DBGDIR)/, $(OBJECTS)))
DBGFLAGS := $(CXXFLAGS) -g -O0

RELEXE := $(RELDIR)/$(TARGET)
RELOBJS := $(addprefix $(OBJDIR)/, $(addprefix $(RELDIR)/, $(OBJECTS)))
RELFLAGS := $(CXXFLAGS) -O3 -fno-math-errno

.PHONY: all clean debug release prep rapid flann yaml install gdip concorde lkh solver trajectories

all: release

############ DEBUG ###############
debug: prep $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CXX) $(DBGFLAGS) $(INCLUDE) -o $(DBGEXE) $^ $(LIBPATH) $(LIBS)

$(OBJDIR)/$(DBGDIR)/%.o: src/%.cpp src/%.h
	$(CXX) $(DBGFLAGS) $(INCLUDE) -c $< -o $@

############ RELEASE ###############
release: prep $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CXX) $(RELFLAGS) $(INCLUDE) -o $(RELEXE) $^ $(LIBPATH) $(LIBS)

$(OBJDIR)/$(RELDIR)/%.o: src/%.cpp src/%.h
	$(CXX) $(RELFLAGS) $(INCLUDE) -c $< -o $@

############# COMMON ##############	
prep:
	@echo "Checking directories..."
	@mkdir -p $(OBJDIR)
	@mkdir -p $(OBJDIR)/$(DBGDIR)
	@mkdir -p $(OBJDIR)/$(RELDIR)
	@mkdir -p $(DBGDIR)
	@mkdir -p $(RELDIR)

install: rapid flann yaml gdip trajectories release

solver: concorde lkh

rapid:
	@echo "Installing RAPID..."
	@$(MAKE) -C lib/rapid-2.01/

flann:
	@echo "Installing FLANN..."
	@cd ./lib/flann; mkdir build; cd build; cmake .. -DCMAKE_INSTALL_PREFIX=../install; make -j4; make install

yaml:
	@echo "Installing YAML-cpp..."
	@cd ./lib/yaml-cpp; mkdir build; cd build; cmake ..; make

gdip:
	@echo "Installing GDIP..."
	@cd ./lib/gdip/gdip; ./install.sh

concorde:
	@echo "Preparing Concorde solver..."
	@cd ./solver/concorde; ./install.sh $(CPLEX_ROOT_DIR)

lkh:
	@echo "Preparing LKH solver..."
	@cd ./solver/lkh; ./install.sh

trajectories:
	@echo "Preparing trajectory generator..."
	@cd ./lib/RapidQuadrocopterTrajectories; make

clean:
	@echo "Cleaning..."
	@rm -rf $(OBJDIR)
	@rm -rf $(RELDIR)
	@rm -rf $(DBGDIR)
