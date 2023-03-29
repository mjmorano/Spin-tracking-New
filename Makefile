##########################################################

# CC compiler options:
#CC=g++
#CC = /opt/rocm-5.2.5/bin/hipcc
CC = g++
#CC_FLAGS= -g -O3 -std=c++17 -fopenmp -target x86_64-pc-linux-gnu -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx1030 -lm
#CC_FLAGS = -g -O3 -std=c++17 -L/opt/rocm-5.2.5/hiprand/lib
CC_FLAGS = -O3 -w -std=c++17 -fopenmp
#CC_LIBS= -L/opt/rocm-5.2.5/hiprand/lib

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

##########################################################

## Make variables ##

# Target executable name:
EXE = run

# Object files:
OBJS = $(OBJ_DIR)/double3.o $(OBJ_DIR)/integrator.o $(OBJ_DIR)/particle.o $(OBJ_DIR)/main.o

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
$(EXE) : $(OBJS)
	$(CC) $(CC_FLAGS) $(OBJS) -o $@ 

# Compile main .cpp file to object files:
$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp include/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE)
