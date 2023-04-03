vpath %.cpp src/
vpath %.h include/
# CC compiler options:

#CC = g++ #uncomment for CPU
#CC = /opt/rocm-5.2.5/bin/hipcc
CC = /usr/local/cuda-11.8/bin/nvcc #uncomment for Nvidia
#CC_FLAGS= -g -w -O3 -std=c++17 -fPIC -fopenmp #uncommend for openmp
SM=86
GENCODE_FLAGS = -gencode arch=compute_$(SM),code=compute_$(SM)
CC_FLAGS= -g -w -O3 -std=c++17 -rdc=true -cudart=shared $(GENCODE_FLAGS)
#GPUCC_FLAGS= $(CC_FLAGS) -fgpu-rdc

#INCLUDES = $(BASEGPUPATH)hiprand/lib
BASEGPUPATH = /opt/nvidia/hpc_sdk/Linux_x86_64/22.11
CC_INCLUDES = -I$(BASEGPUPATH)/math_libs/include
LIBRARY_PATH = -L$(BASEGPUPATH)/math_libs/lib64 -L$(BASEGPUPATH)/cuda/lib64
LIBRARIES = -lcudart -lcurand
## Project file structure ##


MAIN = main
SOURCES = double3.cpp optionsParser.cpp integrator.cpp particle.cpp simulation.cpp
INCLUDES = $(SOURCES:.cpp=.h)
OBJECTS = $(SOURCES:.cpp=.o) $(MAIN).o
BUILD = build/
EXECS =

all: $(MAIN)

$(MAIN): $(OBJECTS)
	$(CC) $(CC_FLAGS) $(CC_INCLUDES) $(LIBRARY_PATH) $(addprefix $(BUILD),$(OBJECTS)) $(LIBRARIES) -o $(EXECS)$@ 

$(MAIN).o : $(MAIN).cpp $(INCLUDES)
	$(CC) $(CC_FLAGS) $(CC_INCLUDES) $(LIBRARY_PATH) -c $< -o $(BUILD)$@ $(LIBRARIES)

%.o : %.cpp %.h
	$(CC) $(CC_FLAGS) $(CC_INCLUDES) $(LIBRARY_PATH) -c $< -o $(BUILD)$@ $(LIBRARIES)

$(shell mkdir -p $(BUILD) $(EXECS))  

clean:
	rm -f bin/* *.o $(MAIN)
