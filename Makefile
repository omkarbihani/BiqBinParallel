#### BiqBin makefile ####

# container image name
IMAGE ?= parallel-biqbin
# container image tag
TAG ?= 1.0.0
DOCKER_BUILD_PARAMS ?=

# Directories
WRAPPER_BUILD_DIR = build/wrapper
C_BUILD_DIR = build/c_build

# Compiler
CC = mpicc
CPP = mpic++

LINALG 	 = -lopenblas -lm 
OPTI     = -O3 -ffast-math -fexceptions -fPIC -fno-common
CPPOPTI  = -O3 -fexceptions -fPIC -fno-common -ffast-math

# Python module (Boost)
PYBOOST_INCLUDES ?= -I$(CONDA_PREFIX)/include/python3.12 -I$(CONDA_PREFIX)/include
PYBOOST_LIBS ?= -L$(CONDA_PREFIX)/lib -lboost_python312 -lboost_numpy312 -lpython3.12

INCLUDES += $(PYBOOST_INCLUDES)
LIB += $(PYBOOST_LIBS)

# Python module (Boost)
PYMOD_OUT = biqbin.so
# C only binary
BINS =  biqbin

# BiqBin objects
C_OBJS = $(C_BUILD_DIR)/bundle.o $(C_BUILD_DIR)/allocate_free.o $(C_BUILD_DIR)/bab_functions.o \
	 	 $(C_BUILD_DIR)/bounding.o $(C_BUILD_DIR)/cutting_planes.o \
         $(C_BUILD_DIR)/evaluate.o $(C_BUILD_DIR)/heap.o $(C_BUILD_DIR)/ipm_mc_pk.o \
         $(C_BUILD_DIR)/heuristic.o $(C_BUILD_DIR)/main.o $(C_BUILD_DIR)/operators.o \
         $(C_BUILD_DIR)/process_input.o $(C_BUILD_DIR)/qap_simulated_annealing.o \
		 $(C_BUILD_DIR)/bqp_data_processing.o

# BiqBin objects
OBJS =   $(WRAPPER_BUILD_DIR)/bundle.o $(WRAPPER_BUILD_DIR)/allocate_free.o $(WRAPPER_BUILD_DIR)/bab_functions.o \
	 	 $(WRAPPER_BUILD_DIR)/bounding.o $(WRAPPER_BUILD_DIR)/cutting_planes.o \
         $(WRAPPER_BUILD_DIR)/evaluate.o $(WRAPPER_BUILD_DIR)/heap.o $(WRAPPER_BUILD_DIR)/ipm_mc_pk.o \
         $(WRAPPER_BUILD_DIR)/heuristic.o $(WRAPPER_BUILD_DIR)/main.o $(WRAPPER_BUILD_DIR)/operators.o \
         $(WRAPPER_BUILD_DIR)/process_input.o $(WRAPPER_BUILD_DIR)/qap_simulated_annealing.o \
		 $(WRAPPER_BUILD_DIR)/bqp_data_processing.o

# All objects

CFLAGS = $(OPTI) -Wall -W -pedantic 
CPPFLAGS = $(CPPOPTI) -Wall -W -pedantic 

#### Rules ####

.PHONY : all clean test tests

# Default rule is to create all binaries #
all: clean $(BINS) $(PYMOD_OUT)

	
clean-output:
	rm -f rudy/*.output*
	rm -f tests/rudy/*.output*
	rm -f tests/qubos/*/*.output*

# Clean rule #
clean: clean-output
	rm -rf $(BINS) $(OBJS) $(PYMOD_OUT) $(C_OBJS) $(WRAPPER_BUILD_DIR)/wrapper.o

# Ensure output directories exist
$(WRAPPER_BUILD_DIR) $(C_BUILD_DIR):
	mkdir -p build
	mkdir -p $@

# Rules for binaries
$(BINS): $(C_OBJS)
	$(CC) -o $@ $^ $(INCLUDES) $(LIB) $(CFLAGS) -DPURE_C $(LINALG)

$(C_BUILD_DIR)/%.o: %.c  | $(C_BUILD_DIR)
	$(CC) $(CFLAGS) -DPURE_C $(INCLUDES) -c -o $@ $<

# BiqBin code rules
$(WRAPPER_BUILD_DIR)/%.o: %.c  | $(WRAPPER_BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(WRAPPER_BUILD_DIR)/wrapper.o: wrapper.cpp  | $(WRAPPER_BUILD_DIR)
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c -o $@ $<

# Python module rule
$(PYMOD_OUT): $(OBJS) build/wrapper/wrapper.o
	$(CPP) -o $@ $^ -shared -fPIC $(INCLUDES) $(LIB) $(LINALG) -Wl,--no-undefined

# Tests
test: clean-output
	tests/test.sh \
	"mpiexec -n 8 ./$(BINS)" tests/rudy/g05_60.0 tests/rudy/g05_60.0-expected_output params

test-all: clean-output
	tests/test_all.sh ./biqbin 60 8
	tests/test_all.sh ./biqbin 80 8
	tests/test_all.sh ./biqbin 100 8

test-python-maxcut: clean-output
	tests/test.sh \
	"mpiexec -n 3 python3 biqbin_maxcut.py" tests/rudy/g05_60.0 tests/rudy/g05_60.0-expected_output params

test-python-maxcut-all: clean-output
	tests/test_all.sh "python3 biqbin_maxcut.py" 60 8
	tests/test_all.sh "python3 biqbin_maxcut.py" 80 8
	tests/test_all.sh "python3 biqbin_maxcut.py" 100 8

test-python-qubo: clean-output
	tests/qubo_test.sh \
	"mpiexec -n 8 python3 -m tests.run_qubo_test" tests/qubos/40/kcluster40_025_10_1.json params

test-python-qubo-all:
	tests/test_all_qubo.sh tests/qubos/40 8
	tests/test_all_qubo.sh tests/qubos/80 8

docker: 
	docker build $(DOCKER_BUILD_PARAMS) --progress=plain -t $(IMAGE):$(TAG)  . 

docker-no-cache: 
	docker build --no-cache $(DOCKER_BUILD_PARAMS) --progress=plain -t $(IMAGE):$(TAG)  . 

docker-clean: 
	docker rmi -f $(IMAGE):$(TAG) 

docker-test:
	docker run --rm $(IMAGE):$(TAG) make test

docker-test-python:
	docker run --rm $(IMAGE):$(TAG) make test-python
	
