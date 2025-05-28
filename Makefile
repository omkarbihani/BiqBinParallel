#### BiqBin makefile ####

# container image name
IMAGE ?= parallel-biqbin
# container image tag
TAG ?= 1.0.0
DOCKER_BUILD_PARAMS ?=

# Directories
OBJ = Obj
OBJC = Objc

# Compiler
CC = mpicc
CPP = mpic++

LINALG 	 = -lopenblas -lm 
OPTI     = -O3 -ffast-math -fexceptions -fPIC -fno-common
CPPOPTI     = -O3 -fexceptions -fPIC -fno-common

PYBOOST ?= -I/home/roman/anaconda3/include/python3.12 \
		  -I/home/roman/anaconda3/include \
          -L/home/roman/anaconda3/lib  -lboost_python312  -lboost_numpy312 -lpython3.12

INCLUDES += $(PYBOOST)
# Python module (Boost)
PYMOD_SRC = wrapper.cpp
PYMOD_OUT = solver.so
# binary
BINS =  biqbin

# BiqBin objects
BBOBJSC = $(OBJC)/bundle.o $(OBJC)/allocate_free.o $(OBJC)/bab_functions.o \
	 	 $(OBJC)/bounding.o $(OBJC)/cutting_planes.o \
         $(OBJC)/evaluate.o $(OBJC)/heap.o $(OBJC)/ipm_mc_pk.o \
         $(OBJC)/heuristic.o $(OBJC)/main.o $(OBJC)/operators.o \
         $(OBJC)/process_input.o $(OBJC)/qap_simulated_annealing.o \
		 $(OBJC)/bqp_data_processing.o

# BiqBin objects
BBOBJS = $(OBJ)/bundle.o $(OBJ)/allocate_free.o $(OBJ)/bab_functions.o \
	 	 $(OBJ)/bounding.o $(OBJ)/cutting_planes.o \
         $(OBJ)/evaluate.o $(OBJ)/heap.o $(OBJ)/ipm_mc_pk.o \
         $(OBJ)/heuristic.o $(OBJ)/main.o $(OBJ)/operators.o \
         $(OBJ)/process_input.o $(OBJ)/qap_simulated_annealing.o \
		 $(OBJ)/bqp_data_processing.o

BBOBJSCPP =	$(OBJ)/wrapper.o

# All objects
OBJS = $(BBOBJS) $(BBOBJSCPP)

CFLAGS = $(OPTI) -Wall -W -pedantic 
CPPFLAGS = $(CPPOPTI) -Wall -W -pedantic 

#### Rules ####

.PHONY : all clean test tests

# Default rule is to create all binaries #
all: clean $(BINS) $(PYMOD_OUT)

	
clean-output:
	rm -f rudy/*.output*
	rm -f tests/rudy/*.output*
	rm -f tests/qubos/*.output*

# Clean rule #
clean: clean-output
	rm -rf $(BINS) $(OBJS) $(PYMOD_OUT) $(BBOBJSC)

# Rules for binaries
$(BINS): $(BBOBJSC)
	$(CC) -o $@ $^ $(INCLUDES) $(LIB) $(CFLAGS) -DPURE_C $(LINALG)

$(OBJC)/%.o: %.c
	$(CC) $(CFLAGS) -DPURE_C $(INCLUDES) -c -o $@ $<
# BiqBin code rules
$(OBJ)/%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c -o $@ $<

# Python module rule
$(PYMOD_OUT): $(OBJS)
	$(CPP) -o $@ $^ -shared -fPIC $(INCLUDES) $(LIB) $(LINALG) -Wl,--no-undefined

# Tests
test: clean-output
	./test.sh \
	"mpiexec -n 8 ./$(BINS)" tests/rudy/g05_60.0 tests/rudy/g05_60.0-expected_output params

test-all: clean-output
	./test_all.sh ./biqbin 60 8
	./test_all.sh ./biqbin 80 8
	./test_all.sh ./biqbin 100 8

test-python-maxcut: clean-output
	./test.sh \
	"mpiexec -n 3 python3 biqbin_maxcut.py" tests/rudy/g05_60.0 tests/rudy/g05_60.0-expected_output params

test-python-maxcut-all: clean-output
	./test_all.sh "python3 biqbin_maxcut.py" 60 8
	./test_all.sh "python3 biqbin_maxcut.py" 80 8
	./test_all.sh "python3 biqbin_maxcut.py" 100 8

test-python-qubo: clean-output
	./qubo_test.sh \
	"mpiexec -n 8 python3 run_qubo_test.py" tests/qubos/40/kcluster40_025_10_1.json params

test-python-qubo-all-small:
	./test_all_qubo.sh tests/qubos/40 8
	./test_all_qubo.sh tests/qubos/80 8

test-python-qubo-all-large:
	./test_all_qubo.sh tests/qubos/100 8
	./test_all_qubo.sh tests/qubos/120 8
	./test_all_qubo.sh tests/qubos/140 8
	./test_all_qubo.sh tests/qubos/160 8

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
	
