#### BiqBin makefile ####

# container image name
IMAGE ?= parallel-biqbin
# container image tag
TAG ?= 1.0.0
DOCKER_BUILD_PARAMS ?=

# Directories
OBJ = Obj

# Compiler
CC = mpicc
CPP = mpic++

LINALG 	 = -lopenblas -lm 
OPTI     = -O3 -ffast-math -fexceptions -fPIC -fno-common

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

#### Rules ####

.PHONY : all clean test tests

# Default rule is to create all binaries #
all: clean $(BINS) $(PYMOD_OUT)


	
# TESTS
TEST = 	./test.sh \
			"mpiexec -n 8 ./$(BINS)" \
			tests/rudy/g05_60.0 \
			tests/rudy/g05_60.0-expected_output \
			params ;\
	done


TEST_ALL_60 = 	for i in $(shell seq 0 9); do \
			./test.sh \
			"mpiexec -n 8 ./$(BINS)" \
			tests/rudy/g05_60.$$i \
			tests/rudy/g05_60.$$i-expected_output \
			params ;\
	done

TEST_ALL_80 = 	for i in $(shell seq 0 9); do \
			./test.sh \
			"mpiexec -n 8 ./$(BINS)" \
			tests/rudy/g05_80.$$i \
			tests/rudy/g05_80.$$i-expected_output \
			params ;\
	done

TEST_ALL_100 = 	for i in $(shell seq 0 9); do \
			./test.sh \
			"mpiexec -n 8 ./$(BINS)" \
			tests/rudy/g05_100.$$i \
			tests/rudy/g05_100.$$i-expected_output \
			params ;\
	done

# TESTS
TEST_PYTHON = ./test.sh \
			"mpiexec -n 3 python3 biqbin.py" \
			tests/rudy/g05_60.0 \
			tests/rudy/g05_60.0-expected_output \
			params \

TEST_PYTHON_HEURISTIC = ./test.sh \
			"mpiexec -n 3 python3 biqbin.py" \
			tests/rudy/g05_60.0 \
			tests/rudy/g05_60.0-expected_output \
			params \

TEST_ALL_60_PYTHON = 	for i in $(shell seq 0 9); do \
			./test.sh \
			"mpiexec -n 8 python3 biqbin.py" \
			tests/rudy/g05_60.$$i \
			tests/rudy/g05_60.$$i-expected_output \
			params ;\
	done

TEST_ALL_80_PYTHON = 	for i in $(shell seq 0 9); do \
			./test.sh \
			"mpiexec -n 8 python3 biqbin.py" \
			tests/rudy/g05_80.$$i \
			tests/rudy/g05_80.$$i-expected_output \
			params ;\
	done

TEST_ALL_100_PYTHON = 	for i in $(shell seq 0 9); do \
			./test.sh \
			"mpiexec -n 8 python3 biqbin.py" \
			tests/rudy/g05_100.$$i \
			tests/rudy/g05_100.$$i-expected_output \
			params ;\
	done

# test-qubos:
# 	@for file in tests/qubo/*.pkl; do \
# 		mpiexec -n 1 python3 run_qubo_tests.py $$file params || exit $$?; \
# 	done

clean-output:
	rm -f rudy/*.output*
	rm -f tests/rudy/*.output*
	rm -f tests/qubo/*.output*

# Clean rule #
clean: clean-output
	rm -rf $(BINS) $(OBJS) $(PYMOD_OUT)


# Rules for binaries
$(BINS): $(OBJS)
	$(CPP) -o $@ $^ $(INCLUDES) $(LIB) $(CFLAGS) $(LINALG)

# BiqBin code rules
$(OBJ)/%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/%.o: %.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Python module rule
$(PYMOD_OUT): $(OBJS)
	$(CPP) -o $@ $^ -shared -fPIC $(INCLUDES) $(LIB) $(CFLAGS) $(LINALG) -Wl,--no-undefined





test-all: clean-output
	n/a
	$(TEST_ALL_60)
	$(TEST_ALL_80)
	$(TEST_ALL_100)


test: clean-output
	n/a
	$(TEST)

test-python: clean-output
	$(TEST_PYTHON)

test-all-python: clean-output
	$(TEST_ALL_60_PYTHON)
	$(TEST_ALL_80_PYTHON)
	$(TEST_ALL_100_PYTHON)

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
	
