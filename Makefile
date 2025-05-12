CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native -fopenmp
LDFLAGS = -fopenmp

# Eigen paths - try different standard locations
EIGEN_FLAGS = -I/usr/include/eigen3 -I/usr/local/include/eigen3 -I$(EIGEN_DIR) -I./eigen

SRCS = gauss.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)
MAIN = main

TEST_SRCS = gauss.cpp gauss_test.cpp
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
TEST = tests

.PHONY: all clean test

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(EIGEN_FLAGS) -o $@ $^ $(LDFLAGS)

tests: $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(EIGEN_FLAGS) -o $@ $^ $(LDFLAGS) -lgtest -pthread

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN_FLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TEST_OBJS) $(MAIN) $(TEST)

test: tests
	./$(TEST) 