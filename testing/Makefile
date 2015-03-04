TARGET = grid
CXX = g++
CXXFLAGS = -O0 -Wall
LIB = -lm

COM = $(CXX) $(CXXFLAGS)

SRC = main.cpp energies.cpp climindex.cpp
OBJ = $(SRC:%.cpp=%.o)

$(TARGET) : $(OBJ)
		$(COM) -o $(TARGET) $(OBJ) $(LIB)

clean 	: rm -f $(TARGET) $(OBJ)

%.o		: %.cpp
		$(COM) -c $< -o $@

