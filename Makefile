# --- Compiler and Linker Settings ---
CXX = g++
CXXFLAGS = -Wall -std=c++17 -O2
LDFLAGS = -lcfitsio
TARGET = app

# --- Directory and Source Definitions ---
SRC_DIR = src
DATA_SRC = $(SRC_DIR)/fits/fits_io.cpp
CATALOG_SRC = $(SRC_DIR)/catalog/catalog.cpp
TRIANGLE_SRC = $(SRC_DIR)/triangle/triangle.cpp
TRIAD_SRC = $(SRC_DIR)/triad/triad.cpp
MAIN_SRC = $(SRC_DIR)/main.cpp

# List all your source files here. Add more as you create them (detector.cpp, solver.cpp)
SRCS = $(MAIN_SRC) $(DATA_SRC) $(CATALOG_SRC) $(TRIANGLE_SRC) $(TRIAD_SRC)
# Convert source file names (.cpp) to object file names (.o)
OBJS = $(SRCS:.cpp=.o)


# --- Primary Build Target ---
all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "Linking executable..."
	$(CXX) $(OBJS) $(LDFLAGS) -o $(TARGET)
	@echo "--- Build successful: $(TARGET) ---"

# --- Rule to compile .cpp files into .o files ---
# This generic rule tells make how to handle any C++ source file.
# $< is the first prerequisite (the .cpp file)
# $@ is the target (the .o file)
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --- Phony Targets ---
.PHONY: all clean

clean:
	@echo "Cleaning up object files and executable..."
	rm -f $(OBJS) $(TARGET)