#Compiler and flagsi
CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pthread -fopenmp
INCLUDES := -Iinclude

# Find all source files in the src/ directory
SRCS     := $(wildcard src/*.cpp)
OBJS     := $(SRCS:.cpp=.o)

# Final executable name
TARGET   := styx_model

# Default target
all: $(TARGET)

# Link object files into the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

# Compile each source file into an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean

