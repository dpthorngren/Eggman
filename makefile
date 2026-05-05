# This file is not used for the normal eggman build.
# It is only here to provide a C++ library if desired.
# You may need to adjust CXX to whatever your compiler is:
CXX = clang++
CXXFLAGS = -Wall -g -I src/include -fPIC
LFLAGS = -shared -lgsl

TARGET = cpp_build/lib_eggman
SRCS=${wildcard src/cpp_source/*.cpp}
OBJS=${patsubst src/cpp_source/%.cpp,cpp_build/%.o,${SRCS}}

$(info Target: $(TARGET))
$(info Sources: $(SRCS))
$(info Objects: $(OBJS))

all: $(TARGET)

$(TARGET): $(OBJS) | cpp_build
	$(info Building executable $<.)
	$(CXX) $(LFLAGS) -o $(TARGET) $(OBJS)

cpp_build/%.o: src/cpp_source/%.cpp | cpp_build
	$(info Compiling $<.)
	$(CXX) $(CXXFLAGS) -c $< -o $@

cpp_build:
	mkdir -p cpp_build

clean:
	rm -rf cpp_build
