
CXX = g++
CXXFLAGS = -std=c++20  -O3  -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT -I./

TARGET = surgery

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 	

clean:
	rm -f $(OBJS) $(TARGET)
