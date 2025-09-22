CXX ?= g++
CXXFLAGS ?= -O3 -Wall
TARGET = BLOOM
SOURCES = mnode.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

clean:
	rm -f $(TARGET)

.PHONY: all clean
