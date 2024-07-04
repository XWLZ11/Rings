CC = g++
CFLAGS = -O3 -Wall -std=c++11 -fPIC `python3 -m pybind11 --includes` -Iinclude 
LDFLAGS = -shared

SUFFIX = $(shell python3-config --extension-suffix)

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = rings$(SUFFIX)

all: $(TARGET) mv_ln

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

mv_ln:
	mkdir -p shared 
	mv $(TARGET) shared/
	rm -rf shared/rings.so && ln -s shared/$(TARGET) shared/rings.so

clean:
	rm -f *.o $(TARGET) 
