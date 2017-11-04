CXX = g++
CXXFLAGS = -Wall -g
SRC= ./src/main.cpp
minhash: ; $(CXX) $(CXXFLAGS) $(SRC) -o ./bin/minhash

clean: ; rm -rfv ./bin/minhash* src/*.o src/*.gch