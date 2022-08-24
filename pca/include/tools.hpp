#include <iostream>
#include <random>
#include <chrono>
#include <string>

#ifndef TOOLS_H
#define TOOLS_H

class RandomNumber {
public:
    RandomNumber();
    int get(int upper_bound);
    
private:
    std::mt19937 gen;
};

class SimpleClock {
public:
    SimpleClock();
    void start();
    double stop();
    void stopAndPrint(std::string message);

private:
    std::chrono::high_resolution_clock::time_point start_time;
};

class FileNameEditor {
public:
    FileNameEditor();
    std::string extractName(std::string file_name);
    std::string addPath(std::string directory, std::string file_name);
};

#endif