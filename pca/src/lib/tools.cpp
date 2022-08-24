#include "tools.hpp"

// Initializes a Mersenne Twister pseudo-random generator
RandomNumber::RandomNumber() {
    std::random_device rd;
    gen = std::mt19937(rd());
}

// Generates a random integer residing between 0 and upper_bound
// Note that 0 and upper_bound are also included in the possible result
int RandomNumber::get(int upper_bound) {
    std::uniform_int_distribution<int> dis(0, upper_bound);
    return dis(gen);
}

SimpleClock::SimpleClock() {}

// Starts a clock
void SimpleClock::start() {
    start_time = std::chrono::high_resolution_clock::now();
    return;
}

// Stops a clock to get the elapsed time in seconds
double SimpleClock::stop() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    return time_span.count();
}

// Stops a clock to get the elapsed time in seconds
// Formatted print message is included
void SimpleClock::stopAndPrint(std::string message) {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << message << " " << time_span.count() << "s" << std::endl;
    return;
}

FileNameEditor::FileNameEditor() {}

// Extracts the name only excluding the directory path and extension
std::string FileNameEditor::extractName(std::string file_name) {
    std::string name = file_name;

    for(int i = name.length()-1; i >= 0; i--) {
        if(name[i] == '.') {
            name = name.substr(0, i);
            break;
        }
    }

    for(int i = name.length()-1; i >= 0; i--) {
        if(name[i] == '/' || name[i] == '\\') {
            name = name.substr(i+1, name.length());
            break;
        }
    }

    return name;
}

// Adds given directory path before the given file name
std::string FileNameEditor::addPath(std::string directory, std::string file_name) {
    if(directory.back() != '/' || directory.back() != '\\') {
        directory.append("/");
    }

    return directory + file_name;
}