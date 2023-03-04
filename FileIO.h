#ifndef _FILEIO_H
#define _FILEIO_H
#include <fstream>
#include <sstream>
#include <vector>

class FileIO {
private:
public:
    std::string load_file(std::string filePath);
    std::vector < std::vector<std::string> > load_and_tokenize_file(
        std::string filePath, char delimiter);
    void save_grid(std::string filePath, std::string data);
    void save_results(std::string filePath, 
        std::vector < std::vector<std::string> > data);
};

#endif
