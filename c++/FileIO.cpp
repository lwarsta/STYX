#include "FileIO.h"

void FileIO::save_grid(std::string filePath, std::string data)
{
    std::ofstream outFile;
    outFile.open(filePath.c_str());
    outFile << data;
    outFile.close();
}

void FileIO::save_results(std::string filePath,
    std::vector < std::vector<std::string> > data)
{
    std::string csvData = "";

    for (size_t i = 0; i < data.size(); i++)
    {
        if (i > 0)
        {
            csvData += "\n";
        }

        for (size_t j = 0; j < data.at(i).size(); j++)
        {
            csvData += data.at(i).at(j);

            if (j < data.at(i).size() - 1)
            {
                csvData += ","; // " " -> ","
            }
        }
    }

    std::ofstream resultsFile;
    resultsFile.open(filePath);
    resultsFile << csvData;
    resultsFile.close();
}

std::vector < std::vector<std::string> > FileIO::load_and_tokenize_file(
    std::string filePath, char delimiter)
{
    std::ifstream fileStream(filePath);
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    fileStream.close();
    std::string line;
    std::vector < std::vector<std::string> > tokens;

    while (std::getline(buffer, line, '\n'))
    {
        std::stringstream ss(line);
        std::vector<std::string> words;
        std::string tmp;

        while (std::getline(ss, tmp, delimiter))
        { // ' '
            words.push_back(tmp);
        }

        tokens.push_back(words);
    }
    return tokens;
}

std::string FileIO::load_file(std::string filePath)
{
    std::ifstream fileStream(filePath);
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    fileStream.close();
    return buffer.str();
}
