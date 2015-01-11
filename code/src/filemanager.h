#pragma once
#include <fstream>
#include <string>
#include <map>

class FileManager
{
private:
    std::map<std::string, std::fstream*> m_openFiles;
public:
    FileManager();
    ~FileManager();
    std::fstream &openFile(std::string filename, std::string key, std::ios_base::openmode mode);
};
