#include "filemanager.h"
#include <iostream>

using std::string;
using std::fstream;
using std::map;
using std::cout;
using std::endl;

FileManager::FileManager()
{

}

FileManager::~FileManager()
{

}

fstream &FileManager::openFile(string filename, string key, std::ios_base::openmode mode)
{
    if(m_openFiles.find(key) == m_openFiles.end()) {
        // File is not already open
        m_openFiles[key] = new fstream(filename, mode);
    }

    if(m_openFiles[key]->is_open()) {
        cout << "Error, file " << filename << " could not be opened. Aborting!" << endl;
        exit(1);
    }

    return *m_openFiles[key];
}

