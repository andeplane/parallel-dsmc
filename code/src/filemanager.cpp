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
    for(auto it : m_openFiles) {
        fstream* file = (fstream*)it.second;
        if(file->is_open()) {
            file->close();
        } else {
            cout << "Warning, file with key " << it.first << " was not open" << endl;
        }
    }

    m_openFiles.clear();
}

void FileManager::closeFile(string key) {
    auto it = m_openFiles.find(key);
    if(it != m_openFiles.end()) {
        fstream *file = (fstream *)it->second;
        if(file->is_open()) {
            file->close();
        } else {
            cout << "Warning, tried to close an non-open file with key " << key << endl;
        }

    }
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

