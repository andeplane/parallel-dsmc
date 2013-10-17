#pragma once
#include <string>
using std::string;

class ProgressBar
{
public:
    ProgressBar(int final_value, string description_);
    int old_progress;
    int final_value;
    int current_value;
    string description;

    void update(int new_value);
};
