#include <progressbar.h>
#include <iomanip>>
#include <cmath>
#include <iostream>
using std::cout;
using std::setfill;
using std::setw;
using std::endl;
ProgressBar::ProgressBar(int final_value_, string description_)
{
    old_progress = -1;
    current_value = 0;
    final_value = final_value_;
    description = description_;
}

void ProgressBar::update(int new_value) {
    current_value = new_value;
    int progress = floor(current_value/double(final_value) / 0.1);
    if (progress > old_progress) {
        cout << description;
        cout << " [";
        cout << setfill('#') << setw(progress) << "";
        cout << setfill('-') << setw(10-progress) << "";
        cout << "]" << endl;
        old_progress = progress;
    }
}
