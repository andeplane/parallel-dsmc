#pragma once
#include <vector>
#include <iostream>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;
using std::max;

template <typename T>
class StatisticalValue {
public:
    vector<T> current_value;
    vector<T> sum;
    vector<T> sum_squared;
    int number_of_samples;
    int number_of_bins;

    StatisticalValue(int number_of_elements = 1)
        : number_of_samples(0)
    {
        number_of_bins = number_of_elements;
        current_value.resize(number_of_elements,0);
        sum.resize(number_of_elements, 0);
        sum_squared.resize(number_of_elements, 0);
    }

    void add_value(T new_value) {
        if(sum.size() > 1) cout << "Warning, appending scalar to vector statistical value." << endl;
        current_value[0] = new_value;
        sum[0] += new_value;
        sum_squared[0] += new_value*new_value;

        number_of_samples++;
    }

    void add_value(vector<T> new_value) {
        current_value = new_value;

        // Assume std vector
        for(int i=0; i<new_value.size(); i++) {
            sum[i] += new_value[i];
            sum_squared[i] += new_value[i]*new_value[i];
        }

        number_of_samples++;
    }

    vector<T> get_current_value() { return current_value; }
    vector<T> get_average() {
        vector<T> average(sum.size());
        for(int i=0; i<sum.size(); i++) {
            average[i] = sum[i] / max(number_of_samples,1); // Divide by at least 1 to make sure we don't get any NaN values
        }
        return average;
    }

    vector<T> get_squared_average() {
        vector<T> squared_average(sum_squared.size());

        for(int i=0; i<sum_squared.size(); i++) {
            squared_average[i] = sum_squared[i] / max(number_of_samples,1); // Divide by at least 1 to make sure we don't get any NaN values
        }

        return squared_average;
    }

    vector<T> get_variance() {
        vector<T> average = get_average();
        vector<T> squared_average = get_squared_average();
        vector<T> variance(average.size());

        for(int i=0; i<average.size(); i++) {
            variance[i] = squared_average[i] - average[i]*average[i];
        }

        return variance;
    }

    vector<T> get_standard_deviation() {
        vector<T> variance = get_variance();
        vector<T> standard_deviation(variance.size());

        for(int i=0; i<variance.size(); i++) {
            standard_deviation[i] = sqrt(variance[i]);
        }

        return standard_deviation;
    }

    void resize(int number_of_bins_) {
        number_of_bins = number_of_bins_;
        current_value.resize(number_of_bins,0);
        sum.resize(number_of_bins,0);;
        sum_squared.resize(number_of_bins,0);;
    }
};
