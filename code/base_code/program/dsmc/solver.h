#pragma once
#include <system.h>
#define VERSION "1.0.08"

class Settings;
class StatisticsSampler;

class Solver {
public:
	int num_processors;
	int myid;
	double t_start;
	Settings *settings;
	StatisticsSampler *sampler;
	System system;
	
	Solver(int num_processors, int myid);
	void step();
	void finalize();
};