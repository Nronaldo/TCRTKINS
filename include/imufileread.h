#ifndef _IMUFILE_READ_H_
#define _IMUFILE_READ_H_
#include <string>
#include <iostream>
#include <stdio.h>
#include "Time_utility.h"
#include <fstream>
#include "Load_Configuration.h"
using namespace std;
int find_IMU_firstepoch(ifstream &imufile, string & buffIMU, const gtime_t &headtime);
void Time_IMU_Extract(const string & buffIMU, vector<double> & dwdv, double &cur_time);

#endif