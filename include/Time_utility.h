# ifndef TIME_UTILITY_H_
# define TIME_UTILITY_H_
#include <stdio.h>
#include <math.h>

#include <iomanip>
#include <string>
#include <sstream>
#include <iostream>
using namespace std;
static const double gpst0[] = { 1980,1, 6,0,0,0 }; /* gps time reference */
static const double bdt0[] = { 2006,1, 1,0,0,0 }; /* beidou time reference */

class gtime_t {        /* time struct */
public:
	gtime_t() {}
	gtime_t(int t, double s) :time(t), sec(s) {}
	gtime_t(double t_s) :time(floor(t_s)), sec(t_s-floor(t_s)) {}
	int time;        /* time (s) expressed by standard time_t */
	double sec;         /* fraction of second under 1 s */
} ; //使用gtime_t的原因是用两个64位数据存取，保证了时间的精度！流下了不学无术的泪水

gtime_t epoch2time(const double *ep);
gtime_t gpst2time(int week,double sec);
gtime_t bdt2time(int week, double sec);
gtime_t bdt2gpst(gtime_t BDST);
void time2str(gtime_t time, string &str,int n);
gtime_t adjweek(gtime_t t, gtime_t t0);
gtime_t string2time(string str);
gtime_t timeadd(gtime_t t, double sec);
double timediff(gtime_t t1, gtime_t t2);
double time2gpst(gtime_t t, int *week);
double time2doy(gtime_t t);
void time2epoch(gtime_t t, double *ep);
double calculate_rcv_time(string bufferR);
#endif // !TIME_UTILITY_H_

#pragma once
