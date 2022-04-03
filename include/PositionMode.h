#ifndef _PositionMODE_H_
#define _PositionMODE_H_
#include "Ephemeris.h"
#include "Load_Configuration.h"
#include "Fusion.h"
#include "imufileread.h"
using namespace std;
void Position_RTK(NAV_RINEX &nav, OBS_RINEX &obs, ifstream &obsfileB, ifstream &obsfileR);
void Position_TC(NAV_RINEX &nav, OBS_RINEX &obs, ifstream &obsfileB, ifstream &obsfileR,const string & imufilename);


#endif