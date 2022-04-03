#ifndef _LOAD_CONFIGURATION_H_
#define	_LOAD_CONFIGURATION_H_
#include <string>
#include <iostream>
#include <fstream>
//#include <stdlib.h>
#include <sstream>
#include <io.h> // windows
#include <direct.h> // windows
// #include <unistd.h> // linux
#include <Rinex.h>
#include <vector>
#include <Eigen\Dense>
#include "matplotlibcpp.h"


int Load_Config_Params(NAV_RINEX & navfile, OBS_RINEX &Bobsfile, ifstream&obsfileB, ifstream& obsfileR,  string & imufile);
//void Position_RTKMode();
//void Position_TCMode();

#endif
