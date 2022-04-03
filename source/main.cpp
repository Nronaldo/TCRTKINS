#include "PositionMode.h"



int main(int argc,char ** argv){ 
	NAV_RINEX nav;
	OBS_RINEX obs;
	ifstream obsfileB;
	ifstream obsfileR;
	string imufilename;
	Load_Config_Params(nav, obs, obsfileB, obsfileR, imufilename);
	//TC_conf.ins_mode = 0;
	if (!TC_conf.ins_mode) {
		Position_RTK(nav, obs, obsfileB, obsfileR);
	}
	else {
		Position_TC(nav, obs, obsfileB, obsfileR, imufilename);
	}
	return 1;
}