#include "imufileread.h"
int find_IMU_firstepoch(ifstream &imufile, string & buffIMU, const gtime_t & headtime) {
	int flag = 0;
	while (getline(imufile, buffIMU)) {
		if (!flag) {
			flag++;
			continue;// the head of the file
		}
		double IMUtime = stod(buffIMU.substr(0, 17));
		double headtimed = (double)(headtime.time) + headtime.sec;
		if (fabs(IMUtime - headtimed) < 0.501 / TC_conf.sample_rate)
			return 1;
	}
	return 0;
}
void Time_IMU_Extract(const string & buffIMU, vector<double> & dwdv,double &cur_time) {
	stringstream iss(buffIMU);
	int i = 0;
	string temp;
	while (getline(iss, temp, ',')) {
		if (i == 0) {
			cur_time = stod(temp);
			i ++;
		}
		else {
			dwdv[i++-1] = stod(temp)*1.0 / TC_conf.sample_rate;
		}
	}
	return;
}
