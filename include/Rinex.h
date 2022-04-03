#ifndef RINEX_H_
#define RINEX_H_



#include <stdio.h>
#include <fstream>
//#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <map>
#include "Time_utility.h"
extern double PI;
extern double DTTOL;
extern double EPS; 
extern double INF;
typedef std::pair<int, vector<double>> SDpair;
typedef std::map <int, vector<double> > SDmp;//存放某卫星ID 和 对应的多普勒

class Conf_struct {
public:
	Conf_struct();
	// GNSS
	int dynamics_on;
	int arfilter;
	std::string outsolfile_name;
	std::vector<double > basepose;
	std::string navsys;
	std::string output_pos_format;
	double elmin_mask;
	double maxinno;
	double maxgdop;
	int armode;
	double sclkstab;
	double varholdamb;
	double Ls_res_threshold;
	double VAR_POS;
	double VAR_VEL;
	double VAR_ACC;
	double tint;
	short maxout;
	short minlock;
	short minfix;
	short optmindropsats;
	short optminholdsats;
	short optminfixsats;
	int exec_mode;
	int nfreq;
	short plot_mode;
	std::vector<double > Error_std;
	std::vector<double > prn;
	std::vector<double > Ini_state_std;
	std::vector<double > thresar;
	std::vector<double > eratio;
	//// INS
	unsigned char ins_mode;
	unsigned char ins_aid_CS;
	unsigned char ins_aid_IGG;
	float sample_rate;
	double psd_gyro;
	double psd_acce;
	double psd_bg;
	double psd_ba;
	std::vector<double> lever;
	std::vector<double> init_att_unc;
	std::vector<double> init_vel_unc;
	std::vector<double> init_pos_unc;
	std::vector<double> init_bg_unc;
	std::vector<double> init_ba_unc;
	short MAX_GNSS_OUTAGE;
	int IGGC0;
	int IGGC1;
	/// out 
	int timef;
	int pos_mode;
	int att_mode; //0: rad  1: deg 
	std::vector<int > out_mode;
};

extern Conf_struct TC_conf;

class NAV_RINEX {
public:
	NAV_RINEX() {};
	NAV_RINEX(string file) { filename = file; }
	~NAV_RINEX() {};
	void NAV_RINEX_CLEAR();
	void setfilename(string file) { filename = file; };
	void display();
	int readfile();
	void solve_epoch(string buff, ifstream &navfile, vector<int> & GNSS_ID, vector<double> & GNSS_NAV_DATA, vector<int> & GNSS_TOC);
	string filename;
	int LEAP_SECONDS;
	vector<string > CORR_TYPE;
	vector<double > ION_CORR;
	vector<double > TIME_CORR; 
	vector<int > GPS_ID;
	vector<int > BDS_ID;
	vector<int > GPS_TOC;
	vector<int > BDS_TOC;
	vector<double> GPS_NAV_DATA;
	vector<double> BDS_NAV_DATA;
	
};

class OBS_RINEX {
public:
	OBS_RINEX() { GPS_L2duplicate_flag = 0; base_interval = 0; rover_interval = 0; };
	OBS_RINEX(string fileB, string fileR) { filenameBase = fileB; filenameRover = fileR; GPS_L2duplicate_flag = 0; base_interval = 0; rover_interval = 0;}
	void setfilesnames(string fileB, string fileR) { filenameBase = fileB; filenameRover = fileR; };
	gtime_t readheader();
	gtime_t obtainheadparas(string file,char mode);
	void findfirstepoch(ifstream &fileB, ifstream &fileR,string &buffB,string &buffR,gtime_t headtime);
	void getepochdata(int,ifstream &fileB,char mode);
	void getepochdatabody(int Frequency_num,int duplicate_flag,string buff, vector<int>&, vector<int>&, vector<double>&);
	void readepoch(ifstream &fileB, ifstream &fileR,string buffB,string buffR, int &Bupdateflag, int &Rupdateflag, int &Btintflag, int &Rtintflag);
	/*与获得观测量有关的函数*/
	void getdoppler();
	void printfdoppler();
	string filenameBase,filenameRover;
	int GPS_Fnum_B;
	int BDS_Fnum_B;
	int GPS_Fnum_R;
	int BDS_Fnum_R;
	int GPS_L2duplicate_flag;
	vector<double> GPS_fre_B;
	vector<double> BDS_fre_B;
	vector<double> GPS_fre_R;
	vector<double> BDS_fre_R;
	int epochflag_B;
	int epochflag_R;
	int cur_sat_num_B;
	int cur_sat_num_R;
	vector<int > GPS_satID_B;
	vector<int > BDS_satID_B;
	vector<int > GPS_flags_B;
	vector<int > BDS_flags_B;
	vector<double > GPS_measurements_B;
	vector<double > BDS_measurements_B;
	gtime_t GPStime_B;
	vector<int > GPS_satID_R;
	vector<int > BDS_satID_R;
	vector<int > GPS_flags_R;
	vector<int > BDS_flags_R;
	vector<double > GPS_measurements_R;
	vector<double > BDS_measurements_R;
	gtime_t GPStime_R;
	/*与提取观测量有关的结构*/
	SDmp sat_doppler;
	/*与文件时间间隔有关系的量*/
	double base_interval;
	double rover_interval;
};

#endif // !RINEX_H_
#pragma once
