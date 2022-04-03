#include "Rinex.h"
double PI = 3.1415926535897932;
double DTTOL = 0.025;
double EPS = 2.220446049e-16F;
double INF = 3.402823466e+30F;
Conf_struct::Conf_struct() {
	// GNSS
	dynamics_on = 1;
	arfilter = 0;
	outsolfile_name = "1";
	navsys = "GC";
	output_pos_format = "xyz";
	nfreq = 1;
	tint = 1;
	plot_mode = 0;
	armode = 0;
	maxinno = 30.0;
	maxgdop = 30.0;
	exec_mode = 0;
	Ls_res_threshold = 100.0;
	elmin_mask = 15 * PI / 180;
	sclkstab = 5E-12;
	VAR_POS = pow(30, 2);
	VAR_VEL = pow(10, 2);
	VAR_ACC = pow(10, 2);
	varholdamb = 1;
	maxout = 3;
	minlock = 5;
	minfix = 5;
	optmindropsats = 10;
	optminholdsats = 5;
	optminfixsats = 2;
	prn = { 1E-4,1E-3,1E-4,0.1,0.01,0.0 };
	Error_std = { 100,0.003,0.003,0,1 };
	Ini_state_std = { 30.0, 0.03, 0.3 };
	thresar = { 3.0,0.999 };
	eratio = { 100.0,100.0,100.0 };
	basepose = { -2168841.004 ,4386629.315 ,4077157.089 };// The lab/
	// IMU
	ins_mode =0;
	ins_aid_CS = 0;
	ins_aid_IGG = 0;
	sample_rate = 50;
	psd_gyro = 5e-9;
	psd_acce = 2e-7;
	psd_bg = 2e-10;
	psd_ba = 1e-7;
	lever = {0,0,0};
	init_att_unc = {0.1,0.1,0.1};
	init_vel_unc = {1,1,1};
	init_pos_unc = {10,10,10};
	init_bg_unc = { 4e-5,4e-5 ,4e-5 };
	init_ba_unc = {9.8e-3,9.8e-3,9.8e-3};
	MAX_GNSS_OUTAGE = 20;
	IGGC0 = 2;
	IGGC1 = 5;
	/// out
	timef = 1;
	pos_mode = 0;
	att_mode = 0;
	out_mode = { 1,0,1,0,1,0,0,0 };//print or not [pos, posvar,vel,velvar,att,attvar,bg,ba]
}
Conf_struct TC_conf;
void NAV_RINEX::display() {}

int NAV_RINEX::readfile() {
	ifstream navfile(filename);
	if (!navfile.is_open()) {
		cout << "Fail to open the navfile" << endl;
		return 0;
	}
	string buff;
	int header_flag = 1;
	stringstream ss;
	while (getline(navfile, buff)) {
		if (header_flag == 1) {
			if (buff.substr(60, 13) == "END OF HEADER") {
				header_flag = 0;
				continue;
			}
			if (buff.substr(60, 12) == "LEAP SECONDS") {
				ss << buff.substr(0, 6);
				ss >> LEAP_SECONDS;
				ss.clear();
				ss.str("");
				continue;
			}
			if(buff.substr(60, 16) == "IONOSPHERIC CORR") {
				CORR_TYPE.push_back(buff.substr(0, 4));
				int i = 0;
				int j = 5;
				for (i = 0; i < 4; i++, j += 12) {
					double temp;
					string temps = buff.substr(j, 13);
					string::iterator it = find(temps.begin(), temps.end(), 'D');
					if (it != temps.end())
						*it = 'E';
					ss << temps;
					ss >> temp;
					ION_CORR.push_back(temp);
					ss.clear();
					ss.str("");
				}				
			}
			if (buff.substr(60, 16) == "TIME SYSTEM CORR") {
				CORR_TYPE.push_back(buff.substr(0, 4));
				double temp;
				string temps;
				temps = buff.substr(5, 17);
				string::iterator it = find(temps.begin(), temps.end(), 'D');
				if (it != temps.end())
					*it = 'E';
				ss << temps;
				ss >> temp;
				TIME_CORR.push_back(temp);
				ss.clear();
				ss.str("");
				temps = buff.substr(22, 16);
				it = find(temps.begin(), temps.end(), 'D');
				if (it != temps.end())
					*it = 'E';
				ss << temps;
				ss >> temp;
				TIME_CORR.push_back(temp);
				ss.clear();
				ss.str("");
				ss << buff.substr(38, 7);
				ss >> temp;
				TIME_CORR.push_back(temp);
				ss.clear();
				ss.str("");
				ss << buff.substr(45, 5);
				ss >> temp;
				TIME_CORR.push_back(temp);
				ss.clear();
				ss.str("");	
			}
		}
		// 读body
		else {
			if (buff[0] == 'G') {
				solve_epoch(buff, navfile, GPS_ID, GPS_NAV_DATA, GPS_TOC);
			}
			else if(buff[0] == 'C'){
				solve_epoch(buff, navfile, BDS_ID, BDS_NAV_DATA, BDS_TOC);
			}
		}
	}
	navfile.close();
	return 1;
}
void NAV_RINEX::solve_epoch(string buff,ifstream &navfile,vector<int> & GNSS_ID, vector<double> & GNSS_NAV_DATA, vector<int> & GNSS_TOC) {
	int i = 1;
	int temp;
	double tempd;
	stringstream ss;
	ss << buff.substr(1, 2);
	ss >> temp;
	ss.clear();
	ss.str("");
	GNSS_ID.push_back(temp);
	gtime_t toc = string2time(buff.substr(4, 19));
	GNSS_TOC.push_back(toc.time);
	for (int k = 0, j = 23; k < 3; k++, j += 19) {
		ss << buff.substr(j, 19);
		ss >> tempd;
		GNSS_NAV_DATA.push_back(tempd);
		ss.clear();
		ss.str("");
	}
	for (i = 1; i < 8;i++) {
		getline(navfile, buff);
		if (i < 7) {
			for (int k = 0, j = 4; k < 4; k++, j += 19) {
				ss << buff.substr(j, 19);
				ss >> tempd;
				GNSS_NAV_DATA.push_back(tempd);
				ss.clear();
				ss.str("");
			}
		}
		else {
			for (int k = 0, j = 4; k < 2; k++, j += 19) {
				ss << buff.substr(j, 19);
				ss >> tempd;
				GNSS_NAV_DATA.push_back(tempd);
				ss.clear();
				ss.str("");
			}
		}
	}
}
void NAV_RINEX::NAV_RINEX_CLEAR() {
	vector<int>().swap(GPS_ID);
	vector<int>().swap(BDS_ID);
	vector<int>().swap(GPS_TOC);
	vector<int>().swap(BDS_TOC);
	vector<string >().swap(CORR_TYPE);
	vector<double >().swap(ION_CORR);
	vector<double >().swap(TIME_CORR);
	vector<double>().swap(GPS_NAV_DATA);
	vector<double>().swap(BDS_NAV_DATA);
}
gtime_t  OBS_RINEX::readheader() {
	//ifstream obsfileB(filenameBase);
	//ifstream obsfileR(filenameRover);
	gtime_t baseheadtime = obtainheadparas(filenameBase,'B');
	gtime_t roverheadtime = obtainheadparas(filenameRover,'R');// 一般来说base 的时间不可能比rover的时间晚
	if (baseheadtime.time == 0 || roverheadtime.time == 0)
		return gtime_t(0, 0);
	return timediff(baseheadtime ,roverheadtime)<0 ? roverheadtime : baseheadtime;
}
gtime_t OBS_RINEX::obtainheadparas(string file, char mode) {
	ifstream obsfile(file);
	if (!obsfile.is_open()) {
		cout << "obsfile 没打开" << endl;
		return gtime_t(0,0.0);
	}
	int GPS_Fnum = 0;
	int BDS_Fnum = 0;
	vector<double > GPS_fre;
	vector<double > BDS_fre;
	string buff;
	stringstream ss;
	int tempi;
	double tempd;
	gtime_t headtime(0,0);
	gtime_t headtime_result(0, 0);
	//int timeflag = 0;
	while (getline(obsfile, buff)) {
		if (buff[0] == '>') {
			
			gtime_t previousheadtime = headtime;
			headtime = string2time(buff.substr(2, 27));
			double timediffv = 0;
			if (previousheadtime.time) {
				timediffv = timediff(headtime, previousheadtime);
				if (timediffv > TC_conf.tint*1.3) {//手机有的时候会给出xxx.7299这样的时间，此处乘个系数是防止对这种文件算不出headtime
					cout << mode << " file does not fit in your tint" << endl;
					return gtime_t(0, 0);
				}
			}
			if (headtime_result.time) {
				if (mode == 'B')
					base_interval = timediffv;
				if (mode == 'R')
					rover_interval = timediffv;
				break;
			}
			//timeflag++;
			//double x = fmod(time2gpst(headtime, NULL) + DTTOL, tint);
			if (fmod(time2gpst(headtime, NULL) + DTTOL, TC_conf.tint) <= DTTOL * 2.0) { //找到的第一个时间必须是足够接近tint的时间才结束
				headtime_result = headtime;
				continue;
			}
			else
				getline(obsfile, buff);
		}
		if (buff.substr(60, 19) == "SYS / # / OBS TYPES") {
			ss << buff.substr(4, 2);
			ss >> tempi;
			ss.clear();
			ss.str("");
			if (buff[0] == 'G') { // RTKLIB 面对C2P C2X 两个，就会选择靠后的那个频点,两个信号只算成是一个频率
				GPS_Fnum = tempi / 4;
				int C2_flag = 0;
				int duplicate_num = 0;
				for (int i = 0; i < GPS_Fnum; i++) {
					if (buff.substr(0 + i * 16 + 7, 3) == "C1C")
						GPS_fre.push_back(1575.42);
					else if (buff.substr(0 + i * 16 + 7, 2) == "C2") {
						if (C2_flag == 0) {
							GPS_fre.push_back(1227.60);
							C2_flag = 1;
							duplicate_num++;
						}
					}	
					else if (buff.substr(0 + i * 16 + 7, 3) == "C5Q" || buff.substr(0 + i * 16 + 7, 3) == "C5X")
						GPS_fre.push_back(1176.45);
				}
				if (C2_flag) {
					GPS_L2duplicate_flag = 1;
					GPS_Fnum -= duplicate_num;
				}
			}
			else if (buff[0] == 'C') {
				BDS_Fnum = tempi / 4;
				for (int i = 0; i < BDS_Fnum; i++) {
					if (buff.substr(0 + i * 16 + 7, 3) == "C2I" || buff.substr(0 + i * 16 + 7, 3) == "C1I")
						BDS_fre.push_back(1561.098);
					else if (buff.substr(0 + i * 16 + 7, 3) == "C7I")
						BDS_fre.push_back(1207.14);
					else if (buff.substr(0 + i * 16 + 7, 3) == "C6I")
						BDS_fre.push_back(1268.52);
				}
			}
		}
	}
	if (mode == 'B') {
		GPS_Fnum_B = GPS_Fnum;
		BDS_Fnum_B = BDS_Fnum;
		GPS_fre_B = GPS_fre;
		BDS_fre_B = BDS_fre;
	}
	else if (mode == 'R') {
		GPS_Fnum_R = GPS_Fnum;
		BDS_Fnum_R = BDS_Fnum;
		GPS_fre_R = GPS_fre;
		BDS_fre_R = BDS_fre;
	}
	else
		cout << "mode is wrong" << endl;
	
	obsfile.close();
	return headtime_result;
}

void OBS_RINEX::findfirstepoch(ifstream &fileB, ifstream &fileR,string & buffB,string &buffR,gtime_t headtime) {
	//外部 readheader 已经确定了接近tint的headtime。0.05.所以此处只需要找到对应base 时间即可
	// 虽然对base 和 rover都进行了寻找，但其实base 的时间不会比rover晚。
	while (getline(fileB, buffB)) {
		if (buffB.substr(60, 13) == "END OF HEADER") {
			break;
		}
	}
	while (getline(fileB, buffB)) {
		if (buffB[0] == '>')
			if (fabs(timediff(headtime,string2time(buffB.substr(2, 27)))) > 2*DTTOL)
				continue;
			else
				break;
	}
	while (getline(fileR, buffR)) {
		if (buffR.substr(60, 13) == "END OF HEADER") {
			break;
		}
	}
	
	while (getline(fileR, buffR)) {
		if (buffR[0] == '>')
			if (fabs(timediff(headtime, string2time(buffR.substr(2, 27)))) > 2 * DTTOL)
				continue;
			else
				break;
	}
	return;
	
}

void OBS_RINEX::readepoch(ifstream &fileB, ifstream &fileR, string buffB, string buffR,int &Bupdateflag,int &Rupdateflag, int &Btintflag, int &Rtintflag) {
	//能进入到此函数，一定是buffB 和 buffR 对应的时间相同
	//并且base rover file中都能找到符合Tint的数据，rover第一组符合tint的数据 对应buffB 和 buffR.
	//正常情况不该出现buffR比buffB 时间提前的情况
	// 但是有的时候rover 一个星也收不到，接收机就不提供当前时刻的时间戳。导致buffR 读取下一历元时间超前了
	// 为了解决这个问题，假设base处于良好opensky ，不至于出现一个星都没有的环境下。
	// 那么只可能是rover 出现这种情况，且rover出现这种情况后，如果时间超前了，
	stringstream ss;
	string temps;
	int flagB = 0;
	int flagR = 0;
	int numB = 0;
	int numR = 0;
	Bupdateflag = 0;
	Rupdateflag = 0;
	gtime_t timestampB, timestampR;

		timestampB = string2time(buffB.substr(2, 27));
		flagB = buffB[31] - '0';
		timestampR = string2time(buffR.substr(2, 27));
		flagR = buffR[31]-'0';
		ss << buffB.substr(32, 3);
		ss >> numB;
		ss.clear();
		ss.str("");
		ss << buffR.substr(32, 3);
		ss >> numR;
		ss.clear();
		ss.str("");
		//先得找到比较接近的rover 和 base 的第一个历元
		// 其实在外部已经找到了，但是手机给出的数据有的时候timeB和 time R差的有点多，所以这里再选择一下。
		//此处应该已经在外部readheader中绕过了了手机的数据可能存在的起始几个历元偏离整s的问题，应该不再需要在此处对时间进行判断。
		/*if (fabs(timediff(timestampB, timestampR)) <= 0.1) {*/ 
		// 两者时间相同，应该都读取 或者 两者时间不同，本站时间不满足tint约束，本站读取
		// base 间隔比 rover 间隔大，rover时间反而超前了，base 需要读取
		if ((fabs(timediff(timestampB, timestampR)) <= 2 * DTTOL) || 
			(fabs(timediff(timestampB, timestampR)) > 2 * DTTOL && fmod(time2gpst(timestampB, NULL) + DTTOL, TC_conf.tint) > DTTOL * 2.0)||
			((base_interval>rover_interval+2*DTTOL)&&(timediff(timestampR, timestampB)>2*DTTOL))) {
			vector<int>().swap(GPS_satID_B);
			vector<int>().swap(BDS_satID_B);
			vector<int>().swap(GPS_flags_B);
			vector<int>().swap(BDS_flags_B);
			vector<double >().swap(GPS_measurements_B);
			vector<double >().swap(BDS_measurements_B);
			getepochdata(numB, fileB, 'B');
			cur_sat_num_B = numB;
			epochflag_B = flagB;
			GPStime_B = timestampB;
			if (fmod(time2gpst(timestampB, NULL) + DTTOL, TC_conf.tint) < DTTOL * 2.0) { //读取完以后时间满足了int指标，才应该置成1
				Btintflag = 1;
			}
			Bupdateflag = 1;
		}
		// rover 间隔比 base 间隔大，base时间反而超前了，rover 需要读取
		if ((fabs(timediff(timestampB, timestampR)) <= 2 * DTTOL) ||
			(fabs(timediff(timestampB, timestampR)) > 2 * DTTOL && fmod(time2gpst(timestampR, NULL) + DTTOL, TC_conf.tint) > DTTOL * 2.0)||
			((rover_interval > base_interval + 2 * DTTOL) && (timediff(timestampB, timestampR) > 2 * DTTOL))) {
			vector<int>().swap(GPS_satID_R);
			vector<int>().swap(BDS_satID_R);
			vector<int>().swap(GPS_flags_R);
			vector<int>().swap(BDS_flags_R);
			vector<double >().swap(GPS_measurements_R);
			vector<double >().swap(BDS_measurements_R);
			getepochdata(numR, fileR, 'R');
			cur_sat_num_R = numR;
			epochflag_R = flagR;
			GPStime_R = timestampR;
			if (fmod(time2gpst(timestampR, NULL) + DTTOL, TC_conf.tint) < DTTOL * 2.0) { //读取完以后时间满足了int指标，才应该置成1
				Rtintflag = 1;
			}
			Rupdateflag = 1;
		}
		return;

	



}
void OBS_RINEX::getepochdata(int num,ifstream &file, char mode) {
	vector<int> GPS_satID;
	vector<int> BDS_satID;
	vector<int> GPS_flags;
	vector<int> BDS_flags;
	vector<double > GPS_measurements;
	vector<double > BDS_measurements;
	int GPS_Fnum, BDS_Fnum;
	string buff;
	if (mode == 'B') {
		GPS_Fnum = GPS_Fnum_B;
		BDS_Fnum = BDS_Fnum_B;
	}
	else if (mode == 'R') {
		GPS_Fnum = GPS_Fnum_R;
		BDS_Fnum = BDS_Fnum_R;
	}
	else
		cout << "no such mode" << endl;
	while (num > 0) {
		getline(file, buff);
		num--;
		if (buff[0] == 'G') {
			getepochdatabody(GPS_Fnum, GPS_L2duplicate_flag,buff,GPS_satID,GPS_flags,GPS_measurements);
		}
		else if(buff[0] == 'C'){
			getepochdatabody(BDS_Fnum, 0,buff, BDS_satID, BDS_flags, BDS_measurements);
		}	
	}
	if (mode == 'B') {
		GPS_satID_B = GPS_satID;
		BDS_satID_B = BDS_satID;
		GPS_flags_B = GPS_flags;
		BDS_flags_B = BDS_flags;
		GPS_measurements_B = GPS_measurements;
		BDS_measurements_B = BDS_measurements;
	}
	else if (mode == 'R') {
		GPS_satID_R = GPS_satID;
		BDS_satID_R = BDS_satID;
		GPS_flags_R = GPS_flags;
		BDS_flags_R = BDS_flags;
		GPS_measurements_R = GPS_measurements;
		BDS_measurements_R = BDS_measurements;
	}
	else
		cout << "no such mode" << endl;
	return;
}
void OBS_RINEX::getepochdatabody(int Frequency_num, int duplicate_flag,string buff, vector<int>&GNSS_ID, vector<int>&GNSS_flags, vector<double>&GNSS_measuremenst) {
	stringstream ss;
	string temps;
	string subbuff;
	int tempint;
	double tempdouble;
	ss << buff.substr(1, 2);
	ss >> tempint;
	GNSS_ID.push_back(tempint);
	ss.clear();
	ss.str("");
	int L2_duplicate_offset = 0;//由于有两个L2信号，按照RTKLIB的做法，只要后一个的
	for (int i = 0; i < Frequency_num; i++) {
		if (duplicate_flag && i == 1) {
			L2_duplicate_offset = 64;
		}
		if (5 + i * 64+ L2_duplicate_offset > buff.size()) {
			GNSS_measuremenst.push_back(0);
			GNSS_measuremenst.push_back(0);
			GNSS_measuremenst.push_back(0);
			GNSS_measuremenst.push_back(0);
			GNSS_flags.push_back(0);
			GNSS_flags.push_back(0);
			GNSS_flags.push_back(0);
		}
		else {
			subbuff = buff.substr(5 + i * 64 + L2_duplicate_offset, 12);
			if (subbuff.substr(0, 4) == "    ") {
				GNSS_measuremenst.push_back(0);
				
			}
			else {
				ss << subbuff;
				ss >> tempdouble;
				GNSS_measuremenst.push_back(tempdouble);
				ss.clear(); ss.str("");
			}
			subbuff = buff.substr(20 + i * 64 + L2_duplicate_offset, 13);
			if (subbuff.substr(0, 4) == "    ") {
				GNSS_measuremenst.push_back(0);
			}
			else {
				ss << subbuff;
				ss >> tempdouble;
				GNSS_measuremenst.push_back(tempdouble);
				ss.clear(); ss.str("");
			}
			subbuff = buff.substr(40 + i * 64 + L2_duplicate_offset, 9);
			if (subbuff.substr(0, 4) == "    ") {
				GNSS_measuremenst.push_back(0);
			}
			else {
				ss << subbuff;
				ss >> tempdouble;
				GNSS_measuremenst.push_back(tempdouble);
				ss.clear(); ss.str("");
			}
			subbuff = buff.substr(59 + i * 64 + L2_duplicate_offset, 6);
			if (subbuff.substr(0, 4) == "    ") {
				GNSS_measuremenst.push_back(0);
			}
			else {
				ss << subbuff;
				ss >> tempdouble;
				GNSS_measuremenst.push_back(tempdouble);
				ss.clear(); ss.str("");
			}
			if (buff[18 + i * 64] != ' ')
				GNSS_flags.push_back(buff[18 + i * 64 + L2_duplicate_offset] - '0');
			else
				GNSS_flags.push_back(0);
			if (buff[34 + i * 64] != ' ')
				GNSS_flags.push_back(buff[34 + i * 64 + L2_duplicate_offset] - '0');
			else
				GNSS_flags.push_back(0);
			if (buff[33 + i * 64] != ' ')
				GNSS_flags.push_back(buff[33 + i * 64 + L2_duplicate_offset] - '0');
			else
				GNSS_flags.push_back(0);	
		}
	}
	return;
}
void OBS_RINEX::getdoppler() {
	for (int i = 0; i < GPS_satID_R.size(); i++) {
		int ID = GPS_satID_R[i];
		double doppler = GPS_measurements_R[i*GPS_Fnum_R * 4 + 2];
		if (sat_doppler.empty() || sat_doppler.find(ID) == sat_doppler.end()) {
			sat_doppler.insert(SDpair(ID, vector<double >(1, 0)));
			sat_doppler[ID][0] = doppler;
		}
		else {
			sat_doppler[ID].push_back(doppler);
		}
	
	}
	for (int i = 0; i < BDS_satID_R.size(); i++) {
		int ID = BDS_satID_R[i]+40;
		double doppler = BDS_measurements_R[i*BDS_Fnum_R * 4 + 2];
		if (sat_doppler.empty() || sat_doppler.find(ID) == sat_doppler.end()) {
			sat_doppler.insert(SDpair(ID, vector<double >(1, 0)));
			sat_doppler[ID][0] = doppler;
		}
		else {
			sat_doppler[ID].push_back(doppler);
		}
	}
}
void OBS_RINEX::printfdoppler() {
	SDmp::iterator it = sat_doppler.begin();
	string RootDir = "D:\\SHMILY\\340_RTK\\RTK_CPP\\doppler\\a";
	for (; it != sat_doppler.end(); it++) {
		int sat = it->first;
		if (sat > 40) {
			RootDir[RootDir.size()-1]='C';
			sat -= 40;
		}
		else
			RootDir[RootDir.size() - 1] = 'G';
		vector<double > doppler_v = it->second;
		string sats = to_string(sat)+".txt";
		string filepath = RootDir + sats;
		ofstream of;
		of.open(filepath, ios::out | ios::app);
		if (!of) {
			cout << "cannot open" << endl;
			return;
		}
		for (int i = 0; i < doppler_v.size(); i++) {
			of << setprecision(7)<<doppler_v[i] << endl;
		}
		//break;
	}
	return;



}