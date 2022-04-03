#include "Receiver.h"

void receiver::GNSS_sever(ofstream &out) {
	gtime_t LASTsoltime;
	int LS_stat = 0;
	sol.set_state(0);
	LASTsoltime = getsoltime();
	LS_stat = LS_pos('N');//最小二乘定位
	settt(LASTsoltime);
	if (LS_stat) {
		LS_vel();//最小二乘定速，只有在最小二乘定位成功时才进行。
	}
	RTK_pos("PL");//不管最小二乘定位成不成功，都要进行RTK的更新
	if(!TC_conf.ins_mode)
		soltofile(LS_stat,out);
	return;
}
void receiver::Set_array_zero() {
	memset(RTK_x,0,nx*sizeof(double));
	memset(RTK_P, 0, nx *nx* sizeof(double));
	memset(outc_G, 0, NSATGPS * sizeof(short));
	memset(slip_G, 0, NSATGPS * sizeof(char));
	memset(half_G, 0, NSATGPS * sizeof(char));
	memset(lock_G, 0, NSATGPS * sizeof(short));
	memset(outc_B, 0, NSATCMP * sizeof(short));
	memset(slip_B, 0, NSATCMP * sizeof(char));
	memset(half_B, 0, NSATCMP * sizeof(char));
	memset(lock_B, 0, NSATCMP * sizeof(short));
	return;

}


void receiver::initializercv(NAV_RINEX &nav_rinex,ofstream & out) {
	//Pos_ecef = vector<double >(3, 0);
	//Pos_llh = vector<double >(3, 0);
	/*与最小二乘有关*/
	Ls_x = Eigen::VectorXd::Zero(5);//此处强行初始化成5维，扩展的时候其实可以根据nav的内容扩展成指定系统数目
	allsateph = All_Sats_EPH(nav_rinex);
	allsateph.satephsolve(nav_rinex);
	DOP = vector<double >(4, 0);
	Ls_res = 0;
	dynamics_on = TC_conf.dynamics_on;
	int pvanum = (dynamics_on ? 9 : 3);
	RTK_fixXpva = Eigen::VectorXd::Zero(pvanum);
	RTK_fixPpva = Eigen::MatrixXd::Zero(pvanum, pvanum);
	holdamb_flag = 0;
	arfilter = TC_conf.arfilter;
	nx = NSATCMP + NSATGPS + pvanum;
	com_bias = 0.0;
	RTK_x = new double[nx]{0};
	RTK_P = new double[nx*nx]{ 0 };//太大了
	outc_G = new short[NSATGPS]{0};
	outc_B = new short[NSATCMP]{0};
	slip_G = new char[NSATGPS]{0};
	half_G = new char[NSATGPS]{0};
	slip_B = new char[NSATCMP]{0};
	half_B = new char[NSATCMP]{0};
	lock_G = new short[NSATGPS] {0};
	lock_B = new short[NSATCMP] {0};
	//basepos.push_back(-2168841.004);//位置设置成实验室天线
	//basepos.push_back(4386629.315);
	//basepos.push_back(4077157.089);
	//basepos.push_back(-3961904.9228);//位置设置成东京数据基站
	//basepos.push_back(3348993.7822);
	//basepos.push_back(3698211.7504);
	basepos.push_back(TC_conf.basepose[0]);
	basepos.push_back(TC_conf.basepose[1]);
	basepos.push_back(TC_conf.basepose[2]);
	SATrefGPSID = 0;
	SATrefBDSID = 0;
	excsatindex = 0;
	nfix = 0;
	sol = solution(2);//强行设置成2系统 G+B
	string solfile_dir = "../solfile/";
	if (_access(solfile_dir.c_str(), 0) == -1)	//如果文件夹不存在
		_mkdir(solfile_dir.c_str());				//则创建
	string filename = solfile_dir + TC_conf.outsolfile_name + ".txt";
	output_pos_format = TC_conf.output_pos_format;
	plot_mode = TC_conf.plot_mode;
	out.open(filename, ios::out);
}
void receiver::updateobs(const OBS_RINEX & obsrinex) {
	if (obsrinex.GPS_satID_B.empty() && obsrinex.BDS_satID_B.empty())
		return;
	//本来我们在这里就进行了选择，选择共同卫星来，抛弃在Base不可见的rover 卫星
	// 但是RTKLIB 和 GINAV都不是这么做的，即使在BASE 不可见，这颗卫星仍然会得到保留以计算LS
	// 在rtkpos里才进行了选择。
	// GINAV 每次tdcp2vel,都利用上个历元的sol.rr和观测量，再算一遍上一历元的距离
	
	if (numiterator >= 2) {
		int aaa = 1;
	}
	distributesats(obsrinex.GPS_satID_B, obsrinex.GPS_satID_R, 'G', obsrinex);
	distributesats(obsrinex.BDS_satID_B, obsrinex.BDS_satID_R, 'B', obsrinex);
	
}
void receiver::distributesats(const vector<int> &GNSS_ID_B, const vector<int> &GNSS_ID_R,char mode, const OBS_RINEX & obsrinex) {
	// 现在来看，似乎不太好，应该是相同的卫星保留，直接更新，多余的卫星再删除
	//但是现在很懒不想这么改，为了载波相位LS定速，在这里打个补丁
	// GINAV 每次tdcp2vel,都利用上个历元的sol.rr和观测量，再算一遍上一历元的距离 只能在这里改，更新cp_t_tu_res
	//注意只更新rover的就可以
	map<int, double> Sat_Cpttures;
	map<int, gtime_t> Sat_gtime;
	
	map<int, OBS_SAT> &mp = mode == 'G'?GPSobs: BDSobs;
	map<int, OBS_SAT> &NoB_R_mp = mode == 'G' ? NoB_R_GPSobs: NoB_R_BDSobs;//choosecommonsats 删除了非BR都可见的星
	//为了测试tdcp2vel函数，人为保留下无B有R的星
		
	map<int, OBS_SAT>::iterator it = mp.begin();
	map<int, OBS_SAT>::iterator it_NoB_R = NoB_R_mp.begin();
	if (numiterator >= 2) {
		int aaa = 1;
	}
	for (; it!=mp.end(); it++) { // 注意，此时是上一历元,这个循环是为了根据上一历元的最终解算位置再算一下cp-r+delta_tu的r
		/*GINAV保存了old_obs 和 cur_obs,但是old_obs中有的星，上一历元可能没参与解算，比如高度角不够
			然而此处都直接参与解算了，没有排除，或许这是一个小问题。这里我们先保持一致*/
		/*if (rover.get_unusedflag() || rover.get_fdemask())
			continue;*/
		it->second.update_cp_t_tu_res(sol.getpos());
		Sat_Cpttures[it->first] = it->second.get_rover_cp_t_tu_res();
		Sat_gtime[it->first] = it->second.getBRtime();
	}
	for (; it_NoB_R != NoB_R_mp.end(); it_NoB_R++) {
		it_NoB_R->second.update_cp_t_tu_res(sol.getpos());
		Sat_Cpttures[it_NoB_R->first] = it_NoB_R->second.get_rover_cp_t_tu_res();
		Sat_gtime[it_NoB_R->first] = it_NoB_R->second.getBRtime();
	}
	mp.clear();//注意这里把上一历元的卫星全清空了
	NoB_R_mp.clear();
	for (int i = 0; i < GNSS_ID_B.size(); i++) {// 此时是当前历元
		vector<int> ::const_iterator itID = find(GNSS_ID_R.begin(), GNSS_ID_R.end(), GNSS_ID_B[i]);
		if(itID!=GNSS_ID_R.end()){ //共同可见星
			int indexR = distance( GNSS_ID_R.begin(),itID);
			OBS_SAT a = OBS_SAT(i,indexR, obsrinex,mode);
			if (Sat_Cpttures.count(GNSS_ID_B[i])) {
				a.set_rover_cp_t_tu_res(Sat_Cpttures[GNSS_ID_B[i]]);
				a.set_rover_last_time(Sat_gtime[GNSS_ID_B[i]]);
			}
				
			mp[GNSS_ID_B[i]] = a;
		}
		else { // base 的在rover中不可见 那么只构建OBS_SAT中的base
			int indexR = -1;
			OBS_SAT a = OBS_SAT(i, indexR, obsrinex, mode);
			mp[GNSS_ID_B[i]] = a;
		}
	}//然后再判断rover 可见 base 不可见的星
	for (int i = 0; i < GNSS_ID_R.size(); i++) {
		vector<int> ::const_iterator itID = find(GNSS_ID_B.begin(), GNSS_ID_B.end(), GNSS_ID_R[i]);
		if (itID == GNSS_ID_B.end()) { // 说明rover可见，base 不可见
			int indexB = -1;
			OBS_SAT a = OBS_SAT(indexB, i, obsrinex, mode);
			if (Sat_Cpttures.count(GNSS_ID_R[i])) {
				a.set_rover_cp_t_tu_res(Sat_Cpttures[GNSS_ID_R[i]]);
				a.set_rover_last_time(Sat_gtime[GNSS_ID_R[i]]);
			}
				
			mp[GNSS_ID_R[i]] = a;
		}
	}
}
void receiver::choosecommonsat() {
	
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	for (; itg != GPSobs.end(); ) {
		map<int, OBS_SAT>::iterator it = itg;
		itg++;
		if (numiterator >= debug_stop && numiterator <= 240 && it->first == 5) {
			GPSobs.erase(it);
			continue;
		}
		if (!it->second.get_roverflag()) {
			GPSobs.erase(it);
		}	
		else {
			if (!it->second.get_baseflag()) {
				NoB_R_GPSobs[it->first] = it->second;
				GPSobs.erase(it);
			}
				
		}
		/*if (numiterator >= debug_stop) {
			int aaaa = 1;
		}*/
		
	}
	for (; itb != BDSobs.end();) {
		map<int, OBS_SAT>::iterator it = itb;
		itb++;
		if (!it->second.get_roverflag())
			BDSobs.erase(it);
		else {
			if (!it->second.get_baseflag()) {
				NoB_R_BDSobs[it->first] = it->second;
				BDSobs.erase(it);
			}
		}
	}
	return;
}
void receiver::displayobs() {
	if (GPSobs.empty()) {
		cout << "no common GPS" << endl;
	}
	else {
		cout << "common GPS measurements:" << endl;
		map<int, OBS_SAT>::iterator it;
		for (it = GPSobs.begin(); it != GPSobs.end(); it++) {
			cout << it->first << "  ";
			it->second.displayOBS_SAT();
		}
	
	}
	if(BDSobs.empty())
		cout << "no common BDS" << endl;
	else {
		cout << "common BDS measurements:" << endl;
		map<int, OBS_SAT>::iterator it;
		for (it = BDSobs.begin(); it != BDSobs.end(); it++) {
			cout << it->first << "  ";
			it->second.displayOBS_SAT();
		}
	
	}
}
void receiver::displaysatpos(int mode) {
	if (mode < 3) {// 打印obs有关
		if (GPSobs.empty()) {
			cout << "no commoon GPS" << endl;
		}
		else {
			cout << "Position of common GPS:" << endl;
			map<int, OBS_SAT>::iterator it;
			for (it = GPSobs.begin(); it != GPSobs.end(); it++) {
				cout << it->first << "  ";
				it->second.displayOBS_SAT_Pos(mode);
			}
		}
		cout << endl;
		if (BDSobs.empty()) {
			cout << "no commoon BDS" << endl;
		}
		else {
			cout << "Position of common BDS:" << endl;
			map<int, OBS_SAT>::iterator it;
			for (it = BDSobs.begin(); it != BDSobs.end(); it++) {
				cout << it->first << "  ";
				it->second.displayOBS_SAT_Pos(mode);
			}
		}
	}
	else if (mode == 3) {//打印位置钟差
		cout << endl << "current time:" << endl;
		cout << setprecision(17) << sol.gettime().time << '.' << sol.gettime().sec;
		cout <<endl<< "current position:" << endl;
		for (int i = 0; i < 3; i++)
			cout <<setprecision(17)<< sol.getpos()[i] << "  ";
		cout << endl << "current velocity:" << endl;
		for (int i = 0; i < 3; i++)
			cout << setprecision(17) << sol.getvel()[i] << "  ";
		/*cout << endl << "当前时刻钟差:" << endl;
		for (int i = 0; i < sol.getdtr().size(); i++)
			cout << setprecision(17) << sol.getdtr()[i] << "  ";
		cout << endl << "当前时刻的协方差:" << endl;
		for (int i = 0; i < sol.getqr().size(); i++)
			cout << setprecision(17) << sol.getqr()[i] << "  ";*/
		//cout << endl << "当前伪距残余:" << endl;
		//cout << v << endl;
	}
	else if (mode == 4) { // print info about DOP and raim
		cout << endl << "TIME:" << endl;
		cout << setprecision(17) << sol.gettime().time << '.' << sol.gettime().sec;
		cout << endl << "DOP:" << endl;
		for (int i = 0; i < DOP.size(); i++)
			cout << setprecision(17) << DOP[i] << "  ";
	}
}
void receiver::displayddres() {
	map<int, double >::iterator itPG =  SatSat_PdresG.begin();//key 值是另外一颗星
	map<int, double >::iterator itLG = SatSat_LdresG.begin();
	map<int, double >::iterator itPB = SatSat_PdresB.begin();//key 值是另外一颗星
	map<int, double >::iterator itLB = SatSat_LdresB.begin();
	cout << "GPS ddres" << endl;
	cout << "ddres of Phase:" << endl;
	for (; itLG != SatSat_LdresG.end(); itLG++) {
		cout << "G" << SATrefGPSID << "-" << "G" << itLG->first << ": " << setprecision(12) << itLG->second << endl;
	}
	
	cout << "ddres of Code:" << endl;
	for (; itPG != SatSat_PdresG.end(); itPG++) {
		cout << "G" << SATrefGPSID << "-" << "G" << itPG->first << ": " << setprecision(12) << itPG->second << endl;
	}
	cout << "BDS ddres" << endl;
	cout << "ddres of Phase:" << endl;
	for (; itLB != SatSat_LdresB.end(); itLB++) {
		cout << "C" << SATrefBDSID << "-" << "C" << itLB->first << ": " << setprecision(12) << itLB->second << endl;
	}
	cout << "ddres of Code:" << endl;
	for (; itPB != SatSat_PdresB.end(); itPB++) {
		cout << "C" << SATrefBDSID << "-" << "C" << itPB->first << ": " << setprecision(12) << itPB->second << endl;
	}

}

void receiver::soltofile(int LS_stat,ofstream & out) {
	gtime_t timestamp = sol.gettime();
	double soltime = timestamp.time + timestamp.sec;
	string stime;
	if (!sol.get_state())
		return;
	if (TC_conf.timef == 1) {
		out << setiosflags(ios::fixed) << setprecision(3) << soltime << "   ";
	}
	else if (TC_conf.timef == 2) { // Week+Seconds of week
		int week;
		double sow;
		//gtime_t timestamp(soltime);
		sow = time2gpst(timestamp, &week);
		out << week << " " << setiosflags(ios::fixed) << setprecision(3) << sow << "   ";
	}
	else { // YMDHMS
		//gtime_t timestamp(soltime);
		string stime;
		time2str(timestamp, stime, 3);
		out << stime << "  ";
	}
	time2str(timestamp,stime,3);
	/*if (!LS_stat)//是不是为了应对之前的什么bug？
		return;*/
	
	if (output_pos_format == "xyz") {
		//out << stime << "  ";
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sol.getpos()[i] << "   ";
		}
		out << sol.get_state() << "  "; 
		for (int i = 0; i <3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqrt(sol.getqr()[i]) << "   ";
		}
		for (int i = 3; i < 6; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(sol.getqr()[i]) << "   ";
		}
	}
	else if (output_pos_format == "llh") {
		//out << stime << "  ";
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sol.getposllh()[i] << "   ";
		}
		out << sol.get_state() << "  ";
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqrt(sol.getqrenu()[i]) << "   ";
		}
		for (int i = 3; i < 6; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(sol.getqrenu()[i]) << "   ";
		}
	}
	else if (output_pos_format == "enu") {
		//double UTCtt = timestamp.time + timestamp.sec-18;
		//out << setiosflags(ios::fixed) << setprecision(4)<< UTCtt << "  ";
		vector<double > posbasellhxyz(3, 0);
		vector<double > posxyz(3,0);
		ecef2pos(basepos, posbasellhxyz);// base 的位置转llh
		vector<double > e_rb(3,0);
		posxyz = sol.getpos();
		for (int i = 0; i < 3; i++)
			e_rb[i] = posxyz[i] - basepos[i];
		vector<double > enu(3, 0);
		ecef2enu(posbasellhxyz, e_rb, enu);
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << enu[i] << "   ";
		}
		out << sol.get_state() << "  ";
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqrt(sol.getqrenu()[i]) << "   ";
		}
		for (int i = 3; i < 6; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(sol.getqrenu()[i]) << "   ";
		}
	
	}
	out << endl;
	/*if (plot_mode) {
		plot_ENU_traj();
	}*/
}
void receiver::updateSatinfo() {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	vector<double > pos_velb(6, 0);
	vector<double > bias_driftb(2,0);
	vector<double > pos_velr(6, 0);
	vector<double > bias_driftr(2, 0);
	double varb,varr;
	varb = 0;
	varr = 0;
	int svh;
	double TGD;
	//double bias_tt;
	/*double Cspeed = 299792458.0;*/
	Station_OBS base, rover;
	int bflag;
	int rflag;
	for (; itg != GPSobs.end(); itg++) {
		pos_velb.assign(6, 0);
		bias_driftb.assign(2, 0);
		pos_velr.assign(6,0);
		bias_driftr.assign(2, 0);
		bflag = itg->second.get_baseflag();
		rflag = itg->second.get_roverflag();
		int tb,tr;
		double sb,sr;
		gtime_t timeb,timer;
		double Pb, Pr;
		int BR_unusedflag = 0;
		if (bflag) {
			base = itg->second.get_base();
			base.gettime(&tb, &sb);
			timeb = gtime_t(tb, sb);
			Pb = itg->second.get_baseP(0);
			timeb = timeadd(timeb, -Pb / Cspeed);
			if (!allsateph.Cal_OneSatinfo(itg->first, 'G', timeb, pos_velb, bias_driftb, varb, svh, TGD))
				BR_unusedflag = 1;
		}
		if (rflag) {
			rover = itg->second.get_rover();
			rover.gettime(&tr, &sr);
			timer = gtime_t(tr, sr);
			Pr = itg->second.get_roverP(0);
			timer = timeadd(timer, -Pr / Cspeed);
			if (!allsateph.Cal_OneSatinfo(itg->first, 'G', timer, pos_velr, bias_driftr, varr, svh, TGD)) {
				BR_unusedflag = 1;
				itg->second.setroverunusedflag(1);
			}
		}
		itg->second.set_BR_unusedflag(BR_unusedflag);
		itg->second.setpvbdst(pos_velb, bias_driftb, varb, pos_velr, bias_driftr, varr, svh, TGD, bflag, rflag);
	}
	for (; itb != BDSobs.end(); itb++) {
		pos_velb.assign(6, 0);
		bias_driftb.assign(2, 0);
		pos_velr.assign(6, 0);
		bias_driftr.assign(2, 0);
		bflag = itb->second.get_baseflag();
		rflag = itb->second.get_roverflag();
		int tb, tr;
		double sb, sr;
		gtime_t timeb, timer;
		double Pb, Pr;
		int BR_unusedflag = 0;
		if (numiterator >= debug_stop && itb->first == 5) {
			int aaaa = 1;
		}
		if (bflag) {
			base = itb->second.get_base();
			base.gettime(&tb, &sb);
			timeb = gtime_t(tb, sb);
			Pb = itb->second.get_baseP(0);
			timeb = timeadd(timeb, -Pb / Cspeed);
			if (!allsateph.Cal_OneSatinfo(itb->first, 'C', timeb, pos_velb, bias_driftb, varb, svh, TGD))
				BR_unusedflag = 1;
		}
		if (rflag) {
			rover = itb->second.get_rover();
			rover.gettime(&tr, &sr);
			timer = gtime_t(tr, sr);
			Pr = itb->second.get_roverP(0);
			timer = timeadd(timer, -Pr / Cspeed);
			/*if ((numiterator >= debug_stop) && (itb->first == 5)) {
				int  aaaa = 1;
			}*/
			if (!allsateph.Cal_OneSatinfo(itb->first, 'C', timer, pos_velr, bias_driftr, varr, svh, TGD)) {
				itb->second.setroverunusedflag(1);
				BR_unusedflag = 1;
			}
		}
		itb->second.set_BR_unusedflag(BR_unusedflag);
		itb->second.setpvbdst(pos_velb, bias_driftb, varb, pos_velr, bias_driftr, varr, svh, TGD, bflag, rflag);
	}
	return;
}
int receiver::DOPcheck() {
	vector<double >azel,tmp;
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	for (; itg != GPSobs.end(); itg++) {
		if (!itg->second.get_rover().get_unusedflag()) {
			tmp = itg->second.get_rover().get_razellos();
			azel.push_back(tmp[1]);
			azel.push_back(tmp[2]);
		}
	}
	for (; itb != BDSobs.end(); itb++) {
		if (!itb->second.get_rover().get_unusedflag()) {
			tmp = itb->second.get_rover().get_razellos();
			azel.push_back(tmp[1]);
			azel.push_back(tmp[2]);
		}
	}
	if (azel.size() < 8) {
		cout << "Satellite is too less， cannot calculate DOP" << endl;
		return 0;
	}
	calculate_DOP(azel);
	if (DOP[0] > TC_conf.maxgdop) {
		cout << "Current GDOP is too big!" << endl;
		return 0;
	}
	return 1;
		
}
void receiver::calculate_DOP(const vector<double > &azel) {
	int num = azel.size() / 2;
	Eigen::MatrixXd H(num, 4);
	Eigen::MatrixXd Q(4, 4);
	double cosel, sinel;
	for (int i = 0; i < num; i++) { // follow the algorithm of RTKLIB but why disregard the clock bias of other constellations
		cosel = cos(azel[1 + i * 2]);
		sinel = sin(azel[1 + i * 2]);
		H.row(i) << cosel * sin(azel[i * 2]), cosel*cos(azel[i * 2]), sinel, 1.0;
	}
	Q = H.transpose()*H;
	Eigen::MatrixXd Q_inv = Q.inverse();
	DOP[0] = sqrt(Q_inv(0, 0) + Q_inv(1, 1) + Q_inv(2, 2) + Q_inv(3, 3));
	DOP[1] = sqrt(Q_inv(0, 0) + Q_inv(1, 1) + Q_inv(2, 2));
	DOP[2] = sqrt(Q_inv(0, 0) + Q_inv(1, 1));
	DOP[3] = sqrt(Q_inv(2, 2));
	return;
}


int receiver::LS_pos(char mode){
	int MAXITR = 10;
	Ls_x.resize(5);
	Ls_x << sol.getpos()[0], sol.getpos()[1], sol.getpos()[2], sol.getdtr()[0], sol.getdtr()[1];
	//cout << Ls_x(0);
	//测试
	/*Ls_x[0] = 0;
	Ls_x[1] = 0;
	Ls_x[2] = 0;*/
	map<int, OBS_SAT> ::iterator itg = GPSobs.begin();
	map<int, OBS_SAT> ::iterator itb = BDSobs.begin();
	int Grover_valid_num = 0;
	int Brover_valid_num = 0;
	for (; itg != GPSobs.end(); itg++) {
		if (itg->second.get_roverflag())
			Grover_valid_num++;
	}
	for (; itb != BDSobs.end(); itb++) {
		if (itb->second.get_roverflag())
			Brover_valid_num++;
	}
	if (Grover_valid_num+Brover_valid_num<=0)//是不是可以改成3
		return 0;
	gtime_t Maybesoltime; 
	if (Grover_valid_num) {
		map<int, OBS_SAT>::iterator itg = GPSobs.begin();
		int t = 0;
		double s = 0;
		itg->second.getBRtime(&t, &s);
		Maybesoltime.time = t;
		Maybesoltime.sec = s;
	}
	if (Brover_valid_num) {
		map<int, OBS_SAT>::iterator itb = BDSobs.begin();
		int t = 0;
		double s = 0;
		itb->second.getBRtime(&t, &s);
		Maybesoltime.time = t;
		Maybesoltime.sec = s;
	}
	sol.settime(Maybesoltime);
	int i = 0;
	int countG = 0;
	int countB = 0;
	for (i = 0; i < MAXITR; i++) {
		static int last_dim;
		if (numiterator >= debug_stop && i ==2) {
			int  aaaa = 1;
		}
		int dim = rescode(allsateph.getGPSIONCORR(), allsateph.getBDSIONCORR(),countG,countB);
		if (mode == 'F' && dim < 6) {
			cout << endl << "less than 6, cannot FDE" << endl;
			return 0;
		}
		if (dim < 5) { //如果少于5 就也可以定
			cout <<endl<< "less than 5，cannot single point" << endl;
			return 0;
		}	
		if (i == 0) {
			//int rows = dim + (!Grover_valid_num || !Brover_valid_num ? 1 : 0);// avoid rank-deficient
			int rows = dim;
			H = Eigen::MatrixXd::Zero(rows,5);//只考虑两星座 TOT
			v = Eigen::VectorXd::Zero(rows, 1);
		}
		else if (dim != last_dim) {
			//int rows = dim + (!Grover_valid_num || !Brover_valid_num ? 1 : 0);// avoid rank-deficient
			int rows = dim;
			H.resize(rows, 5);
			v.resize(rows, 1);
		}
		last_dim = dim;
		
		if (LS_pos_body(mode, countG, countB)) {
			dim = rescode(allsateph.getGPSIONCORR(), allsateph.getBDSIONCORR(),countG,countB);
			break;
		}
			
	}
	if (i >= MAXITR) {
		cout << endl << "最小二乘没收敛" << endl;
		
		return -1;
	}
	
	return 1;
}
int receiver::rescode(const vector<double > &GPS_ION_CORR, const vector<double > &BDS_ION_CORR,int &countG,int &countB) {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	vector<double > Pos_ecef(3, 0);
	vector<double > Pos_llh(3,0);
	vector<double > sysclk_bias(2, 0);

	for (int i = 0; i < 3; i++)
		Pos_ecef[i] = Ls_x(i);
	for(int i = 0;i<2;i++)
		sysclk_bias[i] = Ls_x(i+3);
	
	ecef2pos(Pos_ecef,Pos_llh);
	
	countG = 0;
	countB = 0;
	int count;
	for (; itg != GPSobs.end(); itg++) {
	/*	if (numiterator >= debug_stop && itg->first == 21) {
			int aaa = 1;
		}*/
		//if (!itg->second.get_roverflag()||itg->second.get_rover().get_unusedflag()||itg->second.get_BR_fdemask())// 如果这里unusedflag 拒绝就无法改变LS中反复横跳的星了
		if (!itg->second.get_roverflag() || itg->second.get_BR_fdemask())
			continue;
		itg->second.Ls_update_razellos(sysclk_bias,Pos_ecef, Pos_llh,countG,GPS_ION_CORR);//这里
	
	}
	for (; itb != BDSobs.end(); itb++) {
		if (numiterator >= debug_stop && itb->first == 5) {
			int aaa = 1;
		}
		//if (!itb->second.get_roverflag() || itb->second.get_rover().get_unusedflag()|| itb->second.get_BR_fdemask())
		if (!itb->second.get_roverflag() || itb->second.get_BR_fdemask())
			continue;
		itb->second.Ls_update_razellos(sysclk_bias,Pos_ecef, Pos_llh,countB, GPS_ION_CORR);//用的是GPS的ION参数，然后再用频率比去转
	}
	count = countG + countB;
	if ((!countG && countB) || (!countB && countG))// avoid rank-deficient
		count++;
	return count;
}
int  receiver::LS_pos_body(char mode, int countG, int countB) {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	vector<double > razellos;
	int i = 0;
	double res, res_var;
	int t = 0; double s = 0.0;
	if (countG) { //  have valid GPS measurements
		for (; itg != GPSobs.end(); itg++) {
			/*if (numiterator >= debug_stop) {
				int aaaa = 1;
			}*/
			if (!itg->second.get_roverflag())
				continue;
			Station_OBS rover = itg->second.get_rover();
			if (!rover.get_unusedflag() && !rover.get_fdemask()) {
				if (t == 0)
					rover.gettime(&t, &s);
				razellos = rover.get_razellos();
				res = rover.get_P_zres();
				Ls_res += res * res;
				res_var = sqrt(rover.get_P_zvar());
				v(i) = res / res_var;
				H.row(i) << -razellos[3] / res_var, -razellos[4] / res_var, -razellos[5] / res_var, 1 / res_var, 0;
				//cout << H.row(i)<<endl;
				i++;
			}
		}
	}
	else {
		v(i) = 0;
		H.row(i) << 0, 0, 0, 1, 0;
		//cout << H.row(i)<<endl;
		i++;
	}
	if (countB) {//  have valid BDS measurements
		for (; itb != BDSobs.end(); itb++) {
			/*if (numiterator >= debug_stop && itb->first == 5) {
				int aaaa = 1;
			}*/
			if (!itb->second.get_roverflag())
				continue;
			Station_OBS rover = itb->second.get_rover();
			if (!rover.get_unusedflag() && !rover.get_fdemask()) {
				if (t == 0)
					rover.gettime(&t, &s);
				razellos = rover.get_razellos();
				res = rover.get_P_zres();
				Ls_res += res * res;
				res_var = sqrt(rover.get_P_zvar());
				v(i) = res / res_var;
				
				H.row(i) << -razellos[3] / res_var, -razellos[4] / res_var, -razellos[5] / res_var, 1 / res_var, 1 / res_var;
				i++;
			}
		}
	}
	else {
		v(i) = 0;
		H.row(i) << 0, 0, 0, 0, 1;
		//cout << H.row(i)<<endl;
		i++;
	}
	Ls_res = sqrt(Ls_res / i);
	/*if (numiterator >= debug_stop) {
		matrixtofile(H, "H", 'T');
		matrixtofile(v, "v", 'N');
		int aaaa = 1;
	}*/
	//
	/*cout << endl << "Hv:" << endl;
	cout <<setprecision(12)<< H << endl;
	cout << setprecision(12) << v << endl;*/
	Eigen::MatrixXd HT = H.transpose();
	Eigen::MatrixXd HV = HT * v;
	Eigen::MatrixXd HTH = HT*H;
	//cout << setprecision(12) << HTH << endl;
	/*cout << setprecision(12) << H << endl;
	cout << setprecision(12) << HTH << endl;*/
	/*cout << setprecision(12) << HV << endl;*/
	Eigen::MatrixXd HTH_inv = HTH.inverse();
	//cout << setprecision(12) << HTH_inv << endl;
	Eigen::VectorXd delta_x = HTH_inv *HV;
	gtime_t time(t, s);
	for (int j = 0; j < Ls_x.size(); j++) 
		Ls_x[j] += delta_x[j];
	if (numiterator >= debug_stop) {
		
		int aaaa = 1;
	}
	
	if (delta_x.norm() < 1E-4 && (mode == 'N'||(mode =='F'&&Ls_res< TC_conf.Ls_res_threshold))) {
		/*gtime_t lastsoltime = sol.gettime();*/
		
		sol.settime(timeadd(time, -Ls_x[3] / Cspeed));
		/*if (lastsoltime.time == 0)
			tt = 0.0;
		else
			tt = timediff(sol.gettime(), lastsoltime);*/
		vector<double > dtr(2, 0);
		dtr[0] = Ls_x[3] / Cspeed;
		dtr[1] = Ls_x[4] / Cspeed;
		sol.setdtr(dtr);
		sol.setpos(Ls_x[0], Ls_x[1], Ls_x[2]);
		sol.setvel(0, 0, 0);
		vector<double > qr(6, 0);
		qr[0] = HTH_inv(0, 0);
		qr[1] = HTH_inv(1, 1);
		qr[2] = HTH_inv(2, 2);
		qr[3] = HTH_inv(1, 0);
		qr[4] = HTH_inv(2, 1);
		qr[5] = HTH_inv(2, 0);
		sol.setqr(qr);
		sol.set_state(4);
		sol.PQxyz2PQllh();
		return 1;
	} //没有写validobs函数
	return 0;
}
void receiver::settt(gtime_t lastsoltime) {
	if (lastsoltime.time == 0)
		tt = 0.0;
	else
		tt = timediff(sol.gettime(), lastsoltime);
	return;
}
int receiver::LS_raim() {
	double normv = v.norm();
	double normv2 = normv * normv;
	if (normv2 > chisqr[v.size() - 5 - 1]) { // we mandatorily set the dimension.
		cout << "residual is too big" << endl;
		//return 0;
		//RTKLIB:
		/* threshold too strict for all use cases, report error but continue on */
	}
	return 1;
}
void receiver::LS_fde() {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	for (; itg != GPSobs.end(); itg++) {
		Station_OBS rover = itg->second.get_rover();
		rover.set_fdemask(1);
		Ls_res = 0.0;
		LS_pos('F');
		rover.set_fdemask(0);
	}
	for (; itb != BDSobs.end(); itb++) {
		Station_OBS rover = itb->second.get_rover();
		rover.set_fdemask(1);
		Ls_res = 0.0;
		LS_pos('F');
		rover.set_fdemask(0);
	}
	return;
}
/*与最小二乘多普勒定速有关的函数*/
int receiver::LS_vel() {
	Ls_x.resize(4);
	Ls_x.setZero();
	int MAXITR = 10;
	int i = 0;
	int state;
	for (i = 0; i < MAXITR; i++) {
		int dim = resdop();
		if (dim < 4) {
			cout << "Satellites are too less to estimate velocity " << endl;
			return 0;
		}
		if (i == 0) {
			H = Eigen::MatrixXd::Zero(dim, 4);//只考虑两星座 TOT
			v = Eigen::VectorXd::Zero(dim, 1);
		}
		if (state = LS_vel_body())
			break;
	}
	if (i >= MAXITR) {
		cout << "Ls_vel diverges "<<endl;
		return 0;
	}
	return state;
}
int receiver::resdop() { 
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	Eigen::Matrix3d E;
	vector<double > Pos_ecef(3, 0);
	vector<double > Pos_llh(3, 0);
	//vector<double > sysclk_bias(2, 0);
	for (int i = 0; i < 3; i++)
		Pos_ecef[i] = sol.getpos()[i];
	ecef2pos(Pos_ecef, Pos_llh);
	xyz2enu(Pos_llh,E);
	int count = 0;

	for (; itg != GPSobs.end(); itg++) {
		if (!itg->second.get_roverflag())
			continue;
		Station_OBS rover = itg->second.get_rover();
		itg->second.Ls_updateDres(Pos_ecef,E,Ls_x,count);
		
	}
	for (; itb != BDSobs.end(); itb++) {
		if (!itb->second.get_roverflag())
			continue;
		Station_OBS rover = itb->second.get_rover();
		/*if (itb->first == 4) {
			int aaa = 1;
		}*/
		itb->second.Ls_updateDres(Pos_ecef, E, Ls_x,count);
	}
	return count;
}

int receiver::LS_vel_body() {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	vector<double > razellos;
	int i = 0;
	double res;
	for (; itg != GPSobs.end(); itg++) {
		if (!itg->second.get_roverflag())
			continue;
		Station_OBS rover = itg->second.get_rover();
		if (!rover.get_unusedflag() && !rover.get_fdemask()) {
			razellos = rover.get_razellos();
			res = rover.get_D_zres();
			v(i) = res;
			H.row(i) << -razellos[3], -razellos[4], -razellos[5], 1;
			i++;
		}
		
	}
	for (; itb != BDSobs.end(); itb++) {
		if (!itb->second.get_roverflag())
			continue;
		Station_OBS rover = itb->second.get_rover();
		if (!rover.get_unusedflag() && !rover.get_fdemask()) {
			razellos = rover.get_razellos();
			res = rover.get_D_zres();
			v(i) = res;
			H.row(i) << -razellos[3], -razellos[4], -razellos[5], 1;
			i++;
		}
		
	}
	
	//matrixtofile(H, "H", 'T');
	Eigen::MatrixXd HT = H.transpose();
	Eigen::MatrixXd HV = HT * v;
	Eigen::MatrixXd HTH = HT * H;
	/*cout << "当前历元：" << endl;
	cout << setprecision(12) << H << endl;*/
	/*cout << setprecision(12) << HTH << endl;*/
	/*cout << setprecision(12) << HV << endl;*/
	Eigen::MatrixXd HTH_inv = HTH.inverse();
	//cout << setprecision(12) << HTH_inv << endl;
	Eigen::VectorXd delta_x = HTH_inv * HV;
	for (int j = 0; j < Ls_x.size(); j++)
		Ls_x[j] += delta_x[j];
	if (delta_x.norm() < 1E-6) {
		sol.setvel(Ls_x[0], Ls_x[1], Ls_x[2]);
		sol.setdrift(Ls_x[3]);
		return 1;
	}
	return 0;
}
int receiver::tdcp2vel(Eigen::Vector3d &vn) {
	if (!sol.getpos()[0])
		return 0;
	vector<double> last_sol_pos = sol.getpos();
	map<int, double> Sat_TDCPresG;
	map<int, double> Sat_TDCPresB;
	double dt = 0;
	if (!res_tdcp(last_sol_pos, Sat_TDCPresG, Sat_TDCPresB))
		return 0; 
	Eigen::Vector3d ve;
	vector<double> xllh(3, 0);
	Eigen::Matrix3d E;
	ecef2pos(last_sol_pos, xllh);
	xyz2enu(xllh, E);
	int over_flag = Form_Solve_TDcp(Sat_TDCPresG, Sat_TDCPresB, ve);
	vn = E * ve;
	return over_flag;
}

int receiver::res_tdcp(const vector<double > & last_sol_pos, map<int, double>&Sat_TDCPresG, map<int, double>&Sat_TDCPresB) {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	
	double sum = 0;
	for (; itg != GPSobs.end(); itg++) {
		if (!itg->second.get_roverflag() || itg->second.get_BR_fdemask())
			continue;
		double TDCPresG = itg->second.calculate_TDCPres(last_sol_pos);
		if (!TDCPresG)
			continue;
		Sat_TDCPresG[itg->first] = TDCPresG;
	}
	for (; itb != BDSobs.end(); itb++) {
		if (!itb->second.get_roverflag() || itb->second.get_BR_fdemask())
			continue;
		double TDCPresB = itb->second.calculate_TDCPres(last_sol_pos);//这里
		if (!TDCPresB)
			continue;
		Sat_TDCPresB[itb->first] = TDCPresB;
	}
	int nv = Sat_TDCPresG.size() + Sat_TDCPresB.size();
	if (nv < 4)
		return 0;
	if (!TDcp_robust_exclude(Sat_TDCPresG, Sat_TDCPresB)) 
		return 0;
	return 1;
}
int receiver::TDcp_robust_exclude(map<int, double> &Sat_TDCPresG, map<int, double> &Sat_TDCPresB) {
	int nvG = Sat_TDCPresG.size(),nvB = Sat_TDCPresB.size();
	int nv = nvG + nvB;
	double sum = 0;
	Eigen::VectorXd v(nv);
	Eigen::VectorXd v_ID(nv);
	vector<int> exec;
	double ave_v;
	double ave_distance;
	int index = 0;
	map<int, double> ::iterator it = Sat_TDCPresG.begin();
	for (; it != Sat_TDCPresG.end(); it++) {
		v[index] = it->second;
		v_ID[index++] = it->first;
		sum += fabs(it->second);
	}
	for (it = Sat_TDCPresB.begin(); it != Sat_TDCPresB.end(); it++) {
		v[index] = it->second;
		v_ID[index++] = it->first+NSATGPS;
		sum += fabs(it->second);
	}
	for (int i = 0; i < nv; i++) {
		double vi = fabs(v[i]);
		double ave_v = (sum - vi) / (nv - 1);
		double sum_distance = 0;
		for (int j = 0; j < nv; j++) {
			if (j == i) continue;
			sum_distance += fabs(vi - fabs(v[j]));
		}
		double ave_distance = sum_distance / (nv - 1);
		if (ave_distance > 3 * ave_v || (ave_distance > 100 && ave_distance > 1.2*ave_v)) {
			exec.push_back(v_ID[i]);
		}
	}
	for (int i = 0; i < exec.size(); i++) {
		if (exec[i] < NSATGPS) {
			Sat_TDCPresG.erase(exec[i]);
		}
		else {
			Sat_TDCPresB.erase(exec[i]-NSATGPS);
		}
	}
	if(Sat_TDCPresG.size()+ Sat_TDCPresB.size()<4)
		return 0;
	return 1;
}
int receiver::Form_Solve_TDcp(map<int, double> &Sat_TDCPresG, map<int, double> &Sat_TDCPresB, Eigen::Vector3d & ve) {
	map<int, double>::iterator itSTG = Sat_TDCPresG.begin();
	map<int, double>::iterator itSTB = Sat_TDCPresB.begin();
	int nvG = Sat_TDCPresG.size();
	int nvB = Sat_TDCPresB.size();
	int nv = nvG + nvB;
	Eigen::VectorXd v0(nv);
	Eigen::MatrixXd H0(nv,4);
	Eigen::MatrixXd P0 = Eigen::MatrixXd::Zero(nv,nv);
	vector<double> ral;
	double dt = 0.0;
	int index = 0;
	for (; itSTG != Sat_TDCPresG.end(); itSTG++) {
		v0[index] = itSTG->second;
		ral = GPSobs[itSTG->first].get_rover().get_razellos();
		if (!dt) {
			dt = timediff(GPSobs[itSTG->first].get_rover().gettime(), GPSobs[itSTG->first].get_rover().getlasttime());
		}
		H0.row(index) << -ral[3], -ral[4], -ral[5], 1;
		P0(index, index) = (sin(ral[2])*sin(ral[2])) / 2;
		index++;
	}
	for (; itSTB != Sat_TDCPresB.end(); itSTB++) {
		v0[index] = itSTB->second;
		ral = BDSobs[itSTB->first].get_rover().get_razellos();
		if (!dt) {
			dt = timediff(BDSobs[itSTB->first].get_rover().gettime(),BDSobs[itSTB->first].get_rover().getlasttime());
		}
		H0.row(index) << -ral[3], -ral[4], -ral[5], 1;
		P0(index, index) = (sin(ral[2])*sin(ral[2])) / 2;
		index++;
	}
	//matrixtofile(P0, "P0", 'N');
	Eigen::MatrixXd N = H0.transpose()*P0*H0;
	Eigen::VectorXd dx = N.inverse()*H0.transpose()*P0*v0;
	ve << dx(0) / dt,dx(1)/dt,dx(2)/dt;
	if (ve.dot(ve) > 3)
		return 1;
	return 0;
}

/*与RTK定位有关的函数*/
//RTKLIB 中 每历元RTK之前都会进行一次LS_Pos 
//我们搭配两者使用，如果前面有个LS_Pos,那么读取sol就是Ls之后的
//否则就是上次RTK的。
int receiver::RTK_pos(string RTKmode)
{
	choosecommonsat();
	
	int niter = 1;
	int nv_dim = 0;//观测量双差对数
	int nb = 0;// 参与模糊度解算的双差对数
	double *RTK_xfix = new double[nx] {0};
	vector<int> ix;
	if (RTKmode == "PL") {
		/*if(sol.gettime().time == 1614221569)
			tt = 1.1520002985615525;*/
		/*if (numiterator >= debug_stop) {
			int aaaa = 1;
		}*/
		RCV_zdres('B',NULL);
		if (validsatnum() <= 0) {
			cout << "no common satellites,do not run RTK" << endl;
			return 0;
		}
		
		if (numiterator >= debug_stop) {
			int aaaa = 1;
			/*sol.setpos(-2168786.0354538686, 4386832.4625583608, 4076910.9500581124);
			sol.setvel(-0.023762395795798166, -0.037954004905715009, 0.012549181822045194);*/
		}
		udstate();
		
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		for (int ii = 0; ii < niter; ii++) {
			RCV_zdres('R',NULL);
			nv_dim = RCV_ddres(RTKmode,NULL);//注意是算上伪距 相位 两种观测量的双差对数，
			
			//displayddres();
			if(nv_dim>=1)
				MeasureUpdate(RTKmode);
			
		}
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		if (!solve_validpos(RTKmode,NULL,1)) {
			/*if (numiterator >= 147)
				return 0;*/
			nfix = 0;
			cout << "有效卫星数太少" << endl;
			//sol.set_state(0);//可以相信别的方法，这里先设置成0
			return 0;
		}
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		
		//如果能到这里，至少有浮点解
		sol.set_state(2);
		if (numiterator >= debug_stop) {
			
			//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
			int aaaa = 1;
		}
		
		if ((nb = manage_amb_LAMBDA(RTK_xfix,ix)) > 1) {
			int firstvalidGPSID = SATrefGPSID;
			int  firstvalidBDSID = SATrefBDSID;//注意，这里不能轻易改，因为solve_validpos会改变SATref的值
			if (solve_validpos(RTKmode, RTK_xfix,0)) {
				if (++nfix >= TC_conf.minfix && TC_conf.armode == 3) {
				//if (++nfix >= 1) {
					if (numiterator >= debug_stop) {
						int aaaa = 1;
					}
					holdamb(RTK_xfix,ix,nb, firstvalidGPSID,firstvalidBDSID);
					if (numiterator >= debug_stop) {
						int aaaa = 1;
					}
				}
				sol.set_state(1);
			}
		}
		RTK_tosol();
		sol.PQxyz2PQllh();
		/*if (numiterator == debug_stop) {
			int aaaa = 1;
		}*/
		//lock_solve();
		
	}
	delete RTK_xfix;
	return 0;
}
void receiver::RTK_tosol() {
	if (sol.get_state() == 1) {
		sol.setpos(RTK_fixXpva[0], RTK_fixXpva[1], RTK_fixXpva[2]);
		vector<double > qr(6, 0);
		vector<double > qv(6, 0);
		qr[0] = RTK_fixPpva(0,0);
		qr[1] = RTK_fixPpva(1, 1);
		qr[2] = RTK_fixPpva(2, 2);
		qr[3] = RTK_fixPpva(1, 0);
		qr[4] = RTK_fixPpva(2, 1);
		qr[5] = RTK_fixPpva(2, 0);
		sol.setqr(qr);
		if (dynamics_on) {
			sol.setvel(RTK_fixXpva[3], RTK_fixXpva[4], RTK_fixXpva[5]);
			qv[0] = RTK_fixPpva(3, 3);
			qv[1] = RTK_fixPpva(4, 4);
			qv[2] = RTK_fixPpva(5, 5);
			qv[3] = RTK_fixPpva(4, 3);
			qv[4] = RTK_fixPpva(5, 4);
			qv[5] = RTK_fixPpva(5, 3);
			sol.setqv(qv);
		}
	}
	else {
		sol.setpos(RTK_x[0], RTK_x[1], RTK_x[2]);
		vector<double > qr(6, 0);
		vector<double > qv(6, 0);
		qr[0] = RTK_P[0+ 0*nx];
		qr[1] = RTK_P[1 + 1 * nx];
		qr[2] = RTK_P[2 + 2 * nx];
		qr[3] = RTK_P[1 + 0 * nx];
		qr[4] = RTK_P[1 + 2 * nx];
		qr[5] = RTK_P[2 + 0 * nx];
		sol.setqr(qr);
		if (dynamics_on) {
			sol.setvel(RTK_x[3], RTK_x[4], RTK_x[5]);
			qv[0] = RTK_P[3 + 3 * nx];
			qv[1] = RTK_P[4 + 4 * nx];
			qv[2] = RTK_P[5 + 5 * nx];
			qv[3] = RTK_P[4 + 3 * nx];
			qv[4] = RTK_P[5 + 4 * nx];
			qv[5] = RTK_P[5 + 3 * nx];
			sol.setqv(qv);
		}
		nfix = 0;
	}
	
}
void receiver::lock_solve() {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	if (SATrefGPSID) {
		for (; itg != GPSobs.end(); itg++) {
			if (itg->second.get_BR_fdemask() || itg->second.get_BR_unusedflag())
				continue;
			if (itg->first == SATrefGPSID|| SatSat_PdresG.count(itg->first)|| SatSat_LdresG.count(itg->first))
				lock_G[itg->first - 1]++;
		}
	}
	if (SATrefBDSID) {
		for (; itb != BDSobs.end(); itb++) {
			if (itb->second.get_BR_fdemask() || itb->second.get_BR_unusedflag())
				continue;
			if (itb->first == SATrefBDSID || SatSat_PdresB.count(itb->first) || SatSat_LdresB.count(itb->first))
				lock_B[itb->first - 1]++;
		}
	}
	
	
}
void receiver::out_solve() {
	map<int, double>::iterator itg = SatSat_PdresG.begin();
	map<int, double>::iterator itgL = SatSat_LdresG.begin();
	map<int, double>::iterator itb = SatSat_PdresB.begin();
	map<int, double>::iterator itbL = SatSat_LdresB.begin();
	//int count = 0;
	if (itg != SatSat_PdresG.end() || itgL != SatSat_LdresG.end()) {
		outc_G[SATrefGPSID - 1] = 0;
	}
	for (; itg != SatSat_PdresG.end(); itg++) {
		//count++;
		outc_G[itg->first - 1] = 0;
	}
	for (; itgL != SatSat_LdresG.end(); itgL++) {
		//count++;
		outc_G[itgL->first - 1] = 0;
	}
	if (itb != SatSat_PdresB.end() || itbL != SatSat_LdresB.end()) {
		outc_B[SATrefBDSID - 1] = 0;
	}
	for (; itb != SatSat_PdresB.end(); itb++) {
		//count++;
		outc_B[itb->first - 1] = 0;
	}
	for (; itbL != SatSat_LdresB.end(); itbL++) {
		//count++;
		outc_B[itbL->first - 1] = 0;
	}
}
int receiver::RCV_zdres(char mode, double * state) {
	map<int, OBS_SAT> ::iterator itg = GPSobs.begin();
	map<int, OBS_SAT> ::iterator itb = BDSobs.begin();
	vector<double > posecef;
	double *state_x = (state == NULL ? RTK_x:state);
	if (mode == 'B') {
		 posecef =basepos;
	
	}
	else if (mode == 'R') {
		for(int i = 0;i<3;i++)
			posecef.push_back(state_x[i]);//是状态更新后的当rcv位置，而不是sol
	}
	else {
		cout << "no such zdres mode" << endl;
		return 0;
	}
	if (numiterator >= debug_stop)
		int aaaa = 1;
	vector<double > posllh(3,0);
	ecef2pos(posecef, posllh);//此处没有做Tidecorr
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	for (; itg != GPSobs.end(); itg++) {
		itg->second.BRzdres(mode, posecef, posllh);
	}
	for (; itb != BDSobs.end(); itb++) {
		itb->second.BRzdres(mode, posecef, posllh);
	}
	return 1;
}
int receiver::RCV_ddres(string RTKmode,double * state){
	/*int GPSrefsatID = 0;
	int BDSrefsatID = 0;*/
	double *state_x = (state == NULL ? RTK_x : state);
	/*if (numiterator >= 147) {
		int aaaa = 1;
	}*/
	if (RTKmode == "PL") {// PL 的时候选择伪距相位都有的星，并且俯仰角高于限制，应该增加flag来选择，此处先暂定使用unused_flag来组成双差
		/*if (numiterator == debug_stop) {
			int aaaa = 1;
		}*/
		findrefsat();
		if (SATrefGPSID == 0 && SATrefBDSID == 0) {
			cout << "no valid sat for PL mode" << endl;
			return 0;
		}

		//cout << "G" << GPSrefsatID << "C" << BDSrefsatID << endl;
		calddres(RTKmode,state_x);
		return SatSat_PdresG.size()+ SatSat_LdresG.size()+ SatSat_PdresB.size()+SatSat_PdresB.size();
	}
	return 0;
	
}
void receiver::findrefsat() {
	double min_el = 0;
	double curel = 0;
	map<int, OBS_SAT> ::iterator itg = GPSobs.begin();
	map<int, OBS_SAT> ::iterator itb = BDSobs.begin();
	SATrefGPSID = 0;
	SATrefBDSID = 0;
	for (; itg != GPSobs.end(); itg++) { //选择出来的参考卫星，似乎不会受是否同时提供伪距双差和相位双差的影响。
		if (itg->second.get_BR_unusedflag())
			continue;
		if ( (curel =itg->second.get_rover().get_razellos()[2]) > min_el) {
			min_el = curel;
			SATrefGPSID = itg->first;
		}
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	min_el = 0;
	for (; itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_unusedflag())
			continue;
		if ((curel = itb->second.get_rover().get_razellos()[2]) > min_el) {
			min_el = curel;
			SATrefBDSID = itb->first;
		}
	}
	SatSat_PdresG.clear();
	SatSat_LdresG.clear();
	SatSat_PdresB.clear();
	SatSat_LdresB.clear();
	/*if(*GPSsatID == 0 )
		cout<<"No "*/
	return;
}
void receiver::calddres(string RTKmode,double * state_x) {
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	double refroverP = 0;
	double refroverL = 0;
	double refroverD = 0;
	double refbaseP = 0;
	double refbaseL = 0;
	double refbaseD = 0;
	 double lambdar, lambdab;
	double dl = 0;
	for (int i = 0; i < 3; i++) 
		dl+= pow((basepos[i] - state_x[i]),2);
	dl = sqrt(dl);
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	if (RTKmode == "PL") {
		/*SatSat_PdresG.clear();
		SatSat_LdresG.clear();
		SatSat_PdresB.clear();
		SatSat_LdresB.clear();*/
		if (SATrefGPSID) {
			Station_OBS refrover = GPSobs[SATrefGPSID].get_rover();
			Station_OBS refbase = GPSobs[SATrefGPSID].get_base();
			refroverP = refrover.get_P_zres();
			refroverL = refrover.get_L_zres();
			refbaseP = refbase.get_P_zres();
			refbaseL = refbase.get_L_zres();
			GPSobs[SATrefGPSID].get_lambda(0, &lambdab, &lambdar);
			GPSobs[SATrefGPSID].calsdresvar(dl);
		}

		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		for (; itg != GPSobs.end()&& SATrefGPSID; itg++) {
			if (itg->second.get_BR_unusedflag()||itg->first == SATrefGPSID)
				continue;
			Station_OBS currover = itg->second.get_rover();
			Station_OBS curbase = itg->second.get_base();
			double pseudres = (refroverP - refbaseP) - (currover.get_P_zres()- curbase.get_P_zres());
			double carrierdres = (refroverL - refbaseL) - (currover.get_L_zres() - curbase.get_L_zres());
			int refindex = SATrefGPSID -1 + (dynamics_on ? 9 : 3);
			int curindex = itg->first - 1 + (dynamics_on ? 9 : 3);
			if (numiterator >= debug_stop) {
				int aaaa = 1;
			}
			carrierdres -= lambdar*(state_x[refindex]- state_x[curindex]);
			if (fabs(pseudres) > TC_conf.maxinno && fabs(carrierdres) > TC_conf.maxinno)
				continue;
			if (fabs(pseudres) <= TC_conf.maxinno) {
				int mapkey = itg->first;
				SatSat_PdresG.insert(make_pair(mapkey, pseudres));
			}
			if (fabs(carrierdres) <= TC_conf.maxinno) {
				int mapkey = itg->first;
				SatSat_LdresG.insert(make_pair(mapkey, carrierdres));
			}
			itg->second.calsdresvar(dl);
			
				
		}
		if (SATrefBDSID) {
			Station_OBS refrover = BDSobs[SATrefBDSID].get_rover();
			Station_OBS refbase = BDSobs[SATrefBDSID].get_base();
			refroverP = refrover.get_P_zres();
			refroverL = refrover.get_L_zres();
			refbaseP = refbase.get_P_zres();
			refbaseL = refbase.get_L_zres();
			BDSobs[SATrefBDSID].get_lambda(0, &lambdab, &lambdar);
			BDSobs[SATrefBDSID].calsdresvar(dl);
		}
		for (; itb != BDSobs.end() && SATrefBDSID; itb++) {
			if (itb->second.get_BR_unusedflag() || itb->first == SATrefBDSID)
				continue;
			Station_OBS currover = itb->second.get_rover();
			Station_OBS curbase = itb->second.get_base();
			double pseudres = (refroverP - refbaseP) - (currover.get_P_zres() - curbase.get_P_zres());
			double carrierdres = (refroverL - refbaseL) - (currover.get_L_zres() - curbase.get_L_zres());
			int refindex = SATrefBDSID - 1 +NSATGPS+ (dynamics_on ? 9 : 3);
			int curindex = itb->first - 1 + NSATGPS+(dynamics_on ? 9 : 3);
			carrierdres -= lambdar*(state_x[refindex] - state_x[curindex]);
			if (fabs(pseudres) > TC_conf.maxinno && fabs(carrierdres) > TC_conf.maxinno)
				continue;
			if (fabs(pseudres) <= TC_conf.maxinno) {
				int mapkey = itb->first;
				SatSat_PdresB.insert(make_pair(mapkey, pseudres));
			}
			if (fabs(carrierdres) <= TC_conf.maxinno) {
				int mapkey = itb->first;
				SatSat_LdresB.insert(make_pair(mapkey, carrierdres));
			}
			itb->second.calsdresvar(dl);
		}	
	}
	//if (!SatSat_PdresG.size())//如果没有满足条件（小于optmaxino的值，那就没办法，必须放弃）
	//	SATrefGPSID = 0;
	//if (!SatSat_PdresB.size())
	//	SATrefBDSID = 0;
	if (!SatSat_PdresG.size()&&!SatSat_LdresG.size())//如果没有满足条件（小于optmaxino的值，那就没办法，必须放弃）
		SATrefGPSID = 0;//为什么当时只考虑了伪距的？相位的没考虑
	if (!SatSat_PdresB.size()&& !SatSat_LdresB.size())
		SATrefBDSID = 0;
	return;


}

void receiver::initxP(double xvalue, double pvalue, int index) {
	RTK_x[index] = xvalue;
	for (int j = 0; j < nx; j++) {
		RTK_P[index + j * nx] = RTK_P[j + index * nx] = index == j ? pvalue : 0.0;
	}
	return;
}

void receiver::rcv_detslp_ll() {
	map<int, OBS_SAT> ::iterator itg = GPSobs.begin();
	map<int, OBS_SAT> ::iterator itb = BDSobs.begin();
	for (; itg != GPSobs.end(); itg++) {
		if (itg->second.get_BR_unusedflag())
			continue;
		slip_G[itg->first - 1] &= 0xFC;
		itg->second.detslp_ll(slip_G+itg->first-1,half_G+itg->first-1,tt);
	
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	for (; itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_unusedflag())
			continue;
		slip_B[itb->first - 1] &= 0xFC;
		itb->second.detslp_ll(slip_B + itb->first - 1, half_B + itb->first - 1, tt);
	}
}
void receiver::reset_outageSAT() {
	int reset;
	for (int i = 1; i <= NSATGPS; i++) {
		reset = ++outc_G[i - 1] > TC_conf.maxout;
		int index = i-1+(dynamics_on ? 9:3);
		if (reset) {
			if (RTK_x[index] != 0.0) {
				initxP(0.0, 0.0, index);
				outc_G[i - 1] = 0;
			}
			lock_G[i - 1] = -TC_conf.minlock;
		}
	}
	for (int i = 1; i <= NSATCMP; i++) {
		reset = ++outc_B[i - 1] > TC_conf.maxout;
		int index = NSATGPS + i - 1+(dynamics_on ? 9: 3);
		if (reset) {
			if (RTK_x[index] != 0.0) {
				initxP(0.0, 0.0, index);
				outc_B[i - 1] = 0;
			}
			lock_B[i - 1] = -TC_conf.minlock;
		}
	}
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	/*if (numiterator >= debug_stop) {
		int aaaa = 1;
	}*/
	for (; itg != GPSobs.end(); itg++) {
		int j = itg->first - 1+(dynamics_on ? 9: 3);
		RTK_P[j + j * nx] += TC_conf.prn[0] * TC_conf.prn[0] * tt;
		char slip = slip_G[itg->first - 1];
		if (!(slip & 1))
			continue;
		RTK_x[j] = 0.0;
		lock_G[itg->first - 1] = -TC_conf.minlock;
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	for (; itb != BDSobs.end(); itb++) {
		int j = NSATGPS + itb->first - 1+(dynamics_on ? 9: 3);
		RTK_P[j + j * nx] += TC_conf.prn[0] * TC_conf.prn[0] * tt;
		char slip = slip_B[itb->first - 1];
		if (!(slip & 1))
			continue;
		RTK_x[j] = 0.0;
		lock_B[itb->first - 1] = -TC_conf.minlock;
	}
	//不考虑glonass
}
void receiver::initialize_bias_state() {
	// 目前只考虑单频
	map<int, OBS_SAT >::iterator itg = GPSobs.begin();
	map<int, OBS_SAT >::iterator itb = BDSobs.begin();
	map<int,double > bias_G,bias_B;
	double offset = 0.0;
	double cp = 0.0;
	double pr = 0.0;
	int j = 0;
	for (; itg != GPSobs.end(); itg++) {
		//此处先只考虑P 和 L
		if (itg->second.get_BR_unusedflag())
			continue;
		itg->second.calP_sres();
		itg->second.calL_sres();
		cp = itg->second.getL_sres();
		pr = itg->second.getP_sres();
		double lambdab, lambdar;
		itg->second.get_lambda(0,&lambdab, &lambdar);
		if (fabs(cp) < 1E-14 || fabs(pr) < 1E-14 || fabs(lambdab) <= 1E-14 || fabs(lambdab-lambdar)>=1e-14) {
			itg->second.set_BR_unusedflag(0);
			bias_G.insert(make_pair(itg->first,0.0));
			continue;
		}
		bias_G.insert(make_pair(itg->first, cp - pr / lambdar));
		int index = itg->first - 1+(dynamics_on ? 9:3);
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		if (RTK_x[index] != 0.0) {
			offset += bias_G[itg->first] - RTK_x[index];
			j++;
		}
	}
	for (; itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_unusedflag())
			continue;
		itb->second.calP_sres();
		itb->second.calL_sres();
		cp = itb->second.getL_sres();
		pr = itb->second.getP_sres();
		double lambdab, lambdar;
		itb->second.get_lambda(0, &lambdab, &lambdar);
		if (fabs(cp) < 1E-14 || fabs(pr) < 1E-14 || fabs(lambdab) <= 1E-14 || fabs(lambdab - lambdar) >= 1e-9) {
			itb->second.set_BR_unusedflag(0);
			bias_B.insert(make_pair(itb->first, 0.0));
			continue;
		}
		bias_B.insert(make_pair(itb->first, cp - pr /lambdar));
		int index = itb->first - 1 +NSATGPS+ (dynamics_on ? 9 : 3);
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		if (RTK_x[index] != 0.0) {
			offset += bias_B[itb->first] - RTK_x[index];
			j++;
		}
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	//int j = bias_B.size() + bias_G.size();
	if (j > 0) {
		com_bias = offset / j;
		int pvanum = (dynamics_on ? 9 : 3);
		for (int i = pvanum; i < nx; i++) {
			/*if (i == 13 && numiterator >= debug_stop) {
				int aaaa = 1;
			}*/
			if (RTK_x[i] != 0.0)
				RTK_x[i] += com_bias;
		}
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	//int pvanum = (dynamics_on ? 9 : 3);
	for (itg = GPSobs.begin(); itg != GPSobs.end(); itg++) {
		if (itg->second.get_BR_unusedflag())
			continue;
		int index = itg->first - 1 + (dynamics_on ? 9 : 3);
		if (bias_G[itg->first] == 0.0 || RTK_x[index] != 0.0)
			continue;
		double lambdab, lambdar;
		itg->second.get_lambda(0, &lambdab, &lambdar);
		/*double xxxx = (bias_G[itg->first] - com_bias) / lambdar;
		int xii = index - 9+1;*/
		
		Reset_Sbias_index.push_back(index);
		initxP(bias_G[itg->first], (TC_conf.Ini_state_std[0]* TC_conf.Ini_state_std[0]), index);
		//lock_G[itg->first - 1] = -minlock;
	}
	for (itb = BDSobs.begin(); itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_unusedflag())
			continue;
		int index = itb->first - 1 +NSATGPS+ (dynamics_on ? 9 : 3);
		if (bias_B[itb->first] == 0.0 || RTK_x[index] != 0.0)
			continue;
		double lambdab, lambdar;
		itb->second.get_lambda(0, &lambdab, &lambdar);
		Reset_Sbias_index.push_back(index);
		initxP(bias_B[itb->first], (TC_conf.Ini_state_std[0] * TC_conf.Ini_state_std[0]), index);
		//lock_B[itb->first - 1] = -minlock;
	}


}
void receiver::udstate() {
	
	udpos();
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	udbias();
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
}
void receiver::udpos() {
	double var = 0.0;
	if (norm_array(RTK_x, 3) <= 0.0) {
		for (int i = 0; i < 3; i++) {
			initxP(sol.getpos()[i], TC_conf.VAR_POS,i);
		}
		if (dynamics_on) {
			for (int i = 3; i < 6; i++) {
				initxP(sol.getvel()[i-3], TC_conf.VAR_VEL, i);
			}
			for (int i = 6; i < 9; i++) {
				initxP(1E-6, TC_conf.VAR_ACC, i);
			}
		}
	}
	//非初始化，状态更新不更新位置，把P矩阵对应位置置成无穷大
	if (!dynamics_on) {
		for (int i = 0; i < 3; i++) {
			initxP(sol.getpos()[i], TC_conf.VAR_POS, i);
		}
		return;
	}
	for (int i = 0; i < 3; i++)
		var += RTK_P[i + i * nx];
	var /= 3.0;
	if (var > TC_conf.VAR_POS) {
		/* reset position with large variance */
		for (int i = 0; i < 3; i++)
			initxP(sol.getpos()[i], TC_conf.VAR_POS, i);
		for (int i = 3; i < 6; i++)
			initxP(sol.getvel()[i-3], TC_conf.VAR_VEL, i);
		for (int i = 6; i < 9; i++)
			initxP(1E-6, TC_conf.VAR_ACC, i);
		cout<<"reset rtk position due to large variance: var="<<setprecision(12)<<var<<endl;
		return;
	}
	short *ix = new short[nx];
	int ix_nx = 0;
	for (int i= 0; i < nx; i++) {
		if (i < 9 || (RTK_x[i] != 0.0&&RTK_P[i + i * nx] > 0.0))
			ix[ix_nx++] = i;
	}
	Eigen::MatrixXd F = Eigen::MatrixXd::Identity(ix_nx, ix_nx);
	Eigen::MatrixXd P = Eigen::MatrixXd::Zero(ix_nx,ix_nx);
	Eigen::MatrixXd FP = Eigen::MatrixXd::Zero(ix_nx, ix_nx);
	Eigen::Matrix3d Q = Eigen::Matrix3d::Zero();
	Eigen::VectorXd x = Eigen::VectorXd::Zero(ix_nx);
	Eigen::VectorXd xp = Eigen::VectorXd::Zero(ix_nx);
	
	for (int i = 0; i < 6; i++) {
		F(i,(i + 3)) = tt;
	}
	/*if (numiterator == 16) {
		int aaaaa = 1;
	
	}*/
	//if (var < optthresar[1]) {
		for (int i = 0; i < 3; i++) {
			F(i, (i + 6)) = pow(tt, 2) / 2.0;
		}
	//}
	//else {
	//	//cout << "pos var too high for accel term" << endl;
	//}
	for (int i = 0; i < ix_nx; i++) {//这一步感觉是把所有的非0的状态量都挪到了x这个矩阵的最前端，x的后面全是0
		x[i] =RTK_x[ix[i]];
		for (int j = 0; j < ix_nx; j++) {// 压缩的P矩阵只保留现有的状态对应的协方差
			P(i ,j) = RTK_P[ix[i] + ix[j] * nx];
		}
	}
	/*if (numiterator >= debug_stop) {
		
		matrixtofile(F, "F", 'N');
		matrixtofile(P, "P", 'N');
		matrixtofile(x, "x", 'N');
		int aaaa = 1;
	}*/
	xp = F * x;
	FP = F * P;
	P = FP * F.transpose();
	
	for (int i = 0; i < ix_nx; i++) {//把压缩的数据返回到大矩阵中
		RTK_x[ix[i]] = xp[i];
		for (int j = 0; j < ix_nx; j++) {
			RTK_P[ix[i] + ix[j] * nx] = P(i, j);
		}
	}
	Q(0, 0) = Q(1, 1) = pow(TC_conf.prn[3], 2)*fabs(tt);
	Q(2, 2) = pow(TC_conf.prn[4], 2)*fabs(tt);
	vector<double > posecef = {RTK_x[0],RTK_x[1],RTK_x[2] };
	vector<double > posllh = vector<double > (3,0);
	/*cout << "Q矩阵是:" << endl;
	cout << setprecision(12) << Q << endl;*/
	ecef2pos(posecef,posllh);
	Eigen::Matrix3d Qv;
	covecef(posllh, Q, Qv);	
	/*cout << "Qv矩阵是:" << endl;
	cout << setprecision(12) << Qv << endl;*/
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			RTK_P[i + 6 + (j + 6)*nx] += Qv(i,j);
	}
	delete ix;
	return;
}
void receiver::udbias() { // 返回
	Reset_Sbias_index.clear();
	vector<int>().swap(Reset_Sbias_index);
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	rcv_detslp_ll();//先检查周跳,单频
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	reset_outageSAT();//对应RTKLIB udbias的第二部分，但是不考虑INST模式,重置lock 等
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	initialize_bias_state();// 借助单差初始化状态量
	if (numiterator == debug_stop) {
		int aaaa = 1;
	}
	
	//matrixtofile(RTK_x_vector,"RTK_x",'N');

	/*if (numiterator >= debug_stop) {
		Eigen::MatrixXd RTKP = Eigen::MatrixXd::Zero(nx, nx);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < nx; j++) {
				RTKP(i, j) = RTK_P[i + j * nx];
			}
		}
		matrixtofile(RTKP,"RTKP",'N');
		int aaa = 1;
	}*/
	int aaaa = 1;
}
vector<int> receiver::get_GBSat_PLdres_size() {
	vector<int> v(4, 0);
	v[0] = SatSat_PdresG.size();
	v[1] = SatSat_LdresG.size();
	v[2] = SatSat_PdresB.size();
	v[3] = SatSat_LdresB.size();
	return v;
}	



void receiver::MeasureUpdate(string RTKmode) {
	int ddpairnum_GPSP = 0, ddpairnum_GPSL = 0, ddpairnum_BDSP = 0, ddpairnum_BDSL = 0,ddpairnum = 0;
	ddpairnum_GPSP = SatSat_PdresG.size();
	ddpairnum_GPSL = SatSat_LdresG.size();
	ddpairnum_BDSP = SatSat_PdresB.size();
	ddpairnum_BDSL = SatSat_LdresB.size();
	ddpairnum = (ddpairnum_GPSP+ ddpairnum_GPSL + ddpairnum_BDSP+ ddpairnum_BDSL);
	int numofSbias;
	//int numofstatepva = (dynamics_on ? 9 : 3);
	//numofSbias = (ddpairnum_GPS ? ddpairnum_GPS + 1 : 0) + (ddpairnum_BDS ? ddpairnum_BDS + 1 : 0);
	//int numofvalidstates = numofSbias + numofstatepva;
	//虽然RTK_x中不为0的单差很多，但是只有真正在当前时刻有效的卫星才参与滤波
	// RTKLIB的压缩矩阵包含所有的不为0状态量，但当前无效的卫星对应的单差对应的H矩阵相关元素为0，因此相当于没有参与观测更新。
	if (ddpairnum < 1) {
		cout << "ddres pairs are less than 1" << endl;
		return;//在RTKLIB中即使nv数目很少，也会继续进行观测更新，但是最终若有效卫星数过少，会认定定位无效
	}
	Eigen::MatrixXd H;// 双差对数*状态量维度,会考虑所有当前不为0 的状态量
	Eigen::VectorXd x;
	Eigen::VectorXd v = Eigen::VectorXd::Zero(ddpairnum);//y-h(x)
	Eigen::MatrixXd P;
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(ddpairnum, ddpairnum);
	vector<int > ix;
	FormHxP(RTKmode,H,x,P, ddpairnum_GPSP, ddpairnum_GPSL, ddpairnum_BDSP, ddpairnum_BDSL,ix);
	FormvR(RTKmode,v,R, ddpairnum_GPSP, ddpairnum_GPSL, ddpairnum_BDSP, ddpairnum_BDSL);
	if (numiterator >= debug_stop) {
		/*matrixtofile(H, "MU_H", 'T');
		matrixtofile(P, "MU_P", 'T');
		matrixtofile(R, "MU_R", 'N');
		matrixtofile(v, "MU_v", 'N');
		matrixtofile(x, "MU_x", 'N');*/
		int aaaa = 1;
	}
	filter_unzip(H, x, P, v,R, ix);
	//delete ix;
}
void receiver::FormHxP(string RTKmode, Eigen::MatrixXd &H, Eigen::VectorXd &x, Eigen::MatrixXd &P, int ddpairnum_GPSP, int ddpairnum_GPSL, int ddpairnum_BDSP, int ddpairnum_BDSL,vector<int> &ix) {
	//此处ix 包含了 位速加 所有当前不为0的单差，即使对应卫星不存在双差，仍然会建立，
	map<int, double > ::iterator itddPG = SatSat_PdresG.begin();
	map<int, double > ::iterator itddLG = SatSat_LdresG.begin();
	map<int, double > ::iterator itddPB = SatSat_PdresB.begin();
	map<int, double > ::iterator itddLB = SatSat_LdresB.begin();
	int numofstatepva = (dynamics_on ? 9 : 3);
	int refGPSindex = -1;
	int refBDSindex = -1;
	for (int i = 0; i < nx; i++) { //寻找参考卫星对应在ix中的位置为后面设置H矩阵做准备
		if (i < numofstatepva) {
			ix.push_back(i);
			continue;
		}
		if (RTK_x[i] != 0.0 && RTK_P[i + i * nx] > 0.0) {
			ix.push_back(i);
			if (i - numofstatepva + 1 == SATrefGPSID)
				refGPSindex = ix.size()-1;
			if (i - numofstatepva + 1 - NSATGPS == SATrefBDSID)
				refBDSindex = ix.size() - 1;
		}
	}
	
	int numofstate = ix.size();
	int numofstateSres = numofstate - numofstatepva;
	int ddpairnum_LPLP = ddpairnum_GPSP+ ddpairnum_GPSL + ddpairnum_BDSP+ ddpairnum_BDSL;
	H = Eigen::MatrixXd::Zero(ddpairnum_LPLP, numofstate);// 双差对数*状态量维度
	x = Eigen::VectorXd::Zero(numofstate);
	P = Eigen::MatrixXd::Zero(numofstate, numofstate);
	int refsatID;
	vector<double> refsat_razellos;
	vector<double> cursat_razellos;
	int phaseindex = 0;
	int codeindex = ddpairnum_GPSL;
	
	//int refxindex = numofstatepva;
	double lambdar, lambdab;
	
	//int index = 0;
	if (RTKmode == "PL") { // 如果模式时PL，那么此时PG和PB的索引应该是完全一致的
		for ( int i = 0; i < numofstate; i++) {
			x(i) = RTK_x[ix[i]];
		}
		// 无论如何，真正解算的P是状态量之间的协方差
		// 真正要做的就是把每个状态量自身的自相关和彼此的互相关从大P中提取出来
		// 滤波中 位速加 相对 单差 也会产生协方差，不能省略这部分
		// 没有双差的量，只要单差存在，就也加入滤波，虽然对应H元素是0，对其它状态量影响小，但是会对自身有影响
		// 似乎不这样做，时间久了会影响LD分解
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		
		}

		//int ixindex = 0;
		int ixindex = numofstatepva;//目的是定位对应双差的另一颗GPS星在ix中的位置,哪怕只有一对伪距（相位）双差，也要算
		if (itddPG != SatSat_PdresG.end()|| itddLG != SatSat_LdresG.end()) {
			refsatID = SATrefGPSID;
			refsat_razellos = GPSobs[refsatID].get_rover().get_razellos();
			//P(xindex, xindex) = RTK_P[numofstatepva + refsatID - 1 + (numofstatepva + refsatID - 1)*nx];
			//refxindex = xindex;
			GPSobs[refsatID].get_lambda(0, &lambdab, &lambdar);
			//ixindex++;
			
		}
		//H矩阵，先填GPS 相位
		for (; itddLG != SatSat_LdresG.end(); itddLG++) {
			int cursatID = itddLG->first;
			cursat_razellos = GPSobs[cursatID].get_rover().get_razellos();
			for (int i = 0; i < 3; i++) {
				H(phaseindex, i) = cursat_razellos[i + 3] - refsat_razellos[i + 3];
			}
			H(phaseindex, refGPSindex) = lambdar;
			if (ix[ixindex] + 1 - numofstatepva == SATrefGPSID)//和 map 一样，都是按顺序从小到大存的
				ixindex++;
			while (ix[ixindex] + 1 - numofstatepva != cursatID)//跳过没有观测量的单差状态量
				ixindex++;
			H(phaseindex, ixindex) = -lambdar;
			ixindex++;
			phaseindex++;
		}
		// H矩阵，再填GPS 伪距
		for (; itddPG != SatSat_PdresG.end(); itddPG++) {
			int cursatID = itddPG->first;
			cursat_razellos = GPSobs[cursatID].get_rover().get_razellos();
			for (int i = 0; i < 3; i++) {
				H(codeindex, i) = cursat_razellos[i + 3] - refsat_razellos[i + 3];
			}
			codeindex++;
		}
		phaseindex = ddpairnum_GPSP+ ddpairnum_GPSL;
		codeindex = phaseindex + ddpairnum_BDSL;
		if (itddPB != SatSat_PdresB.end()|| itddLB != SatSat_LdresB.end()) {//目的是定位对应双差的另一颗BDS星在ix中的位置,哪怕只有一对伪距（相位）双差，也要算
			refsatID = SATrefBDSID;
			refsat_razellos = BDSobs[refsatID].get_rover().get_razellos();
			//P(xindex, xindex) = RTK_P[numofstatepva +NSATGPS+ refsatID - 1 + (numofstatepva + NSATGPS  + refsatID - 1)*nx];
			/*refxindex = xindex;*/
			BDSobs[refsatID].get_lambda(0, &lambdab, &lambdar);
			//xindex++;
		}
		// H矩阵，填完GPS 填BDS 相位
		for (; itddLB != SatSat_LdresB.end(); itddLB++) {
			int cursatID = (itddLB->first);
			cursat_razellos = BDSobs[cursatID].get_rover().get_razellos();
			for (int i = 0; i < 3; i++) {
				H(phaseindex, i) = cursat_razellos[i + 3] - refsat_razellos[i + 3];
			}
			H(phaseindex, refBDSindex) = lambdar;
			if (ix[ixindex] + 1 - numofstatepva - NSATGPS == SATrefBDSID)//和 map 一样，都是按顺序从小到大存的
				ixindex++;
			while (ix[ixindex] + 1 - numofstatepva - NSATGPS != cursatID)//跳过没有观测量的单差状态量
				ixindex++;
			H(phaseindex, ixindex) = -lambdar;
			ixindex++;
			phaseindex++;
		}
		// H矩阵，再填BDS伪距
		for (; itddPB != SatSat_PdresB.end(); itddPB++) {
			int cursatID = (itddPB->first);
			cursat_razellos = BDSobs[cursatID].get_rover().get_razellos();
			for (int i = 0; i < 3; i++) {
				H(codeindex, i) = cursat_razellos[i + 3] - refsat_razellos[i + 3];
			}
			codeindex++;
		}
		for (int i = 0; i < ix.size(); i++) {
			for (int j = 0; j < ix.size(); j++) {
				P(i, j) = RTK_P[ix[i] + ix[j] * nx];
			}
		}
	}
	return;
}
void receiver::FormvR(string RTKmode, Eigen::VectorXd &y, Eigen::MatrixXd &R, int ddpairnum_GPSP, int ddpairnum_GPSL, int ddpairnum_BDSP, int ddpairnum_BDSL) {
	map<int, double > ::iterator itddPG = SatSat_PdresG.begin();
	map<int, double > ::iterator itddLG = SatSat_LdresG.begin();
	map<int, double > ::iterator itddPB = SatSat_PdresB.begin();
	map<int, double > ::iterator itddLB = SatSat_LdresB.begin();
	int yindex = 0;
	int refsatID,cursatID;
	int phaseindex = 0;
	int codeindex = ddpairnum_GPSL ;
	double RrefP, RrefL, RcurP, RcurL;
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	// 可能只有双差伪距而没有双差相位？
	if (RTKmode == "PL") { // 如果模式时PL，那么此时PG和PB的索引应该是完全一致的
		if (itddPG != SatSat_PdresG.end()|| itddLG != SatSat_LdresG.end()) {
			refsatID = SATrefGPSID;
			RrefP = GPSobs[refsatID].getP_svar();
			RrefL = GPSobs[refsatID].getL_svar();
		}
		// GPS 相位
		for (; itddLG != SatSat_LdresG.end(); itddLG++) {
			cursatID = (itddLG->first);
			y[phaseindex] = SatSat_LdresG[itddLG->first];
			RcurL = GPSobs[cursatID].getL_svar();
			R(phaseindex, phaseindex) = RrefL + RcurL;
			phaseindex++;
		}
		for (int i = 0; i < ddpairnum_GPSL; i++) {
			for (int j = 0; j < ddpairnum_GPSL; j++) {
				if (i != j)
					R(i, j) = RrefL;
			}
		}
		// GPS 伪距
		for (; itddPG != SatSat_PdresG.end(); itddPG++) {
			cursatID = (itddPG->first);
			y[codeindex] = itddPG->second;
			RcurP = GPSobs[cursatID].getP_svar();
			R(codeindex, codeindex) = RrefP + RcurP;
			codeindex++;
		}
		for (int i = ddpairnum_GPSL; i < ddpairnum_GPSL+ ddpairnum_GPSP; i++) {
			for (int j = ddpairnum_GPSL; j < ddpairnum_GPSL + ddpairnum_GPSP; j++) {
				if (i != j)
					R(i, j) = RrefP;
			}
		}
		phaseindex = ddpairnum_GPSL+ ddpairnum_GPSP;
		codeindex = phaseindex + ddpairnum_BDSL;

		if (itddPB != SatSat_PdresB.end()|| itddLB != SatSat_LdresB.end()) {
			refsatID = SATrefBDSID;
			RrefP = BDSobs[refsatID].getP_svar();
			RrefL = BDSobs[refsatID].getL_svar();

		}
		// BDS 相位
		for (; itddLB != SatSat_LdresB.end(); itddLB++) {
			cursatID = (itddLB->first);
			y[phaseindex] = SatSat_LdresB[itddLB->first];
			RcurL = BDSobs[cursatID].getL_svar();
			R(phaseindex, phaseindex) = RrefL + RcurL;
			phaseindex++;
		}
		for (int i = ddpairnum_GPSP+ ddpairnum_GPSL; i < ddpairnum_GPSP+ ddpairnum_GPSL+ddpairnum_BDSL; i++) {
			for (int j = ddpairnum_GPSP+ ddpairnum_GPSL; j < ddpairnum_GPSP + ddpairnum_GPSL + ddpairnum_BDSL; j++) {
				if (i != j)
					R(i, j) = RrefL;
			}
		}
		// BDS伪距
		for (; itddPB != SatSat_PdresB.end(); itddPB++) {
			cursatID = (itddPB->first);
			y[codeindex] = itddPB->second;
			RcurP = BDSobs[cursatID].getP_svar();
			R(codeindex, codeindex) = RrefP + RcurP;
			codeindex++;
		}
		for (int i = ddpairnum_GPSP + ddpairnum_GPSL + ddpairnum_BDSL; i < ddpairnum_GPSP + ddpairnum_GPSL + ddpairnum_BDSL+ ddpairnum_BDSP; i++) {
			for (int j = ddpairnum_GPSP + ddpairnum_GPSL + ddpairnum_BDSL; j < ddpairnum_GPSP + ddpairnum_GPSL + ddpairnum_BDSL + ddpairnum_BDSP; j++){
				if (i != j)
					R(i, j) = RrefP;
			}
		}
	}
	return;
}
void receiver::filter_unzip(Eigen::MatrixXd &H, Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::VectorXd &v, Eigen::MatrixXd &R, const vector<int > &ix) {
	//注意，在P的更新中，位置速度加速度等状态量 也会和 单差模糊度产生协方差，所以在压缩P矩阵和返回P矩阵的时候，必须连带位速加的状态所在的索引一起。
	int numofstatepva = dynamics_on ? 9 : 3;
	Eigen::MatrixXd PHT,HPHT_add_R,K;
	Eigen::MatrixXd I;
	I.setIdentity(P.rows(), P.cols());
	/*if (numiterator == debug_stop) {
		int aaaa = 1;
	}*/
	PHT = P * H.transpose();
	HPHT_add_R = (H*P*H.transpose()+R);
	K = PHT * HPHT_add_R.inverse();
	/*if (numiterator == debug_stop) {
		matrixtofile(H, "H", 'T');
		matrixtofile(P, "P", 'N');
		matrixtofile(x, "x", 'N');
		matrixtofile(v, "v", 'N');
		int aaaa = 1;
	}*/
	//
	x = x + K * v;
	
	//matrixtofile(P, "P", 'N');
	P = (I - K * H)*P;
	//matrixtofile(P, "P+", 'N');
	/*以下把当前的压缩过的xP填回到对应的大xP中*/
	//int numofstatepva = (dynamics_on ? 9 : 3);
	//for (int i = 0; i < numofstatepva; i++) {
	//	RTK_x[i] = x[i];
	//	for (int j = 0; j < numofstatepva; j++)
	//		RTK_P[i + j * nx] = P(i, j);
	//}
	for (int i = 0; i < ix.size(); i++) {
		RTK_x[ix[i]] = x[i];
		for (int j = 0; j < ix.size(); j++) {
			RTK_P[ix[i] + ix[j] * nx] = P(i,j);
		}
	}
	return;
}
int receiver::solve_validpos(string RTKmode,double * state,int mode) {//RTKLIB 在做完卡尔曼滤波后先要对残余进行一次判定，但并不影响值，只是提示。同时会对有效卫星的outc进行清零操作
	RCV_zdres('R', state);//还包括对卫星数的判断
	int nv_dim = RCV_ddres(RTKmode, state);//注意是算上伪距 相位 两种观测量的双差对数，
	double thres = 4.0;
	
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	if (valpos(thres)) {
		if (mode) {
			lock_solve();
			out_solve();
		}
		//map<int, double>::iterator itg = SatSat_PdresG.begin();
		//map<int, double>::iterator itgL = SatSat_LdresG.begin();
		//map<int, double>::iterator itb = SatSat_PdresB.begin();
		//map<int, double>::iterator itbL = SatSat_LdresB.begin();
		////int count = 0;
		//if (itg != SatSat_PdresG.end()|| itgL!= SatSat_LdresG.end()) {
		//	outc_G[SATrefGPSID - 1] = 0;
		//}
		//for (; itg != SatSat_PdresG.end(); itg++) {
		//	//count++;
		//	outc_G[itg->first - 1] = 0;
		//}
		//for (; itgL != SatSat_LdresG.end(); itgL++) {
		//	//count++;
		//	outc_G[itgL->first - 1] = 0;
		//}
		//if (itb != SatSat_PdresB.end() || itbL!= SatSat_LdresB.end()) {
		//	outc_B[SATrefBDSID - 1] = 0;
		//}
		//for (; itb != SatSat_PdresB.end(); itb++) {
		//	//count++;
		//	outc_B[itb->first - 1] = 0;
		//}
		//for (; itbL != SatSat_LdresB.end(); itbL++) {
		//	//count++;
		//	outc_B[itbL->first - 1] = 0;
		//}
		int Satnum_thres = TC_conf.ins_mode ? 3 : 4;
		if (nv_dim / 2 + (SATrefGPSID ? 1 : 0) + (SATrefBDSID ? 1 : 0) < Satnum_thres) {//能产生双差，且双差残余比较小的
			return 0;
		}
		return 1;

	}
	return 0;
}
int receiver::valpos(double thres) {
	double fact = thres * thres;
	double RrefP, RrefL, RcurP, RcurL;
	map<int, double> ::iterator itPG = SatSat_PdresG.begin();
	//map<int, double> ::iterator itLG = SatSat_LdresG.begin();
	map<int, double> ::iterator itPB = SatSat_PdresB.begin();
	//map<int, double> ::iterator itLB = SatSat_LdresB.begin();
	if (SATrefGPSID) {
		RrefL = GPSobs[SATrefGPSID].getL_svar();
		RrefP = GPSobs[SATrefGPSID].getP_svar();
		for (; itPG != SatSat_PdresG.end(); itPG++) {
			double Pres2 = itPG->second*itPG->second;
			double Lres2 = pow(SatSat_LdresG[itPG->first], 2);
			RcurL = GPSobs[itPG->first].getL_svar();
			RcurP = GPSobs[itPG->first].getP_svar();
			double factRP = (RcurP + RrefP)*fact;
			double factRL = (RcurL + RrefL)*fact;
			if (Pres2 >= factRP) {
				cout << "G" << SATrefGPSID << "-" << "G" << itPG->first << ": ddres of Code is too large;" << endl;
			}
			if (Lres2 >= factRL) {
				cout << "G" << SATrefGPSID << "-" << "G" << itPG->first << ": ddres of Phase is too large;" << endl;
			}
		}
	}
	
	if (SATrefBDSID) {
		RrefL = BDSobs[SATrefBDSID].getL_svar();
		RrefP = BDSobs[SATrefBDSID].getP_svar();
		for (; itPB != SatSat_PdresB.end(); itPB++) {
			double Pres2 = itPB->second*itPB->second;
			double Lres2 = pow(SatSat_LdresB[itPB->first], 2);
			RcurL = BDSobs[itPB->first].getL_svar();
			RcurP = BDSobs[itPB->first].getP_svar();
			double factRP = (RcurP + RrefP)*fact;
			double factRL = (RcurL + RrefL)*fact;
			if (Pres2 >= factRP) {
				cout << "C" << SATrefBDSID << "-" << "C" << itPB->first << ": ddres of Code is too large;" << endl;
			}
			if (Lres2 >= factRL) {
				cout << "C" << SATrefBDSID << "-" << "C" << itPB->first << ": ddres of Phase is too large;" << endl;
			}
		}
	}
	
	return 1;
}
/*与模糊度固定有关的函数*/
int receiver::manage_amb_LAMBDA(double *rtk_fixX,vector<int> & ix) {
	int exc_satsID,exc_sats_lock;//单频
	vector<int > arsatsg, arsatsb;
	int excflag = 0;
	int nb = 0;
	int rerun = 0;
	float ratio1;
	if (numiterator == debug_stop) {

		//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}
	if (TC_conf.exec_mode) {
		excflag = exclude_different_sat(exc_satsID, exc_sats_lock);
	}
	//if (holdamb_flag) {
	if (numiterator >= debug_stop) {
		//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}
	nb = resamb_LAMBDA(rtk_fixX,ix);
	ratio1= sol.get_ratio();
	if (arfilter) {
		rerun = 0;
		if (nb >= 0 && ((sol.get_ratio() < TC_conf.thresar[0] * 1.1 && sol.get_ratio() < sol.get_prev_ratio1() / 2.0) || // 一般情况下，prev_ratio1 和 prev_ratio2 都会设置成
			(sol.get_ratio() < sol.get_thres() && sol.get_prev_ratio2() >= sol.get_thres()))) {
			rerun = exclude_new_sat();
		}
		if (rerun) {
			nb = resamb_LAMBDA(rtk_fixX,ix);
		}	
	}
	sol.set_prev_ratio1(ratio1);
	//}
	/*else {
		ratio1 = 0;
		nb = 0;
	}*/
	float a = ratio1 > 0 ? ratio1 : sol.get_ratio();
	sol.set_prev_ratio1(a);
	sol.set_prev_ratio2(sol.get_ratio());
	if (TC_conf.exec_mode) {
		if (excflag &&sol.get_ratio() < sol.get_thres()) {
			recover_different_sat(exc_satsID, exc_sats_lock);
		}
	}
	return nb;
}
int receiver::exclude_different_sat(int &exc_satsID,int &exc_sats_lock) {
	//可能会有问题因为vsat 对应的是贡献了有效双差残差的卫星
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	int valid_num = 0;
	for (; itg != GPSobs.end(); itg++) {
		if (itg->second.get_BR_fdemask() || itg->second.get_BR_unusedflag() || lock_G[itg->first - 1] <= 0)
			continue;
		valid_num++;
	}
	for (; itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_fdemask() || itb->second.get_BR_unusedflag() || lock_B[itb->first - 1] <= 0)
			continue;
		valid_num++;
	}
	if (sol.get_prev_ratio2() > sol.get_thres() || valid_num < TC_conf.optmindropsats) {
		exc_satsID = 0;
		exc_sats_lock = 0;
		return 0;
	}
	int i = 0;
	for (itg = GPSobs.begin(); itg != GPSobs.end(); itg++) {
		if (itg->second.get_BR_fdemask() || itg->second.get_BR_unusedflag()||lock_G[itg->first-1]<=0)
			continue;
		if (i == excsatindex) {
			exc_satsID = itg->first;
			exc_sats_lock = lock_G[itg->first - 1];
			lock_G[itg->first - 1] = -TC_conf.minfix;
			return 1;
		}
		i++;
	}
	for (itb = BDSobs.begin(); itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_fdemask() || itb->second.get_BR_unusedflag() || lock_B[itb->first - 1] <= 0)
			continue;
		if (i == excsatindex) {
			exc_satsID = itb->first+40;
			exc_sats_lock = lock_B[itb->first - 1];
			lock_B[itb->first - 1] = -TC_conf.minfix;
			return 1;
		}
		i++;
	}
	excsatindex = 0;//超出当前的星数目
	return 1;
}
void receiver::recover_different_sat(int  exc_satsID, int exc_sats_lock) {
	 
	if (exc_satsID < 40) {
		lock_G[exc_satsID - 1] = exc_sats_lock;
	}
	else {
		exc_satsID = exc_satsID % 40;
		lock_B[exc_satsID - 1] = exc_sats_lock;
	}
	
}
int receiver::validsatnum() {
	map<int, OBS_SAT> ::iterator itg = GPSobs.begin();
	map<int, OBS_SAT> ::iterator itb = BDSobs.begin();
	int ns = 0;
	for (; itg != GPSobs.end(); itg++) {
		if (!itg->second.get_BR_unusedflag() && !itg->second.get_BR_fdemask())
			ns++;
	}
	for (; itb != BDSobs.end(); itb++) {
		if (!itb->second.get_BR_unusedflag() && !itb->second.get_BR_fdemask())
			ns++;
	}
	return ns;
}
int receiver::exclude_new_sat() {
	//只考虑单频
	int rerun = 0;
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	int dly = 2;
	for (; itg != GPSobs.end(); itg++) {
		if (itg->second.get_BR_unusedflag() || itg->second.get_BR_fdemask() || itg->second.get_fix() != 2)
			continue;
		if (lock_G[itg->first - 1] == 0) {
			lock_G[itg->first - 1] = -TC_conf.minlock - dly;
			dly += 2;
			rerun = 1;
		}
	}
	for (; itb != BDSobs.end(); itb++) {
		if (itb->second.get_BR_unusedflag() || itb->second.get_BR_fdemask() || itb->second.get_fix() != 2)
			continue;
		if (lock_B[itb->first - 1] == 0) {
			lock_B[itb->first - 1] = -TC_conf.minlock - dly;
			dly += 2;
			rerun = 1;
		}
	}
	return rerun;
}
int receiver::resamb_LAMBDA(double * RTK_fixX,vector<int > & ix) {
	double var = 0.0;
	Eigen::MatrixXd D,P;
	Eigen::VectorXd x;
	//vector<int> ix;
	ix.clear();
	vector<int>().swap(ix);
	int nb;
	if (numiterator>= debug_stop) {

		//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}
	for (int i = 0; i < 3; i++)
		var += RTK_P[i + i * nx];
	var = var / 3.0;
	/*if (numiterator == debug_stop ) {
		arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}*/
	/*if (var > optthresar[1]) {
		cout << "KF var too large, cannot to IAR. var=mean(P.diag(1,2,3)):Var = " << setprecision(12) << var << endl;
		return 0;
	}*/
	if ((nb = ddmat_FormxP(D, x, P, ix)) < TC_conf.optminfixsats - 1) {
		cout << "not enough valid double-differences " << endl;
		return -1;
	}
	/*if (numiterator == debug_stop) {
		arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}*/
	int cny = D.rows();
	int cnx = D.cols();
	int na = dynamics_on ? 9 : 3;
	int solnum = 2;
	int info = 0;
	
	Eigen::VectorXd y(cny);
	Eigen::MatrixXd Qy(cny, cny);
	Eigen::MatrixXd DP(cny, cnx);
	Eigen::MatrixXd b(nb, solnum);
	Eigen::VectorXd db(nb, 1);
	Eigen::MatrixXd Qb(nb, nb);
	Eigen::MatrixXd Qab(na, nb);
	Eigen::MatrixXd QQ(na, nb);
	/*if (numiterator == debug_stop) {
		arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}*/
	y = D * x;
	if (numiterator == debug_stop) {
		//matrixtofile(D,"D",'N');
		int aaaa = 1;
	}
	DP = D * P;
	Qy = DP * D.transpose();
	Qb = Qy.block(na,na,nb,nb);
	Qab = Qy.block(0, na, na, nb);
	if (numiterator == debug_stop) {
		//cout << Qy.block(na, na, nb, nb) << endl;
		//matrixtofile(Qb,"Qb",'N');
		//matrixtofile(Qab, "Qab", 'N');
		int aaaa = 1;
	}
	Eigen::VectorXd yb = y.segment(na,nb);
	vector<double> s(2, 0);
	int i, j;
	if (numiterator == debug_stop) {
		//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}
	if (!(info = lambda(nb, solnum, yb, Qb, b, s))) {
		/*if (numiterator == debug_stop) {
			arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
			int aaaa = 1;
		}*/
		float aa  = s[0]>0?(float)(s[1] / s[0]) : 0.0f;
		if (aa > 999.9)
			aa = 999.9f;
		sol.set_ratio(aa);
		sol.set_thres(TC_conf.thresar[0]);
		if (s[0] <= 0.0 || s[1] / s[0] >= sol.get_thres()) {
			for (i = 0; i < na; i++) {
				RTK_fixXpva[i] = RTK_x[i];//先把当前浮点解作为固定解初值
				for (j = 0; j < na; j++)
					RTK_fixPpva(i,j) = RTK_P[i + j * nx];// 先把当前浮点解协方差作为固定解协方差初值
			}
			//arraymattofile(RTK_P,nx,nx,0,0,na-1,na-1,'F',"P-a");
			//先删去bias
			//matrixtofile(RTK_fixPpva, "Pa", 'N');
			y.segment(na, nb) = y.segment(na, nb)-b.block(0,0,nb,1);// y是原始位置状态量加D矩阵转换的双差状态量，这里是用浮点双差减去整数双差 manual 式 E.7.19
			yb = y.segment(na, nb);
			db = Qb.inverse()*yb;
			RTK_fixXpva = RTK_fixXpva - Qab * db;
			QQ = Qab * Qb.inverse();
			RTK_fixPpva = RTK_fixPpva - QQ * Qab.transpose();
			//if (numiterator == debug_stop) {
			//	/*matrixtofile(QQ,"QQ",'N');
			//	matrixtofile(RTK_fixPpva, "Pa", 'N');*/
			//	int aaaa = 1;
			//}
			restamb(RTK_fixX,ix,b,0);
			if (numiterator == debug_stop) {
				int aaaa = 1;
			}
		}
		else {
			nb = 0;
		}
	}
	else {
		//nb = 0;
		cout << "lambda error info=" << info << endl;
	}
	//outddint(info, ix,D,b);
	return nb;
	 
}
int receiver::ddmat_FormxP(Eigen::MatrixXd & D, Eigen::VectorXd &x, Eigen::MatrixXd & P, vector<int> &ix) {
	int GPSfixddresnum = 0;
	int BDSfixddresnum = 0;
	/*int refG = 0;
	int refB = 0;*/
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	ix = find_joinfixsat(GPSfixddresnum, BDSfixddresnum);//有可能找出来的索引包含
	if (GPSfixddresnum + BDSfixddresnum < TC_conf.optminfixsats - 1)
		return GPSfixddresnum + BDSfixddresnum;
	/*if (!GPSfixddresnum)
		SATrefGPSID = 0;
	if (!BDSfixddresnum)
		SATrefBDSID = 0;*/
	int pvanum = dynamics_on ? 9 : 3;
	int Drows = GPSfixddresnum+ BDSfixddresnum + pvanum;
	int Dcolumns = pvanum+(GPSfixddresnum? GPSfixddresnum + 1:0) + (BDSfixddresnum? BDSfixddresnum + 1:0);

	D = Eigen::MatrixXd::Zero(Drows,Dcolumns);
	x = Eigen::VectorXd::Zero(Dcolumns);
	P= Eigen::MatrixXd::Zero(Dcolumns, Dcolumns);
	
	for (int i = 0; i < pvanum; i++) {
		D(i, i) = 1.0;
		x(i) = RTK_x[i];
		//P(i, i) = RTK_P[i + i * nx];
	}
	int k = 0;
	if (SATrefGPSID) {
		x(pvanum + k) = RTK_x[ix[k]];//RTK_x[pvanum + SATrefGPSID - 1];
		k++;
	}
	for (int i = pvanum; i < pvanum+GPSfixddresnum; i++,k++) {
		D(i, pvanum) = 1.0;
		D(i, pvanum + k) = -1.0;
		x(pvanum + k) = RTK_x[ix[k]];//如果GPSfixddresnum为0 不可能进来
	}
	if (SATrefBDSID) {
		x(pvanum + k) = RTK_x[ix[k]];//RTK_x[pvanum + NSATGPS+SATrefBDSID - 1];
		k++;
	}
	for (int i = pvanum+ GPSfixddresnum; i < pvanum + GPSfixddresnum+ BDSfixddresnum; i++,k++) {
		D(i, pvanum + (GPSfixddresnum? GPSfixddresnum+1:0)) = 1.0;
		D(i, pvanum+k) = -1.0;
		x(pvanum + k) = RTK_x[ix[k]];
	}
	
	vector<int > allstateix;
	if (numiterator == debug_stop) {
		int aaaa = 1;
	}
	for (int i = 0; i < pvanum; i++)
		allstateix.push_back(i);
	allstateix.insert(allstateix.end(),ix.begin(),ix.end());
	for (int i = 0; i < allstateix.size(); i++) {
		for (int j = 0; j < allstateix.size(); j++) {
			P(i, j) = RTK_P[allstateix[i] + allstateix[j] * nx];
		}
	}
	if (numiterator == debug_stop) {
		//matrixtofile(D,"amb_D",'N');
		int aaaa = 1;
	}
	return GPSfixddresnum + BDSfixddresnum;

}
/*有模糊度的星由于lock 的限制，未必参加解算，所以增加此步骤，统计确定能参加模糊度固定的卫星在RTK_x中的位置
同时还要计算基准卫星*/
/*手机可能只有一个GPS卫星，所以没有有效的双差的时候，这个卫星也要丢掉，不对它做双差*/
vector<int> receiver::find_joinfixsat(int &GPSfixddresnum,int &BDSfixddresnum){
	SATrefGPSID = 0;
	SATrefBDSID = 0;//只要第一个满足的就可以
	if (numiterator == debug_stop) {
		int aaaa = 1;
	}
	GPSfixddresnum = 0;
	BDSfixddresnum = 0;
	int pvanum = dynamics_on ? 9 : 3;
	vector<int > ix;
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	if (numiterator == debug_stop) {
		int aaaa = 1;
	}
	for (; itg != GPSobs.end(); itg++) {
		if (RTK_x[itg->first + pvanum - 1] == 0.0 || itg->second.get_BR_unusedflag() || itg->second.get_BR_fdemask())
			continue;
		if (lock_G[itg->first - 1] > 0 && !(slip_G[itg->first - 1] & 2)) {//仰角应该不用判断 在上面unusedflag 就已经判断过了//
			itg->second.set_fix(2);
			if (SATrefGPSID == 0)
				SATrefGPSID = itg->first;
			else 
				GPSfixddresnum++;
			ix.push_back(itg->first - 1 + pvanum);
		}
		else
			itg->second.set_fix(1);
	}
	if (!GPSfixddresnum&&!ix.empty()) {//只有一个有效星
		GPSobs[ix[0] + 1 - pvanum].set_fix(1);//该星也不参与整周解算
		ix.erase(ix.begin());
		SATrefGPSID = 0;
	}
	int GPSixsize = ix.size();
	for (; itb != BDSobs.end(); itb++) {
		if (RTK_x[itb->first + NSATGPS + pvanum - 1] == 0.0 || itb->second.get_BR_unusedflag() || itb->second.get_BR_fdemask())
			continue;
		if (lock_B[itb->first - 1] > 0 && !(slip_B[itb->first - 1] & 2)) {//仰角应该不用判断 在上面unusedflag 就已经判断过了
			itb->second.set_fix(2);
			if (SATrefBDSID == 0)
				SATrefBDSID = itb->first;
			else 	
				BDSfixddresnum++;
			ix.push_back(itb->first - 1 + pvanum + NSATGPS);
		}
		else
			itb->second.set_fix(1);
	}// 寻找第一颗有效星，作为双差的基准
	if (numiterator == debug_stop) {
		int aaa = 1;
	}
	if (!BDSfixddresnum && ix.size()>GPSixsize) {//只有一个有效星
		BDSobs[ix[GPSixsize] + 1 - pvanum-NSATGPS].set_fix(1);//该星也不参与整周解算
		ix.erase(ix.end()-1);
		SATrefBDSID = 0;
	}
	return ix;
}
void receiver::restamb(double *RTK_fixX,const vector<int> & ix,  const Eigen::MatrixXd &b,int mode) { //此处ix是Sbias_ix
	int na = dynamics_on ? 9 : 3;
	int pos_dim;
	if (!mode) {
		
		pos_dim = na;
		//Sbias_index_cor = 0;
	}
	else {
		pos_dim = 6;
		//Sbias_index_cor = -15+ (dynamics_on ? 9 : 3);
	}
	for (int i = 0; i < nx; i++)
		RTK_fixX[i] = RTK_x[i];  /* init all fixed states to float state values */
	for (int i = 0; i < pos_dim; i++)
		RTK_fixX[i] = RTK_fixXpva[i];  /* overwrite non phase-bias states with fixed values */ // 把刚初始化成的浮点解转成固定解
	double refstate = 0;
	int i = 0;
	int nv = 0;
	for (i = 0; i < ix.size()&&ix[i]+1-na<NSATGPS; i++) {
		if (ix[i] + 1 - na == SATrefGPSID) {
			RTK_fixX[ix[i]] = RTK_x[ix[i]];
			refstate = RTK_x[ix[i]];
		}
		else  {
			//RTK_fixX[ix[i]] = RTK_x[ix[i]];
			RTK_fixX[ix[i]] = refstate -b(nv++,0);
		}
	}
	for (; i < ix.size(); i++) {
		if (ix[i] + 1 - na - NSATGPS == SATrefBDSID) {
			RTK_fixX[ix[i]] = RTK_x[ix[i]];
			refstate = RTK_x[ix[i]];
		}
		else {
			//RTK_fixX[ix[i]] = RTK_x[ix[i]];
			RTK_fixX[ix[i]] = refstate  - b(nv++,0);
		}
	}
	return;
}
void receiver::holdamb(double *RTK_xfix, const vector<int > &ddres_stateindex, const int nb, int firstvalidGPSID, int firstvalidBDSID) {
	//虽然观测量只有双差，并且完全可以由单差计算出来，但是构建滤波器的时候仍然要把pva建立出来，哪怕对应的H部分都是0
	// 同时只要不为0的单差状态，哪怕当前该星不可用导致没有对应的双差，仍然要加入滤波器滤波
	//这里保持了和RTKLIB一致
	//firstvalidGPSID = SATrefGPSID;
	//firstvalidBDSID = SATrefBDSID;
	Eigen::MatrixXd P, R, H;
	Eigen::VectorXd x, v;
	vector<int > all_stateindex;
	int nv = holdamb_Form_xPRHv(RTK_xfix, ddres_stateindex, nb, firstvalidGPSID,firstvalidBDSID,all_stateindex, x, H, v, P, R);
	if (nv <= 0) return;
	holdamb_Filter_FB(x,H,v,P,R, all_stateindex);
	return;
}

int receiver::holdamb_Form_xPRHv(double *RTK_xfix, const vector<int > &ddres_stateindex, const int nb, int firstvalidGPSID, int firstvalidBDSID,
	vector<int > &all_stateindex, 
	Eigen::VectorXd &x, Eigen::MatrixXd &H, Eigen::VectorXd & v, Eigen::MatrixXd & P,Eigen::MatrixXd &R) {
	int nv = 0;
	int pvanum = dynamics_on ? 9 : 3;
	int phasestatenum = ddres_stateindex.size();
	//vector<int > all_stateindex;
	/*int firstvalidGPSID = SATrefGPSID;
	int firstvalidBDSID = SATrefBDSID;*/
	int refGPSindex = -1;
	int refBDSindex = -1;
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	for (int i = 0; i < nx; i++) {
		if (i < pvanum) {
			all_stateindex.push_back(i);
			continue;
		}
		if (RTK_x[i] != 0.0 && RTK_P[i + i * nx] > 0.0) {
			/*if (numiterator >= 229 && numiterator <= 240 && i + 1 - pvanum == 5) {
				continue;
			}*/
			all_stateindex.push_back(i);
			if (i + 1 - pvanum == firstvalidGPSID)
				refGPSindex = all_stateindex.size() - 1;
			if (i + 1 - pvanum - NSATGPS == firstvalidBDSID)
				refBDSindex = all_stateindex.size() - 1;
		}
	} // 这样得到的all_stateindex 是包含了所有不为0的状态量，这其中有可能有的星并不会参与构建双差
	//关键要在H矩阵中把指定的设置成0.
	//  ddmat中，有个set_fix的操作，只有满足条件（无周跳，仰角）才会set_fix 3 才能参与构建双差伪观测量
	// 我们的结构中，没有用到set_fix而是把所有参加的都放到ddres_stateindex这个数组中
	//不能直接把所有状态不为0的星都拿来构建伪观测量
	// 但是，在滤完波校正的时候，产生的校正量，必须施加于包括（未参与构建双差伪观测量但当前状态不为0的星）
	// 就像必须施加于位速加上一样
	//但是这意味着，参与滤波的x要保留针对这个 未参与构建双差伪观测量但当前状态不为0 的状态
	// 只不过相应的H矩阵是0
	int all_statenum = all_stateindex.size();
	v = Eigen::VectorXd::Zero(nb);
	H = Eigen::MatrixXd::Zero(nb, all_statenum);
	x = Eigen::VectorXd::Zero(all_statenum);
	//Eigen::MatrixXd H_ = Eigen::MatrixXd::Zero(nb, pvanum+phasestatenum);
	//Eigen::VectorXd x_ = Eigen::VectorXd::Zero(all_statenum);
	R = Eigen::MatrixXd::Zero(nb, nb);
	P = Eigen::MatrixXd::Zero(all_statenum, all_statenum);

	//int na = dynamics_on ? 9 : 3;
	double reffixstate = 0;
	double reffloatstate = 0;
	double curfixstate = 0;
	double curfloatstate = 0;
	int i = 0;
	//int index = 0;
	int k = pvanum;//找到双差索引在单差索引矩阵中的位置
	for (i = 0; i < ddres_stateindex.size() && ddres_stateindex[i] + 1 - pvanum <= NSATGPS; i++) {
		while (all_stateindex[k] != ddres_stateindex[i]) {
			k++;
		}
		if (ddres_stateindex[i] + 1 - pvanum == firstvalidGPSID) {
			reffloatstate = RTK_x[ddres_stateindex[i]];
			reffixstate = RTK_xfix[ddres_stateindex[i]];
			//index = i+ pvanum;
		}
		else {
			curfloatstate = RTK_x[ddres_stateindex[i]];
			curfixstate = RTK_xfix[ddres_stateindex[i]];
			v(nv) = (reffixstate - curfixstate) - (reffloatstate - curfloatstate);
			//H(nv, index) = 1;
			H(nv, refGPSindex) = 1;
			H(nv, k) = -1;
			R(nv, nv) = TC_conf.varholdamb;
			nv++;
		}
		k++;
		GPSobs[ddres_stateindex[i] + 1 - pvanum].set_fix(3);
		//x[i] = RTK_x[ix[i]];
	}
	//matrixtofile(H, "hold_H", 'N');
	k--;
	for (; i < ddres_stateindex.size(); i++) {
		while (all_stateindex[k] != ddres_stateindex[i]) {
			k++;
		}
		if (ddres_stateindex[i] + 1 - pvanum - NSATGPS == firstvalidBDSID) {
			reffloatstate = RTK_x[ddres_stateindex[i]];
			reffixstate = RTK_xfix[ddres_stateindex[i]];
			//index = i+ pvanum;
		}
		else {
			curfloatstate = RTK_x[ddres_stateindex[i]];
			curfixstate = RTK_xfix[ddres_stateindex[i]];
			v(nv) = (reffixstate - curfixstate) - (reffloatstate - curfloatstate);
			//H(nv, index) = 1;
			H(nv, refBDSindex) = 1;
			H(nv, k) = -1;
			R(nv, nv) = TC_conf.varholdamb;
			nv++;
		}
		k++;
		BDSobs[ddres_stateindex[i] + 1 - pvanum - NSATGPS].set_fix(3);
		//x[i] = RTK_x[ix[i]];
	}

	for (int i = 0; i < all_statenum; i++) {
		x(i) = RTK_x[all_stateindex[i]];
		for (int j = 0; j < all_statenum; j++) {
			P(i, j) = RTK_P[all_stateindex[i] + all_stateindex[j] * nx];
		}
	}
	if (numiterator >= 0) {
		/*matrixtofile(H,"hold_H",'T');
		matrixtofile(P, "hold_P", 'N');
		matrixtofile(R, "hold_R", 'N');
		matrixtofile(x, "hold_x", 'N');
		matrixtofile(v, "hold_v", 'N');*/
		int aaaa = 1;
	}
	//matrixtofile(H, "hold_H", 'N');
	/*if (nv < TC_conf.optminholdsats-1) { //新版本取消
		cout << "holdamb: not enough sats to hold ambiguity" << endl;
		return;
	}*/
	if (nv != nb) {
		cout << "why nv is not equal to nb?" << endl;
	}
	return nv;
}
void receiver::holdamb_Filter_FB(const Eigen::VectorXd &x, const Eigen::MatrixXd &H, const Eigen::VectorXd & v, const Eigen::MatrixXd & P, const Eigen::MatrixXd &R,
	const vector<int > 	&all_stateindex) {
	int all_statenum = all_stateindex.size();
	holdamb_flag = 1;
	Eigen::MatrixXd PHT, HPHT_add_R, K,P_MU;
	Eigen::MatrixXd I;
	Eigen::VectorXd x_MU;
	I.setIdentity(P.rows(), P.cols());

	PHT = P * H.transpose();
	//cout << "H" << H << endl;
	HPHT_add_R = (H*P*H.transpose() + R);
	//cout <<"R"<< R << endl;
	K = PHT * HPHT_add_R.inverse();
	//matrixtofile(K, "hold_K", 'N');
	//cout << "K" << K << endl;
	x_MU = x + K * v;
	P_MU = (I - K * H)*P;
	for (int i = 0; i < all_statenum; i++) {
		RTK_x[all_stateindex[i]] = x_MU[i];
		for (int j = 0; j < all_statenum; j++) {
			RTK_P[all_stateindex[i] + all_stateindex[j] * nx] = P_MU(i, j);
		}
	}
}
void receiver::outddint(int info, const vector<int > & ix, const Eigen::MatrixXd & D, const Eigen::MatrixXd &b) {//输出双差模糊度
	static bool firstwrite = true;
	string ddfilename = "../solfile/ddint.txt";
	ofstream ddintout;
	if (firstwrite) {
		ddintout.open(ddfilename, ios::out);
		firstwrite = false;
	}
	else {
		ddintout.open(ddfilename, ios::app);
	}
	gtime_t timestamp = sol.gettime();
	string stime;
	time2str(timestamp, stime, 3);
	ddintout << stime << endl;
	if (info) {
		ddintout << "lambda failed" << endl;
	}
	else {
		int pvanum = dynamics_on ? 9 : 3;
		bool GPSref = false;
		bool BDSref = false;
		string GPSsatstring = "G";
		string BDSsatstring = "C";
		string cursatstring = "";
		int index_b = 0;
		for (int i = 0; i < ix.size(); i++) {
			if (ix[i] < NSATGPS + pvanum) { //第一颗GPS卫星 并且是双差
				if (!GPSref) {
					int satid = ix[i] + 1 - pvanum;
					GPSsatstring += to_string(satid);
					GPSref = true;
				}
				else {
					int satid = ix[i] + 1 - pvanum;
					cursatstring = "G" + to_string(satid);
					ddintout << GPSsatstring << "-" << cursatstring << ": " << b(index_b++,0) << endl;
				}
			}
			else { // 说明是BDS 卫星
				if (!BDSref) {
					int satid = ix[i] + 1 - pvanum-NSATGPS;
					BDSsatstring += to_string(satid);
					BDSref = true;
				}
				else {
					int satid = ix[i] + 1 - pvanum-NSATGPS;
					cursatstring = "C" + to_string(satid);
					ddintout << BDSsatstring << "-" << cursatstring << ": " << b(index_b++,0) << endl;
				}
			}
		}
	}
	ddintout.close();
}
void receiver::saveslips(const OBS_RINEX & obs,int Btintflag,int Rtintflag) {
	int satID = 0;
	int LLIB = 0;
	int LLIR = 0;
	//当前时刻更新了 Btintflag 和 Rtintflag，才应该更新对应的B R 的周跳
	if (Btintflag) {
		for (int i = 0; i < obs.GPS_satID_B.size(); i++) {
			LLIB = obs.GPS_flags_B[i*obs.GPS_Fnum_B * 3 + 2];
			map<int, BRslip>::iterator it = GPSslip.find(obs.GPS_satID_B[i]);
			if (it == GPSslip.end()) { //未找到当前星
				BRslip a = BRslip();
				a.setBslip(LLIB);
				GPSslip[obs.GPS_satID_B[i]] = a;
			}
			else { //已经有当前星
				int original_LLB = GPSslip[obs.GPS_satID_B[i]].getBslip();
				GPSslip[obs.GPS_satID_B[i]].setBslip(LLIB | original_LLB);
			}
		}
		for (int i = 0; i < obs.BDS_satID_B.size(); i++) {
			LLIB = obs.BDS_flags_B[i*obs.BDS_Fnum_B * 3 + 2];
			map<int, BRslip>::iterator it = BDSslip.find(obs.BDS_satID_B[i]);
			if (it == BDSslip.end()) { //未找到当前星
				BRslip a = BRslip();
				a.setBslip(LLIB);
				BDSslip[obs.BDS_satID_B[i]] = a;
			}
			else { //已经有当前星
				int original_LLB = BDSslip[obs.BDS_satID_B[i]].getBslip();
				BDSslip[obs.BDS_satID_B[i]].setBslip(LLIB | original_LLB);
			}
		}
	}
	if (Rtintflag) {
		for (int i = 0; i < obs.GPS_satID_R.size(); i++) {
			LLIR = obs.GPS_flags_R[i*obs.GPS_Fnum_R * 3 + 2];
			map<int, BRslip>::iterator it = GPSslip.find(obs.GPS_satID_R[i]);
			if (it == GPSslip.end()) { //未找到当前星
				BRslip a = BRslip();
				a.setRslip(LLIR);
				GPSslip[obs.GPS_satID_R[i]] = a;
			}
			else { //已经有当前星
				int original_LLR = GPSslip[obs.GPS_satID_R[i]].getRslip();
				GPSslip[obs.GPS_satID_R[i]].setRslip(LLIR | original_LLR);
			}
		}
		for (int i = 0; i < obs.BDS_satID_R.size(); i++) {
			LLIR = obs.BDS_flags_R[i*obs.BDS_Fnum_R * 3 + 2];
			map<int, BRslip>::iterator it = BDSslip.find(obs.BDS_satID_R[i]);
			if (it == BDSslip.end()) { //未找到当前星
				BRslip a = BRslip();
				a.setRslip(LLIR);
				BDSslip[obs.BDS_satID_R[i]] = a;
			}
			else { //已经有当前星
				int original_LLR = BDSslip[obs.BDS_satID_R[i]].getRslip();
				BDSslip[obs.BDS_satID_R[i]].setRslip(LLIR | original_LLR);
			}
		}
	}
	//if (numiterator >= debug_stop - 1) {
	//	int aaaa = 1;
	//}
	
	
	return;
}
void receiver::restoreslips() { //在本函数完成了对LLI的更新。假设tint 1s 时间间隔0.1，那么.3s的LLI也会对1s的LLI产生影响。
	map<int, OBS_SAT>::iterator itg = GPSobs.begin();
	map<int, OBS_SAT>::iterator itb = BDSobs.begin();
	for (; itg != GPSobs.end(); itg++) {
		map<int, BRslip> ::iterator itBR = GPSslip.find(itg->first);
		if (itBR != GPSslip.end()) {
			itg->second.orBRLLI(itBR->second.getBslip(), itBR->second.getRslip());
		}
	}
	for (; itb != BDSobs.end(); itb++) {
		map<int, BRslip> ::iterator itBR = BDSslip.find(itb->first);
		if (itBR != BDSslip.end()) {
			itb->second.orBRLLI(itBR->second.getBslip(), itBR->second.getRslip());
		}
	}
	GPSslip.clear();
	BDSslip.clear();
	return;
}
//作用：读取一个历元的数据并更新 receiver结构体的GPSobs 和 BDSobs
// 比较重要的是要筛选出指定间隔的数据。(可能base 10Hz 更新，rover 1Hz更新，但是10Hz中出了周跳，周跳也应该得到保留)
//rover和base要同步。该间隔由全局变量tint指定
// 必须rover base 都支持指定间隔的更新速率，才能进行运算
// 否则，只能按低频率的那个运算，但是，高频率的在非整s中产生的周跳，会被保存下来，并参与解算
int screent_obs(ifstream &obsfileB, ifstream &obsfileR, string &buffB, string &buffR,receiver & rcv, OBS_RINEX &obs) {
	//第一次进入的时候此处buffB buffR时间对应，且一定是第一组符合tint限制的时间。
	//每一次执行本函数，或者rover 或者base 或者二者都要读取一个历元的数据
	// 两个站都满足tint的时候才往下执行，一个站满足tint的时候就不再读取，另一个站持续读取
	// 另，当接收到某时刻数据的时候，不输出本时刻的时间戳，曾造成本函数陷入死循环的bug
	int Bupdateflag = 0;
	int Rupdateflag = 0; //BR 是否发生一次读取一个历元的数据
	int Bintflag = 0;
	int Rintflag = 0;// BR是否已经都到达了满足tint约束的时间
	while (1) {
		/*if (numiterator >= debug_stop-1) {
			int aaaa = 1;
		}*/
		obs.readepoch(obsfileB, obsfileR, buffB, buffR, Bupdateflag, Rupdateflag, Bintflag, Rintflag);
		rcv.saveslips(obs, Bupdateflag, Rupdateflag);//只要有更新，周跳都要保存
		if (Bintflag && Rintflag)// 都满足tint约束才应该跳出
			break;
		if (Bupdateflag) {//若b未更新，buffB不用变。更新了才应该往下读一行？哪个站更新了，就该往下读
			if (!getline(obsfileB, buffB))
				return 0;
		}
		if (Rupdateflag) {
			if (!getline(obsfileR, buffR))
				return 0;
		}
	}
	/*if (numiterator >= debug_stop - 1) {
		int aaaa = 1;
	}*/
	rcv.updateobs(obs);
	rcv.restoreslips();
	rcv.updateSatinfo();
	
	return 1;
}
void Fusion_set_RTKx_posvel(receiver & rcv, const Eigen::Vector3d & Ext_pos, const Eigen::Vector3d & Ext_vel) {
	for (int i = 0; i < 3; i++)
		rcv.RTK_x[i] = Ext_pos[i];
	if (rcv.dynamics_on) {
		for (int i = 0; i < 3; i++)
			rcv.RTK_x[i+3] = Ext_vel[i];
	}
	return;
}
void Fusion_set_RTKxa_posvel(receiver & rcv, const Eigen::Vector3d & Ext_pos, const Eigen::Vector3d & Ext_vel) {
	for (int i = 0; i < 3; i++)
		rcv.RTK_fixXpva[i] = Ext_pos[i];
	if (rcv.dynamics_on) {
		for (int i = 0; i < 3; i++)
			rcv.RTK_fixXpva[i + 3] = Ext_vel[i];
	}
	return;
}
void Fusion_set_RTKxP_bias(receiver & rcv, const Eigen::VectorXd &x, const Eigen::MatrixXd &P, const vector<int> &TC_ix) {//TC结果直接补偿GNSS的单差状态中
	int start_Sbias_index = rcv.dynamics_on ? 9 : 3;
	int rcv_nx = rcv.get_nx();
	for (int i = 15; i < TC_ix.size(); i++) {
		int ii = TC_ix[i] - 15 + start_Sbias_index;
		rcv.RTK_x[ii] += x[i];
		for (int j = 15; j < TC_ix.size(); j++) {
			int jj = TC_ix[j] - 15 + start_Sbias_index;
			rcv.RTK_P[ii+jj* rcv_nx] = P(i, j);
		}
	}
	return;
}
//void receiver::plot_ENU_traj() {
//	static int flag = 0;
//	static vector<double > posbasellhxyz(3, 0);
//	
//	if (!flag) {
//		plt::clf();
//		plt::title("The trajectory.");
//		ecef2pos(basepos, posbasellhxyz);// base 的位置转llh	
//		flag++;
//	}
//	vector<double > posxyz(3, 0);
//	vector<double > e_rb(3, 0);
//	posxyz = sol.getpos();
//	int stat = sol.get_state();
//	if (plot_mode == 2 && stat != 1 && stat != 2)
//		return;
//	for (int i = 0; i < 3; i++)
//		e_rb[i] = posxyz[i] - basepos[i];
//	vector<double > enu(3, 0);
//	ecef2enu(posbasellhxyz, e_rb, enu);
//	vector<double > e_data,n_data;
//	e_data.push_back(enu[0]);
//	n_data.push_back(enu[1]);
//	plt::scatter(e_data, n_data);
//	plt::pause(0.0001);
//	//plt::draw();
//}