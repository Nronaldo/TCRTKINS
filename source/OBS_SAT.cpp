#include "OBS_SAT.h"

void solution::PQxyz2PQllh() {
	ecef2pos(pos,posllh);
	Eigen::Matrix3d P;
	Eigen::Matrix3d Q = Eigen::Matrix3d::Zero();
	P << qr[0], qr[3], qr[5],
		 qr[3], qr[1], qr[4],
		 qr[5], qr[4], qr[2];
	covenu(posllh, P, Q);
	posllh[0] *= 180 / PI;
	posllh[1] *= 180 / PI;
	for (int i = 0; i < 3; i++)
		qrenu[i] = Q(i, i);
	qrenu[3] = Q(0, 1);
	qrenu[4] = Q(1, 2);
	qrenu[5] = Q(0, 2);
	return;
}
OBS_SAT::OBS_SAT(int satindex, int satindexR, const OBS_RINEX &obsrinex, char mode) {
	base_flag = 0;
	rover_flag = 0;
	BR_unused_flag = 0;
	BR_fde_mask = 0;
	int GNSS_Fnum_B;
	int GNSS_Fnum_R;
	//double Cspeed = 299792458.0;
	vector<double> GNSS_fre_B;
	vector<double> GNSS_fre_R;
	vector<int > GNSS_satID_B;
	vector<int > GNSS_satID_R;
	vector<int > GNSS_flags_B;
	vector<int > GNSS_flags_R;
	vector<double > GNSS_measurements_B;
	vector<double > GNSS_measurements_R;
	
	BRtime = obsrinex.GPStime_R;
	if (mode == 'G') {
		constellation = 'G';
		GNSS_Fnum_B = obsrinex.GPS_Fnum_B;
		GNSS_Fnum_R = obsrinex.GPS_Fnum_R;
		GNSS_fre_B = obsrinex.GPS_fre_B;
		GNSS_fre_R = obsrinex.GPS_fre_R;
		GNSS_satID_B = obsrinex.GPS_satID_B;
		GNSS_satID_R = obsrinex.GPS_satID_R;
		GNSS_flags_B = obsrinex.GPS_flags_B;
		GNSS_flags_R = obsrinex.GPS_flags_R;
		GNSS_measurements_B = obsrinex.GPS_measurements_B;
		GNSS_measurements_R = obsrinex.GPS_measurements_R;
	}
	else {
		constellation = 'C';
		GNSS_Fnum_B = obsrinex.BDS_Fnum_B;
		GNSS_Fnum_R = obsrinex.BDS_Fnum_R;
		GNSS_fre_B = obsrinex.BDS_fre_B;
		GNSS_fre_R = obsrinex.BDS_fre_R;
		GNSS_satID_B = obsrinex.BDS_satID_B;
		GNSS_satID_R = obsrinex.BDS_satID_R;
		GNSS_flags_B = obsrinex.BDS_flags_B;
		GNSS_flags_R = obsrinex.BDS_flags_R;
		GNSS_measurements_B = obsrinex.BDS_measurements_B;
		GNSS_measurements_R = obsrinex.BDS_measurements_R;
	}
	if (satindex >= 0) {
		base = Station_OBS(GNSS_Fnum_B, satindex, GNSS_measurements_B, GNSS_flags_B, GNSS_fre_B, obsrinex.GPStime_B);
		base_flag = 1;//说明当前星对base可见
	}
	if (satindexR >= 0) {
		rover = Station_OBS(GNSS_Fnum_R, satindexR, GNSS_measurements_R, GNSS_flags_R, GNSS_fre_R, obsrinex.GPStime_R);
		rover_flag = 1;
	}
	
}
void OBS_SAT::displayOBS_SAT() {
	cout << endl << "Base:" << endl;
	base.displayOBS_station();
	
	cout << endl << "Rover:" << endl;
	rover.displayOBS_station();
	
	cout << endl;
}

void OBS_SAT::displayOBS_SAT_Pos(int mode) {
	if (mode == 1) {
		cout << "****Base Sat info:" << endl;
		base.displaySAT_station();
		cout << "****Rover Sat info:" << endl;
		rover.displaySAT_station();
		cout << endl << "****SAT SVH FLAG：" << endl;
		cout << SAT_SVH << " ";
		cout << endl;
	}
	else if (mode == 2) {
		int unused_flag = rover.get_unusedflag();
		vector<double > r_azel_los = rover.get_razellos();
		cout << endl << "Unused flag：" << endl;
		cout << unused_flag;
		cout << endl << "distance between user and sat：" << endl;
		cout << setprecision(17) << r_azel_los[0] << endl;
		cout << endl << "azimuth，elevation：" << endl;
		for (int i = 1; i < 3; i++)
			cout << setprecision(17) << r_azel_los[i] << " ";
		cout << endl << "unit Los：" << endl;
		for (int i = 3; i < 6; i++)
			cout << setprecision(17) << r_azel_los[i] << " ";
		cout << endl;
	}


	return;
}
void OBS_SAT::detslp_ll(char * slip, char * half,  double tt) {
	//只考虑单频
	*half = (base.Sget_LLI(0) & 2) ? 0 : (rover.Sget_LLI(0) & 2) ? 0 : 1;// 这个应该是指示是否半周有效
	if (base.Sget_L(0) == 0.0 || rover.Sget_L(0) == 0.0)
		return;
	int rcvpreviousLLI, basepreviousLLI,roverslip,baseslip;
	rcvpreviousLLI = getbitu(slip, 0, 2);
	basepreviousLLI = getbitu(slip, 2, 2);
	if (tt >= 0.0) {
		baseslip = base.Sget_LLI(0);
		roverslip = rover.Sget_LLI(0);
	}
	else {
		baseslip = basepreviousLLI;
		roverslip = rcvpreviousLLI;
	}
	if (((basepreviousLLI & 2) && !(base.Sget_LLI(0) & 2)) || (!(basepreviousLLI & 2) && (base.Sget_LLI(0) & 2))) {
		baseslip |= 1;
	}// slip 低4位 是一个标志，只统计是否发生周跳，不管是整周还是半周。不关心。
	if (((rcvpreviousLLI & 2) && !(rover.Sget_LLI(0) & 2)) || (!(rcvpreviousLLI & 2) && (rover.Sget_LLI(0) & 2))) {
		roverslip |= 1;
	}// slip 低4位 是一个标志，只统计是否发生周跳，不管是整周还是半周。不关
	setbitu(slip, 0, 2, rover.Sget_LLI(0));// 更新rtk->ssat 状态
	setbitu(slip, 2, 2, base.Sget_LLI(0));// 更新rtk->ssat 状态
	*slip|= (char)baseslip;//rtk->ssat[sat-1].slip[f] 包含了rover 和 base的，所以这里要或
	*slip|= (char)roverslip;
	
	return;
}
void OBS_SAT::orBRLLI(int b, int r) {
	//目前只考虑单频
	if (base_flag) {
		int B_originalLLI = base.Sget_LLI(0);
		int B_currentLLI = B_originalLLI | b;
		base.Sset_LLI(0, B_currentLLI);
	}
	if (rover_flag) {
		int R_originalLLI = rover.Sget_LLI(0);
		int R_currentLLI = R_originalLLI | r;
		rover.Sset_LLI(0, R_currentLLI);
	}	
	return;
}

void OBS_SAT::setpvbdst(const vector<double > &pvb, const vector<double > & bdb, double varb,
					  const vector<double > &pvr, const vector<double > & bdr, double varr,
					  int svh, double tgd, int bflag, int rflag) {
	if (bflag) {
		base.Ssetpvbd(pvb, bdb, varb);
	}
	if (rflag) {
		rover.Ssetpvbd(pvr, bdr, varr);
	}
	SAT_SVH = svh;
	TGD = tgd;
	return;
}
double OBS_SAT::calculate_TDCPres(const vector<double > &last_sol_rr) {
	if (SAT_SVH == 1||SAT_SVH<0) 
		return 0;
	return rover.Scalculate_TDCPres(last_sol_rr);
}
void OBS_SAT::Ls_update_razellos(const vector<double > &sysclk_bias, const vector<double> &receiver_pos_ecef,
	const vector<double> &receiver_pos_llh, int &count, const vector<double > &ION_CORR) {
	rover.update_razellos(sysclk_bias, receiver_pos_ecef, receiver_pos_llh,count,ION_CORR,SAT_SVH,constellation,TGD);
	return;
}
void OBS_SAT::Ls_updateDres(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x, int &count) {
	rover.updateDres(receiver_pos_ecef,E,Ls_x,count);
	return;
}

int OBS_SAT::BRzdres(char mode, const vector<double > & posecef, const vector<double > & posllh) {
	if (numiterator >= debug_stop && mode == 'R') {
		int aaa = 1;
	
	}
		
	if (SAT_SVH < 0 || SAT_SVH == 1) {
		BR_unused_flag = 1;
		return 0;
	}
	int state = 0;
	if (mode == 'B') {
		state =  base.zdres(posecef, posllh);
		
	}
	else if (mode == 'R') {
		state =  rover.zdres(posecef, posllh);
	}
	else {
		cout << "no such zdres mode" << endl;
		return 0;
	}
	if(!state ) //base或者rover有一个站计算零差不成功，都不能参与最后的rtk(P+L)解算。
		BR_unused_flag = 1;
	return state;
}
void OBS_SAT::calsdresvar(double bl) {
	double a, b, c, d;
	double el = rover.get_razellos()[2];
	double sinel = sin(el);
	int timeb, timer;
	double secb, secr;
	rover.gettime(&timer, &secr);
	base.gettime(&timeb, &secb);
	double dt = timer-timeb+secr-secb;
	c = TC_conf.Error_std[3] * bl / 1E4;
	d = Cspeed * TC_conf.sclkstab*dt;
	double factL = 1.0,factP = TC_conf.eratio[0];
	a = factL * TC_conf.Error_std[1];
	b = factL * TC_conf.Error_std[2];
	L_svar = 2.0*(a*a + b * b / sinel / sinel + c * c) + d * d;

	a = factP * TC_conf.Error_std[1];
	b = factP * TC_conf.Error_std[2];
	P_svar = 2.0*(a*a + b * b / sinel / sinel + c * c) + d * d;
	return;
}

Station_OBS::Station_OBS(const int GNSS_Fnum, const int satindex, const vector<double > & GNSS_measurements, const vector<int > & GNSS_flags, const vector<double > & GNSS_fre, gtime_t timerinex) {
	unused_flag = 0;
	fde_mask = 0;
	r_azel_los = vector<double>(6, 0);
	time = timerinex;
	SAT_POS_VEL = vector<double>(6, 0);
	SATCLK_BIAS_DRIFT = vector<double >(2, 0);
	cp_t_tu_res = 0;
	for (int i = 0; i < GNSS_Fnum; i++) {
		P.push_back(GNSS_measurements[(satindex*GNSS_Fnum + i) * 4 + 0]);
		L.push_back(GNSS_measurements[(satindex*GNSS_Fnum + i) * 4 + 1]);
		D.push_back(GNSS_measurements[(satindex*GNSS_Fnum + i) * 4 + 2]);
		S.push_back(GNSS_measurements[(satindex*GNSS_Fnum + i) * 4 + 3]);
		qualP.push_back(GNSS_flags[(satindex*GNSS_Fnum + i) * 3 + 0]);
		qualL.push_back(GNSS_flags[(satindex*GNSS_Fnum + i) * 3 + 1]);
		LLI.push_back(GNSS_flags[(satindex*GNSS_Fnum + i) * 3 + 2]);
		lambda.push_back(Cspeed / GNSS_fre[i] / 1e6);
	}
	P_zres = 0;//零差
	L_zres = 0;
	D_zres = 0;
	P_zvar = 0;
}
int Station_OBS::zdres(const vector<double > & posecef, const vector<double > & posllh) {
	if (numiterator >= debug_stop) {
		int aaa = 1;

	}
	if (posecef.empty() || norm_vector(posecef, 0, 3) <= 0.0)
		return 0;

	geodist(posecef);
	satazel(posllh);
	/*if (unused_flag == 1)
		return 0;*/
	unused_flag = 0;
	if (r_azel_los[0] < 0 || r_azel_los[2] < TC_conf.elmin_mask)//
		unused_flag = 1;
	if (unused_flag || fde_mask)
		return 0;
	double r = r_azel_los[0];
	r -= Cspeed * SATCLK_BIAS_DRIFT[0];
	vector<double  > zazel = { 0.0,90.0*D2R };
	double zhd = tropmodel(posllh, zazel, 0.0);
	vector<double  > azel = { r_azel_los[1],r_azel_los[2] };
	double zhd_mapf = tropmapf(time, posllh, azel);
	r += zhd * zhd_mapf;
	//此处没有做antmodel
	//只考虑了单频的情况，后续再扩展iono-free linear combination
	if (numiterator >= debug_stop) {
		int aaa = 1;

	}
	if (lambda[0] == 0.0 || L[0] == 0.0 || P[0] == 0.0) {
		unused_flag = 1;
		return 0;
	}
	L_zres = L[0] * lambda[0] - r;
	P_zres = P[0] - r;
	return 1;
}


void Station_OBS::displayOBS_station() {
	for (int i = 0; i < P.size(); i++)
		cout << setprecision(12) << P[i] << " " << qualP[i] << " "
		<< setprecision(10) << L[i] << " " << LLI[i] << " " << qualL[i]
		<< setprecision(8) << " " << D[i] << " " << setprecision(4) << S[i] << " ";
}
void Station_OBS::displaySAT_station() {
	cout << endl << "Position：" << endl;
	for (int i = 0; i < 3; i++)
		cout << setprecision(17) << SAT_POS_VEL[i] << " ";
	cout << endl << "Velocity：" << endl;
	for (int i = 3; i < 6; i++)
		cout << setprecision(17) << SAT_POS_VEL[i] << " ";
	cout << endl << "Clock Bias：" << endl;
	cout << setprecision(17) << SATCLK_BIAS_DRIFT[0] << " ";
	cout << endl << "Clock drift：" << endl;
	cout << setprecision(17) << SATCLK_BIAS_DRIFT[1] << " ";
	cout << endl << "VAR：" << endl;
	cout << setprecision(17) << SATPOS_VAR << " ";

}

//即使是L1信号，也要dcb调整一下
// 单频接收机，不播发频率差异，用群波延时校正，实际上就是群波延时校正卫星钟差？
double  Station_OBS::cal_P(double *P, string mode, const vector<int> &fre,char constllation,double TGD) {
	if (mode == "single") {
		*P = Sget_P(fre[0]);
		if (fabs(*P) < 1e-14)//有可能当前频点不存在观测量
			return 0;
		double l1, lc;
		l1 = Sget_lambda(0);
		lc = Sget_lambda(fre[0]);
		double gamma = (lc / l1)*(lc / l1);
		if (constllation == 'G') // RTKLIB没有校正北斗TGD 应该改过来，且取负值
			*P = *P - TGD * Cspeed;
		else if (constllation == 'C')
			*P = *P + TGD * Cspeed;
		return gamma;
	}
	return 0.0;
}


void Station_OBS::update_razellos(const vector<double > &sysclk_bias, const vector<double> &receiver_pos_ecef,
	const vector<double> &receiver_pos_llh, int &count, const vector<double > &ION_CORR, double sat_svh, char constellation,double TGD) {
	/*sysclk_bias 是估计出来的接收机钟差*/
	geodist(receiver_pos_ecef);
	satazel(receiver_pos_llh);
	unused_flag = 0;
	//这里有一个反复横跳的过程。。
	// 有一个星很靠近cutoff，LS迭代中前次低于cutoff被舍弃了，后次高于curoff又回来了，但是前次被清零后一直没有回复。
	// 所以出现和RTKIIB 不一致的情况。理论上应该恢复
	if (r_azel_los[0] <= 0 || r_azel_los[2] < TC_conf.elmin_mask || sat_svh == 1 || sat_svh < 0)//
		unused_flag = 1;
	if (unused_flag || fde_mask)
		return;
	count++;
	double P;
	double Fre_rate = cal_P(&P, "single", vector<int>(1, 0),constellation,TGD);
	if (fabs(P) < 1e-14) {
		unused_flag = 1;
		return;
	}
	double vmeas = ERR_CBIAS * ERR_CBIAS;//沿用RTKLIB的方法
	double dion = ion_BRDC(receiver_pos_llh, ION_CORR);
	double vion = (dion*dion)*(ERR_BRDCI*ERR_BRDCI);
	if (constellation == 'G')
		dion *= Fre_rate;//单频时ION模型频率调整
	else {
		double lr;
		lr = Sget_lambda(0);
		dion = dion * lr*lr / lam_carr[0] / lam_carr[0];
	}
	vector<double > azel(2, 0);
	azel[0] = r_azel_los[1];
	azel[1] = r_azel_los[2];
	double dtrop = tropmodel(receiver_pos_llh,azel,REL_HUMI);
	double vtrop = pow(ERR_SAAS / (sin(r_azel_los[2]) + 0.1), 2);
	P_zres = P - (r_azel_los[0] + sysclk_bias[0] - Cspeed * SATCLK_BIAS_DRIFT[0] + dion + dtrop);//SATCLK_BIAS_DRIFT 是卫星钟差，P在原始观测量基础上补偿了TGD
	P_zvar = SATPOS_VAR + vmeas + vion + vtrop + varerr(r_azel_los[2]);
	if (constellation == 'C')
		P_zres -= sysclk_bias[1];
	return;
}

void Station_OBS::Ssetpvbd(const vector<double > &pv, const vector<double > & bd, double var) {
	SAT_POS_VEL[0] = pv[0];
	SAT_POS_VEL[1] = pv[1];
	SAT_POS_VEL[2] = pv[2];
	SATCLK_BIAS_DRIFT[0] = bd[0];
	SAT_POS_VEL[3] = pv[3];
	SAT_POS_VEL[4] = pv[4];
	SAT_POS_VEL[5] = pv[5];
	SATCLK_BIAS_DRIFT[1] = bd[1];
	SATPOS_VAR = var;
	return;
}

void Station_OBS::geodist(const vector<double> &receiver_pos_ecef) {
	double rr = norm_vector(SAT_POS_VEL, 0, 3);
	if (rr < RE_WGS84) {
		r_azel_los[0] = -1.0;
	}
	else {
		for (int i = 0; i < 3; i++)
			r_azel_los[3 + i] = SAT_POS_VEL[i] - receiver_pos_ecef[i];
		double r = norm_vector(r_azel_los, 3, 3);
		divide_vector(r_azel_los, 3, 3, r);//第一个参数是引用，会修改los。
		r_azel_los[0] = r + 7.2921151467E-5*(SAT_POS_VEL[0] * receiver_pos_ecef[1] - SAT_POS_VEL[1] * receiver_pos_ecef[0]) / Cspeed;
	}
	return;
}
void Station_OBS::satazel(const vector<double> &receiver_pos_llh) {
	double az = 0.0, el = PI / 2.0;
	vector<double > e(r_azel_los.end() - 3, r_azel_los.end());
	vector<double > enu(3, 0);
	if (receiver_pos_llh[2] > -RE_WGS84) {
		ecef2enu(receiver_pos_llh, e, enu);
		az = dot(enu, 0, enu, 0, 2) < 1e-12 ? 0.0 : atan2(enu[0], enu[1]);
		if (az < 0.0)
			az += 2 * PI;
		el = asin(enu[2]);
	}
	r_azel_los[1] = az;
	r_azel_los[2] = el;
	return;
}
double Station_OBS::ion_BRDC(const vector<double> &pos, const vector<double> & ION_CORR) {
	const double ion_default[] = { /* 2004/1/1 */
		0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
		0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
	};
	vector<double> ion(8, 0);
	double tt, f, psi, phi, lam, amp, per, x;
	int week;
	if (pos[2] < -1E3 || r_azel_los[2] <= 0) return 0.0;
	if (ION_CORR.empty() || norm_vector(ION_CORR, 0, 8) <= 0.0) {
		for (int i = 0; i < 8; i++)
			ion[i] = ion_default[i];
	}
	else {
		ion = ION_CORR;
	}
	psi = 0.0137 / (r_azel_los[2] / PI + 0.11) - 0.022;
	/* subionospheric latitude/longitude (semi-circle) */
	phi = pos[0] / PI + psi * cos(r_azel_los[1]);
	if (phi > 0.416) phi = 0.416;
	else if (phi < -0.416) phi = -0.416;
	lam = pos[1] / PI + psi * sin(r_azel_los[1]) / cos(phi*PI);
	/* geomagnetic latitude (semi-circle) */
	phi += 0.064*cos((lam - 1.617)*PI);

	/* local time (s) */
	tt = 43200.0*lam + time2gpst(time, &week);
	tt -= floor(tt / 86400.0)*86400.0; /* 0<=tt<86400 */

	/* slant factor */
	f = 1.0 + 16.0*pow(0.53 - r_azel_los[2] / PI, 3.0);
	/* ionospheric delay */
	amp = ion[0] + phi * (ion[1] + phi * (ion[2] + phi * ion[3]));
	per = ion[4] + phi * (ion[5] + phi * (ion[6] + phi * ion[7]));
	amp = amp < 0.0 ? 0.0 : amp;
	per = per < 72000.0 ? 72000.0 : per;
	x = 2.0*PI*(tt - 50400.0) / per;

	return Cspeed * f*(fabs(x) < 1.57 ? 5E-9 + amp*(1.0 + x * x*(-0.5 + x * x / 24.0)) : 5E-9);
}
double Station_OBS::tropmodel(const vector<double> &pos,const vector<double > &azel, double humi) {
	const double temp0 = 15.0; /* temparature at sea level */
	double hgt, pres, temp, e, z, trph, trpw;

	if (pos[2] < -100.0 || 1E4 < pos[2] || azel[1] <= 0) return 0.0;

	/* standard atmosphere */
	hgt = pos[2] < 0.0 ? 0.0 : pos[2];

	pres = 1013.25*pow(1.0 - 2.2557E-5*hgt, 5.2568);
	temp = temp0 - 6.5E-3*hgt + 273.16;
	e = 6.108*humi*exp((17.15*temp - 4684.0) / (temp - 38.45));

	/* saastamoninen model */
	z = PI / 2.0 - azel[1];
	trph = 0.0022768*pres / (1.0 - 0.00266*cos(2.0*pos[0]) - 0.00028*hgt / 1E3) / cos(z);
	trpw = 0.002277*(1255.0 / temp + 0.05)*e / cos(z);
	return trph + trpw;
}
void Station_OBS::updateDres(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x, 
	int &count) {
		//updateDres_Doppler(receiver_pos_ecef, E, Ls_x, count);
	if (unused_flag || fde_mask)
		return;
	double l, doppler;
	doppler = D[0];
	l = lambda[0];
	if (doppler == 0.0 || l == 0 || norm_vector(SAT_POS_VEL, 3, 3) <= 0.0) {
		unused_flag = 1;
		return;
	}
	double cosel = cos(r_azel_los[2]);
	Eigen::Vector3d a, e;
	a << sin(r_azel_los[1])*cosel, cos(r_azel_los[1])*cosel, sin(r_azel_los[2]);
	e = E.transpose() * a;// enu->xyz
	vector<double > vs(3, 0);
	for (int j = 0; j < 3; j++) {
		r_azel_los[3 + j] = e(j);
		vs[j] = SAT_POS_VEL[j + 3] - Ls_x(j);
	}
	double x1 = dot(vs, 0, r_azel_los, 3, 3);
	double x2 = OMGE / Cspeed * (SAT_POS_VEL[4] * receiver_pos_ecef[0] + SAT_POS_VEL[1] * Ls_x(0) -
		SAT_POS_VEL[3] * receiver_pos_ecef[1] - SAT_POS_VEL[0] * Ls_x(1));
	double rate = dot(vs, 0, r_azel_los, 3, 3) + OMGE / Cspeed * (SAT_POS_VEL[4] * receiver_pos_ecef[0] + SAT_POS_VEL[1] * Ls_x(0) -
		SAT_POS_VEL[3] * receiver_pos_ecef[1] - SAT_POS_VEL[0] * Ls_x(1));
	D_zres = -l * doppler - (rate + Ls_x(3) - Cspeed * SATCLK_BIAS_DRIFT[1]);
	count++;
	return;
}

//void Station_OBS::updateDres_Doppler(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x,
//	int &count){
//	
//}
void Station_OBS::Sset_LLI(int Fi, int value) {
	if (Fi >= LLI.size()) {
		cout << "LLI index out of range" << endl;
		return;
	}
	LLI[Fi] = value;
	return;
}

void OBS_SAT::calP_sres() {
	double Pb = base.Sget_P(0);
	double Pr = rover.Sget_P(0);
	if (fabs(Pb) > 1e-15 && fabs(Pr) > 1e-15)
		P_sres = Pr - Pb;
	else
		P_sres = 0;
}
void OBS_SAT::calL_sres() {
	double Lb = base.Sget_L(0);
	double Lr = rover.Sget_L(0);
	if (fabs(Lb) > 1e-15 && fabs(Lr) > 1e-15)
		L_sres = Lr - Lb;
	else
		L_sres = 0;
}
void OBS_SAT::calD_sres() {
	double Db = base.Sget_D(0);
	double Dr = rover.Sget_D(0);
	if (fabs(Db) > 1e-15 && fabs(Dr) > 1e-15)
		D_sres = Dr - Db;
	else
		D_sres = 0;
}
void Station_OBS::Supdate_cp_t_tu_res(const vector<double> &posxyz,int SVH) {// 这里传入的是上个历元的最终解算位置
	/*/*GINAV保存了old_obs 和 cur_obs,但是old_obs中有的星，上一历元可能没参与解算，比如高度角不够
			然而此处都直接参与解算了，没有排除，或许这是一个小问题。这里我们先保持一致*/
	if (!L[0] || !lambda[0]||SVH<0||SVH==1)
		cp_t_tu_res = 0;
	else {
		geodist(posxyz);
		cp_t_tu_res = L[0] * lambda[0] - r_azel_los[0] + Cspeed * SATCLK_BIAS_DRIFT[0];
	}
	return;
}
double  Station_OBS::Scalculate_TDCPres(const vector<double > &last_sol_rr) {
	if (!cp_t_tu_res||unused_flag||fde_mask)//上一时刻没本星，或上一时刻本星不可用，或本时刻本星不可用直接返回
		return 0;
	if (!L[0] || !lambda[0])//本时刻没有载波相位
		return 0;
	geodist(last_sol_rr);
	vector<double > lastposllh(3, 0);
	ecef2pos(last_sol_rr, lastposllh);
	satazel(lastposllh);
	double cur_cp_t_tu_res = L[0] * lambda[0] - r_azel_los[0] + Cspeed * SATCLK_BIAS_DRIFT[0];
	double result = cur_cp_t_tu_res - cp_t_tu_res;
	return result;

	

}