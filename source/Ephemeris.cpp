#include "Ephemeris.h"

int uraindex(double value) {
	int  i = 0;
	for (i = 0; i < 15; i++) 
		if (ura_value[i] >= value) 
			break;
	return i;
}

Sat_EPH::Sat_EPH(int SatID,vector<double>::iterator Satnavmessage) {
	sat = SatID;
	f0 = *Satnavmessage;
	f1 = *(Satnavmessage+1);
	f2 = *(Satnavmessage+2);
	A = *(Satnavmessage+10)* (*(Satnavmessage+10));
	e = *(Satnavmessage+8);
	i0 = *(Satnavmessage+15);
	OMG0 = *(Satnavmessage+13);
	omg = *(Satnavmessage+17);
	M0 = *(Satnavmessage+6);
	deln = *(Satnavmessage+5);
	OMGd = *(Satnavmessage+18);
	idot = *(Satnavmessage+19);
	crc = *(Satnavmessage+16);
	crs = *(Satnavmessage+4);
	cuc = *(Satnavmessage+7);
	cus = *(Satnavmessage+9);
	cic = *(Satnavmessage+12);
	cis = *(Satnavmessage+14);
	svh = (int)(*(Satnavmessage + 24));
	sva = uraindex(*(Satnavmessage + 23));
	week = (int)(*(Satnavmessage + 21));
	toes = *(Satnavmessage + 11);
}
int Sat_EPH::getCRS(double &CRS) {
	CRS = crs;
	return sat;
}
double Sat_EPH::cal_POSVAR() {
	int ura = sva;
	return ura < 0 || 14 < ura ? 6144.0*6144.0 : ura_value[ura]* ura_value[ura];
}

void Sat_EPH::calculate_Satpos_inorbit(gtime_t time) {
	//选择合适卫星的操作在obs里面做，此处只做初始卫星钟差，轨道平面内卫星位置计算
	// 所以如果能进入本函数，一定是选择了距离obs时间最近的卫星了
	// 传入的time 一定是 补偿了伪距的发射时间
	// BDS包含 GEO卫星，因此将GPS BDS计算ECEF坐标放置到各自函数内去做
	double t = timediff(time,toc);
	for (int i = 0; i < 2; i++) { //RTKLIB的做法，这里沿用
		t -= f0 + f1*t + f2*t*t;
	}
	double dt =  f0 + f1*t + f2*t*t;
	time  = timeadd(time,-dt);
	double M, E, Ek, cosE, u, r, i, sin2u, cos2u, x, y;
	int n;
	tk = timediff(time,toe);
	M = M0 + (sqrt(mu / (A*A*A)) + deln)*tk;
	for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > 1E-13 &&n < 30; n++) {
		Ek = E;
		E -= (E - e*sin(E) - M) / (1.0 - e*cos(E));
	}
	sinE = sin(E); 
	cosE = cos(E);
	u = atan2(sqrt(1.0 - e*e)*sinE, cosE - e) + omg;
	r = A*(1.0 - e*cosE);
	i = i0 + idot*tk;
	sin2u = sin(2.0*u); 
	cos2u = cos(2.0*u);
	u += cus*sin2u + cuc*cos2u;
	r += crs*sin2u + crc*cos2u;
	i += cis*sin2u + cic*cos2u;
	x = r * cos(u); 
	y = r * sin(u); 
	//cosi = cos(i);
	x_inorbit = x;
	y_inorbit = y;
	i_inorbit = i;
	return;
}


GPS_Sat_EPH::GPS_Sat_EPH(int SatID, vector<double>::iterator Satnavmessage, int SatToc):Sat_EPH(SatID, Satnavmessage) {
	mu = 3.9860050E14;
	omge = 7.2921151467E-5;
	toc.time = SatToc;
	toc.sec = 0.0;
	iode = (int)(*(Satnavmessage+3));
	iodc = (int)(*(Satnavmessage+26));
	toe = adjweek(gpst2time(week, *(Satnavmessage+11)),toc);
	ttr = adjweek(gpst2time(week, *(Satnavmessage+27)),toc);
	code = (int)(*(Satnavmessage+20));
	flag = (int)(*(Satnavmessage+22));
	tgd = *(Satnavmessage+25);
	fit = *(Satnavmessage+28);
}
void GPS_Sat_EPH:: Cali_GPS_Satinfo(gtime_t obstime, vector<double > &SAT_Pos, double  &CLK_BIAS) {
	double O,cosO,sinO,cosi,sini;
	calculate_Satpos_inorbit(obstime);
	O = OMG0 + (OMGd - omge)*tk - omge * toes;
	sinO = sin(O); 
	cosO = cos(O);
	cosi = cos(i_inorbit);
	sini = sin(i_inorbit);
	SAT_Pos[0] = x_inorbit * cosO - y_inorbit * cosi*sinO;
	SAT_Pos[1] = x_inorbit * sinO + y_inorbit * cosi*cosO;
	SAT_Pos[2] = y_inorbit * sini;
	tk = tk + timediff(toe,toc);
	double CLIGHT = 299792458.0;
	double x = sqrt(mu*A);
	CLK_BIAS = f0 + f1*tk + f2*tk*tk- 2.0*x*e*sinE / CLIGHT/ CLIGHT;
	return;
}
BDS_Sat_EPH::BDS_Sat_EPH(int SatID, vector<double>::iterator Satnavmessage, int SatToc) :Sat_EPH(SatID, Satnavmessage) {
	mu = 3.986004418E14;
	omge= 7.292115E-5;
	toc = bdt2gpst(gtime_t(SatToc,0.0));
	iode = (int)(*(Satnavmessage + 3));
	iodc = (int)(*(Satnavmessage + 28));
	toe = bdt2gpst(bdt2time(week, *(Satnavmessage + 11))); /* bdt -> gpst */
	ttr = bdt2gpst(bdt2time(week, *(Satnavmessage + 27))); /* bdt -> gpst */
	toe = adjweek(toe, toc);
	ttr = adjweek(ttr, toc);

	tgd[0] = *(Satnavmessage+25);
	tgd[1] = *(Satnavmessage+26);
}
void BDS_Sat_EPH::Cali_BDS_Satinfo(int ID,gtime_t obstime, vector<double > &SAT_Pos, double  &CLK_BIAS) {
	double O, cosO, sinO, cosi, sini;
	calculate_Satpos_inorbit(obstime);
	cosi = cos(i_inorbit);
	sini = sin(i_inorbit);
	if (ID <= 5) {
		O = OMG0 + OMGd*tk - omge *toes;
		sinO = sin(O); 
		cosO = cos(O);
		double xg = x_inorbit * cosO - y_inorbit * cosi*sinO;
		double yg = x_inorbit * sinO + y_inorbit * cosi*cosO;
		double zg = y_inorbit * sin(i_inorbit);
		double sino = sin(omge*tk);
		double coso = cos(omge*tk);
		double COS_5 = 0.9961946980917456;
		double SIN_5 = -0.0871557427476582;
		SAT_Pos[0] = xg * coso + yg * sino*COS_5 + zg * sino*SIN_5;
		SAT_Pos[1] = -xg * sino + yg * coso*COS_5 + zg * coso*SIN_5;
		SAT_Pos[2] = -yg * SIN_5 + zg * COS_5;
		
	}
	else { // 非GEO卫星 和 GPS这段代码一致
		O = OMG0 + (OMGd - omge)*tk - omge * toes;
		sinO = sin(O);
		cosO = cos(O);
		cosi = cos(i_inorbit);
		sini = sin(i_inorbit);
		SAT_Pos[0] = x_inorbit * cosO - y_inorbit * cosi*sinO;
		SAT_Pos[1] = x_inorbit * sinO + y_inorbit * cosi*cosO;
		SAT_Pos[2] = y_inorbit * sini;
	}
	tk = tk + timediff(toe, toc);
	double CLIGHT = 299792458.0;
	double x = sqrt(mu*A);
	CLK_BIAS = f0 + f1 * tk + f2 * tk*tk - 2.0*x*e*sinE / CLIGHT / CLIGHT;
	return;
}

All_Sats_EPH::All_Sats_EPH(NAV_RINEX &nav) {
    LEAP_SECONDS = nav.LEAP_SECONDS;
	for (int i = 0; i < nav.CORR_TYPE.size(); i++) {
		if (nav.CORR_TYPE[i] == "GPSA"|| nav.CORR_TYPE[i] == "GPSB") {
			for (int j = 0; j < 4; j++) 
				GPS_ION_CORR.push_back(nav.ION_CORR[4*i+j]);
		}
		else if(nav.CORR_TYPE[i] == "BDSA" || nav.CORR_TYPE[i] == "BDSB") {
			for (int j = 0; j < 4; j++)
				BDS_ION_CORR.push_back(nav.ION_CORR[4*i+j]);
			
		}
		else if (nav.CORR_TYPE[i] == "GPUT") {
			for (int j = 0; j < 4; j++)
				GPS_TIME_CORR.push_back(nav.TIME_CORR[j]);
		}
		else {
			for (int j = 0; j < 4; j++)
				BDS_TIME_CORR.push_back(nav.TIME_CORR[j]);
		}
	}
}
void All_Sats_EPH:: satephsolve(NAV_RINEX & nav) {
	set<string> Gset,Bset;
	for (int i = 0; i < nav.GPS_ID.size(); i++) {
		GPS_Sat_EPH GPSsateph(nav.GPS_ID[i], nav.GPS_NAV_DATA.begin() + i * 29, nav.GPS_TOC[i]);
		int iode;
		double crs;
		int ID;
		iode = GPSsateph.get_iode();
		ID = GPSsateph.getCRS(crs);// 有时候rinex给的电文很多参数都是0，这样的星我们舍掉
		string Gs = to_string(ID) + to_string(iode);
		if ((Gset.empty() || !Gset.count(Gs))&&crs) {
			Gset.insert(Gs);
			Gmmp.insert(gepair(ID, GPSsateph));
		}

	}
	for (int i = 0; i < nav.BDS_ID.size(); i++) { //B1I 和 B3I这些信号和新体制的B2a信号不同，前者采用AODE，并不像后者和GPS的IODE一样6小时不重复。不能依赖这个去重
		BDS_Sat_EPH BDSsateph(nav.BDS_ID[i], nav.BDS_NAV_DATA.begin() + i * 29, nav.BDS_TOC[i]);
		double crs;
		int ID;
		ID = BDSsateph.getCRS(crs);
		double deln = BDSsateph.getdeln();
		string Bs = to_string(ID) + to_string(crs)+to_string(deln);
		if ((Bset.empty() || !Bset.count(Bs))&&crs) {
			Bset.insert(Bs);
			Bmmp.insert(bepair(ID, BDSsateph));
		}
		int a = 1;
	}
	nav.NAV_RINEX_CLEAR();

}
int All_Sats_EPH::Cal_OneSatinfo(int SatID, char mode, gtime_t obs_time, vector<double> & pos_vel, vector<double>  &bias_drift, double &var,int &svh,double &TGD) {
	//肩负着选择卫星的重要责任
	// 没有考虑IODE 而只考虑了时间最接近的卫星
	//GPS 星历有效时间7200 BDS 21600
	vector<double> pos(3, 0);
	vector<double> pos_tt(3, 0);
	double bias;
	double bias_tt;
	double tt = 1e-3;
	gtime_t time_tt = timeadd(obs_time, tt);
	if (mode == 'G') {
		int num = Gmmp.count(SatID);
		//int index = 0;
		int MAXdTOE = 7200;
		double MINdTOE = 7201;
		double dt = 0;
		Gsatephmmp::iterator it = Gmmp.find(SatID);
		Gsatephmmp::iterator it_c = it;
		for (int i = 0; i < num; i++,it++) {
			dt = fabs(timediff(obs_time,it->second.gettoe()));
			if (dt <= MAXdTOE && dt <= MINdTOE) { // 如果出现两个一样的时刻，信任后面的nav文件
				MINdTOE = dt;
				it_c = it;
			}
		}
		if (MINdTOE > (double)MAXdTOE) {
			cout << SatID << "号GPS卫星观测量" << "距离最近的卫星星历也有2小时以上的差异" << endl;
			return 0;
		}
		else {
			svh  = it_c->second.getSVH();
			var = it_c->second.cal_POSVAR();
			TGD = it_c->second.getTGD();
			it_c->second.Cali_GPS_Satinfo(obs_time, pos, bias);
			it_c->second.Cali_GPS_Satinfo(time_tt, pos_tt, bias_tt);
		}
	
	}
	else {
		int num = Bmmp.count(SatID);
		int MAXdTOE = 21600;
		double MINdTOE = 21601;
		double dt = 0;
		Bsatephmmp::iterator it = Bmmp.find(SatID);
		Bsatephmmp::iterator it_b = it;
		for (int i = 0; i < num; i++, it++) {
			dt = fabs(timediff(obs_time, it->second.gettoe()));
			if (dt < MAXdTOE && dt < MINdTOE) {
				MINdTOE = dt;
				it_b = it;
			}
		}
		if (MINdTOE > (double)MAXdTOE) {
			cout << SatID << "号BDS卫星观测量" << "距离最近的卫星星历也有6小时以上的差异" << endl;
			return 0;
		}
		else {
			svh = it_b->second.getSVH();
			var = it_b->second.cal_POSVAR();
			TGD = it_b->second.getTGD();
			/*if (it_b->first == 4) {
				int aaa = 1;
			}*/
			it_b->second.Cali_BDS_Satinfo(SatID, obs_time, pos, bias);
			it_b->second.Cali_BDS_Satinfo(SatID, time_tt, pos_tt, bias_tt);
			
		}
	}
	for (int i = 0; i < 3; i++) {
		pos_vel[i] = pos[i];
		pos_vel[i + 3] = (pos_tt[i] - pos[i]) / tt;
	}
	bias_drift[0] = bias;
	bias_drift[1] = (bias_tt- bias) / tt;
	return 1;



}