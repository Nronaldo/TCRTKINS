#include "sins_TC.h"
Eigen::Matrix3d I33 = Eigen::Matrix3d::Identity();
Eigen::Matrix3d O33 = Eigen::Matrix3d::Zero();

CEarth::CEarth()
{
	Re = 6378137;
	f = 1 / 298.257223563;
	Rp = (1 - f)*Re;
	e1 = sqrt(Re * Re - Rp*Rp) / Re;
	e2 = sqrt(Re * Re - Rp*Rp) / Rp;
	wie = 7.2921151467e-5;
	g0 = 9.7803267714;
}
void CEarth::Update(const Eigen::Vector3d& pos, const Eigen::Vector3d &vel)
{
	double B = pos[0];
	double L = pos[1];
	double h = pos[2];
	double ve = vel[0];
	double vn = vel[1];
	double vu = vel[2];
	RN = Re / sqrt(1 - e1 * e1 * sin(B) * sin(B));
	RM = RN*(1 - e1 * e1) / (1 - e1 * e1 * sin(B) * sin(B));
	secL = 1.0/cos(B);
	tanL  = tan(B);
	sinL = sin(B);
	cosL = cos(B);
	sin2L = sin(2 * B);
	cos2L = cos(2 * B);
	Mpv2 = secL / (RN + h);
	Mpv4 = 1 / (RM + h);
	RNh = RN + h;
	RMh = RM + h;
	wnie = { 0,wie*cosL,wie*sinL };
	wnen = {-vn/RMh,ve/RNh,ve/RNh*tanL};
	wnin = wnie + wnen;
	double sinL2 = sinL * sinL,sinL4 = sinL2*sinL2;
	g = g0*(1 + 5.27094e-3*sinL2 + 2.32718e-5*sinL4) - 3.086e-6*h;
	gn[2] = -(g0*(1 + 5.27094e-3*sinL2 + 2.32718e-5*sinL4) - 3.086e-6*pos[2]);
	gcc = gn - (wnie + wnin)*vn;
	// gcc:gravitational acceleration / coriolis acceleration / centripetal acceleration
	wnien = 2 * wnie + wnen;
	gcc = vel.cross(wnien) + gn;
}
/*CIMU*/
void CIMU::update(const vector<double > &imu) {
	for (int i = 0;i<3;i++){
		dw[i] = imu[i];
		dv[i] = imu[i+3];
	}
	return;
}

/*CSINS*/
void CSINS::ins_init(const Eigen::Vector3d &att0, const Eigen::Vector3d &vel0, const Eigen::Vector3d &pos0) {
	nt = 1.0/TC_conf.sample_rate;
	att = att0;
	vel = vel0;
	pos << pos0[0]/180*PI, pos0[1] / 180 * PI,pos0[2];
	acc =bg=ba=wib=fb=web ={ 0,0,0 };
	Cnb = att2Cnb(att);
	Kg =Ka= Eigen::Matrix3d::Identity();
	tauA = tauG={ INF,INF,INF };
	eth.Update(pos,vel);
	Mpv << 0,		eth.Mpv4,	 0,
		   eth.Mpv2,0,			 0,
		   0,		0,			 1;
	fn = -eth.gn;
	lever <<TC_conf.lever[0], TC_conf.lever[1], TC_conf.lever[2];
	//old_imud = Eigen::VectorXd(6,0);
	/*
	x = Eigen::VectorXd(15,0);
	x<<att,
	   vel,
	   pos,
	   bg,
	   ba;
	*/
	
	Eigen::VectorXd P_diag = Eigen::VectorXd::Zero(15);
	for(int i = 0;i<3;i++){
		P_diag[i] = pow(TC_conf.init_att_unc[i],2);
		P_diag[i+3] = pow(TC_conf.init_vel_unc[i],2);
		P_diag[i+6] = pow(TC_conf.init_pos_unc[i],2);
		P_diag[i+9] = pow(TC_conf.init_bg_unc[i],2);
		P_diag[i+12] = pow(TC_conf.init_ba_unc[i],2);
	}
	P = Eigen::MatrixXd(P_diag.asDiagonal());
	Pa = Eigen::MatrixXd::Zero(15,15);
	xa = Eigen::VectorXd::Zero(15);
	Eigen::VectorXd Q_diag = Eigen::VectorXd::Zero(15);
	for(int i = 0;i<3;i++){
		Q_diag[i] = TC_conf.psd_gyro*nt;
		Q_diag[i+3] = TC_conf.psd_acce*nt;
		Q_diag[i+9] = TC_conf.psd_bg*nt;
		Q_diag[i+12] = TC_conf.psd_ba*nt;
	}
	Q = Eigen::MatrixXd(Q_diag.asDiagonal());
	//matrixtofile(Q, "Q", 'N');
	// xa = Eigen::VectorXd(15,0);
	// Pa = Eigen::MatrixXd(15,15);
	update_trans_mat();
	
	int aaaa = 1;
}
void CSINS::update_trans_mat() { //update_trans_mat
	double tl = eth.tanL, secl = eth.secL, secl2 = secl * secl,
		wN = eth.wnie[1], wU = eth.wnie[2], vE = vel[0], vN = vel[1];
	double f_RMh = 1.0/eth.RMh, f_RNh = 1/eth.RNh, f_clRNh = secl/eth.RNh,
		f_RMh2 = f_RMh * f_RMh, f_RNh2 = f_RNh * f_RNh;
	Eigen::Matrix3d Avn, Mp1, Mp2, Maa, Mav, Map, Mva, Mvv, Mvp, Mpp;
	Eigen::MatrixXd Ft(15,15);
	Avn = askew(vel);
	Mp1<<0, 0, 0, 
		-wU, 0, 0, 
		wN, 0, 0;
	Mp2<<0, 0, vN*f_RMh2, 
		0, 0, -vE * f_RNh2, 
		vE*secl2*f_RNh, 0, -vE * tl*f_RNh2;
	Maa = askew(-eth.wnin);
	Mav << 0, -f_RMh, 0, 
		  f_RNh, 0, 0, 
		tl*f_RNh, 0, 0;
	Map = Mp1 + Mp2;
	Mva = askew(fn);
	Mvv = Avn * Mav - askew(eth.wnie + eth.wnin);
	Eigen::Matrix3d Mvp_pre = Mp1 + Map;
	
	double scl = eth.sinL*eth.cosL;
	Mvp_pre(2, 0) -=  eth.get_g0()*(5.27094e-3 * 2 * scl + 2.32718e-5 * 4 * eth.sin2L*eth.cos2L);
	Mvp_pre(2, 2) += 3.086e-6;
	Mvp = Avn * (Mvp_pre);
	Mpv <<0, f_RMh, 0, 
		f_clRNh, 0, 0, 
		0, 0, 1;
	Mpp<< 0, 0, -vN * f_RMh2, 
		vE*tl*f_clRNh, 0, -vE * secl*f_RNh2,
		0, 0, 0;
	Eigen::Matrix3d Zeromat3 = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d bgmat3,bamat3;
	if (tauG[0] != INF) {
		Eigen::Vector3d tauG_temp;
		for (int i = 0; i < 3; i++) {
			tauG_temp[i] = 1.0 / tauG[i];
		}
		bgmat3 = Eigen::Matrix3d(tauG_temp.asDiagonal());
	}
	else {
		bgmat3 = Eigen::Matrix3d::Zero();
	}
	if (tauA[0] != INF) {
		Eigen::Vector3d tauA_temp;
		for (int i = 0; i < 3; i++) {
			tauA_temp[i] = 1.0 / tauA[i];
		}
		bamat3 = Eigen::Matrix3d(tauA_temp.asDiagonal());
	}
	else {
		bamat3 = Eigen::Matrix3d::Zero();
	}
	//Ft << Maa, Mav, Map, -Cnb, O33;
	Ft << Maa,		Mav,		Map,	-Cnb,		O33,
		Mva,		Mvv,		Mvp,     O33,		Cnb,
		O33,		Mpv,		Mpp,	 O33,		O33,
		O33,		O33,		O33,	 bgmat3,	O33,
		O33,		O33,		O33,	O33,		bamat3;
	Eigen::MatrixXd I_1515 = Eigen::MatrixXd::Identity(15,15);
	Phi = I_1515 + Ft * nt;
}
void CSINS::Update() {
	static Eigen::Vector3d old_dw = { 0,0,0 }; //单子样加前一周期
	static Eigen::Vector3d old_dv = { 0,0,0 };
	Eigen::Vector3d dw0 = imu.get_dw();
	Eigen::Vector3d dv0 = imu.get_dv();
	Eigen::Vector3d dw, dv;
	// correct bias and scaling factor errors
	dw = Kg * dw0 - bg * nt;
	dv = Ka*dv0 - ba*nt;
	//extrapolate velocity and position
	Eigen::Vector3d vel_mid, pos_mid;
	vel_mid = vel + acc*(nt / 2);
	pos_mid = pos + Mpv*(vel + vel_mid) / 2 * nt / 2; //这里不用再÷2吗
	// update the related parameters
	eth.Update(pos_mid, vel_mid);
	wib = dw / nt;
	fb = dv / nt;
	fn = Cnb*fb;
	web = wib - Cnb.transpose()*eth.wnie;
	// update velocity
	Eigen::Vector3d dv_rot, dv_scul, dv_sf, dv_cor, vel_new, pos_new;
	dv_rot = 0.5*dw0.cross(dv0);
	dv_scul = 1 / 12 * (old_dw.cross(dv0) + old_dv.cross(dw0));
	dv_sf = (Eigen::Matrix3d::Identity() - 0.5*nt*askew(eth.wnin))*Cnb*dv + Cnb*(dv_rot + dv_scul);
	dv_cor = eth.gcc*nt;
	vel_new = vel + dv_sf + dv_cor;

	// update position
	Mpv(1,0) = eth.Mpv2;
	Mpv(0,1) = eth.Mpv4;
	pos_new = pos + Mpv*(vel + vel_new) / 2 * nt;

	// update attitude
	Eigen::Vector3d dw_cone, phi_b_ib, phi_n_in, att_new;
	Eigen::Matrix3d Cbb, Cnn, Cnb_new;
	dw_cone = 1 / 12 * old_dw.cross(dw0);
	phi_b_ib = dw + dw_cone;
	phi_n_in = eth.wnin*nt;
	Cbb = rvec2mat(phi_b_ib);
	Cnn = rvec2mat(phi_n_in).transpose();
	Cnb_new = Cnn * Cnb*Cbb;
	att_new = Cnb2att(Cnb_new);
	//update INS result
	Cnb = Cnb_new;
	att = att_new;
	vel = vel_new;
	pos = pos_new;
	//x << att,
	//	 vel,
	//	 pos,
	//	 bg,
	//	 ba;
	old_dw = dw0;
	old_dv = dv0;
}
void CSINS::ins_Error_time_update() {
	// The TC kalman includes [sins error states]15*1 and [SD ambiguity]Num_sat*1
	// The time update procedure contains two parts: ins_error_time_update and SD ambiguity updates
	//The two parts are operated in different rates.
	update_trans_mat();
	Eigen::MatrixXd G = Eigen::MatrixXd::Zero(15,15);
	G.block(0, 0, 3, 3) << -Cnb;
	G.block(3, 3, 3, 3) << Cnb;
	for (int i = 0; i < 6; i++)
		G(9 + i, 9 + i) = 1;   
	Eigen::MatrixXd Q0 = G * Q*G.transpose();// Q0 equals to Q if Q is diagonal
	P = Phi * P*Phi.transpose() + Q0;
	//matrixtofile(Phi,"INS_Phi",'N');
	//matrixtofile(P, "INS_P", 'N');
	int aaaa = 1;
	return;
}
void CSINS::feedback(const Eigen::VectorXd & x_fb, const double * State_P,  int State_P_dimension) {
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	Cnb = (I + askew(x_fb.segment(0, 3)))*Cnb;
	att = Cnb2att(Cnb);
	vel -= x_fb.segment(3,3);
	pos -= x_fb.segment(6, 3);
	bg += x_fb.segment(9, 3);
	ba += x_fb.segment(12, 3);
	for (int i = 0; i < 15; i++) {
		for (int j = 0; j < 15; j++) {
			P(i, j) = State_P[i + j * State_P_dimension];
		}
	}
	//matrixtofile(P,"SINS_P",'N');
	int aaaa = 1;
	return;
}
void  CSINS::feedback(const Eigen::VectorXd & x_fb, const Eigen::MatrixXd &State_P) {
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d Cnb_xa;
	Cnb_xa = (I + askew(x_fb.segment(0, 3)))*Cnb;
	xa.segment(0,3) = Cnb2att(Cnb);
	xa.segment(3, 3) = vel - x_fb.segment(3, 3);
	xa.segment(6, 3) = pos - x_fb.segment(6, 3);
	xa.segment(9, 3) = bg + x_fb.segment(9, 3);
	xa.segment(12, 3) = ba + x_fb.segment(12, 3);
	Pa = State_P;
	//matrixtofile(P,"SINS_P",'N');
	int aaaa = 1;
	return;

}


/*Transformation functions*/
Eigen::Matrix3d att2Cnb(const Eigen::Vector3d &att) {
	double	si = sin(att[0]), ci = cos(att[0]),
		sj = sin(att[1]), cj = cos(att[1]),
		sk = sin(att[2]), ck = cos(att[2]);
	Eigen::Matrix3d Cnb;
	Cnb(0, 0) = cj * ck - si * sj*sk;	Cnb(0, 1) = -ci * sk;	Cnb(0, 2) = sj * ck + si * cj*sk;
	Cnb(1, 0) = cj * sk + si * sj*ck;	Cnb(1, 1) = ci * ck;	Cnb(1, 2) = sj * sk - si * cj*ck;
	Cnb(2, 0) = -ci * sj;				Cnb(2, 1) = si;		Cnb(2, 2) = ci * cj;
	return Cnb;
}
Eigen::Vector3d Cnb2att(const Eigen::Matrix3d &Cnb) {
	Eigen::Vector3d att;
	att[0] = asinEx(Cnb(2, 1));
	att[1] = atan2Ex(-Cnb(2, 0), Cnb(2, 2));
	att[2] = atan2Ex(-Cnb(0, 1), Cnb(1, 1));
	return att;
}
Eigen::Matrix3d askew(const Eigen::Vector3d &v) {
	Eigen::Matrix3d CMat3;
	CMat3<<0, -v[2], v[1],
		v[2], 0.0, -v[0],
		-v[1], v[0], 0;
	return CMat3;

}
Eigen::Matrix3d rvec2mat(const Eigen::Vector3d &rot_vector) {
	double theta = rot_vector.norm();
	Eigen::Vector3d norm_rot_vector = rot_vector / theta;
	Eigen::Matrix3d As_nrv = askew(norm_rot_vector);
	Eigen::Matrix3d	R;
	R = I33 + sin(theta)*As_nrv + (1 - cos(theta))*As_nrv*As_nrv;
	return R;
}
// MATH function
// determine the sign of 'val' with the sensitivity of 'eps'
int sign(double val, double eps)
{
	int s;

	if (val < -eps)
	{
		s = -1;
	}
	else if (val > eps)
	{
		s = 1;
	}
	else
	{
		s = 0;
	}
	return s;
}
// set double value 'val' between range 'minVal' and 'maxVal'
double range(double val, double minVal, double maxVal)
{
	double res;

	if (val < minVal)
	{
		res = minVal;
	}
	else if (val > maxVal)
	{
		res = maxVal;
	}
	else
	{
		res = val;
	}
	return res;
}
double atan2Ex(double y, double x)
{
	double res;

	if ((sign(y) == 0) && (sign(x) == 0))
	{
		res = 0.0;
	}
	else
	{
		res = atan2(y, x);
	}
	return res;
}
double asinEx(double x) {
	return asin(range(x, -1.0, 1.0));
}
double acosEx(double x) {
	return acos(range(x, -1.0, 1.0));
}
/*Test 函数*/
void Test_CEarth() {
	CEarth eth;
	string filename = "../Testfile/eth.txt";
	Eigen::Vector3d pos, vel;
	eth.Load_CEarth(filename,pos,vel);
	eth.Update(pos, vel);
	return;
}
void Test_CSINS() {
	string filenameavp = "../Testfile/avp.txt";
	Eigen::Vector3d att, vel, pos;
	Load_avp(filenameavp,att,vel,pos);
	CSINS sins;
	sins.ins_init(att,vel,pos);// ins_init
	sins.Load_CSINS();
	sins.Update();
	//matrixtofile(sins.Cnb, "Cnb", 'N');
	//matrixtofile(sins.Kg, "Kg", 'N');
	//matrixtofile(sins.Ka, "Ka", 'N');
	int a = 1;
}
void CSINS::Load_CSINS() {
	string filename_sins = "../Testfile/ins.txt";
	string filename_eth = "../Testfile/eth.txt";
	Eigen::Vector3d pos, vel;
	eth.Load_CEarth(filename_eth, pos, vel);
	ifstream fin(filename_sins);
	if (!fin.is_open()) {
		cout << "Fail to open the csins file." << endl;
		return;
	}
	string line, realname;
	int i = 0;
	vector<double> imu_vec(6,0);
	while (getline(fin, line)) {
		stringstream iss(line);
		i = 0;
		string temp_str;
		while (getline(iss, temp_str, ',')) {
			if (i == 0) {
				realname = temp_str;
				i++;
			}
			else {
				if (realname == "imu") {
					imu_vec[i++ - 1] = stod(temp_str);
				}
				else if (realname == "att") {
					att[i++ - 1] = stod(temp_str);
				}
				else if (realname == "vel") {
					vel[i++ - 1] = stod(temp_str);
				}
				else if (realname == "pos") {
					pos[i++ - 1] = stod(temp_str);
				}
				else if (realname == "acc") {
					acc[i++ - 1] = stod(temp_str);
				}
				else if (realname == "bg") {
					bg[i++-1] = stod(temp_str);
				}
				else if (realname == "ba") {
					ba[i++ - 1] = stod(temp_str);
				}
				else if (realname == "tauG") {
					tauG[i++ - 1] = INF;
				}
				else if (realname == "tauA") {
					tauA[i++ - 1] = INF;
				}
				else if (realname == "wib") {
					wib[i++ - 1] = stod(temp_str);
				}
				else if (realname == "fb") {
					fb[i++ - 1] = stod(temp_str);
				}
				else if (realname == "web") {
					web[i++ - 1] = stod(temp_str);
				}
				else if (realname == "fn") {
					fn[i++ - 1] = stod(temp_str);
				}
				else if (realname == "lever") {
					lever[i++ - 1] = stod(temp_str);
				}
				else if (realname == "Cnb") {
					Cnb(i++ - 1) = stod(temp_str);
				}
				else if (realname == "Kg") {
					Kg(i++ - 1) = stod(temp_str);
				}
				else if (realname == "Ka") {
					Ka(i++ - 1) = stod(temp_str);
				}
			}

		}

	}
	Eigen::Matrix3d mat = Cnb.transpose();
	Cnb = mat;
	mat = Kg.transpose();
	Kg = mat;
	mat = Ka.transpose();
	Ka = mat;
	imu.update(imu_vec);
	fin.close();
}
void CSINS::lever_trans(int mode) { // 0 GNSS->INS 1 INS->GNSS
	lever_trans(mode, pos, vel);
	/*int notion = mode ? 1 : -1;
	Eigen::Vector3d pos_modify = notion * Mpv*Cnb*lever;
	Eigen::Vector3d vel_modify = notion *Cnb*askew(web)*lever;
	pos += pos_modify;
	vel += vel_modify;*/
	return;
}
void CSINS::lever_trans(int mode, Eigen::Vector3d &Ext_pos, Eigen::Vector3d &Ext_vel) {
	int notion = mode ? 1 : -1;
	Eigen::Vector3d pos_modify = notion * Mpv*Cnb*lever;
	Eigen::Vector3d vel_modify = notion * Cnb*askew(web)*lever;
	if (mode == 1||mode == 0) {
		Ext_pos = pos + pos_modify;
		Ext_vel = vel + vel_modify;
	}
	else if (mode == 2) {
		Ext_pos = xa.segment(6,3) + pos_modify;
		Ext_vel = xa.segment(3, 3) + vel_modify;
	}
	return;
}
void CEarth::Load_CEarth(const string &filename, Eigen::Vector3d&pos, Eigen::Vector3d&vel) {
	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Fail to open the Earth file." << endl;
		return;
	}
	string line,realname;
	int i = 0;
	while (getline(fin,line)) {
		stringstream iss(line);
		i = 0;
		string temp_str;
		while (getline(iss, temp_str, ',')) {
			if (i == 0) {
				realname = temp_str.substr(3, temp_str.size() - 3);
				i++;
			}
			else {
				if (realname == "posllh") {
					pos[i++ - 1] = stod(temp_str);
				}
				else if (realname == "venu") {
					vel[i++ - 1] = stod(temp_str);
				}
				else if (realname == "RN") {
					RN = stod(temp_str);
				}
				else if (realname == "RM") {
					RM = stod(temp_str);
				}
				else if (realname == "Mpv2") {
					Mpv2 = stod(temp_str);
				}
				else if (realname == "Mpv4") {
					Mpv4 = stod(temp_str);
				}
				else if (realname == "RNh") {
					RNh = stod(temp_str);
				}
				else if (realname == "RMh") {
					RMh = stod(temp_str);
				}
				else if (realname == "tanL") {
					tanL = stod(temp_str);
				}
				else if (realname == "secL") {
					secL = stod(temp_str);
				}
				else if (realname == "sinL") {
					sinL = stod(temp_str);
				}
				else if (realname == "cosL") {
					cosL = stod(temp_str);
				}
				else if (realname == "sin2L") {
					sin2L = stod(temp_str);
				}
				else if (realname == "cos2L") {
					cos2L = stod(temp_str);
				}
				else if (realname == "g") {
					g = stod(temp_str);
				}
				else if (realname == "wnie") {
					wnie[i++ - 1] = stod(temp_str);
				}
				else if (realname == "wnen") {
					wnen[i++ - 1] = stod(temp_str);
				}
				else if (realname == "wnin") {
					wnin[i++ - 1] = stod(temp_str);
				}
				else if (realname == "gn") {
					gn[i++ - 1] = stod(temp_str);
				}
				else if (realname == "wnien") {
					wnien[i++ - 1] = stod(temp_str);
				}
				else if (realname == "gcc") {
					gcc[i++ - 1] = stod(temp_str);
				}
			}
		
		}
	
	}
	fin.close();
}
void Load_avp(const string& filenameavp, Eigen::Vector3d &att, Eigen::Vector3d &vel, Eigen::Vector3d &pos) {
	ifstream fin(filenameavp);
	if (!fin.is_open()) {
		cout << "Fail to open the avp file." << endl;
		return;
	}
	string line;
	while (getline(fin,line)) {
		stringstream iss(line);
		string temp,name;
		int i = 0;
		while (getline(iss, temp, ',')) {
			if (i == 0) {
				name = temp;
				i++;
			}
			else {
				if (name == "att") {
					att[i++ - 1] = stod(temp);
				}
				else if (name == "vel"){
					vel[i++ - 1] = stod(temp);
				}
				else {
					pos[i++ - 1] = stod(temp);
				}
			}
		}
	}
	fin.close();
}
