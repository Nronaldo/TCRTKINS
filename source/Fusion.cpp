#include "Fusion.h"
void Fusion::ins_align(ofstream &out) {

	if (ins_align_flag)//完成对齐了，就不再对齐了，
		return;

	Eigen::Vector3d vel;
	ins_align_flag = rcv.tdcp2vel(vel);
	
	rcv.GNSS_sever(out);
	if (rcv.get_sol_state()&& ins_align_flag) {
		Eigen::Vector3d pos(rcv.get_sol_pos(1).data());//llh
		if (fabs(vel(0)) < 1e-4)
			vel(0) = 1e-4;
		double yaw = -atan2(vel(0), vel(1));
		Eigen::Vector3d att = {0,0,yaw};
		sins.ins_init(att,vel,pos);
		sins.lever_trans(0);
		set_Solstate(2);
	}
	return;
}
void Fusion::ins_mech_ErrorTimeUpdate() {
	if (!ins_align_flag)
		return;
	sins.Update();
	sins.ins_Error_time_update();
	return;
}
int Fusion::rtk_ins() {
	//按目前的设计，进入本函数之前应该就已经把卫星信息更新完了
	numiterator++;
	cout << numiterator << endl;
	
	if (numiterator >= debug_stop) {
		int aaa = 1;
	}
	rcv.choosecommonsat();
	int niter = 1;
	int nv_dim = 0;//观测量双差对数
	rcv.RCV_zdres('B',NULL);
	string RTKmode = "PL";
	if (rcv.validsatnum() <= 0) {
		cout << "no common satellites,do not run RTK" << endl;
		set_ngnsslock(0);
		return 0;
	}
	if (numiterator >= debug_stop) {
		int aaa = 1;
	}
	udstate_rtkins();
	if (numiterator >= debug_stop) {
		
		int aaa = 1;
	}
	// 与rtklib 不同的是，pos vel的更新都由惯导+lever提供
	// 此处已经转换过了
	//bias的更新是否一致？先待定
	double * RTK_fixX = (double *)calloc(rcv.get_nx(), sizeof(double));
	for (int ii = 0; ii < niter; ii++) {
		rcv.RCV_zdres('R',NULL);
		nv_dim = rcv.RCV_ddres(RTKmode, NULL);//注意是算上伪距 相位 两种观测量的双差对数，改过来了吧？单独的伪距或相位也可以
			//displayddres();
		if (nv_dim >= 1)
			rtkins_Measureupdate();//要在里面进行周跳检测和IGGmodel
		if (numiterator >= debug_stop) {
			//matrixtofile(INS_Error_x, "INSEx", 'N');
			int aaaa = 1;
		}
		if (!rcv.solve_validpos(RTKmode, NULL,1)) {
			/*if (numiterator >= 147)
				return 0;*/ 
			rcv.set_nfix(0);
			set_Solstate(3);
			cout << "有效卫星数太少" << endl;
			//sol.set_state(0);//可以相信别的方法，这里先设置成0
			return 0;
		}
		sins.feedback(INS_Error_x,INS_Sbias_P,TC_nx);//如果能过来，就校正惯导的
		//sol.set_state(2);能到这里的,至少是浮点解
		set_Solstate(2);
		ngnsslock++;
		vector<int > Sbias_ix;
		if (numiterator >= debug_stop) {
			//matrixtofile(INS_Error_x, "INSEx", 'N');
			int aaaa = 1;
		}
		int nb = rtkins_resamb(RTK_fixX, Sbias_ix);//内部 restabm 已经把RTK_fixX的位置和速度都对应了惯导的舒尔补解
		if (numiterator >= debug_stop) {
			int aaa = 1;
		}
		if (nb > 1) {
			int firstvalidGPSID = rcv.getSATrefGPSID();
			int  firstvalidBDSID = rcv.getSATrefBDSID();//注意，这里不能轻易改，因为solve_validpos会改变SATref的值

			/*if (numiterator >= 229 && numiterator <= 240) {
				firstvalidGPSID = 15;
			}*/
			
			if (rcv.solve_validpos("PL", RTK_fixX,0)) {
				rcv.add_nfix();
				if (rcv.get_nfix() >= TC_conf.minfix && TC_conf.armode ==3) {
					//if (++nfix >= 1) {
						if (numiterator >= debug_stop) {
							int aaaa = 1;
						}
						rtkins_holdamb(RTK_fixX, Sbias_ix,nb, firstvalidGPSID, firstvalidBDSID);
						if (numiterator >= debug_stop) {
							int aaaa = 1;
						}
				}
				set_Solstate(1);
				//sol.set_state(1);//待会调整
			}
		}
		int aaa = 1;
		delete RTK_fixX;
	}
}
void Fusion::udstate_rtkins() {
	
	Eigen::Vector3d GNSS_pos_llh, GNSS_pos_xyz;
	Eigen::Vector3d GNSS_vel_enu, GNSS_vel_xyz;
	Eigen::Matrix3d Cne;
	sins.lever_trans(1, GNSS_pos_llh, GNSS_vel_enu);// rad
	GNSS_pos_xyz = blh2xyz(GNSS_pos_llh, Cne);
	GNSS_vel_xyz = Cne.transpose()*GNSS_vel_enu;
	Fusion_set_RTKx_posvel(rcv, GNSS_pos_xyz, GNSS_vel_xyz);//vel，RTK需要ECEF，但是在RTKins模式下，不会用到这里传入的速度？
	if (numiterator >= debug_stop) {
		int aaa = 1;
	}
	rcv.udbias();
	Fusion_get_RTK_INS_P();//把rcv的P直接复制过来，Ins的P也复制过来，会影响互相关，如果某一个发生了周跳的话，需要清把和INS状态的互相关都清0
	return;
}
void Fusion::Fusion_get_RTK_INS_P() {
	Eigen::MatrixXd Sins_P = sins.get_P();
	/*把Timeupdate的惯导协方差更新复制过来*/
	//matrixtofile(Sins_P,"Sins_P",'N');
	for (int i = 0; i < 15; i++) {
		for (int j = 0; j < 15; j++) {
			INS_Sbias_P[i + j * TC_nx] = Sins_P(i, j);
		}
	}
	const double* GNSS_RTK_P = rcv.get_RTK_P();
	int GNSS_nx = rcv.get_nx();
	int GNSS_bias_Sindex = TC_conf.dynamics_on ? 9 : 3;
	/*udbias中可能导致单差清零，所以需要把单差的协方差也复制过来*/
	//此处(1)不影响bias和ins的互相关，(2)其实只需要重置相应的项，这样全部赋值一遍大大冗余了，但是，hhhh太懒了，
	//一种情况下会影响互相关，那就是sbias中的一个单差发生周跳，重置
	vector<int> Reset_Sbias_index = rcv.get_Reset_Sbias_index();
	vector<int> Reset_Sbias_index_inTC(Reset_Sbias_index);
	for (int i = 0; i < Reset_Sbias_index.size(); i++) {
		Reset_Sbias_index_inTC[i] = Reset_Sbias_index_inTC[i] + 15 - GNSS_bias_Sindex;
	}
	for (int i = 15; i < TC_nx; i++) {
		for (int j = 15; j < TC_nx; j++) {
			INS_Sbias_P[i + j * TC_nx] = GNSS_RTK_P[i-15+ GNSS_bias_Sindex+(j-15+ GNSS_bias_Sindex)*GNSS_nx];
		}
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
		
	for (int i = 0; i < Reset_Sbias_index_inTC.size(); i++) { // 重置所有的单差周跳互相关的惯导P
		for (int j = 0; j < 15; j++) {
			INS_Sbias_P[j + Reset_Sbias_index_inTC[i] * TC_nx] = 0.0;
			INS_Sbias_P[Reset_Sbias_index_inTC[i]+j * TC_nx] = 0.0;
		}
	}
	/*Eigen::MatrixXd INS_M = Eigen::MatrixXd::Zero(TC_nx, TC_nx);
	for (int i = 0; i < TC_nx; i++) {
		for (int j = 0; j < TC_nx; j++) {
			INS_M(i,j) = INS_Sbias_P[i + j * TC_nx] ;
		}
	}
	matrixtofile(INS_M,"INS_M",'N');*/
	int aaaa = 1;

}
void Fusion::rtkins_Measureupdate() {
	int ddpairnum_GPSP = 0, ddpairnum_GPSL = 0, ddpairnum_BDSP = 0, ddpairnum_BDSL = 0, ddpairnum = 0;
	vector<int > GBSatSat_PLdres_size = rcv.get_GBSat_PLdres_size();
	ddpairnum_GPSP = GBSatSat_PLdres_size[0];
	ddpairnum_GPSL = GBSatSat_PLdres_size[1];
	ddpairnum_BDSP = GBSatSat_PLdres_size[2];
	ddpairnum_BDSL = GBSatSat_PLdres_size[3];
	ddpairnum = (ddpairnum_GPSP + ddpairnum_GPSL + ddpairnum_BDSP + ddpairnum_BDSL);
	int numofSbias;
	if (ddpairnum < 1) {
		cout << "ddres pairs are less than 1" << endl;
		return;//在RTKLIB中即使nv数目很少，也会继续进行观测更新，但是最终若有效卫星数过少，会认定定位无效
	}
	Eigen::MatrixXd H_GNSS,H_TC;// 双差对数*状态量维度,会考虑所有当前不为0 的状态量
	Eigen::VectorXd x_GNSS,x_TC;
	Eigen::VectorXd v_GNSS = Eigen::VectorXd::Zero(ddpairnum);//y-h(x)
	Eigen::MatrixXd P_GNSS,P_TC;
	Eigen::MatrixXd R_GNSS = Eigen::MatrixXd::Zero(ddpairnum, ddpairnum);
	vector<int > ix;
	string mode = "PL";
	if (numiterator >= debug_stop) {
		int aaa = 1;
	}
	rcv.FormHxP(mode, H_GNSS, x_GNSS, P_GNSS, ddpairnum_GPSP, ddpairnum_GPSL, ddpairnum_BDSP, ddpairnum_BDSL, ix);//GNSS的VxP
	rcv.FormvR(mode, v_GNSS, R_GNSS, ddpairnum_GPSP, ddpairnum_GPSL, ddpairnum_BDSP, ddpairnum_BDSL);//GNSS 的vR 应该是和TC 一样的
	if (numiterator >= debug_stop) {
		int aaa = 1;
	}
	vector<int > TC_ix = rtkins_FormHxP(H_GNSS, x_GNSS, P_GNSS,H_TC,x_TC,P_TC,ix);
	if (TC_conf.ins_aid_IGG && ngnsslock > 10) {
		Robust_IGG(H_TC, v_GNSS, R_GNSS,P_TC);
	}
	if (numiterator >= debug_stop) {
		//matrixtofile(H_GNSS,"HGNSS",'N');
		int aaa = 1;
	}
	rtkins_filter_unzip(H_TC, x_TC, P_TC,v_GNSS, R_GNSS,TC_ix);//先调通基本，再加INS周跳和IGG-3
	INS_Error_x = x_TC.segment(0,15);//先保存，不校正，但是对单差的校正已经进行了
	Eigen::Vector3d Extpos_llh = sins.get_sins_pos()-x_TC.segment(6,3);
	Eigen::Matrix3d E;
	Eigen::Vector3d Extpos_xyz = blh2xyz(Extpos_llh, E);
	Eigen::Vector3d Extvel = { 0,0,0 };
	Fusion_set_RTKx_posvel(rcv, Extpos_xyz,Extvel);//因为还要用校正后的位置再算一遍zdres 和 ddres,验证效果（实际感觉目前没啥意义，因为肯定会赋值）
	//
}
vector<int >  Fusion:: rtkins_FormHxP( const Eigen::MatrixXd &H_gnss, const Eigen::VectorXd &x_gnss, const Eigen::MatrixXd &P_gnss,
							  Eigen::MatrixXd &H_TC,  Eigen::VectorXd &x_TC, Eigen::MatrixXd &P_TC,const vector<int > &ix) {
	int ddpairnum_GPSP = 0, ddpairnum_GPSL = 0, ddpairnum_BDSP = 0, ddpairnum_BDSL = 0, ddpairnum = 0;
	vector<int > GBSatSat_PLdres_size = rcv.get_GBSat_PLdres_size();
	ddpairnum_GPSP = GBSatSat_PLdres_size[0];
	ddpairnum_GPSL = GBSatSat_PLdres_size[1];
	ddpairnum_BDSP = GBSatSat_PLdres_size[2];
	ddpairnum_BDSL = GBSatSat_PLdres_size[3];
	int GNSS_Sbias_startindex = 3;//纯RTK的单差开始索引
	if (TC_conf.dynamics_on)
		GNSS_Sbias_startindex = 9;
	int H_TC_rows = H_gnss.rows(); 
	int H_TC_cols = 15 + x_gnss.rows() - GNSS_Sbias_startindex;
	int Sbias_dimension = H_TC_cols - 15;
	H_TC = Eigen::MatrixXd::Zero(H_TC_rows, H_TC_cols);
	x_TC = Eigen::VectorXd::Zero(H_TC_cols);
	P_TC = Eigen::MatrixXd::Zero(H_TC_cols, H_TC_cols);
	vector<double> roverpos_xyz(3,0);
	vector<double> roverpos_llh(3, 0);
	Eigen::Matrix3d Cne;
	for (int i = 0; i < 3; i++) {
		roverpos_xyz[i] = x_gnss[i];
	}
	ecef2pos(roverpos_xyz, roverpos_llh);
	xyz2enu(roverpos_llh, Cne);
	//matrixtofile(Cen, "Cen", 'N');
	Eigen::Vector3d Cnb_lever = sins.get_Cnb_lever();
	Eigen::Matrix3d C_phi = askew(Cne.transpose()*Cnb_lever);
	Eigen::Matrix3d T;
	Dblh2Dxyz(roverpos_llh,T);
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	/*组装H矩阵*/
	/*姿态*/
	H_TC.block(0, 0, ddpairnum_GPSL, 3) = -H_gnss.block(0, 0, ddpairnum_GPSL, 3)*C_phi;
	H_TC.block(ddpairnum_GPSL, 0, ddpairnum_GPSP, 3) = -H_gnss.block(ddpairnum_GPSL, 0, ddpairnum_GPSP, 3)*C_phi;
	H_TC.block(ddpairnum_GPSL+ ddpairnum_GPSP, 0, ddpairnum_BDSL, 3) = -H_gnss.block(ddpairnum_GPSL + ddpairnum_GPSP, 0, ddpairnum_BDSL, 3)*C_phi;
	H_TC.block(ddpairnum_GPSL + ddpairnum_GPSP+ ddpairnum_BDSL, 0, ddpairnum_BDSP, 3) =-H_gnss.block(ddpairnum_GPSL + ddpairnum_GPSP + ddpairnum_BDSL, 0, ddpairnum_BDSP, 3)*C_phi;
	/*位置*/
	H_TC.block(0, 6, ddpairnum_GPSL, 3) = -H_gnss.block(0, 0, ddpairnum_GPSL, 3)*T;
	H_TC.block(ddpairnum_GPSL, 6, ddpairnum_GPSP, 3) = -H_gnss.block(ddpairnum_GPSL, 0, ddpairnum_GPSP, 3)*T;
	H_TC.block(ddpairnum_GPSL + ddpairnum_GPSP, 6, ddpairnum_BDSL, 3) = -H_gnss.block(ddpairnum_GPSL + ddpairnum_GPSP, 0, ddpairnum_BDSL, 3)*T;
	H_TC.block(ddpairnum_GPSL + ddpairnum_GPSP + ddpairnum_BDSL, 6, ddpairnum_BDSP, 3) = -H_gnss.block(ddpairnum_GPSL + ddpairnum_GPSP + ddpairnum_BDSL, 0, ddpairnum_BDSP, 3)*T;
	/*相位*/
	H_TC.block(0, 15, ddpairnum_GPSL, Sbias_dimension) = H_gnss.block(0, GNSS_Sbias_startindex, ddpairnum_GPSL, Sbias_dimension);
	H_TC.block(ddpairnum_GPSL + ddpairnum_GPSP, 15, ddpairnum_BDSL, Sbias_dimension) = H_gnss.block(ddpairnum_GPSL + ddpairnum_GPSP, GNSS_Sbias_startindex, ddpairnum_BDSL, Sbias_dimension);
	//matrixtofile(H_TC,"H_TC",'N');
	/*组装P矩阵*/
	// 注意当前不为0的状态都可能会产生协方差 惯导误差和单差之间的协方差也应该取下来
	vector<int > TC_ix;
	for (int i = 0; i < 15; i++)
		TC_ix.push_back(i);
	for (int i = GNSS_Sbias_startindex; i < ix.size();i++) {
		TC_ix.push_back(ix[i]- GNSS_Sbias_startindex+15);
	}
	for (int i = 0; i < TC_ix.size(); i++) {
		for (int j = 0; j < TC_ix.size(); j++) {
			P_TC(i, j) = INS_Sbias_P[TC_ix[i] + TC_ix[j] * TC_nx];
		}
	}
	
	int aaaa = 1;
	return TC_ix;
}
void Fusion::Robust_IGG(Eigen::MatrixXd & H,Eigen::VectorXd &v, Eigen::MatrixXd &R,const Eigen::MatrixXd &P) {
	//IGG-3 model
	//需要注意这里只针对伪距呢
	double c0 = TC_conf.IGGC0;
	double c1 = TC_conf.IGGC1;
	Eigen::MatrixXd Q = H * P*H.transpose()+R;
	Eigen::VectorXd std_res = Eigen::VectorXd::Zero(v.size());
	Eigen::VectorXd rfact = Eigen::VectorXd::Ones(v.size());
	int ddpairnum_GPSP = 0, ddpairnum_GPSL = 0, ddpairnum_BDSP = 0, ddpairnum_BDSL = 0, ddpairnum = 0;
	vector<int > GBSatSat_PLdres_size = rcv.get_GBSat_PLdres_size();
	ddpairnum_GPSP = GBSatSat_PLdres_size[0];
	ddpairnum_GPSL = GBSatSat_PLdres_size[1];
	ddpairnum_BDSP = GBSatSat_PLdres_size[2];
	ddpairnum_BDSL = GBSatSat_PLdres_size[3];
	for (int i = 0; i < v.size(); i++) {
		std_res(i) = fabs(v(i)) / sqrt(Q(i, i));
	}
	int si = ddpairnum_GPSL;
	int ei = ddpairnum_GPSP + ddpairnum_GPSL;
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	int j = 2;//两星座
	vector<int> Remove_P_index;
	while (j) {
		for (int i = si; i < ei; i++) {
			if (std_res(i) <= c0) {
				continue;
			}
			else if (std_res(i) >= c0 && std_res(i) <= c1) {
				rfact(i) = fabs(std_res(i)) / c0 * pow(((c1 - c0) / (c1 - fabs(std_res(i)))), 2);
			}
			else {
				rfact(i) = -1;
				Remove_P_index.push_back(i);
			}
			R(i, i) *= rfact[i];
		}
		j--;
		si = v.size() - ddpairnum_BDSP;
		ei = v.size();
	}
	int Remove_Pnum = Remove_P_index.size();
	if (!Remove_Pnum) 
		return;
	Eigen::MatrixXd R_pre_copy = R;
	Eigen::MatrixXd H_pre_copy = H;
	Eigen::VectorXd v_pre_copy = v;
	if (numiterator >= debug_stop) {
		//matrixtofile(R_pre_copy, "IGGR_Copy", 'N');
		int aaaa = 1;

	}
	int H_pre_rows = H.rows();
	int H_cols = H.cols();
	int H_aft_rows = H_pre_rows - Remove_Pnum;
	R.resize(H_aft_rows, H_aft_rows);
	R.setZero();
	H.resize(H_aft_rows, H_cols);
	H.setZero();
	v.resize(H_aft_rows);
	v.setZero();
	j = 0;
	for (int i = 0; i < H_pre_rows; i++) {
		if (rfact[i] == -1) {
			j++;
			continue;
		}
		v[i - j] = v_pre_copy[i];
		H.row(i - j) << H_pre_copy.row(i);
	}
	int R_startindex = 0;
	j = 0;
	int dim = 0;
	/*for (int k = 0; k < Remove_Pnum; k++) {
		dim = k?Remove_P_index[k]- Remove_P_index[k-1]-1: Remove_P_index[k];
		R.block(R_startindex, R_startindex, dim, dim) = R_pre_copy.block(R_startindex+j, R_startindex+j, dim, dim);
		R_startindex= R_startindex+dim;
		j++;
	}
	dim = H_aft_rows - R_startindex;
	R.block(R_startindex, R_startindex, dim, dim) = R_pre_copy.block(R_startindex + j, R_startindex + j, dim, dim);*/
	//matrixtofile(R_pre_copy, "IGGR_Copy", 'N');
	int ji = 0;
	int jk = 0;
	for (int i = 0; i < H_pre_rows; i++) {
		jk = 0;
		if (rfact[i] == -1) {
			ji++;
			continue;
		}
		for (int k = 0; k < H_pre_rows; k++) {
			if (rfact[k] == -1) {
				jk++;
				continue;
			}
			R(i - ji, k - jk) = R_pre_copy(i, k);
		}
	}
	//matrixtofile(R, "IGGR", 'N');
	if (numiterator >= debug_stop) {
		//matrixtofile(R, "IGGR", 'N');
		int aaaa = 1;
	
	}
	return;
}
void Fusion::rtkins_filter_unzip(Eigen::MatrixXd &H, Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::VectorXd &v, Eigen::MatrixXd &R, const vector<int > &TC_ix) {
	/*matrixtofile(v,"Tv",'N');
	matrixtofile(P, "TP", 'N');
	matrixtofile(R, "TR", 'N');
	matrixtofile(H, "TH", 'N');*/
	Eigen::MatrixXd PHT, HPHT_add_R, K;
	Eigen::MatrixXd I;
	I.setIdentity(P.rows(), P.cols());
	PHT = P * H.transpose();
	HPHT_add_R = (H*P*H.transpose() + R);
	K = PHT * HPHT_add_R.inverse();
	x = K * v;//每次校正都回零，所以时间更新以后的INS误差先验x是0，单差认为不变，
	//matrixtofile(x, "Tx", 'N');
	//matrixtofile(K, "TK", 'N');
	if (numiterator >= debug_stop) {
		/*matrixtofile(v, "Tv", 'N');
		matrixtofile(x,"Tx",'N');
		matrixtofile(R, "TR", 'N');
		matrixtofile(H, "TH", 'N');*/
		/*matrixtofile(P, "TP", 'N');
		matrixtofile(K, "TK", 'N');
		matrixtofile(PHT, "TPHT", 'N');
		matrixtofile(HPHT_add_R, "THPHT_add_R", 'N');*/
		int aaa = 1;
	}
	P = (I - K * H)*P;
	for (int i = 0; i < TC_ix.size(); i++) { // 此处应该调整rcv的
		for (int j = 0; j < TC_ix.size(); j++) {
			INS_Sbias_P[TC_ix[i] + TC_ix[j] * TC_nx] = P(i, j);
		}
	}
	Fusion_set_RTKxP_bias(rcv,x,P,TC_ix);
	//sins.feedback(x);
	return;
}
/*与固定模糊度有关的*/
int Fusion::rtkins_resamb(double *RTK_fixX,vector<int > &Sbias_ix) {
	Eigen::MatrixXd D_TC,D_GNSS, P_GNSS,P_TC;
	Eigen::VectorXd x_TC,x_GNSS;
	int nb;
	vector<int> TC_ix;
	if ((nb = rcv.ddmat_FormxP(D_GNSS, x_GNSS, P_GNSS, Sbias_ix)) < TC_conf.optminfixsats - 1) {//ddmat_FormxP应该拆成formx 和 formP
		cout << "not enough valid double-differences " << endl;
		return -1;
	}
	int RTK_pvanum = TC_conf.dynamics_on?9:3;
	int Sbias_dimension = Sbias_ix.size();
	int TC_dimension = Sbias_dimension + 15;
	/*以下组建紧组合 舒尔补的 D x P，其中D x 可从D_GNSS x_GNSS推得，是分块的，P必须根据索引依次得到*/
	x_TC = Eigen::VectorXd::Zero(TC_dimension);
	x_TC.segment(0, 3) = sins.get_sins_att();
	x_TC.segment(3, 3) = sins.get_sins_vel();
	x_TC.segment(6, 3) = sins.get_sins_pos();
	x_TC.segment(9, 3) = sins.get_sins_bg();
	x_TC.segment(12, 3) = sins.get_sins_ba();
	x_TC.segment(15, Sbias_dimension) = x_GNSS.segment(RTK_pvanum, Sbias_dimension);
	int  D_Cols = D_GNSS.cols();
	int  D_Rows = D_GNSS.rows();
	D_TC = Eigen::MatrixXd::Zero(D_GNSS.rows() - RTK_pvanum + 15, D_Cols - RTK_pvanum+15);
	D_TC.block(0, 0, 15, 15) = Eigen::MatrixXd::Identity(15,15);
	D_TC.block(15, 15, D_TC.rows()-15, D_TC.cols() - 15) = D_GNSS.block(RTK_pvanum, RTK_pvanum, D_Rows- RTK_pvanum, D_Cols- RTK_pvanum);
	for (int i = 0; i < 15; i++)
		TC_ix.push_back(i);
	for(int i = 0;i< Sbias_dimension;i++)
		TC_ix.push_back(Sbias_ix[i]- RTK_pvanum +15);
	P_TC = Eigen::MatrixXd::Zero(TC_dimension, TC_dimension);
	for (int i = 0; i < TC_dimension; i++) {
		for (int j = 0; j < TC_dimension; j++) {
			P_TC(i, j) = INS_Sbias_P[TC_ix[i]+TC_ix[j]*TC_nx];
		}
	}
	/*开始舒尔补*/
	int cny = D_TC.rows();
	int cnx = D_TC.cols();
	int na = 15;
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
	y = D_TC * x_TC;
	DP = D_TC * P_TC;
	Qy = DP * D_TC.transpose();
	Qb = Qy.block(na, na, nb, nb);
	Qab = Qy.block(0, na, na, nb);
	Eigen::VectorXd yb = y.segment(na, nb);
	vector<double> s(2, 0);
	int i, j;
	if (!(info = lambda(nb, solnum, yb, Qb, b, s))) {
		/*if (numiterator == debug_stop) {
			arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
			int aaaa = 1;
		}*/
		float aa = s[0] > 0 ? (float)(s[1] / s[0]) : 0.0f;
		if (aa > 999.9)
			aa = 999.9f;
		rcv.set_sol_ratio(aa);
		Eigen::VectorXd x_fb(15);	
		Eigen::MatrixXd INS_P_a = sins.get_P();
		if (s[0] <= 0.0 || s[1] / s[0] >= TC_conf.thresar[0]) {
			//arraymattofile(RTK_P,nx,nx,0,0,na-1,na-1,'F',"P-a");
			//先删去bias
			//matrixtofile(RTK_fixPpva, "Pa", 'N');
			y.segment(na, nb) = y.segment(na, nb) - b.block(0, 0, nb, 1);// y是原始位置状态量加D矩阵转换的双差状态量，这里是用浮点双差减去整数双差 manual 式 E.7.19
			yb = y.segment(na, nb);
			db = Qb.inverse()*yb;
			x_fb =  - Qab * db;
			QQ = Qab * Qb.inverse();
			INS_P_a = INS_P_a - QQ * Qab.transpose();
			//if (numiterator == debug_stop) {
			//	/*matrixtofile(QQ,"QQ",'N');
			//	matrixtofile(RTK_fixPpva, "Pa", 'N');*/
			//	int aaaa = 1;
			//}
			//matrixtofile(x_fb, "x_fb", 'N');
			sins.feedback(x_fb, INS_P_a);
			Eigen::Vector3d Extpos_llh = sins.get_sins_pos();
			Eigen::Matrix3d Cne;
			Eigen::Vector3d Extpos_xyz = blh2xyz(Extpos_llh, Cne);
			Eigen::Vector3d Extvelxyz = Cne.transpose()*sins.get_sins_vel();
			Fusion_set_RTKxa_posvel(rcv, Extpos_xyz, Extvelxyz);
			rcv.restamb(RTK_fixX, Sbias_ix, b,1);//不可以省略，ddres中要用到相位单差
			//RTKLIB 和 GINAV都要用整数解算一次zdres 和 ddres然后再判断残差
			// 虽然这个固定解 不会参与kalman滤波，只做判断定位有效性用
		}
		else {
			rcv.set_nfix(0);
			nb = 0;
		}
	}
	else {
		//nb = 0;
		rcv.set_nfix(0);
		cout << "lambda error info=" << info << endl;
	}
	//outddint(info, ix,D,b);
	return nb;
}
void Fusion::rtkins_holdamb(double *RTK_xfix, const vector<int > &Sbias_ix, int nb, int firstvalidGPSID, int firstvalidBDSID) {
	Eigen::MatrixXd P_GNSS, R, H_GNSS;
	Eigen::VectorXd x, v;
	vector<int > all_stateindex;
	//Sbias_ix 只包含实际参与运算的单差状态 但是当前不为0的单差状态有的会有周跳，不需要构建对应H
	//但是仍需要得到校正，这部分单差会在all_stateindex里面
	int nv = rcv.holdamb_Form_xPRHv(RTK_xfix, Sbias_ix, nb, firstvalidGPSID, firstvalidBDSID, all_stateindex, x, H_GNSS, v, P_GNSS, R);
	if (nv <= 0) return;
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	/*根据all_stateindex 组装出rtk_ins的hold的P矩阵*/
	vector<int > TC_ix;
	int GNSS_Sbias_startindex = 3;//纯RTK的单差开始索引
	if (TC_conf.dynamics_on)
		GNSS_Sbias_startindex = 9;
	for (int i = 0; i < 15; i++)
		TC_ix.push_back(i);
	//TC_ix.insert(TC_ix.end(),all_stateindex.begin()+ GNSS_Sbias_startindex,all_stateindex.end());
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	for (int i = GNSS_Sbias_startindex; i < all_stateindex.size(); i++) {
		TC_ix.push_back(all_stateindex[i] - GNSS_Sbias_startindex + 15);
	}
	if (numiterator >= debug_stop) {
		int aaaa = 1;
	}
	int TC_ix_size = TC_ix.size();
	Eigen::MatrixXd P_TC = Eigen::MatrixXd::Zero(TC_ix_size, TC_ix_size);
	for (int i = 0; i < TC_ix_size; i++) {
		for (int j = 0; j < TC_ix_size; j++) {
			P_TC(i, j) = INS_Sbias_P[TC_ix[i] + TC_ix[j] * TC_nx];
		}
	}
	/*根据H_GNSS，在前面补0*/
	int H_GNSS_rows = H_GNSS.rows();
	int H_GNSS_cols = H_GNSS.cols();
	Eigen::MatrixXd H_TC = Eigen::MatrixXd::Zero(H_GNSS_rows, 15 - GNSS_Sbias_startindex+ H_GNSS_cols);
	H_TC << Eigen::MatrixXd::Zero(H_GNSS_rows,15- GNSS_Sbias_startindex),H_GNSS;
	rtkins_holdamb_Filter_FB(H_TC,v,P_TC,R,TC_ix);
}
void Fusion::rtkins_holdamb_Filter_FB(const Eigen::MatrixXd &H, const Eigen::VectorXd & v, const Eigen::MatrixXd & P, const Eigen::MatrixXd &R,
	const vector<int > 	&TC_ix) {
	rcv.set_holdamb_flag(1);
	Eigen::MatrixXd PHT, HPHT_add_R, K, P_hold;
	Eigen::MatrixXd I;
	Eigen::VectorXd x_fb;
	I.setIdentity(P.rows(), P.cols());

	PHT = P * H.transpose();
	//cout << "H" << H << endl;
	HPHT_add_R = (H*P*H.transpose() + R);
	//cout <<"R"<< R << endl;
	K = PHT * HPHT_add_R.inverse();
	//matrixtofile(K, "hold_K", 'N');
	//cout << "K" << K << endl;
	x_fb = K * v;
	P_hold = (I - K * H)*P;
	//matrixtofile(K, "hold_K", 'N');
	for (int i = 0; i < TC_ix.size(); i++) { // 此处应该调整rcv的
		for (int j = 0; j < TC_ix.size(); j++) {
			INS_Sbias_P[TC_ix[i] + TC_ix[j] * TC_nx] = P_hold(i, j);
		}
	}
	Fusion_set_RTKxP_bias(rcv, x_fb, P_hold, TC_ix);
	Eigen::MatrixXd P_INS = P_hold.block(0,0,15,15);
	sins.feedback(x_fb, INS_Sbias_P, TC_nx);//校正惯导的状态
	return;
}
void Fusion::TC_instofile(double soltime, ofstream& out) {
	static int head_flag = 0;
	if (!head_flag) {
		TC_instofile_head(out);
		head_flag = 1;
	}
	TC_instofile_time(soltime, out);
	if(TC_conf.out_mode[0])
		TC_instofile_pos(out);
	out << Sol_state << "  ";
	 out << setiosflags(ios::fixed) << setprecision(1)<< rcv.get_sol_ratio() << "  ";
	 if (rcv.get_sol_ratio() < 10)
		 out << " ";
	 if (rcv.get_sol_ratio() < 100)
		 out << " ";
	if (TC_conf.out_mode[1])
		TC_instofile_posvar(out);
	if (TC_conf.out_mode[2])
		TC_instofile_vel(out);
	if (TC_conf.out_mode[3])
		TC_instofile_velvar(out);
	if (TC_conf.out_mode[4])
		TC_instofile_att(out);
	if (TC_conf.out_mode[5])
		TC_instofile_attvar(out);
	if (TC_conf.out_mode[6])
		TC_instofile_bg(out);
	if (TC_conf.out_mode[7])
		TC_instofile_ba(out);
	out << endl;
}
void Fusion::TC_instofile_head(ofstream& out) {
	out << "% program   : TC_RTK_INS_Cpp v1"<<endl;
	out << "% I am sorry I only support GPS/BDS single Frequency RTK" << endl;
	if(!TC_conf.ins_mode)
		out << "% mode: RTK" << endl;
	else 
		out << "% mode: TC_RTK_INS" << endl;
	if (TC_conf.pos_mode)
		out << "% Posmode: INS" << endl;
	else
		out << "% Posmode: GNSS"  << endl;
	out << "% Outputting items(from the left to the right):"<<endl;
	out << "%";
	if (TC_conf.timef == 1) {
		out<<"    " << "Unix" << "          ";
	}
	else if (TC_conf.timef == 2) {
		out<<"    " << "week+sow" << "       ";
	}
	else { // YMDHMS
		out <<"     "<< "YMDHMS" << "               ";
	}
	if (TC_conf.out_mode[0]) {
		if (TC_conf.output_pos_format == "xyz") {
			out << " " <<"x-ecef(m)"<<"     "<< "y-ecef(m)"<<"      "<<"z - ecef(m)";
		}
		else if (TC_conf.output_pos_format == "llhrad") {
			out << "lat(rad)" << "     " << "longi(rad)" << "    " << "Hei(m)";
		}
		else if (TC_conf.output_pos_format == "llhdeg") {
			out << "lat(deg)" << "    " << "longi(deg)" << "    " << "Hei(m)";
		}
		else { //enu
			out  << "east(m)" << "     " << "north(m)" << "     " << "up(m)";
		}
	}
	out << "    Q";
	out << ","<<"ratio";
	if (TC_conf.out_mode[1]) {
		if (TC_conf.output_pos_format == "xyz") {
			out << "," << "sdx(m)" << "," << "sdy(m)" << "," << "sdz(m)" << "," << "sdxy(m)"<< "," << "sdyz(m)" << "," << "sdxz(m)";
		}
		else if (TC_conf.output_pos_format == "llhrad") {
			out << "," << "sdB(rad)" << "," << "sdL(rad)" << "," << "sdH(m)";
		}
		else if (TC_conf.output_pos_format == "llhdeg") {
			out << "," << "sdB(deg)" << "," << "sdL(deg)" << "," << "sdH(m)";
		}
		else { //enu
			out << "," << "sde(m)" << "," << "sdn(m)" << "," << "sdu(m)" << "," << "sden(m)" << "," << "sdnu(m)" << "," << "sdeu(m)";
		}
	}
	if (TC_conf.out_mode[2]) {
		if (TC_conf.output_pos_format == "xyz") {
			out << "," << "vel-x(m)" << "," << "vel-y(m)" << "," << "vel-z(m)" << ",";
		}
		else { //enu
			out << "," << "vel-e(m)" << "," << "vel-n(m)" << "," << "vel-u(m)";
		}
	}
	if (TC_conf.out_mode[3]) {
		if (TC_conf.output_pos_format == "xyz") {
			out << "," << "vsdx(m)" << "," << "vsdy(m)" << "," << "vsdz(m)" << "," << "vsdxy(m)" << "," << "vsdyz(m)" << "," << "vsdxz(m)";
		}
		else { //enu
			out << "," << "vsde(m)" << "," << "vsdn(m)" << "," << "vsdu(m)" << "," << "vsden(m)" << "," << "vsdnu(m)" << "," << "vsdeu(m)";
		}
	}
	if (TC_conf.out_mode[4]) {
		if (TC_conf.att_mode == 0) {
			out << "," << "pitch(rad)" << "," << "roll(rad)" << "," << "yaw(rad)" ;
		}
		else { //enu
			out << "," << "pitch(deg)" << "," << "roll(deg)" << "," << "yaw(deg)";
		}
	}
	if (TC_conf.out_mode[5]) {
		if (TC_conf.att_mode == 0) {
			out << "," << "pitchsd(rad)" << "," << "rollsd(rad)" << "," << "yawsd(rad)";
		}
		else { //enu
			out << "," << "pitchsd(deg)" << "," << "rollsd(deg)" << "," << "yawsd(deg)";
		}
	}
	if (TC_conf.out_mode[6]) {
			out << "," << "bgx(rad/s)" << "," << "bgy(rad/s)" << "," << "bgz(rad/s)";
	}
	if (TC_conf.out_mode[7]) {
		out << "," << "bax(m/s^2)" << "," << "bay(m/s^2)" << "," << "baz(m/s^2)";
	}
	out << endl;
}
void Fusion::TC_instofile_time(double soltime, ofstream& out) {
	if (TC_conf.timef == 1) {
		out << setiosflags(ios::fixed) << setprecision(3) << soltime << "   ";
	}
	else if (TC_conf.timef == 2) { // Week+Seconds of week
		int week;
		double sow;
		gtime_t timestamp(soltime);
		sow = time2gpst(timestamp, &week);
		out << week << " " << setiosflags(ios::fixed) << setprecision(3) << sow << "   ";
	}
	else { // YMDHMS
		gtime_t timestamp(soltime);
		string stime;
		time2str(timestamp, stime, 3);
		out << stime << "  ";
	}
}
void Fusion:: TC_instofile_pos(ofstream& out) {// 2 3 对应的解之间输出ins的值，1对应的解输出xa的值，xa不参与下一次浮点kalman
	Eigen::Vector3d posllh_rad, GNSS_vel_enu;
	if (Sol_state == 2 || Sol_state == 3) {
		if (TC_conf.pos_mode)
			posllh_rad = sins.get_sins_pos(); //输出惯导的位置
		else
			sins.lever_trans(1, posllh_rad, GNSS_vel_enu);//输出杆臂转换的GNSS的位置
	}
	else { 
		if (TC_conf.pos_mode)
			posllh_rad = sins.get_sins_xapos(); //输出整数解惯导的位置
		else
			sins.lever_trans(2, posllh_rad, GNSS_vel_enu);//输出杆臂转换的固定解GNSS的位置
	}
	
	Eigen::Vector3d posxyz;
	if (TC_conf.output_pos_format == "xyz") {
		Eigen::Matrix3d Cne;
		posxyz = blh2xyz(posllh_rad, Cne);
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << posxyz[i] << "   ";
		}
	}
	else if (TC_conf.output_pos_format == "llhrad") { 
		for (int i = 0; i < 2; i++) {
			out << setiosflags(ios::fixed) << setprecision(9) << posllh_rad[i] << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << posllh_rad[2] << "   ";
	}
	else if (TC_conf.output_pos_format == "llhdeg") {
		for (int i = 0; i < 2; i++) {
			out << setiosflags(ios::fixed) << setprecision(7) << posllh_rad[i]/PI*180 << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << posllh_rad[2] << "   ";
	}
	else { // ENU
		vector<double > posbasellhxyz(3, 0);
		ecef2pos(TC_conf.basepose, posbasellhxyz);// base 的位置转llh
		vector<double > e_rb(3, 0);
		Eigen::Matrix3d Cne;
		posxyz = blh2xyz(posllh_rad, Cne);
		for (int i = 0; i < 3; i++)
			e_rb[i] = posxyz[i] - TC_conf.basepose[i];
		vector<double > enu(3, 0);
		ecef2enu(posbasellhxyz, e_rb, enu);
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << enu[i] << "   ";
		}
	}
	
}
void Fusion::TC_instofile_posvar(ofstream& out) {
	Eigen::Matrix3d ins_pos_var;
	Eigen::Vector3d posllh_rad;
	if (Sol_state == 2 || Sol_state == 3) {
		ins_pos_var = sins.get_P().block(6, 6, 3, 3);// llh rad
		posllh_rad = sins.get_sins_pos();
	}
	else {
		ins_pos_var = sins.get_Pa().block(6, 6, 3, 3);// llh rad
		posllh_rad = sins.get_sins_xapos();
	}
	
	if (TC_conf.output_pos_format == "xyz") {
		Eigen::Matrix3d T;
		vector<double >ins_pos_llh_v(3, 0);
		for (int i = 0; i < 3; i++) {
			ins_pos_llh_v[i] = posllh_rad(i);
		}
		Dblh2Dxyz(ins_pos_llh_v,T);
		Eigen::Matrix3d pos_var_xyz = T * ins_pos_var*T.transpose();
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_xyz(i,i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_xyz(0, 1)) << "   ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_xyz(1, 2))<< "   ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_xyz(0, 2))<< "   ";
	}
	else if (TC_conf.output_pos_format == "llhrad") { 
		for (int i = 0; i < 2; i++) {
			out << setiosflags(ios::fixed) << setprecision(9) << sqvar(ins_pos_var(i, i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_pos_var(2, 2)) << "   ";
	}
	else if (TC_conf.output_pos_format == "llhdeg") {
		for (int i = 0; i < 2; i++) {
			out << setiosflags(ios::fixed) << setprecision(7) << sqvar(ins_pos_var(i, i))/PI*180 << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_pos_var(2, 2)) << "   ";
	}
	else { // ENU
		Eigen::Matrix3d T,Cne;
		vector<double >ins_pos_llh_v(3, 0);
		
		for (int i = 0; i < 3; i++) {
			ins_pos_llh_v[i] = posllh_rad(i);
		}
		Dblh2Dxyz(ins_pos_llh_v, T);
		Eigen::Matrix3d pos_var_xyz = T * ins_pos_var*T.transpose();
		blh2xyz(posllh_rad,Cne);
		Eigen::Matrix3d pos_var_enu = Cne * pos_var_xyz*Cne.transpose();
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_enu(i, i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_enu(0, 1)) << "   ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_enu(1, 2)) << "   ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(pos_var_enu(0, 2)) << "   ";
	}
}
void Fusion::TC_instofile_vel(ofstream& out) {
	Eigen::Vector3d posllh_rad, vel_enu;
	if (Sol_state == 2 || Sol_state == 3) {
		if (TC_conf.pos_mode) {
			posllh_rad = sins.get_sins_pos(); //输出惯导的位置，速度
			vel_enu = sins.get_sins_vel();
		}
		else
			sins.lever_trans(1, posllh_rad, vel_enu);//输出杆臂转换的GNSS的位置,速度
	}
	else {
		if (TC_conf.pos_mode) {
			posllh_rad = sins.get_sins_xapos(); //输出整数解惯导的位置,速度
			vel_enu = sins.get_sins_xavel();
		}
			
		else
			sins.lever_trans(2, posllh_rad, vel_enu);//输出杆臂转换的固定解GNSS的位置
	}
	if (TC_conf.output_pos_format == "xyz") {
		Eigen::Matrix3d Cne;
		blh2xyz(posllh_rad, Cne);
		Eigen::Vector3d vel_xyz = Cne.transpose()*vel_enu;
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(5) << vel_xyz[i] << "   ";
			if (vel_xyz[i] >= 0)
				out << " ";

		}
	}
	else { // ENU//速度，非xyz要求下均输出ENU的速度
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(5) << vel_enu[i] << "   ";
			if (vel_enu[i] >= 0)
				out << " ";
		}
	}

}
void Fusion::TC_instofile_velvar(ofstream& out) {
	Eigen::Matrix3d ins_vel_var;
	Eigen::Vector3d posllh_rad;
	if (Sol_state == 2 || Sol_state == 3) {
		ins_vel_var = sins.get_P().block(3, 3, 3, 3);// llh rad
		posllh_rad = sins.get_sins_pos();
	}
	else {
		ins_vel_var = sins.get_Pa().block(3, 3, 3, 3);// llh rad
		posllh_rad = sins.get_sins_xapos();
	}
	if (TC_conf.output_pos_format == "xyz") {
		Eigen::Matrix3d Cne,ins_vel_var_xyz;
		blh2xyz(posllh_rad, Cne);
		ins_vel_var_xyz = Cne.transpose()*ins_vel_var*Cne;
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var_xyz(i,i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var_xyz(0, 1)) << "   ";
		if (ins_vel_var_xyz(0, 1) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var_xyz(1, 2)) << "   ";
		if (ins_vel_var_xyz(1, 2) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var_xyz(0, 2)) << "   ";
		if (ins_vel_var_xyz(0, 2) >= 0)
			out << " ";
	}
	else { // ENU//速度，非xyz要求下均输出ENU的速度
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var(i,i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var(0, 1)) << "   ";
		if (ins_vel_var(0, 1) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var(1, 2)) << "   ";
		if (ins_vel_var(1, 2) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_vel_var(0, 2)) << "   ";
		if (ins_vel_var(0, 2) >= 0)
			out << " ";
	}
}
void Fusion::TC_instofile_att(ofstream& out) {
	Eigen::Vector3d att_rad;
	if (Sol_state == 2 || Sol_state == 3) {
		att_rad = sins.get_sins_att();
	}
	else {
		att_rad = sins.get_sins_xaatt();
	}
	if (TC_conf.att_mode == 0) { // rad 
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << att_rad(i) << "   ";
			if (att_rad(i) >= 0)
				out << " ";
		}
	}
	else {
		Eigen::Vector3d att_deg = att_rad/PI*180;
		att_deg(2) = -att_deg(2) + (att_deg(2) >= 0 ? 360 : 0);//严老师 北偏西为正 -pi~pi 航向角一般北偏东为正，0~2pi。详细见严老师讲义附录姿态角一章
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << att_deg(i) << "   ";
			if (att_deg(i) >= 0)
				out << " ";
		}
	}
}
void Fusion::TC_instofile_attvar(ofstream& out) {
	Eigen::Matrix3d ins_att_var;
	if (Sol_state == 2 || Sol_state == 3) {
		ins_att_var = sins.get_P().block(0, 0, 3, 3);// llh rad
	}
	else {
		ins_att_var = sins.get_Pa().block(0, 0, 3, 3);// llh rad
	}
	if (TC_conf.att_mode == 0) {
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var(i, i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var(0, 1)) << "   ";
		if (ins_att_var(0, 1) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var(1, 2)) << "   ";
		if (ins_att_var(1, 2) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var(0, 2)) << "   ";
		if (ins_att_var(0, 2) >= 0)
			out << " ";
	}
	else { // 转为deg
		Eigen::Matrix3d ins_att_var_deg = ins_att_var/PI*180/PI*180;
		for (int i = 0; i < 3; i++) {
			out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var_deg(i, i)) << "   ";
		}
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var_deg(0, 1)) << "   ";
		if (ins_att_var_deg(0, 1) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var_deg(1, 2)) << "   ";
		if (ins_att_var_deg(1, 2) >= 0)
			out << " ";
		out << setiosflags(ios::fixed) << setprecision(4) << sqvar(ins_att_var_deg(0, 2)) << "   ";
		if (ins_att_var_deg(0, 2) >= 0)
			out << " ";
	}
}
void Fusion::TC_instofile_bg(ofstream& out) {
	Eigen::Vector3d ins_bg;
	if (Sol_state == 2 || Sol_state == 3) {
		ins_bg = sins.get_sins_bg();// llh rad
	}
	else {
		ins_bg = sins.get_sins_xabg();// llh rad
	}
	for (int i = 0; i < 3; i++) {
		out << setiosflags(ios::fixed) << setprecision(11) << ins_bg(i) << "   ";
		if (ins_bg(i) >= 0)
			out << " ";
	}
}
void Fusion::TC_instofile_ba(ofstream& out) {
	Eigen::Vector3d ins_ba;
	if (Sol_state == 2 || Sol_state == 3) {
		ins_ba = sins.get_sins_ba();// llh rad
	}
	else {
		ins_ba = sins.get_sins_xaba();// llh rad
	}
	for (int i = 0; i < 3; i++) {
		out << setiosflags(ios::fixed) << setprecision(11) << ins_ba(i) << "   ";
		if (ins_ba(i) >= 0)
			out << " ";
	}
}
void Fusion::Fusion_plot() {
	static int flag = 0;
	static vector<double > posbasellhxyz(3, 0);
	
	Eigen::Vector3d posllh_rad, GNSS_vel_enu;
	if (Sol_state == 2 || Sol_state == 3) {
		if (TC_conf.pos_mode)
			posllh_rad = sins.get_sins_pos(); //画惯导的位置
		else
			sins.lever_trans(1, posllh_rad, GNSS_vel_enu);//画杆臂转换的GNSS的位置
	}
	else {
		if (TC_conf.pos_mode)
			posllh_rad = sins.get_sins_xapos(); //画整数解惯导的位置
		else
			sins.lever_trans(2, posllh_rad, GNSS_vel_enu);//画杆臂转换的固定解GNSS的位置
	}
	
	if (!flag) {
		plt::clf();
		plt::title("The trajectory.");
		ecef2pos(TC_conf.basepose, posbasellhxyz);// base 的位置转llh	
		flag++;
	}
	Eigen::Vector3d posxyz;
	vector<double > e_rb(3, 0);
	Eigen::Matrix3d Cne;
	posxyz = blh2xyz(posllh_rad, Cne);
	for (int i = 0; i < 3; i++)
		e_rb[i] = posxyz(i) - TC_conf.basepose[i];
	vector<double > enu(3, 0);
	ecef2enu(posbasellhxyz, e_rb, enu);
	vector<double > e_data,n_data;
	e_data.push_back(enu[0]);
	n_data.push_back(enu[1]);
	plt::scatter(e_data, n_data);
	plt::pause(0.0001);
	plt::draw();


}