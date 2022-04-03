#ifndef _FUSION_H_
#define _FUSION_H_
#include "Receiver.h"
#include "sins_TC.h"

class Fusion {
	/*friend void Fusionmatrixtofile(const Fusion & class_fusion);*/
public:
	Fusion() { ins_align_flag = 0;  TC_nx = 15 + NSATGPS + NSATCMP; INS_Sbias_P = new double[TC_nx*TC_nx] { 0 };
	INS_Error_x = Eigen::VectorXd::Zero(15); Sol_state = 0; ngnsslock = 0;}
	~Fusion() {  free(INS_Sbias_P); };
	/*GNSS related function*/
	//const receiver get_receiver() { return rcv; }
	receiver &get_receiver() { return rcv; }
	void initialze_receiver(NAV_RINEX &nav,ofstream & out) { rcv.initializercv(nav,out); }
	void Clear_rcv_array() { rcv.Set_array_zero(); }
	/*Fuction to set receiver*/

	/*INS related function*/
	void set_sins(const vector<double > &dwdv, double t) { set_sins_time(t); set_sins_IMU(dwdv); }
	void set_sins_time(double t) { sins.set_time(t); }
	void set_sins_IMU(const vector<double > &dwdv) { sins.set_IMU(dwdv); }
	void ins_align(ofstream & out);
	//CSINS get_sins() { return sins; }
	int get_ins_align_flag() { return ins_align_flag; }
	void ins_mech_ErrorTimeUpdate();
	/*Fusion related function*/
	void set_fusion_time(double t) { fusion_time = t; }
	double get_fusion_time() { return fusion_time; }
	int rtk_ins();
	void udstate_rtkins();//为了利用rtk的函数，用ins的位置速度+lever设置RTK_x,然后利用receiver的zdres和ddres计算残差，
	vector<int > rtkins_FormHxP(const Eigen::MatrixXd &H_gnss, const Eigen::VectorXd &x_gnss, const Eigen::MatrixXd &P_gnss,
		Eigen::MatrixXd &H_TC,  Eigen::VectorXd &x_TC,  Eigen::MatrixXd &P_TC, const vector<int > &ix);
	void rtkins_Measureupdate();//再组合vHx，兼顾IGG和周跳检测
	void Fusion_get_RTK_INS_P();
	void rtkins_filter_unzip(Eigen::MatrixXd &H, Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::VectorXd &v, Eigen::MatrixXd &R, const vector<int > &ix);
	/*AR related function*/
	int rtkins_resamb(double *RTK_fixX, vector<int > &Sbias_ix);
	void rtkins_holdamb(double *RTK_xfix, const vector<int > &Sbias_ix, int nb, int firstvalidGPSID, int firstvalidBDSID);
	void rtkins_holdamb_Filter_FB(const Eigen::MatrixXd &H, const Eigen::VectorXd & v, const Eigen::MatrixXd & P, const Eigen::MatrixXd &R,
		const vector<int > 	&allstateindex);
	/*save to file*/
	void TC_instofile(double soltime,ofstream& out);
	void TC_instofile_time(double soltime, ofstream& out);
	void TC_instofile_pos(ofstream& out);
	void TC_instofile_posvar(ofstream& out);
	void TC_instofile_vel(ofstream& out);
	void TC_instofile_velvar(ofstream& out);
	void TC_instofile_att(ofstream& out);
	void TC_instofile_attvar(ofstream& out);
	void TC_instofile_bg(ofstream& out);
	void TC_instofile_ba(ofstream& out);
	void TC_instofile_head(ofstream& out);
	/*与计算状态有关1 固定 2 浮点 3 纯惯导*/
	void set_Solstate(int x) { Sol_state = x; }
	int get_Solstate() { return Sol_state; }
	/*与IGG使用有关*/
	void set_ngnsslock(int a) { ngnsslock = a; }
	void add_ngnsslock() { ngnsslock++; }
	void Robust_IGG(Eigen::MatrixXd & H_TC, Eigen::VectorXd &v_GNSS, Eigen::MatrixXd &R_GNSS,const Eigen::MatrixXd &P);
	/*与画图有关*/
	void Fusion_plot();
private:
	receiver rcv;
	CSINS sins;
	int ins_align_flag;
	int TC_nx;
	//int ins_realign_flag;
	double fusion_time;
	double *INS_Sbias_P;//是一个过渡变量，存放ins 以及各个卫星单差状态的协方差。
	//rcv 的 RTK_P为了保持和RTKLIB一致，是按照列存的。我们这里也按列存
	Eigen::VectorXd INS_Error_x;
	int Sol_state;
	int ngnsslock;
};


#endif