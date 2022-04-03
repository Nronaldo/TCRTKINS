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
	void udstate_rtkins();//Ϊ������rtk�ĺ�������ins��λ���ٶ�+lever����RTK_x,Ȼ������receiver��zdres��ddres����в
	vector<int > rtkins_FormHxP(const Eigen::MatrixXd &H_gnss, const Eigen::VectorXd &x_gnss, const Eigen::MatrixXd &P_gnss,
		Eigen::MatrixXd &H_TC,  Eigen::VectorXd &x_TC,  Eigen::MatrixXd &P_TC, const vector<int > &ix);
	void rtkins_Measureupdate();//�����vHx�����IGG���������
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
	/*�����״̬�й�1 �̶� 2 ���� 3 ���ߵ�*/
	void set_Solstate(int x) { Sol_state = x; }
	int get_Solstate() { return Sol_state; }
	/*��IGGʹ���й�*/
	void set_ngnsslock(int a) { ngnsslock = a; }
	void add_ngnsslock() { ngnsslock++; }
	void Robust_IGG(Eigen::MatrixXd & H_TC, Eigen::VectorXd &v_GNSS, Eigen::MatrixXd &R_GNSS,const Eigen::MatrixXd &P);
	/*�뻭ͼ�й�*/
	void Fusion_plot();
private:
	receiver rcv;
	CSINS sins;
	int ins_align_flag;
	int TC_nx;
	//int ins_realign_flag;
	double fusion_time;
	double *INS_Sbias_P;//��һ�����ɱ��������ins �Լ��������ǵ���״̬��Э���
	//rcv �� RTK_PΪ�˱��ֺ�RTKLIBһ�£��ǰ����д�ġ���������Ҳ���д�
	Eigen::VectorXd INS_Error_x;
	int Sol_state;
	int ngnsslock;
};


#endif