#ifndef OBS_SAT_H_
#define OBS_SAT_H_

#include "Ephemeris.h"
#include <algorithm>

#include "earth_utility.h"


class solution {
public:
	solution() {}
	solution(int sysnum) {
		pos = vector<double>(3, 0); 
		posllh = vector<double>(3, 0);
		vel = vector<double>(3, 0);
		dtr = vector<double >(sysnum, 0);
		qr = vector<double>(6, 0);
		qrenu = vector<double>(6, 0);
		time = gtime_t(0, 0);
		state = 0;
		ratio = 0;
		prev_ratio1 = 0;
		prev_ratio2 = 0;
		thres = 0;
	}
	int get_state() { return state; }
	void set_state(int a) { state = a; }
	vector<double > getpos() { return pos; }
	vector<double > getposllh() { return posllh; }
	vector<double > getvel() { return vel; }
	vector<double > getdtr() { return dtr; }
	vector<double > getqr() { return qr; }
	vector<double > getqrenu() { return qrenu; }
	gtime_t gettime() { return time; }
	void setpos(double x, double y, double z) { pos[0] = x; pos[1] = y; pos[2] = z; }
	void setvel(double x, double y, double z) { vel[0] = x; vel[1] = y; vel[2] = z; }
	void setdtr(const vector<double> &dtr_c) { dtr = dtr_c; }
	void setqr(const vector<double> & qr_c) { qr = qr_c; }
	void setqv(const vector<double> & qr_c) { qv = qr_c; }
	void PQxyz2PQllh();
	void settime(gtime_t t) { time = t; }
	void setdrift(double a) { drift = a; }
	/*获得与求解模糊度有关的解*/
	float get_ratio() { return ratio; }
	float get_prev_ratio1() { return prev_ratio1; }
	float get_prev_ratio2() { return prev_ratio2; }
	float get_thres() {return thres; }
	void set_ratio(float a) { ratio = a; }
	void set_prev_ratio1(float a) { prev_ratio1 = a; }
	void set_prev_ratio2(float a) { prev_ratio2 = a; }
	void set_thres(float a) { thres = a; }
private:
	gtime_t time;
	int state;
	vector<double > pos;
	vector<double > posllh;
	vector<double > vel;
	vector<double > dtr;
	double drift;//clock drift
	vector<double > qr;//6维，xx yy zz xy yz zx
	vector<double > qv;//6维，xx yy zz xy yz zx
	vector<double > qrenu;//6维，ee nn uu en nu eu 对应llh和enu的P是一样的。
	//vector<double > q;//6维，xx yy zz xy yz zx
	/*与载波相位模糊度固定有关的变量*/
	float ratio;        /* AR ratio factor for valiation */
	float prev_ratio1;   /* previous initial AR ratio factor for validation */
	float prev_ratio2;   /* previous final AR ratio factor for validation */
	float thres;        /* AR ratio threshold for valiation */
};
class Station_OBS {
public:
	Station_OBS() {};
	Station_OBS(const int GNSS_Fnum, const int satindex, const vector<double > & GNSS_measurements, const vector<int > & GNSS_flags, const vector<double > & GNSS_fre,gtime_t time);
	~Station_OBS() {};
	void displayOBS_station();
	void displaySAT_station();
	double Sget_lambda(int i) { return lambda[i]; }
	double Sget_P(int i) { return P[i]; }
	double Sget_L(int i) { return L[i]; }
	double Sget_D(int i) { return D[i]; }
	double Sget_S(int i) { return S[i]; }
	int Sget_qualL(int i) { return qualL[i]; }
	int Sget_qualP(int i) { return qualP[i]; }
	int Sget_LLI(int i) { return LLI[i]; }
	void Sset_LLI(int Fi, int value);
	vector<double> get_razellos() { return r_azel_los; }
	double get_P_zres() { return P_zres; }
	double get_P_zvar() { return P_zvar; }
	double get_D_zres() { return D_zres; }
	double get_L_zres() { return L_zres; }
	double Sget_cp_t_tu_res() { return cp_t_tu_res; }
	int get_unusedflag() { return unused_flag; }
	int get_fdemask() { return fde_mask; }
	void set_fdemask(int a) { fde_mask = a; }
	void gettime(int *t, double *s) { *t = time.time; *s = time.sec; }
	gtime_t gettime() { return time; }
	gtime_t getlasttime() { return last_time; }
	void setlasttime(gtime_t t) { last_time = t; }
	void Ssetpvbd(const vector<double > &pvb, const vector<double > & bdb, double varb);
	void setunusedflag(int a) { unused_flag = a; };
	void set_cp_t_tu_res(double a) { cp_t_tu_res = a; }
	/*与单点定位有关的函数*/
	double cal_P(double *P, string mode, const vector<int> &fre,char constellation,double TGD);//mode 反应 单频，还是组合，还是L1C等补偿，这些以后扩展，目前就是单频
	void update_razellos(const vector<double > &sysclk_bias, const vector<double> &receiver_pos_ecef,
						 const vector<double> &receiver_pos_llh, int &count, const vector<double > &ION_CORR,
						 double sat_svh,char constellation,double TGD);
	void geodist(const vector<double> &receiver_pos_ecef);
	void satazel(const vector<double> &receiver_pos_llh);
	double ion_BRDC(const vector<double> &receiver_pos_llh, const vector<double> & ION_CORR);
	double tropmodel(const vector<double> &pos,const vector<double > &azel,double para);

	/*与最小二乘定速有关的函数*/
	void updateDres(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x, 
		int &count);
	/*void updateDres_Doppler(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x,
		int &count);*/
	//void updateDres_TDPhase();
	void Supdate_cp_t_tu_res(const vector<double> &posxyz,int SVH);
	double Scalculate_TDCPres(const vector<double > &last_sol_rr);
	/*与RTK定位有关的函数*/
	int zdres(const vector<double > & posecef, const vector<double > & posllh);//返回0 说明计算零差不成功，可能缺伪距或相位 返回1 零差计算成功
private:
	gtime_t time;
	gtime_t last_time;//为了TD载波相位
	int unused_flag;// 不使用标志位，比如卫星俯仰角过低，此位就置位，在构建H矩阵时就不再考虑此卫星
	int Ls_vel_flag;
	int fde_mask;//raim算法中排除某一个星的标志
	/*因为base rover 观测量的时间不完全一致，把卫星位置信息放到此处*/
	vector<double> SAT_POS_VEL;
	vector<double> SATCLK_BIAS_DRIFT;
	double SATPOS_VAR;
	/*从观测量和卫星信息计算得来的量*/
	double P_zres;//零差
	double L_zres;
	double D_zres;//用于LS定速存多普勒信息，版本2时增加了TD载波相位定速功能
	/*为了利用载波相位定速设置的一个量(cp_m-r_m+clight*delta_tu)*/
	double cp_t_tu_res;
	double P_zvar;
	vector<double > r_azel_los;//base 和 rover 都应该有不同的俯仰角、距离 和 LOS
	/*观测量本身的值*/
	vector<double > lambda;
	vector<double > P;
	vector<double > L;
	vector<double > D;
	vector<double > S;
	vector<int > qualL;
	vector<int > qualP;
	vector<int > LLI;
};
// 存放某一时刻两个站对某星的观测量 以及该时刻该星的相关信息
class OBS_SAT {
public:
	OBS_SAT() { BR_unused_flag = 0; BR_fde_mask = 0; fix = 0; base_flag = 0; rover_flag = 0; };
	OBS_SAT(int indexB, int indexR, const OBS_RINEX &obsrinex, char mode);
	double get_baseP(int i) { return base.Sget_P(i); };
	double get_roverP(int i) { return rover.Sget_P(i); };
	double get_rover_cp_t_tu_res() { return rover.Sget_cp_t_tu_res(); }// base 没必要算这个值
	//void get_P(int i, double * P_B, double *P_R) { *P_B = base.Sget_P(i); *P_R = rover.Sget_P(i); };//直接返回原始伪距，i代表频点 0 1 2 对应 1 2 5
	void get_lambda(int i, double * l_B, double *l_R) { *l_B = base.Sget_lambda(i); *l_R = rover.Sget_lambda(i); };
	void calP_sres();
	void calL_sres();
	void calD_sres();
	double getP_sres() { return P_sres; }
	double getP_svar() { return P_svar; }
	double getL_sres() { return L_sres; }
	double getL_svar() { return L_svar; }
	double getD_sres() { return D_sres; }
	double getD_svar() { return D_svar; }
	int get_baseflag() { return base_flag; }
	int get_roverflag() { return rover_flag; }
	int get_fix() { return fix; }
	void set_fix(int a) { fix = a; }
	int get_BR_unusedflag() { return BR_unused_flag; }
	int get_BR_fdemask() { return BR_fde_mask; }
	void set_BR_fdemask(int a) { BR_fde_mask = a; }
	void set_BR_unusedflag(int a) { BR_unused_flag = a; }
	void setroverunusedflag(int a) { rover.setunusedflag(a);}
	void set_rover_cp_t_tu_res(double a) { rover.set_cp_t_tu_res(a); }
	void set_rover_last_time(gtime_t t) { rover.setlasttime(t); }
	void update_cp_t_tu_res(const vector<double> pos) { rover.Supdate_cp_t_tu_res(pos,SAT_SVH); }
	void orBRLLI(int b, int r);
	Station_OBS get_rover() { return rover; }
	Station_OBS get_base() { return base; }

	void detslp_ll(char * slip,char * half,double tt);
	void setpvbdst(const vector<double > &pvb, const vector<double > & bdb, double varb, 
		         const vector<double > &pvr, const vector<double > & bdr, double varr, 
				 int svh, double TGD,int bflag,int rflag);
	void displayOBS_SAT();
	void displayOBS_SAT_Pos(int mode);
	void getBRtime(int *t, double *s) { *t = BRtime.time; *s =BRtime.sec; }
	gtime_t getBRtime() { return rover.gettime(); }
	///*与单点定位有关的函数*/
	void Ls_update_razellos(const vector<double > &sysclk_bias, const vector<double> &receiver_pos_ecef,
		const vector<double> &receiver_pos_llh, int &count, const vector<double > &ION_CORR);
	/*与最小二乘定速有关的函数*/
	void Ls_updateDres(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x, int &count);
	//void OBS_update_cp_t_tu_res() { rover.Supdate_cp_t_tu_res(); }//只更新rover的，base的我们暂时用不到
	double calculate_TDCPres(const vector<double > &last_sol_rr);
	/*与RTK定位有关的函数*/
	int BRzdres(char mode,const vector<double > &posecef,  const vector<double > & posllh);// B: base R: rover
	void calsdresvar(double bl);
private:
	gtime_t BRtime; 
	char constellation;
	/*卫星的位置、速度、钟差等，只需要时间就可以确定*/
	/*理论上，BR应该处理的是同一时刻的卫星数据，但可惜的是我们后处理的时候获得的时间标往往会有差别
	RTKLIB是分别在各自的时间标计算卫星位置，然后参与解算，但是这实际上是有问题的，因为计算的实际上是不同时刻的零差
	不同时刻的观测量在一起运作了，为了适应实际情况，我们把和位置钟差有关的卫星信息放到了base rover两个结构体中*/
	/*vector<double> SAT_POS_VEL;
	vector<double> SATCLK_BIAS_DRIFT;
	double SATPOS_VAR;*/
	double SAT_SVH;
	double TGD;
	/*卫星是否可用标志，要求BR对本星都可用*/
	int BR_unused_flag;
	int BR_fde_mask;
	
	/*与具体解算有关，当前位置距离卫星的距离r，航向角俯仰角azel，r到s的单位观测向量*/
	double P_sres, L_sres, D_sres;//对应每颗星的单差，可能是相位单差，伪距单差 和 多普勒单差// 注意这个单差只是观测量的单差，还没有和状态量做差和状态量做差仍然需要zrdes
	double P_svar, L_svar, D_svar;//对应ddres 的 Ri Rj
	/*与模糊度解算有关的变量*/
	int fix;
	Station_OBS base;
	Station_OBS rover;
	int base_flag;
	int rover_flag;
	
};
#endif