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
	/*��������ģ�����йصĽ�*/
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
	vector<double > qr;//6ά��xx yy zz xy yz zx
	vector<double > qv;//6ά��xx yy zz xy yz zx
	vector<double > qrenu;//6ά��ee nn uu en nu eu ��Ӧllh��enu��P��һ���ġ�
	//vector<double > q;//6ά��xx yy zz xy yz zx
	/*���ز���λģ���ȹ̶��йصı���*/
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
	/*�뵥�㶨λ�йصĺ���*/
	double cal_P(double *P, string mode, const vector<int> &fre,char constellation,double TGD);//mode ��Ӧ ��Ƶ��������ϣ�����L1C�Ȳ�������Щ�Ժ���չ��Ŀǰ���ǵ�Ƶ
	void update_razellos(const vector<double > &sysclk_bias, const vector<double> &receiver_pos_ecef,
						 const vector<double> &receiver_pos_llh, int &count, const vector<double > &ION_CORR,
						 double sat_svh,char constellation,double TGD);
	void geodist(const vector<double> &receiver_pos_ecef);
	void satazel(const vector<double> &receiver_pos_llh);
	double ion_BRDC(const vector<double> &receiver_pos_llh, const vector<double> & ION_CORR);
	double tropmodel(const vector<double> &pos,const vector<double > &azel,double para);

	/*����С���˶����йصĺ���*/
	void updateDres(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x, 
		int &count);
	/*void updateDres_Doppler(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x,
		int &count);*/
	//void updateDres_TDPhase();
	void Supdate_cp_t_tu_res(const vector<double> &posxyz,int SVH);
	double Scalculate_TDCPres(const vector<double > &last_sol_rr);
	/*��RTK��λ�йصĺ���*/
	int zdres(const vector<double > & posecef, const vector<double > & posllh);//����0 ˵���������ɹ�������ȱα�����λ ����1 ������ɹ�
private:
	gtime_t time;
	gtime_t last_time;//Ϊ��TD�ز���λ
	int unused_flag;// ��ʹ�ñ�־λ���������Ǹ����ǹ��ͣ���λ����λ���ڹ���H����ʱ�Ͳ��ٿ��Ǵ�����
	int Ls_vel_flag;
	int fde_mask;//raim�㷨���ų�ĳһ���ǵı�־
	/*��Ϊbase rover �۲�����ʱ�䲻��ȫһ�£�������λ����Ϣ�ŵ��˴�*/
	vector<double> SAT_POS_VEL;
	vector<double> SATCLK_BIAS_DRIFT;
	double SATPOS_VAR;
	/*�ӹ۲�����������Ϣ�����������*/
	double P_zres;//���
	double L_zres;
	double D_zres;//����LS���ٴ��������Ϣ���汾2ʱ������TD�ز���λ���ٹ���
	/*Ϊ�������ز���λ�������õ�һ����(cp_m-r_m+clight*delta_tu)*/
	double cp_t_tu_res;
	double P_zvar;
	vector<double > r_azel_los;//base �� rover ��Ӧ���в�ͬ�ĸ����ǡ����� �� LOS
	/*�۲��������ֵ*/
	vector<double > lambda;
	vector<double > P;
	vector<double > L;
	vector<double > D;
	vector<double > S;
	vector<int > qualL;
	vector<int > qualP;
	vector<int > LLI;
};
// ���ĳһʱ������վ��ĳ�ǵĹ۲��� �Լ���ʱ�̸��ǵ������Ϣ
class OBS_SAT {
public:
	OBS_SAT() { BR_unused_flag = 0; BR_fde_mask = 0; fix = 0; base_flag = 0; rover_flag = 0; };
	OBS_SAT(int indexB, int indexR, const OBS_RINEX &obsrinex, char mode);
	double get_baseP(int i) { return base.Sget_P(i); };
	double get_roverP(int i) { return rover.Sget_P(i); };
	double get_rover_cp_t_tu_res() { return rover.Sget_cp_t_tu_res(); }// base û��Ҫ�����ֵ
	//void get_P(int i, double * P_B, double *P_R) { *P_B = base.Sget_P(i); *P_R = rover.Sget_P(i); };//ֱ�ӷ���ԭʼα�࣬i����Ƶ�� 0 1 2 ��Ӧ 1 2 5
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
	///*�뵥�㶨λ�йصĺ���*/
	void Ls_update_razellos(const vector<double > &sysclk_bias, const vector<double> &receiver_pos_ecef,
		const vector<double> &receiver_pos_llh, int &count, const vector<double > &ION_CORR);
	/*����С���˶����йصĺ���*/
	void Ls_updateDres(const vector<double> &receiver_pos_ecef, const Eigen::Matrix3d &E, const Eigen::VectorXd &Ls_x, int &count);
	//void OBS_update_cp_t_tu_res() { rover.Supdate_cp_t_tu_res(); }//ֻ����rover�ģ�base��������ʱ�ò���
	double calculate_TDCPres(const vector<double > &last_sol_rr);
	/*��RTK��λ�йصĺ���*/
	int BRzdres(char mode,const vector<double > &posecef,  const vector<double > & posllh);// B: base R: rover
	void calsdresvar(double bl);
private:
	gtime_t BRtime; 
	char constellation;
	/*���ǵ�λ�á��ٶȡ��Ӳ�ȣ�ֻ��Ҫʱ��Ϳ���ȷ��*/
	/*�����ϣ�BRӦ�ô������ͬһʱ�̵��������ݣ�����ϧ�������Ǻ����ʱ���õ�ʱ����������в��
	RTKLIB�Ƿֱ��ڸ��Ե�ʱ����������λ�ã�Ȼ�������㣬������ʵ������������ģ���Ϊ�����ʵ�����ǲ�ͬʱ�̵����
	��ͬʱ�̵Ĺ۲�����һ�������ˣ�Ϊ����Ӧʵ����������ǰѺ�λ���Ӳ��йص�������Ϣ�ŵ���base rover�����ṹ����*/
	/*vector<double> SAT_POS_VEL;
	vector<double> SATCLK_BIAS_DRIFT;
	double SATPOS_VAR;*/
	double SAT_SVH;
	double TGD;
	/*�����Ƿ���ñ�־��Ҫ��BR�Ա��Ƕ�����*/
	int BR_unused_flag;
	int BR_fde_mask;
	
	/*���������йأ���ǰλ�þ������ǵľ���r������Ǹ�����azel��r��s�ĵ�λ�۲�����*/
	double P_sres, L_sres, D_sres;//��Ӧÿ���ǵĵ����������λ���α�൥�� �� �����յ���// ע���������ֻ�ǹ۲����ĵ����û�к�״̬�������״̬��������Ȼ��Ҫzrdes
	double P_svar, L_svar, D_svar;//��Ӧddres �� Ri Rj
	/*��ģ���Ƚ����йصı���*/
	int fix;
	Station_OBS base;
	Station_OBS rover;
	int base_flag;
	int rover_flag;
	
};
#endif