#ifndef RECEIVER_H_
#define RECEIVER_H_

#include "Ephemeris.h"
#include <algorithm>
#include "OBS_SAT.h"
#include "Math_utility.h"
#include "Load_Configuration.h"
namespace plt = matplotlibcpp;
class BRslip {
public:
	BRslip() { Bslip = 0; Rslip = 0; };
	BRslip(int b, int r) { Bslip = b; Rslip = r; };
	int getBslip() { return Bslip; };
	int getRslip() { return Rslip; };
	void setBslip(int a) { Bslip = a; };
	void setRslip(int a) { Rslip = a; };
private:
	int Bslip;
	int Rslip;
};

class receiver {
	/*友元,为了充分利用receiver现有的函数，由Fusion传入参数*/
	friend void Fusion_set_RTKx_posvel(receiver & rcv,const Eigen::Vector3d & Ext_pos, const Eigen::Vector3d & Ext_vel);
	friend void Fusion_set_RTKxa_posvel(receiver & rcv, const Eigen::Vector3d & Ext_pos, const Eigen::Vector3d & Ext_vel);
	friend void Fusion_set_RTKxP_bias(receiver & rcv, const Eigen::VectorXd &v, const Eigen::MatrixXd &P, const vector<int> &TC_ix);
public:
	receiver() {};
	//receiver(NAV_RINEX &nav_rinex,int dynamics_on,int arfilter,string outsolname);
	~receiver() {delete RTK_x; delete RTK_P; delete outc_G; delete outc_B; delete slip_G;
	delete slip_B; delete half_G; delete half_B; delete lock_G; delete lock_B; 
	//out.close();
	}
	void initializercv(NAV_RINEX &nav_rinex,ofstream & out);
	/*与读入观测量和更新卫星状态有关的函数*/
	void updateobs(const OBS_RINEX & obsrinex);
	void distributesats(const vector<int> &GNSS_ID_B, const vector<int> &GNSS_ID_R, char mode, const OBS_RINEX & obsrinex);//把观测量分配到两个结构中
	void choosecommonsat();//删除上个函数产生的 只有base 或者只有rover的
	void updateSatinfo();
	void Set_array_zero();
	/*与显示有关的函数*/
	void displayobs();
	void displayddres();
	void displaysatpos(int mode);
	void soltofile(int LS_stat,ofstream & out);
	void closefile(ofstream & out) { out.close(); };
	void plot_ENU_traj();
	/*与定位有关的函数*/
	int DOPcheck();//1 normal 0 gdop too big
	void calculate_DOP(const vector<double > &azel);
	int validsatnum();
	/*与最小二乘单点定位有关的函数*/
	int LS_pos(char mode);// s如果是 N 代表普通定位，可以更新sol，如果F 代表 FDE，此时不在LS_pos中更新sol，需要在外层的fde函数中判断一下是否更新
	//返回值：1 正常 0 卫星数不足导致无法FDE 或者定位，-1 卫星数够，但是发散了
	int rescode(const vector<double > &GPS_ION_CORR, const vector<double > &BDS_ION_CORR,int &countG,int &countB);
	
	int LS_pos_body(char mode,int countG,int countB);
	int LS_raim();// 1 normal 0 residual too big	
	void LS_fde();
	
	/*与最小二乘多普勒定速有关的函数*/
	int LS_vel();
	int resdop();
	int LS_vel_body();
	/*与TD载波相位定速有关的函数*/
	int tdcp2vel(Eigen::Vector3d &vn);
	int res_tdcp(const vector<double > & last_sol_pos, map<int, double>&Sat_TDCPresG, map<int, double>&Sat_TDCPresB);
	int TDcp_robust_exclude(map<int, double> &Sat_TDCPresG, map<int, double> &Sat_TDCPresB);
	int Form_Solve_TDcp(map<int, double> &Sat_TDCPresG, map<int, double> &Sat_TDCPresB, Eigen::Vector3d & vn);
	/*与RTK定位有关的函数*/
	int RTK_pos(string RTKmode);// mode：P L D 分别对应 DGPS RTK 加入多普勒的RTK。设置mode的目的是不重复计算。先从伪距加相位算起
	void RTK_tosol();
	void lock_solve();
	void out_solve();
	/*总定位函数*/
	void GNSS_sever(ofstream &out);
	/*1与计算零差有关*/
	int RCV_zdres(char mode, double *);
	/*2 与计算双差有关*/
	int RCV_ddres(string RTKmode, double *);
	void findrefsat();
	void calddres(string RTKmode, double *);
	vector<int> get_GBSat_PLdres_size();
	/*3与状态更新有关*/
	void initxP( double xvalue, double pvalue, int index);
	void rcv_detslp_ll();
	void reset_outageSAT();
	void initialize_bias_state();
	void udstate();
	void udpos();
	void udbias();// 返回值，是因为周跳需要重置的INS单差状态在RTK_X的索引，搭配INS_TC重置互相关使用，注意索引的转换
	const double * get_RTK_P() { return RTK_P; }
	int get_nx() { return nx; }
	/*4 与观测更新有关*/
	void MeasureUpdate(string RTKmode);
	void FormHxP(string RTKmode,Eigen::MatrixXd &H, Eigen::VectorXd &x,Eigen::MatrixXd &P, int ddpairnum_GPSP, int ddpairnum_GPSL, int ddpairnum_BDSP, int ddpairnum_BDSL, vector<int> &ix);
	void FormvR(string RTKmode, Eigen::VectorXd &y,Eigen::MatrixXd &R, int ddpairnum_GPSP, int ddpairnum_GPSL, int ddpairnum_BDSP, int ddpairnum_BDSL);
	void filter_unzip(Eigen::MatrixXd &H, Eigen::VectorXd &x, Eigen::MatrixXd &P, Eigen::VectorXd &v, Eigen::MatrixXd &R,const vector<int > &ix);
	int solve_validpos(string RTKmode,double *,int mode);// mode:0 只判断 为1  进行lock 和 out的更改
	int valpos(double t);
	/*5 与固定模糊度有关*/
	void set_holdamb_flag(int a) { holdamb_flag = a; }
	int getSATrefGPSID() {return SATrefGPSID;}
	int getSATrefBDSID() { return SATrefBDSID; }
	int manage_amb_LAMBDA(double *rtk_fixX,vector<int > &ix);//主要作用是进行lambda，除此之外设计了一些策略，比如说一颗星出现一定时间之后再进行模糊度固定，或者进行部分模糊度固定
	int exclude_different_sat(int &exc_satsID, int &exc_sats_lock);//每个历元排除一个不同的卫星（实际上策略过于简单，只是数组索引递增，并不能保证每次排除的都是不同的卫星）
	void recover_different_sat(int,int);//排除和恢复的方法都是对lock进行操作。但恢复也不会重新计算模糊度，
	int exclude_new_sat();//排除当前历元新参与计算模糊度的卫星。
	int resamb_LAMBDA(double * a, vector<int > &ix);
	void restamb(double *RTK_fixX, const vector<int> & ix,const Eigen::MatrixXd &b,int mode);//mode 0 : GNSS 1: TC
	void holdamb(double *RTK_xfix, const vector<int > &ddres_stateindex,int nb, int firstvalidGPSID, int firstvalidBDSID);
	int holdamb_Form_xPRHv(double *RTK_xfix, const vector<int > &ddres_stateindex, const int nb, int firstvalidGPSID, int firstvalidBDSID,
		vector<int > &allstateindex,
		Eigen::VectorXd &x, Eigen::MatrixXd &H, Eigen::VectorXd & v, Eigen::MatrixXd & P, Eigen::MatrixXd &R);
	void holdamb_Filter_FB(const Eigen::VectorXd &x, const Eigen::MatrixXd &H, const Eigen::VectorXd & v, const Eigen::MatrixXd & P, const Eigen::MatrixXd &R,
	const vector<int > 	&allstateindex);
	int ddmat_FormxP(Eigen::MatrixXd & D,Eigen::VectorXd &x,Eigen::MatrixXd & P,vector<int> &ix);
	vector<int> find_joinfixsat(int &nb1,int &nb2);
	void outddint(int info, const vector<int > & ix, const Eigen::MatrixXd & D, const Eigen::MatrixXd &b);//输出双差模糊度
	void set_sol_ratio(double a) { sol.set_ratio(a); }
	double get_sol_ratio() { return sol.get_ratio(); }
	void add_nfix() { nfix++; }
	int get_nfix() { return nfix; }
	void set_nfix(int a) { nfix = a; }
	vector<int > get_Reset_Sbias_index() { return Reset_Sbias_index; }
	/*void rcv_set_allsateph();*/
	/*6 与获得solution有关的接口*/
	gtime_t getsoltime() { return sol.gettime(); };
	void setsoltime(gtime_t t) { sol.settime(t); }
	void settt(gtime_t LASTsoltime);
	int get_sol_state() { return sol.get_state(); }
	vector<double> get_sol_pos(int mode) { return !mode?sol.getpos():sol.getposllh();}
	//solution get_sol() { return sol; }
	/*7 与统计非整s历元的LLI标志有关*/
	void saveslips(const OBS_RINEX & obs, int Btintflag, int Rtintflag);
	void restoreslips();
	/*8 与画图有关的*/
	short plot_mode;
	
private:
	/*与定位模式有关的变量*/
	//int Pos_stat;//0 无法定位 4 单点定位 1固定RTK 2 浮点RTK  
	/*与最小二乘单点定位有关的中间量*/
	Eigen::VectorXd Ls_x;
	Eigen::VectorXd v;
	Eigen::MatrixXd H;
	double Ls_res;
	/*与RTK定位有关的中间量*/
	Eigen::VectorXd RTK_fixXpva;
	Eigen::MatrixXd RTK_fixPpva;
	double tt;//当前时刻和上一时刻的时间差。若上一时刻还未开始定位，tt值应为0,RTKLIB是不管能不能定位都给了值
	// 我们仍然沿用这种做法
	//Eigen::MatrixXd RTK_Fpva;
	//Eigen::MatrixXd RTK_Q;
	int dynamics_on;
	double com_bias;
	int nx;
	double *RTK_x;
	double *RTK_P;//是一个过渡变量，存放理论上各个卫星状态的协方差。//实在是有点大。。 为了和RTKLIB一致好算我们按列存
	/*与固定模糊度有关的变量*/
	int nfix;
	short * outc_G;//GPS某星丢失时间
	short * outc_B;//BDS某星丢失时间
	short * lock_G;//GPS模糊度锁定历元
	short * lock_B;//BDS 模糊度锁定历元
	char * slip_G;//GPS 某星 周跳情况
	char * half_G;
	char * slip_B;//BDS 某星 周跳情况,主要原因在于RTKLIB 在LLI检测周跳的时候还要和上一时刻的周跳情况联合判断，所以必须得有一个全局数组保留，上一时刻的
	char * half_B;
	int excsatindex;
	int holdamb_flag;
	int arfilter;
	vector<int> Reset_Sbias_index; //发生周跳需重置的卫星单差状态在RTK_x中的位置
	// 双差量，按照RTKLIB的做法，两个星座各自组成双差
	int SATrefGPSID, SATrefBDSID;
	map<int, double > SatSat_PdresG;//key 值是另外一颗星
	map<int, double > SatSat_LdresG;
	map<int, double > SatSat_DdresG;
	map<int, double > SatSat_PdresB;//key 值是另外一颗星
	map<int, double > SatSat_LdresB;
	map<int, double > SatSat_DdresB;
	vector<double> DOP;// GDOP PDOp HDOP VDOP
	/*与定位有关的星历观测量*/
	All_Sats_EPH allsateph;
	map<int, OBS_SAT> GPSobs;
	map<int, OBS_SAT> BDSobs;
	map<int, OBS_SAT> NoB_R_GPSobs;//为了TD载波相位定速而引入的量，choosecomonsat 把非BR共存的obs擦除了
	map<int, OBS_SAT> NoB_R_BDSobs;//这里把无B有R的星保留下来，为了TD载波相位定速的时候用
	vector<double > basepos;//这里我们都使用提前标定好的位置
	/*与保存slip有关的量*/
	map<int, BRslip> GPSslip;//目前只考虑了单频
	map<int, BRslip> BDSslip;
	solution sol;
	//ofstream out;
	/*与输出模式有关的量*/
	string output_pos_format;
	
};
int screent_obs(ifstream &fileB, ifstream &fileR, string &buffB, string &buffR,receiver & rcv, OBS_RINEX &obs);
#endif
