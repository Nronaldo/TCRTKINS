# ifndef _SINS_TC_H_
# define _SINS_TC_H_
#include "Load_Configuration.h"
#include "Time_utility.h"
#include "Math_utility.h"
using namespace std;
// class definition
extern Eigen::Matrix3d I33, O33;
class CEarth;
class CEarth{
public:
	double RN,RM,Mpv2,Mpv4,RNh,RMh,sinL,tanL,secL,cosL,sin2L,cos2L,g;
	Eigen::Vector3d wnie, wnen, wnin,wnien, gn, gcc;

	CEarth();
	void Update(const Eigen::Vector3d &pos, const Eigen::Vector3d &vn = Eigen::Vector3d(0, 0, 0));
	double get_g0() { return g0; }
	/*Test*/
	void Load_CEarth(const string & filename, Eigen::Vector3d&pos, Eigen::Vector3d&vel);
private:
	double Re, Rp;
	double f, e1, e2;
	double wie, g0;
};
class CIMU
{
public:
	CIMU() { };
	void update(const vector<double > &imu);
	//gtime_t get_time() const { return time; }
	Eigen::Vector3d get_dw() const { return dw; }
	Eigen::Vector3d get_dv() const { return dv; }
	//CIMU & operator=(const CIMU & imu_o);
private:
	//gtime_t time;
	Eigen::Vector3d  dw, dv;
};
class CSINS
{
public:
	
	//Eigen::VectorXd old_imud;
	//Eigen::VectorXd x,xa;
	
	
	//Eigen::MatrixXd Pa;
	// Eigen::Vector3d oldpos;
	//struct oldobsr,oldobsb;
	
	
	//Eigen::Matrix3d Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp;	// for etm
	
	/*INS更新函数*/
	CSINS() {};
	void ins_init(const Eigen::Vector3d &att0, const Eigen::Vector3d &vn0, const Eigen::Vector3d &pos0);    // initialization using quat attitude, velocity & position
	//void SetTauGA(CVect3 &tauG, CVect3 &tauA);
	void Update();		// SINS update using Gyro&Acc samples
	//void Extrap(const CVect3 &wm, const CVect3 &vm, double ts);						// SINS fast extrapolation using 1 Gyro&Acc sample
	void lever_trans(int mode);		// lever arm 0: GNSS->INS 1 INS->GNSS
	void lever_trans(int mode,Eigen::Vector3d &pos, Eigen::Vector3d &vel);
	Eigen::Vector3d get_Cnb_lever() { return Cnb * lever; }
	Eigen::MatrixXd get_P() { return P; }
	Eigen::MatrixXd get_Pa() { return Pa; }
	Eigen::Vector3d get_sins_pos() { return pos; }
	Eigen::Vector3d get_sins_xapos() { return xa.segment(6,3); }
	Eigen::Vector3d get_sins_xavel() { return xa.segment(3, 3); }
	Eigen::Vector3d get_sins_att() { return att; }
	Eigen::Vector3d get_sins_xaatt() { return xa.segment(0, 3); }
	Eigen::Vector3d get_sins_vel() { return vel; }
	Eigen::Vector3d get_sins_bg() { return bg; }
	Eigen::Vector3d get_sins_xabg() { return xa.segment(9, 3);}
	Eigen::Vector3d get_sins_xaba() { return xa.segment(12, 3); }
	Eigen::Vector3d get_sins_ba() { return ba; }
	void set_time(double t) { time = t; }
	double get_time() { return time; }
	void set_IMU(const vector<double> & dwdv) { imu.update(dwdv); }
	/*ins related error states update function*/
	void update_trans_mat();							// SINS error transform matrix coefficients
	void ins_Error_time_update();// SINS error time update (15*1 state) only update P because delta_x will be set as zero after correction.
	void feedback(const Eigen::VectorXd & x,const double * State_P,int State_P_dimension);//浮点解，参与kalman
	void feedback(const Eigen::VectorXd & x, const Eigen::MatrixXd &State_P);//重载，固定解，只做本历元输出，不参与下一历元kalman
	/*Test*/
	void Load_CSINS();
private:
	double nt;
	CEarth eth;
	double time;
	Eigen::MatrixXd P, Q,Pa;
	Eigen::VectorXd xa;
	Eigen::MatrixXd Phi;
	CIMU imu;
	Eigen::Vector3d att, vel, pos, acc, bg, ba, tauG, tauA, wib, web, fb, fn, lever;
	Eigen::Matrix3d Cnb, Kg, Ka;
	Eigen::Matrix3d Mpv;	// for etm
	
};
Eigen::Matrix3d att2Cnb(const Eigen::Vector3d &att);
Eigen::Vector3d Cnb2att(const Eigen::Matrix3d &Cnb);
Eigen::Matrix3d rvec2mat(const Eigen::Vector3d &rot_vector);
Eigen::Matrix3d askew(const Eigen::Vector3d &v);
/*MATH Function*/
int		sign(double val, double eps = EPS);
double	range(double val, double minVal, double maxVal);
double	atan2Ex(double y, double x);
double asinEx(double x);
double acosEx(double x);
/*Test function*/
void Test_CSINS();
void Test_CEarth();
void Load_avp(const string& filenameavp, Eigen::Vector3d &att, Eigen::Vector3d &vel, Eigen::Vector3d &pos);
#endif
