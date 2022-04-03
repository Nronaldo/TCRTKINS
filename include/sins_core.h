#ifndef _SINS_CORE_H_
#define _SINS_CORE_H_
//从严老师的C++代码修改过来的
//后来感觉大多数用不到的，只要能完成更新就可以
//所以最终没用这个文件，基于Eigen写的sins_TC.h
#include "Load_Configuration.h"
using namespace std;
// class definition
class CGLV;
class CVect3; class CMat3; class CQuat;
class CEarth; class CIMU; class CSINS;

extern const CVect3 I31, O31;
extern const CQuat	qI;
extern const CMat3	I33, O33;
extern const CGLV	glv;
int		sign(double val, double eps = EPS);
double	range(double val, double minVal, double maxVal);
double	atan2Ex(double y, double x);
double asinEx(double x);
double acosEx(double x);
class CGLV
{
public:
	double Re, f, g0, wie;											// the Earth's parameters
	double e, e2;
	double mg, ug, deg, min, sec, hur, ppm, ppmpsh;					// commonly used units
	double dps, dph, dpsh, dphpsh, ugpsh, ugpsHz, mpsh, mpspsh, secpsh;
	CGLV();
};
class CEarth
{
public:
	double a, b;
	double f, e, e2;
	double wie;

	double sl, sl2, sl4, cl, tl, RMh, RNh, clRNh, f_RMh, f_RNh, f_clRNh;
	CVect3 pos, vn, wnie, wnen, wnin, gn, gcc;

	CEarth(double a0 = glv.Re, double f0 = glv.f, double g0 = glv.g0);
	void Update(const CVect3 &pos, const CVect3 &vn = CVect3(0,0,0));
	CVect3 vn2dpos(const CVect3 &vn, double ts = 1.0) const;
};
class CVect3 {
	friend ostream & operator<<(ostream & os, const CVect3 &cv3);
	friend istream & operator>>(istream & is,CVect3 &cv3);
	friend CVect3 operator +(const CVect3 &lhs, const CVect3 &rhs);
	friend CVect3 operator- (const CVect3 & lhs,const CVect3 &rhs);//subtract
	friend CVect3 operator*(const CVect3 & lhs, const CVect3 &rhs);// cross 
	friend CVect3 operator*(const CVect3 & lhs, double f);//*double
	friend CVect3 operator*(double f,const CVect3 & lhs );//*double
	friend CVect3 operator/(const CVect3 & lhs, double f);//*double
	friend CVect3 operator-(const CVect3 &cv3);// minus
	friend double norm(const CVect3 &v);// 
	friend double normXY(const CVect3 &v);
	friend CVect3 sqrt(const CVect3 &v);
	friend double dot(const CVect3 & v1,const CVect3 & v2);
public:	
	CVect3() {};
	~CVect3() {};
	CVect3(double xx, double yy, double zz) { elements[0] = xx; elements[1] = yy; elements[2] = zz;}
	double &operator[](int n);//可修改
	const double & operator[](int n) const;// 返回值不可修改，且不可修改任何成员
	CVect3 & operator =(const CVect3 & cv3);
	CVect3 & operator +=(const CVect3 & cv3);
	CVect3 & operator -=(const CVect3 & cv3);
	CVect3 & operator *=(double f);
	CVect3 & operator /=(double f);
	bool IsZero() const;
	bool IsZeroXY() const;
private:
	double elements[3];
};
class CQuat {
	friend ostream & operator<< (ostream &os,CQuat q);
	friend istream & operator>> (istream &is, CQuat q);
	friend CVect3 operator-(const CQuat & nominalq, const CQuat & trueq);// 四元数与真实四元数作差，转为姿态误差
	friend CQuat operator-(const CQuat & nominalq, const CVect3 & v_phi);//四元数中校正姿态
	friend CQuat operator*(const CQuat & q1, const CQuat & q2);//四元数乘法
	friend CVect3 operator*(const CQuat & q1, const CVect3 & v2);//矢量经过四元数旋转 q*v*(~q)
	friend CQuat operator~(const CQuat &q);//共轭
	
public:
	CQuat() {};
	CQuat(double w,double x,double y,double z) { elements[0] = w; elements[1] = x; elements[2] = y; elements[3] = z;}
	double &operator[](int n);
	const double &operator[](int n) const;
	CQuat & operator=(const CQuat & phi);
	CQuat & operator-=(const CVect3 & v_phi);
	CQuat & operator*=(const CQuat & q);
	void normlize();
private:
	double elements[4];


};
class CMat3 {
	friend ostream & operator<<(ostream &os,CMat3 & m);
	friend istream & operator<<(istream &os, CMat3 & m);
	friend CMat3 operator+(const CMat3 &m1,const CMat3 & m2);
	friend CMat3 operator-(const CMat3 &m1, const CMat3 & m2);
	friend CMat3 operator-(const CMat3 &m);
	friend CMat3 operator*(const CMat3 &m1, const CMat3 & m2);
	friend CMat3 operator*(const CMat3 &m,double f);
	friend CMat3 operator*( double fconst, CMat3 &m);
	friend CVect3 operator*(const CMat3 &m, const CVect3 & v);
	friend CMat3 operator~(const CMat3 &m);					// matirx transposition
	friend double det(const CMat3 &m);
	friend CMat3 inv(const CMat3 &m);						// matirx inverse
	friend CVect3 diag(const CMat3 &m);
	friend CMat3 diag(const CVect3 &v);
public:
	CMat3() {};
	CMat3(double xx, double xy, double xz,
		double yx, double yy, double yz,
		double zx, double zy, double zz);
	CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2);  // M = [v0; v1; v2]);0
	double &operator()(int r,int c);
	const double &operator()(int r, int c)const;
	CMat3&operator=(const CMat3 & m);
	CMat3&operator+=(const CMat3 & m);
	CMat3&operator-=(const CMat3 & m);
private:
	double elements[9];


};
class CIMU
{
public:
	int prefirst;
	CVect3 phim, dvbm, wm_1, vm_1;
	CIMU() { prefirst = 1; }
	void Update(const vector<CVect3> &wm, const vector<CVect3>&vm);// wm.size need equal to vm.size,  
};
class CSINS
{
public:
	double nts, tk;
	CEarth eth;
	CIMU imu;
	CQuat qnb;
	CMat3 Cnb, Cnb0, Cbn, Kg, Ka;
	CVect3 wib, fb, fn, an, web, wnb, att, vn, vb, pos, eb, db, _tauGyro, _tauAcc;
	CMat3 Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp;	// for etm
	CVect3 vnL, posL; CMat3 CW, MpvCnb;		// for lever
	CQuat qnbE; CVect3 vnE, posE;			// for extrapolation

	CSINS(const CQuat &qnb0 = qI, const CVect3 &vn0 = O31, const CVect3 &pos0 = O31, double tk0 = 0.0);    // initialization using quat attitude, velocity & position
	//void SetTauGA(CVect3 &tauG, CVect3 &tauA);
	void Update(const vector<CVect3> &wm, const vector<CVect3> &vm, double ts);		// SINS update using Gyro&Acc samples
	//void Extrap(const CVect3 &wm, const CVect3 &vm, double ts);						// SINS fast extrapolation using 1 Gyro&Acc sample
	void lever(const CVect3 &dL = O31);		// lever arm
	void etm(void);							// SINS error transform matrix coefficients
};

CMat3 a2mat(const CVect3 &att);
CQuat a2qua(double pitch, double roll, double yaw);	// Euler angles to quaternion
CQuat a2qua(const CVect3 &att);					// Euler angles to quaternion
CQuat rv2q(const CVect3 &rv);					// rotation vector to quaternion
CMat3 askew(const CVect3 &v);					// askew matrix;
CMat3 pos2Cen(const CVect3 &pos);				// to geographical position matrix llg->Cen
CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts = 1.0, CEarth *pEth = NULL);  // position difference to velocity
CVect3 q2att(const CQuat &qnb);		// quaternion to Euler angles 
CMat3  q2mat(const CQuat &qnb);		// quaternion to DCM
CVect3 q2rv(const CQuat &q);			// quaternion to rotation vector
CVect3 m2att(const CMat3 &Cnb);					// DCM to Euler angles 
CQuat  m2qua(const CMat3 &Cnb);					// DCM to quaternion
#endif