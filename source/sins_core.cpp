#include "sins_core.h"
const CVect3 I31(1, 1, 1), O31(0, 0, 0);
const CQuat  qI(1.0, 0, 0, 0);
const CMat3  I33(1, 0, 0, 0, 1, 0, 0, 0, 1), O33(0, 0, 0, 0, 0, 0, 0, 0, 0);
const CGLV   glv;
/*MATH function*/
int		sign(double val, double eps=EPS);
double	range(double val, double minVal, double maxVal);
double	atan2Ex(double y, double x);

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
double asinEx(double x){
	return asin(range(x, -1.0, 1.0));
}	
double acosEx(double x) {
	return acos(range(x, -1.0, 1.0));
}	
/***************************  class CGLV  *********************************/
CGLV::CGLV()
{
	Re = 6378137.0;
	f = 1.0 / 298.257;  
	wie = 7.2921151467e-5;
	g0 = 9.7803267714;
	e = sqrt(2 * f - f * f); e2 = e * e;
	mg = g0 / 1000.0;
	ug = mg / 1000.0;
	deg = PI / 180.0;
	min = deg / 60.0;
	sec = min / 60.0;
	ppm = 1.0e-6;
	hur = 3600.0;
	dps = deg / 1.0;
	dph = deg / hur;
	dpsh = deg / sqrt(hur);
	dphpsh = dph / sqrt(hur);
	ugpsHz = ug / sqrt(1.0);
	ugpsh = ug / sqrt(hur);
	mpsh = 1 / sqrt(hur);
	mpspsh = 1 / 1 / sqrt(hur);
	ppmpsh = ppm / sqrt(hur);
	secpsh = sec / sqrt(hur);
}
/***************************  class CEarth  *********************************/
CEarth::CEarth(double a0, double f0, double g0)
{
	a = a0;	f = f0; wie = glv.wie;
	b = (1 - f)*a;
	e = sqrt(a*a - b * b) / a;	e2 = e * e;
	gn = CVect3(0, 0, -g0);
}
void CEarth::Update(const CVect3 &pos, const CVect3 &vn)
{
	this->pos = pos;  this->vn = vn;
	sl = sin(pos[0]), cl = cos(pos[0]), tl = sl / cl;
	double sq = 1 - e2 * sl*sl, sq2 = sqrt(sq);
	RMh = a * (1 - e2) / sq / sq2 + pos[2];	f_RMh = 1.0 / RMh;
	RNh = a / sq2 + pos[2];    clRNh = cl * RNh;  f_RNh = 1.0 / RNh; f_clRNh = 1.0 / clRNh;
	wnie[0] = 0.0, wnie[1] = wie * cl, wnie[2] = wie * sl;
	wnen[0] = -vn[1]*f_RMh, wnen[1] = vn[0]*f_RNh, wnen[2] = wnen[1]*tl;
	wnin = wnie + wnen;
	sl2 = sl * sl, sl4 = sl2 * sl2;
	gn[2] = -(glv.g0*(1 + 5.27094e-3*sl2 + 2.32718e-5*sl4) - 3.086e-6*pos[2]);
	gcc = gn - (wnie + wnin)*vn;
}
CVect3 CEarth::vn2dpos(const CVect3 &vn, double ts) const
{
	return CVect3(vn[1]*f_RMh, vn[0]*f_clRNh, vn[2])*ts;// p_dot = M_pv * vn = [0 1/RMh 0; 1/(cosL*RNh),0,0; 0, 0,1]* vn. V_n+0.5 = (V_n+V_n+1)/2
}
/*CVect3 友元函数*/
istream & operator>>(istream & is, CVect3 &cv3) {
	is >> cv3[0]>> cv3[1]>> cv3[2];
	return is;
}
ostream & operator<<(ostream & os, const CVect3 & cv3) {
	os << setprecision(14) << cv3[0]<<endl;
	os << setprecision(14) << cv3[1]<< endl;
	os << setprecision(14) << cv3[2]<< endl;
	return os;
}	
CVect3 operator+(const CVect3 &lhs, const CVect3 &rhs) {
	CVect3 r = lhs;
	r += rhs;
	return r;
}
CVect3 operator-(const CVect3 &lhs, const CVect3 &rhs) {
	CVect3 r = lhs;
	r -= rhs;
	return r;
}
CVect3 operator*(const CVect3 &lhs, const CVect3 &rhs) {
	CVect3 r(lhs[1]*rhs[2]- lhs[2] * rhs[1], lhs[2] * rhs[0] - lhs[0] * rhs[2], lhs[0] * rhs[1] - lhs[1] * rhs[0]);
	return r;
}
CVect3 operator*(const CVect3 &lhs, double f) {
	CVect3 r = lhs;
	r *= f;
	return r;
}
CVect3 operator/(const CVect3 & lhs, double f) {
	CVect3 r = lhs;
	r /= f;
	return r;
}
CVect3 operator*(double f, const CVect3 & rhs) {
	CVect3 r = rhs;
	r *= f;
	return r;
}
CVect3 operator-(const CVect3 &cv3) {
	CVect3 r(-cv3[0], -cv3[1], -cv3[2]);
	return r;
}
double norm(const CVect3 &v) {
	return sqrt(dot(v,v));
}
double normXY(const CVect3 &v) {
	return sqrt(v[0] * v[0] + v[1] * v[1]);
}
CVect3 sqrt(const CVect3 &v) {
	return CVect3(sqrt(v[0]), sqrt(v[1]), sqrt(v[2]));
}
double dot(const CVect3 &v1,const CVect3 &v2) {
	return(v1[0]*v2[0]+ v1[1] * v2[1]+ v1[2] * v2[2]);
}
/*CVect3 成员函数*/
double &CVect3::operator[](int n) {
	if (n >= 3 || n < 0) {
		cout << "The CVect3 index is out of range." << endl;
		return elements[0];
	}
	return *(this->elements + n);
}
const double & CVect3::operator[](int n) const {
	if (n >= 3 || n < 0) {
		cout << "The CVect3 index is out of range." << endl;
		return elements[0];
	}
	return *(this->elements + n);
}
CVect3 &CVect3::operator=(const CVect3 & cv3) {
	for (int i = 0; i < 3; i++) {
		this->elements[i] = cv3[i];
	}
	return *this;
}
CVect3 & CVect3::operator+=(const CVect3 & cv3) {
	for (int i = 0; i < 3; i++) {
		this->elements[i] += cv3[i];
	}
	return *this;
}
CVect3 & CVect3::operator-=(const CVect3& cv3) {
	for (int i = 0; i < 3; i++) {
		this->elements[i] -= cv3[i];
	}
	return *this;
}
CVect3 & CVect3::operator*=(double f) {
	for (int i = 0; i < 3; i++) {
		this->elements[i] *= f;
	}
	return *this;
}
CVect3 & CVect3::operator/=(double f) {
	for (int i = 0; i < 3; i++) {
		this->elements[i] /= f;
	}
	return *this;
}
bool CVect3::IsZero()const {
	for (int i = 0; i < 3; i++) {
		if (this->elements[i]<-EPS || this->elements[i] > EPS) {
			return false;
		}
	}
	return true;
}
bool CVect3::IsZeroXY()const {
	for (int i = 0; i < 2; i++) {
		if (this->elements[i]<-EPS || this->elements[i] > EPS) {
			return false;
		}
	}
	return true;
}

/*CQuat 友元*/
ostream & operator<< (ostream &os, CQuat q) {
	os << setprecision(14) << q[0] << endl;
	os << setprecision(14) << q[1] << endl;
	os << setprecision(14) << q[2] << endl;
	os << setprecision(14) << q[3] << endl;
	return os;
}
istream & operator>>(istream & is, CQuat q) {
	is >> q[0] >> q[1] >> q[2]>>q[3];
	return is;
}
CQuat operator*(const CQuat & q1, const CQuat & q2) {
	CQuat q;
	q[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
	q[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
	q[2] = q1[0] * q2[2] + q1[2] * q2[0] + q1[3] * q2[1] - q1[1] * q2[3];
	q[3] = q1[0] * q2[3] + q1[3] * q2[0] + q1[1] * q2[2] - q1[2] * q2[1];
	return q;
}
CQuat operator-(const CQuat & lhs, const CVect3 & v_phi) {
	CQuat temp = lhs;
	temp -= v_phi;
	return temp;
}
CVect3 operator-(const CQuat & q_nominal, const CQuat & q_true) {
	CQuat dq = q_true* (~q_nominal);
	return q2rv(dq);
}
CVect3 operator*(const CQuat & q1, const CVect3 & v2) {
	CQuat  q_v2(0.0,v2[0],v2[1],v2[2]);
	CQuat q_temp = q1 * q_v2*(~q1);
	CVect3 result(q_temp[1],q_temp[2],q_temp[3]);
	return result;
}
CQuat operator~(const CQuat &q) {
	return CQuat(q[0], -q[1], -q[2], -q[3]);
}

/*CQuat 成员*/
double &CQuat::operator[](int n) {
	if (n >= 4 || n < 0) {
		cout << "The CQuat index is out of range." << endl;
		return elements[0];
	}
	return *(this->elements + n);
}
const double &CQuat::operator[](int n) const {
	if (n >= 4 || n < 0) {
		cout << "The CQuat index is out of range." << endl;
		return elements[0];
	}
	return *(this->elements + n);
}
CQuat & CQuat::operator=(const CQuat & q){
	for (int i = 0; i < 4; i++)
		*(this->elements + i) = q[i];
	return *this;
}
CQuat & CQuat::operator*=(const CQuat & q) {
	return (*this = *this * q);
}
CQuat & CQuat::operator-=(const CVect3 & v_phi) {
	CQuat q = rv2q(v_phi);
	return (*this=q*(*this));
}
void CQuat:: normlize() {
	double sum = 0;
	for (int i = 0; i < 4; i++)
		sum += elements[i];
	sum=sqrt(sum);
	for (int i = 0; i < 4; i++)
		elements[i] /= sum;
	return;
}

/*CMat3 */
ostream & operator<<(ostream &os, CMat3 & m) {
	os << setprecision(14) << m.elements[0]<<',' << setprecision(14) << m.elements[1]<<',' << setprecision(14) << m.elements[2] << endl;
	os << setprecision(14) << m.elements[3] << ',' << setprecision(14) << m.elements[4] << ',' << setprecision(14) << m.elements[5] << endl;
	os << setprecision(14) << m.elements[6] << ',' << setprecision(14) << m.elements[7] << ',' << setprecision(14) << m.elements[8] << endl;
	return os;
}
istream & operator<<(istream &is, CMat3 & m) {
	is >> m.elements[0] >> m.elements[1] >> m.elements[2] >> m.elements[3] >> m.elements[4] >> m.elements[5] >> m.elements[6] >> m.elements[7] >> m.elements[8];
	return is;
}
CMat3 operator+(const CMat3 &m1, const CMat3 & m2) {
	CMat3 r = m1;
	r += m2;
	return r;
}
CMat3 operator-(const CMat3 &m1, const CMat3 & m2) {
	CMat3 r = m1;
	r -= m2;
	return r;
}
CMat3 operator-(const CMat3 &m) {
	CMat3 r;
	for (int i = 0; i < 9; i++)
		r.elements[i] = -m.elements[i];
	return r;
}
CMat3 operator*(const CMat3 &m1, const CMat3 & m2) {
	CMat3 r;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k <3; k++) {
				r(i,j) += m1(i,k) * m2(k,j);
			}
		}
	}
	return r;
}
CMat3 operator*(const CMat3 &m, double f) {
	CMat3 r;
	for (int i = 0; i < 9; i++) {
		r.elements[i] = m.elements[i] * f;
	}
	return r;
}
CMat3 operator*(double f, CMat3 &m) {
	CMat3 r;
	for (int i = 0; i < 9; i++) {
		r.elements[i] = m.elements[i] * f;
	}
	return r;
}
CVect3 operator*(const CMat3 &m, const CVect3 & v) {
	CVect3 v_r(0,0,0);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			v_r[i] += m.elements[i*3+j] * v[j];
		}
	}
	return v_r;
}
CMat3 operator~(const CMat3 &m) {
	return CMat3(m(0,0), m(1,0), m(2,0), m(0,1), m(1,1), m(2,1), m(0,2), m(1,2), m(2,2));
}
double det(const CMat3 &m) {
	return m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1)) - m(0,1)*(m(1,0)*m(2,2) - m(1,2)*m(2,0)) + m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));
}
CMat3 inv(const CMat3 &m) {
	double nm;
	nm = m(0,0)*(m(1,1)*m(2,2) - m(1, 2)*m(2,1)) - m(0, 1)*(m(1, 0)*m(2,2) - m(1,2)*m(2, 0)) + m(0, 2)*(m(1, 0)*m(2, 1) - m(1, 1)*m(2, 0));
	CMat3 mtmp;
	mtmp(0, 0) = (m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)) / nm;
	mtmp(1,0) = -(m(1, 0)*m(2, 2) - m(1, 2)*m(2, 0)) / nm;
	mtmp(2,0) = (m(1, 0)*m(2, 1) - m(1, 1)*m(2, 0)) / nm;
	mtmp(0, 1) = -(m(0, 1)*m(2, 2) - m(0, 2)*m(2, 1)) / nm;
	mtmp(1, 1) = (m(0, 0)*m(2, 2) - m(0, 2)*m(2, 0)) / nm;
	mtmp(2, 1) = -(m(0, 0)*m(2, 1) - m(0, 1)*m(2, 0)) / nm;
	mtmp(0,2) = (m(0, 1)*m(1, 2) - m(0, 2)*m(1, 1)) / nm;
	mtmp(1, 2) = -(m(0, 0)*m(1, 2) - m(0, 2)*m(1, 0)) / nm;
	mtmp(2, 2) = (m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0)) / nm;
	return mtmp;
}
CVect3 diag(const CMat3 &m) {
	return CVect3(m(0,0), m(1,1), m(2,2));
}
CMat3 diag(const CVect3 &v) {
	return CMat3(v[0], 0, 0, 0, v[1], 0, 0, 0, v[2]);
}
/*成员*/
CMat3::CMat3(double xx, double xy, double xz,
	double yx, double yy, double yz,
	double zx, double zy, double zz) {
	elements[0] = xx;
	elements[1] = xy;
	elements[2] = xz;
	elements[3] = yx;
	elements[4] = yy;
	elements[5] = yz;
	elements[6] = zx;
	elements[7] = zy;
	elements[8] = zz;
}
CMat3::CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2) {
	for (int i = 0; i < 3; i++) {
		elements[i * 3 + 0] = v0[i];
		elements[i * 3 + 1] = v1[i];
		elements[i * 3 + 2] = v2[i];
	}
}
double &CMat3::operator()(int r, int c) {
	if (r < 0 || c < 0 || r >= 3 || c >= 3) {
		cout << "CMat3 index out of range." << endl;
		return elements[0];
	}
	return *(this->elements + r * 3 + c);
}
const double &CMat3::operator()(int r, int c)const {
	if (r < 0 || c < 0 || r >= 3 || c >= 3) {
		cout << "CMat3 index out of range." << endl;
		return elements[0];
	}
	return *(this->elements + r * 3 + c);
}
CMat3& CMat3::operator=(const CMat3 & m) {
	for (int i = 0; i < 9; i++)
		this->elements[i] = m.elements[i];
	return *this;
}
CMat3& CMat3::operator+=(const CMat3 & m) {
	for (int i = 0; i < 9; i++)
		this->elements[i] += m.elements[i];
	return *this;
}
CMat3& CMat3::operator-=(const CMat3 & m) {
	for (int i = 0; i < 9; i++)
		this->elements[i] -= m.elements[i];
	return *this;
}
/*CIMU*/
void CIMU::Update(const vector<CVect3> &wm, const vector<CVect3>&vm)
{// 单子样，除非指定，否则选取单子样+前一周期
	// 多子样，选择圆锥误差修正，而非多项式修正
	//优化方式：用圆锥运动，来补偿，多项式方式：用多项式拟合，前者适合正弦角振动，比后者强
	// 因为有噪声，会导致多项式很难拟合
	// 归根结底两者的目的都是为了补偿不可交换
	if (wm.empty() || vm.empty() || wm.size() != vm.size()) {
		cout << "The subsamples go wrong." << endl;
		return;
	}
	int nSamples = wm.size();
	static double conefactors[5][4] = {				// coning coefficients
		{2. / 3},										// 2
		{9. / 20, 27. / 20},							// 3
		{54. / 105, 92. / 105, 214. / 105},				// 4
		{250. / 504, 525. / 504, 650. / 504, 1375. / 504}	// 5
	};
	int i;
	double *pcf = conefactors[nSamples - 2];
	CVect3 cm(0, 0, 0), sm(0, 0, 0), wmm(0, 0, 0), vmm(0, 0, 0);

	if (nSamples == 1)  // one-plus-previous sample
	{
		if (prefirst == 1) { 
			wm_1 = wm[0]; 
			vm_1 = vm[0];
			prefirst = 0; 
		}
		cm = 1.0 / 12 * wm_1;
		wm_1 = wm[0];
		sm = 1.0 / 12 * vm_1;
		vm_1 = vm[0];
	}
	if (nSamples > 1) 
		prefirst = 1;
	for (i = 0; i < nSamples - 1; i++)
	{
		cm += pcf[i] * wm[i];
		sm += pcf[i] * vm[i];
		wmm += wm[i];
		vmm += vm[i];
	}
	wmm += wm[i];
	vmm += vm[i];
	phim = wmm + cm * wm[i];
	dvbm = vmm + 1.0 / 2 * wmm*vmm + (cm*vm[i] + sm * wm[i]);
}
/***************************  class CSINS  *********************************/
CSINS::CSINS(const CQuat &qnb0, const CVect3 &vn0, const CVect3 &pos0, double tk0)
{
	tk = tk0;  nts = 0.0;
	qnb = qnb0;	vn = vn0, pos = pos0;  	qnbE = qnb, vnE = vn, posE = pos;
	eth.Update(pos0, vn0);
	Cnb = q2mat(qnb); att = m2att(Cnb); Cnb0 = Cnb; Cbn = ~Cnb; vb = Cbn * vn;
	Kg = Ka = I33; eb = db = _tauGyro = _tauAcc = O31;
	wib = fb = fn = an = wnb = web = O31;
	Maa = Mav = Map = Mva = Mvv = Mvp = Mpv = Mpp = O33;
	etm(); lever();
}
void CSINS::Update(const vector<CVect3> &wm, const vector<CVect3> &vm, double ts)
{
	if (wm.empty() || vm.empty() || wm.size() != vm.size()) {
		cout << "The subsamples go wrong." << endl;
		return;
	}
	int nSamples = wm.size();
	nts = nSamples * ts;	
	tk += nts;
	double nts2 = nts / 2;
	imu.Update(wm, vm);
	imu.phim = Kg * imu.phim - eb * nts; imu.dvbm = Ka * imu.dvbm - db * nts;  // IMU calibration
	CVect3 vn01 = vn + an * nts2, pos01 = pos + eth.vn2dpos(vn01, nts2);
	eth.Update(pos01, vn01);
	wib = imu.phim / nts; fb = imu.dvbm / nts;
	web = wib - Cbn * eth.wnie;
	//wnb = wib - (qnb*rv2q(imu.phim / 2))*eth.wnin;
	//wnb = wib - (Cnb*rv2m(imu.phim / 2))*eth.wnin;
	wnb = wib - Cbn * eth.wnin;
	fn = qnb * fb;
	an = rv2q(-eth.wnin*nts2)*fn + eth.gcc;
	CVect3 vn1 = vn + an * nts;
	pos = pos + eth.vn2dpos(vn + vn1, nts2);	vn = vn1;
	Cnb0 = Cnb;
	qnb = rv2q(-eth.wnin*nts)*qnb*rv2q(imu.phim);
	Cnb = q2mat(qnb); att = m2att(Cnb); Cbn = ~Cnb; vb = Cbn * vn;
}
void CSINS::lever(const CVect3 &dL)
{
	Mpv = CMat3(0, eth.f_RMh, 0, eth.f_clRNh, 0, 0, 0, 0, 1);
	CW = Cnb * askew(web), MpvCnb = Mpv * Cnb;
	vnL = vn + CW * dL; posL = pos + MpvCnb * dL;
}

void CSINS::etm(void)
{
	double tl = eth.tl, secl = 1.0 / eth.cl, secl2 = secl * secl,
		wN = eth.wnie[1], wU = eth.wnie[2], vE = vn[0], vN = vn[1];
	double f_RMh = eth.f_RMh, f_RNh = eth.f_RNh, f_clRNh = eth.f_clRNh,
		f_RMh2 = f_RMh * f_RMh, f_RNh2 = f_RNh * f_RNh;
	CMat3 Avn = askew(vn),
		Mp1(0, 0, 0, -wU, 0, 0, wN, 0, 0),
		Mp2(0, 0, vN*f_RMh2, 0, 0, -vE * f_RNh2, vE*secl2*f_RNh, 0, -vE * tl*f_RNh2);
	Maa = askew(-eth.wnin);
	Mav = CMat3(0, -f_RMh, 0, f_RNh, 0, 0, tl*f_RNh, 0, 0);
	Map = Mp1 + Mp2;
	Mva = askew(fn);
	Mvv = Avn * Mav - askew(eth.wnie + eth.wnin);
	Mvp = Avn * (Mp1 + Map);
	double scl = eth.sl*eth.cl;
	Mvp(2,0) = Mvp(2, 0) - glv.g0*(5.27094e-3 * 2 * scl + 2.32718e-5 * 4 * eth.sl2*scl); 
	Mvp(2,2) = Mvp(2,2) + 3.086e-6;
	Mpv = CMat3(0, f_RMh, 0, f_clRNh, 0, 0, 0, 0, 1);
	Mpp = CMat3(0, 0, -vN * f_RMh2, vE*tl*f_clRNh, 0, -vE * secl*f_RNh2, 0, 0, 0);
}
/*转换函数*/
CMat3 a2mat(const CVect3 &att) {
	double	si = sin(att[0]), ci = cos(att[0]),
			sj = sin(att[1]), cj = cos(att[1]),
			sk = sin(att[2]), ck = cos(att[2]);
	CMat3 Cnb;
	Cnb(0,0) = cj * ck - si * sj*sk;	Cnb(0,1) = -ci * sk;	Cnb(0,2) = sj * ck + si * cj*sk;
	Cnb(1,0) = cj * sk + si * sj*ck;	Cnb(1,1) = ci * ck;	Cnb(1,2) = sj * sk - si * cj*ck;
	Cnb(2,0) = -ci * sj;				Cnb(2,1) = si;		Cnb(2,2) = ci * cj;
	return Cnb;
}
CQuat a2qua(double pitch, double roll, double yaw) {	// Euler angles to quaternion// 类似姿态矩阵，同样是三个轴，分别列出四元数，从左到右依次相乘
	pitch /= 2.0, roll /= 2.0, yaw /= 2.0;
	double	sp = sin(pitch), sr = sin(roll), sy = sin(yaw),
		cp = cos(pitch), cr = cos(roll), cy = cos(yaw);
	CQuat qnb;
	qnb[0] = cp * cr*cy - sp * sr*sy;
	qnb[1] = sp * cr*cy - cp * sr*sy;
	qnb[2] = cp * sr*cy + sp * cr*sy;
	qnb[3] = cp * cr*sy + sp * sr*cy;
	return qnb;
}
CQuat a2qua(const CVect3 &att) {					// Euler angles to quaternion
	return a2qua(att[0], att[1], att[2]);
}
CQuat rv2q(const CVect3 &rv) { // rotation vector to quaternion
	double F1 = 2 * 1;		// define: Fk=2^k*k! 
	double F2 = F1 * 2 * 2;
	double F3 = F2 * 2 * 3;
	double F4 = F3 * 2 * 4;
	double F5 = F3 * 2 * 5;
	double n2 = rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2], c, f;
	if (n2 < (PI / 180.0*PI / 180.0))	// 0.017^2 
	{
		double n4 = n2 * n2;
		c = 1.0 - n2 * (1.0 / F2) + n4 * (1.0 / F4);
		f = 0.5 - n2 * (1.0 / F3) + n4 * (1.0 / F5);
	}
	else
	{
		double n_2 = sqrt(n2) / 2.0;
		c = cos(n_2);
		f = sin(n_2) / n_2 * 0.5;
	}
	return CQuat(c, f*rv[0], f*rv[1], f*rv[2]);
}					
CMat3 askew(const CVect3 &v) {					// askew matrix;
	return CMat3(0, -v[2], v[1],
				v[2], 0.0, -v[0],
				-v[1], v[0], 0);
}
CMat3 pos2Cen(const CVect3 &pos) {				// to geographical position matrix
	double si = sin(pos[0]), ci = cos(pos[0]), sj = sin(pos[1]), cj = cos(pos[1]);
	return CMat3(-sj, -si * cj, ci*cj,
		cj, -si * sj, ci*sj,
		0, ci, si);	//Cen
}
CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts = 1.0, CEarth *pEth = NULL) {  // position difference to velocity
	double sl, cl, sl2, sq, sq2, RMh, RNh, clRNh;
	if (pEth)
	{
		RMh = pEth->RMh; clRNh = pEth->clRNh;
	}
	else
	{
		sl = sin(pos0[0]); cl = cos(pos0[0]); sl2 = sl * sl;
		sq = 1 - glv.e2*sl2; sq2 = sqrt(sq);
		RMh = glv.Re*(1 - glv.e2) / sq / sq2 + pos0[2];
		RNh = glv.Re / sq2 + pos0[2];    clRNh = cl * RNh;
	}
	CVect3 vn = pos1 - pos0;
	return CVect3(vn[1]*clRNh / ts, vn[0]*RMh / ts, vn[2] / ts);
}
CVect3 q2att(const CQuat &qnb) {		// quaternion to Euler angles 
	double	q11 = qnb[0]*qnb[0], q12 = qnb[0] *qnb[1], q13 = qnb[0] *qnb[2], q14 = qnb[0] *qnb[3],
	q22 = qnb[1] *qnb[1], q23 = qnb[1] *qnb[2], q24 = qnb[1]*qnb[3],
	q33 = qnb[2] *qnb[2], q34 = qnb[2] *qnb[3],
	q44 = qnb[3] *qnb[3];
	CVect3 att;
	att[0] = asinEx(2 * (q34 + q12));
	att[1] = atan2Ex(-2 * (q24 - q13), q11 - q22 - q33 + q44);
	att[2] = atan2Ex(-2 * (q23 - q14), q11 - q22 + q33 - q44);
	return att;
}
CMat3  q2mat(const CQuat &qnb) {		// quaternion to DCM
	double	q11 = qnb[0]*qnb[0], q12 = qnb[0] *qnb[1], q13 = qnb[0] *qnb[2], q14 = qnb[0] *qnb[3],
		q22 = qnb[1] *qnb[1], q23 = qnb[1] *qnb[2], q24 = qnb[1] *qnb[3],
		q33 = qnb[2] *qnb[2], q34 = qnb[2] *qnb[3],
		q44 = qnb[3] *qnb[3];
	CMat3 Cnb;
	Cnb(0,0) = q11 + q22 - q33 - q44, Cnb(0, 1) = 2 * (q23 - q14), Cnb(0, 2) = 2 * (q24 + q13),
		Cnb(1, 0) = 2 * (q23 + q14), Cnb(1, 1) = q11 - q22 + q33 - q44, Cnb(1, 2) = 2 * (q34 - q12),
		Cnb(2, 0) = 2 * (q24 - q13), Cnb(2, 1) = 2 * (q34 + q12), Cnb(2, 2) = q11 - q22 - q33 + q44;
	return Cnb;
}
CVect3 q2rv(const CQuat &q) { // quaternion to rotation vector
	CQuat dq;
	dq = q;
	if (dq[0] < 0) { for (int i = 0; i < 4; i++) dq[i] = -dq[i]; }
	if (dq[0] > 1.0) dq[0] = 1.0;
	double n2 = acos(dq[0]), f;
	if (n2 > 1.0e-20)
	{
		f = 2.0 / (sin(n2) / n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq[1], dq[2], dq[3])*f;
}	
CVect3 m2att(const CMat3 &Cnb)
{
	CVect3 att;
	att[0] = asinEx(Cnb(2,1));
	att[1] = atan2Ex(-Cnb(2,0), Cnb(2,2));
	att[2] = atan2Ex(-Cnb(0,1), Cnb(1,1));
	return att;
}
CQuat m2qua(const CMat3 &Cnb)
{
	double q0, q1, q2, q3, qq4;
	if (Cnb(0,0) >= Cnb(1,1) + Cnb(2,2))
	{
		q1 = 0.5*sqrt(1 + Cnb(0, 0) - Cnb(1, 1) - Cnb(2, 2));  qq4 = 4 * q1;
		q0 = (Cnb(2, 1) - Cnb(1, 2)) / qq4; q2 = (Cnb(0, 1) + Cnb(1, 0)) / qq4; q3 = (Cnb(0,2) + Cnb(2,0)) / qq4;
	}
	else if (Cnb(1, 1) >= Cnb(0, 0) + Cnb(2, 2))
	{
		q2 = 0.5*sqrt(1 - Cnb(0, 0) + Cnb(1, 1) - Cnb(2, 2));  qq4 = 4 * q2;
		q0 = (Cnb(0, 2) - Cnb(2, 0)) / qq4; q1 = (Cnb(0, 1) - Cnb(1, 0)) / qq4; q3 = (Cnb(1, 2) + Cnb(2, 1)) / qq4;
	}
	else if (Cnb(2,2) >= Cnb(0, 0) + Cnb(1, 1))
	{
		q3 = 0.5*sqrt(1 - Cnb(0, 0) - Cnb(1, 1) + Cnb(2, 2));  qq4 = 4 * q3;
		q0 = (Cnb(1, 0) - Cnb(0, 1)) / qq4; q1 = (Cnb(0, 2) + Cnb(2, 0)) / qq4; q2 = (Cnb(1, 2) + Cnb(2, 1)) / qq4;
	}
	else
	{
		q0 = 0.5*sqrt(1 + Cnb(0, 0) + Cnb(1, 1) + Cnb(2, 2));  qq4 = 4 * q0;
		q1 = (Cnb(2, 1) - Cnb(1, 2)) / qq4; q2 = (Cnb(0, 2) - Cnb(2, 0)) / qq4; q3 = (Cnb(1, 0) - Cnb(0, 1)) / qq4;
	}
	double nq = sqrt(q0*q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 /= nq; q1 /= nq; q2 /= nq; q3 /= nq;
	return CQuat(q0, q1, q2, q3);
}