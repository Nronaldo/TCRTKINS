#include "earth_utility.h"

double Cspeed = 299792458.0;
double D2R = PI / 180.0;
double R2D = 1 / D2R;
double RE_WGS84 = 6378137.0;
double ECC_WGS84 = 0.0818191908426215;
double FE_WGS84 = 1.0 / 298.257223563;
double OMGE = 7.2921151467E-5;
double ERR_BRDCI = 0.5;
double ERR_CBIAS = 0.3;
double ERR_SAAS = 0.3;
double REL_HUMI = 0.7;
const double lam_carr[3] = { /* carrier wave length (m) */
	Cspeed / 1.57542E9,Cspeed / 1.22760E9,Cspeed / 1.17645E9
};
const double chisqr[100] = {
	10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
	31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
	46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
	61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
	74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
	88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
	101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
	113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
	126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
	138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};

const int NSATGPS = 32;
const int NSATCMP = 45;//should change  35 45 63


/*坐标转换相关函数*/
void ecef2pos(const std::vector<double > &xecef, std::vector<double > &xllh) {
	double e2 = FE_WGS84 * (2.0 - FE_WGS84);
	double r2 = dot(xecef, 0, xecef, 0, 2), z, zk, v = RE_WGS84, sinp;
	for (z = xecef[2], zk = 0.0; fabs(z - zk) >= 1E-4;) {
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = RE_WGS84 / sqrt(1.0 - e2 * sinp*sinp);
		z = xecef[2] + v * e2*sinp;
	}
	xllh[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (xecef[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
	xllh[1] = r2 > 1E-12 ? atan2(xecef[1], xecef[0]) : 0.0;
	xllh[2] = sqrt(r2 + z * z) - v;
	return;
}
void ecef2enu(const std::vector<double > &xllh, const std::vector<double > &los, std::vector<double > &enu) {
	Eigen::Matrix3d E;
	xyz2enu(xllh, E);
	Eigen::Vector3d xecef_eigen(los.data());
	//Eigen::Vector3d enu_eigen;
	Eigen::Vector3d enu_eigen = E * xecef_eigen;
	double a;
	for (int i = 0; i < 3; i++)
		enu[i] = enu_eigen(i);
	return;
}
void Dblh2Dxyz(const std::vector<double > &blh, Eigen::Matrix3d & T) {
	double B, L, H, sinB, cosB, sinL, cosL, N, NH, e2;
	B = blh[0]; L = blh[1]; H = blh[2];
	sinB = sin(B); cosB = cos(B); sinL = sin(L); cosL = cos(L);
	N = RE_WGS84 / sqrt(1 - pow(ECC_WGS84,2) * pow(sinB,2)); 
	NH = N + H;
	e2 = ECC_WGS84* ECC_WGS84;
	T << -NH * sinB*cosL, -NH * cosB*sinL, cosB*cosL,
		-NH * sinB*sinL, NH*cosB*cosL, cosB*sinL,
		(NH - e2 * N)*cosB, 0, sinB;
	return;

}
void covecef(const std::vector<double > &posllh, const Eigen::Matrix3d &Q, Eigen::Matrix3d &Qv) {
	Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d EQ;
	xyz2enu(posllh, E);
	EQ = E.transpose()*Q;
	Qv = EQ * E;
	return;
}
void covenu(const std::vector<double > &xllh, const Eigen::Matrix3d &P, Eigen::Matrix3d &Q) {
	Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
	xyz2enu(xllh, E);
	Q = E * P*E.transpose();
	return;
}
void xyz2enu(const std::vector<double > &xllh, Eigen::Matrix3d &E) {
	double sinp = sin(xllh[0]), cosp = cos(xllh[0]), sinl = sin(xllh[1]), cosl = cos(xllh[1]);
	E << -sinl, cosl, 0.0,
		-sinp * cosl, -sinp * sinl, cosp,
		cosp * cosl, cosp * sinl, sinp;
	return;
}
Eigen::Vector3d blh2xyz(const Eigen::Vector3d &blh, Eigen::Matrix3d &E) {
	double B, L, H, N, X, Y, Z,sinB,cosB,sinL,cosL;
	B = blh[0]; 
	L = blh[1]; 
	H = blh(2);
	vector<double > xllh = {B,L,H};
	sinB = sin(B); cosB = cos(B); sinL = sin(L); cosL = cos(L);
	N = RE_WGS84 / sqrt(1 - pow(ECC_WGS84,2) * pow(sinB,2));
	X = (N + H)*cosB*cosL;
	Y = (N + H)*cosB*sinL;
	Z = (N*(1 - pow(ECC_WGS84,2)) + H)*sinB;
	Eigen::Vector3d xyz = { X, Y, Z };
	xyz2enu(xllh, E);
	return xyz;
}
double varerr(double el) {
	double fact, varr;
	fact = 1.0;
	varr = pow(TC_conf.Error_std[0], 2)*(pow(TC_conf.Error_std[1], 2) + pow(TC_conf.Error_std[2], 2) / sin(el));
	//先不考虑双频的事，只考虑单频
	return pow(fact, 2)*varr;
}
double interpc(const double coef[], double lat) {
	int i = (int)(lat / 15.0);
	if (i < 1) return coef[0]; else if (i > 4) return coef[4];
	return coef[i - 1] * (1.0 - lat / 15.0 + i) + coef[i] * (lat / 15.0 - i);
}
double tropmapf(gtime_t time, const std::vector<double > &posllh, vector<double > azel) {
	if (posllh[2]<-1000.0 || posllh[2]>20000.0) {
		return 0.0;
	}
	return nmf(time, posllh, azel);
}
double nmf(gtime_t time, const std::vector<double > &pos, vector<double > azel) {
	/* ref [5] table 3 */
	/* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
	const double coef[][5] = {
		{ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
		{ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
		{ 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

		{ 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
		{ 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
		{ 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

		{ 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
		{ 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
		{ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
	};
	const double aht[] = { 2.53E-5, 5.49E-3, 1.14E-3 }; /* height correction */
	double y, cosy, ah[3], aw[3], dm, el = azel[1], lat = pos[0] * R2D, hgt = pos[2];
	int i;
	if (el <= 0.0) {
		return 0.0;
	}
	y = (time2doy(time) - 28.0) / 365.25 + (lat < 0.0 ? 0.5 : 0.0);

	cosy = cos(2.0*PI*y);
	lat = fabs(lat);

	for (i = 0; i < 3; i++) {
		ah[i] = interpc(coef[i], lat) - interpc(coef[i + 3], lat)*cosy;
		aw[i] = interpc(coef[i + 6], lat);
	}
	/* ellipsoidal height is used instead of height above sea level */
	dm = (1.0 / sin(el) - mapf(el, aht[0], aht[1], aht[2]))*hgt / 1E3;
	

	return mapf(el, ah[0], ah[1], ah[2]) + dm;
}

double mapf(double el, double a, double b, double c) {
	double sinel = sin(el);
	return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))));
}