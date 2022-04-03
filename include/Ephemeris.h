#ifndef EPHEMERIS_H_
#define EPHEMERIS_H_
//#include "Time_utility.h"
#include "Rinex.h"
#include <vector>
//#include <map>
#include <set>
#include <string>

const double ura_value[] = {              /* ura max values */
	2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
	3072.0,6144.0
};
int uraindex(double value);

class Sat_EPH {
public:
	Sat_EPH(int SatID, vector<double> ::iterator Satnavmessage);
	//virtual void Sat_EPH_GBparaset(const double * Satnavmessage) = 0;
	void calculate_Satpos_inorbit(gtime_t time);//����˴���time�ǲ�����α��ʱ���ķ���ʱ��
	int getCRS(double & CRS);
	int getSVH() { return svh; }
	double cal_POSVAR();
	gtime_t gettoe() { return toe; }
	~Sat_EPH() {};
protected:
	/*��������λ���йز���*/
	double mu; //Ϊ�˼�������λ�õĲ���
	double omge;// RTKILIB���ú��ʾ��GPSBDS ��һ��
	double sinE;//�м����
	double tk;//����ʱ����toc֮��
	double x_inorbit;//���ƽ��λ�ã���ͨ�����ǽ�һ����������λ��
	double y_inorbit;
	double i_inorbit;
	/*��������*/
	int sat;            /* satellite number */
	int iode, iodc;      /* IODE,IODC */
	int sva;            /* SV accuracy (URA index) */
	int svh;            /* SV health (0:ok) */
	int week;           /* GPS/QZS: gps week, GAL: galileo week */
	//int code;           /* GPS/QZS: code on L2, GAL/CMP: data sources */
	gtime_t toe, toc, ttr; /* Toe,Toc,T_trans */
						/* SV orbit parameters */
	double A, e, i0, OMG0, omg, M0, deln, OMGd, idot;
	double crc, crs, cuc, cus, cic, cis;
	double toes;        /* Toe (s) in week */
	double f0, f1, f2;    /* SV clock parameters (af0,af1,af2) */
	//double tgd[2];      /* group delay parameters */
						/* GPS/QZS:tgd[0]=TGD */
						/* GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E1 */
						/* CMP    :tgd[0]=BGD1,tgd[1]=BGD2 */
	//double Adot, ndot;   /* Adot,ndot for CNAV */
};
class GPS_Sat_EPH :public Sat_EPH {
public:
	GPS_Sat_EPH(int SatID, vector<double>::iterator Satnavmessage, int SatToc);
	int get_iode() { return iode; }
	~GPS_Sat_EPH() {}
	void Cali_GPS_Satinfo(gtime_t obstime, vector<double > &SAT_Pos, double  &CLK_BIAS);
	double getTGD() { return tgd; }
private:
	int flag;
	double fit;
	int code;
	double tgd;
};
class BDS_Sat_EPH :public Sat_EPH {
public:
	BDS_Sat_EPH(int SatID, vector<double>::iterator Satnavmessage, int SatToc);
	double getdeln() { return deln; }
	~BDS_Sat_EPH() {}
	void Cali_BDS_Satinfo(int ID,gtime_t obstime, vector<double > &SAT_Pos, double  &CLK_BIAS);
	double getTGD() { return tgd[0]; }
	//void Cali_BDS_Satinfo(double obstime, vector<double > &SAT_Pos_Vel, vector<double> &CLK_BIAS_DRIFT);
private:
	double tgd[2];
};
typedef std::pair<int, GPS_Sat_EPH> gepair;
typedef std::pair<int, BDS_Sat_EPH> bepair;
typedef std::multimap <int, GPS_Sat_EPH> Gsatephmmp;
typedef std::multimap <int, BDS_Sat_EPH> Bsatephmmp;



class All_Sats_EPH {
public:
	All_Sats_EPH() {};
	All_Sats_EPH( NAV_RINEX &nav);
	void satephsolve( NAV_RINEX & nav);
	vector<double > getGPSIONCORR() { return GPS_ION_CORR; }
	vector<double > getBDSIONCORR() { return BDS_ION_CORR; }
	int Cal_OneSatinfo(int SatID,char mode, gtime_t obs_time,vector<double> & pos_vel, vector<double> &bias_drift, double &var,int &svh,double &TGD);
private:
	int LEAP_SECONDS;
	vector<double > GPS_ION_CORR;
	vector<double > BDS_ION_CORR;
	vector<double > GPS_TIME_CORR;
	vector<double > BDS_TIME_CORR;
	Gsatephmmp Gmmp;
	Bsatephmmp Bmmp;
};

#endif
