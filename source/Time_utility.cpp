#include "Time_utility.h"

gtime_t epoch2time(const double *ep) {
	const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
	gtime_t time(0,0);
	int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];
	if (year < 1970 || 2099 < year || mon < 1 || 12 < mon) return time;
	days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
	sec = (int)floor(ep[5]);
	time.time = (int)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
	time.sec = ep[5] - sec;
	return time;
}
gtime_t gpst2time(int week, double sec)
{
	gtime_t t = epoch2time(gpst0);

	if (sec < -1E9 || 1E9 < sec) sec = 0.0;
	t.time += (int)86400 * 7 * week + (int)sec;
	t.sec = sec - (int)sec;
	return t;
}

gtime_t bdt2time(int week, double sec)
{
	gtime_t t = epoch2time(bdt0);

	if (sec < -1E9 || 1E9 < sec) sec = 0.0;
	t.time += (int)86400 * 7 * week + (int)sec;
	t.sec = sec - (int)sec;
	return t;
}
gtime_t bdt2gpst(gtime_t t)
{
	return timeadd(t, 14.0);
}

gtime_t adjweek(gtime_t t, gtime_t t0)
{
	double tt = timediff(t, t0);
	if (tt < -302400.0) return timeadd(t, 604800.0);
	if (tt > 302400.0) return timeadd(t, -604800.0);
	return t;
}


gtime_t timeadd(gtime_t t, double sec)
{
	double tt;
	t.sec += sec; 
	tt = floor(t.sec); 
	t.time += (int)tt; 
	t.sec -= tt;
	return t;
}
double timediff(gtime_t t1, gtime_t t2)
{
	return t1.time-t2.time + t1.sec - t2.sec;
}

double time2gpst(gtime_t t, int *week) {
	gtime_t t0(315964800, 0);
	int sec = t.time - t0.time;
	int w = (int)(sec / 86400 / 7);
	if (week)
		*week = w;
	return (double)(sec - (double)w * 86400 * 7) + t.sec;
}

gtime_t string2time(std::string str) {
	double ymdhmsp[6];
	stringstream ss(str);
	string splitstr;
	gtime_t result (0,0);
	int i = 0;
	while (ss >> splitstr) {
		stringstream ss_temp;
		ss_temp << splitstr;
		ss_temp >> ymdhmsp[i];
		i++;
	}
	int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
	if (ymdhmsp[0] < 1970 or ymdhmsp[0]> 2099 or ymdhmsp[1] < 1 or ymdhmsp[1]>12) {
		cout << "time format is wrong" << endl;
		return result;
	}
	else {
		int days = (ymdhmsp[0] - 1970) * 365 + (ymdhmsp[0] - 1969) / 4 + doy[(int)ymdhmsp[1] - 1] + ymdhmsp[2] - 2 + ((int)ymdhmsp[0] % 4 == 0 and ymdhmsp[1] >= 3 ? 1 : 0);
		int secs = (int)floor(ymdhmsp[5]);
		result.time = (int)days * 86400 + (int)ymdhmsp[3] * 3600 + (int)ymdhmsp[4] * 60 + secs;
		result.sec = ymdhmsp[5] - secs;
		return result;
	}
}
double time2doy(gtime_t t) {
	double ep[6];
	time2epoch(t, ep);
	ep[1] = ep[2] = 1.0; ep[3] = ep[4] = ep[5] = 0.0;
	return timediff(t, epoch2time(ep)) / 86400.0 + 1.0;
}
void time2epoch(gtime_t t, double *ep) {
	const int mday[] = { /* # of days in a month */
	   31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
	   31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
	};
	int days, sec, mon, day;

	/* leap year if year%4==0 in 1901-2099 */
	days = (int)(t.time / 86400);
	sec = (int)(t.time - (time_t)days * 86400);
	for (day = days % 1461, mon = 0; mon < 48; mon++) {
		if (day >= mday[mon]) day -= mday[mon]; else break;
	}
	ep[0] = 1970 + days / 1461 * 4 + mon / 12; ep[1] = mon % 12 + 1; ep[2] = day + 1;
	ep[3] = sec / 3600; ep[4] = sec % 3600 / 60; ep[5] = sec % 60 + t.sec;
}

void time2str(gtime_t t, string &str,int n) {
	if (1.0 - t.sec < 0.5 / pow(10.0, n)) { 
		t.time++; 
		t.sec = 0.0; 
	}
	double ep[6];
	time2epoch(t,ep);
	stringstream ss;
	ss << setw(4) << ep[0] << '/';
	ss << setfill('0') << setw(2) << ep[1] << '/';
	ss << setfill('0') << setw(2) << ep[2] << ' ';
	ss << setfill('0') << setw(2) << ep[3] << ':';
	ss << setfill('0') << setw(2) << ep[4] << ':';
	ss << setfill('0') << setw(2)<<(int)ep[5]<<'.';
	int a = (ep[5] - (int)ep[5])*pow(10,n);
	ss << setfill('0') << setw(n) << a;

	str = ss.str();
	return;
}
double calculate_rcv_time(string bufferR) {
	gtime_t times = string2time(bufferR.substr(2, 27));
	return (times.time + times.sec); 
}