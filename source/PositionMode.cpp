#include "PositionMode.h"
int numiterator = 0;
int debug_stop = 229;
void Position_RTK(NAV_RINEX &nav, OBS_RINEX &obs, ifstream &obsfileB, ifstream &obsfileR) {
	receiver rcv;
	ofstream out;
	rcv.initializercv(nav, out);
	gtime_t headtime = obs.readheader();
	
	if (headtime.time == 0) {//找不到headtime，不往下进行。
		cout << "Something wrong with your obs file" << endl;
		return;
	}
	string buffB, buffR;
	obs.findfirstepoch(obsfileB, obsfileR, buffB, buffR, headtime);//此处已经找到了和rover第一个历元对应的base的时间
	do {
		numiterator++;
		
		if (!screent_obs(obsfileB, obsfileR, buffB, buffR, rcv, obs))
			continue;
		
		cout << "numiterator:" << numiterator << endl;
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		//rcv.updateSatinfo();
		//rcv.tdcp2vel();//TD载波相位最小二乘
		rcv.GNSS_sever(out);
	} while (getline(obsfileB, buffB) && getline(obsfileR, buffR));//readepoch中的getepochdata向下推进了obsfileB
	obsfileB.close();
	obsfileR.close();
	rcv.closefile(out);
	return;

}
void Position_TC(NAV_RINEX &nav, OBS_RINEX &obs, ifstream &obsfileB, ifstream &obsfileR,const string &imufilename) {
	gtime_t headtime = obs.readheader();
	 ofstream out;
	gtime_t LASTsoltime;
	int LS_stat = 0;
	if (headtime.time == 0) {//找不到headtime，不往下进行。
		cout << "Something wrong with your obs file" << endl;
		return;
	}
	string buffB, buffR,buffIMU;
	obs.findfirstepoch(obsfileB, obsfileR, buffB, buffR, headtime);//此处已经找到了和rover第一个历元对应的base的时间
	Fusion TCfusion;
	TCfusion.initialze_receiver(nav,out);
	ifstream imufile(imufilename);
	if (!imufile.is_open()) {
		cout << "Fail to open the imu file." << endl;
		return;
	}
	if (!find_IMU_firstepoch(imufile, buffIMU, headtime)) {
		cout << "IMU file time cannot match with the GNSS obs file" << endl;
		return;
	}
	vector<double > dwdv(6, 0);
	double imu_time,GNSS_time;
	int Fusion_Epoch_Update = 0;//是否应该往下读取一个历元的GNSS文件，为融合做准备，但是不一定能完成融合，因为可能没有GNSS数据
	int IMU_Continue_Flag = 0;

	do {
		
		if (!Fusion_Epoch_Update) {
			while (!screent_obs(obsfileB, obsfileR, buffB, buffR, TCfusion.get_receiver(), obs));
			GNSS_time = calculate_rcv_time(buffR);
			Fusion_Epoch_Update = 1;// 读取了新一个历元的GNSS数据，等待IMU追上，二者融合
			if (!(getline(obsfileB, buffB) && getline(obsfileR, buffR))) {//搭配screent_obs使用，往下读一行
				//如果无法读取，说明obsfile结束了，对于RTK/INS紧组合来说也可以结束了
				// 但是要把最后一组算完,依靠短路特性
				IMU_Continue_Flag = 1;
			}
		}
		Time_IMU_Extract(buffIMU, dwdv, imu_time);
		
		TCfusion.set_sins(dwdv,imu_time);
		TCfusion.ins_mech_ErrorTimeUpdate();//每一个IMU读取都应该进行sins解算和紧组合时间更新,此处并非TC kalman全部状态，只更新与ins有关的状态
		if (TCfusion.get_ins_align_flag()) {
				int aaaa = 1;
			
		}
		//内部维护，虽然代码在ins_align之前，但是ins_align未完成前，不会进行任何更新。
		if (Fusion_Epoch_Update) {
			double Thres_interval = 0.501 / TC_conf.sample_rate;
			
			if (fabs(imu_time - GNSS_time) > Thres_interval + 0.5 && fabs(imu_time - round(imu_time)) < Thres_interval&&TCfusion.get_ins_align_flag()
				&& fabs(imu_time-TCfusion.get_fusion_time())> Thres_interval + 0.5) {
				// Rinex Tokyo rover 4:05:50 没有星，对应没有时间戳
				//此时的GNSS_time应该是读到了 应有的整数s之后的的某个时间
				//同时IMU取值已经到了整数s 此时应该输出结果
				//ins2sol();
				TCfusion.set_ngnsslock(0);
				TCfusion.set_Solstate(3);
				TCfusion.TC_instofile(imu_time,out);
				if (TC_conf.plot_mode) {
					TCfusion.Fusion_plot();
				}
			}
			else if (fabs(imu_time - GNSS_time) < Thres_interval) {//到了IMU和GNSS重合的时间
				/*if (fabs(GNSS_time - 1545192630) < Thres_interval)
					int aaaaa = 1;*/
				
				if (!TCfusion.get_ins_align_flag()) {
					TCfusion.ins_align(out);//取消realign架构，ins_align只在gnss元进行 且只要进行完了就不会再进行
					if (TCfusion.get_ins_align_flag()) {
						TCfusion.Clear_rcv_array();
						TCfusion.TC_instofile(GNSS_time, out);
						if (TC_conf.plot_mode) {
							TCfusion.Fusion_plot();
						}
					}
					Fusion_Epoch_Update = 0;//这组数据对准完就应该往下再读取一组GNSS数据,即使对准了，也是下一历元再开始组合
					continue;
				}
				/*if (imu_time > GNSS_time) {
					cout << "hello" << endl;
				}*/
				if (numiterator >= debug_stop-1) {
					int aaaa = 1;
				}
				TCfusion.set_fusion_time(GNSS_time);
				TCfusion.rtk_ins();//只包含观测更新，里面也应加上，如果ins_align没完成，直接返回
				TCfusion.TC_instofile(GNSS_time, out);
				if (TC_conf.plot_mode) {
					TCfusion.Fusion_plot();
				}
				Fusion_Epoch_Update = 0;//更新完之后或者初始对准完就应该往下再读取一组GNSS数据	
				if (IMU_Continue_Flag) {
					//IMU_Continue_Flag = 0;
					break;
				}
			}
		}
	} while (getline(imufile, buffIMU));
	
	obsfileB.close();
	obsfileR.close();
	out.close();
	return;
}