#include "PositionMode.h"
int numiterator = 0;
int debug_stop = 229;
void Position_RTK(NAV_RINEX &nav, OBS_RINEX &obs, ifstream &obsfileB, ifstream &obsfileR) {
	receiver rcv;
	ofstream out;
	rcv.initializercv(nav, out);
	gtime_t headtime = obs.readheader();
	
	if (headtime.time == 0) {//�Ҳ���headtime�������½��С�
		cout << "Something wrong with your obs file" << endl;
		return;
	}
	string buffB, buffR;
	obs.findfirstepoch(obsfileB, obsfileR, buffB, buffR, headtime);//�˴��Ѿ��ҵ��˺�rover��һ����Ԫ��Ӧ��base��ʱ��
	do {
		numiterator++;
		
		if (!screent_obs(obsfileB, obsfileR, buffB, buffR, rcv, obs))
			continue;
		
		cout << "numiterator:" << numiterator << endl;
		if (numiterator >= debug_stop) {
			int aaaa = 1;
		}
		//rcv.updateSatinfo();
		//rcv.tdcp2vel();//TD�ز���λ��С����
		rcv.GNSS_sever(out);
	} while (getline(obsfileB, buffB) && getline(obsfileR, buffR));//readepoch�е�getepochdata�����ƽ���obsfileB
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
	if (headtime.time == 0) {//�Ҳ���headtime�������½��С�
		cout << "Something wrong with your obs file" << endl;
		return;
	}
	string buffB, buffR,buffIMU;
	obs.findfirstepoch(obsfileB, obsfileR, buffB, buffR, headtime);//�˴��Ѿ��ҵ��˺�rover��һ����Ԫ��Ӧ��base��ʱ��
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
	int Fusion_Epoch_Update = 0;//�Ƿ�Ӧ�����¶�ȡһ����Ԫ��GNSS�ļ���Ϊ�ں���׼�������ǲ�һ��������ںϣ���Ϊ����û��GNSS����
	int IMU_Continue_Flag = 0;

	do {
		
		if (!Fusion_Epoch_Update) {
			while (!screent_obs(obsfileB, obsfileR, buffB, buffR, TCfusion.get_receiver(), obs));
			GNSS_time = calculate_rcv_time(buffR);
			Fusion_Epoch_Update = 1;// ��ȡ����һ����Ԫ��GNSS���ݣ��ȴ�IMU׷�ϣ������ں�
			if (!(getline(obsfileB, buffB) && getline(obsfileR, buffR))) {//����screent_obsʹ�ã����¶�һ��
				//����޷���ȡ��˵��obsfile�����ˣ�����RTK/INS�������˵Ҳ���Խ�����
				// ����Ҫ�����һ������,������·����
				IMU_Continue_Flag = 1;
			}
		}
		Time_IMU_Extract(buffIMU, dwdv, imu_time);
		
		TCfusion.set_sins(dwdv,imu_time);
		TCfusion.ins_mech_ErrorTimeUpdate();//ÿһ��IMU��ȡ��Ӧ�ý���sins����ͽ����ʱ�����,�˴�����TC kalmanȫ��״̬��ֻ������ins�йص�״̬
		if (TCfusion.get_ins_align_flag()) {
				int aaaa = 1;
			
		}
		//�ڲ�ά������Ȼ������ins_align֮ǰ������ins_alignδ���ǰ����������κθ��¡�
		if (Fusion_Epoch_Update) {
			double Thres_interval = 0.501 / TC_conf.sample_rate;
			
			if (fabs(imu_time - GNSS_time) > Thres_interval + 0.5 && fabs(imu_time - round(imu_time)) < Thres_interval&&TCfusion.get_ins_align_flag()
				&& fabs(imu_time-TCfusion.get_fusion_time())> Thres_interval + 0.5) {
				// Rinex Tokyo rover 4:05:50 û���ǣ���Ӧû��ʱ���
				//��ʱ��GNSS_timeӦ���Ƕ����� Ӧ�е�����s֮��ĵ�ĳ��ʱ��
				//ͬʱIMUȡֵ�Ѿ���������s ��ʱӦ��������
				//ins2sol();
				TCfusion.set_ngnsslock(0);
				TCfusion.set_Solstate(3);
				TCfusion.TC_instofile(imu_time,out);
				if (TC_conf.plot_mode) {
					TCfusion.Fusion_plot();
				}
			}
			else if (fabs(imu_time - GNSS_time) < Thres_interval) {//����IMU��GNSS�غϵ�ʱ��
				/*if (fabs(GNSS_time - 1545192630) < Thres_interval)
					int aaaaa = 1;*/
				
				if (!TCfusion.get_ins_align_flag()) {
					TCfusion.ins_align(out);//ȡ��realign�ܹ���ins_alignֻ��gnssԪ���� ��ֻҪ�������˾Ͳ����ٽ���
					if (TCfusion.get_ins_align_flag()) {
						TCfusion.Clear_rcv_array();
						TCfusion.TC_instofile(GNSS_time, out);
						if (TC_conf.plot_mode) {
							TCfusion.Fusion_plot();
						}
					}
					Fusion_Epoch_Update = 0;//�������ݶ�׼���Ӧ�������ٶ�ȡһ��GNSS����,��ʹ��׼�ˣ�Ҳ����һ��Ԫ�ٿ�ʼ���
					continue;
				}
				/*if (imu_time > GNSS_time) {
					cout << "hello" << endl;
				}*/
				if (numiterator >= debug_stop-1) {
					int aaaa = 1;
				}
				TCfusion.set_fusion_time(GNSS_time);
				TCfusion.rtk_ins();//ֻ�����۲���£�����ҲӦ���ϣ����ins_alignû��ɣ�ֱ�ӷ���
				TCfusion.TC_instofile(GNSS_time, out);
				if (TC_conf.plot_mode) {
					TCfusion.Fusion_plot();
				}
				Fusion_Epoch_Update = 0;//������֮����߳�ʼ��׼���Ӧ�������ٶ�ȡһ��GNSS����	
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