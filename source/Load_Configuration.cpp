#include "Load_Configuration.h"


int Load_Config_Params(NAV_RINEX & nav, OBS_RINEX &obs,ifstream&obsfileB, ifstream& obsfileR,string & imufilename) {
	/*This function aims to load the configuration parameters from a file*/
	std::string navfile, Bobsfile, Robsfile;
	std::string confdir = "../conf/";
	std::string conffile = confdir+"TC_conf.ini";
	std::ifstream conf_stream(conffile);
	if (!conf_stream.is_open()) {
		std::cout << "Fail to open the Configuration file" << std::endl;
		return 0;
	}
	std::string line;
	while (getline(conf_stream, line)){
		if (line[0] == '#')
			continue;
		std::stringstream iss(line);
		std::string line_part[2];
		int i = 0;
		while (i<2&&getline(iss, line_part[i], '=')) { // conf 必须以= 分割，但不一定有空格，先分割成两部分
			if ( line_part[i].empty()) {
				std::cout << "Configuration file's format is wrong" << std::endl;
				return 0;
			}
			if (i == 0) 
				line_part[0].erase(line_part[0].find_last_not_of(" ") + 1);
			else {
				line_part[1].erase(0, line_part[1].find_first_not_of(" "));
				if(line_part[1].find("%")!= line_part[1].npos) // 可能会有注释
					line_part[1].erase(line_part[1].find("%"));	
				line_part[1].erase(line_part[1].find_last_not_of(" ") + 1);
			}	
			i++;
		}
		if (line_part[0] == "data_path") {
			navfile = line_part[1] + navfile;
			Bobsfile = line_part[1] + Bobsfile;
			Robsfile = line_part[1] + Robsfile;
			imufilename = line_part[1] + imufilename;
		}
		else if (line_part[0] == "data_name") {
			navfile += (line_part[1] + "\\base.nav");
			Bobsfile += line_part[1] + "\\base.obs";
			Robsfile += line_part[1] + "\\rover.obs";
			imufilename+= line_part[1] + "\\imu.txt";
		}
		else if (line_part[0] == "dynamics_on") {
			TC_conf.dynamics_on = stoi(line_part[1].c_str()) ? 1 : 0;
		}
		else if (line_part[0] == "pos_mode") {
			TC_conf.pos_mode = stoi(line_part[1]);
		}
		else if (line_part[0] == "arfilter") {
			TC_conf.arfilter = stoi(line_part[1].c_str()) ? 1 : 0;
		}
		else if (line_part[0] == "outsolfile_name") {
			TC_conf.outsolfile_name = line_part[1];
		}
		else if (line_part[0] == "base_position") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss,temp_s,',')) {
				TC_conf.basepose[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "navsys"|| line_part[0] == "nfreq") {
			std::cout << "I am sorry I can only support GPS/BDS single frequency now." << endl;
			std::cout << "Now the configuration option has no influence on me." << endl;
			std::cout << "I wish I could become as wondeful as RTKLIB and GINAV one day." << endl;
		}
		else if (line_part[0] == "elmin") {
			TC_conf.elmin_mask = stod(line_part[1])*PI / 180;
		}
		else if (line_part[0] == "Error_std") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.Error_std[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "prn") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.prn[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "Ini_state_std") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.Ini_state_std[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "out_mode") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.out_mode[j++] = stoi(temp_s);
			}
		}
		else if (line_part[0] == "thresar") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.thresar[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "eratio") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.eratio[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "maxinno") {
			TC_conf.maxinno = stod(line_part[1]);
		}
		else if (line_part[0] == "maxout") {
			TC_conf.maxout = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "minlock") {
			TC_conf.minlock = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "plot_mode") {
			TC_conf.plot_mode = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "MAX_GNSS_OUTAGE") {
			TC_conf.MAX_GNSS_OUTAGE = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "IGGC0") {
			TC_conf.IGGC0 = stoi(line_part[1]);
		}
		else if (line_part[0] == "IGGC1") {
			TC_conf.IGGC1 = stoi(line_part[1]);
		}
		else if (line_part[0] == "minfix") {
			TC_conf.minfix = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "optminfixsats") {
			TC_conf.optminfixsats = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "optminholdsats") {
			TC_conf.optminholdsats = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "timef") {
		TC_conf.timef = stoi(line_part[1]);
		}
		else if (line_part[0] == "optmindropsats") {
			TC_conf.optmindropsats = (short)stoi(line_part[1]);
		}
		else if (line_part[0] == "sclkstab") {
			TC_conf.sclkstab = stod(line_part[1]);
		}
		else if (line_part[0] == "varholdamb") {
			TC_conf.varholdamb = stod(line_part[1]);
		}
		else if (line_part[0] == "Ls_res_threshold") {
			TC_conf.Ls_res_threshold = stod(line_part[1]);
		}
		else if (line_part[0] == "var_pos_std") {
			TC_conf.VAR_POS = pow(stod(line_part[1]),2);
		}
		else if (line_part[0] == "var_vel_std") {
			TC_conf.VAR_VEL = pow(stod(line_part[1]), 2);
		}
		else if (line_part[0] == "var_acc_std") {
			TC_conf.VAR_ACC = pow(stod(line_part[1]), 2);
		}
		else if (line_part[0] == "posf") {
			TC_conf.output_pos_format = line_part[1];
		}
		else if (line_part[0] == "t_interval") {
			TC_conf.tint = stod(line_part[1]);
		}
		else if (line_part[0] == "ins_mode") {
			TC_conf.ins_mode = (unsigned char)stoi(line_part[1]);
			if (!TC_conf.ins_mode)
				TC_conf.pos_mode = 0;
		} 
		else if (line_part[0] == "att_mode") {
			TC_conf.att_mode = stoi(line_part[1]);
		}
		else if (line_part[0] == "ins_aid_IGG") {
			TC_conf.ins_aid_IGG = (unsigned char)stoi(line_part[1]);
		}
		else if (line_part[0] == "sample_rate") {
			TC_conf.sample_rate = (float)stod(line_part[1]);
		}
		else if (line_part[0] == "armode") {
			TC_conf.armode = stoi(line_part[1]); 
		}
		else if (line_part[0] == "exec_mode") {
			TC_conf.exec_mode = stoi(line_part[1]);
		}
		else if (line_part[0] == "lever") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.lever[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "init_att_unc") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.init_att_unc[j++] = stod(temp_s)/180*PI;
			}
		}
		else if (line_part[0] == "init_vel_unc") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.init_vel_unc[j++] = stod(temp_s);
			}
		}
		else if (line_part[0] == "init_pos_unc") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.init_pos_unc[j++] = stod(temp_s);
			}
			TC_conf.init_pos_unc[0] /= 6378137.0;
			TC_conf.init_pos_unc[1] /= 6378137.0;
		}
		else if (line_part[0] == "init_bg_unc") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.init_bg_unc[j++] = stod(temp_s);
			}
			if (j == 1) {
				TC_conf.init_bg_unc[1] = TC_conf.init_bg_unc[0];
				TC_conf.init_bg_unc[2] = TC_conf.init_bg_unc[0];
			}
		}
		else if (line_part[0] == "init_ba_unc") {
			std::stringstream temp_ss(line_part[1]);
			std::string temp_s;
			int j = 0;
			while (getline(temp_ss, temp_s, ',')) {
				TC_conf.init_ba_unc[j++] = stod(temp_s);
			}
			if (j == 1) {
				TC_conf.init_ba_unc[1] = TC_conf.init_ba_unc[0];
				TC_conf.init_ba_unc[2] = TC_conf.init_ba_unc[0];
			}
		}
		else if (line_part[0] == "psd_gyro") {
			TC_conf.psd_gyro = stod(line_part[1]);
		}
		else if (line_part[0] == "psd_acce") {
			TC_conf.psd_acce = stod(line_part[1]);
		}
		else if (line_part[0] == "psd_bg") {
			TC_conf.psd_bg = stod(line_part[1]);
		}
		else if (line_part[0] == "psd_ba") {
			TC_conf.psd_ba = stod(line_part[1]);
		}
	}
	conf_stream.close();
	nav.setfilename(navfile);
	if (!nav.readfile()) {
		return 0;
	}
	obs.setfilesnames(Bobsfile, Robsfile);
	obsfileB.open(Bobsfile);
	if (!obsfileB.is_open()) {
		std::cout << "Fail to open the base obs file" << endl;
		return 0;
	}
	obsfileR.open(Robsfile);
	if (!obsfileR.is_open()) {
		std::cout << "Fail to open the rover obs file" << endl;
		return 0;
	}
	return 1;
}

