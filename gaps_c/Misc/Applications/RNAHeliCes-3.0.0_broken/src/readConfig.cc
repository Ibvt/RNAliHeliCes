#include "RNAStack.hh"
/***********************************************************
                         Config
***********************************************************/
void Config::setConfig(const std::string & configfile) {
	if(!fexists(configfile.c_str())) {
		std::cerr<<"The specified config file does not exist."<<std::endl;
		exit(1);
	}
	std::ifstream configin(configfile.c_str());
	char buf[MAXCHARNUM];
	while(!configin.eof()){
		configin.getline(buf, MAXCHARNUM);
		std::string sline = std::string(buf);
		if(sline[0] == '#') continue;
		if(sz(sline)==0) continue;
		int indx = sline.find("=");
		if(indx != std::string::npos) {
			std::string name = sline.substr(0,indx);
			std::string val = sline.substr(indx+1);
			name = trim(name);
			val = trim(val);
			std::istringstream iss(val);
			if(name.compare("ALLOWGU")==0){
				if(val.compare("TRUE")==0 || val.compare("True")==0 || val.compare("true")==0) {
					ALLOWGU = true;
				} else if(val.compare("FALSE")==0 || val.compare("False")==0 || val.compare("false")==0) {
					ALLOWGU = false;
				} else {
					std::cerr<<"Can not recognize ALLOWGU = "<<val<<std::endl<<"Please specify TRUE or FALSE in your configuration file (by default: config)."<<std::endl;
					exit(1);
				}
			}else if(name.compare("MINSTACKLEN")==0){
				iss>>MINSTACKLEN;
			}else if(name.compare("MINLOOPLEN")==0){
				iss>>MINLOOPLEN;
			}else if(name.compare("MINHYDROGEN")==0){
				iss>>MINHYDROGEN;
			}else if(name.compare("MAXSEQLEN")==0){
				iss>>MAXSEQLEN;
			}else if(name.compare("MAXMATRIXSZ")==0){
				iss>>MAXMATRIXSZ;
			}else if(name.compare("BARRIERCUTOFF")==0){
				iss>>BARRIERCUTOFF;
			}else if(name.compare("BARRIERCUTOFFDELTA")==0){
				iss>>BARRIERCUTOFFDELTA;
			}else if(name.compare("MAXINWARDOVERLAP")==0) {
				iss>>MAXINWARDOVERLAP;
			}else {
				std::cerr<<"Can not analyze "<<sline<<"in your configuration file:"<<configfile<<std::endl<<"Please update it."<<std::endl;
				exit(1);
			}
		}
	}
}

std::string Config::trim(const std::string & in_s) {
	std::string ret = in_s;
	while(ret[0] == ' ' || ret[0] == '\t'){
		ret = ret.substr(1);
	}
	while(ret[sz(ret)-1] == ' ' || ret[sz(ret)-1] == '\t' 
	   || ret[sz(ret)-1] == '\n') {
		ret = ret.substr(0, sz(ret)-1);
	}
	return ret;
}

void Config::show(std::ostream & fout) {
	fout<<"ALLOWGU = ";
	if(ALLOWGU)  fout<<"TRUE"<<std::endl;
	else fout<<"FALSE"<<std::endl;
	
	fout<<"MINSTACKLEN = "<<MINSTACKLEN<<std::endl;
	fout<<"MINLOOPLEN = "<<MINLOOPLEN<<std::endl;
	fout<<"MINHYDROGEN = "<<MINHYDROGEN<<std::endl;
	fout<<"MAXSEQLEN = "<<MAXSEQLEN<<std::endl;
	fout<<"MAXMATRIXSZ = "<<MAXMATRIXSZ<<std::endl;
	fout<<"BARRIERCUTOFF = "<<BARRIERCUTOFF<<std::endl;
	fout<<"BARRIERCUTOFFDELTA = "<<BARRIERCUTOFFDELTA<<std::endl;
}

void Config::increaseBARRIERCUTOFF() {
	BARRIERCUTOFF += BARRIERCUTOFFDELTA;
}

void Config::resetBARRIERCUTOFF(const int & num ) {
	if(num <= 5000) {
	} else if(num <= 10000){ 
		BARRIERCUTOFF += BARRIERCUTOFFDELTA * 1;
	} else if(num <= 20000) {
		BARRIERCUTOFF += BARRIERCUTOFFDELTA * 2;
	} else if(num <= 30000) {
		BARRIERCUTOFF += BARRIERCUTOFFDELTA * 3;
	} else if(num <= 50000) {
		BARRIERCUTOFF += BARRIERCUTOFFDELTA * 4;
	} else if(num <= 70000) {
		BARRIERCUTOFF += BARRIERCUTOFFDELTA * 5;
	} else {
		BARRIERCUTOFF += BARRIERCUTOFFDELTA * 6;
	}
}
