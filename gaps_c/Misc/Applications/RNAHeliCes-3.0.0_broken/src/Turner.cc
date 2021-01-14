#include "Turner.hh"

std::string trimspace(const std::string & str) {
//	cout<<"input std::string ["<<str<<"]"<<std::endl;
	int s = 0,  e = sz(str)-1; 
	bool stop = false;
	FOR(i,0,sz(str)) {
		if(str[i] == ' ' || str[i] == '\t' ||  str[i] == '\n') {
		} else {
			s = i;
			stop = true;
			break;
		}
	}
	if(!stop) return "";
	
	for(int i= sz(str)-1; i>=s; i--) {
		if(str[i] == ' ' || str[i] == '\t' ||  str[i] == '\n') {
		} else {
			e = i;
			break;
		}
	}
	std::string ret = str.substr(s, (e-s+1));
//	cout<<"output std::string ["<<ret<<"]"<<std::endl;
	return ret;
}

void TurnerE::init(const std::string & LEfile, const std::string & SEfile, const std::string & 
  TerminalSMHEfile, const std::string & TerminalSMIEfile, const std::string & 
  SingleMIEfile, const std::string & OneTwoMIEfile, const std::string & TetraLoopEfile, 
  const std::string & TwoTwoMIEfile, const std::string & d3Efile, const std::string & d5Efile) {
	init_LE(LEfile);
	init_SE(SEfile);
	init_TerminalSMHE(TerminalSMHEfile);
	init_TerminalSMIE(TerminalSMIEfile);
	init_SingleMIE(SingleMIEfile);
	init_OneTwoMIE(OneTwoMIEfile);
	init_TetraLoopE(TetraLoopEfile);
	init_TwoTwoMIE(TwoTwoMIEfile);
	init_d3E(d3Efile);
	init_d5E(d5Efile);
}

void TurnerE::init_LE(const std::string & file)
{
//	cout<<"LoopEfile="<<file<<std::endl;
//format: #Size INTERNAL BULGE	HAIRPIN                                       \
1	.	3.80	.                                                                  \
2	.	2.80	.                                                                  \
3	.	3.20	5.70                                                               \
4	1.70	3.60	5.60                                                            \
5	1.80	4.00	5.60                                                            \
6	2.00	4.40	5.40                                                            \
7	2.20	4.60	5.90

	std::ifstream fin(file.c_str());
	while(fin.good() && !fin.eof()){
		char buf[MAXCHARNUM];
		fin.getline(buf, MAXCHARNUM);
		std::string str(buf);
		str = trimspace(str);
		if(str[0] == '#' || sz(str)==0) {
			//comments, don't process
			continue;
		}
		std::stringstream ss(str); 
		std::string s1, s2, s3, s4;
		ss>>s1>>s2>>s3>>s4;
//		cout<<s1<<"	"<<s2<<"	"<<s3<<"	"<<s4<<"	"<<std::endl;
		int i1;
		int i2, i3, i4;
		std::stringstream ss1(s1); ss1>>i1;
		if(s2.compare(".") != 0) {
			i2 = atoi(s2.c_str());
//			cout<<d2<<" ";
			INTERNAL[i1] = i2;
		}
		if(s3.compare(".") != 0) {
			i3 = atoi(s3.c_str());
//			cout<<d3<<" ";
			BULGE[i1] = i3;
		}
		if(s4.compare(".") != 0) {
			i4 = atoi(s4.c_str());
//			cout<<d4<<std::endl;
			HAIRPIN[i1] = i4;
		}
	}
	fin.close();
		
	Mh = -100;//energy for every hydrogen bond
	Mb = 10;//penalty for an unpaired nt in a dangling region
	Mc = 340;//penalty for opening a multiloop
	Mphi = 40;//penalty for every branch in a multiloop
	Mbb = 0; //penalty for every unpaired nt in multiploop
	Mninio_m = 50;//penalty for difference of lengths between two gaps of an interior loop
	Mninio_max = 300; //max penalty for difference of lengths between two gaps of an interior loop
	MtermianAU = 50;//penalty for having AU in the terminal
	Mlongloop = 107.9;//For internal/bulge/hairpin loops > 30, dS(T) = ds(30)+Mlongloop*ln(n/30), Mlongloop = 1.079
	Mgubous = -220;
	Mcslope = 30;
	Mcintercept = 160;
	Mc3 = 140;
	/*
	std::map<int, int>::iterator it;
	cout<<"INTERNAL"<<std::endl;
	for(it=INTERNAL.begin(); it!=INTERNAL.end(); it++){
		cout<<(*it).first<<" "<<(*it).second<<std::endl;
		if((*it).second<10 && (*it).second>-10 && (*it).second!=0) {
			std::cerr<<"error!"<<std::endl;
		}
	}
	
	cout<<"BULGE"<<std::endl;
	for(it=BULGE.begin(); it!=BULGE.end(); it++){
		cout<<(*it).first<<" "<<(*it).second<<std::endl;
		if((*it).second<10 && (*it).second>-10 && (*it).second!=0) {
			std::cerr<<"error!"<<std::endl;
		}
	}
	
	cout<<"HAIRPIN"<<std::endl;
	for(it=HAIRPIN.begin(); it!=HAIRPIN.end(); it++){
		cout<<(*it).first<<" "<<(*it).second<<std::endl;
		if((*it).second<10 && (*it).second>-10 && (*it).second!=0) {
			std::cerr<<"error!"<<std::endl;
		}
	}*/
}

void TurnerE::readfiletype(const std::string & file, std::map <std::string, int> & table, 
  const int & pos1, const int & pos2) {
//format:                                                                      \
#         UX                                                                   \
#         GY                                                                   \
>UXYG                                                                          \
  .     .     .   -100                                                         \
  .     .   -150   .                                                           \
  .   -140   .    30                                                           \
-60     .   -50   .                                                            \
assure name[pos1]=='X' and name[pos2]=='Y'
	char alphabet[4];
	bool setalphabet = false;
	std::ifstream fin(file.c_str());
	while(fin.good() && !fin.eof()){
		char buf[MAXCHARNUM];
		fin.getline(buf, MAXCHARNUM);
		std::string str(buf);
		str =  trimspace(str);
//		cout<<"readline: ["<<str<<"]"<<std::endl;
		if(sz(str)==0 || str[0] == '#') {
			//comments, don't process
			continue;
		}
		std::stringstream ss(str);
		std::string s1, s2;
		if(str.find("alphabet") != std::string::npos) {
			str = trimspace(str.substr(str.find("=")+1));
			FOR(i,0,4)
				alphabet[i] = str[i];
//			cout<<"alphabet="<<alphabet[0]<<alphabet[1]<<alphabet[2]<<alphabet[3]<<std::endl;
			setalphabet=true;
			continue;
		}
		std::string name;
		std::vector<std::string> vval;
		FOR(i,0,4) vval.PB("");
		
		if(str[0] == '>') {
			if(!setalphabet) {
				std::cerr<<"file:"<<file<<" does not specify alphabet, please add alphabet=ACGU to the first line."<<std::endl;
				exit(0);
			}
			name = str.substr(1);
//			cout<<"name is ["<<name<<"]"<<std::endl;
			if(pos2>=0)	assert(name[pos1]='X' && name[pos2]=='Y');
			else assert(name[pos1]=='X');
			if(pos2>=0) {
				FOR(i,0,4) {			
					//row i
					fin.getline(buf, MAXCHARNUM);
					std::string str(buf);
					str = trimspace(str);
					std::stringstream ss2(str);
					FOR(j,0,4) {
						//column j
						ss2>>vval[j];
		//				cout<<"vval[j]="<<vval[j]<<std::endl;
						if(vval[j].compare(".")!=0) {
							std::string newname = name;
							newname[pos1] = alphabet[i];
							newname[pos2] = alphabet[j];
							std::stringstream ss3(vval[j]);
							int ival;
							ss3>>ival;
		//					cout<<"newname=["<<newname<<"], dval="<<dval<<std::endl;
		//					cout<<"table["<<newname<<"]"<<table[newname]<<std::endl;
							table[newname] = ival;
		//					cout<<newname<<" = "<<dval<<std::endl;
						}
					}
				}
			} else {
				fin.getline(buf, MAXCHARNUM);
				std::string str(buf);
				str = trimspace(str);
				std::stringstream ss2(str);
				FOR(j,0,4) {
					//colum j;
					std::string newname = name;
					newname[pos1] = alphabet[j];
					int ival;
					ss2 >> ival;
					table[newname] = ival;
				}
			}
		}
	}
	/*
	std::map<std::string, int>::iterator it;
	for(it=table.begin(); it!=table.end(); it++){
		cout<<(*it).first<<" "<<(*it).second<<std::endl;
		if((*it).second<10 && (*it).second>-10 && (*it).second!=0) {
			std::cerr<<"error!"<<std::endl;
		}
	}*/
}

void TurnerE::init_SE(const std::string & file) {
//	cout<<"init_SE"<<std::endl;
	readfiletype(file, STACKING, 1, 2);
}
void TurnerE::init_TerminalSMHE(const std::string & file) {
//	cout<<"init_TerminalSMHE"<<std::endl;
	readfiletype(file, TerminalSMHE, 1, 2);
}
void TurnerE::init_TerminalSMIE(const std::string & file) {
//	cout<<"init_TerminalSMIE"<<std::endl;
	readfiletype(file, TerminalSMIE, 1, 2);
}
void TurnerE::init_SingleMIE(const std::string & file) {
//	cout<<"init_SingleMIE"<<std::endl;
	readfiletype(file, SingleMIE, 1, 4);
}
void TurnerE::init_OneTwoMIE(const std::string & file) {
//	cout<<"init_OneTwoMIE"<<std::endl;
	readfiletype(file, OneTwoMIE, 1, 5);
}

void TurnerE::init_TetraLoopE(const std::string & file) {
	std::ifstream fin(file.c_str());
	while(fin.good() && !fin.eof()){
		char buf[MAXCHARNUM];
		fin.getline(buf, MAXCHARNUM);
		std::string str(buf);
		str =  trimspace(str);
//		cout<<"readline: ["<<str<<"]"<<std::endl;
		if(sz(str)==0 || str[0] == '#') {
			//comments, don't process
			continue;
		}
		std::stringstream ss(str);
		std::string s1, s2;
		ss>>s1>>s2;
		s1 = trimspace(s1);
		assert(sz(s1)==6);
		int ival;
		std::stringstream ss2(s2);
		ss2>>ival;
		TetraLoopE[s1] = ival;
	}
	fin.close();
	/*
	std::map<std::string, int>::iterator it;
	for(it=TetraLoopE.begin(); it!=TetraLoopE.end(); it++){
		cout<<(*it).first<<" "<<(*it).second<<std::endl;
		if((*it).second<10 && (*it).second>-10 && (*it).second!=0) {
			std::cerr<<"error!"<<std::endl;
		}
	}*/
}

void TurnerE::init_TwoTwoMIE(const std::string & file) {
	char alphabet[4];//by default, ACGU
	std::vector<std::string> alphabet2;//by default, AA AC AG AU CA CC CG CU GA GC GG GU UA UC UG UU
	bool setalphabet = false;
//format: 
//>AXYAUQPU
	int posX = 1; int posY = 2; int posP = 6; int posQ = 5;
	std::ifstream fin(file.c_str());
	while(fin.good() && !fin.eof()){
		char buf[MAXCHARNUM];
		fin.getline(buf, MAXCHARNUM);
		std::string str(buf);
		str =  trimspace(str);
		if(sz(str)==0 || str[0] == '#') {
			//comments, don't process
			continue;
		}
		std::stringstream ss(str);
		std::string s1, s2;
		if(str.find("alphabet") != std::string::npos) {
			str = trimspace(str.substr(str.find("=")+1));
			FOR(i,0,4)
				alphabet[i] = str[i];
//			cout<<"alphabet="<<alphabet[0]<<alphabet[1]<<alphabet[2]<<alphabet[3]<<std::endl;
//			cout<<"alphabet2=";
			FOR(i,0,4) FOR(j,0,4) {
				std::string tmp = ""; tmp += alphabet[i]; tmp += alphabet[j];
				alphabet2.PB(tmp);
//				cout<<alphabet2[sz(alphabet2)-1]<<" ";
			}
//			cout<<std::endl;
			setalphabet=true;
			continue;
		}

		std::string name;
		std::vector<std::string> vval;
		FOR(i,0,sz(alphabet2)) vval.PB("");
		
		if(str[0] == '>') {
			if(!setalphabet) {
				std::cerr<<"file:"<<file<<" does not specify alphabet, please add alphabet=ACGU to the first line."<<std::endl;
				exit(0);
			}
			name = str.substr(1);
			assert(name[posX]=='X' && name[posY]=='Y' && name[posP] == 'P' && name[posQ] == 'Q');
		}		
		FOR(i,0,sz(alphabet2)) {
			//row i, alphabet2[i] = (XP)
			fin.getline(buf, MAXCHARNUM);
			std::string str(buf);
			str = trimspace(str);
			std::stringstream ss2(str);
			FOR(j,0,sz(alphabet2)) {
				//column j, alphabet2[j] = (YQ)
//				cout<<"row["<<i<<"]"<<alphabet2[i]<<", col["<<j<<"]"<<alphabet2[j]<<std::endl;
				ss2>>vval[j];
//				cout<<"vval[j]="<<vval[j]<<std::endl;
				if(vval[j].compare(".")!=0) {
					std::string newname = name;
					newname[posX] = alphabet2[i][0];
					newname[posP] = alphabet2[i][1];
					newname[posY] = alphabet2[j][0];
					newname[posQ] = alphabet2[j][1];
					std::stringstream ss3(vval[j]);
					int ival; ss3>>ival;
					TwoTwoMIE[newname] = ival;
//					cout<<"newname=["<<newname<<"], ival="<<ival<<std::endl;
//					cout<<"table["<<newname<<"]="<<TwoTwoMIE[newname]<<std::endl;
				}
			}
		}
	}
	/*
	std::map<std::string, int>::iterator it;
	for(it=TwoTwoMIE.begin(); it!=TwoTwoMIE.end(); it++){
		cout<<(*it).first<<" "<<(*it).second<<std::endl;
		if((*it).second<10 && (*it).second>-10 && (*it).second!=0) {
			std::cerr<<"error!"<<std::endl;
		}
	}*/
}


void TurnerE::init_d3E(const std::string & file) {
	readfiletype(file, d3E, 1, -1);
}

void TurnerE::init_d5E(const std::string & file) {
	readfiletype(file, d5E, 1, -1);
}
















