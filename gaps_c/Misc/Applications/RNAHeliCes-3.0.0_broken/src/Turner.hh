#ifndef _TURNER_HH_
#define _TURNER_HH_
#include <map>
#include "RNAStack.hh"
#include <math.h>

#define MAXLOOPLEN 30 //the length of the longest loop in a bulge, an interior \
loop or a hairpin loop
//#define MAXINTERIORLENDIFF 10//the difference between lengths of two loops in  \
an interior loop

class TurnerE{
public:
	std::map<int, int> INTERNAL;//Destablizing energy of interior loops
	std::map<int, int> BULGE;//Destablizing energy of bulges
	std::map<int, int> HAIRPIN;//Destablizing energy of hairpin loops
	
	std::map<std::string, int> STACKING;//stacking energy
	std::map<std::string, int> TerminalSMHE;//terminal stacking mismatch hairpin energy
	std::map<std::string, int> TerminalSMIE;//terminal stacking mismatch interior energy
	std::map<std::string, int> SingleMIE;//single mismatch interior energy
	std::map<std::string, int> OneTwoMIE;//1x2 mismatch interior energy
	std::map<std::string, int> TetraLoopE;//specially stable tetra loop energy;
	std::map<std::string, int> TwoTwoMIE;//2x2 mismatch interior energy
	
	std::map<std::string, int> d3E;//dangling single base stacking energy at 3' end
	std::map<std::string, int> d5E;//dangling single base stacking energy at 5' end
	
	int Mh;//energy for every hydrogen bond
	int Mb;//penalty for an unpaired nt in a dangling region
	int Mc;//penalty for opening a multiloop
	int Mphi;//penalty for having every branch in a multiloop
	int Mbb;//panlaty for an unpaired nt in a multiloop
	int Mninio_m;//penalty for difference of lengths between two gaps of an interior loop
	int Mninio_max; //max penalty for difference of lengths between two gaps of an interior loop
	int MtermianAU;//penalty for having AU in the terminal;
	double Mlongloop; //For internal/bulge/hairpin loops > 30, dS(T) = ds(30)+Mlongloop*ln(n/30), Mlongloop = 1.079
	int Mgubous; //For GU enclosure preceded by GG
	int Mcslope; //c hairpin slope 
	int Mcintercept; //c hairpin intercept 
	int Mc3; //c hairpin of 3
	
	TurnerE(){}
	TurnerE(const std::string & LEfile, const std::string & SEfile, const std::string & 
	  TerminalSMHEfile, const std::string & TerminalSMIEfile, const std::string & 
	  SingleMIEfile, const std::string & OneTwoMIEfile, const std::string & TetraLoopEfile,
	  const std::string & TwoTwoMIEfile, const std::string & d3Efile, const std::string & d5Efile)
	{
		init(LEfile, SEfile, TerminalSMHEfile, TerminalSMIEfile, SingleMIEfile, OneTwoMIEfile, TetraLoopEfile, TwoTwoMIEfile, d3Efile, d5Efile);
	}
	void init(const std::string & LEfile, const std::string & SEfile, const std::string & 
	  TerminalSMHEfile, const std::string & TerminalSMIEfile, const std::string & 
	  SingleMIEfile, const std::string & OneTwoMIEfile, const std::string & TetraLoopEfile,
	  const std::string & TwoTwoMIEfile, const std::string & d3Efile, const std::string & d5Efile);
	void readfiletype(const std::string & file, std::map <std::string, int> & table, 
	  const int & pos1, const int & pos2);
	void init_LE(const std::string & LEfile);
	void init_SE(const std::string & SEfile);
	void init_TerminalSMHE(const std::string & TerminalSMHEfile);
	void init_TerminalSMIE(const std::string & TerminalSMIEfile);
	void init_SingleMIE(const std::string & SingleMIEfile);
	void init_OneTwoMIE(const std::string & OneTwoMIEfile);
	void init_TetraLoopE(const std::string& TetraLoopEfile);
	void init_TwoTwoMIE(const std::string & TwoTwoMIEfile);
	void init_d3E(const std::string & d3Efile);
	void init_d5E(const std::string & d5Efile);
};

#endif

