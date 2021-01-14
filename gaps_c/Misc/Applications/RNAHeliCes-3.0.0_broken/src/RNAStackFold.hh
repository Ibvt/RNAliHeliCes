#ifndef _RNASTACKFOLD_HH_
#define _RNASTACKFOLD_HH_
#include "Turner.hh"

const std::string LEfile = "TurnerModel/LoopE.dat";
const std::string SEfile = "TurnerModel/StackingE.dat";
const std::string TerminalSMHEfile = "TurnerModel/TerminalStackingMismatchHairpinE.dat";
const std::string TerminalSMIEfile = "TurnerModel/TerminalStackingMismatchInteriorE.dat";
const std::string SingleMIEfile = "TurnerModel/SingleMismatchInteriorE.dat";
const std::string OneTwoMIEfile = "TurnerModel/OneTwoMismatchInteriorE.dat";
const std::string TetraLoopEfile = "TurnerModel/TetraLoopE.dat";
const std::string TwoTwoMIEfile = "TurnerModel/TwoTwoMismatchInteriorE.dat";
const std::string d3Efile = "TurnerModel/d3E.dat";
const std::string d5Efile = "TurnerModel/d5E.dat";

const int OOE = 100000; //const double OOENERGY = 1000.00;
const int OUTWIDTH = 8;

//Energy types 
enum ETYPE{
	etype_unknown = 0,
	etype_F = 1,
	etype_C = 2,
	etype_FM = 3,
	etype_M1 = 4,
	etype_M2 = 5,
	etype_M4 = 6,
	etype_Fn = 7,
	etype_Fd = 8,
	etype_F2n = 9,
	etype_F2d = 10
};
//Case types
enum CTYPE {
	ctype_unknown = -1,
	ctype_0 = 0,
	ctype_1 = 1,
	ctype_2 = 2,
	ctype_3 = 3,
	ctype_4 = 4,
	ctype_1a = 5,
	ctype_1b = 6,
	ctype_1c = 7,
	ctype_2a = 8,
	ctype_2b = 9,
	ctype_2c = 10,	
	ctype_3a = 11,
	ctype_3b = 12,
	ctype_3c = 13,
	ctype_3d = 14,
	ctype_3e = 15,
	ctype_3f = 16,
	ctype_3g = 17,
	ctype_3h = 18,
	ctype_3i = 19,
	ctype_4a = 20,
	ctype_4b = 21,
	ctype_5a = 22,
	ctype_5b = 23
};

//define classes
class PartStructT:public PartStruct{
	public:
	std::vector<ETYPE> sgt;//function type F, C, FM
	int elphi ;//energy of known structures;
	
	std::string log;
	PartStructT(){elphi=0; log="";}
	PartStructT(const seglist & in_sgl, const std::vector<ETYPE> in_sgt, const stemlist & in_stl, const int & in_e=0, const std::string & in_s="");
	void clear();
	PartStructT(const PartStructT & phi);
	bool isLocOpt(const stemlist & in_stl, const std::vector<std::vector<int> > & H, const std::vector<std::vector<int> > & C);
	void show(std::ostream & fout=std::cout);
};

class StackFoldAuxiliary{
	public:
	CTYPE ctype;
	stem st;
	
	StackFoldAuxiliary();
	StackFoldAuxiliary(const CTYPE & in_ctype, const stem & in_t);
	StackFoldAuxiliary(const StackFoldAuxiliary & in_sfa);
//	~StackFoldAuxiliary();
	void set(const CTYPE & in_ctype, const stem & t);
	void show(std::ostream & fout = std::cout);
};


std::string dna2rna(const std::string & seq);
double max(const double & d1, const double & d2);
double min(const double & d1, const double & d2);
int max(const int & d1, const int & d2);
int min(const int & d1, const int & d2);

void init_vi(std::vector<int> & vi, const int & n, const int & val);
void init_vvi(std::vector<std::vector<int> > & vvi, const int & n, const int & m, const int & val);
void init_vstem(std::vector<stem> & vstem, const int & n, stem & val);
void init_vvstem(std::vector<std::vector<stem> > & vvstem, const int & n, const int & m, stem & val);
void init_vsfa(std::vector<StackFoldAuxiliary> & vsfa, const int & n, const StackFoldAuxiliary & val);
void init_vvsfa(std::vector<std::vector<StackFoldAuxiliary> > & vvsfa, const int & n, const int & m, const StackFoldAuxiliary & val);


void showETYPE(const ETYPE & etype, std::ostream & fout = std::cout);
void showCTYPE(const CTYPE & ctype, std::ostream & fout = std::cout);
void show(const std::vector<int>    & vi, std::ostream & fout = std::cout);
void show(const std::vector<ETYPE>  & vet,std::ostream & fout = std::cout);
void show(const std::vector<double> & vd, std::ostream & fout = std::cout);
void show(const std::vector<std::string> & vs, std::ostream & fout = std::cout);
void show(const std::vector<std::vector<int> >    & vvi, std::ostream & fout = std::cout);
void print2file(const std::vector<std::vector<int> > & vvi, const std::string & file);
void show(const std::vector<std::vector<double> > & vvd, std::ostream & fout = std::cout);
void show(const std::vector<std::vector<std::string> > & vvs, std::ostream & fout = std::cout);
void showd(TurnerE & TE, const std::string & seq, const int & a, const int & b, const int & c, std::ostream & fout=std::cout);
void showmind(TurnerE & TE, const std::string & seq, const int & a, const int & b, const int & c, const int & e, const int & f, const int & g, std::ostream & fout=std::cout);
void showH(TurnerE & TE, const std::string & seq, stem & st, std::ostream & fout=std::cout);
void showS(TurnerE & TE, const std::string & seq, stem & st, std::ostream & fout=std::cout);
void showphi(TurnerE & TE, const std::string & seq, stem & st, stem & stk, std::ostream & fout=std::cout);
void showMphi(TurnerE & TE, std::ostream & fout=std::cout);
void showMc(TurnerE & TE, std::ostream & fout=std::cout);
void showTP(TurnerE & TE, const std::string & seq, const int & x, const int & y, std::ostream & fout=std::cout);



//int TerminalPenalty(const char & x, const char & y);
int TerminalPenalty(TurnerE & TE, const char & x, const char & y);
//int d(const std::string & seq, const int & i1, const int & i2, const int & i3);
int d(TurnerE & TE, const std::string & seq, const int & i1, const int & i2, const int & i3);
int d5(TurnerE & TE, const std::string & seq, const int & i1, const int & i2, const int & i3);
int d3(TurnerE & TE, const std::string & seq, const int & i1, const int & i2, const int & i3);
int computeEnergy(TurnerE & TE, const std::string & seq, RStack & in_stl); //given a secondary structure, predict its energy.
int computeH(TurnerE & TE, const std::string & seq, stem & st);
int computeBulge(TurnerE & TE, const std::string & seq, const int & x, const int & y, const int & p, const int & q);
int computeInterior(TurnerE & TE, const std::string & seq, const int & x, const int & y, const int & p, const int & q);
int computeS(TurnerE & TE, const std::string & seq, const int & x, const int & y, const int & p, const int & q);
int computeS(TurnerE & TE, const std::string & seq, stem & s, stem & t);
int computeS(TurnerE & TE, const std::string & seq, stem & st);
int computephi(TurnerE & TE, const std::string & seq, stem & s, stem & t);
stem getUsedPartOfStem(stem & s, stem & t);
stem getUsedPartOfStem(stemlist & in_stl, const int & indx_s, const int & indx_t); //s.enclose(t) or t.inward_extend(s), return the used part of S


//functions for folding and generating all possible partial structures
//int fold(const std::string & seq); //Given putative stacsk in an RNA sequence , fold the minimum base pair structure
//int fold(const std::string & seq, stemlist & in_stl); //Given putative stacsk in an RNA sequence , fold the minimum base pair structure
int generate_all_RStacks(const std::string & seq, std::ostream & fout, const int delta_val);//Given an RNA sequence, the output stream, the delta val, generate all possible structures that are within (global maximum - deltaval).

void generate_all_local_RStacks(const std::string & seq, std::ostream & fout, const int delta_val, std::vector<RStack> & vRStack, int & num_struct, int & num_local_struct);//generate all possible local mimimum structures
int generate_all_local_RStacks(const std::string & seq, std::ostream & fout, const int delta_val);//generate all possible local minimum structures 
int calscore(PartStruct pst) ;//Given a parital structure, return the maximum score possible


class StackFold{
	public:
	std::string seq;
	stemlist stl;
	std::vector<std::vector<int> > H; //H[i][j] = energy of a segment [i,j] provided that there exists a hairpin stem with outer=(i,j) in the segment
	std::vector<std::vector<int> > S; //S[i][j] = stacking energy of s=stl.vstems[i] followed by t=stl.vstems[j], where s.enclose(t) or t.inward_extend(s)
	std::vector<std::vector<int> > phi;//phi[i][j] = energy of interior loop between s=stl.vstem[i] and t=stl.vstem[j], where s.enclose(t) or t.inward_extend(s)

	std::vector<int> F; //F[j] = energy of segment [0,j]
	std::vector<int> Fn;
	std::vector<int> Fd;
	std::vector<int> F2n;
	std::vector<int> F2d;
	
	std::vector<std::vector<int> > C; //C[i][j] = energy of a segment [i,j] provided that the segment is enclosed by a stem with outer=(i,j)
	std::vector<std::vector<int> > FM; //FM[i][j] = energy of segment [i,j], which is a part of a multi-loop and it contains at least one stem
	std::vector<std::vector<int> > M1; //M1[i][j] = energy of segment [i,j], which is a part of a multi-loop and it contains at least one stem starting at i
	std::vector<std::vector<int> > M2; //M2[i][j] = energy of segment [i,j], which is a part of a multi-loop and it contains at least one stem ending at j
	std::vector<std::vector<int> > M4; //M4[i][j] = energy of segment [i,j], which is a part of a multi-loop and it contains at least one stem starting at i and one stem ending at j (the stem starting at i and the stem ending at j can be identical)

	std::vector<StackFoldAuxiliary> F_sfa; //Auxiliary
	std::vector<StackFoldAuxiliary> Fn_sfa; //Auxiliary  
	std::vector<StackFoldAuxiliary> Fd_sfa; //Auxiliary  
	std::vector<StackFoldAuxiliary> F2n_sfa; //Auxiliary  
	std::vector<StackFoldAuxiliary> F2d_sfa; //Auxiliary  
	
	std::vector<std::vector<StackFoldAuxiliary> > C_sfa;
	std::vector<std::vector<StackFoldAuxiliary> > FM_sfa;
	std::vector<std::vector<StackFoldAuxiliary> > M1_sfa;
	std::vector<std::vector<StackFoldAuxiliary> > M2_sfa;
	std::vector<std::vector<StackFoldAuxiliary> > M4_sfa;

	std::vector<stem> F_S; //the right most stem in F_S
	std::vector<stem> Fn_S;
	std::vector<stem> Fd_S;
	std::vector<stem> F2n_S;
	std::vector<stem> F2d_S;
		
	std::vector<std::vector<stem> > FM_LS; //the left most stem in FM[i,j]
	std::vector<std::vector<stem> > FM_RS; //the right most stem in FM[i,j]
	std::vector<std::vector<stem> > M1_LS; //the left most stem in M1[i,j]
	std::vector<std::vector<stem> > M1_RS; //the right most stem in M1[i,j]	
	std::vector<std::vector<stem> > M2_LS; //the left most stem in M2[i,j]
	std::vector<std::vector<stem> > M2_RS; //the righStackFoldt most stem in M2[i,j]	
	std::vector<std::vector<stem> > M4_LS; //the left most stem in M4[i,j]
	std::vector<std::vector<stem> > M4_RS; //the right most stem in M4[i,j]
	
	StackFold(){}
	StackFold(const std::string & in_seq);
	
	int fold(TurnerE & TE);
	std::string traceback(TurnerE & TE);
	//generate all possible secondary structures with energy <= Emin+delta;
	//std::vector<std::string> generate_subopt(std::ostream & fout, const int & e_ub, bool & reqLocOpt, TurnerE & TE);	
   std::vector<RNA> generate_subopt(std::ostream & fout, const int & e_ub, bool & reqLocOpt, TurnerE & TE);
   
private:	
	void init(const int & n, const int & m); // n is the length of the sequence, m is the number of stems in stemlist
	void fillH (TurnerE & TE);
	void fillS (TurnerE & TE);
	void fillphi(TurnerE & TE);
	void fillC (TurnerE & TE, const int & i, const int & j);
	void fillFM(TurnerE & TE, const int & i, const int & j);
	void fillM1(TurnerE & TE, const int & i, const int & j);
	void fillM2(TurnerE & TE, const int & i, const int & j);
	void fillM4(TurnerE & TE, const int & i, const int & j);
	void fillF (TurnerE & TE, const int & j);
	void fillFnFd(TurnerE & TE, const int & j);
	void fillF2nF2d(TurnerE & TE, const int & j);
	
	std::vector<stem> tracebackF(const int & j, TurnerE & TE);
	std::vector<stem> tracebackFn(const int & j, TurnerE & TE);
	std::vector<stem> tracebackFd(const int & j, TurnerE & TE);
	std::vector<stem> tracebackFnFd(const int & j, TurnerE & TE, StackFoldAuxiliary & sfa);
	std::vector<stem> tracebackF2n(const int & j, TurnerE & TE);
	std::vector<stem> tracebackF2d(const int & j, TurnerE & TE);
	std::vector<stem> tracebackF2nF2d(const int & j, TurnerE & TE, StackFoldAuxiliary & sfa);
	
	std::vector<stem> tracebackC(const int & i, const int & j, TurnerE & TE);
	std::vector<stem> tracebackFM(const int & i, const int & j, TurnerE & TE);
	std::vector<stem> tracebackM1(const int & i, const int & j, TurnerE & TE);
	std::vector<stem> tracebackM2(const int & i, const int & j, TurnerE & TE);
	std::vector<stem> tracebackM4(const int & i, const int & j, TurnerE & TE);
	
	
	int getE(const ETYPE & etype, const int & i=-1, const int & j=-1);
	stem getStem(const std::vector<std::vector<stem> > & vvstem, const int & i, const int & j);
	stem getStem(const std::vector<stem> & vstem, const int & j);
	int getsegment_esum(seglist & in_sgl, std::vector<ETYPE> in_sgt);
	

   bool suboptF (std::vector<PartStructT> & R, const segment & in_sg, const PartStructT & in_phi, const int & energy_ub, const bool & reqLocOpt, TurnerE & TE);
   bool suboptC (std::vector<PartStructT> & R, const segment & in_sg, const PartStructT & in_phi, const int & energy_ub, const bool & reqLocOpt, TurnerE & TE);
   bool suboptFM(std::vector<PartStructT> & R, const segment & in_sg, const PartStructT & in_phi, const int & energy_ub, const bool & reqLocOpt, TurnerE & TE);
//   bool subopt_backtrackM (std::vector<PartStructT> & R, const segment & in_sg, const PartStructT & phi, const int & energy_ub);
   int get_esum(const seglist & in_sgl, const std::vector<ETYPE> in_sgt, stemlist & in_stl) ;//return the sum of minimum energy of segments in in_sgl
};

#endif
