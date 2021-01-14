#ifndef _RNASTACK_HH_
#define _RNASTACK_HH_

#include <cstdlib>
#include <cstdio>
#include <string>
#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <stack>
#include <assert.h>

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h> //for checking if a file exists
#include <time.h>	    //for rand() and initializing seed
//#include <hash_map>//for hashing table

//NOTE: For TR1 random numbers: #include <tr1/random> and use std::tr1::normal_distribution.
//NOTE: For c++11 random numbers: compile with flag -std=c++0x or -std=gnu++0x, then #include <random> and use std::normal_distribution.
//#include <tr1/unordered_set>
#include <algorithm> //for using abs
#include <iomanip>

//using namespace std;


#define sz(x) (int)x.size()
//##
#define FOR(i,a,b) for(int i=a;i<b;++i)
//##PB==push_back
#define PB push_back
#define MAXCHARNUM 20000
const double PRECISION = 1e-6;
//#define ALLOWGU 1
//#define MAXSEQLEN 500
//#define MINSTACKLEN 3 /* minimum stack length */
//#define MINLOOPLEN 3 /* minimum number nt in a loop; unpaired bases between any pairing bases */
//#define MINHYDROGEN 8


//const std::string RNAevalPath = "/home/li/Software/ViennaRNA-1.8.4/bin/RNAeval -d1 ";//the parameters are exactly the same as PClotes' RNAPathFinder

/* prototype of classes */
class BP;
class Action;
class RStruct;
class stem;
class SD;

/* functions used in common */
bool eq (double d1, double d2);
int myrand(int n); // Random number generator
int getindx(const std::vector<int>& vi, const int& query, const int & startpos=0, const int & direction=1);
int getindx(const std::vector<BP> & in_vbp, const BP & in_bp, const int & startpos=0, const int & direction=1);//return the index of in_bp in in_vbp
int getindx(const std::vector<Action> & in_va, const Action & in_a, const int & startpos=0, const int & direction=1);//return the index of the first in_a in in_va
int getindx(const std::vector<stem> & in_vst, const stem & in_st, const int & startpos=0, const int & direction=1);
std::vector<BP> s2v(const std::string & str);//given a std::string representation of a strucuture, return vbp
std::string v2s(const int & l, const std::vector<BP> & in_vbp);//given indice of the left/right parenthesis and the length of the RStruct, return s;
std::vector<int> getunp(const int & l, const std::vector<BP> & in_vbp);//return indices of unpaired bases
std::vector<int> getunp(const std::string & str);//return indices of unpaired bases
bool fexists(const char *filename);//check if a file exists
bool match(const char & x , const char & y);
char getcomplementchar(const char & x); 
bool descendingsort(double d1, double d2);
bool descendingsortBP(BP b1, BP b2);
bool descendingsortST(stem st1, stem st2);
bool descendingsortSD(SD st1, SD st2);
bool ascendingsort(double d1, double d2);
bool ascendingsortBP(BP b1, BP b2);
bool ascendingsortST(stem st1, stem st2);
bool ascendingsortSD(SD st1, SD st2);
void initialize_fg(int val=0);
bool mycheckclash(const BP & in_bp, const std::vector<BP>& in_vbp);

class Config{
public:
	bool ALLOWGU;
	int MINSTACKLEN ;//3
	int MINLOOPLEN  ;//3
	int MINHYDROGEN ;//8
	int MAXSEQLEN   ;//500
	int MAXMATRIXSZ ;//1000
	double BARRIERCUTOFF; //1
	double BARRIERCUTOFFDELTA; //3
	double MAXINWARDOVERLAP;//0.25
	Config () {
		ALLOWGU = true;
		MINSTACKLEN = 4; /* minimum stack length */
		MINLOOPLEN  = 3; /* minimum number nt in a loop; unpaired bases between any pairing bases */
		MINHYDROGEN = 8;
		MAXSEQLEN = 400;
		MAXMATRIXSZ = 1000; //the maximum size (maximum number of elements in a row) of a n*n/2 double matrix, the maximum number of structures considered.
		BARRIERCUTOFF = 12; //the algorithm reduces the number of structures considered by representing a set of structures that are close to each other by the one with minimum free energy.
		BARRIERCUTOFFDELTA = 1; //if the number of structures to consider is greater than MAXMATRIXSZ, increase DISTANCECUTOFF by DISTANCECUTOFFDELTA, so as to reduce the size of the distance matrix.
		MAXINWARDOVERLAP = 0.20;
	}
	Config(const std::string & configfile) { setConfig(configfile);}
	void setConfig(const std::string & configfile);
	std::string trim(const std::string & in_s) ;
	void show(std::ostream & fout=std::cout);
	void resetBARRIERCUTOFF(const int & num);
	void increaseBARRIERCUTOFF();
};

class BP{//base pair
	public: 
	int lp; //index of the left parenthesis
	int rp; //index of the right parenthesis
	BP() {lp=rp=0;}
	BP(const BP & in_bp) {
		lp = in_bp.lp; rp = in_bp.rp;
	}
	BP(int i1, int i2) {
		lp=i1; rp=i2;
	}
	void show() {//std::cout<<"("<<lp<<","<<rp<<")";}
		std::cout<<"("<<lp<<" "<<rp<<")";
	}
	/*
	bool operator == (const BP & in_bp) {
		return ( (lp==in_bp.lp && rp==in_bp.rp) || (lp==in_bp.rp && rp==in_bp.lp));
	}
	bool operator != (const BP & in_bp) {
		return ( lp != in_bp.lp ||  rp != in_bp.rp);
	}*/
	bool conflictwithBP(const BP & in_bp) {
		if(in_bp.lp <= lp && lp <= in_bp.rp && in_bp.rp <= rp) return true;
		if(lp <= in_bp.lp && in_bp.lp <= rp && rp <= in_bp.rp) return true;
		return false;
	}
};
bool operator == (const BP & in_bp1, const BP & in_bp2);
bool operator != (const BP & in_bp1, const BP & in_bp2);


/* define classes */
enum Act{
	del=1,
	add=2
};

class Action{
	public:	
	Act a;
	BP  b;
	Action(){b=BP();}
	Action(const Action & in_a) {a=in_a.a; b=in_a.b; check();}
	Action(const Act & input1, const int & input2, const int & input3){
		a = input1; b.lp = input2; b.rp = input3; check();
	}
	Action(const Act & input1, const BP & in_bp) {
		a = input1; b = in_bp; check();
	}
	Action(const std::string & input1, const BP & in_bp) {
		if(input1 == "add") a = add;
		if(input1 == "del") a = del; 
		check();
		b = in_bp;
	}
	Action(const std::string & input1, const int & input2, const int & input3) {
		if(input1 == "add") a = add;
		if(input1 == "del") a = del;
		check();
		b.lp = input2; b.rp = input3;
	}
	bool operator == (const Action & in_a) {
		return (a==in_a.a && b.lp==in_a.b.lp && b.rp==in_a.b.rp);
	}
	bool operator != (const Action & in_a) {
		return (a!=in_a.a || b.lp!=in_a.b.lp || b.rp!=in_a.b.rp);
	}
	void show() {
		if(a==del) std::cout<<"del ";
		else if(a==add) std::cout<<"add ";
		check();
		b.show();
	}
	void check() {if(a!=add && a!=del) {std::cerr<<"Action = "<<a<<std::endl;exit(0);}}
};


//## stem
class stem {
	public:
	BP outer; //for short (x=outer.lp, y=outer.rp, p=inner.lp, q=inner.rp)
	BP inner;
	stem(){}
	stem(const int & x, const int & y, const int & p, const int & q) {
		outer.lp = x; outer.rp = y;
		inner.lp = p; inner.rp = q;
	}
	stem(BP bp1, BP bp2) {
		outer = bp1;
		inner = bp2;
	}
	stem(const stem & st) {
		outer = st.outer;
		inner = st.inner;
	}
	
	bool valid() {
		if(outer.lp<0 || outer.rp<0 || 
		   inner.lp<0 || inner.rp<0 || 
		   outer.lp>=outer.rp || inner.lp>=inner.rp || outer.lp>inner.lp)
			return false;
		return true;
	}
	bool enclose(const stem & st) {//check if this stem totally enclose stem st \
i.e. ((((......(((((.....))))).........)))), the outer stem enclose the inner  \
stem
		return (inner.lp < st.outer.lp && inner.rp > st.outer.rp);
	}
	bool enclosedby(const stem & st) {
		return (st.inner.lp < outer.lp && outer.rp < st.inner.rp);
	}
	bool inward_extend(const stem & s);
	int inwhichseg(const int & query);
	int pairwithwhich(const int i) {
		if(i<outer.lp || i>outer.rp || (i>inner.lp && i<inner.rp)) 
			return -1;
		else 
			return (outer.lp+outer.rp-i);
	}
//	int len() {return (outer.rp-outer.lp+1);}
	int width() {return outer.rp-inner.rp+1;}
//	bool ispartof (const stem & st){//check if this stem is part of another stem
//		return ();
//	}
	//check if this stem is conflicting with in_bp?
	bool conflictwithBP(const BP & in_bp) {
		if (in_bp.rp < outer.lp || in_bp.lp > outer.rp) return false;
		if (in_bp.lp > inner.lp && in_bp.rp < inner.rp) return false;
		if (in_bp.lp < outer.lp && in_bp.rp > outer.rp) return false;
		return true;
	}
	bool conflictwith(stem st){
		if(st.outer.lp > outer.rp || st.outer.rp < outer.lp || 
		   (st.outer.lp > inner.lp && st.outer.rp < inner.rp) ||
		   (st.inner.lp < outer.lp && st.inner.rp > outer.rp))
		   return false;
		else return true;
	}
	bool containBP(const BP & in_bp) {
		if (in_bp.lp >= outer.lp && in_bp.rp <= outer.rp 
		 && in_bp.lp <= inner.lp && in_bp.rp >= inner.rp 
		 && in_bp.lp-outer.lp == outer.rp-in_bp.rp)
		 return true;
		return false;
	}
	bool containNT(const int & in_nt) {
		if (in_nt >= outer.lp && in_nt <= inner.lp) return true;
		if (in_nt >= inner.rp && in_nt <= outer.rp) return true;
		return false;
	}
	BP getBPcontainNT(const int & in_nt) {
		//return BP(in_nt, x) or BP(x, in_nt)
		if (in_nt>=outer.lp && in_nt<=inner.lp) return BP(in_nt, (outer.rp+outer.lp-in_nt));
		if (in_nt>=inner.rp && in_nt<=outer.rp) return BP((outer.rp+outer.lp-in_nt), in_nt);
	}
	void show() {
		//std::cout<<"["<<interval_s<<","<<interval_e<<"], 
		std::cout<<"("<<outer.lp<<","<<outer.rp<<") ~ ("<<inner.lp<<","<<inner.rp<<")";
	}
	int gethydrobonds(const std::string & seq) {
		int ret = 0; 
		int i = outer.lp; int j=outer.rp;
		while(i<=inner.lp) {
			if((seq[i] == 'A' && seq[j] == 'U') || (seq[i] == 'U' && seq[j] == 'A')) {
				ret += 2;
			}
			if((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'C' && seq[j] == 'G')) {
				ret += 3;
			}
			if((seq[i] == 'G' && seq[j] == 'U') || (seq[i] == 'U' && seq[j] == 'G')) {
				ret += 1;
			}
			i++; j--;
		}
		return ret;
	}
};

bool operator == (const stem & st1, const stem & st2);
bool operator != (const stem & st1, const stem & st2);
bool operator <  (const stem & st1, const stem & st2);
bool operator >  (const stem & st1, const stem & st2);

class stemlist {
	public:
	std::vector<stem> vstems;
	std::string s;
	
	stemlist(){}
	stemlist(const std::string & seq, const std::set<stem> & vstem_set);
	stemlist(const std::string & seq);
	void clear() {vstems.clear(); }
	void show();
//	void sort();
	int size() {return sz(vstems);}
	std::vector<stem> getstems_startfrom(int i);
	std::vector<stem> getstems_startfrom_within(int in_s, int in_e);
	std::vector<stem> getstems_endat(int in_e);
	std::vector<stem> getstems_startfrom_endat(int in_s, int in_e);
	std::vector<stem> getstems_within(int in_s, int in_e);

	std::vector<int> getstemindice_startfrom(int i);
	std::vector<int> getstemindice_startfrom_within(int in_s, int in_e);
	std::vector<int> getstemindice_endat(int in_e);
	std::vector<int> getstemindice_endat_within(int in_s, int in_e);
	std::vector<int> getstemindice_startfrom_endat(int in_s, int in_e);
	std::vector<int> getstemindice_enclosedby(int indx_s);
	std::vector<int> getstemindice_within(int in_s, int in_e);
	std::vector<int> getstemindice_inward_extend(int indx_s);
	std::vector<int> getstemindice_enclosedby_or_inward_extend(int indx_s);

	bool conflictwith(stem st);
	bool containST(stem st);
	
	void insert(stem & in_st);
	void insert(stemlist & in_stl);
	void PB(stem & in_st);
	void PB(stemlist & in_stl);
	stem & operator [] (const int & i);
	stem pop() {stem st = vstems[0]; vstems.erase(vstems.begin()); return st;}
	void add(stem & in_st);
	stem findLS(const stem & in_st);
	stem findRS(const stem & in_st);
};


class segment{
	public: 
	int lp;  //index of the left parenthesis
	int rp;	//index of the right parenthesis
	
	segment(int i1, int i2) {
		if(i2<i1) {std::cout<<"("<<i1<<","<<i2<<")"<<std::endl;}
		assert(i2>=i1);		
		lp=i1; rp=i2;
	}
	segment(BP bp) {assert(bp.lp<bp.rp); lp=bp.lp; rp=bp.rp;}
	bool overlap(stem st) {
		if(rp < st.outer.lp || lp > st.outer.rp ||
		   (st.inner.lp < lp && st.inner.rp > rp) )
		   return false;
		return true;
	}
	bool enclose(const stem & st) {
		if(lp <= st.outer.lp && rp >= st.outer.rp) 
			return true;
		return false;
	}
	bool valid() {
		if(lp>=0 && rp>=0 && lp<=rp) return true;
		else return false;
	}
	bool operator == (const segment & sg) {
		return (lp ==  sg.lp && rp == sg.rp);
	}
	bool operator != (const segment & sg) {
		return (lp != sg.lp || rp != sg.rp);
	}
	bool operator < (const segment & sg) {
		return (lp < sg.lp || (lp==sg.rp && rp>sg.rp));
	}
	bool operator > (const segment & sg) {
		return (lp > sg.lp || (lp==sg.lp && rp<sg.rp));
	}
	void show() {
		std::cout<<"("<<lp<<","<<rp<<") ";
	}
};

class seglist{
	public:
	std::vector<segment> vsegs;
	seglist(){}
	seglist(const std::vector<segment> & in_vsg) {
		vsegs = in_vsg;
	}
	int size() {return sz(vsegs);}
	void clear() {vsegs.clear();}
	bool overlap(stem st) {
		FOR(i,0,sz(vsegs)) {
			if(vsegs[i].overlap(st)) return true;
		}
		return false;
	}
	segment & operator [] (const int & i) {return vsegs[i];}
	void insert(segment & in_sg) {vsegs.insert(vsegs.begin(), in_sg);}
	void insert(seglist & in_sgl) {vsegs.insert(vsegs.begin(), 
	     in_sgl.vsegs.begin(), in_sgl.vsegs.end());}
	void PB(segment in_sg) {vsegs.PB(in_sg);}
	void PB(seglist in_sgl) {vsegs.insert(vsegs.end(), 
	     in_sgl.vsegs.begin(), in_sgl.vsegs.end());}
	segment pop() {segment sg = vsegs[0]; vsegs.erase(vsegs.begin()); return sg;}
	segment findRS(const segment & in_sg);
	void show() {
		FOR(i,0,sz(vsegs)){
			std::cout<<"segs ["<<i<<"/"<<sz(vsegs)-1<<"] :";
			vsegs[i].show();
			std::cout<<std::endl;
		}
	}
};

//##
class RStruct {//using base pairs
	public: 
	std::string      s;   //std::string representation of a RStruct
	std::vector<BP>  vbp; //std::vector of base pairs
	std::vector<int> unp; //indice of the unpaired parenthesis
	int         len; //length of the RStruct
	
	RStruct(){}
	RStruct(const RStruct & in_rs) {s=in_rs.s; vbp = in_rs.vbp; unp=in_rs.unp; len=in_rs.len;}
	RStruct(const std::string & str) {s=str; len=sz(s); vbp=s2v(str); unp=getunp(str);}
	RStruct(const int & l, const std::vector<BP> & in_vbp) {len=l; vbp=in_vbp; s=v2s(len, vbp); unp=getunp(len,vbp);}
	~RStruct(){len=0;s="";vbp.clear(); unp.clear();}
	
	bool operator == (const RStruct & rs);
	bool operator != (const RStruct & rs);
	
	void   show();
//	int    dist(const RStruct & rs);        //compute the distance between this and rs
	int    pairwithwhich(const int & query);//determine with which the query integer (bp) pairs in this structure; -1 not pair with any bp
	bool   del_bp_with_indx(const int & i); //delete base pair with index i from vbp and add lp and rp to unp; false failed; true succeeded
	bool   del_bp(const BP & in_bp);        //delete base pair in_bp from vbp and add lp and rp to unp; false failed; true succeed
	bool   del_bp_contain(const int & i);   //delete base pair with either left or right parenthesis equals to i; false failed; true succeed
	bool   add_bp(const int & i, const int & j); //add base pair (i,j) to vbp and delete them from unp; false failed; true succeed
	bool   add_bp(const BP & in_bp);        //add base pair in_bp to vbp and delete in_bp.lp and in_bp.rp from unp; false failed; true succeed
	bool   applyAction(const Action & in_a);        //apply action to this, but don't check clash
	bool   checkclash(const Action & in_a);        //check if in_ac leads to clash if applied
	std::vector<BP> intersectBP(const RStruct & rs); //get paired based that are in both this and rs
	std::vector<BP> getconflictingBP(const BP & in_bp); //get BPs in this RStruct that are conflicting with in_bp
	std::vector<stem> create_stems();	//identify stems in the secondary structure
};

class RStack{
	public:
	std::vector<stem> vstems; //all stacks in the structure
	int len; //the length of the sequence 
	
	RStack(){len=0;vstems.clear();}
	RStack(const std::vector<stem> in_vst, const int & in_len) {vstems = in_vst;len = in_len;}
	RStack(const std::string & str) {
		RStruct rs(str);
		vstems = rs.create_stems();
		len = sz(str);
		sort(vstems.begin(), vstems.end(), ascendingsortST);
	}

	bool ismatch(const std::string & seq);
	std::string tostring();
	int size() {return sz(vstems);}
	void show(){std::cout<<tostring()<<std::endl;}
	void show2(){
		std::cout<<"The stems in RStack is "<<std::endl;
		FOR(i,0,sz(vstems)) {
			vstems[i].show();
			std::cout<<std::endl;
		}
	}
	bool containST(const stem & st);//whether this instance contains stack st?
	int  inwhichST(const int & query); //in which stack is query.
	bool isLocOpt(const stemlist & in_stl);
	bool isLocOpt(const stemlist & in_stl, const std::vector<std::vector<int> > & H, const std::vector<std::vector<int> > & C);
	std::vector<RStack> getneighbor_delstack();
	std::vector<int> findConflictingST(stem & sti);//return indices of stems that conflict with sti
};



class PartStruct{
//partial structure: PartStruct = (seglist; stemlist)
	public: 
	seglist sgl;
	stemlist stl;
	
	PartStruct(){}
	PartStruct(const PartStruct & phi) {
		clear();
		sgl.vsegs.insert(sgl.vsegs.begin(), phi.sgl.vsegs.begin(), phi.sgl.vsegs.end());
		stl.vstems.insert(stl.vstems.begin(), phi.stl.vstems.begin(), phi.stl.vstems.end());
	}
	PartStruct(const seglist & in_sgl, const stemlist & in_stl){
		sgl = in_sgl; stl = in_stl;
	}
	PartStruct(const std::vector<segment> & in_vsg, const std::vector<stem> & in_vst) {
		sgl.vsegs = in_vsg; stl.vstems= in_vst;
	}
	void clear() {sgl.clear(); stl.clear();}
	void show() {
		std::cout<<"PartStruct segments:";
		FOR(i,0,sz(sgl))
			std::cout<<"["<<sgl.vsegs[i].lp<<", "<<sgl.vsegs[i].rp<<"],";
		std::cout<<std::endl<<"PartStruct stacks:";
		FOR(i,0,sz(stl)){
			stl.vstems[i].show();
		}
		std::cout<<std::endl;
	}
	bool isLocOpt(const stemlist & in_stl);
};

std::string traceback(const int & n, stemlist & in_stl); //trace back for constructing the best structure with minimum base pairs 


class SD{
	public:
	std::string str;
	double cd;
	SD(const std::string & in_s, const double & in_d) {
		str = in_s; cd = in_d;
	}
	bool operator < (const SD & in_sd) {
		return (cd < in_sd.cd);
	}
	bool operator > (const SD & in_sd) {
		return (cd > in_sd.cd);
	}
	bool operator == (const SD & in_sd) {
		return (str == in_sd.str && cd == in_sd.cd);
	}
	bool operator != (const SD & in_sd) {
		return (str != in_sd.str || cd != in_sd.cd);
	}
};

class RNA{
	public:
	std::string str;
	int e;
 	std::string hishape;
	RNA(const std::string & in_str, const int & in_e) {
		str = in_str; e = in_e;
	}
 	RNA(const std::string & in_str, const int & in_e, const std::string & in_hishape) {
 		str = in_str; e = in_e; hishape = in_hishape;
 	}
};

bool descendingsortRNA(RNA rna1, RNA rna2);
bool ascendingsortRNA(RNA rna1, RNA rna2);
#endif

