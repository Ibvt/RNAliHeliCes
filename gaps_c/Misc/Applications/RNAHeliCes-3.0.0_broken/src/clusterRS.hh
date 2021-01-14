#ifndef _CLUSTERRS_HH_
#define _CLUSTERRS_HH_
#include "RNAStackFold.hh"
#include <utility>

const double samedirect = 0.5;

//##
class stemstem {
public: 
	stem sti;
	stem stj;
	stemstem(){}
	stemstem(stem & in_sti, stem & in_stj) {
		sti = in_sti; stj = in_stj;
	}
};


bool operator == (const stemstem & stst1, const stemstem & stst2);
bool operator != (const stemstem & stst1, const stemstem & stst2);
bool operator <  (const stemstem & stst1, const stemstem & stst2);
bool operator >  (const stemstem & stst1, const stemstem & stst2);

std::vector<int> str2vectorint(const std::string & str);
std::vector<int> RStack2vectorint(const RStack & rs) ;

/* vi[.]= -1: value is not set;
   vi[.]= 0 : set to 0 for sure;
   vi[.]= 1 : set to 1 for sure;
   vi[.]= 2 : set to samedirect or greater; */
/*
std::vector<int> compute_conflictingdegree_between_stems(stem & sti, stem & stj);
std::vector<int> get_conflictingdegree_between_stems(std::map<stemstem, std::vector<int> > & CDSTmap, stem & sti, stem & stj);
std::vector<double> compute_conflictingdegree(const std::string & seq, RStack & rs1, RStack & rs2, std::map<stemstem, std::vector<int> > & CDSTmap);
std::vector<double> compute_conflictingdegree(const std::string & seq, const std::string & s1, const std::string & s2, std::map<stemstem, std::vector<int> > & CDSTmap);
double compute_pairwise_conflictingdegree(TurnerE & TE, const std::string & seq, 
       RStack & rs1, RStack & rs2, std::map<stem, int> & stemEmap, std::map<stemstem, std::vector<int> > & CDSTmap);

*/
std::vector<int> compute_conflictingdegree_between_stems(stem & sti, stem & stj);
std::vector<int> get_conflictingdegree_between_stems(std::map<stemstem, std::vector<int> > & CDSTmap, stem & sti, stem & stj);
std::vector<double> compute_conflictingdegree(const std::string & seq, RStack & rs1, RStack & rs2, std::map<stemstem, std::vector<int> > & CDSTmap);
std::vector<double> compute_conflictingdegree(const std::string & seq, const std::string & s1, const std::string & s2, std::map<stemstem, std::vector<int> > & CDSTmap);
double compute_pairwise_conflictingdegree(TurnerE & TE, const std::string & seq, 
       RStack & rs1, RStack & rs2, std::map<stem, int> & stemEmap, std::map<stemstem, std::vector<int> > & CDSTmap);


double compute_distance(const std::string & seq, const std::string & s1, const std::string & s2);
//computing the pairwise distance (approximated barrier) between rs1 and rs2.
std::map<stem, int> computeStemE(TurnerE & TE, const std::string & seq, std::vector<stem> & vst);//compute stacking energy of all stems in vst and return a std::map std::mapping stems to their energies (which are integers);
int getStemE(TurnerE & TE, const std::string & seq, stem & st, std::map<stem, int> & stemEmap);//obtain the stacking energy of stem st from std::map stemEmap, if it is not available, call computeS and add it to stemEmap.


#endif
