#include "RNAStackFold.hh"

extern TurnerE TE;
extern Config cf;

int TerminalPenalty(const char & x, const char & y)
{
	return TerminalPenalty(TE, x, y);
}

int TerminalPenalty(TurnerE & TE, const char & x, const char & y) {
	if ((x=='U') || (y=='U'))
		return TE.MtermianAU;
	return 0;
}

int d(const std::string & seq, const int & i1, const int & i2, const int & i3){
	return d(TE, seq, i1, i2, i3);
}

int d(TurnerE & TE, const std::string & seq, const int & i1, const int & i2, const int & i3) {
/*
# 5'-->3'
# d5(h,k-1,k)   d3(i,i+1,j)
# .....\k.....i/-----\
#       |||||||      |
#       h.....j      |
# ...../       \_____/
# d3(h,h+1,k)   d5(i,j-1,j)
# 3'<--5' 
*/
//	cout<<"i1="<<i1<<", i2="<<i2<<", i3="<<i3<<std::endl;
	assert(i1>=0 && i1<i2 && i2<=sz(seq)-1);
	assert(match(seq[i1], seq[i2]));
	if(i3==-1 || i3==sz(seq)) return 0;
	if(i3==i1-1)	{//k=i1,h=i2;
		int k=i1, h=i2;
		return d5(TE, seq, h,k-1,k);
	} else if(i3==i1+1) {
		int i=i1, j=i2;
		return d3(TE, seq, i,i+1,j);
	} else if(i3==i2+1) {
		int k=i1, h=i2;
		return d3(TE, seq, h,h+1,k);
	} else if(i3==i2-1){
		int i=i1, j=i2;
		return d5(TE, seq, i,j-1,j);
	} else {
		assert(false);
	}
}

int d5(TurnerE & TE, const std::string & seq, const int & i1, const int & i2, const int & i3) {
	std::string name = "XYZ";
	name[0] = seq[i1]; name[1] = seq[i2]; name[2] = seq[i3];
//	cout<<"d5("<<name<<")="<<TE.d5E[name]<<std::endl;
	int ret = TE.d5E[name];
	return ret;
}


//seq[i1] pairs with seq[i2]
//#5'-->3'      5'-->3'
//# 1 2             3       ---\ d3(123)
//# 3             2 1       ---/
//#3'<--5'      3'<--5'
int d3(TurnerE & TE, const std::string & seq, const int & i1, const int & i2, const int & i3) {
	std::string name = "XYZ";
	name[0] = seq[i1]; name[1] = seq[i2]; name[2] = seq[i3];
//	cout<<"d3("<<name<<")="<<TE.d3E[name]<<std::endl;
	int ret = TE.d3E[name];
	return ret;
}

int computeH(TurnerE & TE, const std::string & seq, stem & st) {
	int n = sz(seq);
	int x = st.outer.lp; int y = st.outer.rp;
	int p = st.inner.lp; int q = st.inner.rp;
	int Si = computeS(TE, seq, x, y, p, q);
	
	int hairpinlen = q-p-1;
	int e_loop = 0;
//internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30)           \
--> 1.079, Ref:http://www.bioinfo.rpi.edu/zukerm/cgi-bin/efiles-3.0.cgi#ASINT1x2
	if(hairpinlen>MAXLOOPLEN) {
		e_loop = TE.HAIRPIN[MAXLOOPLEN] + TE.Mlongloop*log(double(hairpinlen)/30);
	} else {
		e_loop = TE.HAIRPIN[hairpinlen];
	}
	std::string tmp;
	FOR(k,0,n){
		if(k>=x && k<=p) tmp += '(';
		else if(k>=q && k<=y) tmp += ')';
		else tmp += '.';
	}
#ifdef MAINDEBUG
	cout<<"structure: "<<std::endl<<tmp<<std::endl;
	cout<<"\thairpin length = "<<hairpinlen<<", loop destablizing energy = "<<e_loop<<std::endl;
#endif
	int energy = 0; 
	int e_terminalsmhe = 0; 
	int e_AUpenalty = 0; //DONT count left AU penalty, which should be          \
	calculated when computing H[][] or C[][] or FM[][] or FM1[][]
	int e_tetraLoopE = 0;
	assert(hairpinlen >= cf.MINLOOPLEN);
	
	if(hairpinlen > cf.MINLOOPLEN){
		//if length of the hairpin loop is greater or equal to 4, only count     \
		the terminal mismatch stacking energy, dont count terminal AU GU         \
		penalty (right end)
		std::string name = "";
		name += seq[p]; name += seq[p+1];
		name += seq[q-1]; name += seq[q];
	 	e_terminalsmhe = TE.TerminalSMHE[name];
#ifdef MAINDEBUG
	 	cout<<"\tterminal single mismatch stack = "<<name<<", energy = "<<e_terminalsmhe<<std::endl;
#endif
	 	if(hairpinlen == 4) {
	 		//check tetra loop energy;
	 		name = seq.substr(p, 6);
	 		std::map<std::string, int>::iterator it;
	 		it = TE.TetraLoopE.find(name);
	 		if(it != TE.TetraLoopE.end()) {
	 			e_tetraLoopE = (*it).second;
#ifdef MAINDEBUG
	 			cout<<"\tspecial tetra loop = "<<name<<", energy:"<<e_tetraLoopE<<std::endl;
#endif
			}
			energy = Si+e_terminalsmhe+e_loop+e_tetraLoopE;
#ifdef MAINDEBUG
			cout<<"Energy [hairpinlen="<<hairpinlen<<"] = StackingE+TerminalSMHE+LoopE+TetraLoopE = "<<Si<<"+"<<e_terminalsmhe<<"+"<<e_loop<<"+"<<e_tetraLoopE<<" = "<<energy<<std::endl;
#endif
		} else {
			energy = Si+e_terminalsmhe+e_loop;
#ifdef MAINDEBUG
			cout<<"Energy [hairpinlen="<<hairpinlen<<"] = StackingE+TerminalSMHE+LoopE = "<<"+"<<Si<<"+"<<e_terminalsmhe<<"+"<<e_loop<<" = "<<energy<<std::endl;
#endif
		}
	} else if (hairpinlen==cf.MINLOOPLEN){
		//if length of the hairpin loop is 3, only count the (right) terminal    \
		AU/GU penalty, do not count terminal mismatch stacking energy.
		e_AUpenalty = TerminalPenalty(TE, seq[p], seq[q]);
		energy = e_AUpenalty+Si+e_loop;
#ifdef MAINDEBUG
		cout<<"Energy [hairpinlen="<<hairpinlen<<"] = AUpenalty(inner)+StackingE+LoopE = "<<
		  e_AUpenalty<<"+"<<Si<<"+"<<e_loop<<" = "<<energy<<std::endl;
#endif
	} else {
		energy = OOE;
	}
	
	//check for GU closeure preceded by GG
	if(seq[p] == 'G' && seq[q] == 'U') {
		if(p >= 2) {
			if(seq[p-1] == 'G' && seq[p-2] == 'G')
				energy += TE.Mgubous;
		}
	}
	
	//check for a poly-c loop
	bool polyc_loop = true;
	FOR(k,p+1,q) {
		if(seq[k] != 'C') { 
			polyc_loop = false;
			break;
		}
	}
	if(polyc_loop) {
		if(hairpinlen==3) {
			energy += TE.Mc3;
		} else {
			energy += TE.Mcintercept + TE.Mcslope * hairpinlen;
		}
	}
	return energy;
}

int computeBulge(TurnerE & TE, const std::string & seq, const int & x, const int & y, 
   const int & p, const int & q) {
//x pairs with y and p pairs with q
//...xp....
//||||\\\\\\
//...y q.....
//   |  \
//   |__/

//   ______
//   |    /
//   |   /
//...x  p.....
//|||| //////
//...yq.....
	int ret = 0;
	assert(x<p && p<q && q<y);
	int leftgap = p-x-1;
	int rightgap = y-q-1;
	int longgap = (leftgap>rightgap)?(leftgap):(rightgap);
	int shortgap = (leftgap<rightgap)?(leftgap):(rightgap);
	assert(shortgap==0);	
	if(longgap==1) {
		//according to mfold 3.5 algorithm.cpp line 3349
		std::string name = "";
		name += seq[x]; name += seq[p]; name += seq[q]; name += seq[y];
		ret = TE.BULGE[longgap] + TE.STACKING[name];
	} else if (longgap>30) {
		ret = TE.BULGE[MAXLOOPLEN]+ TE.Mlongloop*log(longgap/MAXLOOPLEN) + 
		      TerminalPenalty(TE, seq[x], seq[y]) + 
		      TerminalPenalty(TE, seq[p], seq[q]);
	} else {
		ret = TE.BULGE[longgap] + 
				TerminalPenalty(TE, seq[x], seq[y]) +
				TerminalPenalty(TE, seq[p], seq[q]);
	}
	return ret;
}

int computeInterior(TurnerE & TE, const std::string & seq, const int & x, const int & y, 
   const int & p, const int & q) {
   int ret = 0;
	assert(x<p && p<q && q<y);
	int leftgap = p-x-1;
	int rightgap = y-q-1;
	int longgap = (leftgap>rightgap)?(leftgap):(rightgap);
	int shortgap = (leftgap<rightgap)?(leftgap):(rightgap);
	assert(shortgap!=0 && longgap!=0);	
	int allgap = shortgap+longgap;
	if(leftgap==1 && rightgap==1) {
		std::string name = "";
		name += seq[x]; name += seq[x+1]; name += seq[x+2];
		name += seq[y-2]; name += seq[y-1]; name += seq[y];
		ret = TE.SingleMIE[name];
	} else if(leftgap==1 && rightgap==2) {
		std::string name = "";
		name += seq[x]; name += seq[x+1]; name += seq[x+2];
		name += seq[y-3]; name += seq[y-2]; name += seq[y-1]; name += seq[y];
		ret = TE.OneTwoMIE[name];
	} else if(leftgap==2 && rightgap==1) {
		std::string name = "";
		name += seq[y-2]; name += seq[y-1]; name += seq[y];
		name += seq[x]; name += seq[x+1]; name += seq[x+2]; name += seq[x+3];
		ret = TE.OneTwoMIE[name];
	} else if (leftgap==2 && rightgap==2) {
		std::string name = "";
		name = seq.substr(x, 4);
		name += seq[y-3]; name += seq[y-2]; name += seq[y-1]; name += seq[y];
//		cout<<"2x2 interior "<<name<<", energy="<<TE.TwoTwoMIE[name]<<std::endl;
		ret = TE.TwoTwoMIE[name];
	}

	if(allgap > MAXLOOPLEN) {
		ret = TE.INTERNAL[MAXLOOPLEN]+ TE.Mlongloop*log(allgap/MAXLOOPLEN);
		ret += (TE.Mninio_max < (TE.Mninio_m*(longgap-shortgap)))? (TE.Mninio_max): (TE.Mninio_m*(longgap-shortgap));
	} else if(allgap > 4 || (allgap==4 && longgap==3 && shortgap==1)) {
		ret = TE.INTERNAL[allgap];

//		cout<<"shortgap="<<shortgap<<", longgap="<<longgap<<std::endl;
//		cout<<"Internal["<<allgap<<"]="<<TE.INTERNAL[allgap]<<std::endl;
		ret += min(TE.Mninio_max, (TE.Mninio_m*(longgap-shortgap)));
//		cout<<"Mninio="<<min(TE.Mninio_max, (TE.Mninio_m*(longgap-shortgap)))<<std::endl;
	}
	if(shortgap==1 && longgap > 2) {
	//according to mfold, assume we had a loop composed of AA between
		std::string name = "";
		name += seq[x]; name += 'A'; name += 'A'; name += seq[y];
		ret  += TE.TerminalSMIE[name];
		name = ""; 
		name += seq[q]; name += 'A'; name += 'A'; name += seq[p];
		ret  += TE.TerminalSMIE[name];
	}
	if(shortgap>=2 && longgap > 2) {
		std::string name = "";
		name += seq[x]; name += seq[x+1]; name += seq[y-1]; name += seq[y];
		ret  += TE.TerminalSMIE[name];
//		cout<<"SMIE["<<name<<"]="<<TE.TerminalSMIE[name]<<std::endl;
		
		name = "";
		name += seq[q]; name += seq[q+1]; name += seq[p-1]; name += seq[p];
		ret  += TE.TerminalSMIE[name];
//		cout<<"SMIE["<<name<<"]="<<TE.TerminalSMIE[name]<<std::endl;
	}
//	cout<<"ret="<<ret;
	return ret;
}

int computeS(TurnerE & TE, const std::string & seq, stem & st) {
	return computeS(TE, seq, st.outer.lp, st.outer.rp, st.inner.lp, st.inner.rp);
}

int computeS(TurnerE & TE, const std::string & seq, const int & x, const int & y, const int & p, const int & q) {
	int e_stacking = 0;
	FOR(t,0,p-x) {
		int i1 = x+t; int i2 = x+t+1;
		int i3 = y-(t+1); int i4 = y-t;
		//i1,i2,i3,i4 (5'--->3')
		std::string name = "";
		name += seq[i1]; name += seq[i2];
		name += seq[i3]; name += seq[i4];
#ifdef MAINDEBUG
//		cout<<"\tstack = "<<name<<", stacking energy = "<<TE.STACKING[name]<<std::endl;
#endif
		e_stacking += TE.STACKING[name];
	}
#ifdef MAINDEBUG
//	cout<<"\tStacking energy in all = "<<e_stacking<<std::endl;
#endif
	return e_stacking;
}

int computeS(TurnerE & TE, const std::string & seq, stem & s, stem & t) {
//(((((((((((((.........)))))))))))))))))
//...........((((.....)))))..............
	int sx = s.outer.lp; int sy = s.outer.rp;
	int sp = s.inner.lp; int sq = s.inner.rp;
	int tx = t.outer.lp; int ty = t.outer.rp;
	int tp = t.inner.lp; int tq = t.inner.rp;
	
	assert(s.enclose(t) || t.inward_extend(s));
	if(s.enclose(t)) {
		return computeS(TE, seq, sx, sy, sp, sq);
	}
	if(t.inward_extend(s)) {
		assert(ty-sq!=sp-tx);
		if( (sp<tx && tq<=sq && sq<=ty && ty<sy) ||
			 (ty>=sq && sp>=tx && ty-sq > sp-tx)) {
			int star = sx+sy-ty-1;
			return computeS(TE, seq, sx, sy, star, ty+1);
		}
		if( (ty<sq && sx<tx && tx<=sp && sp<=tp) ||
			 (ty>=sq && sp>=tx && ty-sq<sp-tx) ){
			int star = sx+sy-tx+1;
			return computeS(TE, seq, sx, sy, tx-1, star);
		}
	}
	assert(false);
}

int computephi(TurnerE & TE, const std::string & seq, stem & s, stem & t) {
	int ret = 0;
	assert(s.enclose(t) || t.inward_extend(s));
	int sx = s.outer.lp; int sy = s.outer.rp;
	int sp = s.inner.lp; int sq = s.inner.rp; 
	int tx = t.outer.lp; int ty = t.outer.rp;
	int tp = t.inner.lp; int tq = t.inner.rp;
#ifdef MIANDEBUG
	cout<<"s="<<sx<<","<<sp<<","<<sq<<","<<sy<<std::endl;
	cout<<"t="<<tx<<","<<tp<<","<<tq<<","<<ty<<std::endl;
#endif	
	int leftgap  = tx-sp-1; //t.x-s.p-1
	int rightgap = sq-ty-1; //s.q-t.y-1
	int shortgap = (leftgap>rightgap)?(rightgap):(leftgap);
	int longgap  = (leftgap>rightgap)?(leftgap):(rightgap);
#ifdef MIANDEBUG
	cout<<"longgap="<<longgap<<", shortgap="<<shortgap<<std::endl;
#endif
	if(shortgap==0) { //bulge
		ret = computeBulge(TE, seq, sp, sq, tx, ty);
	} else if(shortgap>0) {
		ret = computeInterior(TE, seq, sp, sq, tx, ty);
	} else if(shortgap<0) {
		assert(ty-sq!=sp-tx);
		if( (sp<tx && tq<=sq && sq<=ty && ty<sy) ||
			 (ty>=sq && sp>=tx && ty-sq > sp-tx)) {
			int star = sx+sy-ty-1;
			ret = computeBulge(TE, seq, star, ty+1, tx, ty);
		}
		if( (ty<sq && sx<tx && tx<=sp && sp<=tp) ||
			 (ty>=sq && sp>=tx && ty-sq<sp-tx) ){
			int star = sx+sy-tx+1;
			ret = computeBulge(TE, seq, tx-1, star, tx, ty);
		}
	}
	return ret;
}

stem getUsedPartOfStem(stem & s, stem & t) {
#ifdef MAINDEBUG	
	cout<<"in getUsedPartOfStem"<<std::endl;
	cout<<"s is "; s.show();
	cout<<"t is "; t.show();
#endif	
	int sx = s.outer.lp; int sy = s.outer.rp;
	int sp = s.inner.lp; int sq = s.inner.rp;
	int tx = t.outer.lp; int ty = t.outer.rp;
	int tp = t.inner.lp; int tq = t.inner.rp;
	assert(s.enclose(t) || t.inward_extend(s));
	if(s.enclose(t)) {
#ifdef MAINDEBUG	
		cout<<"s enclose t, return s"<<std::endl;
#endif
		return s;
	}
	if(t.inward_extend(s)) {
		assert(ty-sq!=sp-tx);
		if( (sp<tx && tq<=sq && sq<=ty && ty<sy) ||
			 (ty>=sq && sp>=tx && ty-sq > sp-tx)) {
			int star = sx+sy-ty-1;
			// computeBulge(seq, star, ty+1, tx, ty);
			// computeS(seq, sx, sy, star, ty+1);
			stem ret(BP(sx, sy), BP(star, ty+1));
#ifdef MAINDEBUG
			cout<<"case 1"<<std::endl;
			cout<<"star="<<star<<std::endl;
			cout<<"return "; ret.show();
#endif
			assert(ret.outer.lp <= ret.inner.lp && ret.inner.lp < ret.inner.rp && ret.inner.rp <= ret.outer.rp);
			return ret;
		}
		if( (ty<sq && sx<tx && tx<=sp && sp<=tp) ||
			 (ty>=sq && sp>=tx && ty-sq<sp-tx) ) {
			int star = sx+sy-tx+1;
			//computeBulge(seq, tx-1, star, tx, ty);
			//computeS(seq, sx, sy, tx-1, star);
			stem ret(BP(sx, sy), BP(tx-1, star));
#ifdef MAINDEBUG
			cout<<"star="<<star<<std::endl;
			cout<<"case 2"<<std::endl;
			cout<<"return "; ret.show();
#endif
			assert(ret.outer.lp <= ret.inner.lp && ret.inner.lp < ret.inner.rp && ret.inner.rp <= ret.outer.rp);
			return ret;
		}
	}
	assert(false);
}

stem getUsedPartOfStem(stemlist & in_stl, const int & indx_s, const int & indx_t) {
	stem s = in_stl.vstems[indx_s];
	stem t = in_stl.vstems[indx_t];
	return getUsedPartOfStem(s, t);
}

void showd(TurnerE & TE, const std::string & seq, const int & a, const int & b, const int & c, std::ostream & fout) {
	fout<<"d("<<a<<","<<b<<","<<c<<")="<<d(TE,seq,a,b,c)<<std::endl;
}

void showmind(TurnerE & TE, const std::string & seq, const int & a, const int & b, const int & c, const int & e, const int & f, const int & g, std::ostream & fout) {
	fout<<"min(d("<<a<<","<<b<<","<<c<<"), d("<<e<<","<<f<<","<<g<<"))="<<min(d(TE,seq,a,b,c), d(TE,seq,e,f,g))<<std::endl;
}

void showH(TurnerE & TE, const std::string & seq, stem & st, std::ostream & fout) {
	int sx = st.outer.lp; int sy = st.outer.rp;
	int sp = st.inner.lp; int sq = st.inner.rp;
	
	fout<<"H[("<<sx<<","<<sy<<")~("<<sp<<","<<sq<<")]="<<computeH(TE,seq,st)<<std::endl;
}

void showS(TurnerE & TE, const std::string & seq, stem & st, std::ostream & fout) {
	int sx = st.outer.lp; int sy = st.outer.rp;
	int sp = st.inner.lp; int sq = st.inner.rp;

	fout<<"S[("<<sx<<","<<sy<<")~("<<sp<<","<<sq<<")]="<<computeS(TE,seq,st)<<std::endl;
}

void showphi(TurnerE & TE, const std::string & seq, stem & st, stem & stk, std::ostream & fout) {
	int sx = st.outer.lp; int sy = st.outer.rp;
	int sp = st.inner.lp; int sq = st.inner.rp;
	int tx = stk.outer.lp; int ty = stk.outer.rp;
	int tp = stk.inner.lp; int tq = stk.inner.rp;
	
	fout<<"phi[("<<sx<<","<<sy<<")~("<<sp<<","<<sq<<")][("<<tx<<","<<ty<<")~("<<tp<<","<<tq<<")]="<<computephi(TE,seq,st,stk)<<std::endl;
}

void showMphi(TurnerE & TE, std::ostream & fout) {
	fout<<"Mphi="<<TE.Mphi<<std::endl;
}

void showMc(TurnerE & TE, std::ostream & fout) {
	fout<<"Mc="<<TE.Mc<<std::endl;
}

void showTP(TurnerE & TE, const std::string & seq, const int & x, const int & y, std::ostream & fout) {
	int res = TerminalPenalty(TE, seq[x], seq[y]);
	if(res != 0) {
		fout<<"AU("<<x<<","<<y<<")="<<res<<std::endl;
	}
}

int computeEnergy(TurnerE & TE, const std::string & seq, RStack & in_stl) {
	if(sz(in_stl)==0) return 0; 
	int ret = 0;
	std::vector<int> vinside;
	//vinside[i] contains a 
	//0 - Hairpin ((((....))))
	//1 - Interior loop or bulge ((((..((((....))))..))))
	//n, n>=2 - Multiple loops ((((....((((....))))....((((....)))))....))))
	//where n is the number of branches

	std::vector<int> voutside;
	//0 - not enclosed by any stem
	//1 - enclosed by an interior loop
	//2 - enclosed by a multi loop

	std::vector<std::vector<int> > relations;
	
	FOR(i,0,sz(in_stl.vstems))	{
		vinside.PB(0);
		voutside.PB(0);
	}
//	cout<<"show stl"<<std::endl;
//	FOR(i,0,sz(in_stl.vstems)) {
//		in_stl.vstems[i].show(); cout<<std::endl;
//	}
//	cout<<"end show"<<std::endl;
	//relations[i][j] has the relationships between i and j
	//0  - no relationship
	//1  - sti enclose stj and there does not exist any other stem stk such that stk enclose stj and sti enclose stk
	//-1  - sti is enclosed by stj, and there does not exist any other stem stk such that stk sti is enclosed by stk and stk is enclosed by stj 
   //2  - sti enclose stj and there exist a stem stk such that ...
   //-2 - sti is enclosed by stj, and there exist a stem stk such that ...	
	FOR(i,0,sz(in_stl.vstems)) {
		std::vector<int> tmpvi;
		relations.PB(tmpvi);
		FOR(j,0,sz(in_stl.vstems)){
			relations[i].PB(0);
		}
	}
	
	FOR(i,0,sz(in_stl.vstems)) {
		stem sti = in_stl.vstems[i];
		FOR(j,0,sz(in_stl.vstems)) {
			stem stj = in_stl.vstems[j];
			if(i==j) {
				relations[i][j] = 0; 
				continue;
			}
			if(sti.enclose(stj)) {
				bool exist = false;
				FOR(k,0,sz(in_stl.vstems)){
					if(k==i || k==j) continue;
					stem stk = in_stl.vstems[k];
					if(sti.enclose(stk) && stk.enclose(stj)){
						exist = true;
						break;
					}
				}
				if(exist) relations[i][j] = 2;
				else relations[i][j] = 1;
			} else if(stj.enclose(sti)){
				bool exist = false;
				FOR(k,0,sz(in_stl.vstems)){
					if(k==i || k==j) continue;
					stem stk = in_stl.vstems[k];
					if(stj.enclose(stk) && stk.enclose(sti)){
						exist = true;
						break;
					}
				}
				if(exist) relations[i][j] = -2;
				else relations[i][j] = -1;
			} else {
				relations[i][j] = 0;
			}
		}
	}

	FOR(i,0,sz(in_stl.vstems)) {
		stem sti = in_stl.vstems[i];
		FOR(j,0,sz(in_stl.vstems)) {
			if(relations[i][j] == 1) {
				vinside[i]++;
			}
		}
	}
	
	FOR(i,0,sz(in_stl.vstems)) {
		stem sti = in_stl.vstems[i];
		FOR(j,0,sz(in_stl.vstems)) {
			if(relations[i][j] == -1) {
				if(vinside[j]==0) assert(false);
				else if(vinside[j]==1) {
					voutside[i] = 1;
				} else if(vinside[j]>=2) {
					voutside[i] = 2;
				}
			}
		}
	}
	
	std::vector<int> vstindx;
	vstindx.PB(0);
	int prei = -1;
	while(sz(vstindx)>0){
		int i = vstindx[0];
		assert(i>prei);
		prei = i;
		stem st = in_stl.vstems[i];
//		FOR(k,0,sz(vstindx)) {
//			in_stl.vstems[vstindx[k]].show(); cout<<std::endl;
//		}
//		cout<<"input vstems"<<std::endl;
//		cout<<"i="<<i<<", stem = "<<std::endl;
//		st.show(); cout<<std::endl;

		vstindx.erase(vstindx.begin());
		
		if(vinside[i] == 0) {
			ret += computeH(TE, seq, st);
//			showH(TE,seq,st);
		} else if (vinside[i] == 1){
			for(int k=sz(in_stl.vstems)-1; k>=0; k--) {
				if(relations[i][k] == 1) {
					stem stk = in_stl.vstems[k];
					ret += computeS(TE, seq, st)+computephi(TE, seq, st, stk);
//					showS(TE,seq,st);
//					showphi(TE,seq,st,stk);
					vstindx.insert(vstindx.begin(),k);
					break;
				}
			}
		} else if (vinside[i] >= 2){
			ret += computeS(TE, seq, st) + TE.Mphi + TE.Mc;
//			showS(TE,seq,st);
//			showMphi(TE);
//			showMc(TE);
			
			for(int k=sz(in_stl.vstems)-1; k>=0; k--) {
				if(relations[i][k] == 1) {
					ret += TE.Mphi;
//					showMphi(TE);
					vstindx.insert(vstindx.begin(), k);
				}
			}
		} else assert(false);
		
		if(sz(vstindx)==0 && i<sz(in_stl.vstems)-1) {
			vstindx.PB(i+1);
		}
	}
	FOR(i,0,sz(in_stl.vstems)) {
		stem st = in_stl.vstems[i];
		int x = st.outer.lp; int y = st.outer.rp;
		int p = st.inner.lp; int q = st.inner.rp;
		
		if(voutside[i]==1) {
			//inside an interior loop, no terminal AU penalty for x,y
		} else {
			ret += TerminalPenalty(TE, seq[x], seq[y]);
//			showTP(TE,seq,x,y);
		}
		
		if(vinside[i]==1 || vinside[i]==0) {
			//enclose an interior loop, no terminal AU penalty for p,q
			//a hairpin loop, no temrinalAU penalty for p,q
		} else {
			ret += TerminalPenalty(TE, seq[p], seq[q]);
//			showTP(TE,seq,p,q);
		}
		
		int k = x-1; int h = x-2;
		int a = y+1; int b = y+2;
		int c = p+1; int m = p+2;
		int e = q-1; int f = q-2;
		int k1, h1, a1, b1, c1, m1, e1, f1;
		k1 = h1 = a1 = b1 = c1 = m1 = e1 = f1 = -1;
		FOR(j,0,sz(in_stl.vstems)) {
			if(i==j) continue;
			stem stj = in_stl.vstems[j];
			if(stj.pairwithwhich(k)!=-1) k1 = stj.pairwithwhich(k);
			if(stj.pairwithwhich(h)!=-1) h1 = stj.pairwithwhich(h);
			if(stj.pairwithwhich(a)!=-1) a1 = stj.pairwithwhich(a);
			if(stj.pairwithwhich(b)!=-1) b1 = stj.pairwithwhich(b);
			if(stj.pairwithwhich(c)!=-1) c1 = stj.pairwithwhich(c);
			if(stj.pairwithwhich(m)!=-1) m1 = stj.pairwithwhich(m);
			if(stj.pairwithwhich(e)!=-1) e1 = stj.pairwithwhich(e);
			if(stj.pairwithwhich(f)!=-1) f1 = stj.pairwithwhich(f);
		}
		
		if(k1!=-1) {
			//seq[k] is pairing with seq[k1], don't add d(x,y,k)
		} else if(k1==-1 && h1!=-1) {
			if(voutside[i]==1) {
			//interior loop, no dangle energy contributions
			} else {
				if(h1>h) {
					//min(d(x,y,k), d(h1,h,k)) should have been included when processing a stem with inner=(h,h1)
				} else {
					//min(d(x,y,k), d(h1,h,k)) should have been included when processing a stem enclosed by (h1,h)
				}
			}
		} else if(k1==-1 && h1==-1) {
			if(voutside[i]==1) {
			//interior loop, no dangle energy contributions
			} else {
				ret += d(TE, seq, x, y, k);
//				showd(TE,seq,x,y,k);
			}
		}
		if(a1!=-1) {
			//seq[a] is pairing with seq[a1], don't add d(x,y,a)
		} else if(a1==-1 && b1!=-1) {
			if(voutside[i]==1) {
				//interior loop, no dangle energy contributions
			} else {
				if(b1>b) {
					ret += min(d(TE, seq, x, y, a), d(TE, seq, b, b1, a));
//					showmind(TE,seq,x,y,a,b,b1,a);
				} else {
					ret += min(d(TE, seq, x, y, a), d(TE, seq, b1, b, a));
//					showmind(TE,seq,x,y,a,b1,b,a);
				}
			}
		} else if(a1==-1 && b1==-1) {
			if(voutside[i]==1) {
				//interior loop, no dangle energy contributions
			} else {
				ret += d(TE, seq, x, y, a);
//				showd(TE,seq,x,y,a);
			}
		}
		
		if(c1!=-1) {
			//dont add d(p,q,c)
		} else if(c1==-1 && m1!=-1) {
			if(vinside[i]==0) {//Hairpin
				assert(false);
			} else if(vinside[i]==1) {//interior or bulge
			} else if(vinside[i]>=2) {//multiloop
				ret += min(d(TE, seq, p, q, c), d(TE, seq, m, m1, c));
//				showmind(TE,seq,p,q,c,m,m1,c);
			}
		} else if(c1==-1 && m1==-1) {
			if(vinside[i]==0) {//Hairpin
			} else if(vinside[i]==1) {//interior or bulge
			} else if(vinside[i]>=2) {//multiloop
				ret += d(TE, seq, p, q, c);
//				showd(TE,seq,p,q,c);
			}
		}

		if(e1!=-1){
		} else if(e1==-1 && f1!=-1) {
			if(vinside[i]==0) {//Hairpin
				assert(false);
			} else if(vinside[i]==1) {//interior or bulge
			} else if(vinside[i]>=2) {//multiloop
				//min(d(p,q,e), d(f1, f, e)) should have been included when processing a stem enclosed by (f1, f);
			}
		} else if(e1==-1 && f1==-1) {
			if(vinside[i]==0) {//Hairpin
			} else if(vinside[i]==1) {//interior or bulge
			} else if(vinside[i]>=2) {//multiloop
				ret += d(TE, seq, p, q, e);
//				showd(TE,seq,p,q,e);
			}
		}
	}
	return ret;
}



