#include "RNAStackFold.hh"
extern Config cf;

/******************************************************************************
                              PartStructT
 ******************************************************************************/
PartStructT::PartStructT(const seglist & in_sgl, const std::vector<ETYPE> in_sgt, 
   const stemlist & in_stl, const int & in_e, const std::string & in_s)
{
	sgl = in_sgl;
	sgt = in_sgt;
	stl = in_stl;
	elphi = in_e;
	log = in_s;
}

PartStructT::PartStructT(const PartStructT & phi) 
{
	clear();
	seglist in_sgl = phi.sgl;
	stemlist in_stl= phi.stl;
	sgl.insert(in_sgl);
	sgt.insert(sgt.begin(), phi.sgt.begin(), phi.sgt.end());
	stl.insert(in_stl);
	elphi = phi.elphi;
	log = phi.log;
}

void PartStructT::clear() 
{
	sgl.clear(); sgt.clear(); stl.clear(); elphi = 0; log = "";
}

bool PartStructT::isLocOpt(const stemlist & in_stl, const std::vector<std::vector<int> > & 
   H, const std::vector<std::vector<int> > & C)
{
	FOR(i,0,sz(in_stl.vstems)){
		stem sti = in_stl.vstems[i];
		if(stl.containST(sti)) continue;
		if(!stl.conflictwith(sti) && !sgl.overlap(sti))
			return false;
	}
	return true;
}


void PartStructT::show(std::ostream & fout)
{
	fout<<"PartStruct segments:"<<std::endl;
	FOR(i,0,sz(sgl)) {
		fout<<"["<<sgl[i].lp<<", "<<sgl[i].rp<<"]_"; 
		showETYPE(sgt[i],fout); fout<<std::endl;
	}
	fout<<"PartStruct stacks:";
	FOR(i,0,sz(stl)){
		stl.vstems[i].show();
		fout<<std::endl;
	}
	fout<<"elphi="<<elphi<<std::endl;
}

/******************************************************************************
                            StackFoldAuxiliary
 ******************************************************************************/
StackFoldAuxiliary::StackFoldAuxiliary() 
{
	ctype = ctype_unknown;
	st = stem(-1,-1,-1,-1);
}

StackFoldAuxiliary::StackFoldAuxiliary(const CTYPE & in_ctype, const stem & in_t) 
{
	set(in_ctype, in_t);
}

StackFoldAuxiliary::StackFoldAuxiliary(const StackFoldAuxiliary & in_sfa) {
	set(in_sfa.ctype, in_sfa.st);
}

/*StackFoldAuxiliary::~StackFoldAuxiliary()
{
	ctype = ctype_unknown;
	t.clear();
}*/

void StackFoldAuxiliary::set(const CTYPE & in_ctype, const stem & in_t) {
	ctype = in_ctype;
	st = in_t;
}

void StackFoldAuxiliary::show(std::ostream & fout) {
	showCTYPE(ctype, fout);
	fout<<"\t("<<st.outer.lp<<","<<st.outer.rp<<")~("<<st.inner.lp<<","<<st.inner.rp<<")"<<std::endl;
}
/******************************************************************************
	                           StackFold
 ******************************************************************************/
StackFold::StackFold(const std::string & in_seq)
{
	seq = in_seq;
	stl = stemlist(seq);
	init(sz(stl.s), sz(stl.vstems));	
}


void StackFold::init(const int & n, const int & m)
{
	init_vvi(H,n,n,OOE);
	init_vvi(S,m,m,OOE);
	init_vvi(phi,m,m,OOE);
	init_vi(F,n,OOE);
	init_vi(Fn,n,OOE);
	init_vi(Fd,n,OOE);
	init_vi(F2n,n,OOE);
	init_vi(F2d,n,OOE);
	
	init_vvi(C,n,n,OOE);
	init_vvi(FM,n,n,OOE);
	init_vvi(M1,n,n,OOE);
	init_vvi(M2,n,n,OOE);
	init_vvi(M4,n,n,OOE);
	
	StackFoldAuxiliary sfa;
	init_vsfa(F_sfa,n,sfa);
	init_vsfa(Fn_sfa,n,sfa);
	init_vsfa(Fd_sfa,n,sfa);
	init_vsfa(F2n_sfa,n,sfa);
	init_vsfa(F2d_sfa,n,sfa);
	
	init_vvsfa(C_sfa,n,n,sfa);
	init_vvsfa(FM_sfa,n,n,sfa);
	init_vvsfa(M1_sfa,n,n,sfa);
	init_vvsfa(M2_sfa,n,n,sfa);
	init_vvsfa(M4_sfa,n,n,sfa);
	
	stem s(-1,-1,-1,-1);
	init_vstem(F_S,n,s);
	init_vstem(Fn_S,n,s);
	init_vstem(Fd_S,n,s);	
	init_vstem(F2n_S,n,s);
	init_vstem(F2d_S,n,s);
			
	init_vvstem(FM_LS,n,n,s);
	init_vvstem(FM_RS,n,n,s);
	init_vvstem(M1_LS,n,n,s);
	init_vvstem(M1_RS,n,n,s);
	init_vvstem(M2_LS,n,n,s);
	init_vvstem(M2_RS,n,n,s);
	init_vvstem(M4_LS,n,n,s);
	init_vvstem(M4_RS,n,n,s);
}

//Given all possible putative stacks in an RNA sequence, fold the minimum free \
energy structure
int StackFold::fold(TurnerE & TE) {
	//n=sz(seq), m=sz(stems);
	int n = sz(seq);
	int m = sz(stl.vstems);
	
	//compute H[s]
//	std::cout<<"fill H"<<std::endl;
	fillH(TE);

	//compute S[0..m-1][0..m-1]
//	std::cout<<"fill S"<<std::endl;
	fillS(TE);
	
	//compute fillphi[0..m-1][0..m-1]
//	std::cout<<"fill phi"<<std::endl;
	fillphi(TE);
	
	FOR(k,1,n) {
		int i = 0; int j = i+k;
		if(j<0 || j>=n) break;
		while(i>=0 && i<=n-1 && j>=0 && j<=n-1) {
#ifdef DEBUG
			std::cout<<"(i,j)=("<<i<<","<<j<<"), fill C, FM, M1, M2, M4"<<std::endl;
#endif
			fillC (TE, i, j);
			fillFM(TE, i, j);			
			fillM1(TE, i, j);
			fillM2(TE, i, j);
			fillM4(TE, i, j);
			i++; j++;
		}
	}
	
	FOR(i,0,n) {
		fillF2nF2d(TE, i);
		fillFnFd(TE, i);
		fillF(TE, i);
	}

/*	print2file(C,  "Output/C");
	print2file(FM, "Output/FM");
	print2file(M1, "Output/M1");
	print2file(M2, "Output/M2");	
	print2file(M4, "Output/M4");*/

	return F[n-1];
}

void StackFold::fillH(TurnerE & TE) {
	int n = sz(seq);
	int m = sz(stl.vstems);
	FOR(i,0,m) {
		stem s = stl.vstems[i];
		int sx = s.outer.lp; int sy = s.outer.rp;
		H[sx][sy] = computeH(TE, seq, s);
	}
}

void StackFold::fillS(TurnerE & TE) {
	int m = sz(stl.vstems);
	FOR(i,0,m) {
		stem s = stl.vstems[i];
		FOR(j,i+1,m) {
			stem t = stl.vstems[j];
			if(s.enclose(t) || t.inward_extend(s)) {
				S[i][j] = computeS(TE, seq, s, t);
			}
		}
	}
}

void StackFold::fillphi(TurnerE & TE) {
	int n = sz(seq);
	int m = sz(stl.vstems);
	FOR(i,0,m) {
		stem s = stl.vstems[i];
		FOR(j,i+1,m) {
			stem t = stl.vstems[j];
//   s       t 
// -----...----.....
// |||||   ||||
// -----...----.....
			if(s.enclose(t) || t.inward_extend(s)) {
				phi[i][j] = computephi(TE, seq, s, t);
			}
		}
	}
}

void StackFold::fillC(TurnerE & TE, const int & i, const int & j) {
	std::vector<int> vindx = stl.getstemindice_startfrom_endat(i, j);
	if(sz(vindx)==0) return;
	assert(sz(vindx)==1);
	int indx_s = vindx[0];
	stem s = stl[indx_s];
//	std::cout<<"fillC["<<i<<"]["<<j<<"]"<<std::endl;
	//case C.(1)
	int e = H[i][j];
	C_sfa[i][j].ctype = ctype_1;
	
	vindx = stl.getstemindice_enclosedby_or_inward_extend(indx_s);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl[indx_t];
		stem s2 = getUsedPartOfStem(s, t);
		
		int sx = s2.outer.lp; int sy = s2.outer.rp;
		int sp = s2.inner.lp; int sq = s2.inner.rp;
		int tx = t.outer.lp; int ty = t.outer.rp;
		int tp = t.inner.lp; int tq = t.inner.rp;
		
		//case C.(2)
		//interior loop, neither AU(tx, ty) nor AU(sp, sq) are considered
		int e1 = S[indx_s][indx_t] + phi[indx_s][indx_t] + C[tx][ty];
		if(e1 < e) {
			e = e1;
			C_sfa[i][j].set(ctype_2, t);
		}
#ifdef DEBUG		
		std::cout<<"t="; t.show();
#endif		
		int delta = S[indx_s][indx_t] + TE.Mc + 2*TE.Mphi +
		            TerminalPenalty(TE,seq[sp],seq[sq]) +
		            C[tx][ty] + TerminalPenalty(TE, seq[tx], seq[ty]);
#ifdef DEBUG
		std::cout<<"delta="<<delta<<std::endl;
#endif
		if(ty==sq-1) {
		} else if(ty==sq-2) {
			delta += min(d(TE, seq, tx, ty, ty+1), d(TE, seq, sp, sq, sq-1));
		} else if(ty<=sq-3) {
			delta += d(TE, seq, tx, ty, ty+1) + d(TE, seq, sp, sq, sq-1);
		}
		int i1=sp+1, i2=sp+2, i3=sp+3;
		int j1=tx-1, j2=tx-2, j3=tx-3;
		
		//case C.(3.a)
		int eij = getE(etype_M4, i1, j1);
		int e3a = delta + eij;
#ifdef DEBUG
		std::cout<<"e3a="<<e3a<<", e="<<e<<std::endl;
#endif
		if(e3a < e && eij != OOE) {
			e = e3a;
			C_sfa[i][j].set(ctype_3a, t);
		}
		
		//case C.(3.b)
		stem w = getStem(M4_RS, i1, j2); //M4_RS[i1][j2];
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M4, i1, j2);
			int e3b = delta + eij
			          + min(d(TE, seq, wx, wy, wy+1), d(TE, seq, tx, ty, tx-1));
#ifdef DEBUG
			std::cout<<"e3b="<<e3b<<", e="<<e<<std::endl;
#endif
			if(e3b < e && eij != OOE) {
				e = e3b;
				C_sfa[i][j].set(ctype_3b, t);
			}
  		}
  		
  		//case C.(3.c)
  		w = getStem(M4_RS, i1, j3); //M4_RS[i1][j3];
  		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M4, i1, j3);
			int e3c = delta + eij
			          + d(TE, seq, wx, wy, wy+1) + d(TE, seq, tx, ty, tx-1);
#ifdef DEBUG
			std::cout<<"e3c="<<e3c<<", e="<<e<<std::endl;
#endif
			if(e3c < e && eij != OOE) {
				e = e3c;
				C_sfa[i][j].set(ctype_3c, t);
			}
		}
		
		//case C.(3.d)
		stem v = getStem(M4_LS, i2, j1); //M4_LS[i2][j1];
		if(v.valid()) {
			int vx = v.outer.lp; int vy = v.outer.rp;
			eij = getE(etype_M4, i2, j1);
			int e3d = delta+ eij
			          + min(d(TE, seq, sp, sq, sp+1), d(TE, seq, vx, vy, vx-1));
#ifdef DEBUG
			std::cout<<"e3d="<<e3d<<", e="<<e<<std::endl;
#endif
			if(e3d < e && eij != OOE) {
				e = e3d;
				C_sfa[i][j].set(ctype_3d, t);
			}
		}

		//case C.(3.e)
		v = getStem(M4_LS, i2, j2); //M4_LS[i2][j2];
		w = getStem(M4_RS, i2, j2); //M4_RS[i2][j2];
		if(v.valid() && w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			int vx = v.outer.lp; int vy = v.outer.rp;
			eij = getE(etype_M4, i2, j2);
			int e3e = delta + eij
			          + min(d(TE, seq, sp, sq, sp+1), d(TE, seq, vx, vy, vx-1))
			          + min(d(TE, seq, wx, wy, wy+1), d(TE, seq, tx, ty, tx-1));
#ifdef DEBUG
			std::cout<<"e3e="<<e3e<<", e="<<e<<std::endl;
#endif
			if(e3e < e && eij != OOE) {
				e = e3e;
				C_sfa[i][j].set(ctype_3e, t);
			}
		}
		
		//case C.(3.f)
		v = getStem(M4_LS, i2, j3); //M4_LS[i2][j3];
		w = getStem(M4_RS, i2, j3); //M4_RS[i2][j3];
		if(v.valid() && w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			int vx = v.outer.lp; int vy = v.outer.rp;
			eij = getE(etype_M4, i2, j3);
			int e3f = delta + eij
			          + min(d(TE, seq, sp, sq, sp+1), d(TE, seq, vx, vy, vx-1)) 
			          + d(TE, seq, wx, wy, wy+1) + d(TE, seq, tx, ty, tx-1);
#ifdef DEBUG
			std::cout<<"e3f="<<e3f<<", e="<<e<<std::endl;
#endif
			if(e3f < e && eij != OOE) {
				e = e3f;
				C_sfa[i][j].set(ctype_3f, t);
			}
		}
		
		//case C.(3.g)
		v = getStem(M4_LS, i3, j1); //M4_LS[i3][j1];
		if(v.valid()) {
			int vx = v.outer.lp; int vy = v.outer.rp;
			eij = getE(etype_M4, i3, j1);
			int e3g = delta + eij
			          + d(TE, seq, sp, sq, sp+1) + d(TE, seq, vx, vy, vx-1);
#ifdef DEBUG
			std::cout<<"e3g="<<e3g<<", e="<<e<<std::endl;
#endif
			if(e3g < e && eij != OOE) {
				e = e3g;
				C_sfa[i][j].set(ctype_3g, t);
			}
		}
		
		//case C.(3.h)
		w = getStem(M4_RS, i3, j2); //M4_RS[i3][j2];
		v = getStem(M4_LS, i3, j2); //M4_LS[i3][j2];
		if(v.valid() && w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			int vx = v.outer.lp; int vy = v.outer.rp;
			eij = getE(etype_M4, i3, j2);
			int e3h = delta + eij
			          + d(TE, seq, sp, sq, sp+1) + d(TE, seq, vx, vy, vx-1)
			          + min(d(TE, seq, wx, wy, wy+1), d(TE, seq, tx, ty, tx-1));
#ifdef DEBUG
          std::cout<<"e3h="<<e3h<<", e="<<e<<std::endl;
#endif
			if(e3h < e && eij != OOE) {
				e = e3h;
				C_sfa[i][j].set(ctype_3h, t);
			}
		}
		
		//case C.(3.i)
		w = getStem(FM_RS, i3, j3); //FM_RS[i3][j3];
		v = getStem(FM_LS, i3, j3); //FM_LS[i3][j3];
#ifdef DEBUG
		std::cout<<"C.(3.i)"<<std::endl;
		std::cout<<"FM_LS["<<i3<<"]["<<j3<<"]"; v.show();
		std::cout<<"FM_RS["<<i3<<"]["<<j3<<"]"; w.show();
#endif		
		if(v.valid() && w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			int vx = v.outer.lp; int vy = v.outer.rp;
			eij = getE(etype_FM, i3, j3);
			int e3i = delta + eij
			          + d(TE, seq, sp, sq, sp+1) + d(TE, seq, vx, vy, vx-1)
			          + d(TE, seq, wx, wy, wy+1) + d(TE, seq, tx, ty, tx-1);
#ifdef DEBUG
			std::cout<<"e3i="<<e3i<<", e="<<e<<", delta="<<delta<<", FM[][]"<<FM[i3][j3]<<std::endl;
#endif
//			showd(TE,seq,sp,sq,sp+1);
//			showd(TE,seq,vx,vy,vx-1);
//			showd(TE,seq,wx,wy,wy+1);
//			showd(TE,seq,tx,ty,tx-1);
			if(e3i < e && eij != OOE) {
				e = e3i;
				C_sfa[i][j].set(ctype_3i, t);
			}
		}
	}
	C[i][j] = e;
#ifdef DEBUG
	std::cout<<"C["<<i<<"]["<<j<<"]="<<C[i][j]<<std::endl;
#endif
}

void StackFold::fillFM(TurnerE & TE, const int & i, const int & j) {
	//case FM.(0)
	int e = OOE;
	FM_sfa[i][j].ctype = ctype_0;
	
	std::vector<int> vindx = stl.getstemindice_within(i, j);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl[indx_t];
		int tx = t.outer.lp; int ty = t.outer.rp;
		assert(tx>=i && ty<=j);
		int delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
		//case FM.(1)
		int e1 = delta;
		if(e1 < e) {
			e = e1;
			FM_sfa[i][j].set(ctype_1, t);
			FM_LS[i][j] = t;
			FM_RS[i][j] = t;
		}
		
		//case FM.(2)
		int eij = getE(etype_M1, ty+1, j);
		int e2 = delta + eij;
		if(e2 < e && eij != OOE) {
			e = e2;
			FM_sfa[i][j].set(ctype_2, t);
			FM_LS[i][j] = t;
			FM_RS[i][j] = M1_RS[ty+1][j];
		}
		
		//case FM.(3)
		stem w = getStem(M1_RS, ty+2, j);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M1, ty+2, j);
			int e3 = delta+ eij 
			         + min(d(TE, seq, tx, ty, ty+1), d(TE, seq, wx, wy, wx-1));
			if(e3 < e && eij != OOE) {
				e = e3;
				FM_sfa[i][j].set(ctype_3, t);
				FM_LS[i][j] = t;
				FM_RS[i][j] = M1_RS[ty+2][j];
			}
		}
		
		//case FM.(4)
		w = getStem(FM_RS, ty+3, j);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_FM, ty+3, j);
			int e4 = delta + eij + 
			         d(TE, seq, tx, ty, ty+1) + d(TE, seq, wx, wy, wx-1);
			if(e4 < e && eij != OOE) {
				e = e4;
				FM_sfa[i][j].set(ctype_4, t);
				FM_LS[i][j] = t;
				FM_RS[i][j] = FM_RS[ty+3][j];
			}
		}
	}
	FM[i][j] = e;
#ifdef DEBUG
	std::cout<<"FM["<<i<<"]["<<j<<"]="<<FM[i][j]<<std::endl;
	FM_sfa[i][j].show();
	std::cout<<"LS=";FM_LS[i][j].show();
	std::cout<<"RS=";FM_RS[i][j].show();
	std::cout<<std::endl;
#endif
}

void StackFold::fillM1(TurnerE & TE, const int & i, const int & j) {
	//case M1.(0)
	int e = OOE;
	M1_sfa[i][j].ctype = ctype_0;
	
	std::vector<int> vindx = stl.getstemindice_startfrom_within(i, j);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl[indx_t];
		int tx = t.outer.lp; int ty = t.outer.rp;
		assert(tx==i && ty<=j);
		int delta  = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
		//case M1.(1)
		int e1 = delta;
		if(e1 < e) {
			e = e1;
			M1_sfa[i][j].set(ctype_1, t);
			M1_LS[i][j] = t;
			M1_RS[i][j] = t;
		}
		
		//case M1.(2)
		int eij = getE(etype_M1, ty+1, j);
		int e2 = delta + eij;
		if(e2 < e && eij != OOE) {
			e = e2;
			M1_sfa[i][j].set(ctype_2, t);
			M1_LS[i][j] = t;
			M1_RS[i][j] = M1_RS[ty+1][j];
		}
		
		//case M1.(3)
		stem w = getStem(M1_LS, ty+2, j);//M1_LS[ty+2][j];
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M1, ty+2, j);
			int e3 = delta + eij
				      + min(d(TE, seq, tx, ty, ty+1), d(TE, seq, wx, wy, wx-1));
			if(e3 < e && eij != OOE) {
				e = e3;
				M1_sfa[i][j].set(ctype_3, t);
				M1_LS[i][j] = t;
				M1_RS[i][j] = M1_RS[ty+2][j];
			}
		}
		
		//case M1.(4)
		w = getStem(FM_LS, ty+3, j);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_FM, ty+3, j);
			int e4 = delta + eij 
			         + d(TE, seq, tx, ty, ty+1) + d(TE, seq, wx, wy, wx-1);
			if(e4 < e && eij != OOE) {
				e = e4; 
				M1_sfa[i][j].set(ctype_4, t);
				M1_LS[i][j] = t;
				M1_RS[i][j] = FM_RS[ty+3][j];
			}
		}
	}
	M1[i][j] = e; 
#ifdef DEBUG
	std::cout<<"M1["<<i<<"]["<<j<<"]="<<M1[i][j]<<std::endl;
	M1_sfa[i][j].show();
#endif
}

void StackFold::fillM2(TurnerE & TE, const int & i, const int & j) {
	//case M2.(0);
	int e = OOE;
	M2_sfa[i][j].ctype = ctype_0;
	//M2_LS[i][j] remains to be (-1,-1,-1,-1)
	//M2_RS[i][j] remains to be (-1,-1,-1,-1)
	
	std::vector<int> vindx = stl.getstemindice_endat_within(i, j);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl.vstems[indx_t];
		int tx = t.outer.lp; int ty = t.outer.rp;
		assert(tx>=i && ty==j);
		int delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
		//case M2.(1)
		int e1 = delta;
		if(e1 < e) {
			e = e1; 
			M2_sfa[i][j].set(ctype_1, t);
			M2_LS[i][j] = t;
			M2_RS[i][j] = t;
		}
		
		//case M2.(2)
		int eij = getE(etype_M2, i, tx-1);
		int e2 = delta + eij;
		if(e1 < e && eij != OOE) {
			e = e1; 
			M2_sfa[i][j].set(ctype_2, t);
			M2_LS[i][j] = M2_LS[i][tx-1];
			M2_RS[i][j] = t;
		}
		
		//case M2.(3)
		stem w = getStem(M2_RS, i, tx-2);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M2, i, tx-2);
			int e3 = delta + eij
			         + min(d(TE, seq, wx, wy, wy+1), d(TE, seq, tx, ty, tx-1));
			if(e3 < e && eij != OOE) {
				M2_sfa[i][j].set(ctype_3, t);
				M2_LS[i][j] = M2_LS[i][tx-2];
				M2_RS[i][j] = t;
			}
		}
		
		//case M2.(4)
		w = getStem(FM_LS, i, tx-3);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_FM, i, tx-3);
			int e4 = delta +  eij 
			         + d(TE, seq, wx, wy, wy+1) + d(TE, seq, tx, ty, tx-1);
			if(e4 < e && eij != OOE) {
				e = e4;
				M2_sfa[i][j].set(ctype_4, t);
				M2_LS[i][j] = FM_LS[i][tx-3];
				M2_RS[i][j] = t;
			}
		}
	}
	M2[i][j] = e;
#ifdef DEBUG
	std::cout<<"M2["<<i<<"]["<<j<<"]="<<M2[i][j]<<std::endl;
	M2_sfa[i][j].show();
#endif
}

void StackFold::fillM4(TurnerE & TE, const int & i, const int & j) {
	//case M4.(0)
	int e = OOE;
	M4_sfa[i][j].ctype = ctype_0;
	//M4_LS[i][j] remains to be (-1,-1,-1,-1)
	//M4_RS[i][j] remains to be (-1,-1,-1,-1)
	
	std::vector<int> vindx = stl.getstemindice_startfrom_within(i, j);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl.vstems[indx_t];
		int tx = t.outer.lp; int ty = t.outer.rp;
		assert(tx==i && ty<=j);
		int delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
		int e1 = delta;
		//case M4.(1)
		if(e1 < e && ty==j) {
			e = e1;
			M4_sfa[i][j].set(ctype_1, t);
			M4_LS[i][j] = t;
			M4_RS[i][j] = t;
		}
		
		//case M4.(2)
		int eij = getE(etype_M4, ty+1, j);
		int e2 = delta + eij;
		if(e2 < e && eij != OOE) {
			e = e2;
			M4_sfa[i][j].set(ctype_2, t);
			M4_LS[i][j] = t;
			M4_RS[i][j] = M4_RS[ty+1][j];
		}
		
		//case M4.(3)
		stem w = getStem(M4_LS, ty+2, j);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M4, ty+2, j);
			int e3 = delta + min(d(TE, seq, tx, ty, ty+1), d(TE, seq, wx, wy, wx-1)) + eij;
			if(e3 < e && eij != OOE) {
				e = e3;
				M4_sfa[i][j].set(ctype_3, t);
				M4_LS[i][j] = t;
				M4_RS[i][j] = M4_RS[ty+2][j];
			}
		}
				
		//case M4.(4)
		w = getStem(M2_LS, ty+3, j);
		if(w.valid()) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			eij = getE(etype_M2, ty+3, j);
			int e4 = delta + d(TE, seq, tx, ty, ty+1) + d(TE, seq, wx, wy, wx-1) + eij;
			if(e4 < e && eij != OOE) {
				e = e4;
				M4_sfa[i][j].set(ctype_4, t);
				M4_LS[i][j] = t;
				M4_RS[i][j] = M2_RS[ty+3][j];
			}
		}
	}
	M4[i][j] = e;
#ifdef DEBUG
	std::cout<<"M4["<<i<<"]["<<j<<"]="<<M4[i][j]<<std::endl;
	M4_sfa[i][j].show();
#endif
}

void StackFold::fillFnFd(TurnerE & TE, const int & j) {
	//compute Fn[j] and Fd[j]
	//case FnFd.(0)
	if(j==0) {
		Fn[j] = Fd[j] = 0;
		Fn_sfa[j].ctype = Fd_sfa[j].ctype = ctype_0;
		return;
	}
	
	stem w = getStem(Fn_S, j-1);
	//case FnFd.(1)
	int efn = Fn[j-1];
	int efd = Fd[j-1];
	Fn_S[j] = Fn_S[j-1];
	Fd_S[j] = Fd_S[j-1];
	Fn_sfa[j].ctype = Fd_sfa[j].ctype = ctype_1;
	
	
	std::vector<int> vindx = stl.getstemindice_endat(j);
	FOR(k,0,sz(vindx)) {
		int indx_s = vindx[k];
		stem s = stl[indx_s];
		int sx = s.outer.lp; int sy = s.outer.rp;
		int delta = C[sx][sy] + TerminalPenalty(TE, seq[sx], seq[sy]);
		int d0 = d(TE,seq,sx,sy,sy+1);
		
		//case FnFd.(2)
		stem t = getStem(F2n_S, sx-1);
		int ej = getE(etype_F2n, sx-1);
		if(t.valid() && ej != OOE) {
			int tx = t.outer.lp; int ty = t.outer.rp;
			int efn2 = ej + delta;
			int efd2 = efn2 + d0;
			if(efn2 < efn) {
				efn = efn2;
				Fn_S[j] = s;
				Fn_sfa[j].set(ctype_2, s);
			}
			if(efd2 < efd) {
				efd = efd2;
				Fd_S[j] = s;
				Fd_sfa[j].set(ctype_2, s);
			}
		}
		
		//case FnFd.(3)
		stem w = getStem(F2d_S, sx-2);
		ej = getE(etype_F2d, sx-2);
		if(w.valid() && ej != OOE) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			int efn3 = ej + delta - d(TE,seq,wx,wy,wy+1) + 
			           min(d(TE,seq,wx,wy,wy+1), d(TE,seq,sx,sy,sx-1));
			int efd3 = efn3 +  d0;
			if(efn3 < efn) {
				efn = efn3;
				Fn_S[j] = s;
				Fn_sfa[j].set(ctype_3, s);
			}
			if(efd3 < efd) {
				efd = efd3;
				Fd_S[j] = s;
				Fd_sfa[j].set(ctype_3, s);
			}
		}
		
		//case FnFd.(4)
		stem v = getStem(Fd_S, sx-3);
		ej = getE(etype_Fd, sx-3);
		if(v.valid() && ej!=OOE) {
			int efn4 = ej + delta + d(TE,seq,sx,sy,sx-1);
			int efd4 = efn4 + d0;
			if(efn4 < efn) {
				efn = efn4;
				Fn_S[j] = s;
				Fn_sfa[j].set(ctype_4, s);
			}
			if(efd4 < efd) {
				efd = efd4;
				Fd_S[j] = s;
				Fd_sfa[j].set(ctype_4, s);
			}
		}
		
		//case FnFd.(5)
		if(sx==0) {
			//case FnFd.(5.a)
			int efn5a = delta;
			int efd5a = efn5a + d0;
			if(efn5a < efn) {
				efn = efn5a;
				Fn_S[j] = s;
				Fn_sfa[j].set(ctype_5a, s);
			}
			if(efd5a < efd) {
				efd = efd5a;
				Fd_S[j] = s;
				Fd_sfa[j].set(ctype_5a, s);
			}
		} else if(sx>=1) {
			//case FnFd.(5.b)
			int efn5b = delta + d(TE,seq,sx,sy,sx-1);
			int efd5b = efn5b + d0;
			if(efn5b < efn) {
				efn = efn5b;
				Fn_S[j] = s;
				Fn_sfa[j].set(ctype_5b, s);
			}
			if(efd5b < efd) {
				efd = efd5b;
				Fd_S[j] = s;
				Fd_sfa[j].set(ctype_5b, s);
			}
		}
	}
	
	Fn[j] = efn;
	Fd[j] = efd;

#ifdef DEBUG	
	std::cout<<"Fn["<<j<<"]="<<Fn[j]<<std::endl;
	std::cout<<"Fd["<<j<<"]="<<Fd[j]<<std::endl;
#endif
}

void StackFold::fillF(TurnerE & TE, const int & j) {
	if(j==0) {
		F[j] = 0;
		F_sfa[j].ctype = ctype_0;
		return;
	}
	
	int ef = OOE;
	int efn = Fn[j];
	int efd = Fd[j];
	stem w = getStem(Fn_S, j);
	int wx = w.outer.lp; int wy = w.outer.rp;
	
	stem v = getStem(Fd_S, j);
	int vx = v.outer.lp; int vy = v.outer.rp;
	
	if(w.valid()) {
		if(wy==j) {
			//	efn += 0;
		} else {
			efn += d(TE,seq,wx,wy,wy+1);
		}
	}
	if(v.valid()) {
		if(vy==j) {
			efd -= d(TE,seq,vx,vy,vy+1);
		} else {
			// efd += 0;
		}
	}	
	if(efn < efd) {
		F[j] = efn;
		F_S[j] = w;
		F_sfa[j].set(ctype_1, w);
	} else {
		F[j] = efd;
		F_S[j] = v;
		F_sfa[j].set(ctype_2, v);
	}
#ifdef DEBUG
	std::cout<<"F["<<j<<"]="<<F[j]<<std::endl;
#endif
}


void StackFold::fillF2nF2d(TurnerE & TE, const int & j) {
	//case (0)
	int ef2n = OOE; int ef2d = OOE;
	F2n_sfa[j].ctype = F2d_sfa[j].ctype = ctype_0;
	
	std::vector<int> vindx = stl.getstemindice_endat(j);
	FOR(k,0,sz(vindx)) {
		int indx_s = vindx[k];
		stem s = stl[indx_s];
		int sx = s.outer.lp; int sy = s.outer.rp;
		
		int delta = C[sx][sy] + TerminalPenalty(TE,seq[sx],seq[sy]);
		int d0 = d(TE,seq,sx,sy,sy+1);
		
		stem t = getStem(F2n_S, sx-1);
		int ej = getE(etype_F2n, sx-1);
		//case (1)
		if(t.valid() && ej!=OOE) {
			int ef2n1 = delta + ej;
			int ef2d1 = ef2n1 + d0;
			if(ef2n1 < ef2n) {
				ef2n = ef2n1;
				F2n_S[j] = s;
				F2n_sfa[j].set(ctype_1, s);
			}
			if(ef2d1 < ef2d) {
				ef2d = ef2d1;
				F2d_S[j] = s;
				F2d_sfa[j].set(ctype_1, s);
			}
		}
		
		//case (2)
		stem w = getStem(F2d_S, sx-2);
		ej = getE(etype_F2d, sx-2);
		if(w.valid() && ej!=OOE) {
			int wx = w.outer.lp; int wy = w.outer.rp;
			int ef2n2 = ej + delta - d(TE,seq,wx,wy,wy+1) + 
			            min(d(TE,seq,wx,wy,wy+1), d(TE,seq,sx,sy,sx-1));
			int ef2d2 = ef2n2 + d0;
			if(ef2n2 < ef2n) {
				ef2n = ef2n2;
				F2n_S[j] = s;
				F2n_sfa[j].set(ctype_2, s);
			}
			if(ef2d2 < ef2d) {
				ef2d = ef2d2;
				F2d_S[j] = s;
				F2d_sfa[j].set(ctype_2, s);
			}
		}
		
		//case (3)
		stem v = getStem(Fd_S, sx-3);
		ej = getE(etype_Fd, sx-3);
		if(v.valid() && ej!=OOE) {
			int ef2n3 = delta + ej + d(TE,seq,sx,sy,sx-1);
			int ef2d3 = ef2n3 + d0;
			if(ef2n3 < ef2n) {
				ef2n = ef2n3;
				F2n_S[j] = s; 
				F2n_sfa[j].set(ctype_3, s);
			}
			if(ef2d3 < ef2d) {
				ef2d = ef2d3;
				F2d_S[j] = s;
				F2d_sfa[j].set(ctype_3, s);
			}
		}
		
		//case (4)
		if(sx==0) {
			//case (4.a)
			int ef2n4a = delta;
			int ef2d4a = ef2n4a + d0;
			if(ef2n4a < ef2n) {
				ef2n = ef2n4a;
				F2n_S[j] = s;
				F2n_sfa[j].set(ctype_4a, s);
			}
			if(ef2d4a < ef2d) {
				ef2d = ef2d4a;
				F2d_S[j] = s;
				F2d_sfa[j].set(ctype_4a,s);
			}
		} else if(sx>=1) {
			//case (4.b)
			int ef2n4b = delta + d(TE,seq,sx,sy,sx-1);
			int ef2d4b = ef2n4b + d0;
			if(ef2n4b < ef2n) {
				ef2n = ef2n4b;
				F2n_S[j] = s;
				F2n_sfa[j].set(ctype_4b, s);
			}
			if(ef2d4b < ef2d) {
				ef2d = ef2d4b;
				F2d_S[j] = s;
				F2d_sfa[j].set(ctype_4b,s);
			}
		}
	}
	
	F2n[j] = ef2n;
	F2d[j] = ef2d;
#ifdef DEBUG
	std::cout<<"F2n["<<j<<"]="<<F2n[j]<<std::endl;
	std::cout<<"F2d["<<j<<"]="<<F2d[j]<<std::endl;
#endif
}


std::string StackFold::traceback(TurnerE & TE) {
	int n = sz(seq);
	std::string ret = "";
	FOR(i,0,n) ret += '.';
	
	std::vector<stem> vst = tracebackF(n-1, TE);
//	std::cout<<"sz(traceback)="<<sz(vst)<<std::endl;
	FOR(i,0,sz(vst)) {
		stem s = vst[i];
		assert(s.valid());
		FOR(j,i+1,sz(vst)) {
			stem t = vst[j];
//			std::cout<<"comparing s and t"<<std::endl<<"s=";s.show();std::cout<<"t=";t.show();
			assert(!s.conflictwith(t));
		}
	}
	
	FOR(i,0,sz(vst)) {
		stem s = vst[i];
		FOR(j,s.outer.lp, s.inner.lp+1) {
			ret[j] = '(';
		}
		FOR(j,s.inner.rp, s.outer.rp+1) {
			ret[j] = ')';
		}
	}
	return ret;
}

std::vector<stem> StackFold::tracebackF(const int & j, TurnerE & TE) {
//	std::cout<<"traceback F["<<j<<"]="<<F[j]<<std::endl;
	int n = sz(seq);
	assert(j>=0 && j<n);
	std::vector<stem> ret;
	
	CTYPE ctypej = F_sfa[j].ctype;
	stem s = F_sfa[j].st;
	int sx = s.outer.lp; int sy = s.outer.rp;
	
	switch(ctypej) {
		case ctype_0: {
			//do nothing
			break;
		} case ctype_1: {
			ret = tracebackFn(j,TE);
			break;
		} case ctype_2: {
			ret = tracebackFd(j,TE);
			break;
		}
	}
	return ret;
}

std::vector<stem> StackFold::tracebackFn(const int & j, TurnerE & TE) {
//	std::cout<<"traceback Fn["<<j<<"]="<<Fn[j]<<std::endl;
	
	int n = sz(seq);
	assert(j>=0 && j<n);
	CTYPE ctypej = Fn_sfa[j].ctype;
	stem s = Fn_sfa[j].st;
	int sx = s.outer.lp; int sy = s.outer.rp;
	
	std::vector<stem> vstC;
	std::vector<stem> vstF;
	
	switch(ctypej) {
		case ctype_0: {
			//do nothing
			break;
		} case ctype_1: {
			vstF = tracebackFn(j-1,TE);
			break;
		} case ctype_2: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackF2n(sx-1,TE);
			break;
		} case ctype_3: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackF2d(sx-2,TE);
			break;
		} case ctype_4: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackFd(sx-3,TE);
			break;
		} case ctype_5a: {
		} case ctype_5b: {
			vstC = tracebackC(sx,sy,TE);
			break;
		} default: {
			assert(false);
		}
	}
	
	std::vector<stem> ret;
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstF.begin(), vstF.end());
	return ret;
}

std::vector<stem> StackFold::tracebackFd(const int & j, TurnerE & TE) {
//	std::cout<<"traceback Fd["<<j<<"]="<<Fd[j]<<std::endl;
	int n = sz(seq);
	assert(j>=0 && j<n);
	CTYPE ctypej = Fd_sfa[j].ctype;
	stem s = Fd_sfa[j].st;
	int sx = s.outer.lp; int sy = s.outer.rp;
	
	std::vector<stem> vstC;
	std::vector<stem> vstF;
	
	switch(ctypej) {
		case ctype_0: {
			//do nothing
			break;
		} case ctype_1: {
			vstF = tracebackFd(j-1,TE);
			break;
		} case ctype_2: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackF2n(sx-1,TE);
			break;
		} case ctype_3: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackF2d(sx-2,TE);
			break;
		} case ctype_4: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackFd(sx-3,TE);
			break;
		} case ctype_5a: {
		} case ctype_5b: {
			vstC = tracebackC(sx,sy,TE);
			break;
		} default: {
			assert(false);
		}
	}
	
	std::vector<stem> ret;
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstF.begin(), vstF.end());
	return ret;
}

std::vector<stem> StackFold::tracebackFnFd(const int & j, TurnerE & TE, StackFoldAuxiliary & sfa) {

}

std::vector<stem> StackFold::tracebackF2n(const int & j, TurnerE & TE) {
	int n = sz(seq);
	assert(j>=0 && j<n);
//	std::cout<<"traceback F2n["<<j<<"]="<<F2n[j]<<std::endl;
//	F2n_sfa[j].show();
	return tracebackF2nF2d(j, TE, F2n_sfa[j]);
}

std::vector<stem> StackFold::tracebackF2d(const int & j, TurnerE & TE) {
	int n = sz(seq);
	assert(j>=0 && j<n);
//	std::cout<<"traceback F2d["<<j<<"]="<<F2d[j]<<std::endl;
//	F2d_sfa[j].show();
	return tracebackF2nF2d(j, TE, F2d_sfa[j]);
}

std::vector<stem> StackFold::tracebackF2nF2d(const int & j, TurnerE & TE, StackFoldAuxiliary & sfa) {
	int n = sz(seq);
	assert(j>=0 && j<n);
	CTYPE ctypej = sfa.ctype;
	stem s = sfa.st;
	int sx = s.outer.lp; int sy = s.outer.rp;
	
	std::vector<stem> vstC;
	std::vector<stem> vstF;
	
	switch(ctypej) {
		case ctype_0: {
			assert(false);
			break;
		} case ctype_1: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackF2n(sx-1,TE);
			break;
		} case ctype_2: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackF2d(sx-2,TE);
			break;
		} case ctype_3: {
			vstC = tracebackC(sx,sy,TE);
			vstF = tracebackFd(sx-3,TE);
			break;
		} case ctype_4a: {
		} case ctype_4b: {
			vstC = tracebackC(sx,sy,TE);
			break;
		} default: {
			assert(false);
		}
	}
	
	std::vector<stem> ret;
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstF.begin(), vstF.end());
	return ret;
}


std::vector<stem> StackFold::tracebackC(const int & i, const int & j, TurnerE & TE) {
	int n = sz(seq);
//	std::cout<<"traceback C["<<i<<"]["<<j<<"]="<<C[i][j]<<std::endl;
	assert(i>=0 && i<n && j>=0 && j<n && i<j);
	std::vector<stem> ret;
	std::vector<stem> vstC;
	std::vector<stem> vstM;
	
	//get stem s
	std::vector<int> vindx = stl.getstemindice_startfrom_endat(i, j);
	assert(sz(vindx)==1);
	stem s = stl.vstems[vindx[0]];
	
	//get stem t
	CTYPE ctypeij = C_sfa[i][j].ctype;
	stem t = C_sfa[i][j].st;
	int tx = t.outer.lp; int ty = t.outer.rp;
	int tp = t.inner.lp; int tq = t.inner.rp;
	//get stem s2
	stem s2 = s;
	if(ctypeij != ctype_1) 
		s2 = getUsedPartOfStem(s, t);
	
	int sx = s2.outer.lp; int sy = s2.outer.rp;
	int sp = s2.inner.lp; int sq = s2.inner.rp;
	ret.PB(s2);
	
	//tracebackC(t)
	if(ctypeij != ctype_1) {
		vstC = tracebackC(tx,ty, TE);
	}
	
	//get i1, i2, i3, j1, j2, j3
	int i1 = sp+1, i2 = sp+2, i3 = sp+3;
	int j1 = tx-1, j2 = tx-2, j3 = tx-3;
	switch(ctypeij) {
		case ctype_1: {
			//do nothing
			break;
		} case ctype_2: {
			break;
		} case ctype_3a: {
			vstM = tracebackM4(i1,j1, TE);
			break;
		} case ctype_3b: {
			vstM = tracebackM4(i1,j2, TE);
			break;
		} case ctype_3c: {
			vstM = tracebackM4(i1,j3, TE);
			break;
		} case ctype_3d: {
			vstM = tracebackM4(i2,j1, TE);
			break;
		} case ctype_3e: {
			vstM = tracebackM4(i2,j2, TE);
			break;
		} case ctype_3f: {
			vstM = tracebackM4(i2,j3, TE);
			break;
		} case ctype_3g: {
			vstM = tracebackM4(i3,j1, TE);
			break;
		} case ctype_3h: {
			vstM = tracebackM4(i3,j2, TE);
			break;
		} case ctype_3i: {
			vstM = tracebackFM(i3,j3, TE);
			break;
		} default :{
			assert(false);
		}
	}
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstM.begin(), vstM.end());
#ifdef DEBUG
	int indx_s = getindx(stl.vstems, s);
	int indx_t = getindx(stl.vstems, t);
	
	int delta = 0;
	
	if(ctypeij!=ctype_1) {
		delta = S[indx_s][indx_t] + TE.Mc + 2*TE.Mphi + 
		        TerminalPenalty(TE, seq[sp], seq[sq]) +
		        C[tx][ty] + TerminalPenalty(TE,seq[tx],seq[ty]);
		if(ty==sq-1)	delta += 0;
		else if(ty==sq-2) delta += min(d(TE, seq, tx, ty, ty+1), d(TE, seq, sp, sq, sq-1));
		else if(ty<=sq-3) delta += d(TE,seq,tx,ty,ty+1) + d(TE,seq,sp,sq,sq-1);
	}

	switch(ctypeij) {
		case ctype_1: {
			assert(C[i][j]==H[i][j]);
			break;
		} case ctype_2: {
			assert(t.valid());
			assert(C[i][j]==S[indx_s][indx_t] + phi[indx_s][indx_t] + C[tx][ty]);
			break;
		} case ctype_3a: {
			assert(t.valid());
			assert(C[i][j]==delta + M4[i1][j1]);
			break;
		} case ctype_3b: {
			stem w = M4_RS[i1][j2];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(w.valid());
			assert(C[i][j]==delta + M4[i1][j2] + min(d(TE,seq,wx,wy,wy+1), d(TE,seq,tx,ty,tx-1)));
			break;
		} case ctype_3c: {
			stem w = M4_RS[i1][j3];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(w.valid());
			assert(C[i][j]==delta + M4[i1][j3] + d(TE,seq,wx,wy,wy+1) + d(TE,seq,tx,ty,tx-1));
			break;
		} case ctype_3d: {
			stem v =  M4_LS[i2][j1];
			int vx = v.outer.lp; int vy = v.outer.rp;
			assert(v.valid());
			assert(C[i][j]==delta + M4[i2][j1] + min(d(TE,seq,sp,sq,sp+1), d(TE,seq,vx,vy,vx-1)));
			break;
		} case ctype_3e: {
			stem v = M4_LS[i2][j2];
			int vx = v.outer.lp; int vy = v.outer.rp;
			stem w = M4_RS[i2][j2];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(v.valid() && w.valid());
			assert(C[i][j]==delta + M4[i2][j2] + min(d(TE,seq,sp,sq,sp+1), d(TE,seq,vx,vy,vx-1)) + min(d(TE,seq,wx,wy,wy+1), d(TE,seq,tx,ty,tx-1)));
			break;
		} case ctype_3f: {
			stem v = M4_LS[i2][j3];
			int vx = v.outer.lp; int vy = v.outer.rp;
			stem w = M4_RS[i2][j3];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(v.valid() && w.valid());
			assert(C[i][j]==delta + M4[i2][j3] + min(d(TE,seq,sp,sq,sp+1), d(TE,seq,vx,vy,vx-1)) + d(TE,seq,wx,wy,wy+1) + d(TE,seq,tx,ty,tx-1));
			break;
		} case ctype_3g: {
			stem v = M4_LS[i3][j1];
			int vx = v.outer.lp; int vy = v.outer.rp;
			stem w = M4_RS[i3][j1];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(v.valid() && w.valid());
			assert(C[i][j]==delta + M4[i3][j1] + d(TE,seq,sp,sq,sp+1) + d(TE,seq,vx,vy,vx-1));
			break;
		} case ctype_3h: {
			stem v = M4_LS[i3][j2];
			int vx = v.outer.lp; int vy = v.outer.rp;
			stem w = M4_RS[i3][j2];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(v.valid() && w.valid());
			assert(C[i][j]==delta + M4[i3][j2] + d(TE,seq,sp,sq,sp+1) + d(TE,seq,vx,vy,vx-1) + min(d(TE,seq,wx,wy,wy+1), d(TE,seq,tx,ty,tx-1)));
			break;
		} case ctype_3i: {
			stem v = FM_LS[i3][j3];
			int vx = v.outer.lp; int vy = v.outer.rp;
			stem w = FM_RS[i3][j3];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(v.valid() && w.valid());
			assert(C[i][j]==delta + FM[i3][j3] + d(TE,seq,sp,sq,sp+1) + d(TE,seq,vx,vy,vx-1) + d(TE,seq,wx,wy,wy+1) + d(TE,seq,tx,ty,tx-1));
			break;
		} default :{
			assert(false);
		}
	}
#endif
	return ret;
}

std::vector<stem> StackFold::tracebackFM(const int & i, const int & j, TurnerE & TE) {
	int n = sz(seq);
	assert(i>=0 && i<n && j>=0 && j<n && i<j);
	std::vector<stem> ret;
	std::vector<stem> vstC; 
	std::vector<stem> vstM;
	
//	std::cout<<"FM["<<i<<"]["<<j<<"]="<<FM[i][j]<<std::endl;
	CTYPE ctypeij = FM_sfa[i][j].ctype;
	stem t = FM_sfa[i][j].st;
	int tx = t.outer.lp; int ty = t.outer.rp;
	switch(ctypeij){
		case ctype_0: {
			assert(false);
			break; 
		} case ctype_1: {
			vstC = tracebackC(tx,ty, TE);
			break; 
		} case ctype_2: {
			vstC = tracebackC(tx,ty, TE);
			vstM = tracebackM1(ty+1,j, TE);
			break;
		} case ctype_3: {
			vstC = tracebackC(tx,ty, TE);
			vstM = tracebackM1(ty+2,j, TE);
			break; 
		} case ctype_4: {
			vstC = tracebackC(tx,ty, TE);
			vstM = tracebackFM(ty+3,j, TE);
			break; 
		} default:{
			assert(false);
		}
	}
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstM.begin(), vstM.end());
#ifdef DEBUG
	int delta = 0;
	if(ctypeij!=ctype_0) {
		delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
		assert(t.valid());
	}
	switch(ctypeij){
		case ctype_0: {
			assert(FM[i][j]==OOE);
			break; 
		} case ctype_1: {
			assert(FM[i][j]==delta);
			break; 
		} case ctype_2: {
			assert(FM[i][j]==delta + M1[ty+1][j]);
			break;
		} case ctype_3: {
			stem w = M1_RS[ty+2][j];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(FM[i][j]==delta + M1[ty+2][j] + min(d(TE,seq,tx,ty,ty+1), d(TE,seq,wx,wy,wx-1)));
			break; 
		} case ctype_4: {
			stem w = FM_RS[ty+3][j];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(FM[i][j]==delta + FM[ty+3][j] + d(TE,seq,tx,ty,ty+1) + d(TE,seq,wx,wy,wx-1));
			break; 
		} default:{
			assert(false);
		}
	}
#endif	
	return ret;
}

std::vector<stem> StackFold::tracebackM1(const int & i, const int & j, TurnerE & TE) {
	int n = sz(seq);
	assert(i>=0 && i<n && j>=0 && j<n && i<j);
	std::vector<stem> ret;
	
	CTYPE ctypeij = M1_sfa[i][j].ctype;
	stem t = M1_sfa[i][j].st;
	int tx = t.outer.lp; int ty = t.outer.rp;
	
//	std::cout<<"M1["<<i<<"]["<<j<<"]="<<M1[i][j]<<std::endl;
	std::vector<stem> vstC;
	std::vector<stem> vstM;
	
	switch(ctypeij){
		case ctype_0: {
			assert(false);
			break; 
		} case ctype_1: {
			vstC = tracebackC(tx,ty, TE);
			break; 
		} case ctype_2: {
			vstC = tracebackC(tx,ty, TE);
			vstM = tracebackM1(ty+1,j, TE);
			break; 
		} case ctype_3: {
			vstC = tracebackC(tx,ty, TE);
			vstM = tracebackM1(ty+2,j, TE);
			break; 
		} case ctype_4: {
			vstC = tracebackC(tx,ty, TE);
			vstM = tracebackFM(ty+3,j, TE);
			break; 
		} default:{
			assert(false);
		}
	}
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstM.begin(), vstM.end());
#ifdef DEBUG
	int delta = 0;
	if(ctypeij!=ctype_0) {
		assert(t.valid());
		delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
	}
	switch(ctypeij){
		case ctype_0: {
			assert(M1[i][j]==OOE);
			break; 
		} case ctype_1: {
			assert(M1[i][j]==delta);
			break; 
		} case ctype_2: {
			assert(M1[i][j]==delta+M1[ty+1][j]);
			break; 
		} case ctype_3: {
			stem w = M1_LS[ty+2][j];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(M1[i][j]==delta + M1[ty+2][j] + min(d(TE,seq,tx,ty,ty+1), d(TE,seq,wx,wy,wx-1)));
			break; 
		} case ctype_4: {
			stem w = M1_LS[ty+3][j];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(M1[i][j]==delta + FM[ty+3][j] + d(TE,seq,tx,ty,ty+1) + d(TE,seq,wx,wy,wx-1));
			break; 
		} default:{
			assert(false);
		}
	}
#endif
	return ret;
}

std::vector<stem> StackFold::tracebackM2(const int & i, const int & j, TurnerE & TE) {
	int n = sz(seq);
	assert(i>=0 && i<n && j>=0 && j<n && i<j);
	std::vector<stem> ret;
	CTYPE ctypeij = M2_sfa[i][j].ctype;
	stem t = M2_sfa[i][j].st;
	int tx = t.outer.lp; int ty = t.outer.rp;
	
//	std::cout<<"M2["<<i<<"]["<<j<<"]="<<M2[i][j]<<std::endl;
	std::vector<stem> vstC;
	std::vector<stem> vstM;
	switch(ctypeij){
		case ctype_0: {
			assert(false);
			break; 
		} case ctype_1: {
			vstC = tracebackC(tx, ty, TE);
			break; 
		} case ctype_2: {
			vstC = tracebackC(tx, ty, TE);
			vstM = tracebackM2(i,tx-1, TE);
			break; 
		} case ctype_3: {
			vstC = tracebackC(tx, ty, TE);
			vstM = tracebackM2(i,tx-2, TE);
			break;
		} case ctype_4: {
			vstC = tracebackC(tx, ty, TE);
			vstM = tracebackFM(i,tx-3, TE);
			break; 
		} default:{
			assert(false);
		}
	}
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstM.begin(), vstM.end());
#ifdef DEBUG
	int delta = 0;
	if(ctypeij!=ctype_0) {
		delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE,seq[tx],seq[ty]);
		assert(t.valid());
	}
	switch(ctypeij){
		case ctype_0: {
			assert(M2[i][j]==OOE);
			break; 
		} case ctype_1: {
			assert(M2[i][j]==delta);
			break; 
		} case ctype_2: {
			assert(M2[i][j]==delta + M2[i][tx-1]);
			break;
		} case ctype_3: {
			stem w = M2_RS[i][tx-2];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(M2[i][j]==delta + M2[i][tx-2]+ min(d(TE,seq,wx,wy,wy+1), d(TE,seq,tx,ty,tx-1)));
			break;
		} case ctype_4: {
			stem w = M2_RS[i][tx-3];
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(M2[i][j]==delta + M2[i][tx-3]+ d(TE,seq,wx,wy,wy+1) + d(TE,seq,tx,ty,tx-1));
			break; 
		} default:{
			assert(false);
		}
	}
#endif
	return ret;
}

std::vector<stem> StackFold::tracebackM4(const int & i, const int & j, TurnerE & TE) {
	int n = sz(seq);
	assert(i>=0 && i<n && j>=0 && j<n && i<j);
	std::vector<stem> ret;
	CTYPE ctypeij = M4_sfa[i][j].ctype;
	stem t = M4_sfa[i][j].st;
	int tx = t.outer.lp; int ty = t.outer.rp;
	
//	std::cout<<"M4["<<i<<"]["<<j<<"]="<<M4[i][j]<<std::endl;
	
	std::vector<stem> vstC;
	std::vector<stem> vstM;
	
	switch(ctypeij){
		case ctype_0: {
			assert(false);
			break; 
		} case ctype_1: {
			vstC = tracebackC(tx, ty, TE);
			break; 
		} case ctype_2: {
			vstC = tracebackC(tx, ty, TE);
			vstM = tracebackM4(ty+1,j, TE);
			break; 
		} case ctype_3: {
			vstC = tracebackC(tx, ty, TE);
			vstM = tracebackM4(ty+2,j, TE);
			break; 
		} case ctype_4: {
			vstC = tracebackC(tx, ty, TE);
			vstM = tracebackM2(ty+3,j, TE);
			break; 
		} default:{
			assert(false);
		}
	}
	ret.insert(ret.begin(), vstC.begin(), vstC.end());
	ret.insert(ret.begin(), vstM.begin(), vstM.end());
#ifdef DEBUG
	int delta = 0;
	if(ctypeij!=ctype_0) {
		delta = C[tx][ty] + TE.Mphi + TerminalPenalty(TE, seq[tx], seq[ty]);
		assert(t.valid());
	}
	switch(ctypeij){
		case ctype_0: {
			assert(M4[i][j]==OOE);
			break; 
		} case ctype_1: {
			assert(M4[i][j]==delta);
			break; 
		} case ctype_2: {
			assert(M4[i][j]==delta + M4[ty+1][j]);
			break; 
		} case ctype_3: {
			stem w = M4_LS[ty+2][j]; 
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(M4[i][j]==delta + min(d(TE,seq,tx,ty,ty+1), d(TE,seq,wx,wy,wx-1)) + M4[ty+2][j]);
			break; 
		} case ctype_4: {
			stem w = M2_LS[ty+3][j]; 
			int wx = w.outer.lp; int wy = w.outer.rp;
			assert(M4[i][j]==delta + d(TE,seq,tx,ty,ty+1) + d(TE,seq,wx,wy,wx-1) + M2[ty+3][j]);
			break; 
		} default:{
			assert(false);
		}
	}
#endif
	return ret;
}
	

int StackFold::getE(const ETYPE & etype, const int & i, const int & j) {
	int ret = OOE;
	switch(etype) {
		case etype_F: { 
			if(!(i<0 || i>=sz(seq))) ret = F[i];
			break;
		} case etype_Fn: {
			if(!(i<0 || i>=sz(seq))) ret = Fn[i];
			break;
		} case etype_Fd: {
			if(!(i<0 || i>=sz(seq))) ret = Fd[i];
			break;
		} case etype_F2n: {
			if(!(i<0 || i>=sz(seq))) ret = F2n[i];
			break;
		} case etype_F2d: {
			if(!(i<0 || i>=sz(seq))) ret = F2d[i];
			break;
		} case etype_C : {
			if(!(i<0 || i>=sz(seq) || j<0 || j>=sz(seq) || (i>j) || (j<i+4)))
				ret = C[i][j];
			break;
		} case etype_FM: {
			if(!(i<0 || i>=sz(seq) || j<0 || j>=sz(seq) || (i>j)))
				ret = FM[i][j];
			break;
		} case etype_M1: {
			if(!(i<0 || i>=sz(seq) || j<0 || j>=sz(seq) || (i>j)))
				ret = M1[i][j];
			break;
		} case etype_M2: {
			if(!(i<0 || i>=sz(seq) || j<0 || j>=sz(seq) || (i>j)))
				ret = M2[i][j];
			break;
		} case etype_M4: {
			if(!(i<0 || i>=sz(seq) || j<0 || j>=sz(seq) || (i>j)))
				ret = M4[i][j];
			break;
		} default: {
			showETYPE(etype);
			assert(false);
		}
	}
	ret = min(OOE, ret);
	return ret;
}

stem StackFold::getStem(const std::vector<std::vector<stem> > & vvstem, const int & i, const int & j) {
	stem ret(-1,-1,-1,-1);
	if(!(i<0 || i>=sz(seq) || j<0 || j>=sz(seq) || (i>j))) {
		ret = vvstem[i][j];
	}
	return ret;
}

stem StackFold::getStem(const std::vector<stem> & vstem, const int & j) {
	stem ret(-1,-1,-1,-1);
	if(!(j<0 || j>=sz(seq))) {
		ret = vstem[j];
	}
	return ret;
}

int StackFold::getsegment_esum(seglist & in_sgl, std::vector<ETYPE> in_sgt) {
	assert(sz(in_sgl)==sz(in_sgt));
	int ret = 0;
	FOR(k,0,sz(in_sgl)) {
		segment sg = in_sgl[k];
		int i = sg.lp; int j = sg.rp;
		ETYPE type = in_sgt[k];
		switch(type) {
			case  etype_F: {ret += F[j]; assert(i==0); break;}
			case  etype_C: {ret += C[i][j]; break;}
			case  etype_FM: {ret += FM[i][j]; break;}
		}
	}
	return ret;
}


//#
std::vector<RNA> StackFold::generate_subopt(std::ostream & fout, const int & e_ub, bool & reqLocOpt, TurnerE & TE) {
	std::vector<RNA> ret;
	int n = sz(seq);
	int m = sz(stl);
	
	std::vector<PartStructT> R;
	
	PartStructT phi;
	segment sg_n(0, n-1);
	phi.sgl.PB(sg_n);
	phi.sgt.PB(etype_F);
	phi.elphi = 0;
	//phi <- ({[0..n-1]}, {F}, {}, 0)
	

	R.PB(phi); //R.PB(phi);
	int indx = 0;
	int num = 0;
	while(sz(R)!=0) {
		bool happen = false;
		PartStructT phi_0 = R[0];
//		phi_0.show(std::cout);
		//(theta, thetype, P, e) <- phi_0
		seglist theta = phi_0.sgl;
		std::vector<ETYPE> thetype = phi_0.sgt;
		stemlist P = phi_0.stl;
		int e = phi_0.elphi;
		std::string log = phi_0.log;
		R.erase(R.begin());
		if(sz(theta)==0) {
			//print the secondary structure of phi
			RStack rst(P.vstems, n);
			int eprim = computeEnergy(TE,seq,rst);
			if(eprim <= e_ub) {
				num++;
				if(num%500==0) {
					std::cout<<"The algorithm has generated "<<num<<" local optimal stack configurations"<<std::endl;
				}
				ret.PB(RNA(rst.tostring(), eprim));
			}
			continue;
		}
		
		//sg_0 =  thetha.pop(); type_0 = thetype.pop();
		segment sg_0 = theta.pop();
		ETYPE   type_0 = thetype[0]; thetype.erase(thetype.begin());
		
		PartStructT pprim(theta, thetype, P, e, log);
		switch(type_0) {
			case etype_F: {happen = suboptF(R, sg_0, pprim, e_ub, reqLocOpt, TE); break;}
			case etype_C: {happen = suboptC(R, sg_0, pprim, e_ub, reqLocOpt, TE); break;}
			case etype_FM:{happen = suboptFM(R, sg_0, pprim, e_ub, reqLocOpt, TE);break;}
			default : {std::cerr<<"ETYPE="<<type_0<<" is unknown."<<std::endl; exit(0);}
		}
		
/*		if(!happen && type_0 != eetype_C) {
			if(type_0 == etype_C) {
				int sx = sg_0.lp; int sy = sg_0.rp;
				if(sy<sz(seq)-1) 
					pprim.elphi -= (d(TE,seq,sx,sy,sy+1)+TerminalPenalty(TE,seq[sx],seq[sy]));
			}
			int esum = getsegment_esum(pprim.sgl, pprim.sgt);
			int emin = esum + pprim.elphi - 170;
			if((!reqLocOpt || (reqLocOpt && pprim.isLocOpt(stl,C,H))) && emin < e_ub) {
				R.PB(pprim);
				std::cout<<"Push back 4"; pprim.show();
			}
		}*/
		indx++;
	}
	return ret;
}

bool StackFold::suboptF(std::vector<PartStructT> & R, 
     const segment & in_sg, const PartStructT & in_phi, const int & e_ub, 
     const bool & reqLocOpt, TurnerE & TE) {
	bool happen = false;
	int i = in_sg.lp;
	int j = in_sg.rp;
	assert(i==0);
	std::vector<int> vindx = stl.getstemindice_endat(j);
		
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl[indx_t];
		int tx = t.outer.lp; int ty = t.outer.rp;
		int tp = t.inner.lp; int tq = t.inner.rp;
		PartStructT pprim (in_phi);
		segment w = pprim.sgl.findRS(segment(tx,ty));
		int wx = w.lp; int wy = w.rp;
		assert(!w.valid() || (w.valid() && ty<wx));
		
		int delta = 0;
		if(w.valid()) {
			if(ty==wx-1) 
				delta += 0;
			else if(ty==wx-2)
				delta += min(d(TE,seq,tx,ty,ty+1), d(TE,seq,wx,wy,wx-1));
			else if(ty<=wx-3)
				delta += d(TE,seq,tx,ty,ty+1) + d(TE,seq,wx,wy,wx-1);
		} else {//!w.valid()
			if(j==sz(seq)-1) {
				delta += 0;
			} else if(j<=sz(seq)-2) {
				delta += d(TE,seq,tx,ty,ty+1);
			}
		}
		
		int delta2 = 0;
		if(tx==0)
			delta2 += 0;
		else if(tx>=1) 
			delta2 += d(TE,seq,tx,ty,tx-1);

		segment sg_t(tx,ty);
		pprim.sgl.insert(sg_t);
		pprim.sgt.insert(pprim.sgt.begin(), etype_C);
		pprim.elphi += delta + delta2 + TerminalPenalty(TE,seq[tx],seq[ty]);
		int esum = getsegment_esum(pprim.sgl, pprim.sgt);
		int emin = esum + pprim.elphi;
		
		if((!reqLocOpt || (reqLocOpt && pprim.isLocOpt(stl,C,H))) && emin < e_ub) {
			std::stringstream ss;
			ss<<"subopt F["<<i<<"]["<<j<<"], elphi="<<pprim.elphi<<"=elphi+dR("<<delta<<")+dL("<<delta2<<")+AU("<<TerminalPenalty(TE,seq[tx],seq[ty])<<")"<<std::endl;
			std::string newlog=ss.str();
			pprim.log += newlog;
			R.PB(pprim);
#ifdef DEBUG
			std::cout<<"Push back F.1, "; pprim.show();
#endif
			happen = true;
		}
		if(tx<2*cf.MINSTACKLEN+cf.MINLOOPLEN) 
			continue;
		
		PartStructT pprim2(in_phi);
		segment sg_f(0,tx-1);
		pprim2.sgl.insert(sg_t);
		pprim2.sgl.insert(sg_f);
		pprim2.sgt.insert(pprim2.sgt.begin(), etype_C);
		pprim2.sgt.insert(pprim2.sgt.begin(), etype_F);
		pprim2.elphi += delta + TerminalPenalty(TE,seq[tx],seq[ty]);
		esum = getsegment_esum(pprim2.sgl, pprim2.sgt);
		emin = esum + pprim2.elphi + (-170*2);
		if((!reqLocOpt || (reqLocOpt && pprim2.isLocOpt(stl,C,H))) && emin < e_ub) 
		{
			std::stringstream ss;
			ss<<"subopt F["<<i<<"]["<<j<<"], elphi="<<pprim2.elphi<<"=elphi+dR("<<delta<<")+AU("<<TerminalPenalty(TE,seq[tx],seq[ty])<<")"<<std::endl;
			std::string newlog=ss.str();
			pprim2.log += newlog;
			
			R.PB(pprim2);
#ifdef DEBUG
			std::cout<<"Push back F.2, "; pprim2.show();
#endif
			happen = true;
		}
	}
	if(j>=2*cf.MINSTACKLEN+cf.MINLOOPLEN) {
		segment sg_f(0,j-1);
		PartStructT pprim3 (in_phi);
		pprim3.sgl.insert(sg_f);
		pprim3.sgt.insert(pprim3.sgt.begin(), etype_F);
		int esum = getsegment_esum(pprim3.sgl, pprim3.sgt);
		int emin = esum + pprim3.elphi + (-170*2);
		if((!reqLocOpt || (reqLocOpt && pprim3.isLocOpt(stl,C,H))) && emin < e_ub) 
		{
			std::stringstream ss;
			ss<<"subopt F["<<i<<"]["<<j<<"], elphi="<<pprim3.elphi<<std::endl;
			std::string newlog=ss.str();
			pprim3.log += newlog;
			R.PB(pprim3);
#ifdef DEBUG
			std::cout<<"Push back F.3, "; pprim3.show();
#endif
			happen = true;
		}
	}
	return happen;
}
bool StackFold::suboptFM(std::vector<PartStructT> & R, 
     const segment & in_sg, const PartStructT & in_phi, const int & e_ub, 
     const bool & reqLocOpt, TurnerE & TE) {
	bool happen = false;
	int i = in_sg.lp;
	int j = in_sg.rp;
	
	std::vector<int> vindx = stl.getstemindice_within(i,j);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl[indx_t];
		int tx = t.outer.lp; int ty = t.outer.rp;
		int tp = t.inner.lp; int tq = t.inner.rp;
		
		PartStructT pprim (in_phi);
		stem v = pprim.stl.findLS(t);
		segment w = pprim.sgl.findRS(segment(tx,ty));

#ifdef DEBUG		
		std::cout<<"pprim.stl is"<<std::endl;
		pprim.stl.show();
		std::cout<<"pprim.sgl is"<<std::endl;
		pprim.sgl.show();

		std::cout<<"t= ";t.show();std::cout<<std::endl;
		std::cout<<"v=stl.findLS(t)= ";v.show();std::cout<<std::endl;
		std::cout<<"w=sgl.findRS(t)= ";w.show();std::cout<<std::endl;
#endif		
		int vx = v.outer.lp; int vy = v.outer.rp;
		int vp = v.inner.lp; int vq = v.inner.rp;
		int wx = w.lp; int wy = w.rp;
		assert(v.valid() && w.valid() && v.enclose(t) && vp<=wx && vq>=wy && ty<wx);

		int delta = 0;
		if(ty==wx-1) 
			delta += 0;
		else if(ty==wx-2)
			delta += min(d(TE,seq,tx,ty,ty+1), d(TE,seq,wx,wy,wx-1));
		else if(ty<=wx-3)
			delta += d(TE,seq,tx,ty,ty+1) + d(TE,seq,wx,wy,wx-1);
		int delta2 = 0;
		if(tx==vp+1)
			delta2 += 0;
		else if(tx==vp+2)
			delta2 += min(d(TE,seq,vp,vq,vp+1), d(TE,seq,tx,ty,tx-1));
		else if(tx>=vp+3)
			delta2 += d(TE,seq,vp,vq,vp+1) + d(TE,seq,tx,ty,tx-1);
				
		segment sg_t(tx,ty);
		pprim.sgl.insert(sg_t);
		pprim.sgt.insert(pprim.sgt.begin(), etype_C);
		pprim.elphi += delta + TerminalPenalty(TE,seq[tx],seq[ty]) + delta2 + 
		               TE.Mphi;
		int esum = getsegment_esum(pprim.sgl, pprim.sgt);
		int emin = esum + pprim.elphi;
		
		if((!reqLocOpt || (reqLocOpt && pprim.isLocOpt(stl,C,H))) && emin < e_ub) {
			std::stringstream ss;
			ss<<"subopt FM["<<i<<"]["<<j<<"], elphi="<<pprim.elphi<<"=elphi+dR("<<delta<<")+dL("<<delta2<<")+AU("<<TerminalPenalty(TE,seq[tx],seq[ty])<<")+Mphi(40)"<<std::endl;
			std::string newlog=ss.str();
			pprim.log += newlog;
		
			R.PB(pprim);
#ifdef DEBUG
			std::cout<<"Push back FM.1, "; pprim.show();
#endif
			happen = true;
		}
		
		if(tx-vp<2*cf.MINSTACKLEN+cf.MINLOOPLEN+1)	{
			//it is impossible to find a stem which is enclosed by v and to the 
			//left of t
			continue;
		}

		PartStructT pprim2 (in_phi);
		segment sg_fm(vp+1,tx-1);
		pprim2.sgl.insert(sg_t);
		pprim2.sgl.insert(sg_fm);
		pprim2.sgt.insert(pprim2.sgt.begin(), etype_C);
		pprim2.sgt.insert(pprim2.sgt.begin(), etype_FM);
		pprim2.elphi += delta + TerminalPenalty(TE,seq[tx],seq[ty]) + TE.Mphi;

		esum = getsegment_esum(pprim2.sgl, pprim2.sgt);
		emin = esum + pprim2.elphi + (-170)*4;
		if((!reqLocOpt || (reqLocOpt && pprim2.isLocOpt(stl,C,H))) && emin < e_ub) 
		{
			std::stringstream ss;
			ss<<"subopt FM["<<i<<"]["<<j<<"], elphi="<<pprim2.elphi<<"=elphi+dR("<<delta<<")+AU("<<TerminalPenalty(TE,seq[tx],seq[ty])<<")"<<std::endl;
			std::string newlog=ss.str();
			pprim2.log += newlog;

		
			R.PB(pprim2);
#ifdef DEBUG
			std::cout<<"Push back FM.2, "; pprim2.show();
#endif
			happen = true;
		}
	}
	return happen;
}

bool StackFold::suboptC(std::vector<PartStructT> & R, 
     const segment & in_sg, const PartStructT & in_phi, const int & e_ub, 
     const bool & reqLocOpt, TurnerE & TE) {
	bool happen = false;
	int i = in_sg.lp;
	int j = in_sg.rp;
	
	std::vector<int> vindx = stl.getstemindice_startfrom_endat(i,j);
	assert(sz(vindx)==1); 
	int indx_s = vindx[0];
	stem s = stl[indx_s];

	PartStructT pprim (in_phi);
	pprim.stl.add(s);
	pprim.elphi += H[i][j];
	int esum = getsegment_esum(pprim.sgl, pprim.sgt);

	int emin = esum + pprim.elphi;
	if((!reqLocOpt || (reqLocOpt && pprim.isLocOpt(stl,C,H))) && emin < e_ub) {
		std::stringstream ss;
		ss<<"subopt C["<<i<<"]["<<j<<"], elphi="<<pprim.elphi<<"=H["<<i<<"]["<<j<<"]("<<H[i][j]<<")"<<std::endl;
		std::string newlog=ss.str();
		pprim.log += newlog;

		R.PB(pprim);
#ifdef DEBUG
		std::cout<<"Push back C.1, "; pprim.show();
#endif
		happen = true;
	}
	
	vindx = stl.getstemindice_enclosedby_or_inward_extend(indx_s);
	FOR(k,0,sz(vindx)) {
		int indx_t = vindx[k];
		stem t = stl[indx_t];
		stem s2 = getUsedPartOfStem(s, t);
		int tx = t.outer.lp; int ty = t.outer.rp;
		int tp = t.inner.lp; int tq = t.inner.rp;
		int sx = s2.outer.lp; int sy = s2.outer.rp;
		int sp = s2.inner.lp; int sq = s2.inner.rp;
		
		//interior loop
		PartStructT pprim (in_phi);
		segment sg_t(tx,ty);
		pprim.sgl.insert(sg_t);
		pprim.sgt.insert(pprim.sgt.begin(), etype_C);
		pprim.stl.add(s2);
		pprim.elphi += S[indx_s][indx_t]+phi[indx_s][indx_t];
		
		esum = getsegment_esum(pprim.sgl, pprim.sgt);
		emin = esum + pprim.elphi;
		if((!reqLocOpt || (reqLocOpt && pprim.isLocOpt(stl,C,H))) && emin < e_ub) 
		{
			std::stringstream ss;
			ss<<"subopt C["<<i<<"]["<<j<<"], elphi="<<pprim.elphi<<"=S["<<i<<"]["<<j<<"]("<<S[indx_s][indx_t]<<")+phi("<<phi[indx_s][indx_t]<<")"<<std::endl;
			std::string newlog=ss.str();
			pprim.log += newlog;
		
			R.PB(pprim);
#ifdef DEBUG
			std::cout<<"Push back C.2, "; pprim.show();
#endif
			happen = true;
		}
		
		if(tx-1-(sp+1)+1 < 2*cf.MINSTACKLEN+cf.MINLOOPLEN) 
			continue;
		//multi-loop
		int delta = 0;
		if(ty==sq-1) {
			delta += 0;
		} else if(ty==sq-2) {
			delta += min(d(TE,seq,tx,ty,ty+1), d(TE,seq,sp,sq,sq-1));
		} else {
			delta += d(TE,seq,tx,ty,ty+1) + d(TE,seq,sp,sq,sq-1);
		}
		segment sg_fm(sp+1,tx-1);
		PartStructT pprim2(in_phi);
		pprim2.sgl.insert(sg_t);
		pprim2.sgl.insert(sg_fm);
		pprim2.sgt.insert(pprim2.sgt.begin(), etype_C);
		pprim2.sgt.insert(pprim2.sgt.begin(), etype_FM);
		pprim2.stl.add(s2);
		pprim2.elphi += S[indx_s][indx_t]+TE.Mc+2*TE.Mphi+
		                TerminalPenalty(TE,seq[sp],seq[sq])+delta+
		                TerminalPenalty(TE,seq[tx],seq[ty]);
      esum = getsegment_esum(pprim2.sgl, pprim2.sgt);
		emin = esum + pprim2.elphi + (-170*4);
		
		if((!reqLocOpt || (reqLocOpt && pprim2.isLocOpt(stl,C,H))) && emin < e_ub)
		{
			std::stringstream ss;
			ss<<"subopt C["<<i<<"]["<<j<<"], s=("<<sx<<","<<sp<<")~("<<sq<<","<<sy<<"), t=("<<tx<<","<<tp<<")~("<<tq<<","<<ty<<")"<<std::endl<<"elphi="<<pprim2.elphi<<"=S["<<i<<"]["<<j<<"]("<<S[indx_s][indx_t]<<")+Mc+2*Mphi+"<<"dR("<<delta<<")+AUspsq("<<TerminalPenalty(TE,seq[sp],seq[sq])<<")+AUtxty("<<TerminalPenalty(TE,seq[tx],seq[ty])<<")"<<std::endl;
			std::string newlog=ss.str();
			pprim2.log += newlog;
		
			R.PB(pprim2);
#ifdef DEBUG			
			std::cout<<"Push back C.3, "; pprim2.show();
#endif
			happen = true;
		} else {
#ifdef DEBUG
			std::cout<<"No push back C.3, "; pprim2.show();
#endif
		}
	}
	
	return happen;
}

































