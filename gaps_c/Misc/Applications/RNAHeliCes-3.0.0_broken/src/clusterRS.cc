#include "clusterRS.hh"
bool operator == (const stemstem & stst1, const stemstem & stst2) {
	if(stst1.sti == stst2.sti && stst1.stj == stst2.stj) 
		return true; 
	else return false;
}

bool operator != (const stemstem & stst1, const stemstem & stst2){
	if(stst1.sti != stst2.sti || stst1.stj != stst2.stj) 
		return true; 
	else return false;
}

int StemCompareHelp(const stem & A, const stem & B)
{
	if(A.outer.lp != B.outer.lp) return (A.outer.lp - B.outer.lp);
	if(A.outer.rp != B.outer.rp) return (A.outer.rp - B.outer.rp);
	if(A.inner.lp != B.inner.lp) return (A.inner.lp - B.inner.lp);
	if(A.inner.rp != B.inner.rp) return (A.inner.rp - B.inner.rp);
	return 0;
}

bool operator < (const stemstem & stst1, const stemstem & stst2) {
	stem sti1 = stst1.sti; stem stj1 = stst1.stj;
	stem sti2 = stst2.sti; stem stj2 = stst2.stj;
	
	int c = StemCompareHelp(sti1, sti2);
	if(c != 0) return c < 0;
	
	return (StemCompareHelp(stj1, stj2) < 0);
}

bool operator > (const stemstem & stst1, const stemstem & stst2) {
	return stst2 < stst1;
}

std::vector<int> RStack2vectorint(const RStack & rs) {	
	RStack rs1(rs);
	std::vector<int> ret;
	FOR(i,0,rs.len) ret.PB(-1);
	FOR(i,0,sz(rs.vstems)) {
		stem s = rs.vstems[i];
		int x = s.outer.lp; int y = s.outer.rp; 
		int p = s.inner.lp; int q = s.inner.rp;
		FOR(k,x,p+1) {
			assert(ret[k] == -1 && ret[x+y-k] == -1);
			ret[k] = x+y-k;
			ret[x+y-k] = k;
		}
	}
	return ret;
}

std::vector<int> str2vectorint(const std::string & str) {
	return RStack2vectorint(RStack(str));
}

//##-2_alternative
std::vector<double> compute_conflictingdegree(const std::string & seq, const std::string & s1, 
               const std::string & s2, std::map<stemstem, std::vector<int> > & CDSTmap) {
	RStack rs1(s1);
	RStack rs2(s2);
	return compute_conflictingdegree(seq, rs1, rs2, CDSTmap);
}

//##-4
std::vector<int> compute_conflictingdegree_between_stems(stem & sti, stem & stj) {
	//SLM_DEBUG std::cout<<"the conflicting degree between "; sti.show(); std::cout<<" "; stj.show(); std::cout<<std::endl;
	int ix = sti.outer.lp; int iy = sti.outer.rp;
	int ip = sti.inner.lp; int iq = sti.inner.rp;

	int jx = stj.outer.lp; int jy = stj.outer.rp;
	int jp = stj.inner.lp; int jq = stj.inner.rp;
	std::vector<int> ret; 
	//if sti and stj are identical 
	if(sti==stj) {
		FOR(i,0,ip-ix+1) {
			ret.PB(0);
		}
//		FOR(i,0,ip-ix+1) 
//			std::cout<<ret[i]<<"\t";
//		std::cout<<std::endl;
		return ret;
	}

	//case -1 : not set;	
	//case 0  : set to 0 for sure;
	//case 1  : set to 1 for sure;
	//case 2  : set to samedirect or greater;
	FOR(i,0,ip-ix+1) {
	//initiallization, -1 means not set 
		ret.PB(-1);
	}
		
	if(jy<ix || iy<jx || !sti.conflictwith(stj)) {
		//since sti and stj do not conflict with each other, return an array ret of size 0.
//		FOR(i,0,ip-ix+1) 
//			std::cout<<ret[i]<<"\t";
//		std::cout<<std::endl;
		return ret;
	}
		
	FOR(i,ix,ip+1) {          // 0123456789012
		int k = ix+iy-i;  // ((((....))))  i=2, ix=0, iy=12 => k=10
		int val;
		int segi = stj.inwhichseg(i);
		int segk = stj.inwhichseg(k);
		     // .....(((....)))................
		     // ...................((....))....
		if( (segi==1 && segk==1) ||        // compatible
		     // .....(((....))).......
		     // ..((............))....
		    (segi==3 && segk==3) ||
		     // ...................((....))....  
		     // .....(((....)))................
		    (segi==5 && segk==5) ||
		    // ..((............))....
		    // .....(((....))).......
		    (segi==1 && segk==5) ) {
		    val = -1; //not set;
		} else if(                          // conflict
		     // .....(((....))).........
		     // ............((....))....
		    (segi==1 && segk==2) ||
		     // .....(((....))).........
		     // .((...))................
		    (segi==4 && segk==5) ||
		    (segi==2 && segk==2) ||
		    (segi==4 && segk==4) ||
		    (segi==1 && segk==3) ||
		    (segi==3 && segk==5) ) {
		    val = 1; //set to 1
		} else if(                          // partially-conflict
		    // .....(((....))).........
		    // ......((........))......
		    (segi==2 && segk==3) ||
		    // .....(((....))).........
		    // ...((........)).........
		    (segi==3 && segk==4) ||
		    // .....(((..........)))........
		    // ........((........)).........    
		    (segi==1 && segk==4) ||
		    // .....(((..........)))...
		    // ......((........))......
		    (segi==2 && segk==5) ) {
		    val = 2; //set to samedirect or greater.
		} else if(                          

		    (segi==2 && segk==4) ) {  
		      // .....(((..........)))........   
		      // ......((..........)).........
		      if(ix+iy == jx+jy) {    // consistent
		    	val = 0;
		      // .....(((..........)))........   
		      // .......((..........)).........
		      } else {
			val = samedirect;
		      }
		} else assert(false);
		ret[i-ix] = val;
	}
	//SLM_DEBUG FOR(i,0,ip-ix+1) 
	//SLM_DEBUG 	std::cout<<ret[i]<<"\t";
	//SLM_DEBUG std::cout<<std::endl;
	return ret;
}


//##-3: the barrier energy between 2 stems, std::vector<int> stores the bp-relationship between two stems
std::vector<int> get_conflictingdegree_between_stems(std::map<stemstem, std::vector<int> > & CDSTmap, stem & sti, stem & stj) {
//	std::cout<<"get conflictingdegree between sti and stj\nsti = "; sti.show(); std::cout<<std::endl<<"stj = "; stj.show(); std::cout<<std::endl;
	int length = sti.inner.lp-sti.outer.lp+1;
	stemstem stst(sti, stj);
	std::map<stemstem, std::vector<int> >::iterator it;
	it = CDSTmap.find(stst);
	if(it!=CDSTmap.end()) {
		return (*it).second;
	} else {
		std::vector<int> value = compute_conflictingdegree_between_stems(sti, stj);
		assert(sz(value) == length);
		CDSTmap[stst] = value;
		return value;
	}
}

//##-2: calculate potential barrier energy of stems and store them in a std::vector
std::vector<double> compute_conflictingdegree(const std::string & seq, RStack & rs1, RStack & rs2, std::map<stemstem, std::vector<int> > & CDSTmap) {
//std::vector<double> compute_conflictingdegree(const std::string & seq, RStack & rs1, RStack & rs2, std::map<stemstem, std::vector<int> > & CDSTmap) {
	std::vector<double> ret;
	std::vector<int> vi1 = RStack2vectorint(rs1);
	std::vector<int> vi2 = RStack2vectorint(rs2);
	
	FOR(i,0,sz(rs1.vstems)) {
		stem sti = rs1.vstems[i];
		int ix = sti.outer.lp; int iy = sti.outer.rp;
		int ip = sti.inner.lp; int iq = sti.inner.rp;
		int leni = ip-ix+1;
		
		if(rs2.containST(sti)) {
			ret.PB(0);
			continue;
		}
		std::vector<int> vcdi;
		FOR(k,0,leni) {
/* vcdi[.]= -1: value is not set;
   vcdi[.]= 0 : set to 0 for sure;
   vcdi[.]= 1 : set to 1 for sure;
   vcdi[.]= 2 : set to samedirect or greater; */
			vcdi.PB(-1);
		}

		std::vector<int> vindx2 = rs2.findConflictingST(sti);
		FOR(indxj,0,sz(vindx2)) {
			int j = vindx2[indxj];
			stem stj = rs2.vstems[j];
			assert(sti.conflictwith(stj));
//			std::vector<int> vcdij = get_conflictingdegree_between_stems(CDSTmap, sti, stj);
			std::vector<int> vcdij = get_conflictingdegree_between_stems(CDSTmap, sti, stj);
			//compute_conflictingdegree_between_stems(sti, stj);
			if(sz(vcdi) != sz(vcdij)) {
				std::cout<<"\nsti is "<<std::endl;
				sti.show();
				std::cout<<"\nstj is "<<std::endl;
				stj.show();
				std::cout<<"size of vcdi is "<<sz(vcdi)<<std::endl;
				std::cout<<"size of vcdij is "<<sz(vcdij)<<std::endl;
				std::cout<<"length of sti is "<<(sti.inner.lp-sti.outer.lp+1)<<std::endl;
			}
			assert(sz(vcdi) == sz(vcdij));
			bool stop = true;
			FOR(k,0,sz(vcdij)) {
				if(vcdi[k]==-1) {
					vcdi[k] = vcdij[k];
					stop = false;
				} else if(vcdi[k]==0) {
					//set to 0 for sure
					assert(vcdij[k]==-1 || vcdij[k]==0);
				} else if(vcdi[k]==1) {
					assert(vcdij[k]!=0);
				} else if(vcdi[k]==2) {
					assert(vcdij[k]!=0);
					if(vcdij[k]==1) {
						vcdi[k] = 1;
					} else {
						stop = false;
					}
					//else if(vcdij[k]==-1) vcdi[k] = 2;
					//else if(vcdij[k]==2 ) vcdi[k] = 2;
				}
			}
			if(stop) break;
		}

		double sum = 0;
		FOR(k,0,sz(vcdi)) {
			if(vcdi[k]==-1) {
				sum += 1; //if never been set, the position counts up to 1
			} else if(vcdi[k]==0) {
				//if this position counts 0
				//sum += 0;
			} else if(vcdi[k]==1) {
				//if this position counts 1
				sum += 1;
			} else if(vcdi[k]==2) {
				//if this position counts samedirect or greater
				sum += samedirect;
			}
		}
		ret.PB(sum/(double)(leni));
	}
	return ret;
}


//distance = the number of different base pairs between a pair of structures
double compute_distance(const std::string & seq, const std::string & s1, const std::string & s2) {
	double ret = 0;

	std::cout<<"si="<<s1<<std::endl;
	std::cout<<"sj="<<s2<<std::endl;
	
	RStack rs1(s1);
	RStack rs2(s2);
	int crossnum = 0;
	int bpnum = 0;
	FOR(i,0,sz(rs1.vstems)) {
		stem sti = rs1.vstems[i];
		int ix = sti.outer.lp; int iy = sti.outer.rp;
		int ip = sti.inner.lp; int iq = sti.inner.rp;
		FOR(k, ix, iy+1) {
			int ilp = k; int irp = ix+iy-k;
			bpnum++;
			FOR(j,0,sz(rs2.vstems)) {
				stem stj = rs2.vstems[j];
				int jx = stj.outer.lp; int jy = stj.outer.rp;
				int jp = stj.inner.lp; int jq = stj.inner.rp;
				FOR(h, jx, jy+1) {
					int jlp = h; int jrp = jx+jy-h;
					if( (ilp<jlp && jlp<irp && irp<jrp) ||  
					    (jlp<ilp && ilp<jrp && jrp<irp) ) {
					    	crossnum++;
					}
				}
			}
		}
	}
	if(bpnum==0) {
		ret = sz(seq);
	} else {
		ret = (double)crossnum/(double)bpnum;
	}
	return ret;
}

//compute stacking energy of all stems in vst and return a std::map std::mapping stems to 
//their energies (which are integers);
std::map<stem, int> computeStemE(TurnerE & TE, const std::string & seq, std::vector<stem> & vst) {
	std::map<stem, int> ret;
	FOR(i,0,sz(vst)) {
		stem sti = vst[i];
		int e = computeS(TE, seq, sti); //computeS(TE, seq, rs1.vstems[h])
		ret[sti] = e;
	}
	return ret;
}

//obtain the stacking energy of stem st from std::map stemEmap, if it is not 
//available, call computeS and add it to stemEmap.
int getStemE(TurnerE & TE, const std::string & seq, stem & st, std::map<stem, int> & stemEmap) {
	std::map<stem,int>::iterator it;
	it = stemEmap.find(st);
	if(it!=stemEmap.end()) {
		return (*it).second;
	} else {
		int e = computeS(TE, seq, st);
		stemEmap[st] = e;
		return e;
	}
}

//##-1
double compute_pairwise_conflictingdegree(TurnerE & TE, const std::string & seq, 
      RStack & rs1, RStack & rs2, std::map<stem, int> & stemEmap, 
      std::map<stemstem, std::vector<int> > & CDSTmap) {
//computing the pairwise distance (approximated barrier) between rs1 and rs2.
//double compute_pairwise_conflictingdegree(TurnerE & TE, const std::string & seq, 
//      RStack & rs1, RStack & rs2, std::map<stem, int> & stemEmap, 
//      std::map<stemstem, std::vector<int> > & CDSTmap) {
  
        /////////////////////////////////////////////
        // compute the result saved in a std::vector vd //
	std::vector<double> vd = compute_conflictingdegree(seq, rs1, rs2, CDSTmap);
	///////////////////////////////
	// sum up all values in vd ////
	double sum = 0;
	FOR(h,0,sz(vd)) {
		int e = getStemE(TE, seq, rs1.vstems[h], stemEmap);
		sum += (vd[h] * (double)e)/(double)100;
	}
	return sum ;
}


/*
std::vector<int> compute_conflictingdegree_between_stems(stem & sti, stem & stj) {
	int ix = sti.outer.lp; int iy = sti.outer.rp;
	int ip = sti.inner.lp; int iq = sti.inner.rp;
	double leni = ip-ix+1;

	std::vector<int> ret; 
	FOR(i,0,ip-ix+1) {
		ret.PB(-1);
	}
	//case -1 : not set;
	//case 0  : set to 0 for sure;
	//case 1  : set to 1 for sure;
	//case 2  : set to samedirect for sure;
	//case 3  : if not set to samedirect, set to 1;
	//case 4  : if not set to 1, set to samedirect;
	//case 5  : 

	int jx = stj.outer.lp; int jy = stj.outer.rp;
	int jp = stj.inner.lp; int jq = stj.inner.rp;
	double lenj = jp-jx+1;
	if(jy<ix || iy<jx || !sti.conflictwith(stj) || sti==stj) {
		return ret;
	}
	
//structure : .....((((((......))))))).......
//segment id: 1     2     3      4     5
	int segjx = sti.inwhichseg(jx); int segjy = sti.inwhichseg(jy);
	int segjp = sti.inwhichseg(jp); int segjq = sti.inwhichseg(jq);
	if(segjy==1) assert(false);
	if(segjy==2) {//?,?,?,2
		if(segjq==1) {//?,?,1,2
			//1,1,1,2
			FOR(i,ix,jy+1) {
				ret[i-ix] = 1; //set to 1 for sure
			}
//			ret = jy-ix+1;
		} else if(segjq==2) {//?,?,2,2
			if(segjp==1) {//?,1,2,2,
				//1,1,2,2
//				ret = jy-ix+1;
				FOR(i,ix,jq) {
					ret[i-ix] = 3;//if not set to samedirect, set to 1;
				}
				FOR(i,jq,jy+1) {
					ret[i-ix] = 1;//set to 1 for sure
				}
			} else if(segjp==2) {//?,2,2,2
				if(segjx==1) {//1,2,2,2
//					ret = (jy-jp) + samedirect*(jp-ix+1);
					FOR(i,ix,jp+1) ret[i-ix] = 2; //set to samedirect for sure. 
					FOR(i,jp+1,jq) ret[i-ix] = 3; //if not set to samedirect, set to 1;
					FOR(i,jq,jy+1) ret[i-ix] = 1; //set to 1 for sure.
				} else if(segjx==2) {//2,2,2,2
//					ret = (jy-jp) + samedirect*lenj;
					FOR(i,jx,jp+1) ret[i-ix] = 2; //set to samedirect for sure.
					FOR(i,jp+1,jq) ret[i-ix] = 3; //if not set to samedirect, set to 1;
					FOR(i,jq,jy+1) ret[i-ix] = 1; //set to 1 for sure
				}
			}
		}
	} else if(segjy==3) {//?,?,?,3
		if(segjq==1) { //?,?,1,3
			//1,1,1,3
//			ret = leni;			
			FOR(i,ix,ip+1) ret[i-ix] = 1; //set to 1 for sure.
		} else if(segjq==2) {//?,?,2,3
			if(segjp==1) {//?,1,2,3
				//1,1,2,3
//				ret = leni;
				FOR(i,ix,jq)   ret[i-ix] = 3; //if not set to samedirect, set to 1;
				FOR(i,jq,ip+1) ret[i-ix] = 1; //set to 1 for sure.
			} else if(segjp==2) {//?,2,2,3
				if(segjx==1) {//1,2,2,3
//					ret = (ip-jp) + samedirect*(jp-ix+1);
					FOR(i,ix,jp+1) ret[i-ix] = 2;//set to samedirect for sure.
					FOR(i,jp+1,jq) ret[i-ix] = 3;//if not set to samedirect, set to 1;
					FOR(i,jq,iy+1) ret[i-ix] = 1;//set to 1 for sure.
				} else if(segjx==2) {//2,2,2,3
//					ret = (ip-jp) + samedirect*lenj;
					FOR(i,jx,jp+1) ret[i-ix] = 2;//set to samedirect for sure.
					FOR(i,jp+1,jq) ret[i-ix] = 3;//if not set to samedirect, set to 1;
					FOR(i,jq,iq+1) ret[i-ix] = 1;//set to 1 for sure.
				}
			}
		} else if(segjq==3) {//?,?,3,3
			if(segjp==1) {//1,1,3,3
//				ret = leni;
				FOR(i,ix,ip+1) ret[i-ix] = 1; //set to 1 for sure.
			} else if(segjp==2) {//?,2,3,3
				if(segjx==1) {//1,2,3,3
//					ret = (ip-jp) + samedirect*(jp-ix+1);
					FOR(i,ix,jp+1)   ret[i-ix] = 2;//set to samedirect for sure
					FOR(i,jp+1,ip+1) ret[i-ix] = 3;//if not set to samedirect, set to 1
				} else if(segjx==2) {//2,2,3,3
//					ret = (ip-jp) + samedirect*lenj;
//					FOR(i,ix,jx)     ret[i-ix] not set;
					FOR(i,jx,jp+1)   ret[i-ix] = 2;//set to samedirect for sure
					FOR(i,jp+1,ip+1) ret[i-ix] = 3;//if not set to samedirect, set to 1
				}
			} else if(segjp==3) {//?,3,3,3
				if(segjx==1) {//1,3,3,3
//					ret = samedirect*leni;
					FOR(i,ix,ip+1)   ret[i-ix] = 2;//set to samedirect for sure.
				} else if(segjx==2) {//2,3,3,3
//					ret = samedirect*(ip-jx+1);
					FOR(i,jx,ip+1)   ret[i-ix] = 2;//set to samedirect for sure.
				} else if(segjx==3) {//3,3,3,3
					assert(false);//!sti.conflictwith(stj);
				}
			}
		}
	} else if(segjy==4) {//?,?,?,4
		if(segjq==1) {//?,?,1,4
			//1,1,1,4
//			ret = leni;
			FOR(i,ix,ip+1) ret[i-ix] = 1;//set to 1 for sure.
		} else if(segjq==2) {//?,?,2,4
			if(segjp==1) {//?,1,2,4
				//1,1,2,4
//				ret = leni;
				FOR(i,ix,jq)   ret[i-ix] = 3;//if not set to samedirect, set to 1
				FOR(i,jq,iq+1) ret[i-ix] = 1;//set to 1 for sure.
			} else if(segjp==2) {//?,2,2,4
				if(segjx==1) {//1,2,2,4
//					ret = (ip-jp) + samedirect*(jp-ix+1);
					???//it seems there are three different cases. 	
				} else if(segjx==2) {//2,2,2,4
					ret = (ip-jp) + samedirect*lenj;
					???
				}
			}
		} else if(segjq==3) {//?,?,3,4
			if(segjp==1) {//?,1,3,4
				//1,1,3,4
//				ret = leni;
				FOR(i,ix,ix+iy-jy+1)   ret[i-ix] = 1; //set to 1 for sure
				FOR(i,ix+iy-jy+1,ip+1) ret[i-ix] = 2;//set to samedirect for sure
			} else if(segjp==2) {//?,2,3,4
				if(segjx==1) {//1,2,3,4
//					ret = (ip-jp) + samedirect*(jp-ix+1);
					int jystar = ix+iy-jy;
					int break1, break2;
					if(jp<jystar) {
						break1 = jp+1; break2 = jystar;
					} else {
						break1 = jystar; break2 = jp+1;
					}
					FOR(i,ix,break1)     ret[i-ix] = 4;//if not set to 1, set to samedirect
					FOR(i,break1,break2) ret[i-ix] = 1;//set to 1 for sure
					FOR(i,break2,ip+1)   ret[i-ix] = 4;//if not set to 1, set to samedirect
				} else if(segjx==2) {//2,2,3,4
//					ret = (ip-jp) + samedirect*lenj;
					
				}
			} else if(segjp==3) {//?,3,3,4
				if(segjx==1) {//1,3,3,4
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,3,4
					ret = max((ip-jx+1), (jy-iq+1));
				} else if(segjx==3) {//3,3,3,4
					ret = samedirect*(jy-iq+1);
				}
			}
		} else if(segjq==4) {//?,?,4,4
			if(segjp==1) {//?,1,4,4
				//1,1,4,4
				ret = samedirect*lenj;
			} else if(segjp==2) {//?,2,4,4
				if(segjx==1) {//1,2,4,4
					ret = samedirect*lenj;
				} else if(segjx==2) {//2,2,4,4
					if(ix+iy==jx+jy) ret = 0;
					else ret = samedirect*lenj;
				} 
			} else if(segjp==3) {//?,3,4,4
				if(segjx==1) {//1,3,4,4
					assert(false);
				} else if(segjx==2) {//2,3,4,4
					ret = (jq-iq) + samedirect*lenj;
				} else if(segjx==3) {//3,3,4,4
					ret = (jq-iq) + samedirect*lenj;
				}
			} else if(segjp==4) {//?,4,4,4
				if(segjx==1) {//1,4,4,4
					assert(false);
				} else if(segjx==2) {//2,4,4,4
					ret = (jq-iq) + samedirect*lenj;
				} else if(segjx==3) {//3,4,4,4
					ret = (jq-iq) + samedirect*lenj;
				} else if(segjx==4) {//4,4,4,4
					ret = (jq-jx) + samedirect*lenj;
				}
			}
		}
	} else if(segjy==5) {//?,?,?,5
		if(segjq==1) {//?,?,1,5
			//1,1,1,5
			ret = leni;
		} else if(segjq==2) {//?,?,2,5
			if(segjp==1) {//?,1,2,5
				//1,1,2,5
				ret = samedirect*leni;
			} else if(segjp==2) {//?,2,2,5
				if(segjx==1) {//1,2,2,5
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,2,5
					assert(false);
				}
			}
		} else if(segjq==3) {//?,?,3,5
			if(segjp==1) {//?,1,3,5
				//1,1,3,5
				ret = samedirect*leni;
			} else if(segjp==2) {//?,2,3,5
				if(segjx==1) {//1,2,3,5
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,3,5
					assert(false);
				}
			} else if(segjp==3) {//?,3,3,5
				if(segjx==1) {//1,3,3,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,3,5
					ret = samedirect*leni;
				} else if(segjx==3) {//3,3,3,5
					ret = samedirect*leni;
				}
			}
		} else if(segjq==4) {//?,?,4,5
			if(segjp==1) {//?,1,4,5
				//1,1,4,5
			} else if(segjp==2) {//?,2,4,5
				if(segjx==1) {//1,2,4,5
					ret = max((jp-ix+1), (iy-jq+1));
				} else if(segjx==2) {//2,2,4,5
					ret = (jx-ix) + samedirect*lenj;
				}
			} else if(segjp==3) {//?,3,4,5
				if(segjx==1) {//1,3,4,5
					ret = samedirect*lenj;
				} else if(segjx==2) {//2,3,4,5
					ret = min((jq-iq) + samedirect*(iy-jq+1), (jx-ix) + samedirect*(ip-jx+1));
				} else if(segjx==3) {//3,3,4,5
					ret = (jq-iq) + samedirect*(iy-jq+1);
				}
			} else if(segjp==4) {//?,4,4,5
				if(segjx==1) {//1,4,4,5
					ret = (jq-iq) + samedirect*(iy-jq+1);
				} else if(segjx==2) {//2,4,4,5
					ret = min((jq-iq) + samedirect*(iy-jq+1), (jx-ix) + samedirect*(ip-jx+1));
				} else if(segjx==3) {//3,4,4,5
					ret = (jq-iq) + samedirect*(iy-jq+1);
				} else if(segjx==4) {//4,4,4,5
					ret = iy-jx+1;
				}
			}
		} else if(segjq==5) {//?,?,5,5
			if(segjp==1) {//?,1,5,5
				//1,1,5,5
				assert(false); //!sti.conflictwith(stj);
			} else if(segjp==2) {//?,2,5,5
				if(segjx==1) {//1,2,5,5
					ret = samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,5,5
					ret = (jx-ix) + samedirect*lenj;
				}
			} else if(segjp==3) {//?,3,5,5
				if(segjx==1) {//1,3,5,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,5,5
					ret = (jx-ix) + samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,3,5,5
					ret = leni;
				}
			} else if(segjp==4) {//?,4,5,5
				if(segjx==1) {//1,4,5,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,4,5,5
					ret = (jx-ix) + samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,4,5,5
					ret = leni;
				} else if(segjx==4) {//4,4,5,5
					ret = iy-jx+1;
				}
			} else if(segjp==5) {//?,5,5,5
				if(segjx==1) {//1,5,5,5
					ret = leni;
				} else if(segjx==2) {//2,5,5,5
					ret = (jx-ix) + samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,5,5,5
					ret = leni;
				} else if(segjx==4) {//4,5,5,5
					ret = iy-jx+1;
				} else if(segjx==5) {//5,5,5,5
					assert(false); //!sti.conflictwith(stj);
				}
			}
		}
	}
	return ret/leni; 
} */
/*
double compute_conflictingdegree_between_stems(stem & sti, stem & stj)  {
	double ret = 0.0;
	std::vector<double> ret; 
	int ix = sti.outer.lp; int iy = sti.outer.rp;
	int ip = sti.inner.lp; int iq = sti.inner.rp;
	double leni = ip-ix+1;

	int jx = stj.outer.lp; int jy = stj.outer.rp;
	int jp = stj.inner.lp; int jq = stj.inner.rp;
	double lenj = jp-jx+1;
	if(jy<ix || iy<jx || !sti.conflictwith(stj) || sti==stj)
		return 0.0;
//structure : .....((((((......))))))).......
//segment id: 1     2     3      4     5
	int segjx = sti.inwhichseg(jx); int segjy = sti.inwhichseg(jy);
	int segjp = sti.inwhichseg(jp); int segjq = sti.inwhichseg(jq);
	//if (segjy==1) return 0;
	
	if(segjy==2) {//?,?,?,2
		if(segjq==1) {//?,?,1,2
			//1,1,1,2
			ret = jy-ix+1;	
		} else if(segjq==2) {//?,?,2,2
			if(segjp==1) {//?,1,2,2,
				//1,1,2,2
				ret = jy-ix+1;
			} else if(segjp==2) {//?,2,2,2
				if(segjx==1) {//1,2,2,2
					ret = (jy-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,2,2
					ret = (jy-jp) + samedirect*lenj;
				}
			}
		}
	} else if(segjy==3) {//?,?,?,3
		if(segjq==1) { //?,?,1,3
			//1,1,1,3
			ret = leni;
		} else if(segjq==2) {//?,?,2,3
			if(segjp==1) {//?,1,2,3
				//1,1,2,3
				ret = leni;
			} else if(segjp==2) {//?,2,2,3
				if(segjx==1) {//1,2,2,3
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,2,3
					ret = (ip-jp) + samedirect*lenj;
				}
			}
		} else if(segjq==3) {//?,?,3,3
			if(segjp==1) {//1,1,3,3
				ret = leni;
			} else if(segjp==2) {//?,2,3,3
				if(segjx==1) {//1,2,3,3
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,3,3
					ret = (ip-jp) + samedirect*lenj;
				}
			} else if(segjp==3) {//?,3,3,3
				if(segjx==1) {//1,3,3,3
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,3,3
					ret = samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,3,3,3
					assert(false);//!sti.conflictwith(stj);
				}
			}
		}
	} else if(segjy==4) {//?,?,?,4
		if(segjq==1) {//?,?,1,4
			//1,1,1,4
			ret = leni;
		} else if(segjq==2) {//?,?,2,4
			if(segjp==1) {//?,1,2,4
				//1,1,2,4
				ret = leni;
			} else if(segjp==2) {//?,2,2,4
				if(segjx==1) {//1,2,2,4
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,2,4
					ret = (ip-jp) + samedirect*lenj;
				}
			}
		} else if(segjq==3) {//?,?,3,4
			if(segjp==1) {//?,1,3,4
				//1,1,3,4
				ret = leni;
			} else if(segjp==2) {//?,2,3,4
				if(segjx==1) {//1,2,3,4
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,3,4
					ret = (ip-jp) + samedirect*lenj;
				}
			} else if(segjp==3) {//?,3,3,4
				if(segjx==1) {//1,3,3,4
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,3,4
					ret = max((ip-jx+1), (jy-iq+1));
				} else if(segjx==3) {//3,3,3,4
					ret = samedirect*(jy-iq+1);
				}
			}
		} else if(segjq==4) {//?,?,4,4
			if(segjp==1) {//?,1,4,4
				//1,1,4,4
				ret = samedirect*lenj;
			} else if(segjp==2) {//?,2,4,4
				if(segjx==1) {//1,2,4,4
					ret = samedirect*lenj;
				} else if(segjx==2) {//2,2,4,4
					if(ix+iy==jx+jy) ret = 0;
					else ret = samedirect*lenj;
				} 
			} else if(segjp==3) {//?,3,4,4
				if(segjx==1) {//1,3,4,4
					assert(false);
				} else if(segjx==2) {//2,3,4,4
					ret = (jq-iq) + samedirect*lenj;
				} else if(segjx==3) {//3,3,4,4
					ret = (jq-iq) + samedirect*lenj;
				}
			} else if(segjp==4) {//?,4,4,4
				if(segjx==1) {//1,4,4,4
					assert(false);
				} else if(segjx==2) {//2,4,4,4
					ret = (jq-iq) + samedirect*lenj;
				} else if(segjx==3) {//3,4,4,4
					ret = (jq-iq) + samedirect*lenj;
				} else if(segjx==4) {//4,4,4,4
					ret = (jq-jx) + samedirect*lenj;
				}
			}
		}
	} else if(segjy==5) {//?,?,?,5
		if(segjq==1) {//?,?,1,5
			//1,1,1,5
			ret = leni;
		} else if(segjq==2) {//?,?,2,5
			if(segjp==1) {//?,1,2,5
				//1,1,2,5
				ret = samedirect*leni;
			} else if(segjp==2) {//?,2,2,5
				if(segjx==1) {//1,2,2,5
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,2,5
					assert(false);
				}
			}
		} else if(segjq==3) {//?,?,3,5
			if(segjp==1) {//?,1,3,5
				//1,1,3,5
				ret = samedirect*leni;
			} else if(segjp==2) {//?,2,3,5
				if(segjx==1) {//1,2,3,5
					ret = (ip-jp) + samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,3,5
					assert(false);
				}
			} else if(segjp==3) {//?,3,3,5
				if(segjx==1) {//1,3,3,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,3,5
					ret = samedirect*leni;
				} else if(segjx==3) {//3,3,3,5
					ret = samedirect*leni;
				}
			}
		} else if(segjq==4) {//?,?,4,5
			if(segjp==1) {//?,1,4,5
				//1,1,4,5
			} else if(segjp==2) {//?,2,4,5
				if(segjx==1) {//1,2,4,5
					ret = max((jp-ix+1), (iy-jq+1));
				} else if(segjx==2) {//2,2,4,5
					ret = (jx-ix) + samedirect*lenj;
				}
			} else if(segjp==3) {//?,3,4,5
				if(segjx==1) {//1,3,4,5
					ret = samedirect*lenj;
				} else if(segjx==2) {//2,3,4,5
					ret = min((jq-iq) + samedirect*(iy-jq+1), (jx-ix) + samedirect*(ip-jx+1));
				} else if(segjx==3) {//3,3,4,5
					ret = (jq-iq) + samedirect*(iy-jq+1);
				}
			} else if(segjp==4) {//?,4,4,5
				if(segjx==1) {//1,4,4,5
					ret = (jq-iq) + samedirect*(iy-jq+1);
				} else if(segjx==2) {//2,4,4,5
					ret = min((jq-iq) + samedirect*(iy-jq+1), (jx-ix) + samedirect*(ip-jx+1));
				} else if(segjx==3) {//3,4,4,5
					ret = (jq-iq) + samedirect*(iy-jq+1);
				} else if(segjx==4) {//4,4,4,5
					ret = iy-jx+1;
				}
			}
		} else if(segjq==5) {//?,?,5,5
			if(segjp==1) {//?,1,5,5
				//1,1,5,5
				assert(false); //!sti.conflictwith(stj);
			} else if(segjp==2) {//?,2,5,5
				if(segjx==1) {//1,2,5,5
					ret = samedirect*(jp-ix+1);
				} else if(segjx==2) {//2,2,5,5
					ret = (jx-ix) + samedirect*lenj;
				}
			} else if(segjp==3) {//?,3,5,5
				if(segjx==1) {//1,3,5,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,3,5,5
					ret = (jx-ix) + samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,3,5,5
					ret = leni;
				}
			} else if(segjp==4) {//?,4,5,5
				if(segjx==1) {//1,4,5,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,4,5,5
					ret = (jx-ix) + samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,4,5,5
					ret = leni;
				} else if(segjx==4) {//4,4,5,5
					ret = iy-jx+1;
				}
			} else if(segjp==5) {//?,5,5,5
				if(segjx==1) {//1,5,5,5
					ret = samedirect*leni;
				} else if(segjx==2) {//2,5,5,5
					ret = (jx-ix) + samedirect*(ip-jx+1);
				} else if(segjx==3) {//3,5,5,5
					ret = leni;
				} else if(segjx==4) {//4,5,5,5
					ret = iy-jx+1;
				} else if(segjx==5) {//5,5,5,5
					assert(false); //!sti.conflictwith(stj);
				}
			}
		}
	}
	return ret/leni; 
} */
/*
std::vector<double> compute_conflictingdegree(const std::string & seq, RStack & rs1, RStack & rs2) {

	std::vector<double> ret;
	
	std::vector<int> vi1 = RStack2vectorint(rs1);
	std::vector<int> vi2 = RStack2vectorint(rs2);
#ifdef DEBUG
	std::cout<<"rs1: "; rs1.show(); std::cout<<std::endl;
	std::cout<<"rs2: "; rs2.show(); std::cout<<std::endl;
#endif
	
	FOR(i,0,sz(rs1.vstems)) {
		stem sti = rs1.vstems[i];

		int ix = sti.outer.lp; int iy = sti.outer.rp;
		int ip = sti.inner.lp; int iq = sti.inner.rp;
		int leni = ip-ix+1;

		double max_cdi = 0;
		FOR(j,0,sz(rs2.vstems)) {

			double max_cdij = 0;
			stem stj = rs2.vstems[j];
			int jx = stj.outer.lp; int jy = stj.outer.rp;
			int jp = stj.inner.lp; int jq = stj.inner.rp;
			int lenj = jp-jx+1;
#ifdef DEBUG
			std::cout<<"si.stem["<<i<<"]=";sti.show();std::cout<<std::endl;
			std::cout<<"sj.stem["<<j<<"]=";stj.show();std::cout<<std::endl;
#endif
			if(jx>iy) {
				break;
			}
			if(!sti.conflictwith(stj)){
				max_cdij = 0;
				continue;
			} else if(ix==jx && ip==jp && iq==jq && iy==jy) {
				max_cdij = 0;
				continue;
			} else {
				double cdleft = 0; 
				double cdright = 0;
				bool stop = false;
				double tmp = 0.0;
				if(ix>=jx && ip<=jy && iq>jy) max_cdij = 1;
				else if(ip<jx && iq>=jx && iy<=jy) max_cdij = 1;
				else {
					FOR(k,ix,ip+1) {
						//k pairs with k1 in stem sti, k pairs with k2 in stem stj,
						int k1 = vi1[k]; int k2 = vi2[k];
						assert(k1 > k);
						if(stop) break;
						if(k>=jx && k<=jp) {
							if(k1==k2) {
								//do nothing
							} else if(k2 <= ip) {
								tmp = (k2-k+1)/(double)leni;//cdleft
							} else if(k2>ip && k2<iq) {
								tmp = (ip-k+1)/(double)leni;//cdleft
							} else if(k2>=iq && k2<=iy) {
								tmp = samedirect*(abs(k1-k2)+1)/(double)leni;//cdleft, cdright
							} else if(vi2[k]>iy) {
								//note that k-ix+1 == iy-k1+1, so tmp == tmp
//								tmp = samedirect*(k-ix+1)/(double)leni;
								tmp = samedirect*(iy-k1+1)/(double)leni;//cdleft, cdright
							}
						} else if(k>jp && k<jq) {
							if(jq<=ip) {
								if(jy<=ip) {
									tmp = (iy-k+1)/(double)leni;
									cdleft = max(cdleft, tmp);
								} else if(jy>ip) {
									tmp = (ip-k+1)/(double)leni;
									cdleft = max(cdleft, tmp);
								}
							} else if(jq>ip && jy<iq) {
								if(ix<=jx) {
									tmp = samedirect*(jp-jx+1)/(double)leni + (ip-jp)/(double)leni;
								} else if(jx<ix) {
									if(jp<ix) {
										max_cdij = 1;
										stop = true;
									} else if(ix<=jp && jp<ip){
										tmp = (ip-jp)/(double)leni + samedirect*(jp-ix+1)/(double)leni;
									} else {assert(false);}
								}
							
								tmp = (ip-k+1)/(double)leni;
								cdleft = max(cdleft, tmp);
							}
						} else if(k>=jq && k<=jy) {
							if(jx>=ix) {
								cdleft = max(cdleft, (k-jx+1)/(double)leni);
							} else {
								cdleft = max(cdleft, (k-ix+1)/(double)leni);
							}
						} else {//k>jy or k<jx
							//do nothing
						}
						cdleft = max(cdleft, tmp);
						cdright = max(cdright, tmp);
					}
					FOR(k,ip+1,iq) {
//						assert(vi1[k] == -1);
						if(stop) break;
						//do nothing
						if(k>=jx && k<=jp) {
							if(vi2[k]<iq) {
								//do nothing
							} else if(vi2[k]>=iq && vi2[k]<=iy) {
								tmp = samedirect*(vi2[k]-iq+1)/(double)leni;
								cdright = (cdright<tmp)?(tmp):(cdright);
							} else if(vi2[k]>iy){
								tmp = 0.5;
								cdright = (cdright<tmp)?(tmp):(cdright);
							}
						} else if(k>jp && k<jq) {
							//do nothing
						} else if(k>=jq && k<=jy) {
							//do nothing
							if(vi2[k]<ix) {
								cdleft = 1; stop = true;
							} else if(vi2[k]>=ix && vi2[k]<=ip) {
								tmp = (ip-vi2[k]+1)/(double)leni;
								cdleft = (cdleft<tmp)?(tmp):(cdleft);
							} else if(vi2[k]>ip  && vi2[k]>iq) {
								//do nothing
							} //else impossible
						} else {//k<jx or k>jy
							//do nothing
						}
						cdleft = max(cdleft, tmp);
						cdright = max(cdright, tmp);
					}
					FOR(k,iq,iy+1) {
						assert(vi1[k] < k);						
						if(stop) break;
						if(k>=jx && k<=jp) {
							if(vi2[k]<=iy) {
								tmp = (vi2[k]-k+1)/(double)leni;
								cdright = (cdright<tmp)?(tmp):(cdright);
							} else if(vi2[k]>iy) {
								tmp = (iy-k+1)/(double)leni;
								cdright = (cdright<tmp)?(tmp):(cdright);
							}
						} else if(k>jp && k<jq) {
							if(jp+1<ix) {
								cdleft = 1; stop = true;
							} else if(jp+1>=ix && jp+1<=ip) {
								tmp = (ip-jp)/(double)leni;
								if(jx<ix){
									tmp += samedirect*(jp-ix+1)/(double)leni;
								} else if(jx>=ix) {
									tmp += samedirect*(jp-jx+1)/(double)leni;
								}
								cdleft = (cdleft<tmp)?(tmp):(cdleft);
							}
						} else if(k>=jq && k<=jy) {
							if(jq>ip && jq<iq) {
								tmp = samedirect*(k-iq+1)/(double)leni;
								cdright = (cdright<tmp)?(tmp):(cdright);
							} else if(jq>=iq && jq<=iy) {
								tmp = samedirect*(k-jq+1)/(double)leni;
								cdright = (cdright<tmp)?(tmp):(cdright);
							} else if(jq>=ix && jq<=ip) {
								tmp = (ip-jq+1)/(double)leni;
								cdleft = (cdleft<tmp)?(tmp):(cdleft);
								tmp = samedirect*(k-iq+1)/(double)leni;
								cdright = (cdright<tmp)?(tmp):(cdright);
							} else if(jq<ix) {
								cdleft = 1; stop = true;
							} else if(jq>iy) { //impossible
							}
						} else if(k<jx|| k>jy) {
							//do nothing
						}
						cdleft = max(cdleft, tmp);
						cdright = max(cdright, tmp);
					}
				}
				if(max_cdij <= cdleft) max_cdij = cdleft;
				if(max_cdij <= cdright) max_cdij = cdright;

#ifdef DEBUG
				std::cout<<"conflicting degree between si.stem["<<i<<"] and sj.stem["<<j<<"] is "<<max_cdij<<std::endl;
#endif
			}
			if(max_cdi <= max_cdij) max_cdi = max_cdij;
		}
#ifdef DEBUG
		std::cout<<std::endl<<"conflicting degree between si.stem["<<i<<"] and sj is "<<max_cdi<<std::endl<<std::endl;
#endif
		ret.PB((int(max_cdi *100))/(double)100);
	}

	return ret;
} */

/*

std::vector<std::vector<double> > compute_pairwiseconflictingdegree(const tring & seq, const std::vector<stem> & vstem) {
	//compute conflicting degree between pairwise stems
	std::vector<std::vector<double> > ret; 
	FOR(i,0,sz(vstem)) {
		std::vector<double> vd;
		ret.PB(vd);
		
		stem sti = vstem[i];
		int ix = sti.outer.lp; int iy = sti.outer.rp;
		int ip = sti.inner.lp; int iq = sti.inner.rp;
		int leni = ip-ix+1;
		
		FOR(j,0,sz(vst)) {
			if(i==j) ret[i].PB(0);
			
			stem stj = vstems[j];
			int jx = stj.outer.lp; int jy = stj.outer.rp;
			int jp = stj.inner.lp; int jq = stj.inner.rp;
			int lenj = jp-jx+1;
			if(!sti.conflictwith(stj))
				ret[i].PB(0);
			
			
			if(ix>=jx && ip<=jy && iq>jy) 
				ret[i].PB();
			else if(ip<jx && iq>=jx && iy<=jy) max_cdij = 1;
			
			FOR(k,ix,ip+1) {
				assert(vi1[k] > k);
				if(stop) break;
				if(k>=jx && k<=jp) {
					if(vi1[k] == vi2[k]) {
						//do nothing
					} else if(vi2[k] <= ip) {
						tmp = (vi2[k]-k+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
					} else if(vi2[k]>ip && vi2[k]<iq ) {
						tmp = (ip-k+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
					} else if(vi2[k]>= iq && vi2[k]<=iy) {
						tmp = samedirect*(abs(vi1[k]-vi2[k])+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
						cdright = (cdright<tmp)?(tmp):(cdright);
					} else if(vi2[k]>iy) {
						tmp = samedirect*(ip-k+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
						tmp = samedirect*(iy-vi1[k]+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					}
				} else if(k>jp && k<jq) {
					if(jq<=ip) {
						tmp = (ip-k+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
					} else if(jq>ip && jq<iq) {
						tmp = (ip-k+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
					} else if(jq>=iq && jq<=iy) {
						//do nothing
					} else if(jq>iy) {
						//do nothing
					}
				} else if(k>=jq && k<=jy) {
					tmp = (k-ix+1)/(double)leni;
					cdleft = (cdleft<tmp)?(tmp):(cdleft);
				} else {//k>jy or k<jx
					//do nothing
				}
			}
			
			FOR(k,ip+1,iq) {
				if(stop) break;
				//do nothing
				if(k>=jx && k<=jp) {
					if(vi2[k]<iq) {
						//do nothing
					} else if(vi2[k]>=iq && vi2[k]<=iy) {
						tmp = samedirect*(vi2[k]-iq+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					} else if(vi2[k]>iy){
						tmp = 0.5;
						cdright = (cdright<tmp)?(tmp):(cdright);
					}
				} else if(k>jp && k<jq) {
					//do nothing
				} else if(k>=jq && k<=jy) {
					//do nothing
					if(vi2[k]<ix) {
						cdleft = 1; stop = true;
					} else if(vi2[k]>=ix && vi2[k]<=ip) {
						tmp = (ip-vi2[k]+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
					} else if(vi2[k]>ip  && vi2[k]>iq) {
						//do nothing
					} //else impossible
				} else {//k<jx or k>jy
					//do nothing
				}
			}
			
			FOR(k,iq,iy+1) {
				assert(vi1[k] < k);						
				if(stop) break;
				if(k>=jx && k<=jp) {
					if(vi2[k]<=iy) {
						tmp = (vi2[k]-k+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					} else if(vi2[k]>iy) {
						tmp = (iy-k+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					}
				} else if(k>jp && k<jq) {
					if(jp+1<ix) {
						cdleft = 1; stop = true;
					} else if(jp+1>=ix && jp+1<=ip) {
						tmp = (ip-jp)/(double)leni;
						if(jx<ix){
							tmp += samedirect*(jp-ix+1)/(double)leni;
						} else if(jx>=ix) {
							tmp += samedirect*(jp-jx+1)/(double)leni;
						}
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
					}
				} else if(k>=jq && k<=jy) {
					if(jq>ip && jq<iq) {
						tmp = samedirect*(k-iq+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					} else if(jq>=iq && jq<=iy) {
						tmp = samedirect*(k-jq+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					} else if(jq>=ix && jq<=ip) {
						tmp = (ip-jq+1)/(double)leni;
						cdleft = (cdleft<tmp)?(tmp):(cdleft);
						tmp = samedirect*(k-iq+1)/(double)leni;
						cdright = (cdright<tmp)?(tmp):(cdright);
					} else if(jq<ix) {
						cdleft = 1; stop = true;
					} else if(jq>iy) { //impossible
					}
				} else if(k<jx|| k>jy) {
					//do nothing
				}
			}
		}
	}
	return ret; 
}
*/

/*
std::vector<double> compute_conflictingdegree(const std::string & seq, RStack & rs1, RStack & rs2) {
	std::vector<double> ret;
	std::vector<int> vi1 = RStack2vectorint(rs1);
	std::vector<int> vi2 = RStack2vectorint(rs2);
	std::cout<<std::endl;
	std::cout<<rs1.tostd::string()<<std::endl;
	std::cout<<rs2.tostd::string()<<std::endl;
	FOR(i,0,sz(rs1.vstems)) {
		stem sti = rs1.vstems[i];
		int ix = sti.outer.lp; int iy = sti.outer.rp;
		int ip = sti.inner.lp; int iq = sti.inner.rp;
		double leni = ip-ix+1;
		std::cout<<"stem["<<i<<"/"<<sz(rs1.vstems)<<"]=("<<ix<<","<<iy<<") ~ ("<<ip<<","<<iq<<")"<<std::endl;
	}

	std::cout<<std::endl;
	FOR(i,0,sz(rs1.vstems)) {
		double cdi = 0;
		bool stop = false;
		stem sti = rs1.vstems[i];
		int ix = sti.outer.lp; int iy = sti.outer.rp;
		int ip = sti.inner.lp; int iq = sti.inner.rp;
		double leni = ip-ix+1;
		std::cout<<"stem["<<i<<"/"<<sz(rs1.vstems)<<"]=("<<ix<<","<<iy<<") ~ ("<<ip<<","<<iq<<")"<<std::endl;
		
		if(rs2.containST(sti)) {
			cdi = 0;
		} else {
			std::vector<double> vcd;
			FOR(k,ix,ip+1) vcd.PB(0);
			assert(sz(vcd)==ip-ix+1);
	

			FOR(k,ix,ip+1) {
				int h  = vi1[k]; //in S1, k~h and k<h
				int k1 = vi2[k]; //in S2, k~k1
				int h1 = vi2[h]; //in S2, h~h1
				double tmp = 0;
				if(k1==h) {
					//do nothing. 
				} else if(k1==-1 && h1==-1) {
					tmp = 1;
				} else if(k1==-1 && h1!=-1) {
					if(h1<k) tmp = samedirect;
					else if(k<h1 && h1<h) tmp = samedirect;
					else if(h<h1) tmp = 1;
				} else if(k1!=-1 && h1==-1) {
					if(k1<k) tmp = 1;
					else if(k<k1 && k1<h) tmp = samedirect; //vcd[k-ix] = samedirect;
					else if(k1>h) tmp = samedirect; //vcd[k-ix] = samedirect;
				} else {//k1!=-1, k2!=-1
					double cdi_l, cdi_r;
					if(k1<k) cdi_l = 1;
					else cdi_l = samedirect;				
					if(h1<h) cdi_r = samedirect;
					else cdi_r = 1;
					tmp = max(cdi_l, cdi_r);
				}
				vcd[k-ix] = tmp;
			}

			std::cout<<"Phase 1 "<<std::endl;
			FOR(k,0,sz(vcd)) std::cout<<vcd[k]<<"\t";
			std::cout<<std::endl;


			FOR(k,ip+1,iq) {
				//in S1, k is unpaired.
				int k1 = vi2[k];//in S2, k~k1
				if(k1==-1) {
					//do nothing
				} else {
					int segk1 = sti.inwhichseg(k1);
					if(segk1==1) {
						int indxj = rs2.inwhichST(k1);
						if(indxj == -1) assert(false);
						stem stj = rs2.vstems[indxj];
						int jx = stj.outer.lp; int jy = stj.outer.rp;
						int jp = stj.inner.lp; int jq = stj.inner.rp;
						if(jp<ix) {
							cdi = 1;
							stop = true;
							break;
						} else {
							FOR(j,jp+1, ip+1) {
								vcd[j-ix] = 1;
							}
						}
					} else if(segk1==2) {
						//vcd[k1-ix] = samedirect; already set
						FOR(j,k1+1,ip+1) {
							if(eq(vcd[j-ix], 0)) {
								vcd[j-ix] = 1;
							}
						}
					} else if(segk1==3) {
						//do nothing
					} else if(segk1==4) {
						for(int j=k1-1; j>=iq; j--) {
							if(eq(vcd[iy-j],0)) {
								vcd[iy-j] = 1;
							}
						}
					} else if(segk1==5) {
						int indxj = rs2.inwhichST(k1);
						if(indxj == -1) assert(false);
						stem stj = rs2.vstems[indxj];
						int jx = stj.outer.lp; int jy = stj.outer.rp;
						int jp = stj.inner.lp; int jq = stj.inner.rp;
						if(iy<jq) {
							cdi = 1;
							stop = true;
							break;
						} else {
							for(int j=jq-1; j>=iq; j--) {
								vcd[iy-j] = 1;
							}
						}
					} else {
						assert(false);
					}
				}
			}

			std::cout<<"Phase 2 "<<std::endl;
			FOR(k,0,sz(vcd)) std::cout<<vcd[k]<<"\t";
			std::cout<<std::endl;

			if(!stop) {
				cdi = 0;
				FOR(k,0,sz(vcd)) {
					cdi += vcd[k];
				}
				cdi = cdi/leni;
			}
		}
		std::cout<<"cd["<<i<<"]="<<cdi<<std::endl;
		ret.PB(cdi);
	}
	return ret;
}
*/
