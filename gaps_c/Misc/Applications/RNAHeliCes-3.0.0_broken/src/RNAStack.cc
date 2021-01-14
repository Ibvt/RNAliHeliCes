#include "RNAStack.hh"
extern Config cf;

/***********************************************************
                         overload operator
***********************************************************/
bool operator == (const BP & in_bp1, const BP & in_bp2) {
	return ( (in_bp1.lp==in_bp2.lp && in_bp1.rp==in_bp2.rp) 
//	||       (in_bp1.lp==in_bp2.rp && in_bp1.rp==in_bp2.lp) 
	);
}
bool operator != (const BP & in_bp1, const BP & in_bp2) {
	return ( in_bp1.lp != in_bp2.lp ||  in_bp1.rp != in_bp2.rp);
}

bool operator == (const stem & st1, const stem & st2) {
	return (st1.outer.lp == st2.outer.lp && st1.outer.rp == st2.outer.rp
	     && st1.inner.lp == st2.inner.lp && st1.inner.rp == st2.inner.rp);
}
bool operator != (const stem & st1, const stem & st2) {
	return (st1.outer.lp != st2.outer.lp || st1.outer.rp != st2.outer.rp 
	     || st1.inner.lp != st2.inner.lp || st1.inner.rp != st2.inner.rp);
}

bool operator < (const stem & st1, const stem & st2) {
	return (st1.outer.lp < st2.outer.lp || 
	       (st1.outer.lp == st2.outer.lp && st1.outer.rp > st2.outer.rp));
}
bool operator > (const stem & st1, const stem & st2) {
	return (st1.outer.lp > st2.outer.lp ||
          (st1.outer.lp == st2.outer.lp && st1.outer.rp < st2.outer.rp));
}


/***********************************************************
                         Common
***********************************************************/
//random number generator
int myrand(int n)
{
	if(n==0) return 0;
   return rand() % n;
}

bool eq(double d1, double d2)
{
	return (d1 >= d2) ? ((d1-d2) <= PRECISION) : ((d2-d1) <= PRECISION);
}

//return the index of query if query is in vi
//otherwise return -1
//if direction=1, search in [startpos...sz(vi)-1], increasing order
//if direction=-1, search in from [startpos..0], decreasing order
int getindx(const std::vector<int> & vi, const int & query, const int & startpos, const int & direction){
	if(direction==1) {
		int i = (startpos>=0)?(startpos):0;
		for(; i<sz(vi); i++) {
			if(vi[i]==query) return i;
		}
		return -1;
	} else if(direction==-1) {
		if(startpos<0) return -1;
		int i = (startpos>=sz(vi)-1)?(sz(vi)-1):(startpos);
		for(; i>=0; i--) {
			if(vi[i]==query) return i;
		}
		return -1;
	} else {
		std::cerr<<"direction can either be 1 or -1 while it is "<<direction<<std::endl; exit(0);
	}
}

//return the index of in_bp in in_vbp
//if direction=1, search in [startpos...sz(in_vbp)-1], increasing order
//if direction=-1, search in from [startpos..0], decreasing order
int getindx(const std::vector<BP> & in_vbp, const BP & in_bp, const int & startpos, const int & direction) {
	if(direction==1) {
		int i = (startpos>=0)?(startpos):0;
		for(; i<sz(in_vbp); i++) {
			if(in_vbp[i].lp == in_bp.lp && in_vbp[i].rp==in_bp.rp) return i;
		}
		return -1;
	} else if(direction==-1) {
		if(startpos<0) return -1;
		int i = (startpos>=sz(in_vbp)-1)?(sz(in_vbp)-1):(startpos);
		for(; i>=0; i--) {
			if(in_vbp[i].lp == in_bp.lp && in_vbp[i].rp==in_bp.rp) return i;
		}
		return -1;
	} else {
		std::cerr<<"direction can either be 1 or -1 while it is "<<direction<<std::endl; exit(0);
	}
}
//return the index of in_a in in_va
//if direction=1, search in [startpos...sz(in_va)-1], increasing order
//if direction=-1, search in from [startpos..0], decreasing order
int getindx(const std::vector<Action> & in_va, const Action & in_a, const int & startpos, const int & direction){
	if(direction==1) {
		int i = (startpos>=0)?(startpos):0;
		for(;i<sz(in_va);i++){
			if(in_va[i].a == in_a.a && in_va[i].b.lp==in_a.b.lp && in_va[i].b.rp==in_a.b.rp) return i;
		}
		return -1;
	} else if(direction==-1){
		if(startpos<0) return -1;
		int i = (startpos>=sz(in_va)-1)?(sz(in_va)-1):(startpos);
		for(; i>=0; i--) {
			if(in_va[i].a == in_a.a && in_va[i].b.lp==in_a.b.lp && in_va[i].b.rp==in_a.b.rp) return i;
		}
		return -1;
	} else {
		std::cerr<<"direction can either be 1 or -1 while it is "<<direction<<std::endl; exit(0);
	}
}

int getindx(const std::vector<stem> & in_vst, const stem & in_st, const int & startpos, const int & direction) {
	stem myst = in_st;
	if(direction==1) {
		int i = (startpos>=0)?(startpos):0;
		for(;i<sz(in_vst);i++){
			if(myst == in_vst[i]) return i;
		}
		return -1;
	} else if(direction==-1){
		if(startpos<0) return -1;
		int i = (startpos>=sz(in_vst)-1)?(sz(in_vst)-1):(startpos);
		for(; i>=0; i--){
			if(myst == in_vst[i]) return i;
		}
		return -1;
	} else {
		std::cerr<<"direction can either be 1 or -1 while it is "<<direction<<std::endl; exit(0);
	}
}

//given a std::string representation of a strucuture, return vbp
std::vector<BP> s2v(const std::string & str)
{
#ifdef DEBUG
	std::cout<<"in s2v, str="<<str<<", sz(str)="<<sz(str)<<std::endl;
#endif
	std::stack<int> sta;
	std::vector<BP> ret;
	FOR(i,0,sz(str)) {
#ifdef DEBUG
//		std::cout<<"i="<<i<<", s[i]="<<s[i]<<std::endl;
#endif
		if (str[i] == '(') {
			sta.push(i);
		} else if(str[i] == ')') {
			assert(sz(sta)!=0);
			ret.PB(BP(sta.top(), i));
			sta.pop();
		}
#ifdef DEBUG
//		std::cout<<"the size of std::stack is "<<sz(sta)<<std::endl;
#endif
	}
	assert(sz(sta)==0);
	return ret;
}

//given indice of the left/right parenthesis and the length of the RStruct, return s
std::string v2s(const int & l, const std::vector<BP> & in_vbp)
{
#ifdef DEBUG
//	std::cout<<"in v2s, len="<<l<<std::endl;
#endif 

	std::string str = "";	
	FOR(i,0,l)
		str += " ";
	
	FOR(i,0,sz(in_vbp)) {
		str[in_vbp[i].lp] = '(';
		str[in_vbp[i].rp] = ')';
	}
	return str;
}


//return indices of unpaired bases
std::vector<int> getunp(const int & l, const std::vector<BP> & in_vbp) {
	std::vector<int> ret; 
	
	int tag_unp[l]; //tag_unp[i] represents if s[i] is an unpaired base
	FOR(i,0,l) {
		tag_unp[i] = 1; //assume every index is in unp
	}
	
	FOR(i,0,sz(in_vbp)) {
		tag_unp[in_vbp[i].lp] = 0;
		tag_unp[in_vbp[i].rp] = 0;
	}
	
	FOR(i,0,l) 
		if(tag_unp[i]==1)
			ret.PB(i);

	return ret;
}

//return indices of unpaired bases
std::vector<int> getunp(const std::string & str) {
	std::vector<int> ret; 
	FOR(i,0,sz(str)) 
		if(str[i]=='(' ||  str[i]==')'){}
		else	ret.PB(i);
	return ret;
}

//check if a file exists or not
bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile;
}

bool match(const char & x , const char & y) {
	if(x=='T' && y=='A') return true;
	if(x=='U' && y=='A') return true;
	if(x=='G' && y=='C') return true;
	if(x=='C' && y=='G') return true;
	if(x=='A' && y=='T') return true;
	if(x=='A' && y=='U') return true;
	if(cf.ALLOWGU==1) {
		if(x=='G' && y=='U') return true;
		if(x=='U' && y=='G') return true;
	}
	return false;
}

char getcomplementchar(const char & x) {
	switch(x) {
	case 'A': return 'U';
	case 'U':
	case 'T': return 'A';
	case 'G': return 'C';
	case 'C': return 'G';
	default : {std::cerr<<"The input nucleotide is "<<x<<std::endl; exit(0);}
	}
	return 'A';
}


bool descendingsort(double d1, double d2) {
	return (d1 > d2);
}
bool descendingsortBP(BP bp1, BP bp2) {
	return (bp1.lp > bp2.lp);
}
bool descendingsortST(stem st1, stem st2) {
	if (st1.outer.lp == st2.outer.lp) {
		return (st1.outer.rp > st2.outer.rp);
		//why not compare (st1.outer.rp == st2.outer.rp)? Because any two std::stacks \
		can not share exactly the same outer.lp and outer.rp
	} else
		return (st1.outer.lp > st2.outer.lp);
}
bool descendingsortSD(SD st1, SD st2) {
	return (st1 > st2);
}
bool descendingsortRNA(RNA rna1, RNA rna2) {
	if (rna1.e == rna2.e) 
		return (rna1.str > rna2.str);
	return (rna1.e > rna2.e);
}


bool ascendingsort(double d1, double d2) {
	return (d1 < d2);
}
bool ascendingsortBP(BP bp1, BP bp2) {
	return (bp1.lp < bp2. lp);
}
bool ascendingsortST(stem st1, stem st2) {
	if (st1.outer.lp == st2.outer.lp) {
		return (st1.outer.rp < st2.outer.rp);
		//why not compare (st1.outer.rp == st2.outer.rp)? Because any two std::stacks \
		can not share exactly the same outer.lp and outer.rp
	} else
		return (st1.outer.lp < st2.outer.lp);
}
bool ascendingsortSD(SD st1, SD st2) {
	return (st1 < st2);
}
bool ascendingsortRNA(RNA rna1, RNA rna2) {
	if (rna1.e == rna2.e) 
		return (rna1.str < rna2.str);
	return (rna1.e < rna2.e);
}

//check if in_bp clash with any base pair in in_vbp
bool mycheckclash(const BP & in_bp, const std::vector<BP>& in_vbp) {
	int x=in_bp.lp; int y=in_bp.rp;
	
	FOR(k,0,sz(in_vbp)) {
		int i = in_vbp[k].lp; int j = in_vbp[k].rp;
		//clash A: i ... x ... j ... y
		//clash B; x ... i ... y ... j
		if (i < x && x < j && j < y) return true;
		if (x < i && i < y && y < j) return true;
		if (x==i || x==j || y==i || y==j) return true;
	}
	return false;	
}

/***********************************************************
                         stem
***********************************************************/
bool stem::inward_extend(const stem & s) {
//(((((((((............)))))))))) --> stem s
//.....((((((.......))))))....... --> this stem
//or 
//(((((((((............)))))))))) --> stem s
//.........((((((....))))))...... --> this stem
//or
//(((((((((............)))))))))) --> stem s
//....((((((....))))))........... --> this stem
	int sx = s.outer.lp; int sy = s.outer.rp;
	int sp = s.inner.lp; int sq = s.inner.rp;
	int tx = outer.lp; int ty = outer.rp;
	int tp = inner.lp; int tq = inner.rp;
	bool ret = false;
	if( (sx<tx && tx<=sp && sp<=tp && tq<=sq && ty<sy && !enclosedby(s)) ||
	    (sx<tx && sp<=tp && tq<=sq && sq<=ty && ty<sy && !enclosedby(s)) )
	{
		double leftoverlap = -1; double rightoverlap = -1; 
		leftoverlap = std::max((double)(sp-tx+1)/(double)(sp-sx+1), 
		                  (double)(sp-tx+1)/(double)(tp-tx+1));
		rightoverlap = std::max((double)(ty-sq+1)/(double)(sp-sx+1),
		                   (double)(ty-sq+1)/(double)(tp-tx+1));
		if((leftoverlap < cf.MAXINWARDOVERLAP || eq(leftoverlap, cf.MAXINWARDOVERLAP)) && 
		   (rightoverlap <= cf.MAXINWARDOVERLAP || eq(rightoverlap, cf.MAXINWARDOVERLAP)) ) {
		   return true;
		} else {
			return false;
		}
	} else return false;
}

//##-5
int stem::inwhichseg(const int & query) {
	//stem st (x,p) ~ (q,y) is composed of five segments, including 
	//segment 1 : [0..x-1]
	//segment 2 : [x..p]
	//segment 3 : [p+1,q-1]
	//segment 4 : [q..y]
	//segment 5 : [y+1..+OO]
	//return in which segment is query
	int x = outer.lp; int y = outer.rp;
	int p = inner.lp; int q = inner.rp;
	
	if(query<x) {
		return 1; 
	} else if(x<=query && query<=p) {
		return 2;
	} else if(p<query  && query<q ) {
		return 3;
	} else if(q<=query && query<=y) {
		return 4;
	} else if(query>y) {
		return 5;
	} else {
		assert(false);
	}
}
/***********************************************************
                         stemlist
***********************************************************/
stemlist::stemlist(const std::string & seq, const std::set<stem> & vstem_set){
	s = seq;
	std::copy(vstem_set.begin(), vstem_set.end(), inserter(vstems, vstems.begin())); 
	sort(vstems.begin(), vstems.end(), ascendingsortST);
}
stemlist::stemlist(const std::string & seq){
//	cf.show();
	s = seq;
	int n = sz(seq);
	int itable[n][n];
	FOR(i,0,n) FOR(j,0,n) itable[i][j] = 0;
	FOR(i,0,n) 
		FOR(j,0,n) {
			if(match(seq[i], seq[n-1-j]) && (n-1-j)-i>cf.MINLOOPLEN) 
				itable[i][j] = 1;
		}
	FOR(i, 0, n)
		FOR(j, 0, n) {
			//start checking
			if(itable[i][j] == 1 && (i==0 || j==0 || (itable[i-1][j-1] == 0))) {
				int start_i = i; int start_j = j;
				while(true) {
					if(start_i+1>=n || start_j+1>=n || 
					   itable[start_i+1][start_j+1]==0) break;
					start_i++; start_j++;
				}
				stem newst(BP(i,  n-1-j), BP(start_i, n-1-start_j));
				if(start_i-i+1>=cf.MINSTACKLEN && newst.gethydrobonds(seq) >= cf.MINHYDROGEN) {
						vstems.PB(newst);
				} else {
					FOR(k, i, start_i+1) {
						itable[i][j] = 0;
					}
				}
			}
		}
	sort(vstems.begin(), vstems.end(), ascendingsortST);
}

void stemlist::show() {
	FOR(i,0,sz(vstems)){
		std::cout<<"stem ["<<i<<"/"<<sz(vstems)-1<<"] :";
		vstems[i].show();
		std::cout<<std::endl;
	}
	
	int n = sz(s);
	int itable[n][n];
	FOR(i,0,n) FOR(j,0,n) itable[i][j] = 0;
	FOR(i,0,n) 
		FOR(j,0,n) {
			if(match(s[i], s[n-1-j]) && (n-1-j)-i>cf.MINLOOPLEN) 
				itable[i][j] = 1;
		}
	FOR(i, 0, n)
		FOR(j, 0, n) {
			//start checking
			if(itable[i][j] == 1 && (i==0 || j==0 || (itable[i-1][j-1] == 0))) {
				int start_i = i; int start_j = j;
				while(true) {
					if(start_i+1>=n || start_j+1>=n || itable[start_i+1][start_j+1]==0) break;
					start_i++; start_j++;
				}
				if(!(start_i-i+1>=cf.MINSTACKLEN)) {
					FOR(k, i, start_i+1) {
						itable[i][j] = 0;
					}
				}
			}
		}
	std::cout<<"  ";
	for(int i = sz(s)-1; i>=0; i--) {
		std::cout<<s[i];
	}
	std::cout<<std::endl;
	FOR(i,0,sz(s)) {
		std::cout<<s[i]<<" ";
		FOR(j,0,n) if(itable[i][j]==1) std::cout<<"."; else std::cout<<" ";
		std::cout<<std::endl;
	}
}
//return std::stacks starting from s
std::vector<stem> stemlist::getstems_startfrom(int in_s) {
	std::vector<stem> ret;
	int i=0; int j=sz(vstems)-1;
	int find = -1;
	while(i<=j) {
		int mid = (i+j)/2;
		if(vstems[mid].outer.lp == in_s) {
			find = mid;
			break;
		} else if (vstems[mid].outer.lp > in_s) {
			j = mid-1;
		} else {
			i = mid+1;
		}
	}
	if(find == -1) return ret;
//	ret.PB(vstems[find]);
	int cursor = find-1;
	int start = find;
	int end = find;
	while(cursor>=0) {
		if(vstems[cursor].outer.lp == in_s) start=cursor;//ret.PB(vstems[cursor]);
		else break;
		cursor--;
	}
	
	cursor = find+1;
	while(cursor<sz(vstems)) {
		if(vstems[cursor].outer.lp ==in_s) end=cursor;
		else break;
		cursor++;
	}
	FOR(i,start, end+1) ret.PB(vstems[i]);
	
	FOR(i,0,sz(ret)) {
		assert(ret[i].outer.lp == in_s);
	}
	return ret;
}


//return std::stacks starting from s
std::vector<int> stemlist::getstemindice_startfrom(int in_s) {
	std::vector<int> ret;
	int i=0; int j=sz(vstems)-1;
	int find = -1;
	while(i<=j) {
		int mid = (i+j)/2;
		if(vstems[mid].outer.lp == in_s) {
			find = mid;
			break;
		} else if (vstems[mid].outer.lp > in_s) {
			j = mid-1;
		} else {
			i = mid+1;
		}
	}
	if(find == -1) return ret;
//	ret.PB(vstems[find]);
	int cursor = find-1;
	int start = find;
	int end = find;
	while(cursor>=0) {
		if(vstems[cursor].outer.lp == in_s) start=cursor;//ret.PB(vstems[cursor]);
		else break;
		cursor--;
	}
	
	cursor = find+1;
	while(cursor<sz(vstems)) {
		if(vstems[cursor].outer.lp == in_s) end=cursor;
		else break;
		cursor++;
	}
	FOR(i,start, end+1) ret.PB(i);
	
	FOR(i,0,sz(ret)) {
		assert(vstems[ret[i]].outer.lp == in_s);
	}
	return ret;
}

//return std::stacks ending at e
std::vector<stem> stemlist::getstems_endat(int in_e) {
	std::vector<stem> ret;
	if(in_e<0|| in_e>=sz(s)) return ret;
	FOR(i,0,sz(vstems)) {
		if(vstems[i].outer.rp == in_e) {
			ret.PB(vstems[i]);
		}
	}
	FOR(i,0,sz(ret)) {
		assert(ret[i].outer.rp==in_e);
	}
	return ret;
}

std::vector<int> stemlist::getstemindice_endat(int in_e) {
	std::vector<int> ret;
	if(in_e<0|| in_e>=sz(s)) return ret;
	FOR(i,0,sz(vstems)) {
		if(vstems[i].outer.rp == in_e) {
			ret.PB(i);
		}
	}
	FOR(i,0,sz(ret)) {
		assert(vstems[ret[i]].outer.rp==in_e);
	}
	return ret;
}

std::vector<stem> stemlist::getstems_startfrom_endat(int in_s, int in_e){
	std::vector<stem> ret;
	if(in_s<0 || in_s>=sz(s) || in_e<0|| in_e>=sz(s)) return ret;
	std::vector<stem> vst = getstems_startfrom(in_s);
	FOR(i,0,sz(vst)){
		if(vst[i].outer.rp == in_e) {
			ret.PB(vst[i]);
		}
	}
	FOR(i,0,sz(ret)) {
		assert(ret[i].outer.lp==in_s && ret[i].outer.rp==in_e);
	}
	return ret;
}

	
std::vector<int> stemlist::getstemindice_startfrom_endat(int in_s, int in_e){
	std::vector<int> ret;
	if(in_s<0 || in_s>=sz(s) || in_e<0|| in_e>=sz(s)) return ret;
	std::vector<int> vi = getstemindice_startfrom(in_s);
	FOR(i,0,sz(vi)){
		if(vstems[vi[i]].outer.rp == in_e) {
			ret.PB(vi[i]);
		}
	}
	FOR(i,0,sz(ret)) {
		assert(vstems[ret[i]].outer.lp == in_s && vstems[ret[i]].outer.rp == in_e);
	}
	return ret;
}


std::vector<stem> stemlist::getstems_within(int in_s, int in_e){
	std::vector<stem> ret;
	FOR(i,0,sz(vstems)){
		if(vstems[i].outer.lp >= in_s && vstems[i].outer.rp <= in_e){
			ret.PB(vstems[i]);
		}
	}
	FOR(i,0,sz(ret)) {
		assert(ret[i].outer.lp >= in_s && ret[i].outer.rp <= in_e);
	}
	return ret;
}

std::vector<int> stemlist::getstemindice_within(int in_s, int in_e){
	std::vector<int> ret;
	FOR(i,0,sz(vstems)){
		if(vstems[i].outer.lp >= in_s && vstems[i].outer.rp <= in_e){
			ret.PB(i);
		}
	}
	FOR(i,0,sz(ret)) {
		assert(vstems[ret[i]].outer.lp >= in_s && vstems[ret[i]].outer.rp <= in_e);
	}
	return ret;
}

std::vector<stem> stemlist::getstems_startfrom_within(int in_s, int in_e) {
	std::vector<stem> ret;
	std::vector<stem> vst;
	vst = getstems_startfrom(in_s); 
	FOR(i,0,sz(vst)) {
		if(vst[i].outer.rp <= in_e) ret.PB(vst[i]);
	}
	FOR(i,0,sz(ret)) {
		assert(ret[i].outer.lp == in_s && ret[i].outer.rp <= in_e);
	}
	return ret;
}

std::vector<int> stemlist::getstemindice_startfrom_within(int in_s, int in_e) {
	std::vector<int> ret;
	std::vector<int> vi = getstemindice_startfrom(in_s);
	FOR(i,0,sz(vi)) {
		if(vstems[vi[i]].outer.rp <= in_e) ret.PB(vi[i]);
	}
	FOR(i,0,sz(ret)) {
		assert(vstems[ret[i]].outer.lp == in_s && vstems[ret[i]].outer.rp <= in_e);
	}
	return ret;
}

std::vector<int> stemlist::getstemindice_endat_within(int in_s, int in_e) {
	std::vector<int> ret;
	std::vector<int> vi = getstemindice_endat(in_e);
	FOR(i,0,sz(vi)) {
		if(vstems[vi[i]].outer.lp >= in_s) ret.PB(vi[i]);
	}
	FOR(i,0,sz(ret)) {
		assert(vstems[ret[i]].outer.lp >= in_s && vstems[ret[i]].outer.rp == in_e);
	}
	return ret;
}

std::vector<int> stemlist::getstemindice_enclosedby(int indx_s) {
	std::vector<int> ret;
	stem s = vstems[indx_s];
	FOR(i,indx_s+1,sz(vstems)){
		stem sti = vstems[i];
		if(sti.outer.lp >= s.outer.rp) break;
		if(s.enclose(sti)) ret.PB(i);
	}
	return ret;
}

//	(((((((.......)))))))  -->stem s
//    ((((((....))))))    -->stem t

// (((((((.......)))))))  -->stem s
//    ((((.....))))       -->stem t

// (((((((.......)))))))  -->stem s
//       ((((...))))      -->stem t

// (((((((.........)))))))-->stem s
//       ((((....))))     -->stem t
std::vector<int> stemlist::getstemindice_inward_extend(int indx_s) {
	std::vector<int> ret;
	stem s = vstems[indx_s];
	FOR(i, indx_s+1, sz(vstems)) {
		stem t = vstems[i];
		if(t.inward_extend(s)){
			ret.PB(i);
		}
	}
	return ret;
}

std::vector<int> stemlist::getstemindice_enclosedby_or_inward_extend(int indx_s){
	std::vector<int> ret;
	stem s = vstems[indx_s];
	FOR(i, indx_s+1, sz(vstems)) {
		stem t = vstems[i];
		if(s.enclose(t) || t.inward_extend(s)){
			ret.PB(i);
		}
	}
	return ret;
}

bool stemlist::conflictwith(stem st) {
	FOR(i,0,sz(vstems)){
		if(vstems[i].conflictwith(st))
			return true;
	}
	return false;
}

bool stemlist::containST(stem st) {
	FOR(i,0,sz(vstems)){
		if(vstems[i] == st) return true;
	}
	return false;
}

void stemlist::insert(stem & in_st) {
	vstems.insert(vstems.begin(), in_st);
}
void stemlist::insert(stemlist & in_stl) {
	vstems.insert(vstems.begin(), in_stl.vstems.begin(), in_stl.vstems.end());
}
void stemlist::PB(stem & in_st) {
	vstems.PB(in_st);
}
void stemlist::PB(stemlist & in_stl) {
	vstems.insert(vstems.end(), in_stl.vstems.begin(), in_stl.vstems.end());
}

stem & stemlist::operator [] (const int & i) {
	return vstems[i];
}
void stemlist::add(stem & in_st) {
	if(sz(vstems)==0) {
		vstems.PB(in_st);
		return;
	}
	int i = 0; int j = sz(vstems)-1;
	
	while(i<=j) {
	  //std::cout << "i=" << i << "j=" << j << std::endl;
		int mid = (i+j)/2;
		stem m = vstems[mid];
		assert(m!=in_st);
		if(i==j) {
			if(m<in_st) {
				vstems.insert(vstems.begin()+i+1, in_st);
			} else {//m>in_st
				vstems.insert(vstems.begin()+i, in_st);
			}
			return ;
		}
		if(m<in_st) {
			i = mid+1;
		} else if(m>in_st) {
			j = mid-1;
		}
	}
	assert(i>j && i==j+1);
	vstems.insert(vstems.begin()+i, in_st);
}

stem stemlist::findLS(const stem & in_st) {
	int indx = -1;
	
	if(sz(vstems)!=0) {
		int i = 0; int j = sz(vstems)-1;
		while(i<=j) {
			int mid = (i+j)/2;
			stem m = vstems[mid];
			if(i==j) {
				if(m<in_st) {
					indx = mid;
				} else if(m>in_st) {
					indx = mid-1;
				} else if(m==in_st) {
					indx = mid-1;
				}
				break;
			}
			if(m<in_st) {
				i = mid+1;
			} else if(m>in_st) {
				j = mid-1;
			} else if(m==in_st) {
				indx = mid-1;
				break;
			}
		}
		if(i>j) indx = j;
	}
	
	if(indx<0 || indx>=sz(vstems))
		return stem(-1,-1,-1,-1);
	else 
		return vstems[indx];
}

stem stemlist::findRS(const stem & in_st) {
	int indx = -1;
	
	if(sz(vstems)!=0) {
		int i = 0; int j = sz(vstems)-1;
		while(i<=j) {
			int mid = (i+j)/2;
			stem m = vstems[mid];
			if(i==j) {
				if(m<in_st) {
					indx = mid+1;
				} else if(m>in_st) {
					indx = mid;
				} else if(m==in_st) {
					indx = mid+1;
				}
				break;
			}
			if(m<in_st) {
				i = mid+1;
			} else if(m>in_st) {
				j = mid-1;
			} else if(m==in_st) {
				indx = mid-1;
				break;
			}
		}
		if(i>j) indx = i;
	}
	
	if(indx<0 || indx>=sz(vstems))
		return stem(-1,-1,-1,-1);
	else 
		return vstems[indx];
	
}

/***********************************************************
                         seglist
***********************************************************/
segment seglist::findRS(const segment & in_sg) {
	int indx = -1;
	
	if(sz(vsegs)!=0) {
		int i = 0; int j = sz(vsegs)-1;
		while(i<=j) {
			int mid = (i+j)/2;
			segment m = vsegs[mid];
			if(i==j) {
				if(m<in_sg) {
					indx = mid+1;
				} else if(m>in_sg) {
					indx = mid;
				} else if(m==in_sg) {
					indx = mid+1;
				}
				break;
			}
			if(m<in_sg) {
				i = mid+1;
			} else if(m>in_sg) {
				j = mid-1;
			} else if(m==in_sg) {
				indx = mid-1;
				break;
			}
		}
		if(i>j) indx = i;
	}
	
	if(indx<0 || indx>=sz(vsegs))
		return segment(-1,-1);
	else 
		return vsegs[indx];
	
}

/***********************************************************
                         RStack
***********************************************************/
std::string RStack::tostring(){
	std::string ret;
	FOR(i,0,len) ret += ".";
	FOR(i,0,sz(vstems)){
//l1<=l2<r2<r1
		int l1 = vstems[i].outer.lp;
		int l2 = vstems[i].inner.lp;
		int r1 = vstems[i].outer.rp;
		int r2 = vstems[i].inner.rp;
		FOR(j,l1,l2+1){
			ret[j] = '('; ret[l1+r1-j] = ')';
		}
	}
	return ret;
}
//check whether this instance contains std::stack st?
bool RStack::containST(const stem & st) {
	FOR(i,0,sz(vstems)) {
		if(vstems[i] ==  st) return true;
	}
	return false;
}

int RStack::inwhichST(const int & query) { 
	//in which std::stack is query.
	FOR(i,0,sz(vstems)) {
		stem st = vstems[i];
		int x = st.outer.lp; int y = st.outer.rp;
		int p = st.inner.lp; int q = st.inner.rp;
		if((x<=query && query<=p) || (q<=query && query<=y)) 
			return i;
	}
	return -1;
}

//isLocOpt checks whether this instance is local minimum or not.
bool RStack::isLocOpt(const stemlist & in_stl) {
	stemlist rs_stl;
	rs_stl.vstems = vstems;
	FOR(i,0,sz(in_stl.vstems)){
		stem sti = in_stl.vstems[i];
		if(containST(sti)) continue;
		if(!rs_stl.conflictwith(sti)){
			return false;
		}
	}
	return true;
}

bool RStack::isLocOpt(const stemlist & in_stl, const std::vector<std::vector<int> > & H, 
const std::vector<std::vector<int> > & C) {
/*	stemlist rs_stl;
	rs_stl.vstems = vstems;
	FOR(i,0,sz(in_stl.vstems)){
		stem sti = in_stl.vstems[i];
		if(containST(sti)) continue;
		if(!rs_stl.conflictwith(sti)){
			if(H[sti.outer.lp][sti.outer.rp] <= 0 ||
			   C[sti.outer.lp][sti.outer.rp] <= 0 ) {
				return false;
			}
		}
	}
	return true; */
	return isLocOpt(in_stl);
}

std::vector<RStack> RStack::getneighbor_delstack() {
	std::vector<RStack> ret;
	FOR(i,0,sz(vstems)) {
		std::vector<stem> vst;
		FOR(j,0,sz(vstems)) {
			if(i==j) continue;
			vst.PB(vstems[j]);
		}
		RStack tmp;
		tmp.vstems = vst;
		ret.PB(tmp);
	}
	return ret;
}

bool RStack::ismatch(const std::string & seq) {
	FOR(i,0,sz(vstems)) {
		stem s = vstems[i];
		int sx = s.outer.lp; int sy = s.outer.rp;
		int sp = s.inner.lp; int sq = s.inner.rp;
		FOR(k,sx,sp+1) {
			char x = seq[k]; char y = seq[sx+sy-k];
			if(!match(x,y)) {
				std::cerr<<"seq["<<k<<"]="<<seq[k]<<" and seq["<<sx+sy-k<<"]="<<seq[sx+sy-k]<<" does not pair."<<std::endl;
				return false;
			}
		}
	}
	return true;
}

std::vector<int> RStack::findConflictingST(stem & st) {
	//return indices of stems that conflict with sti
	std::vector<int> ret;
	FOR(i,0,sz(vstems)) {
		if(st.conflictwith(vstems[i])) {
			ret.PB(i);
		}
	}
	return ret; 
}

/***********************************************************
                         RStruct
***********************************************************/
void RStruct::show(){
	std::cout<<s<<std::endl;
	FOR(i,0,sz(vbp)) {vbp[i].show();std::cout<<std::endl;}	
}

//compute the distance between this and rs
//int RStruct::dist(const RStruct & rs)
//{
//	return sz(convertfrom(rs));
//}

//determine with which the query integer (bp) pairs in this structure; -1 not pair with any bp
int RStruct::pairwithwhich(const int & query) 
{
	FOR(i,0,sz(vbp)) {
		if(vbp[i].lp==query) {
			return vbp[i].rp;
		}
		if(vbp[i].rp==query) {
			return vbp[i].lp;
		}
	}
	return -1;
}


//delete base pair with index i from vbp and add lp and rp to unp
//bool false fail; true succeed;
bool RStruct::del_bp_with_indx(const int & i) {
	assert(sz(vbp) > i);
//	std::cout<<"\tin del_bp_with_indx, i="<<i<<std::endl;
	if(sz(vbp) <= i) return false;
	//if i >= sz(l), base pair with index i does not exist
	unp.PB(vbp[i].lp); unp.PB(vbp[i].rp);
	//update s (this should be done before update vbp)
	s[vbp[i].lp] = '.'; s[vbp[i].rp] = '.';

	//update vbp
	vbp.erase(vbp.begin()+i);
//	std::cout<<"del bp with index "<<i<<"\t"<<", std::string="<<s<<std::endl;
	return true;
}

//delete base pair with index i from vbp and add lp and rp to unp
//bool false fail; true succeed;
bool RStruct::del_bp(const BP & in_bp) {
	int indx = getindx(vbp, in_bp);
//	std::cout<<"\tin_bp=("<<in_bp.lp<<","<<in_bp.rp<<")"<<std::endl;
//	std::cout<<"\tindx of in_bp is "<<indx<<", vbp[indx]=("<<vbp[indx].lp<<","<<vbp[indx].rp<<")"<<std::endl;
//	std::cout<<"\tstructure is "<<s<<std::endl;
	if(indx==-1) return false;
	return del_bp_with_indx(indx);
}

//delete base pair with either left or right parenthesis equals to i
//bool false fail; true succeed;
bool RStruct::del_bp_contain(const int & i){
	FOR(k,0,sz(vbp)) 
		if(vbp[k].lp==i || vbp[k].rp==i)
			return del_bp_with_indx(k);
//	std::cout<<"del bp contain "<<i<<"\t"<<", std::string="<<s<<std::endl;
	return true;
}

//add base pair (i,j) to vbp and delete them from unp
//bool false fail; true succeed;
bool RStruct::add_bp(const int & i, const int & j) {
	BP x(i,j);
	int indxi = getindx(unp, i);
	if(indxi == -1) return false;
	unp.erase(unp.begin()+indxi);

	int indxj = getindx(unp, j);
	if(indxj == -1) return false;
	unp.erase(unp.begin()+indxj);

	vbp.PB(x);
	
	//update s
	s[i] = '('; s[j] = ')';

//	std::cout<<"add bp "<<i<<","<<j<<", std::string="<<s<<std::endl;
	return true;
}

//add base pair in_bp to vbp and delete in_bp.lp and in_bp.rp from unp 
bool RStruct::add_bp(const BP & in_bp) {
	return add_bp(in_bp.lp, in_bp.rp);
}

bool RStruct:: operator == (const RStruct & rs) {
	if (rs.len != len ) return false;
	if (rs.unp != unp ) return false;
	if (sz(rs.vbp) != sz(vbp) ) return false;
	FOR(i,0,sz(vbp)) {
		if(getindx(rs.vbp, vbp[i]) == -1) return false;
	}
	return true;
}

bool RStruct:: operator != (const RStruct & rs) {
	if (rs.len != len ) return true;
	if (rs.unp != unp ) return true;
	if (sz(rs.vbp) != sz(vbp) ) return true;	
	FOR(i,0,sz(vbp)) {
		if(getindx(rs.vbp, vbp[i]) == -1) return true;
	}
	return false;
}

//apply action to this, but don't check clash
bool RStruct::applyAction(const Action &  in_a){ 
//	std::cout<<"in applyAction(Action)"<<std::endl;

	if(in_a.a == add) {
//		std::cout<<"action is add "<<in_a.b.lp<<","<<in_a.b.rp<<std::endl;
		if(!add_bp(in_a.b)) 
			return false;
	}
	else if(in_a.a == del) {
//		std::cout<<"action is del "<<in_a.b.lp<<","<<in_a.b.rp<<std::endl;
		if(!del_bp(in_a.b))
			return false;
	}
	else {
		std::cerr<<"Action = "<<in_a.a<<std::endl;
		std::cerr<<"The action is neither add nor del"<<std::endl;exit(0);
	}
	return true;
}
//check if in_a leads to clash if applied
//true clash; false not clash
bool RStruct::checkclash(const Action & in_a){
	if(in_a.a == del) return false;
	if (in_a.a == add) 
			return mycheckclash(in_a.b, vbp); 
	else {
		std::cerr<<"Action = "<<in_a.a<<std::endl;
		std::cerr<<"This action is neither add nor del"<<std::endl;exit(0);
	}
	
}
//get paired bases that are in both this and rs
std::vector<BP> RStruct::intersectBP(const RStruct & rs){
	std::vector<BP> ret;
	FOR (i,0,sz(vbp)) {
		bool find = false;
		FOR (j,0,sz(rs.vbp)) {
			if(vbp[i] == rs.vbp[j]) {
				find = true; break;
			}
		}
		if(find) 
			ret.PB(vbp[i]);
	}
	return ret;
}

//get BPs in this RStruct that are conflicting with in_bp
std::vector<BP> RStruct::getconflictingBP(const BP & in_bp) {
	std::vector<BP> ret;
	int i=in_bp.lp, j=in_bp.rp;
	FOR(t,0,sz(vbp)) {
		int h = vbp[t].lp, k = vbp[t].rp;
		if( (i<h && h<j && j<k) //i   h   j  k   || h  i  k  j
		  ||(h<i && i<k && k<j) ) {
			//conflicting
			ret.PB(vbp[t]);
		}
	}
	return ret;
}

//identify stems in the secondary structure
std::vector<stem> RStruct::create_stems() {
	std::vector<stem> ret;
	if (sz(vbp)==0) return ret;
	sort(vbp.begin(), vbp.end(), ascendingsortBP);
//	FOR(i,0,sz(vbp)){generate_all_RStacks
//		vbp[i].show(); std::cout<<" ";
//	}
	BP outer = vbp[0];
	BP inner = vbp[0];
	int i=1;
	while(i<sz(vbp)) {
		BP bpx = vbp[i-1];
		BP bpy = vbp[i];
		if (bpx.lp==bpy.lp-1 && bpx.rp==bpy.rp+1) {
			//continue to read the std::stack
		} else {
			ret.PB(stem(outer, inner));
			outer = bpy;
		}
		inner = bpy;
		i++;
	}
	ret.PB(stem(outer, inner));
	return ret;
}

/***********************************************************
                         PartStruct
***********************************************************/
bool PartStruct::isLocOpt(const stemlist & in_stl){
	FOR(i,0,sz(in_stl.vstems)){
		stem sti = in_stl.vstems[i];
		if(stl.containST(sti)) continue;
		if(!stl.conflictwith(sti) && !sgl.overlap(sti))
			return false;
	}
	return true;
}

















