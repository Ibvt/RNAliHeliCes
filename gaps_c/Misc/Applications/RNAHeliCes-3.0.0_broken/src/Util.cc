#include "RNAStackFold.hh"

std::string dna2rna(const std::string & seq) {
	std::string ret = "";
	FOR(i,0,sz(seq)){
		if(seq[i] == 'T'|| seq[i] == 't' || seq[i] == 'U' || seq[i] == 'u') {
			ret += "U";
		} else if (seq[i] == 'A' || seq[i] == 'a') {
			ret += "A";
		} else if (seq[i] == 'G' || seq[i] == 'g') {
			ret += "G";
		} else if (seq[i] == 'C' || seq[i] == 'c') {
			ret += "C";
		} else {
			std::cerr<<"Unknown character "<<seq[i]<<" found in an RNA sequence: "<< 
			  seq<< std::endl;
			exit(0);
		}
	}
	return ret;
}

double max(const double & d1, const double & d2){
	return (d1>d2)?(d1):(d2);
}

double min(const double & d1, const double & d2){
	return (d1>d2)?(d2):(d1);
}

int max(const int & d1, const int & d2){
	return (d1>d2)?(d1):(d2);
}

int min(const int & d1, const int & d2){
	return (d1>d2)?(d2):(d1);
}

void init_vi(std::vector<int> & vi, const int & n, const int & val) {
	vi.clear();
	FOR(i,0,n)
		vi.PB(val);
}

void init_vvi(std::vector<std::vector<int> > & vvi, const int & n, const int & m, const int & val) {
	vvi.clear();
	FOR(i,0,n) {
		std::vector<int> vi;
		vvi.PB(vi);
		FOR(j,0,m)
			vvi[i].PB(val);
	}
}

void init_vstem(std::vector<stem> & vstem, const int & n, stem & val) {
	vstem.clear();
	FOR(i,0,n)
		vstem.PB(val);
}

void init_vvstem(std::vector<std::vector<stem> > & vvstem, const int & n, const int & m, stem & val) {
	vvstem.clear();
	FOR(i,0,n) {
		std::vector<stem> vstem;
		vvstem.PB(vstem);
		FOR(j,0,m)
			vvstem[i].PB(val);
	}
}

void init_vsfa(std::vector<StackFoldAuxiliary> & vsfa, const int & n, const StackFoldAuxiliary & val) {
	vsfa.clear();
	FOR(i,0,n)
		vsfa.PB(val);
}

void init_vvsfa(std::vector<std::vector<StackFoldAuxiliary> > & vvsfa, const int & n, const int & m, const StackFoldAuxiliary & val) {
	vvsfa.clear();
	FOR(i,0,n) {
		std::vector<StackFoldAuxiliary> vsfa;
		vvsfa.PB(vsfa);
		FOR(j,0,m)
			vvsfa[i].PB(val);
	}
}

void showETYPE(const ETYPE & etype, std::ostream & fout) {
/*	etype_unknown = 0,
	etype_F = 1,
	etype_C = 2,
	etype_FM = 3,
	etype_M1 = 4,
	etype_M2 = 5,
	etype_M4 = 6,
	etype_Fn = 7,
	etype_Fd = 8,
	etype_F2n = 9,
	etype_F2d = 10*/
	switch(etype) {
		case etype_unknown: {fout<<"unknown";break;}
		case etype_F:{fout<<"F ";break;}
		case etype_C:{fout<<"C ";break;}
		case etype_FM:{fout<<"FM ";break;}
		case etype_M1:{fout<<"M1 ";break;}
		case etype_M2:{fout<<"M2 ";break;}
		case etype_M4:{fout<<"M4 ";break;}
		case etype_Fn:{fout<<"Fn ";break;}
		case etype_Fd:{fout<<"Fd ";break;}
		case etype_F2n:{fout<<"F2n ";break;}
		case etype_F2d:{fout<<"F2d ";break;}
		default:{std::cerr<<"unknown ETYPE="<<etype<<std::endl;exit(0);}
	}
}

void showCTYPE(const CTYPE & ctype, std::ostream & fout) {
/*	ctype_unknown = -1,
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
	ctype_5b = 23 */
	switch(ctype) {
		case ctype_unknown: {fout<<"unknown";break;}
		case ctype_0: {fout<<"0";break;}
		case ctype_1: {fout<<"1";break;}
		case ctype_2: {fout<<"2";break;}
		case ctype_3: {fout<<"3";break;}
		case ctype_4: {fout<<"4";break;}
		case ctype_1a: {fout<<"1.a";break;}
		case ctype_1b: {fout<<"1.b";break;}
		case ctype_1c: {fout<<"1.c";break;}
		case ctype_2a: {fout<<"2.a";break;}
		case ctype_2b: {fout<<"2.b";break;}
		case ctype_2c: {fout<<"2.c";break;}
		case ctype_3a: {fout<<"3.a";break;}
		case ctype_3b: {fout<<"3.b";break;}
		case ctype_3c: {fout<<"3.c";break;}
		case ctype_3d: {fout<<"3.d";break;}
		case ctype_3e: {fout<<"3.e";break;}
		case ctype_3f: {fout<<"3.f";break;}
		case ctype_3g: {fout<<"3.g";break;}
		case ctype_3h: {fout<<"3.h";break;}
		case ctype_3i: {fout<<"3.i";break;}
		case ctype_4a: {fout<<"4.a";break;}
		case ctype_4b: {fout<<"4.b";break;}
		case ctype_5a: {fout<<"5.a";break;}
		case ctype_5b: {fout<<"5.b";break;}
		default:{std::cerr<<"unknown CTYPE="<<ctype<<std::endl;exit(0);}
	}
}

void show(const std::vector<ETYPE> & vet, std::ostream & fout){
	FOR(i,0,sz(vet)) {
		showETYPE(vet[i], fout);
	}
	fout<<std::endl;
}

void show(const std::vector<int> & vi, std::ostream & fout){
	FOR(i,0,sz(vi)) {
		if(vi[i] != OOE) {
			std::cout.width(OUTWIDTH);std::cout<<vi[i];
		} else {
			std::cout.width(OUTWIDTH);std::cout<<"INF";
		}
	}
	std::cout<<std::endl;
}

void show(const std::vector<double> & vd, std::ostream & fout){
	FOR(i,0,sz(vd)) {
		if(!eq(vd[i],OOE)) {
			std::cout.width(OUTWIDTH);std::cout<<vd[i];
		} else {
			std::cout.width(OUTWIDTH);std::cout<<"INF";
		}
	}
	std::cout<<std::endl;
}


void show(const std::vector<std::string> & vs, std::ostream & fout){
	FOR(i,0,sz(vs)) {
		std::cout.width(OUTWIDTH);std::cout<<vs[i];
	}
	std::cout<<std::endl;
}

void print2file(const std::vector<std::vector<int> > & vvi, const std::string & file) {
	std::ofstream fout(file.c_str());
	show(vvi, fout);
	fout.close();
}

void show(const std::vector<std::vector<int> > & vvi, std::ostream & fout){
	FOR(i,0,sz(vvi)) {
		FOR(j,0,sz(vvi)){
			if(vvi[i][j] != OOE) {
//				std::cout.width(OUTWIDTH);std::cout<<vvi[i][j];
				fout<<vvi[i][j];
			} else {
//				std::cout.width(OUTWIDTH);std::cout<<"INF";
				fout<<"X";
			}
			if(j!=sz(vvi)-1) fout<<"\t";
		}
		fout<<std::endl;
	}
	fout<<std::endl;
}

void show(const std::vector<std::vector<double> > & vvd, std::ostream & fout){
	FOR(i,0,sz(vvd)) {
		FOR(j,0,sz(vvd)){
			if(!eq(vvd[i][j], OOE)) {
				std::cout.width(OUTWIDTH);std::cout<<vvd[i][j];
			} else {
				std::cout.width(OUTWIDTH);std::cout<<"INF";
			}
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

void show(const std::vector<std::vector<std::string> > & vvs, std::ostream & fout){
	FOR(i,0,sz(vvs)) {
		FOR(j,0,sz(vvs)){
			std::string str = "";
			FOR(k,0,OUTWIDTH-sz(vvs[i][j])){
				str += " ";
			}
			std::cout<<str<<vvs[i][j];
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}






















