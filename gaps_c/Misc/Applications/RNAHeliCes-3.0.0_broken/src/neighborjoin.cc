#include "neighborjoin.hh"
/*******************************************************************************
                             mytable
*******************************************************************************/
mytable::mytable(std::vector<std::vector<double> > & in_D) {
	original_num = valid_num = sz(in_D);
	int n = sz(in_D);
	D = in_D;
	FOR(i,0,n) {
		assert(sz(in_D[i]) == i);
	}
	
	FOR(i,0,n) {
		valid.PB(true);
		tree.PB(TNode());
	}
}




void mytable::pickminD(int & indx_i, int & indx_j) 
{	//Pick a pair i,j for which Q[i][j] is the smallest among all the cells that are valid.
	int n = sz(D);
	double minval = 0;
	indx_i = indx_j = -1;
	bool start = false;
	FOR(i,0,n) {
		if(!valid[i]) continue;
		FOR(j,0,i) {
			if(!valid[j]) continue;
			if(!start) {
				start = true;
				minval = D[i][j];
				indx_i = i;
				indx_j = j;
			} else {
				if(minval > D[i][j]) {
					minval = D[i][j];
					indx_i = i;
					indx_j = j;
				}
			}
		}
	}
	if(!start) {
		std::cerr<<"The overall matrix has been processed."<<std::endl;
		exit(0);
	}
}


//join neighbors according to the distance between structures.
std::vector<int> mytable::DJ(const int & knum, const std::vector<RNA> & vrna, std::ostream & fout, std::vector<int> & vorder) {
	int delta = 0; 
	int n = original_num;
	vorder.clear();
	if(n==1) {
		vorder.PB(0);
	}
	while(n>1) {
		int indx_i = -1;
		int indx_j = -1; 
		pickminD(indx_i, indx_j);
		
		//for debugging DJ
//		if(n<=10) {
//			fout<<"n = "<< n << ", display the current structures."<<std::endl;
//			FOR(i,0,sz(D)) {
//				if(!valid[i]) continue;
//				fout<<i<<"\t"<<vrna[i].str<<std::endl;
//			}
//			fout<<"Display the current matrix."<<std::endl;
//			FOR(i,0,sz(D)) {
//				if(!valid[i]) continue;
//				FOR(j,0,sz(D[i])) {
//					if(!valid[j]) continue;
//					fout.precision(4);
//					fout<<D[i][j]<<"\t";
//				}
//				fout<<std::endl;
//			}
//		}
//		fout<<"D["<<indx_i<<"]["<<indx_j<<"]="<<D[indx_i][indx_j]<<std::endl;
//		fout<<"join structure["<<indx_i<<"] and structure["<<indx_j<<"], the distance between them is "<<D[indx_i][indx_j]<<std::endl;
//		fout<<vrna[indx_i].str<<std::endl;
//		fout<<vrna[indx_j].str<<std::endl;
		assert(indx_i>indx_j);
		
		int indx_rm = -1;
		if(vrna[indx_i].e < vrna[indx_j].e) {
			tree[indx_j].pt = indx_i;
			valid[indx_j] = false;
//			fout<<"structure["<<indx_j<<"]-->structure["<<indx_i<<"]"<<std::endl;		
			indx_rm = indx_j;
		} else {
			tree[indx_i].pt = indx_j;
			valid[indx_i] = false;
//			fout<<"structure["<<indx_i<<"]-->structure["<<indx_j<<"]"<<std::endl;		
			indx_rm = indx_i;
		}
		//## 删除小的那个记录
		vorder.PB(indx_rm);
		if(n==2) {
			vorder.PB(0);
		}
//		fout<<"index="<<indx_rm<<", structure="<<vrna[indx_rm].str<<std::endl;
		n--;
	}

	//## 颠倒 vector, 如是有k值，取k值
	std::vector<int> ret;
	FOR(k,0,std::min(knum, sz(vorder))) {
		ret.PB(vorder[sz(vorder)-1-k]);
	}

	return ret;
} 
	

void mytable::show() {
	std::cout<<"The original num is "<<original_num<<std::endl;
	std::cout<<"The valid num is "<<valid_num<<std::endl;
	//print out D
	std::cout<<"D"<<std::endl;
	FOR(i,0,sz(D)) {
		std::cout<<std::setw(6)<<i<<" ";
		FOR(j,0,sz(D[i])) {
			std::cout.precision(2);
			std::cout<<std::setw(6)<<D[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;

/*
	//print out Q
	std::cout<<"Q"<<std::endl;
	FOR(i,0,sz(Q)) {
		std::cout<<setw(8)<<i<<" ";
		FOR(j,0,sz(Q[i])) {
			std::cout.precision(3);
			std::cout<<setw(10)<<Q[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl; */

	std::cout<<"valid"<<std::endl;
	//print out valid
	FOR(i,0,sz(valid)) {
		std::cout<<"["<<i<<"] = ";
		if(valid[i]) 
			std::cout<<"T"<<std::endl;
		else 
			std::cout<<"F"<<std::endl;
	} 
	std::cout<<std::endl;
	//print out tree
	std::cout<<"Tree"<<std::endl;
	showtree(sz(tree)-1, 0);
	std::cout<<std::endl;	
}

void mytable::showtree(const int & n, const int & indent) {

}

/*
void mytable::initQ() 
{	//update Matrix Q based on Matrix D
	int n = original_num;
	//construct Q
	FOR(i,0,n) {
		std::vector<double> vd;
		Q.PB(vd);
		double drowi = 0;
		FOR(j,0,n) {
			if(j<i) {
				drowi += D[i][j];
				Q[i].PB(0);
			} else if(i==j) {
				//do nothing
			} else {
				drowi += D[j][i];
			}
		}
		drow.PB(drowi);
	}
	
	FOR(i,0,n) {
		FOR(j,0,i) {
			Q[i][j] = (original_num-2)*D[i][j] - drow[i] - drow[j];
		}
	}
}


void mytable::updateQdrow(const int & indx_rm)
{	//Only one of i and j are removed from the matrix,	update Q and drow accordingly
	int u = indx_rm;
	
	std::cout<<"Update Q"<<std::endl;
	FOR(k,0,sz(Q)) {
		if(!valid[k]) continue;
		FOR(h,0,k) {
			if(!valid[h]) continue;
			Q[k][h] += drow[k] + drow[h] - D[k][h];
		}
	}
	

	FOR(k,0,sz(drow)) {
		if(!valid[k]) continue;
		drow[k] += - ((k<u)?(D[u][k]):(D[k][u]));
		//drow[k] += -D[k][u];
	}

		
	FOR(k,0,sz(Q)) {
		if(!valid[k]) continue;
		FOR(h,0,k) {
			if(!valid[h])
				Q[k][h] = 0;
			else 
				Q[k][h] += -drow[k] - drow[h];
		}
	}
} */
/*
void mytable::pickminQ(int & indx_i, int & indx_j) 
{	//Pick a pair i,j for which Q[i][j] is the smallest among all the cells that are valid.
	int n = sz(Q);
	double minval = 0;
	indx_i = indx_j = -1;
	bool start = false;
	FOR(i,0,n) {
		if(!valid[i]) continue;
		FOR(j,0,i) {
			if(!valid[j]) continue;
			if(!start) {
				start = true;
				minval = Q[i][j];
				indx_i = i;
				indx_j = j;
			} else {
				if(minval > Q[i][j]) {
					minval = Q[i][j];
					indx_i = i;
					indx_j = j;
				}
			}
		}
	}
	if(!start) {
		std::cerr<<"The overall matrix has been processed."<<std::endl;
		exit(0);
	}
}
*/
//This is a variant of the neighbor joining algorithm of Saitou 
//Saitou's NJ algorithm:
//(1) compute a low triangular matrix Q based on the low triangular distance matrix D
//(2) pick up a pair i,j for which Q[i][j] is the smallest
//(3) create a new node u and infer the distance
//    u represents the consensus of i and j
//(4) compute diu, the distance between i and u 
//            dju, the distance between j and u
//    for all the other nodes k, duk, the distance between k and u
//    and update the distance matrix D
//(5) repeat (1)-(4)
//The variant of the Saitou algorithm
//(1) compute a low triangular matrix Q based on the low triangular distance matrix D
//(2) pick up a pair i,j for which Q[i][j] is the smallest
//(3) choose from i and j the one with the lowest free energy and delete the other one from the matrix D
//(4) repeat (1)-(3)
/*
std::vector<int> mytable::NJ(const int & knum, const std::vector<RNA> & vrna, ofstream & fout)
{
	int n = original_num;
	std::cout<<"original number is "<<n<<std::endl;
	while(n>=knum) {
		std::cout<<"n = "<<n<<std::endl;
		//pick a pair i,j for which Q[i][j] is the smallest 
		int indx_i = -1;
		int indx_j = -1;
		pickminQ(indx_i, indx_j);
		fout<<"Q["<<indx_i<<"]["<<indx_j<<"]="<<Q[indx_i][indx_j]<<std::endl;
		fout<<"join structure["<<indx_i<<"] and structure["<<indx_j<<"], the distance between them is "<<D[indx_i][indx_j]<<std::endl;
		fout<<vrna[indx_i].str<<std::endl;
		fout<<vrna[indx_j].str<<std::endl;
	
		assert(indx_i>indx_j);
		
		show();
		//(3) choose from i and j the one with the lowest free energy and delete the other one from the matrix D
		int indx_rm = -1;
		if(vrna[indx_i].e < vrna[indx_j].e) {
			tree[indx_j].pt = indx_i;
			valid[indx_j] = false;
			fout<<"structure["<<indx_j<<"]-->structure["<<indx_i<<"]"<<std::endl;		
			indx_rm = indx_j;
		} else {
			tree[indx_i].pt = indx_j;
			valid[indx_i] = false;
			fout<<"structure["<<indx_i<<"]-->structure["<<indx_j<<"]"<<std::endl;		
			indx_rm = indx_i;
		}

		//update n
		n--;
		
		//update Q and drow
		std::cout<<"update Q and drow"<<std::endl;
		if(n>=knum) 
			updateQdrow(indx_rm);
	}
	
	std::vector<int> ret;
	FOR(k,0,sz(valid)) {
		if(valid[k]) 
			ret.PB(k);
	}
	return ret;
} */


