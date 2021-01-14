#ifndef _NEIGHBORJOIN_HH_
#define _NEIGHBORJOIN_HH_
#include "RNAStack.hh"


class TNode {
public:
	int pt; //point to 
	TNode() {pt = -1;}
	TNode(const int & in_pt) {pt = in_pt;}
	~TNode() {}
	void reset() {pt = -1;} 
};

class mytable {
public:
	int original_num;  //the number of original OTUs
	int valid_num;     //the number of OTUs remain valid
	std::vector<bool> valid;//indicate whether column x and row x are valid, e.g. valid[i]==false --> matrix[i][.] and maxtrix[j][.] are invalid
	std::vector<TNode> tree;
	std::vector<std::vector<double> > D; //the distance matrix
//	std::vector<std::vector<double> > Q; //the Q matrix, computed based on D
	std::vector<double> drow;      //the sum of distances in each row
	
	mytable() {original_num=0;}
	mytable(std::vector<std::vector<double> > & in_D);
	~mytable() {}
//	std::vector<int>   NJ(const int & knum, const std::vector<RNA> & vrna, ofstream & fout); //A variant of the neighbor joining algorithm of Saitou
//	std::vector<int>   DJ(const int & knum, const std::vector<RNA> & vrna, ofstream & fout); //join neighbors according to the distance between structures.
	std::vector<int>   DJ(const int & knum, const std::vector<RNA> & vrna, std::ostream & cout, std::vector<int> & vorder); //join neighbors according to the distance between structures.
	void   show();
	void   showtree(const int & n, const int & indent);
	std::vector<int> cluster2k(const int & k);//split the whole tree into k parts
		
private:
//	void   initQ();//compute matrix Q based on matrix D
//	void   updateQdrow(const int & u); //Only one of i and j are removed from the matrix,	update Q and drow accordingly
	void   pickminD(int & indx_i, int & indx_j);
//	void   pickminQ(int & indx_i, int & indx_j);//Pick a pair i,j for which matrix[i][j] is the smallest among all the cells that are valid.	
};
#endif
