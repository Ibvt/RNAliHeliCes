#ifndef ALIFOLD_HH
#define ALIFOLD_HH

#include "rnaoptions_defaults.hh"
#include "alignment.hh"

//static const float cfactor = 1.0; //Set the weight of the covariance term in the energy function (default=`1.0')
//static const float nfactor = 1.0; //Set the penalty for non-compatible sequences in the covariance term of the energy function (default=`1.0')
//static const int MINPSCORE = -200;

template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
  if (j<=i+1) {
    return false;
  }
  Basic_Subsequence<alphabet, pos_type> sub(seq, i, j);
  return int(covscore(sub, int(i), int(j)-1)*-1*rows(sub)) >= int(getAlifold_minscore_basepair()); //convert to int, because float differences are very small, but will have a hugh impact on small changes of nfactor or cfactor!
}

//basepair filter for an un-interrupted stem, as they appear in pseudoknots for alpha, beta and gamma helices.
inline bool regionpair(int i, int j, int len) {
	return true;
}

template<typename alphabet, typename pos_type, typename T>
inline bool unpaired(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
	return true;
}


template<typename alphabet, typename pos_type, typename T>
inline bool fuzzy_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
    T i, T j, std::vector<std::string>& tokens)
{
    if (j<=i+1) {
        return false;
    }
    Basic_Subsequence<alphabet, pos_type> sub(seq, i, j);
    
    if (tokens.empty())
    {
        return int(covscore(sub, int(i), int(j)-1)*-1*rows(sub)) >= int(getAlifold_minscore_basepair());
    }
    else
    {
	std::vector<std::string>::iterator it;
	float helix_center_float = (i+j+1)/2.0f;

	bool found = false;
	float fuzzy = 0.0f;
	for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
	    fuzzy = atof((*it).c_str());
	    if ((helix_center_float >= (fuzzy-1.5f)) && (helix_center_float <= (fuzzy+1.5f))) {
		found = true;
		break;
	    }
	}
	return (int(covscore(sub, int(i), int(j)-1)*-1*rows(sub)) >= int(getAlifold_minscore_basepair())) && found;
    }
}

template<typename alphabet, typename pos_type, typename T>
inline bool exact_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
    T i, T j, std::vector<std::string>& tokens)
{
    if (j<=i+1) {
        return false;
    }
    Basic_Subsequence<alphabet, pos_type> sub(seq, i, j);
    
    if (tokens.empty())
    {
        return int(covscore(sub, int(i), int(j)-1)*-1*rows(sub)) >= int(getAlifold_minscore_basepair());
    }
    else
    {
	std::vector<std::string>::iterator it;
	float helix_center_float = (i+j+1)/2.0f;

	bool found = false;
	float fuzzy = 0.0f;
	for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
	    fuzzy = atof((*it).c_str());
	    if ((helix_center_float >= (fuzzy-0.1f)) && (helix_center_float <= (fuzzy+0.1f))) {
		found = true;
		break;
	    }
	}
	return (int(covscore(sub, int(i), int(j)-1)*-1*rows(sub)) >= int(getAlifold_minscore_basepair())) && found;
    }
}

#endif


