#ifndef FILTER_ADD_HH
#define FILTER_ADD_HH

//#include "rtlib/rope.hh"
#include "Extensions/alifold.hh"

//#include "string.hh"
//#include "sequence.hh"

//#include <cassert>
#include <iostream>
#include <vector>


// template<typename alphabet, typename pos_type, typename T>
// inline bool exact_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
//     T i, T j, const std::string &str)      // modify line 1: match_str ==> str
// {
//     //std::cout << "line" << std::endl;
//     ////std::vector<std::string> tokens;  // comparing with using function mytokenize comment line 1
//     ////mytokenize(match_str, tokens);    // comment line 2
//     
//     std::vector<std::string> tokens;      // add line 1
//     const std::string delimiters = ",";   // add line 2
//     // Skip delimiters at beginning.
//     std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
//     // Find first "non-delimiter".
//     std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
// 
//     while (std::string::npos != pos || std::string::npos != lastPos)
//     {
//         // Found a token, add it to the vector.
//         tokens.push_back(str.substr(lastPos, pos - lastPos));
//         // Skip delimiters.  Note the "not_of"
//         lastPos = str.find_first_not_of(delimiters, pos);
//         // Find next "non-delimiter"
//         pos = str.find_first_of(delimiters, lastPos);
//     }
// 
// 
//     std::vector<std::string>::iterator it;
//     //std::cout << tokens.size() << std::endl;
//     for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
//       std::string tmp = *it;
//       if (tmp[tmp.length()-1] == 'm' || tmp[tmp.length()-1] == 'b') {
//           tokens.erase(it);
// 	  it -= 1;
//       }
//      // std::cout << *it << std::endl;
//       //(*it).at((*it).size()-1);
//     }
// 
//     //std::copy(tokens.begin(), tokens.end(), std::ostream_iterator<std::string>(std::cout, ":yangyang: "));
//     ////11  std::cout << "ok:reference_hishape_=" << reference_hishape_ << "  x.first=" << x.first << "  x.second" << x.second << std::endl;
//     
//     //TODO: assert(sth.)
//     //TODO: improve the runtime in the following piece of code! using string may be a good idea.
// 	        Rope helix_index;
// 		int helix_center_int;
// 		helix_center_int = (i+j+1)/2;
// 		if ( helix_center_int*2 > (i+j+1) ) helix_center_int = helix_center_int - 1;  
// 		append(helix_index, helix_center_int);
// 		if ( helix_center_int*2 != (i+j+1) ) append(helix_index, ".5", 2);
//       std::ostringstream outs;    // Declare an output string stream.
//       outs << helix_index;            // Convert value into a string.
//       std::string helix_index_str = outs.str();     // Get the created string from the output stream.
// 
//     bool found = false;
//     for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
//         if (helix_index_str == *it) {
// 	    found = true;
// 	    break;
// 	}
//     }
//     
//     return (i+3 < j) && basepairing(seq, i, j) && basepairing(seq, i+1, j-1) && found;  //((i+j+1)/2 == l);
// }

// template<typename alphabet, typename pos_type, typename T>
// inline bool exact_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
//     T i, T j, std::vector<std::string>& tokens)      // modify line 1: match_str ==> str
// {
//     std::vector<std::string>::iterator it;
//   
//     Rope helix_index;
//     int helix_center_int;
//     helix_center_int = (i+j+1)/2;
//     if ( helix_center_int*2 > (i+j+1) ) helix_center_int = helix_center_int - 1;  
//     append(helix_index, helix_center_int);
//     if ( helix_center_int*2 != (i+j+1) ) append(helix_index, ".5", 2);
//     std::ostringstream outs;    // Declare an output string stream.
//     outs << helix_index;            // Convert value into a string.
//     std::string helix_index_str = outs.str();     // Get the created string from the output stream.
// 
//     bool found = false;
//     for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
//         if (helix_index_str == *it) {
// 	    found = true;
// 	    break;
// 	}
//     }
//     
//     return (i+3 < j) && basepairing(seq, i, j) && basepairing(seq, i+1, j-1) && found;  //((i+j+1)/2 == l);
// }



template<typename alphabet, typename pos_type, typename T>
inline bool exact_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
    T i, T j, std::vector<std::string>& tokens)
{
    if (tokens.empty())
    {
        return basepairing(seq, i, j);
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
	
	return basepairing(seq, i, j) && found;
    }
}

template<typename alphabet, typename pos_type, typename T>
inline bool fuzzy_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
    T i, T j, std::vector<std::string>& tokens)
{
    if (tokens.empty())
    {
        return basepairing(seq, i, j);
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
	
	return basepairing(seq, i, j) && found;
    }
}

#endif
