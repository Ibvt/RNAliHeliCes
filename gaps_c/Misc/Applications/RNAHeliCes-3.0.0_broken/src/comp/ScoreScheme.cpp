//*****************************************************************************
#include "SimpleIni.h"
#include "ScoreScheme.hpp"
#include "DebugTools.hpp"
#include <sstream>
#include <limits>

namespace hishape
{
namespace comp
{
  
    // ################ function #################
    // the third parameter of from_string() should be 
    // one of std::hex, std::dec or std::oct
    template <class T>
    bool from_string(T& t, 
		    const std::string& s, 
		    std::ios_base& (*f)(std::ios_base&))
    {
      std::istringstream iss(s);
      return !(iss >> f >> t).fail();
    }
    
    
CScoreScheme::CScoreScheme()
{
    // LOADING DATA

    // load from a data file
    CSimpleIniA ini(true, true, true);
    //SI_Error rc = ini.LoadFile("ScoreScheme.ini");
    //if (rc < 0) return false;
    //assert(rc>=0);
    
    // load from a string
    /*    
	std::string strData;
	rc = ini.LoadData(strData.c_str(), strData.size());
	//if (rc < 0) return false;
	assert(rc>=0);
    */

    float fB2iValue = 3.0f;
    float fH2iValue = 8.0f;
    float fH2bValue = 8.0f;
    float fM2iValue = 8.0f;
    float fM2bValue = 8.0f;
    float fM2hValue = 8.0f;
    float fIndel2iValue = 5.0f;
    float fIndel2bValue = 5.0f;
    float fIndel2hValue = 100.0f;
    float fIndel2mValue = 75.0f;    
    if (ini.LoadFile("ScoreScheme.ini")>=0 || ini.LoadFile("../ScoreScheme.ini")>=0)
    {
	// GETTING SECTIONS AND KEYS

	// get all sections
	//CSimpleIniA::TNamesDepend sections;
	//ini.GetAllSections(sections);

	// get all keys in a section
	//CSimpleIniA::TNamesDepend keys;
	//ini.GetAllKeys("integer", keys);

	
	
	// GETTING VALUES

	// get the value of a key
	const char * b2iValue = ini.GetValue("integer",
	    "B_I", NULL /*default*/);
	std::string strB2iValue(b2iValue);
	from_string<float>(fB2iValue, strB2iValue, std::dec);
	
	const char * h2iValue = ini.GetValue("integer",
	    "H_I", NULL /*default*/);
	std::string strH2iValue(h2iValue);
	from_string<float>(fH2iValue, strH2iValue, std::dec);
	
	const char * h2bValue = ini.GetValue("integer",
	    "H_B", NULL /*default*/);
	std::string strH2bValue(h2bValue);
	from_string<float>(fH2bValue, strH2bValue, std::dec);
	    
	const char * m2iValue = ini.GetValue("integer",
	    "M_I", NULL /*default*/);
	std::string strM2iValue(m2iValue);
	from_string<float>(fM2iValue, strM2iValue, std::dec);
	
	const char * m2bValue = ini.GetValue("integer",
	    "M_B", NULL /*default*/);
	std::string strM2bValue(m2bValue);
	from_string<float>(fM2bValue, strM2bValue, std::dec);
	
	const char * m2hValue = ini.GetValue("integer",
	    "M_H", NULL /*default*/);
	std::string strM2hValue(m2hValue);
	from_string<float>(fM2hValue, strM2hValue, std::dec);
	
	const char * indel2iValue = ini.GetValue("integer",
	    "Indel_I", NULL /*default*/);
	std::string strIndel2iValue(indel2iValue);
	from_string<float>(fIndel2iValue, strIndel2iValue, std::dec);
	
	const char * indel2bValue = ini.GetValue("integer",
	    "Indel_B", NULL /*default*/);
	std::string strIndel2bValue(indel2bValue);
	from_string<float>(fIndel2bValue, strIndel2bValue, std::dec);
	
	const char * indel2hValue = ini.GetValue("integer",
	    "Indel_H", NULL /*default*/);
	std::string strIndel2hValue(indel2hValue);
	from_string<float>(fIndel2hValue, strIndel2hValue, std::dec);
	
	const char * indel2mValue = ini.GetValue("integer",
	    "Indel_M", NULL /*default*/);
	std::string strIndel2mValue(indel2mValue);
	from_string<float>(fIndel2mValue, strIndel2mValue, std::dec);
     
    }
    
    float fI2bValue = fB2iValue;
    float fI2hValue = fH2iValue;
    float fB2hValue = fH2bValue;
    float fI2mValue = fM2iValue;
    float fB2mValue = fM2bValue;
    float fH2mValue = fM2hValue;
    float fI2indelValue = fIndel2iValue;
    float fB2indelValue = fIndel2bValue;
    float fH2indelValue = fIndel2hValue;
    float fM2indelValue = fIndel2mValue; 
    
  // std::numeric_limits<int>::max()

    m_iScoreMatrixWidth = 6;
    m_iScoreMatrixHeight = 6;
    /*
    float standardScoreMatrix[6][6]
    = { {0.0f,      std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity()},  // N
        {std::numeric_limits<float>::infinity(),0.0f,      3.0f,      8.0f,      8.0f,      5.0f      },  // I
        {std::numeric_limits<float>::infinity(),3.0f,      0.0f,      8.0f,      8.0f,      5.0f      },  // B
        {std::numeric_limits<float>::infinity(),8.0f,      8.0f,      0.0f,      8.0f,      100.0f    },  // H
        {std::numeric_limits<float>::infinity(),8.0f,      8.0f,      8.0f,      0.0f,      75.0f     },  // M       
        {std::numeric_limits<float>::infinity(),5.0f,      5.0f,      100.0f,    75.0f,     0.0f      }  };   // Indel
        */
    float standardScoreMatrix[6][6]
    = { {0.0f,std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity(),std::numeric_limits<float>::infinity()},  // N
        {std::numeric_limits<float>::infinity(),0.0f,          fB2iValue,     fH2iValue,     fM2iValue,     fIndel2iValue    },  // I
        {std::numeric_limits<float>::infinity(),fI2bValue,     0.0f,          fH2bValue,     fM2bValue,     fIndel2bValue    },  // B
        {std::numeric_limits<float>::infinity(),fI2hValue,     fB2hValue,     0.0f,          fM2hValue,     fIndel2hValue    },  // H
        {std::numeric_limits<float>::infinity(),fI2mValue,     fB2mValue,     fH2mValue,     0.0f,          fIndel2mValue    },  // M       
        {std::numeric_limits<float>::infinity(),fI2indelValue, fB2indelValue, fH2indelValue, fM2indelValue, -std::numeric_limits<float>::infinity()}  };   // Indel
        m_tScoreMatrix.resize(6);
	for (int i = 0; i < 6; i++)
	{
		m_tScoreMatrix[i].resize(6);
		for (int j = 0; j < 6; j++)
		{
			m_tScoreMatrix[i][j] = standardScoreMatrix[i][j];
		}
	}

/*
    m_iScoreMatrixWidth = 6;
    m_iScoreMatrixHeight = 6;
    int standardScoreMatrix[6][6]
    = { {0,1,1,1,1,1},  // N
        {1,0,1,1,1,1},  // I
        {1,1,0,1,1,1},  // B
        {1,1,1,0,1,1},  // H
        {1,1,1,1,0,1},  // M       
        {1,1,1,1,1,0}  };   // Indel
        m_tScoreMatrix.resize(6);
	for (int i = 0; i < 6; i++)
	{
		m_tScoreMatrix[i].resize(6);
		for (int j = 0; j < 6; j++)
		{
			m_tScoreMatrix[i][j] = standardScoreMatrix[i][j];
		}
	}
	*/
}


/**
Default constructor
*******************
Creates a @c iWidth x @c iWidth score matrix with each cell initialized to 0.

@param const int &iWidth:
score scheme width.

@param const int &iHeight:
score scheme height.

@return CScoreScheme:
 a new instance of a score scheme object.
*/
CScoreScheme::CScoreScheme(const int &iWidth, const int &iHeight)
: m_iScoreMatrixWidth(iWidth)
, m_iScoreMatrixHeight(iHeight)
{
	TRACEI("CScoreScheme::CScoreScheme(const int &iWidth, const int &iHeight)");
	// inits score matrix
	m_tScoreMatrix.resize(iWidth);
	for (int i = 0; i < iWidth; i++)
	{
		m_tScoreMatrix[i].resize(iHeight);
		for (int j = 0; j < iHeight; j++)
		{
			m_tScoreMatrix[i][j] = 0;
		}
	}
}


/**
Copy constructor
****************
Copy constructor.

@param const CScoreScheme &tScoreScheme:
score scheme to copy.

@return CScoreScheme:
 a new copy of tScoreScheme.
*/
CScoreScheme::CScoreScheme(const CScoreScheme &tScoreScheme)
: m_tScoreMatrix(tScoreScheme.m_tScoreMatrix)
, m_iScoreMatrixWidth(tScoreScheme.m_iScoreMatrixWidth)
, m_iScoreMatrixHeight(tScoreScheme.m_iScoreMatrixHeight)
{
	TRACEI("CScoreScheme::CScoreScheme(const CScoreScheme &tScoreScheme)");
}


/**
Destructor
**********
Standard virtual destructor.
*/
CScoreScheme::~CScoreScheme(void)
{
	TRACEI("CScoreScheme::~CScoreScheme(void)");
}


/**
@fn CScoreScheme &CScoreScheme::operator=(const CScoreScheme &tRightScoreScheme)
****************
Equal operator.

@param const CScoreScheme &tScoreScheme:
score scheme to copy.

@return CScoreScheme:
 the score scheme copy.
*/
CScoreScheme &CScoreScheme::operator=(const CScoreScheme &tRightScoreScheme)
{
	TRACEI("CScoreScheme &CScoreScheme::operator=(const CScoreScheme &tRightScoreScheme)");
	if (this == &tRightScoreScheme)
	{
		return *this;
	}
	m_tScoreMatrix = tRightScoreScheme.m_tScoreMatrix;
	m_iScoreMatrixWidth = tRightScoreScheme.m_iScoreMatrixWidth;
	m_iScoreMatrixHeight = tRightScoreScheme.m_iScoreMatrixHeight;
	return *this;
}


/**
@fn void CScoreScheme::SetScore(const int &iFirstElement, const int &iSecondElement, const float &fScore)
*******************************************************************************************************
*/
void CScoreScheme::SetScore(const int &iFirstElement, const int &iSecondElement, const float &fScore)
{
	TRACEI("void CScoreScheme::SetScore(const int &iFirstElement, const int &iSecondElement, const float &fScore)");
	m_tScoreMatrix[iFirstElement][iSecondElement] = fScore;
}


/**
@fn const float& CScoreScheme::GetScore(const int &iFirstElement, const int &iSecondElement) const
为什么 int 传进来的时候也用 & 指针？
************************************************************************************************
*/
const float& CScoreScheme::GetScore(const int &iFirstElement, const int &iSecondElement) const
{
	TRACEI("const int& CScoreScheme::GetScore(const int &iFirstElement, const int &iSecondElement) const");
	return m_tScoreMatrix[iFirstElement][iSecondElement];
}


/**
@fn std::string CScoreScheme::ToString() const
**********************************************
Returns a string that contains score scheme specifications.

@return std::string:
 the score scheme details.
*/
std::string CScoreScheme::ToString() const
{
	TRACEI("std::string CScoreScheme::ToString() const");
	std::ostringstream tOutputStream;

	for (int j = 0; j < m_iScoreMatrixHeight; ++j)
	{
		for (int i=0; i < m_iScoreMatrixWidth; ++i)
		{
			tOutputStream << m_tScoreMatrix[i][j] << "\t";
		}
		tOutputStream << "\n";
	}
	return tOutputStream.str();
}

}; // namespace comp
}; // namespace hishape
