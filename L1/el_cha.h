#pragma once
#undef _DEBUG

#include <stdlib.h>
//#include <boost/numeric/mtl/mtl.hpp>
//#include <boost/numeric/itl/itl.hpp>

#include <cmath>

#include "el_cha.h"
#include "usefull.h"
#include "py_sympy.h"

#define _DEBUG

#define EL_I 3
#define EL_U 4
#define EL_L 5
#define EL_C 6

//using namespace mtl;
//using namespace itl;
using namespace std;





std::string myreplace(std::string &s,
	const std::string &toReplace,
	const std::string &replaceWith);
class L_S
{
public:
	string s;
	L_S()
	{

	}
	L_S(double val)
	{
		s = to_string(val);
	}
	L_S(string val)
	{
		s = val;
	}
	friend L_S operator +(L_S a, L_S b)
	{
		L_S res;
		res.s = sympy_sim("(" + a.s + ")+(" + b.s + ")");
		return res;
	}
	friend L_S operator -(L_S a, L_S b)
	{
		L_S res;
		res.s = sympy_sim("(" + a.s + ")-(" + b.s + ")");
		return res;
	}
	friend L_S operator -(L_S a)
	{
		L_S res;
		res.s = sympy_sim("-(" + a.s + ")");
		return res;
	}
	friend L_S operator -=(L_S &a, L_S b)
	{
		a.s = sympy_sim("(" + a.s + ")-(" + b.s + ")");
		return a;
	}
	friend L_S operator +=(L_S &a, L_S b)
	{
		a.s = sympy_sim("(" + a.s + ")+(" + b.s + ")");
		return a;
	}
	friend bool operator ==(L_S &a, double zz)
	{
		if (zz == 0 && a.s == "0")
			return 1;
		return 0;
	}
	friend L_S operator /(L_S a, L_S b)
	{
		L_S res;
		res.s = sympy_sim("(" + a.s + ")/(" + b.s + ")");
		return res;
	}
	friend L_S operator *(L_S a, L_S b)
	{
		L_S res;
		res.s = sympy_sim("(" + a.s + ")*(" + b.s + ")");
		return res;
	}
	friend ostream& operator<<(ostream& os, const  L_S& dt);
};



class ELEM_L
{
public:
	int t;//0-connect 1-no connect 2-R 3- ->; 4- +-; 5-L 6-C
	int p1, p2;
	L_S zn, U, I;
	ELEM_L()
	{
		t = 0;
		zn = 0;
		U = 0;
		I = 0;
		p1 = 0;
		p2 = 0;
	}
};

class POI_L
{
public:
	L_S U;
	POI_L()
	{
		U = 0;
	}
};

class UZL_P_L
{
public:
	int id_uz, id_ss;
	double val_add;
};


class EL_CHAIN_L
{
public:
	vector<POI_L> po;
	vector<ELEM_L> el;

	string h1_l, h1_sympy;
	vector<complex<double>> v_polus, v_zero;

	//UR_SOST ur_so;
	//vector<PER_S> per_s;
	EL_CHAIN_L();
	EL_CHAIN_L(string _data);
	void comp_h1(int id_res);
	void comp_par_1_iu_uns();


};

