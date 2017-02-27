#pragma once

#undef _DEBUG
#include <python.h>
#pragma comment( lib, "python36.lib" )

#include <stdlib.h>
#include "usefull.h"
//#include <boost/numeric/mtl/mtl.hpp>
//#include <boost/numeric/itl/itl.hpp>

#include <cmath>



#define _DEBUG

#define EL_I 3
#define EL_U 4
#define EL_L 5
#define EL_C 6

//using namespace mtl;
//using namespace itl;
using namespace std;

void sympy_init();
string sympy_sim(string funct);
string sympy_lap(string funct);
string sympy_alap(string f1);
string sympy_eva(string f1, string f2, string f3);
string sympy_dif(string f1, string f2);
string sympy_sep(string f1);
string sympy_roo(string f1);
string sympy_im(string f1);
string sympy_re(string f1);

class NEOPR_COEF;


class POLY
{
public:
	vector<double> c;

	void set_coef(int n, double a);
	string str();
	friend POLY operator *(const POLY& a, const POLY& b);
};


class DROB
{
public:
	int t;//1 - A/x+a   2 - Ax+B/x2+ax+b
	int ste;
	double a, b;
	double A, B;
	int is_max_st;

};

class NEOPR_COEF
{
public:
	POLY p1, p2;
	vector<DROB> dr;
	NEOPR_COEF(POLY _p1, POLY _p2);

	string calc();
};

POLY get_pol_coef(string pol);
string crazy_alap(string funct);
vector<complex<double>> parse_roots(string raw);