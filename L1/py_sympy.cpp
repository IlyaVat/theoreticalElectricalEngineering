#undef _DEBUG
#include <python.h>
#pragma comment( lib, "python36.lib" )

#include <stdlib.h>
//#include <boost/numeric/mtl/mtl.hpp>
//#include <boost/numeric/itl/itl.hpp>

#include <cmath>
#include "el_cha.h"
#include "usefull.h"

//using namespace mtl;
//using namespace itl;
using namespace std;


#define _DEBUG



PyObject *myM, *sympyM;
PyObject *pName, *pDict, *pFunc, *pValue, *sympy_simpyfy, *lapl, *py_print;
PyObject *Fhello;
PyObject *Fsim, *Flap, *Feva, *Fdif, *Falap;
PyObject *Fsep;
PyObject *Froo;
PyObject *Fim;
PyObject *Fre;
PyObject *Flim;

PyObject *moD;

string sympy_sep(string f1);
string sympy_roo(string f1);

string sympy_sim(string funct)
{
	myreplace(funct, ",", ".");
	if (PyCallable_Check(Fsim))
	{
		auto res = PyObject_CallObject(Fsim, Py_BuildValue("(s)", funct.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		if (!res)
		{
			cout << "\nsympy sim error at:\n"<<funct;
			system("pause");
		}
		char *s;
		s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}
string sympy_lap(string funct)
{
	myreplace(funct, ",", ".");
	if (PyCallable_Check(Flap))
	{
		auto res = PyObject_CallObject(Flap, Py_BuildValue("(s)", funct.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s;
		s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}



void POLY::set_coef(int n, double a)
{
	c.resize(max((int)c.size(), n + 1));

	c[n] = a;
}
string POLY::str()
{
	string res = "0";

	for (int i = 0; i < c.size(); i++)
	{
		res += "+";
		res += "s**";
		res += ftos(i);
		res += "*(";
		res += ftos(c[i]);
		res += ")";
	}
	return res;
}
POLY operator *(const POLY& a, const POLY& b)
{
	POLY res;
	res.c.resize(a.c.size() + b.c.size() - 1);
	for (int i = 0; i < a.c.size(); i++)
	for (int r = 0; r < b.c.size(); r++)
	{
		res.c[i + r] += a.c[i]*b.c[r];
	}
	return res;
}

POLY get_pol_coef(string pol)
{
	vector<double> res;

	string part;

	POLY pol1;

	pol += " + ";

	double mno = 1;

	for (int i = 0; i < pol.length(); i++)
	{
		if (pol[i] == ' ')
		{

			string pr[3];
			string part1;
			part += "*";
			for (int r = 0, tt = 0; r < part.length(); r++)
			{
				if (part[r] == '*')
				{
					if (part1 == "s" && tt == 0)
					{
						pr[0] = "1";
						tt = 1;
					}
					if (part1 == "-s" && tt == 0)
					{
						pr[0] = "-1";
						tt = 1;
					}
					pr[tt] = part1;



					tt = 1;
					if (part[r + 1] == '*')
					{
						r++;
						tt = 2;
					}
					part1 = "";
				}
				else
					part1 += part[r];
			}


			if (pr[1] == "")
			{
				pol1.set_coef(0, atof1(pr[0].data())*mno);
			}
			else
			{
				int v2 = atof1(pr[2].data());
				v2 = atof1(pr[0].data());
				if (pr[2] == "")
					pol1.set_coef(1, atof1(pr[0].data())*mno);
				else
					pol1.set_coef((int)atof1(pr[2].data()), atof1(pr[0].data())*mno);

			}

			if (i + 1 < pol.length())
			{
				if (pol[i + 1] == '-')
					mno = -1;
				if (pol[i + 1] == '+')
					mno = 1;
			}

			i += 2;

			part = "";
		}
		else
			part += pol[i];
	}
	return pol1;
}

vector<complex<double>> parse_roots(string raw)
{
	vector<complex<double>> res;

	myreplace(raw, " ", "");
	//myreplace(raw, "", "");

	string part;

	for (int i = 1; i < raw.length(); i++)
	{
		if (raw[i] == ',' || raw[i] == ']')
		{
			complex<double> roo;
			string pr[3];
			string part1;
			//part += "*";
			if (part[part.length() - 1] == 'I')
			{

				string inum_s;
				for (int r = part.length() - 3; r >= 0; r--)
				{
					if (part[r] == '-')
					{
						inum_s = '-' + inum_s;
						r = 0;
					}
					else if (part[r] == '+')
					{
						r = 0;
					}
					else
						inum_s = part[r] + inum_s;
				}
				roo._Val[1] = atof1(inum_s.data());
			}

			string num_s;
			num_s += part[0];
			for (int r = 1; r <part.length(); r++)
			{
				if (part[r] == '-' || part[r] == '+' || r ==part.length()-1)
				{
					if (r == part.length() - 1 && !(part[r] == '-' || part[r] == '+'))
						num_s = num_s + part[r];
					roo._Val[0] = atof1(num_s.data());
					r = part.length();
				}
				else
					num_s = num_s + part[r];
			}
			
			res.resize(res.size()+1);
			res[res.size() - 1] = roo;
			part = "";
		}
		else
			part += raw[i];
	}

	return res;
}



NEOPR_COEF::NEOPR_COEF(POLY _p1, POLY _p2)
{
	p1 = _p1;
	p2 = _p2;
	//calc();
}

long double fact(int N)
{
	if (N < 0) // если пользователь ввел отрицательное число
		return 0; // возвращаем ноль
	if (N == 0) // если пользователь ввел ноль,
		return 1; // возвращаем факториал от нул€ - не удивл€етесь, но это 1 =)
	else // ¬о всех остальных случа€х
		return N * fact(N - 1); // делаем рекурсию.
}

string NEOPR_COEF::calc()
{
	auto root = parse_roots(sympy_roo(p2.str()));
	dr.resize(0);
	for (int i = 0; i < root.size(); i++)
	{
		if (root[i].imag()==0 || abs(root[i].imag() / root[i].real())<0.0000001)//условие нулегого комплексного
		{
			int count_dr = 0;
			for (int r = 0; r < dr.size(); r++)
			{
				if (dr[r].t == 1)
				if (
					dr[r].a == root[i].real() || 
					(min(abs(dr[r].a), abs(root[i].real())) != 0) && 
					abs(dr[r].a + root[i].real()) / min(abs(dr[r].a), abs(root[i].real()))<0.000001
					)//условие малого отличи€
				{
					count_dr++;
					dr[r].is_max_st = 0;
				}
			}
			count_dr++;
			dr.resize(dr.size() + 1);

			dr[dr.size() - 1].t = 1;
			dr[dr.size() - 1].is_max_st = 1;
			dr[dr.size() - 1].ste = count_dr;
			dr[dr.size() - 1].a = -root[i].real();

		}
		else
		{
			if (root[i].imag() > 0)//комплексники всегда парные
			{
				double a, b;

				a = -root[i].real() * 2;//из вольфрама:(x+a+bi)*(x+a-bi)
				b = root[i].real()*root[i].real() + root[i].imag()*root[i].imag();


				int count_dr = 0;
				for (int r = 0; r < dr.size(); r++)
				{
					if (dr[r].t == 2)
					if (abs(dr[r].a - a) / min(abs(dr[r].a), abs(a))<0.000001)//условие малого отличи€
					if (abs(dr[r].b - b) / min(abs(dr[r].b), abs(b))<0.000001)//условие малого отличи€
					{
						dr[r].is_max_st = 0;
						count_dr++;
						a = dr[r].a;
						b = dr[r].b;
					}
				}

				count_dr++;
				dr.resize(dr.size() + 1);

				dr[dr.size() - 1].is_max_st = 1;
				dr[dr.size() - 1].t = 2;
				dr[dr.size() - 1].ste = count_dr;
				dr[dr.size() - 1].a = a;
				dr[dr.size() - 1].b = b;

			}
			
		}

		//
	}

	cout << "\ndr:\n";

	for (int i = 0; i < dr.size(); i++)
	{
		cout << dr[i].t << " " << dr[i].a << " "<< dr[i].b<<endl;
	}

	vector <double> coefs;
	vector <POLY> pols;

	Matrix m(root.size(), root.size());

	for (int o = 0; o < root.size(); o++)
		for (int j = 0; j < root.size(); j++)
			m[o][j] = 0;

	//cout << m;

	coefs.resize(root.size());
	pols.resize(root.size());

	for (int i = 0, r = 0; i < dr.size(); i++)
	{
		if (dr[i].t == 1)
		{
			POLY mn;
			mn.set_coef(0, 1);
			for (int o = 0; o < dr.size(); o++)
				if (o != i)
				{
					POLY sub;
					if (dr[o].t == 1)
					{
						sub.set_coef(1, 1);
						sub.set_coef(0, dr[o].a);
					}
					if (dr[o].t == 2)
					{
						sub.set_coef(2, 1);
						sub.set_coef(1, dr[o].a);
						sub.set_coef(0, dr[o].b);
					}
					if (dr[o].is_max_st == 1)
					{
						int is_sam = 0;
						if (dr[o].t == dr[i].t && dr[o].a == dr[i].a && dr[o].b == dr[i].b)
						{
							is_sam = 1;
						}



						if (is_sam==0)
						for (int j = 0; j < dr[o].ste; j++)
							mn = mn * sub;
						if (is_sam != 0)
						{
							for (int j = 0; j < dr[o].ste - dr[i].ste; j++)
								mn = mn * sub;

						}
					}
				}

			pols[r] = mn;

			for (int o = 0; o < pols[r].c.size(); o++)
				m[o][r] = pols[r].c[o];

			r++;
		}
		if (dr[i].t == 2)
		{
			if (dr[i].ste > 1)
				return "";
			POLY mn;
			mn.set_coef(0, 1);
			for (int o = 0; o < dr.size(); o++)
				if (o != i)
				{
					POLY sub;
					if (dr[o].t == 1)
					{
						sub.set_coef(1, 1);
						sub.set_coef(0, dr[o].a);
					}
					if (dr[o].t == 2)
					{
						sub.set_coef(2, 1);
						sub.set_coef(1, dr[o].a);
						sub.set_coef(0, dr[o].b);
					}
					for (int j = 0; j < dr[o].ste; j++)
						mn = mn * sub;
				}

			POLY ps;
			ps.set_coef(1, 1);

			pols[r+1] = mn;
			pols[r] = mn*ps;

			//cout << m;
			for (int o = 0; o < pols[r].c.size() && o<root.size(); o++)
				m[o][r] = pols[r].c[o];
			//cout << m;
			for (int o = 0; o < pols[r+1].c.size() && o<root.size(); o++)
				m[o][r + 1] = pols[r + 1].c[o];
			//cout << m;

			r++;
			r++;
		}
			
	}
		
	cout << m;

	Vec a(root.size());

	for (int o = 0; o < root.size(); o++)
	{
		if (o<p1.c.size())
			a[o] = p1.c[o];
		else
			a[o] = 0;

	}
	//cout << a;

	auto nko = solveAXB(m, a);
	//cout << "\nnko:\n" << nko;

	//коэффициенты получены
	//примен€ем таблицу преобразовани€ лапласса
	//
	//n!/s^n = t^(n-1)
	//1=b(t)... вр€тли б будет в f2
	//1/(s+a)=exp(-a*t)
	//n!/(s+a)^(n+1)=t^n*exp(-a*t)
	//b/((s+a)^2+b^2)=exp(-a*t)*sin(b*t)
	//(s-a)/((s+a)^2+b^2)=exp(-a*t)*cos(b*t)

	string full_res="0";

	for (int i = 0, r = 0; i < dr.size(); i++)
	{
		if (dr[i].t == 1)
		{
			dr[i].A = nko[r];

			//A/(s+a)^ste

			///1/(s+a)^(n+1)=t^n*exp(-a*t)/n!

			full_res += "+";
			full_res += "(" + ftos(dr[i].A/fact(dr[i].ste-1)) + ")*t**(" + ftos(dr[i].ste-1) + ")*exp(t*(-(" + ftos(dr[i].a) + ")))";

			if (dr[i].a != 0)
			{
				//full_res += "+";
				//full_res += "(" + ftos(dr[i].A) + ")*exp(t*(-(" + ftos(dr[i].a) + ")))";
			}
			else
			{
				//full_res += "+";
				//full_res += "(" + ftos(dr[i].A) + ")";
			}
			r++;
		}
		if (dr[i].t == 2)
		{
			dr[i].A = nko[r];
			dr[i].B = nko[r+1];

			//(Ax+B)/(s^2+a*s+b)
			//b/((s+a)^2+b^2)=exp(-a*t)*sin(b*t)
			//(s-a)/((s+a)^2+b^2)=exp(-a*t)*cos(b*t)
			//(s+a)^2+b^2=s^2+2a*s+a^2+b^2
			//a=a/2

			double a = dr[i].a / 2;
			double b = sqrt(dr[i].b - a*a);
			double c = dr[i].A;
			double d = dr[i].B - dr[i].A*a;


			//(c(s+a)+d)/((s+a)^2+b^2)
			//d*(1)/((s+a)^2+b^2)
			//c*(s+a)/((s+a)^2+b^2)
			full_res += "+";
			full_res += "(" + ftos(d / b) + ")*exp(t*(-(" + ftos(a) + ")))*sin(" + ftos(b) + "*t)";

			full_res += "+";
			full_res += "(" + ftos(c) + ")*exp(t*(-(" + ftos(a) + ")))*cos(" + ftos(b) + "*t)";


			r++;
			r++;
		}

	}

	cout << "\ndr:\n";

	for (int i = 0; i < dr.size(); i++)
	{
		cout << dr[i].t << " " << dr[i].a << " " << dr[i].b << " " << dr[i].A << " " << dr[i].B << endl;
	}

	cout << "\nRES:\n" << full_res;

	return "("+full_res+")/("+ftos(p2.c[p2.c.size()-1])+")";
}


string crazy_alap(string f)
{
	string res;
	string po1, po2;
	
	for (int i = 0, tt = 0; i < f.length(); i++)
	{
		if (tt == 0)
		{
			if (f[i] == '/')
			{
				tt = 1;
			}
			else
				po1 += f[i];
		}
		else
		{
			po2 += f[i];
		}
	}



	po1 = sympy_sep(po1);
	po2 = sympy_sep(po2);


	cout << "\np1:\n" << po1;
	cout << "\np2:\n" << po2;

	POLY p1, p2;
	p1 = get_pol_coef(po1);
	p2 = get_pol_coef(po2);


	cout << "\np1:\n" << p1.str();
	cout << "\np2:\n" << p2.str();


	auto root=parse_roots(sympy_roo(po2));

	cout << "\nroots:\n";

	//надо разобрать корни по парам
	NEOPR_COEF calc(p1,p2);
	res = calc.calc();

	return "("+res+")*Heaviside(t)";
}


string sympy_alap(string funct)
{
	myreplace(funct, ",", ".");
	if (PyCallable_Check(Falap))
	{
		auto res1 = PyObject_CallObject(Falap, Py_BuildValue("(s)", funct.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s;
		s = PyUnicode_AsUTF8(res1);
		return s;
	}
	else
	{
		PyErr_Print();
	}
	return "";
}
string sympy_eva(string f1, string f2, string f3)
{
	myreplace(f1, ",", ".");
	myreplace(f2, ",", ".");
	myreplace(f3, ",", ".");
	if (PyCallable_Check(Feva))
	{
		auto res = PyObject_CallObject(Feva, Py_BuildValue("(s,s,s)", f1.data(), f2.data(), f3.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		string s="-9999999";
		if (!res)
		{
			cout << "\nsympy eva error at:\n" << f1 << "  ||  " << f2 << "  ||  " << f3;
			system("pause");
		}
		if (res)
		{
			s = PyUnicode_AsUTF8(res);
		}
		return s;
	}
	else
	{
		PyErr_Print();
	}
}
string sympy_dif(string f1, string f2)
{
	myreplace(f1, ",", ".");
	myreplace(f2, ",", ".");

	if (PyCallable_Check(Fdif))
	{
		auto res = PyObject_CallObject(Fdif, Py_BuildValue("(s,s)", f1.data(), f2.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s;
		s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}

string sympy_sep(string f1)
{
	myreplace(f1, ",", ".");

	if (PyCallable_Check(Fsep))
	{
		auto res = PyObject_CallObject(Fsep, Py_BuildValue("(s)", f1.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s = "";
		if (!res)
		{
			cout << "\nsympy sep error at:\n" << f1;
			system("pause");
		}
		if (res != nullptr)
			s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}

string sympy_roo(string f1)
{
	int co = 0;
	for (int i = 0; i < f1.size(); i++)
	{
		if (f1[i] == 's')
			co++;
	}
	if (co == 0)
		return "";

	if (PyCallable_Check(Froo))
	{
		auto res = PyObject_CallObject(Froo, Py_BuildValue("(s)", f1.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s = "";
		if (res != nullptr)
			s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}

string sympy_im(string f1)
{
	myreplace(f1, ",", ".");

	if (PyCallable_Check(Fim))
	{
		auto res = PyObject_CallObject(Fim, Py_BuildValue("(s)", f1.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s = "";
		if (res != nullptr)
			s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}

string sympy_re(string f1)
{
	myreplace(f1, ",", ".");

	if (PyCallable_Check(Fre))
	{
		auto res = PyObject_CallObject(Fre, Py_BuildValue("(s)", f1.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s = "";
		if (res != nullptr)
			s = PyUnicode_AsUTF8(res);
		return s;
	}
	else
	{
		PyErr_Print();
	}
}

string sympy_lim(string f1, string f2, string f3, string f4)
{
	myreplace(f1, ",", ".");

	if (PyCallable_Check(Flim))
	{
		auto res = PyObject_CallObject(Flim, Py_BuildValue("(s,s,s,s)", f1.data(), f2.data(), f3.data(), f4.data()));
		//PyObject_CallObject(pFunc, Py_BuildValue("(O)", res));
		char *s = "";
		if (res != nullptr)
			s = PyUnicode_AsUTF8(res);
		string ss = s;
		myreplace(ss, "inf", "oo");
		return ss;
	}
	else
	{
		PyErr_Print();
	}
}

void sympy_init()
{
	Py_Initialize();

	PyRun_SimpleString("import os, sys");
	PyRun_SimpleString("sys.path.append('.')");
	myM = PyImport_ImportModule("mo");
	sympyM = PyImport_ImportModule("sympy");

	pDict = PyModule_GetDict(sympyM);
	moD = PyModule_GetDict(myM);
	// pFunc is also a borrowed reference
	pFunc = PyDict_GetItemString(pDict, "pprint");
	sympy_simpyfy = PyDict_GetItemString(pDict, "simplify");
	lapl = PyDict_GetItemString(pDict, "laplace_transform");
	//py_print = PyDict_GetItemString(NULL, "print");
	Fsim = PyDict_GetItemString(moD, "sim");
	Flap = PyDict_GetItemString(moD, "lap");
	Falap = PyDict_GetItemString(moD, "alap");
	Feva = PyDict_GetItemString(moD, "eva");
	Fdif = PyDict_GetItemString(moD, "dif");
	Fsep = PyDict_GetItemString(moD, "sep");
	Froo = PyDict_GetItemString(moD, "roo");
	Fim = PyDict_GetItemString(moD, "im");
	Fre = PyDict_GetItemString(moD, "re");
	Flim = PyDict_GetItemString(moD, "lim");

}











