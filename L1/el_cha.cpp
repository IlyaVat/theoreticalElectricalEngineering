#pragma once
#undef _DEBUG

#include <stdlib.h>
//#include <boost/numeric/mtl/mtl.hpp>
//#include <boost/numeric/itl/itl.hpp>

#include <cmath>

#include "el_cha.h"
#include "py_sympy.h"
#include "usefull.h"

#define _DEBUG

#define EL_K 0
#define EL_H 1
#define EL_R 2
#define EL_I 3
#define EL_U 4
#define EL_L 5
#define EL_C 6

//using namespace mtl;
//using namespace itl;
using namespace std;


typedef MATRIX<L_S>            MatrixL;
typedef vector<L_S>   VecL;


std::string myreplace(std::string &s,
	const std::string &toReplace,
	const std::string &replaceWith)
{
	size_t fi;
	do
	{
		fi = s.find(toReplace);
		if (fi != -1)
			s.replace(fi, toReplace.length(), "_78_5_f##@@");
	} while (fi != -1);
	do
	{
		fi = s.find("_78_5_f##@@");
		if (fi != -1)
			s.replace(fi, string("_78_5_f##@@").length(), replaceWith);
	} while (fi != -1);
	return s;
}




vector < string > split_command1(string str)
{
	vector < string > res;
	for (int i = 0; i < str.length(); i++)
	{
		if (i == 0 && str[i] == ' ')
		{
			str.erase(i, 1);
			i--;
		}
		if (i>0 && i == str.length() - 1 && str[i] == ' ')
		{
			str.erase(i, 1);
			i -= 2;
		}
		if (i > 0 && str[i] == ' ' && str[i - 1] == ' ')
		{
			str.erase(i, 1);
			i--;
		}
	}

	str += ' ';

	string resres;
	for (int i = 0; i < str.length(); i++)
	{
		if (str[i] == ' ')
		{
			res.push_back(resres);
			resres = "";
		}
		else
		{
			resres += str[i];
		}
	}


	return res;
}

ostream& operator<<(ostream& os, const  L_S& dt)
{
	os << dt.s;
	return os;
}

int get_r_uzl(vector<UZL_P_L> vec, int uzl)
{
	bool b;
	do
	{
		b = 1;
		for (int i = 0; i < vec.size(); i++)
			if (vec[i].id_uz == uzl)
			{
				uzl = vec[i].id_ss;
				b = 0;
				if (vec[i].id_ss == vec[i].id_uz)
					b = 1;
			}

	} while (b == 0);
	return uzl;
}



EL_CHAIN_L::EL_CHAIN_L()
{

}

EL_CHAIN_L::EL_CHAIN_L(string _data)
{
	auto raw = split_command1(_data);

	po.resize(100);
	el.resize(100);

	int max_po = 0;
	int max_el = 0;


	for (int i = 0; i < raw.size(); i += 5)
	{
		int n_e = atoi(raw[i + 0].data());
		int n_p1 = atoi(raw[i + 1].data());
		int n_p2 = atoi(raw[i + 2].data());
		char cht = raw[i + 3][0];
		L_S val = raw[i + 4].data();
		if (cht == 'R')
			el[n_e].t = 2;
		if (cht == 'U')
			el[n_e].t = 4;
		if (cht == 'I')
			el[n_e].t = 3;
		if (cht == 'L')
			el[n_e].t = 5;
		if (cht == 'C')
			el[n_e].t = 6;
		if (cht == 'K')
			el[n_e].t = 0;

		el[n_e].p1 = n_p1;
		el[n_e].p2 = n_p2;
		el[n_e].zn = val;

		max_po = max(max(max_po, n_p1), n_p2);
		max_el = max(max_el, n_e);
	}

	po.resize(max_po + 1);
	el.resize(max_el + 1);


}



void EL_CHAIN_L::comp_h1(int id_res)
{
	comp_par_1_iu_uns();
	h1_l = sympy_sim("(" + el[id_res].I.s + ")/s");
}

void EL_CHAIN_L::comp_par_1_iu_uns()
{

	L_S a,b,c;
	a = (string)"100";
	b = (string)"s+1/s";
	c = a + b;
	c += c;

	vector<int> id_uzl;
	for (int i = 0; i < el.size(); i++)
	{
		if (el[i].p1 != el[i].p2)
		{
			int b = 1;
			for (int r = 0; r < id_uzl.size(); r++)
			{
				if (el[i].p1 == id_uzl[r])
					b = 0;
			}
			if (b)
			{
				id_uzl.push_back(el[i].p1);
			}

			b = 1;
			for (int r = 0; r < id_uzl.size(); r++)
			{
				if (el[i].p2 == id_uzl[r])
					b = 0;
			}
			if (b)
			{
				id_uzl.push_back(el[i].p2);
			}
		}

	}

	//id_uzl-массив номеров существующих узлов
	//последний узел имеет потенциал 0 а с остальными мы и будем работать


	vector<UZL_P_L> prop_uzl;

	for (int i = 0; i < el.size(); i++)
	{
		if (el[i].p1 != el[i].p2)
		{
			if (el[i].t == 0)
			{
				int p1 = get_r_uzl(prop_uzl, el[i].p1);
				int p2 = get_r_uzl(prop_uzl, el[i].p2);

				UZL_P_L puz;
				puz.id_uz = p2;
				puz.id_ss = p1;
				puz.val_add = 0;

				prop_uzl.push_back(puz);
				int o = 0;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p2)
						o = r;
				id_uzl.erase(id_uzl.begin() + o);
				o += 2;
			}
		}

	}

	//id_uzl-массив номеров существующих узлов
	//prop_uzl-массив номеров МНИМЫХ узлов
	//get_r_uzl даёт реальный узел

	int baz_uzl;
	int skip_uzl;

	for (int i = 0; i < el.size(); i++)
	{
		if (el[i].t == 3 || el[i].t == 4)
		{
			baz_uzl = get_r_uzl(prop_uzl, el[i].p1);

			UZL_P_L puz;
			puz.id_uz = el[i].p1;
			puz.id_ss = el[i].p1;
			puz.val_add = 0;

			prop_uzl.push_back(puz);
			int o = 0;
			for (int r = 0; r < id_uzl.size(); r++)
				if (id_uzl[r] == el[i].p1)
					o = r;
			id_uzl.erase(id_uzl.begin() + o);

		}

	}

	for (int i = 0; i < el.size(); i++)
	{
		if (el[i].t == 4)
		{
			baz_uzl = get_r_uzl(prop_uzl, el[i].p2);
			skip_uzl = get_r_uzl(prop_uzl, el[i].p1);

			int o = 0;
			for (int r = 0; r < id_uzl.size(); r++)
				if (id_uzl[r] == get_r_uzl(prop_uzl, el[i].p2))
					o = r;
			id_uzl.erase(id_uzl.begin() + o);

		}

	}

	MatrixL A(id_uzl.size(), id_uzl.size());
	VecL B(id_uzl.size());

	for (int i = 0; i < id_uzl.size(); i++)
		for (int r = 0; r < id_uzl.size(); r++)
		{
			A[i][r] = 0;
			B[i] = 0;
		}

	int check_id = -1;
	L_S check_zn;

	for (int i = 0; i < el.size(); i++)
	{
		int p1 = get_r_uzl(prop_uzl, el[i].p1);
		int p2 = get_r_uzl(prop_uzl, el[i].p2);

		if (p1 != p2)
		{
			if (el[i].t == 4)
			{
				int o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p1)
						o = r;

				int p1arr = o;

				o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p2)
						o = r;

				int p2arr = o;

				check_id = p2;
				check_zn = -el[i].zn;






			}
			if (el[i].t == 3)
			{
				int o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p1)
						o = r;

				int p1arr = o;

				o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p2)
						o = r;

				int p2arr = o;

				if (p1arr != -1)
				{
					B[p1arr] -= el[i].zn;
				}
				if (p2arr != -1)
				{
					B[p2arr] += el[i].zn;
					//cout << B;
				}

			}
			if (el[i].t == 2)
			{

				int o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p1)
						o = r;

				int p1arr = o;

				o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == p2)
						o = r;

				int p2arr = o;

				if (p1arr != -1)
				{
					L_S t = 1 / el[i].zn;
					A[p1arr][p1arr] += 1 / el[i].zn;
				}
				if (p2arr != -1)
				{
					A[p2arr][p2arr] += 1 / el[i].zn;
					if (p1arr != -1)
					{
						A[p1arr][p2arr] -= 1 / el[i].zn;
						A[p2arr][p1arr] -= 1 / el[i].zn;
					}
				}

			}

		}
	}

	for (int i = 0; i < el.size(); i++)
	{
		int p1 = get_r_uzl(prop_uzl, el[i].p1);
		int p2 = get_r_uzl(prop_uzl, el[i].p2);

		if (p1 != p2)
		{
			if (p1 == check_id || p2 == check_id)
			{
				if (el[i].t == 2)
				{
					int o = -1;
					for (int r = 0; r < id_uzl.size(); r++)
						if (id_uzl[r] == p1)
							o = r;

					int p1arr = o;

					o = -1;
					for (int r = 0; r < id_uzl.size(); r++)
						if (id_uzl[r] == p2)
							o = r;

					int p2arr = o;

					if (p1 != check_id)
					{
						if (p1arr != -1)
							B[p1arr] += check_zn * 1 / el[i].zn;
					}
					if (p2 != check_id)
					{
						if (p2arr != -1)
							B[p2arr] += check_zn * 1 / el[i].zn;

					}
				}
			}

		}
	}

	cout << "------------------------------" << endl;
	cout << "------------------------------" << endl;
	cout << endl << "Moon matrix" << endl << A;
	//cout << endl << "Moon vector" << endl << B;
	cout << endl << "Moon vert vector" << endl;
	for (int i = 0; i < id_uzl.size(); i++)
		cout << get_r_uzl(prop_uzl, id_uzl[i]) << " ";

	cout << "------------------------------" << endl;
	cout << "------------------------------" << endl;
	auto res = solveAXB(A, B);

	cout << "------------------------------" << endl;
	//cout << endl << "Res vector" << endl << res;

	for (int i = 0; i < id_uzl.size(); i++)
	{
		po[id_uzl[i]].U = res[i];
	}
	if (check_id != -1)
		po[check_id].U = check_zn;
	for (int i = 0; i < prop_uzl.size(); i++)
		if (prop_uzl[i].id_ss == prop_uzl[i].id_uz)
		{
			po[prop_uzl[i].id_ss].U = 0;
			prop_uzl.erase(prop_uzl.begin() + i);
		}

	int reps = prop_uzl.size();

	for (int i = 0; i < reps; i++)
		for (int r = 0; r < prop_uzl.size(); r++)
			if (prop_uzl[r].id_ss != prop_uzl[r].id_uz)
			{
				int b = 1;
				for (int o = 0; o < prop_uzl.size(); o++)
					if (prop_uzl[o].id_uz == prop_uzl[r].id_ss)
						b = 0;
				if (b == 1)
				{
					po[prop_uzl[r].id_uz].U = po[prop_uzl[r].id_ss].U;
					prop_uzl.erase(prop_uzl.begin() + r);
				}
			}


	cout << "------------------------------" << endl;
	//cout << endl << "Uzl potents" << endl << res;

	for (int i = 0; i < po.size(); i++)
		cout << po[i].U << " ";

	//I p1 to p2
	//U p1+  p2-

	for (int i = 0; i < el.size(); i++)
	{
		if (el[i].t == 0)
		{
			el[i].U = 0;
		}
		if (el[i].t == 1)
		{
			el[i].U = po[el[i].p1].U - po[el[i].p2].U;
			el[i].I = 0;
		}
		if (el[i].t == 2)
		{
			el[i].U = po[el[i].p1].U - po[el[i].p2].U;
			el[i].I = el[i].U / el[i].zn;
		}
		if (el[i].t == 3)
		{
			el[i].I = el[i].zn;
			el[i].U = po[el[i].p1].U - po[el[i].p2].U;
		}
		if (el[i].t == 4)
		{
			el[i].U = el[i].zn;
		}
	}

	vector<int> comp(el.size());

	for (int i = 0; i < el.size(); i++)
		comp[i] = 0;

	for (int rep = 0; rep < po.size(); rep++)
		for (int i = 0; i < po.size(); i++)
		{
			L_S i_u = 0;
			int unk_el = -1;
			int count = 0;
			for (int r = 0; r < el.size(); r++)
				if (el[r].p1 != el[r].p2)
					if (el[r].p1 == i || el[r].p2 == i)
					{

						if (comp[r] == 0 && (el[r].t == EL_K || el[r].t == EL_U))
						{
							unk_el = r;
							count++;
						}

						if (el[r].t == EL_R || el[r].t == EL_I || comp[r] == 1 && (el[r].t == EL_K || el[r].t == EL_U))
						{
							if (el[r].p1 == i)
							{
								i_u += el[r].I;
							}
							if (el[r].p2 == i)
							{
								i_u -= el[r].I;
							}
						}

					}

			if (count == 1)
			{
				comp[unk_el] = 1;

				if (el[unk_el].p2 == i)
				{
					el[unk_el].I = i_u;
				}
				if (el[unk_el].p1 == i)
				{
					el[unk_el].I = -i_u;
				}
			}
		}


	cout << "------------------------------" << endl;
//	cout << endl << "Out info" << endl << res;

	for (int i = 0; i < el.size(); i++)
	{
		if (el[i].t == 2)
			cout << "R" << " ";
		if (el[i].t == 0)
			cout << "K" << " ";
		if (el[i].t == 1)
			cout << "H" << " ";
		if (el[i].t == 3)
			cout << "I" << " ";
		if (el[i].t == 4)
			cout << "U" << " ";
		cout << el[i].p1 << " " << el[i].p2 << "   zn:" << el[i].zn << " " << "U:" << el[i].U << " I:" << el[i].I << " " << endl;
	}

	int i = 5;
}
