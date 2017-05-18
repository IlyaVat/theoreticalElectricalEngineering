#undef _DEBUG
#include <python.h>
#pragma comment( lib, "python36.lib" )

#define _USE_MATH_DEFINES
#include <cmath>
#include "usefull.h"
#include <stdlib.h>

#include "el_cha.h"
#include "py_sympy.h"
#include "samrrr_bibl.h"

#define _DEBUG

#define EL_K 0
#define EL_H 1
#define EL_R 2
#define EL_I 3
#define EL_U 4
#define EL_L 5
#define EL_C 6

#define SIGN_T 5
#define SIGN_A 1
#define SIGN_B 2
#define SIGN_V 3
#define SIGN_G 4
#define SIGN_D 5
#define SIGN_Z 100

//using namespace itl;
using namespace std;

string sympy_sim(string funct);
string sympy_lap(string funct);
string sympy_eva(string f1, string f2, string f3);
string sympy_dif(string f1, string f2);


#define cout std::cout


VecC get_l(Matrix A)
{
	if (A.size_x() == 2)
	{
		VecC rr(2);

		double a, b, c;
		a = 1;
		b = -A[1][1] - A[0][0];
		c = A[1][1] * A[0][0] - A[1][0] * A[0][1];

		double d = b*b - 4 * a*c;

		if (d >= 0)
		{
			rr[0] = (-b + sqrt(d)) / (2 * a);
			rr[1] = (-b - sqrt(d)) / (2 * a);
		}
		if (d < 0)
		{
			rr[0] = (-b + sqrt((complex<double>)d)) / (2 * a);
			rr[1] = (-b - sqrt((complex<double>)d)) / (2 * a);
		}
		return rr;
	}

	cout << "111";
	system("pause");

	return {0};
}


class EL_CHAIN;

void outche(EL_CHAIN& cha);


vector < string > split_command(string str)
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


class UR_SOST
{
public:

	Matrix A;
	Matrix B;
	vector<int> id_cl;//id C and L
	vector<int> id_ui;//id U and I

};


class ELEM
{
public:
	int t;//0-connect 1-no connect 2-R 3- ->; 4- +-; 5-L 6-C
	int p1, p2;
	double zn, U, I;
	ELEM()
	{
		t = 0;
		zn = 0;
		U = 0;
		I = 0;
		p1 = 0;
		p2 = 0;
	}
};

class POI
{
public:
	double U;
	POI()
	{
		U = 0;
	}
};

class UZL_P
{
public:
	int id_uz, id_ss;
	double val_add;
};

class PER_S
{
public:
	vector<double> coef;
	vector<int> type;//0-base 1-sin 2-cos
	vector<int> lev;
	vector<complex<double>> lamd;
	string cool_str;
	string sym_str;

	double inf;
	double f_at(double x)
	{
		double res = 0;
		if (x > 0)
		{
			res += inf;
			int n = lamd.size();
			for (int i = 0; i < n; i++)
			{
				if (type[i] == 0)
				{
					res += exp(x*lamd[i].real())*pow(x, lev[i]);
				}
				if (type[i] == 1)
				{
					res += sin(x*lamd[i].imag())*exp(x*lamd[i].real())*pow(x, lev[i]);
				}
				if (type[i] == 2)
				{
					res += cos(x*lamd[i].imag())*exp(x*lamd[i].real())*pow(x, lev[i]);
				}
			}
		}
		return res;
	}
};

bool in_arr(vector<int> vec, int val)
{
	for (int i = 0; i < vec.size(); i++)
		if (vec[i] == val)
			return true;
	return false;
}

int get_r_uzl(vector<UZL_P> vec, int uzl)
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

class UR_H1_1
{
public:
	vector<double> k_per_s;
	double k_u;
	string sym_str;
};

class UR_H1_2
{
public:
	vector<double> k_per_s;
	double k_u;
	string sym_str;
	vector<complex<double>> v_polus,v_zero;
	L_S str_l;
};


class EL_CHAIN
{
public:
	vector<POI> po;
	vector<ELEM> el;
	

	UR_SOST ur_so;
	vector<PER_S> per_s;
	string obs_vid;
	UR_H1_1 h1_1;
	UR_H1_2 h1_2;
	string f2_s, f2_t;
	string f1_s, f1_t;

	EL_CHAIN()
	{

	}

	EL_CHAIN(string _data)
	{
		auto raw = split_command(_data);

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
			double val = atof1(raw[i + 4].data());
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

	void replace_el(int el1, int elnew, int zn=0)
	{
		for (int i = 0; i < el.size(); i++)
		{
			if (el[i].t == el1)
			{
				el[i].t = elnew;
				el[i].zn = zn;
			}
		}
	}

	void comp_ur_so()
	{
		ur_so.id_cl.resize(0);
		ur_so.id_ui.resize(0);
		for (int i = 0; i < el.size(); i++)
		{
			if (el[i].t == 3 || el[i].t == 4)
				ur_so.id_ui.push_back(i);
			if (el[i].t == 5 || el[i].t == 6)
				ur_so.id_cl.push_back(i);
		}

		ur_so.A.resize(ur_so.id_cl.size(), ur_so.id_cl.size());
		ur_so.B.resize(ur_so.id_cl.size(), ur_so.id_ui.size());

		for (int i = 0; i < ur_so.id_cl.size(); i++)
			for (int r = 0; r < ur_so.id_cl.size(); r++)
				ur_so.A[i][r] = 0;

		for (int i = 0; i < ur_so.id_cl.size(); i++)
		{
			EL_CHAIN podc = *this;


			for (int r = 0; r < el.size(); r++)
			{

				if (podc.el[r].t == 5)
				{
					podc.el[r].t = 3;
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 6)
				{
					podc.el[r].t = 4;
					podc.el[r].zn = 1;
				}
				if (r != ur_so.id_cl[i])
				{
					if (podc.el[r].t == 4)
						podc.el[r].t = 0;

					if (podc.el[r].t == 3)
					{
						podc.el[r].t = 1;
					}
				}

			}

			int id = ur_so.id_cl[i];



			podc.comp_par_1_iu_uns();


			//outche(podc);

			//getch();

			//в цепи подц есть только 1 источник с id_cl[i]

			for (int r = 0; r < ur_so.id_cl.size(); r++)
			{
				int id2 = ur_so.id_cl[r];

				if (el[id2].t == 5)//L
				{
					ur_so.A[r][i] = podc.el[id2].U / el[id2].zn;
				}
				if (el[id2].t == 6)//C
				{
					ur_so.A[r][i] = podc.el[id2].I / el[id2].zn;
				}
			}

			cout << ur_so.A;
		}

		for (int i = 0; i < ur_so.id_ui.size(); i++)
		{
			EL_CHAIN podc = *this;


			for (int r = 0; r < el.size(); r++)
			{

				if (podc.el[r].t == 5)
				{
					podc.el[r].t = 3;
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 6)
				{
					podc.el[r].t = 4;
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 3)
				{
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 4)
				{
					podc.el[r].zn = 1;
				}
				if (r != ur_so.id_ui[i])
				{
					if (podc.el[r].t == 4)
						podc.el[r].t = 0;

					if (podc.el[r].t == 3)
					{
						podc.el[r].t = 1;
					}
				}

			}

			int id = ur_so.id_ui[i];

			podc.comp_par_1_iu_uns();


			//в цепи подц есть только 1 источник с id_cl[i]

			for (int r = 0; r < ur_so.id_cl.size(); r++)
			{
				int id2 = ur_so.id_cl[r];
				//int id = ur_so.id_cl[i];//unused

				if (el[id2].t == 5)//I
				{
					double f = podc.el[id2].U;
					ur_so.B[r][i] = podc.el[id2].U / el[id2].zn;
				}
				if (el[id2].t == 6)//U
				{
					double f = podc.el[id2].I;
					ur_so.B[r][i] = podc.el[id2].I / el[id2].zn;
				}
			}

			cout << ur_so.B;
		}

		cout << ur_so.A;

		auto l_s = get_l(ur_so.A);

		//cout << l_s;
		cout << ur_so.B;

		cout << endl << "ur_sost A:" << ur_so.A;
		cout << endl << "ur_sost B:" << ur_so.B;





	}

	string comp_h1_l(int id_res)
	{
		EL_CHAIN_L cha1;
		cha1.el.resize(el.size());
		cha1.po.resize(po.size());

		for (int i = 0; i < el.size(); i++)
		{
			cha1.el[i].p1 = el[i].p1;
			cha1.el[i].p2 = el[i].p2;
			cha1.el[i].t = el[i].t;
			cha1.el[i].zn = el[i].zn;
			if (el[i].t == EL_U)
			{
				cha1.el[i].zn = 1;
			}
			if (el[i].t == EL_I)
			{
				cha1.el[i].zn = 1;
			}
			if (el[i].t == EL_L)
			{
				cha1.el[i].t = EL_R;
				cha1.el[i].zn = "(" + to_string(el[i].zn) + ")*s";
			}
			if (el[i].t == EL_C)
			{
				cha1.el[i].t = EL_R;
				cha1.el[i].zn = "1/((" + to_string(el[i].zn) + ")*s)";
			}
		}
		cha1.comp_h1(id_res);

		return cha1.h1_l;

	}

	void comp_h1(int id_res)
	{
		

		//h1 (от ступеньки)

		//A+B=f0
		//Al+bl=f'0

		Vec B;
		Matrix A;
		int n = ur_so.A.size_x();
		
		A.resize(n,n);
		B.resize(n);


		auto l_s = get_l(ur_so.A);

		Vec col_B,tv;
		col_B.resize(ur_so.B.size_x());
		for (int o = 0; o < ur_so.B.size_x(); o++)
			col_B[o] = ur_so.B[o][0];

		tv = col_B;
		for (int i = 0; i < n; i++)
			tv[i] *= -1;

		auto F_inf = solveAXB(ur_so.A, tv);
		for (int i = 0; i < n; i++)
		{
			//F_inf[i] *= el[ur_so.id_cl[i]].zn;
		}
		//cout << F_inf;

		vector<Vec> vec_pe_so_lkoef;
		vec_pe_so_lkoef.resize(ur_so.A.size_x());


		vector<Vec> vec_pr_zero;
		vec_pr_zero.resize(ur_so.A.size_x());
		for (int r = 0; r < vec_pr_zero.size(); r++)
		{
			vec_pr_zero[r].resize(ur_so.id_cl.size());
			//как хорошо что есть закон коммутации
			if (r == 0)
			{
				for (int o = 0; o < vec_pr_zero.size(); o++)
				{
					vec_pr_zero[r][o] = 0;
				}
			}
			if (r == 1)
			{
				vec_pr_zero[r] = col_B;
			}
			if (r >= 2)
			{
				//vec_pr_zero[r] = ur_so.A*vec_pr_zero[r-1];//ну наверно это работает но никакой информации в интернете я не нашёл
			}
			

		}

		//l_s[0]._Val[0] = 10;
		//l_s[0]._Val[1] = 15;
		//l_s[1]._Val[0] = 10;
		//l_s[1]._Val[1] = -15;

		vector<int>equal_lev(n);
		vector<int>ignoer(n);
		int is_need_symp = 0;

		//надо бы сначала комплексники по парам разбить а то в цепи 4 порядка имеется небольшой шанс багульки
		for (int i = 0; i < n; i++)
		{
			if (abs(l_s[i].imag()) != 0)
				is_need_symp = 1;

			for (int r = 0; r < i; r++)
			{
				if (ignoer[r] == 0)
					if (abs(l_s[r] - l_s[i])<abs(l_s[r])/100000.0)
					{
						ignoer[r] = 1;
						equal_lev[i] = equal_lev[r] + 1;
						is_need_symp = 1;
					}
			}
		}

		vector<string> l_arr(n);
		vector<string> l_arr_b(n);

		//if(is_need_symp == 1)
		{
			obs_vid = "f(inf)";
			for (int i = 0; i < n; i++)
			{
				string s;
				s = s + "t**(" + to_string(equal_lev[i]) + ")*";

				if (abs(l_s[i].imag()) == 0)
				{
					s = s + "exp(t*(" + to_string(l_s[i].real()) + "))";
				}
				if (abs(l_s[i].imag()) != 0)
				{
					s = s + "exp(t*(" + to_string(l_s[i].real()) + "))";
					s = s + "*" + (l_s[i].imag()>0 ? "sin" : "cos") + "(t*(" + to_string(abs(l_s[i].imag())) + "))";
				}
				cout << "L: " << s << endl;
				l_arr[i] = s;
				l_arr_b[i] = s;
				obs_vid += "+A"+ftos(i+1)+"*";
				obs_vid += s;

			}
		}


		per_s.resize(ur_so.id_cl.size());

		for (int num_cl = 0; num_cl < ur_so.id_cl.size(); num_cl++)
		{
			l_arr = l_arr_b;
			//A-суб лямда матрица
			//B-вектор производных в нуле
			for (int i = 0; i < n; i++)
			{
				B[i] = vec_pr_zero[i][num_cl];
				cout << B[i];
				if (i == 0)
				{
					B[i] -= F_inf[num_cl];
				}
			}
			for (int i = 0; i < n; i++)
				for (int r = 0; r < n; r++)
				{
					if (is_need_symp == 0)
						A[i][r] = pow(l_s[r].real(), i);
					else
					{
						A[i][r] = atof1(sympy_eva(l_arr[r],"t","0").data());
						cout << "\n larr=" << l_arr[r];
						cout << "\n larr(0)=" << A[i][r];
						l_arr[r] = sympy_dif(l_arr[r], "t");
						cout << "\n larr(diff)=" << l_arr[r];
					}
				}

			cout << "\nA\n" << A;
			//cout << "\nB\n" << B;

			auto res = solveAXB(A,B);

			//cout << "\nx\n" << res;

			if (is_need_symp == 1)
			{
				for (int i = 0; i < n; i++)
					cout << "\n" << l_arr_b[i];
			}

			string sres;

			if (el[ur_so.id_cl[num_cl]].t == EL_L)
			{
				sres+= "L";
			}
			if (el[ur_so.id_cl[num_cl]].t == EL_C)
			{
				sres += "C";
			}
			sres += to_string(ur_so.id_cl[num_cl]);
			sres += "(t)=(";

			string sym_r;

			sym_r += to_string(F_inf[num_cl]) + "+";

			for (int i = 0; i < n; i++)
			{
				if (i > 0)
					sym_r += "+";
				sym_r += "(" + l_arr_b[i] + ")*(" + to_string(res[i]) + ")";

			}

			myreplace(sym_r, "exp", "secretthing");
			sym_r = sympy_sim(sym_r);
			myreplace(sym_r, "secretthing", "exp");

			sres += sym_r;
			sres += ")*b1(t)";

			per_s[num_cl].cool_str = sres;
			per_s[num_cl].sym_str = sympy_sim(sym_r);
			per_s[num_cl].coef.resize(n);
			per_s[num_cl].lamd.resize(n);
			per_s[num_cl].lev.resize(n);
			per_s[num_cl].type.resize(n);
			per_s[num_cl].inf=F_inf[num_cl];
			for (int i = 0; i < n; i++)
			{
				per_s[num_cl].coef[i] = res[i];
				per_s[num_cl].lamd[i] = l_s[i];
				per_s[num_cl].lev[i]=equal_lev[i];
				if (l_s[i].imag() == 0)
					per_s[num_cl].type[i] = 0;
				if (l_s[i].imag() > 0)
					per_s[num_cl].type[i] = 1;
				if (l_s[i].imag() < 0)
					per_s[num_cl].type[i] = 2;
			}

		}

		//теперь переменные состояния посчитаны но ежели попадутся 2 одинаковые лямды....
		for (int num_cl = 0; num_cl < ur_so.id_cl.size(); num_cl++)
		{
			cout << "\n" << per_s[num_cl].cool_str;
			cout << "\n" << per_s[num_cl].sym_str;
		}

		h1_1.k_per_s.resize(ur_so.id_cl.size());

		for (int i = 0; i < ur_so.id_cl.size(); i++)
		{
			EL_CHAIN podc = *this;


			for (int r = 0; r < el.size(); r++)
			{

				if (podc.el[r].t == 5)
				{
					podc.el[r].t = 3;
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 6)
				{
					podc.el[r].t = 4;
					podc.el[r].zn = 1;
				}
				if (r != ur_so.id_cl[i])
				{
					if (podc.el[r].t == 4)
						podc.el[r].t = 0;

					if (podc.el[r].t == 3)
					{
						podc.el[r].t = 1;
					}
				}

			}

			int id = ur_so.id_cl[i];

			podc.comp_par_1_iu_uns();


			//в цепи подц есть только 1 источник с id_cl[i]


			//h1_1.k_per_s[i] = podc.el[ur_so.id_ui[i]].I;
			h1_1.k_per_s[i] = podc.el[id_res].I;

		}

		for (int i = 0; i < ur_so.id_ui.size(); i++)
		{
			EL_CHAIN podc = *this;


			for (int r = 0; r < el.size(); r++)
			{

				if (podc.el[r].t == 5)
				{
					podc.el[r].t = 3;
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 6)
				{
					podc.el[r].t = 4;
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 3)
				{
					podc.el[r].zn = 1;
				}
				if (podc.el[r].t == 4)
				{
					podc.el[r].zn = 1;
				}
				if (r != ur_so.id_ui[i])
				{
					if (podc.el[r].t == 4)
						podc.el[r].t = 0;

					if (podc.el[r].t == 3)
					{
						podc.el[r].t = 1;
					}
				}

			}

			int id = ur_so.id_ui[i];

			podc.comp_par_1_iu_uns();


			//в цепи подц есть только 1 источник с id_cl[i]

			//h1_1.k_per_s[i] = podc.el[ur_so.id_ui[i]].I;
			h1_1.k_u = podc.el[id_res].I;
			if (i > 0)
			{
				cout << "so many U and I only 1 allowed\n";
				system("pause");
			}

		}

		cout << "\nh1:";

		string h_s;
		h_s += to_string(h1_1.k_u) + "+";
		
		for (int i = 0; i < ur_so.id_cl.size(); i++)
		{
			if (i > 0)
				h_s += "+";
			h_s = h_s + "(" + to_string(h1_1.k_per_s[i]) + ")*("+per_s[i].sym_str+")";
		}

		cout << h_s<<'\n';
		cout << sympy_sim(h_s)<< '\n';

		h1_1.sym_str=sympy_sim(h_s);

	}

	void comp_par_1_iu_uns()
	{
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
		/*
		проводимость есть: 1/R
		проводимость хх=0
		проводимость кз-значит ошибка в пропуске узлов
		кз это же ИН....



		*/

		vector<UZL_P> prop_uzl;

		for (int i = 0; i < el.size(); i++)
		{
			if (el[i].p1 != el[i].p2)
			{
				if (el[i].t == 0)
				{
					int p1 = get_r_uzl(prop_uzl, el[i].p1);
					int p2 = get_r_uzl(prop_uzl, el[i].p2);

					UZL_P puz;
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

				UZL_P puz;
				puz.id_uz = get_r_uzl(prop_uzl, el[i].p1);
				puz.id_ss = get_r_uzl(prop_uzl, el[i].p1);
				puz.val_add = 0;

				prop_uzl.push_back(puz);
				int o = -1;
				for (int r = 0; r < id_uzl.size(); r++)
					if (id_uzl[r] == get_r_uzl(prop_uzl, el[i].p1))
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

		Matrix A(id_uzl.size(), id_uzl.size());
		Vec B(id_uzl.size());

		for (int i = 0; i < id_uzl.size(); i++)
			for (int r = 0; r < id_uzl.size(); r++)
			{
				A[i][r] = 0;
				B[i] = 0;
			}

		int check_id = -1;
		double check_zn;

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
		
		Vec res;
		if(id_uzl.size()>0)
			res= solveAXB(A, B);

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
				//prop_uzl.erase(prop_uzl.begin() + i);
				//i--;
			}

		int reps = prop_uzl.size();

		for (int i = 0; i < reps; i++)
			for (int r = 0; r < prop_uzl.size(); r++)
				if (prop_uzl[r].id_ss != prop_uzl[r].id_uz)
				{
					po[prop_uzl[i].id_uz].U = po[get_r_uzl(prop_uzl, prop_uzl[i].id_uz)].U;
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
				double i_u = 0;
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
		//cout << endl << "Out info" << endl << res;

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
			cout << el[i].p1 << " " << el[i].p2 << " " << "U:" << el[i].U << " I:" << el[i].I << " " << endl;
		}

		int i = 5;
	}

};






















struct W_MOUSE
{
	int x, y;
	int mchl, mchr, mdol, mdor;
};

typedef struct													// Create A Structure
{
	GLubyte	*imageData;											// Image Data (Up To 32 Bits)
	GLuint	bpp;												// Image Color Depth In Bits Per Pixel.
	GLuint	x;													// Image Width
	GLuint	y;													// Image Height
	GLuint	texID;												// Texture ID Used To Select A Texture
} TextureImage;													// Structure Name

bool LoadTGA(TextureImage *texture, char *filename, int fff = GL_LINEAR/*GL_NEAREST*/)				// Loads A TGA File Into Memory
{
	GLubyte		TGAheader[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };		// Uncompressed TGA Header
	GLubyte		TGAcompare[12];									// Used To Compare TGA Header
	GLubyte		header[6];										// First 6 Useful Bytes From The Header
	GLuint		bytesPerPixel;									// Holds Number Of Bytes Per Pixel Used In The TGA File
	GLuint		imageSize;										// Used To Store The Image Size When Setting Aside Ram
	GLuint		temp;											// Temporary Variable
	GLuint		type = GL_RGBA;									// Set The Default GL Mode To RBGA (32 BPP)

	glEnable(GL_TEXTURE_2D);

	FILE *file = fopen(filename, "rb");							// Open The TGA File

	if (file == NULL ||											// Does File Even Exist?
		fread(TGAcompare, 1, sizeof(TGAcompare), file) != sizeof(TGAcompare) ||	// Are There 12 Bytes To Read?
		memcmp(TGAheader, TGAcompare, sizeof(TGAheader)) != 0 ||	// Does The Header Match What We Want?
		fread(header, 1, sizeof(header), file) != sizeof(header))				// If So Read Next 6 Header Bytes
	{
		if (file == NULL)										// Did The File Even Exist? *Added Jim Strong*
			return false;										// Return false
		else													// Otherwise
		{
			fclose(file);										// If Anything Failed, Close The File
			return false;										// Return false
		}
	}

	texture->x = header[1] * 256 + header[0];				// Determine The TGA x	(highbyte*256+lowbyte)
	texture->y = header[3] * 256 + header[2];				// Determine The TGA y	(highbyte*256+lowbyte)

	if (texture->x <= 0 ||									// Is The x Less Than Or Equal To Zero
		texture->y <= 0 ||									// Is The y Less Than Or Equal To Zero
		(header[4] != 24 && header[4] != 32))						// Is The TGA 24 or 32 Bit?
	{
		fclose(file);											// If Anything Failed, Close The File
		return false;											// Return false
	}

	texture->bpp = header[4];								// Grab The TGA's Bits Per Pixel (24 or 32)
	bytesPerPixel = texture->bpp / 8;							// Divide By 8 To Get The Bytes Per Pixel
	imageSize = (texture->x)*(texture->y)*bytesPerPixel;	// Calculate The Memory Required For The TGA Data
	//imageSize	=128*128*128;	
	texture->imageData = (GLubyte *)malloc(imageSize);			// Reserve Memory To Hold The TGA Data

	if (texture->imageData == NULL ||								// Does The Storage Memory Exist?
		fread(texture->imageData, 1, imageSize, file) != imageSize)	// Does The Image Size Match The Memory Reserved?
	{
		if (texture->imageData != NULL)							// Was Image Data Loaded
			free(texture->imageData);							// If So, Release The Image Data

		fclose(file);											// Close The File
		return false;											// Return false
	}

	for (GLuint i = 0; i<int(imageSize); i += bytesPerPixel)			// Loop Through The Image Data
	{															// Swaps The 1st And 3rd Bytes ('R'ed and 'B'lue)
		temp = texture->imageData[i];								// Temporarily Store The Value At Image Data 'i'
		texture->imageData[i] = texture->imageData[i + 2];		// Set The 1st Byte To The Value Of The 3rd Byte
		texture->imageData[i + 2] = temp;						// Set The 3rd Byte To The Value In 'temp' (1st Byte Value)
	}

	fclose(file);												// Close The File

	// Build A Texture From The Data
	glGenTextures(1, &texture[0].texID);						// Generate OpenGL texture IDs
	fff = GL_LINEAR;
	glBindTexture(GL_TEXTURE_2D, texture[0].texID);				// Bind Our Texture
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, fff);	// Linear Filtered
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, fff);	// Linear Filtered

	if (texture[0].bpp == 24)										// Was The TGA 24 Bits
	{
		type = GL_RGB;											// If So Set The 'type' To GL_RGB
	}

	glTexImage2D(GL_TEXTURE_2D, 0, type, texture[0].x, texture[0].y, 0, type, GL_UNSIGNED_BYTE, texture[0].imageData);

	return true;												// Texture Building Went Ok, Return true
}

class GL_PTINT_TEXT
{
	GLuint base;
	TextureImage tex;
public:
	~GL_PTINT_TEXT()										// Delete The Font List
	{
		free(tex.imageData);
	}

	void init()									// Build Our Bitmap Font
	{
		LoadTGA(&tex, "font.tga");
	}


	void operator ()(const char *fmt, ...)					// Custom GL "Print" Routine
	{
		char		text[512];								// Holds Our String
		va_list		ap;										// Pointer To List Of Arguments

		if (fmt == NULL)									// If There's No Text
			return;											// Do Nothing

		va_start(ap, fmt);									// Parses The String For Variables
		vsprintf(text, fmt, ap);							// And Converts Symbols To Actual Numbers
		va_end(ap);											// Results Are Stored In Text

		int i;
		char *s;
		for (s = text, i = 0; *s; s++, i++)
		{
			int x, y;
			x = *s % 16;
			y = *s / 16;

			glBindTexture(GL_TEXTURE_2D, tex.texID);
			glBegin(GL_POLYGON);
			glTexCoord2f((0 + x) / 16.0, (0 + y) / 16.0); glVertex3f(0 + i*0.8, 0, 0);
			glTexCoord2f((0 + x) / 16.0, (1 + y) / 16.0); glVertex3f(0 + i*0.8, -1, 0);
			glTexCoord2f((1 + x) / 16.0, (1 + y) / 16.0); glVertex3f(1 + i*0.8, -1, 0);
			glTexCoord2f((1 + x) / 16.0, (0 + y) / 16.0); glVertex3f(1 + i*0.8, 0, 0);
			glEnd();

		}


	}
};

GL_PTINT_TEXT glPrint;


#define LOWORD(f) ((f)&0xffff) 
#define HIWORD(f) (((f)>>16)&0xffff) 

class OPENGL_WINDOW
{
private:

	HDC			hDC = NULL;		// Private GDI Device Context
	HGLRC		hRC = NULL;		// Permanent Rendering Context
	HWND		hWnd = NULL;		// Holds Our Window Handle
	HINSTANCE	hInstance;		// Holds The Instance Of The Application

	bool	keys[256];			// Array Used For The Keyboard Routine
	bool	active = true;		// Window Active Flag Set To true By Default
	bool	fullscreen = true;	// Fullscreen Flag Set To Fullscreen Mode By Default

	bool enabled;

	bool initialise;

	int razokx;
	int razoky;
	float pers_angle;
	W_MOUSE mouse;

	//static LRESULT CALLBACK WinMessage(HWND _window, unsigned int _message, WPARAM _wParam, LPARAM _lParam);	
	//LRESULT CALLBACK WinMessage(HWND _window, unsigned int _message, WPARAM _wParam, LPARAM _lParam);

	static LRESULT CALLBACK WndProcS(HWND _window, unsigned int _message, WPARAM _wParam, LPARAM _lParam)
	{
		OPENGL_WINDOW* application = 0;

		if (_message == WM_NCCREATE)
			application = (OPENGL_WINDOW*)(_lParam);
		else
			application = (OPENGL_WINDOW*)GetWindowLongPtr(_window, GWLP_USERDATA);

		//application = (OPENGL_WINDOW*)((long)application & 0x00000000ffffffff);

		return (application)->WndProc(_window, _message, _wParam, _lParam);
	}
	LRESULT CALLBACK WndProc(HWND	hWnd,			// Handle For This Window
		UINT	uMsg,			// Message For This Window
		WPARAM	wParam,			// Additional Message Information
		LPARAM	lParam)			// Additional Message Information
	{
		if (this)
		switch (uMsg)									// Check For Windows Messages
		{
		case WM_ACTIVATE:							// Watch For Window Activate Message
		{
			
			if (!HIWORD(wParam))					// Check Minimization State
			{
				active = true;						// Program Is Active
			}
			else
			{
				//disable();
				//enable();
				active = false;						// Program Is No Longer Active
			}

			return 0;
		}

		case WM_SYSCOMMAND:							// Intercept System Commands
		{
			switch (wParam)							// Check System Calls
			{
			case SC_SCREENSAVE:					// Screensaver Trying To Start?
			case SC_MONITORPOWER:				// Monitor Trying To Enter Powersave?
				//rly this your problem....
				return 0;
			}
			break;
		}

		case WM_CLOSE:								// Did We Receive A Close Message?
		{
			disable();
			return 0;
		}

		case WM_KEYDOWN:							// Is A Key Being Held Down?
		{
			if (wParam >= 0 && wParam < 256)
			keys[wParam] = true;					// If So, Mark It As true
			if (wParam == 910)
			{
				//disable();
				//enable();
			}
			return 0;

		}
		case WM_ACTIVATEAPP:
		{
			break;
		}

		case WM_KEYUP:								// Has A Key Been Released?
		{
			keys[wParam] = false;					// If So, Mark It As false
			return 0;
		}

		case WM_SIZE:								// Resize The OpenGL Window
		{
			int g = lParam;
			razokx = LOWORD(lParam);
			razoky = HIWORD(lParam);
			ReSizeGLScene(LOWORD(lParam), HIWORD(lParam));  // LoWord=Width, HiWord=Height
			break;
		}

		case WM_MOUSEMOVE:
			mouse.x = LOWORD(lParam);
			mouse.y = HIWORD(lParam);
			break;
		case WM_LBUTTONDOWN:
			mouse.mdol = 1;
			mouse.mchl = 1;
			break;
		case WM_LBUTTONUP:
			mouse.mchl = 0;
			break;
		case WM_RBUTTONDOWN:
			mouse.mdor = 1;
			mouse.mchr = 1;
			break;
		case WM_RBUTTONUP:
			mouse.mchr = 0;
			break;

		}

		// Pass All Unhandled Messages To DefWindowProc
		return DefWindowProc(hWnd, uMsg, wParam, lParam);
	}

	GLvoid ReSizeGLScene(GLsizei width, GLsizei height)		// Resize And Initialize The GL Window
	{
		if (height == 0)										// Prevent A Divide By Zero By
		{
			height = 1;										// Making Height Equal One
		}

		glViewport(0, 0, width, height);						// Reset The Current Viewport

		glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
		glLoadIdentity();									// Reset The Projection Matrix

		// Calculate The Aspect Ratio Of The Window
		gluPerspective(90.0f, (GLfloat)width / (GLfloat)height, 0.1, 10000.0f);

		glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
		glLoadIdentity();									// Reset The Modelview Matrix
	}

	int InitGL(GLvoid)										// All Setup For OpenGL Goes Here
	{
		for (int i = 0; i < 256;i++)
			keys[i] = 0;

		glShadeModel(GL_SMOOTH);							// Enable Smooth Shading
		glClearColor(1.0f, 0.0f, 0.0f, 0.5f);				// Black Background
		glClearDepth(1.0f);									// Depth Buffer Setup
		glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
		glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glAlphaFunc(GL_GREATER, 0.99f);
		glDisable(GL_ALPHA_TEST);

		return true;										// Initialization Went OK
	}

	int DrawGLScene(GLvoid)									// Here's Where We Do All The Drawing
	{

		return true;										// Everything Went OK
	}

	GLvoid KillGLWindow(GLvoid)								// Properly Kill The Window
	{
		if (fullscreen)										// Are We In Fullscreen Mode?
		{
			ChangeDisplaySettings(NULL, 0);					// If So Switch Back To The Desktop
			ShowCursor(true);								// Show Mouse Pointer
		}

		if (hRC)											// Do We Have A Rendering Context?
		{
			if (!wglMakeCurrent(NULL, NULL))					// Are We Able To Release The DC And RC Contexts?
			{
				MessageBox(NULL, L"Release Of DC And RC Failed.", L"SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			}

			if (!wglDeleteContext(hRC))						// Are We Able To Delete The RC?
			{
				MessageBox(NULL, L"Release Rendering Context Failed.", L"SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			}
			hRC = NULL;										// Set RC To NULL
		}

		if (hDC && !ReleaseDC(hWnd, hDC))					// Are We Able To Release The DC
		{
			MessageBox(NULL, L"Release Device Context Failed.", L"SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			hDC = NULL;										// Set DC To NULL
		}

		if (hWnd && !DestroyWindow(hWnd))					// Are We Able To Destroy The Window?
		{
			MessageBox(NULL, L"Could Not Release hWnd.", L"SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			hWnd = NULL;										// Set hWnd To NULL
		}

		if (!UnregisterClass(L"OpenGL", hInstance))			// Are We Able To Unregister Class
		{
			MessageBox(NULL, L"Could Not Unregister Class.", L"SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
			hInstance = NULL;									// Set hInstance To NULL
		}
	}

	/*	This Code Creates Our OpenGL Window.  Parameters Are:					*
	*	title			- Title To Appear At The Top Of The Window				*
	*	width			- Width Of The GL Window Or Fullscreen Mode				*
	*	height			- Height Of The GL Window Or Fullscreen Mode			*
	*	bits			- Number Of Bits To Use For Color (8/16/24/32)			*
	*	fullscreenflag	- Use Fullscreen Mode (true) Or Windowed Mode (false)	*/

	BOOL CreateGLWindow(char* title, int width, int height, int bits, bool fullscreenflag)
	{
		GLuint		PixelFormat;			// Holds The Results After Searching For A Match
		WNDCLASS	wc;						// Windows Class Structure
		DWORD		dwExStyle;				// Window Extended Style
		DWORD		dwStyle;				// Window Style
		RECT		WindowRect;				// Grabs Rectangle Upper Left / Lower Right Values
		WindowRect.left = (long)0;			// Set Left Value To 0
		WindowRect.right = (long)width;		// Set Right Value To Requested Width
		WindowRect.top = (long)0;				// Set Top Value To 0
		WindowRect.bottom = (long)height;		// Set Bottom Value To Requested Height

		fullscreen = fullscreenflag;			// Set The Global Fullscreen Flag

		hInstance = GetModuleHandle(NULL);				// Grab An Instance For Our Window
		wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
		wc.lpfnWndProc = (WNDPROC)WndProcS;					// WndProc Handles Messages
		wc.cbClsExtra = 0;									// No Extra Window Data
		wc.cbWndExtra = 0;									// No Extra Window Data
		wc.hInstance = hInstance;							// Set The Instance
		wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);			// Load The Default Icon
		wc.hCursor = LoadCursor(NULL, IDC_ARROW);			// Load The Arrow Pointer
		wc.hbrBackground = NULL;									// No Background Required For GL
		wc.lpszMenuName = NULL;									// We Don't Want A Menu
		wc.lpszClassName = L"OpenGL";								// Set The Class Name

		if (!RegisterClass(&wc))									// Attempt To Register The Window Class
		{
			MessageBox(NULL, L"Failed To Register The Window Class.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;											// Return false
		}

		if (fullscreen)												// Attempt Fullscreen Mode?
		{
			DEVMODE dmScreenSettings;								// Device Mode
			memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));	// Makes Sure Memory's Cleared
			dmScreenSettings.dmSize = sizeof(dmScreenSettings);		// Size Of The Devmode Structure
			dmScreenSettings.dmPelsWidth = width;				// Selected Screen Width
			dmScreenSettings.dmPelsHeight = height;				// Selected Screen Height
			dmScreenSettings.dmBitsPerPel = bits;					// Selected Bits Per Pixel
			dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;

			// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
			if (ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL)
			{
				// If The Mode Fails, Offer Two Options.  Quit Or Use Windowed Mode.
				if (MessageBox(NULL, L"The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?", L"NeHe GL", MB_YESNO | MB_ICONEXCLAMATION) == IDYES)
				{
					fullscreen = false;		// Windowed Mode Selected.  Fullscreen = false
				}
				else
				{
					// Pop Up A Message Box Letting User Know The Program Is Closing.
					MessageBox(NULL, L"Program Will Now Close.", L"ERROR", MB_OK | MB_ICONSTOP);
					return false;									// Return false
				}
			}
		}

		if (fullscreen)												// Are We Still In Fullscreen Mode?
		{
			dwExStyle = WS_EX_APPWINDOW;								// Window Extended Style
			dwStyle = WS_POPUP;										// Windows Style
			ShowCursor(false);										// Hide Mouse Pointer
		}
		else
		{
			dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
			dwStyle = WS_OVERLAPPEDWINDOW;							// Windows Style
		}

		AdjustWindowRectEx(&WindowRect, dwStyle, false, dwExStyle);		// Adjust Window To true Requested Size

		// Create The Window
		if (!(hWnd = CreateWindowEx(dwExStyle,							// Extended Style For The Window
			L"OpenGL",							// Class Name
			L"title not work",								// Window Title
			dwStyle |							// Defined Window Style
			WS_CLIPSIBLINGS |					// Required Window Style
			WS_CLIPCHILDREN,					// Required Window Style
			0, 0,								// Window Position
			WindowRect.right - WindowRect.left,	// Calculate Window Width
			WindowRect.bottom - WindowRect.top,	// Calculate Window Height
			NULL,								// No Parent Window
			NULL,								// No Menu
			hInstance,							// Instance
			this)))
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Window Creation Error.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		SetWindowLongPtr(hWnd, GWLP_USERDATA, (long long)this);

		static	PIXELFORMATDESCRIPTOR pfd =				// pfd Tells Windows How We Want Things To Be
		{
			sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
			1,											// Version Number
			PFD_DRAW_TO_WINDOW |						// Format Must Support Window
			PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
			PFD_DOUBLEBUFFER,							// Must Support Double Buffering
			PFD_TYPE_RGBA,								// Request An RGBA Format
			bits,										// Select Our Color Depth
			0, 0, 0, 0, 0, 0,							// Color Bits Ignored
			0,											// No Alpha Buffer
			0,											// Shift Bit Ignored
			0,											// No Accumulation Buffer
			0, 0, 0, 0,									// Accumulation Bits Ignored
			16,											// 16Bit Z-Buffer (Depth Buffer)  
			0,											// No Stencil Buffer
			0,											// No Auxiliary Buffer
			PFD_MAIN_PLANE,								// Main Drawing Layer
			0,											// Reserved
			0, 0, 0										// Layer Masks Ignored
		};

		if (!(hDC = GetDC(hWnd)))							// Did We Get A Device Context?
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Can't Create A GL Device Context.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		if (!(PixelFormat = ChoosePixelFormat(hDC, &pfd)))	// Did Windows Find A Matching Pixel Format?
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Can't Find A Suitable PixelFormat.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		if (!SetPixelFormat(hDC, PixelFormat, &pfd))		// Are We Able To Set The Pixel Format?
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Can't Set The PixelFormat.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		if (!(hRC = wglCreateContext(hDC)))				// Are We Able To Get A Rendering Context?
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Can't Create A GL Rendering Context.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		if (!wglMakeCurrent(hDC, hRC))					// Try To Activate The Rendering Context
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Can't Activate The GL Rendering Context.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		ShowWindow(hWnd, SW_SHOW);						// Show The Window
		SetForegroundWindow(hWnd);						// Slightly Higher Priority
		SetFocus(hWnd);									// Sets Keyboard Focus To The Window
		ReSizeGLScene(width, height);					// Set Up Our Perspective GL Screen

		if (!InitGL())									// Initialize Our Newly Created GL Window
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL, L"Initialization Failed.", L"ERROR", MB_OK | MB_ICONEXCLAMATION);
			return false;								// Return false
		}

		return true;									// Success
	}




public:
	OPENGL_WINDOW()
	{
		enabled = 0;
	}
	~OPENGL_WINDOW()
	{
		disable();
	}

	void draw()
	{
		DrawGLScene();
		SwapBuffers(hDC);
	}
	void upd()
	{
		mouse.mdol = 0;
		mouse.mdor = 0;
		MSG msg;
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}

	bool enable()
	{
		if (enabled)
			return 0;

		enabled = 1;

		hDC = NULL;
		hRC = NULL;
		hWnd = NULL;

		active = true;

		fullscreen = false;
		CreateGLWindow("NeHe's OpenGL Framework", 640, 480, 16, fullscreen);
		razokx = 640; razoky = 480;

		//fullscreen = true;
		//CreateGLWindow("NeHe's OpenGL Framework", GetSystemMetrics(SM_CXSCREEN), GetSystemMetrics(SM_CYSCREEN), 16, fullscreen);


		return 1;
	}
	bool is_enabled()
	{
		return enabled;
	}
	bool disable()
	{
		if (!enabled)
			return 0;
		enabled = 0;
		KillGLWindow();
		return 1;
	}
	int get_wx()
	{
		return razokx;
	}
	int get_wy()
	{
		return razoky;
	}
	float get_angle()
	{
		return 60;
	}
	const bool* get_keys()
	{
		return keys;
	}
	W_MOUSE get_mouse()
	{
		return mouse;
	}

};

struct COLOR
{
	float r, g, b, a;
};




void put_icon()
{
	glDisable(GL_DEPTH_TEST);
	static float angle = 1;


	angle += 1;

	glPushMatrix();
	glRotatef(angle, 0, 1, 0);

		glColor3f((rand() % 10) / 10.0, (rand() % 10) / 10.0, (rand() % 10) / 10.0);
		//треугольник
		glBegin(GL_POLYGON);
		glVertex3f(0, 1, 0);
		glVertex3f(-1, -1, 0);
		glVertex3f(1, -1, 0);
		glEnd();

		glColor3f((rand() % 10) / 10.0, (rand() % 10) / 10.0, (rand() % 10) / 10.0);
		for (int i = -90; i < 90; i += 10)
		{
			glBegin(GL_POLYGON);
			glVertex3f(i / 180.0, cos(i / 180.0*3.1415926) / 4, 0);
			glVertex3f(i / 180.0, -cos(i / 180.0 * 3.1415926) / 4, 0);
			glVertex3f((i + 10) / 180.0, -cos((i + 10) / 180.0 * 3.1415926) / 4, 0);
			glVertex3f((i + 10) / 180.0, cos((i + 10) / 180.0 * 3.1415926) / 4, 0);
			glEnd();
		}

		glColor3f((rand() % 10) / 10.0, (rand() % 10) / 10.0, (rand() % 10) / 10.0);
		for (int i = 0; i < 360; i += 10)
		{
			glBegin(GL_POLYGON);
			glVertex3f(sin(i / 180.0*3.1415926) / 4, cos(i / 180.0*3.1415926) / 4, 0);
			glVertex3f(sin((i + 10) / 180.0*3.1415926) / 4, cos((i + 10) / 180.0 * 3.1415926) / 4, 0);
			glVertex3f(0, 0, 0);
			glEnd();
		}

	glPopMatrix();
}



class EASY_TEX
{
private:
	vector<byte> data;
	GLuint texID;
	int init;
	int x, y;
public:
	EASY_TEX()
	{
		init = 0;
		x = 0; y = 0; texID = 0;
	}
	GLuint ID()
	{
		return texID;
	}
	int gx()
	{
		return x;
	}
	int gy()
	{
		return y;
	}
	void resize(int _x, int _y)
	{
		if (_x > 100000)_x = 100000;
		if (_y > 100000)_y = 100000;
		int xx = 1;
		int yy = 1;
		while (xx<_x)
			xx = xx << 1;
		while (yy < _y)
			yy = yy << 1;
		x = xx;
		y = yy;

		data.resize(x * y * 3);

		for (int i = 0; i < data.size(); i++)
		{
			data[i] = ((i/3+i/3/x)%2)<<7;//шахматное поле
		}
	}
	const vector<byte> gdata()
	{
		return data;
	}
	void setpixel(int _x, int _y, byte _r, byte _g, byte _b)
	{
		if (_x < 0 || _y < 0 || _x >= x || _y >= y)
		{
			//cout << "draw error ";
		}
		else
		{
			data[(_x + _y*x) * 3+0] = _r;
			data[(_x + _y*x) * 3+1] = _g;
			data[(_x + _y*x) * 3+2] = _b;
			stringstream ss;
			int a = data[(_x + _y*x) * 3 + 0];
			int b= data[(_x + _y*x) * 3 + 1];
			int c = data[(_x + _y*x) * 3 +2];

		}
	}

	void numbers(int px, int py, string s)
	{
		vector<int> n_d[12][5];
		string ru_d[40][5];

		int i;

		i = 0;

		ru_d[i][0] = " 111 ";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "11111";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "11111";
		ru_d[i][1] = "1    ";
		ru_d[i][2] = "1111 ";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1111 ";

		i++;
		ru_d[i][0] = "1111 ";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "1111 ";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1111 ";

		i++;
		ru_d[i][0] = "11111";
		ru_d[i][1] = "1    ";
		ru_d[i][2] = "1    ";
		ru_d[i][3] = "1    ";
		ru_d[i][4] = "1    ";

		i++;
		ru_d[i][0] = " 111 ";
		ru_d[i][1] = " 1 1 ";
		ru_d[i][2] = " 1 1 ";
		ru_d[i][3] = "11111";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "11111";
		ru_d[i][1] = "1    ";
		ru_d[i][2] = "11111";
		ru_d[i][3] = "1    ";
		ru_d[i][4] = "11111";

		i++;
		ru_d[i][0] = "11111";
		ru_d[i][1] = "1    ";
		ru_d[i][2] = "11111";
		ru_d[i][3] = "1    ";
		ru_d[i][4] = "11111";

		i++;
		ru_d[i][0] = "1 1 1";
		ru_d[i][1] = " 111 ";
		ru_d[i][2] = "  1  ";
		ru_d[i][3] = " 111 ";
		ru_d[i][4] = "1 1 1";

		i++;
		ru_d[i][0] = "1111 ";
		ru_d[i][1] = "    1";
		ru_d[i][2] = " 111 ";
		ru_d[i][3] = "    1";
		ru_d[i][4] = "1111 ";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1  11";
		ru_d[i][2] = "1 1 1";
		ru_d[i][3] = "11  1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1  11";
		ru_d[i][2] = "1 1 1";
		ru_d[i][3] = "11  1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1 11 ";
		ru_d[i][2] = "11   ";
		ru_d[i][3] = "1 11 ";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = " 1111";
		ru_d[i][1] = " 1  1";
		ru_d[i][2] = " 1  1";
		ru_d[i][3] = " 1  1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "11 11";
		ru_d[i][2] = "1 1 1";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "11111";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = " 111 ";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "1   1";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = " 111 ";

		i++;
		ru_d[i][0] = "11111";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "1   1";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "1111 ";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "1111 ";
		ru_d[i][3] = "1    ";
		ru_d[i][4] = "1    ";

		i++;
		ru_d[i][0] = " 1111";
		ru_d[i][1] = "1    ";
		ru_d[i][2] = "1    ";
		ru_d[i][3] = "1    ";
		ru_d[i][4] = " 1111";

		i++;
		ru_d[i][0] = "11111";
		ru_d[i][1] = "  1  ";
		ru_d[i][2] = "  1  ";
		ru_d[i][3] = "  1  ";
		ru_d[i][4] = "  1  ";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = " 1111";
		ru_d[i][3] = "    1";
		ru_d[i][4] = "  11 ";

		i++;
		ru_d[i][0] = "  1  ";
		ru_d[i][1] = " 111 ";
		ru_d[i][2] = "1 1 1";
		ru_d[i][3] = " 111 ";
		ru_d[i][4] = "  1  ";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = " 1 1 ";
		ru_d[i][2] = "  1  ";
		ru_d[i][3] = " 1 1 ";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "1  1 ";
		ru_d[i][1] = "1  1 ";
		ru_d[i][2] = "1  1 ";
		ru_d[i][3] = "11111";
		ru_d[i][4] = "    1";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "11111";
		ru_d[i][3] = "    1";
		ru_d[i][4] = "    1";

		i++;
		ru_d[i][0] = "1 1 1";
		ru_d[i][1] = "1 1 1";
		ru_d[i][2] = "1 1 1";
		ru_d[i][3] = "1 1 1";
		ru_d[i][4] = "11111";

		i++;
		ru_d[i][0] = "1 1 1";
		ru_d[i][1] = "1 1 1";
		ru_d[i][2] = "1 1 1";
		ru_d[i][3] = "11111";
		ru_d[i][4] = "    1";

		i++;
		ru_d[i][0] = "1    ";
		ru_d[i][1] = "1    ";
		ru_d[i][2] = "1111 ";
		ru_d[i][3] = "1   1";
		ru_d[i][4] = "1111 ";

		i++;
		ru_d[i][0] = "1   1";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = "11  1";
		ru_d[i][3] = "1 1 1";
		ru_d[i][4] = "11  1";

		i++;
		ru_d[i][0] = "11   ";
		ru_d[i][1] = " 1   ";
		ru_d[i][2] = " 111 ";
		ru_d[i][3] = " 1  1";
		ru_d[i][4] = " 111 ";

		i++;
		ru_d[i][0] = "1111 ";
		ru_d[i][1] = "    1";
		ru_d[i][2] = "11111";
		ru_d[i][3] = "    1";
		ru_d[i][4] = "1111 ";

		i++;
		ru_d[i][0] = "1  1 ";
		ru_d[i][1] = "1 1 1";
		ru_d[i][2] = "111 1";
		ru_d[i][3] = "1 1 1";
		ru_d[i][4] = "1  1 ";

		i++;
		ru_d[i][0] = " 1111";
		ru_d[i][1] = "1   1";
		ru_d[i][2] = " 1111";
		ru_d[i][3] = " 1  1";
		ru_d[i][4] = "1   1";

		i++;
		ru_d[i][0] = "     ";
		ru_d[i][1] = "     ";
		ru_d[i][2] = "     ";
		ru_d[i][3] = "     ";
		ru_d[i][4] = "     ";


		i = 0;
		n_d[i][0] = { 0, 1, 0 };
		n_d[i][1] = { 1, 0, 1 };
		n_d[i][2] = { 1, 0, 1 };
		n_d[i][3] = { 1, 0, 1 };
		n_d[i][4] = { 0, 1, 0 };

		i++;
		n_d[i][0] = { 0, 1, 0 };
		n_d[i][1] = { 1, 1, 0 };
		n_d[i][2] = { 0, 1, 0 };
		n_d[i][3] = { 0, 1, 0 };
		n_d[i][4] = { 1, 1, 1 };

		i++;
		n_d[i][0] = { 1, 1, 0 };
		n_d[i][1] = { 0, 0, 1 };
		n_d[i][2] = { 0, 0, 1 };
		n_d[i][3] = { 0, 1, 0 };
		n_d[i][4] = { 1, 1, 1 };

		i++;
		n_d[i][0] = { 1, 1, 0 };
		n_d[i][1] = { 0, 0, 1 };
		n_d[i][2] = { 0, 1, 0 };
		n_d[i][3] = { 0, 0, 1 };
		n_d[i][4] = { 1, 1, 0 };

		i++;
		n_d[i][0] = { 1, 0, 1 };
		n_d[i][1] = { 1, 0, 1 };
		n_d[i][2] = { 1, 1, 1 };
		n_d[i][3] = { 0, 0, 1 };
		n_d[i][4] = { 0, 0, 1 };

		i++;
		n_d[i][0] = { 1, 1, 1 };
		n_d[i][1] = { 1, 0, 0 };
		n_d[i][2] = { 1, 1, 0 };
		n_d[i][3] = { 0, 0, 1 };
		n_d[i][4] = { 1, 1, 0 };

		i++;
		n_d[i][0] = { 0, 1, 0 };
		n_d[i][1] = { 1, 0, 0 };
		n_d[i][2] = { 1, 1, 0 };
		n_d[i][3] = { 1, 0, 1 };
		n_d[i][4] = { 0, 1, 0 };

		i++;
		n_d[i][0] = { 1, 1, 1 };
		n_d[i][1] = { 0, 0, 1 };
		n_d[i][2] = { 0, 0, 1 };
		n_d[i][3] = { 0, 1, 0 };
		n_d[i][4] = { 0, 1, 0 };

		i++;
		n_d[i][0] = { 0, 1, 0 };
		n_d[i][1] = { 1, 0, 1 };
		n_d[i][2] = { 0, 1, 0 };
		n_d[i][3] = { 1, 0, 1 };
		n_d[i][4] = { 0, 1, 0 };

		i++;
		n_d[i][0] = { 0, 1, 0 };
		n_d[i][1] = { 1, 0, 1 };
		n_d[i][2] = { 0, 1, 1 };
		n_d[i][3] = { 0, 0, 1 };
		n_d[i][4] = { 0, 1, 0 };

		i++;
		n_d[i][0] = { 0, 0, 0 };
		n_d[i][1] = { 0, 0, 0 };
		n_d[i][2] = { 0, 0, 0 };
		n_d[i][3] = { 0, 0, 0 };
		n_d[i][4] = { 0, 0, 1 };

		i++;
		n_d[i][0] = { 0, 0, 0 };
		n_d[i][1] = { 0, 0, 0 };
		n_d[i][2] = { 1, 1, 1 };
		n_d[i][3] = { 0, 0, 0 };
		n_d[i][4] = { 0, 0, 0 };

		int offs = 0;

		for (int i = 0; i < s.length(); i++)
		{
			int id = s[i]-'0';
			if (id < 0 || id >= 10)
			{
				id = -1;
				if (s[i] == '.' || s[i] == ',')
					id = 10,offs+=2;
				if (s[i] == '-')
					id = 11;
			}
			if (id >= 0 && id <= 11)
			{
				if (id >= 0)
					for (int r = 0; r < 10; r++)
						for (int o = 0; o < 6; o++)
							if (n_d[id][r / 2][o / 2])
							{
								setpixel(px + i * 4 * 2 + o - offs * 2, py - r, 255, 255, 255);
							}
			}
			else
			{
				//ё= -72
				switch (s[i])
				{
				case 'а':id = 0; break;
				case 'б':id = 1; break;
				case 'в':id = 2; break;
				case 'г':id = 3; break;
				case 'д':id = 4; break;
				case 'е':id = 5; break;
				case 'ё':id = 6; break;
				case 'ж':id = 7; break;
				case 'з':id = 8; break;
				case 'и':id = 9; break;
				case 'й':id = 10; break;
				case 'к':id = 11; break;
				case 'л':id = 12; break;
				case 'м':id = 13; break;
				case 'н':id = 14; break;
				case 'о':id = 15; break;
				case 'п':id = 16; break;
				case 'р':id = 17; break;
				case 'с':id = 18; break;
				case 'т':id = 19; break;
				case 'у':id = 20; break;
				case 'ф':id = 21; break;
				case 'х':id = 22; break;
				case 'ц':id = 23; break;
				case 'ч':id = 24; break;
				case 'ш':id = 25; break;
				case 'щ':id = 26; break;
				case 'ь':id = 27; break;
				case 'ы':id = 28; break;
				case 'ъ':id = 29; break;
				case 'э':id = 30; break;
				case 'ю':id = 31; break;
				case 'я':id = 32; break;
				}


				if (id >= 0 && id <= 32)
				{
					for (int r = 0; r < 10; r++)
						for (int o = 0; o < 10; o++)
							if (ru_d[id][r / 2][o / 2]!=' ')
							{
								setpixel(px + i * 4 * 2 + o - offs * 2, py - r, 255, 255, 255);
							}
					if (id == 6)
					{
						for (int r = -4; r < -2; r++)
							for (int o = 2; o < 4; o++)
								setpixel(px + i * 4 * 2 + o - offs * 2, py - r, 255, 255, 255);
						for (int r = -4; r < -2; r++)
							for (int o = 6; o < 8; o++)
								setpixel(px + i * 4 * 2 + o - offs * 2, py - r, 255, 255, 255);
					}
					if (id == 10)
					{
						for (int r = -4; r < -2; r++)
							for (int o = 2; o < 8; o++)
								setpixel(px + i * 4 * 2 + o - offs * 2, py - r, 255, 255, 255);
					}
					offs -= 2;
				}
			}
		}


	}

	void line(double x1, double y1, double x2, double y2, int _r=255, int _g=255, int _b=255)
	{
		if ((y1 - y2) * 0 != 0 || (x1 - x2) * 0 != 0)
			return;

		if (abs(x1 - x2) > abs(y1 - y2))//горизонтальная линия
		{
			double a = (y2-y1)/(x2-x1);
			double b = y1-a*x1;



			//y=ax+b

			if (x1 < 0)
			{
				x1 = 0;
				y1 = b;
			}
			if (x1 > x)
			{
				x1 = x;
				y1 = a*x+b;
			}


			if (y1 < 0 && y2 < 0 || y1 > y && y2 > y)
				return;
			int maxi = abs(x1 - x2)+1;
			for (int i = 0; i <= maxi; i++)
			{
				setpixel(x1*i / maxi + x2*(maxi - i) / maxi, y1*i / maxi + y2*(maxi - i) / maxi, _r, _g, _b);
			}

		}
		else
		{
			double a = (x2 - x1) / (y2 - y1);
			double b = x1 - a*y1;



			//x=ay+b

			if (y1 < 0)
			{
				y1 = 0;
				x1 = b;
			}
			if (y1 > y)
			{
				y1 = y;
				x1 = a*y + b;
			}


			if (x1 < 0 && x2 < 0 || x1 > x && x2 > x)
				return;
			int maxi = abs(y1 - y2)+1;
			for (int i = 0; i <= maxi; i++)
			{
				setpixel(x1*i / maxi + x2*(maxi - i) / maxi, y1*i / maxi + y2*(maxi - i) / maxi, _r,_g,_b);
			}

		}
	}

	void gen_broken_tex()
	{
		init = 0;
		gentex();
	}
	void gentex()
	{
		if (init)
			return;
		init = 1;
		if (texID)
		{
			glDeleteTextures(1, &texID);
		}
		// Build A Texture From The Data
		glGenTextures(1, &texID);						// Generate OpenGL texture IDs
		int fff = GL_LINEAR;
		glBindTexture(GL_TEXTURE_2D, texID);				// Bind Our Texture
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, fff);	// Linear Filtered
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, fff);	// Linear Filtered


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, x, y, 0, GL_RGB, GL_UNSIGNED_BYTE, &data[0]);

	}


};

mutex mtx;

class OUTPUT_IMAGES
{
public:
	vector < EASY_TEX > a;
	void add(EASY_TEX& tt)
	{
		a.push_back(tt);
	}
};

void print_chem(EASY_TEX &tex, const EL_CHAIN &cha);

OUTPUT_IMAGES img;

float pow1(double x, double y)
{
	return exp(y* log(x));
}
double pow10(int y)
{
	if (y >= 0)
	{
		double b = 1;
		for (int i = 0; i < y; i++)
			b=b*10;
		return b;
	}
	else
	{
		return 1 / pow10(-y);
	}
}

EASY_TEX create_double_plot(int tx, int ty, vector<double> vx, vector<double> vy, vector<double> vz)
{
	double min_x, max_x, min_y, max_y;
	max_x = vx[0];
	min_x = vx[0];
	max_y = vy[0];
	min_y = vy[0];
	for (int i = 0; i < vx.size(); i++)
	{
		if (vx[i] < min_x)
			min_x = vx[i];
		if (vx[i] > max_x)
			max_x = vx[i];

		if (vy[i] < min_y)
			min_y = vy[i];
		if (vy[i] > max_y)
			max_y = vy[i];
		if (vz[i] < min_y)
			min_y = vz[i];
		if (vz[i] > max_y)
			max_y = vz[i];
	}

	min_x -= (max_x - min_x)*0.1;
	min_y -= (max_y - min_y)*0.1;
	max_x += (max_x - min_x) * 2 / tx;
	max_y += (max_y - min_y) * 2 / ty;

	EASY_TEX tex;
	tex.resize(tx, ty);

	for (int i = 0; i < tx; i++)
		for (int r = 0; r < ty; r++)
			tex.setpixel(i, r, 0, 0, 0);


	//теперь выведем циферки и оси

	double sx, sy;
	int porx, pory;
	sx = abs(max_x - min_x);
	sy = abs(max_y - min_y);

	porx = 10;//10^10
	while (sx / pow1((double)10, porx) < 3 && porx < 100)
		porx--;
	pory = 10;//10^10
	while (sy / pow1((double)10, pory) < 3 && pory < 100)
		pory--;



	float mnox = pow10(porx);
	float mnoy = pow10(pory);

	if (mnox != 0)
		for (int i = min_x / mnox; i < max_x / mnox; i++)
		{
			tex.line((i * mnox - min_x) / sx * tx, 40, (i * mnox - min_x) / sx * tx, ty, 50, 50, 50);

			if ((i * mnox - min_x) / sx * tx>20 && (i * mnox - min_x) / sx * tx < tx - 30)
				tex.numbers((i * mnox - min_x) / sx * tx - 3, 20, ftos(i*mnox));
		}
	if (mnoy != 0)
		for (int i = min_y / mnoy; i < max_y / mnoy; i++)
		{
			tex.line(20, (i * mnoy - min_y) / sy * ty, tx, (i * mnoy - min_y) / sy * ty, 50, 50, 50);

			if ((i * mnoy - min_y) / sy * ty>20 && (i * mnoy - min_y) / sy * ty < ty - 30)
				tex.numbers(4, (i * mnoy - min_y) / sy * ty + 5, ftos(i*mnoy));
		}
	//tex.numbers(30, 20, "111111");


	for (int i = 1; i <vx.size(); i++)
	{
		double rx[2], ry[2];
		double px[2], py[2];

		for (int r = 0; r < 2; r++)
		{
			rx[r] = vx[i - 1 + r];
			ry[r] = vy[i - 1 + r];
			px[r] = (rx[r] - min_x) / (max_x - min_x)*tx;
			py[r] = (ry[r] - min_y) / (max_y - min_y)*ty;
		}
		tex.line(px[0], py[0], px[1], py[1]);

		for (int r = 0; r < 2; r++)
		{
			rx[r] = vx[i - 1 + r];
			ry[r] = vz[i - 1 + r];
			px[r] = (rx[r] - min_x) / (max_x - min_x)*tx;
			py[r] = (ry[r] - min_y) / (max_y - min_y)*ty;
		}
		tex.line(px[0], py[0], px[1], py[1]);
		//tex.numbers(px[0], py[0],ftos(py[0]));
	}

	return tex;

}

EASY_TEX create_plot(int tx, int ty, vector<double> vx, vector<double> vy, double min_x, double max_x, double min_y, double max_y)
{


	if (min_x == min_y)
	{
		max_x = vx[0];
		min_x = vx[0];
		max_y = vy[0];
		min_y = vy[0];
		for (int i = 0; i < vx.size(); i++)
		{
			if (vx[i] < min_x)
				min_x = vx[i];
			if (vx[i] > max_x)
				max_x = vx[i];

			if (vy[i] < min_y)
				min_y = vy[i];
			if (vy[i] > max_y)
				max_y = vy[i];
		}

		min_x -= (max_x - min_x)*0.1;
		min_y -= (max_y - min_y)*0.1;
		max_x += (max_x - min_x) * 2 / tx;
		max_y += (max_y - min_y) * 2 / ty;

	}


	EASY_TEX tex;
	tex.resize(tx, ty);

	for (int i = 0; i < tx; i++)
		for (int r = 0; r < ty; r++)
			tex.setpixel(i, r, 0, 0, 0);


	//теперь выведем циферки и оси

	double sx, sy;
	int porx, pory;
	sx = abs(max_x - min_x);
	sy = abs(max_y - min_y);

	porx = 10;//10^10
	while (sx / pow1((double)10, porx) < 3 && porx < 100)
		porx--;
	pory = 10;//10^10
	while (sy / pow1((double)10, pory) < 3 && pory < 100)
		pory--;



	float mnox = pow10(porx);
	float mnoy = pow10(pory);
	if (max_x / mnox<1000)
	if (mnox != 0)
		for (int i = min_x / mnox; i < max_x / mnox; i++)
		{
			tex.line((i * mnox - min_x) / sx * tx, 40, (i * mnox - min_x) / sx * tx, ty, 50, 50, 50);

			if ((i * mnox - min_x) / sx * tx>20 && (i * mnox - min_x) / sx * tx < tx - 30)
				tex.numbers((i * mnox - min_x) / sx * tx - 3, 20, ftos(i*mnox));
		}
	if (max_y / mnoy<1000)
	if (mnoy != 0)
		for (int i = min_y / mnoy; i < max_y / mnoy; i++)
		{
			tex.line(20, (i * mnoy - min_y) / sy * ty, tx, (i * mnoy - min_y) / sy * ty, 50, 50, 50);

			if ((i * mnoy - min_y) / sy * ty>20 && (i * mnoy - min_y) / sy * ty < ty - 30)
				tex.numbers(4, (i * mnoy - min_y) / sy * ty + 5, ftos(i*mnoy));
		}
	//tex.numbers(30, 20, "111111");


	for (int i = 1; i <vx.size(); i++)
	{
		double rx[2], ry[2];
		double px[2], py[2];

		for (int r = 0; r < 2; r++)
		{
			rx[r] = vx[i - 1 + r];
			ry[r] = vy[i - 1 + r];
			px[r] = (rx[r] - min_x) / (max_x - min_x)*tx;
			py[r] = (ry[r] - min_y) / (max_y - min_y)*ty;
		}
		tex.line(px[0], py[0], px[1], py[1]);
		//tex.numbers(px[0], py[0],ftos(py[0]));
	}

	return tex;
}


EASY_TEX create_diskr_plot(int tx, int ty, vector<double> vx, vector<double> vy, double min_x, double max_x, double min_y, double max_y)
{


	if (min_x == min_y)
	{
		max_x = vx[0];
		min_x = vx[0];
		max_y = vy[0];
		min_y = vy[0];
		for (int i = 0; i < vx.size(); i++)
		{
			if (vx[i] < min_x)
				min_x = vx[i];
			if (vx[i] > max_x)
				max_x = vx[i];

			if (vy[i] < min_y)
				min_y = vy[i];
			if (vy[i] > max_y)
				max_y = vy[i];
		}

		if (min_y > 0)
			min_y = 0;

		max_x += 0.5;
		min_x -= (max_x - min_x)*0.1;
		min_y -= (max_y - min_y)*0.1;
		max_x += (max_x - min_x) * 0.05;
		max_y += (max_y - min_y) * 0.05;


	}


	EASY_TEX tex;
	tex.resize(tx, ty);

	for (int i = 0; i < tx; i++)
		for (int r = 0; r < ty; r++)
			tex.setpixel(i, r, 0, 0, 0);


	//теперь выведем циферки и оси

	double sx, sy;
	int porx, pory;
	sx = abs(max_x - min_x);
	sy = abs(max_y - min_y);

	porx = 10;//10^10
	while (sx / pow1((double)10, porx) < 3 && porx < 100)
		porx--;
	pory = 10;//10^10
	while (sy / pow1((double)10, pory) < 3 && pory < 100)
		pory--;



	float mnox = pow10(porx);
	float mnoy = pow10(pory);

	if (mnox != 0)
		for (int i = min_x / mnox; i < max_x / mnox; i++)
		{
			tex.line((i * mnox - min_x) / sx * tx, 40, (i * mnox - min_x) / sx * tx, ty, 50, 50, 50);

			if ((i * mnox - min_x) / sx * tx>20 && (i * mnox - min_x) / sx * tx < tx - 30)
				tex.numbers((i * mnox - min_x) / sx * tx - 3, 20, ftos(i*mnox));
		}
	if (mnoy != 0)
		for (int i = min_y / mnoy; i < max_y / mnoy; i++)
		{
			tex.line(20, (i * mnoy - min_y) / sy * ty, tx, (i * mnoy - min_y) / sy * ty, 50, 50, 50);

			if ((i * mnoy - min_y) / sy * ty>20 && (i * mnoy - min_y) / sy * ty < ty - 30)
				tex.numbers(4, (i * mnoy - min_y) / sy * ty + 5, ftos(i*mnoy));
		}
	//tex.numbers(30, 20, "111111");


	for (int i = 0; i <vx.size(); i++)
	{
		double rx[2], ry[2];
		double px[2], py[2];

		for (int r = 1; r < 2; r++)
		{
			rx[r] = vx[i - 1 + r];
			ry[r] = vy[i - 1 + r];
			px[r] = (rx[r] - min_x) / (max_x - min_x)*tx;
			py[r] = (ry[r] - min_y) / (max_y - min_y)*ty;
		}
		px[0] = (vx[i] - min_x) / (max_x - min_x)*tx;
		py[0] = (0 - min_y) / (max_y - min_y)*ty;
		tex.line(px[0], py[0], px[1], py[1]);
		for (int r = 0; r < 12; r++)
		{
			tex.line(px[1] + sin(r / 12.0 * M_PI * 2) * 3, py[1] + cos(r / 12.0 * M_PI * 2) * 3,
				px[1] + sin((r + 1) / 12.0 * M_PI * 2) * 3, py[1] + cos((r + 1) / 12.0 * M_PI * 2) * 3);
		}
		//tex.numbers(px[0], py[0],ftos(py[0]));
	}

	return tex;
}

void task_graphix()
{
	OPENGL_WINDOW w;



	w.enable();
	//tex = create_plot(512, 512, {0,1,2,3}, {2,0,1,4},-1,4,-1,5);


	double avg_time = 0;
	double scroll = 0;
	int tex_ssy = 0;
	int lasttime = time(nullptr);

	int imgz = 0;

	while (1)
	{

		avg_time = avg_time*0.9 + (GetTickCount() - lasttime)*0.1;
		lasttime = GetTickCount();
		glColor3f(0.3, 0.3, 0.3);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glLoadIdentity();
		glClearColor(0, 0, 0, 0);




		glColor3f(1,1,1);

		mtx.lock();

		for (; imgz < img.a.size(); imgz++)
		{
			string ttt;
			ttt = "img" + ftos(imgz) + ".bmp";
			vector<byte> v,v1;
			v1=v = img.a[imgz].gdata();
			for (int r = 0; r < img.a[imgz].gy(); r++)
				for (int i = 0; i < img.a[imgz].gx(); i++)
				{
					v[(i + r*img.a[imgz].gx()) * 3] = v1[(i + (img.a[imgz].gy() - r - 1)*img.a[imgz].gx()) * 3];
					v[(i + r*img.a[imgz].gx()) * 3+1] = v1[(i + (img.a[imgz].gy() - r - 1)*img.a[imgz].gx()) * 3+1];
					v[(i + r*img.a[imgz].gx()) * 3+2] = v1[(i + (img.a[imgz].gy() - r - 1)*img.a[imgz].gx()) * 3+2];
				}
			SaveArrFile(wstring(ttt.begin(), ttt.end()).c_str(), (int *)v.data(), img.a[imgz].gx(), img.a[imgz].gy(),24);
		}

		int yy = 0;
		for (int i = 0; i < img.a.size(); i++)
		{
			img.a[i].gentex();
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, img.a[i].ID());

			int wx = w.get_wx(), wy = w.get_wy();

			int tx, ty;
			tx = img.a[i].gx();
			ty = img.a[i].gy();

			glBegin(GL_POLYGON);
			glTexCoord2d(0, 1); glVertex3f(-wx + 000000, wy - 000000 - yy * 2 + scroll * 2, -wy);
			glTexCoord2d(1, 1); glVertex3f(-wx + tx * 2, wy - 000000 - yy * 2 + scroll * 2, -wy);
			glTexCoord2d(1, 0); glVertex3f(-wx + tx * 2, wy - ty * 2 - yy * 2 + scroll * 2, -wy);
			glTexCoord2d(0, 0); glVertex3f(-wx + 000000, wy - ty * 2 - yy * 2 + scroll * 2, -wy);
			glEnd();

			yy += img.a[i].gy()+10;
			
		}
		tex_ssy = yy-10;
		mtx.unlock();
		/*
		glBegin(GL_POLYGON);
		glTexCoord2d(0, 0); glVertex3f(-1, -1, -1);
		glTexCoord2d(1, 0); glVertex3f(1, -1, -1);
		glTexCoord2d(1, 1); glVertex3f(1, 1, -1);
		glTexCoord2d(0, 1); glVertex3f(-1, 1, -1);
		glEnd();
		*/
		glDisable(GL_TEXTURE_2D);

		glPushMatrix();
		glTranslatef(0, 0, -4);


		//put_icon();
		glPopMatrix();

		Sleep(10);
		w.draw();
		w.upd();
		const bool *keys=w.get_keys();

		if (keys[VK_UP])
		{
			scroll -= avg_time;
			if (scroll < 0)
				scroll = 0;
		}
		if (keys[VK_DOWN])
		{
			scroll += avg_time;
			if (scroll + w.get_wy()>tex_ssy)
				scroll = tex_ssy - w.get_wy();
			if (scroll < 0)
				scroll = 0;
		}

		if (!w.is_enabled())
		{
			w.enable();
			for (int i = 0; i < img.a.size(); i++)
				img.a[i].gen_broken_tex();

		}
	}
}

class COMP_2_RES
{
public:
	string f1_s, f1_t, f2_s, f2_t;
	string h1_s, h1_t;
	string s0, s9, s0s, s9s;
	string h1_0, h1_9, h1_0s, h1_9s;
	vector<complex<double>> v_polus, v_zero;
};

COMP_2_RES comp_2(EL_CHAIN cha,int id_res,int t_s, double im, double ti)
{
	COMP_2_RES res;
	EASY_TEX tex;

	res.h1_s = cha.comp_h1_l(id_res);

	EL_CHAIN cha_temp;

	cha_temp = cha;
	cha_temp.replace_el(EL_I, EL_I,1);
	cha_temp.replace_el(EL_U, EL_U,1);
	cha_temp.replace_el(EL_L, EL_K);
	cha_temp.replace_el(EL_C, EL_H);
	cha_temp.comp_par_1_iu_uns();

	tex.resize(512, 256);
	print_chem(tex, cha_temp);
	tex.numbers(10, tex.gy() - 20, "это схема при эс равному нулю");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	double HS_0 = (cha_temp.el[id_res].I);


	cha_temp = cha;
	cha_temp.replace_el(EL_I, EL_I, 1);
	cha_temp.replace_el(EL_U, EL_U, 1);
	cha_temp.replace_el(EL_L, EL_H);
	cha_temp.replace_el(EL_C, EL_K);
	cha_temp.comp_par_1_iu_uns();

	tex.resize(512, 256);
	print_chem(tex, cha_temp);
	tex.numbers(10, tex.gy() - 20, "это схема при эс стремящемуся");
	tex.numbers(10, tex.gy() - 20-18, "к бесконечности");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	double HS_9 = (cha_temp.el[id_res].I);

	string HS_0S = sympy_lim("(" + res.h1_s + ")*s", "s", "0","+");
	string HS_9S = sympy_lim("(" + res.h1_s + ")*s", "s", "oo","-");



	if (HS_0S.find("nan") != -1 || HS_0S.find("oo") != -1)
	{
		cout << "HS_0S:" << HS_0S << endl;
		system("pause");
	}
	if (HS_0S.find("nan") != -1 || HS_0S.find("oo") != -1)
	{
		cout << "HS_9S:" << HS_9S << endl;
		system("pause");
	}

	if (HS_0 == 0 && (atof1(HS_0S) - HS_0)>0.000001 || (atof1(HS_0S) - HS_0) / HS_0>0.001)
	{
		cout << "err: h0!=H(0)" << endl;
		cout << "h1:" << HS_0 << endl;
		cout << "hs:" << HS_0S << endl;
		system("pause");
	}
	if (HS_9 == 0 && (atof1(HS_9S) - HS_9)>0.000001 || (atof1(HS_9S) - HS_9) / HS_9>0.001)
	{
		cout << "err: h9!=H(9)" << endl;
		cout << "h1:" << HS_9 << endl;
		cout << "hs:" << HS_9S << endl;
		system("pause");
	}

	res.s0 = ftos(HS_0);
	res.s9 = ftos(HS_9);
	res.s0s = HS_0S;
	res.s9s = HS_9S;

	res.h1_t = crazy_alap(res.h1_s);
	if (res.h1_t == "")
		res.h1_t = sympy_alap(res.h1_s);



	string po1, po2;
	string f = sympy_sim("("+res.h1_s+")*s");

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


	auto roo1 = parse_roots(sympy_roo(po1));
	auto roo2 = parse_roots(sympy_roo(po2));



	res.v_polus = roo2;
	res.v_zero = roo1;



	vector<string> a_t, a_l, a_res;
	vector<double> t_sdv;



	if (t_s == SIGN_D)//2 анти ступени
	{
		a_t.resize(3);
		a_l.resize(3);
		a_res.resize(3);
		t_sdv.resize(3);
		//сдвиг по времени на программном уровне для более красивого вывода
		a_t[0] = "Heaviside(t)*(" + ftos(im) + ")";
		t_sdv[0] = 0;
		a_t[1] = "-Heaviside(t) * 2 * (" + ftos(im) + ")";
		t_sdv[1] = ti / 2;
		a_t[2] = "Heaviside(t)*(" + ftos(im) + ")";
		t_sdv[2] = ti;
		//cout << "\n" << b_s;



		//b_l = sympy_lap(b_s);
		//cout << "\n" << b_l;
	}

	if (t_s == SIGN_V)//пирамидка
	{
		a_t.resize(3);
		a_l.resize(3);
		a_res.resize(3);
		t_sdv.resize(3);
		//сдвиг по времени на программном уровне для более красивого вывода
		a_t[0] = "Heaviside(t)*t*(" + ftos(im / ti * 2) + ")";
		t_sdv[0] = 0;
		a_t[1] = "-Heaviside(t) * t * 2 * (" + ftos(im / ti * 2) + ")";
		t_sdv[1] = ti / 2;
		a_t[2] = "Heaviside(t)*t*(" + ftos(im / ti * 2) + ")";
		t_sdv[2] = ti;
		//cout << "\n" << b_s;



		//b_l = sympy_lap(b_s);
		//cout << "\n" << b_l;
	}

	if (t_s == SIGN_Z)//cтупенька
	{
		a_t.resize(2);
		a_l.resize(2);
		a_res.resize(2);
		t_sdv.resize(2);
		//сдвиг по времени на программном уровне для более красивого вывода
		a_t[0] = "Heaviside(t)*(" + ftos(im ) + ")";
		t_sdv[0] = 0;
		a_t[1] = "-Heaviside(t)* (" + ftos(im) + ")";
		t_sdv[1] = ti;
		//cout << "\n" << b_s;



		//b_l = sympy_lap(b_s);
		//cout << "\n" << b_l;
	}


	for (int i = 0; i < a_t.size(); i++)
	{
		a_l[i] = sympy_lap(a_t[i]);
		string ts;
		if (t_sdv[i] != 0)
			ts = "(" + a_l[i] + ")*exp(-s*(" + ftos(t_sdv[i]) + "))";
		else
			ts = "(" + a_l[i] + ")";


	}

	res.f1_t = "0";
	res.f1_s = "0";
	for (int i = 0; i < a_t.size(); i++)
	{
		string temp;

		temp = a_t[i];
		myreplace(temp, "t", "(t-(" + ftos(t_sdv[i]) + "))");
		res.f1_t = res.f1_t + "+(" + temp + ")";

		temp = a_l[i];
		res.f1_s = res.f1_s + "+(" + temp + ")*exp(-s*(" + ftos(t_sdv[i]) + "))";
	}

	L_S H;
	H = res.h1_s*(L_S)"s";


	res.f2_s = "";
	for (int i = 0; i < a_t.size(); i++)
	{
		L_S temp;
		temp = (a_l[i] * H).s;

		//Формула сдвига лапласса
		//f(t-a)=exp(-a*s)*F(s)

		if (i > 0)
			res.f2_s += "+";


		if (t_sdv[i] != 0)
			res.f2_s += "(" + temp.s + ")*exp(-s*(" + ftos(t_sdv[i]) + "))";
		else
			res.f2_s += "(" + temp.s + ")";
		//
		string rr = crazy_alap(temp.s);
		if (rr == "")
			rr = sympy_alap(temp.s);
		a_res[i] = rr;
	}



	res.f2_t = "";

	for (int i = 0; i < a_t.size(); i++)
	{
		string temp = a_res[i];
		myreplace(temp, "exp", "hidden_exp");//этот симпи не должен изгаживать экспоненты
		temp = sympy_sim(temp);
		myreplace(temp, "hidden_exp", "exp");

		//тэшку только в конце переписываем а то скобки раскроются
		if (t_sdv[i] != 0)
			myreplace(temp, "t", "(t-" + ftos(t_sdv[i]) + ")");

		myreplace(temp, "((t-" + ftos(t_sdv[i]) + "))", "(t-" + ftos(t_sdv[i]) + ")");

		//myreplace(temp, "Heaviside", "b1");

		if (i > 0)
			res.f2_t += "+";
		res.f2_t += "(" + temp + ")";


		myreplace(a_res[i], "t", "(t-" + ftos(t_sdv[i]) + ")");
	}

	//Определить переходную ( ) 1 h t характеристику цепи, сравнить с найденной в п. 1.2 задания. 
	//Проверить ( ) 1 0h и ( ) 1h ? по аналитическому выражению ( ) 1 h t и непосредственно по схеме цепи. 


	cha_temp = cha;

	cha_temp.replace_el(EL_L, EL_H);
	cha_temp.replace_el(EL_C, EL_K);
	cha_temp.replace_el(EL_U, EL_U, 1);
	cha_temp.replace_el(EL_I, EL_I, 1);

	tex.resize(512, 256);
	print_chem(tex, cha_temp);
	tex.numbers(10, tex.gy()-20, "это схема для тэ равному 0 с плюсом");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	cha_temp.comp_par_1_iu_uns();
	double h1_0=cha_temp.el[id_res].I;



	cha_temp = cha;
	cha_temp.replace_el(EL_L, EL_K);
	cha_temp.replace_el(EL_C, EL_H);
	cha_temp.replace_el(EL_U, EL_U, 1);
	cha_temp.replace_el(EL_I, EL_I, 1);

	cha_temp.comp_par_1_iu_uns();
	double h1_9 = cha_temp.el[id_res].I;

	tex.resize(512, 256);
	print_chem(tex, cha_temp);
	tex.numbers(10, tex.gy() - 20, "это схема при тэ стремящемуся");
	tex.numbers(10, tex.gy() - 20 - 18, "к бесконечности");
	mtx.lock();
	img.add(tex);
	mtx.unlock();


	string h1_0s = sympy_lim(res.h1_t, "t", "0", "+");
	string h1_9s = sympy_lim(res.h1_t, "t", "oo", "-");


	res.h1_0 = ftos(h1_0);
	res.h1_9 = ftos(h1_9);
	res.h1_0s = h1_0s;
	res.h1_9s = h1_9s;


	vector<double> vx(200);
	vector<double> vy;
	vector<double> vz;

	vz.resize(vx.size());
	vy.resize(vx.size());

	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*ti*2/vx.size();
		string s;
		s = sympy_eva(res.f1_t, "t", ftos(vx[i]));
		vy[i] = atof1(s);
		s = sympy_eva(res.f2_t, "t", ftos(vx[i]));
		vz[i] = atof1(s);
	}

	tex = create_double_plot(1024, 512, vx, vy, vz);
	tex.numbers(100, tex.gy() - 20, "ток реакции и сигнал на 1 графике");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	return res;
}






class POl_PROPUSK
{
public:
	int id; //0 - низких_ч, 1 - высоких_ч, 2 - полосовой, 3 - заграждающий
	double w0, w1; //-1 - бесконечность
	void show()const
	{
		double dw;
		switch (id)
		{
		case 0:
			cout << "Фильтр низких частот";
			dw = w1 - w0;
			cout << "\nПолоса пропускания dw = " << dw << " [ " << w0 << ", " << w1 << " ]\n";
			break;
		case 1:
			cout << "Фильтр высоких частот";
			dw = -1;
			cout << "\nПолоса пропускания dw = inf [ " << w0 << ", +inf ]\n";
			break;
		case 2:
			cout << "Полосовой фильтр";
			dw = w1 - w0;
			cout << "\nПолоса пропускания dw = " << dw << " [ " << w0 << ", " << w1 << " ]\n";
			break;
		case 3:
			cout << "Заграждающий фильтр";
			dw = w1 - w0;
			cout << "\nПолоса пропускания [ 0 , " << w0 << " U " << w1 << ", +inf ]\n";
			break;
		}
	}
};

class COMP_3_RES
{
public:
	string ACH_H, FCH_H;
	double max_ACH;
	double ACH_H_0, ACH_H_inf;
	POl_PROPUSK polosa;
	string AS_F1, FS_F1;
	double dw_AS_F1, mid_persent, isk_persent, AS_F2_0;
	string p34;
};

COMP_3_RES comp_3(L_S H_S, string F1_T, string F1_S, int sgn_t, double sgn_time)
{
	COMP_3_RES RES;
	//H_S это переходная функция тока в лапласе

	string H_JW, H_IM, H_RE;

	//3.1. Используя найденное в 2.1 выражение HU(s) или HI(s), вычислить и построить графики АЧХ и ФЧХ 
	//функций передачи цепи HU(jw) или HI(jw). Произвести проверку АЧХ при w = 0 и w -> inf.  
	H_JW = sympy_eva(H_S.s, "s", "I*w");
	H_IM = sympy_im(H_JW);
	H_RE = sympy_re(H_JW);

	string ACH_H = "sqrt((" + H_IM + ")**2+(" + H_RE + ")**2)";
	string FCH_H = "atan((" + H_IM + ")/(" + H_RE + "))-Heaviside(-(" + H_RE + "))*pi";


	cout << "\n\nACH_H = " << ACH_H; //АЧХ
	cout << "\n\nFCH_H = " << FCH_H; //ФЧХ

	double ACH_H_0 = atof1(sympy_eva(ACH_H, "w", "0.0000001")), ACH_H_inf = atof1(sympy_lim(ACH_H, "w", "oo", "-"));
	//проверка при w = 0
	cout << "\n\nACH_H [w = 0] = " << ACH_H_0;

	//проверка при w -> inf
	cout << "\nACH_H [w -> inf] = " << ACH_H_inf << "\n";

	RES.ACH_H = ACH_H;
	RES.FCH_H = FCH_H;
	RES.ACH_H_0 = ACH_H_0;
	RES.ACH_H_inf = ACH_H_inf;



	//3.2. Определить полосу пропускания цепи по уровню 0,707|H(jw)|[max].
	double first, last, min = 1000000, max = 0, temp;

	double x[193], x_max = 0.001, x_min = 0.001;
	for (int i = 0; i < 193; i++)
	{
		if (i == 0)
		{
			x[i] = 0.001;
			first = min = max = atof1(sympy_eva(ACH_H, "w", ftos(x[i])));
		}
		else
		{
			x[i] = x[i - 1] * 1.1;
			temp = atof1(sympy_eva(ACH_H, "w", ftos(x[i])));
			if (min > temp) { min = temp; x_min = x[i]; }
			if (max < temp) { max = temp; x_max = x[i]; }
		}

		if (i == 192)
			last = atof1(sympy_lim(ACH_H, "w", "oo", "-"));
	}

	//точное определение max и min
	//определение max
	double shag = 0.01;
	for (; shag > 0.0001; shag /= 2)
		if (atof1(sympy_eva(ACH_H, "w", ftos(x_max))) < atof1(sympy_eva(ACH_H, "w", ftos(x_max + shag))))
			x_max += shag;
		else
			if ((x_max - shag) > 0) x_max -= shag;

	max = atof1(sympy_eva(ACH_H, "w", ftos(x_max)));

	//определение min
	/*for(shag = 0.01; atof1(sympy_eva(ACH, "w", ftos(x_min))) > atof1(sympy_eva(ACH, "w", ftos(x_min + shag))); shag = x_min - x_min/1.1)
	x_min += shag;*/
	shag = 0.01;
	for (; shag > 0.0001; shag /= 2)
		if (atof1(sympy_eva(ACH_H, "w", ftos(x_min))) > atof1(sympy_eva(ACH_H, "w", ftos(x_min + shag))))
			x_min += shag;
		else
			if ((x_min - shag) > 0) x_min -= shag;

	min = atof1(sympy_eva(ACH_H, "w", ftos(x_min)));

	//cout << "\nx_max = " << x_max << " x_min = " << x_min << "\n";

	/*if (first > last && (last - pogr < min && min < last + pogr) && (first - pogr < max && max < first + pogr))
	POLOSA.id = 0;
	if (first < last && (first - pogr < min && min < first + pogr) && (last - pogr < max && max < last + pogr))
	POLOSA.id = 1;
	if (first < max && (first - pogr < last && last < first + pogr) && ((first - pogr < min && min < first + pogr) || (last - pogr < min && min < last + pogr)))
	POLOSA.id = 2;
	if (first > min && (first - pogr < last && last < first + pogr) && ((first - pogr < max && max < first + pogr) || (last - pogr < max && max < last + pogr)))
	POLOSA.id = 3;*/

	POl_PROPUSK POLOSA;
	bool b1 = (max / sqrt(2)) > first;
	bool b2 = (max / sqrt(2)) > last;
	if (b1 && b2) POLOSA.id = 2;
	if (!b1 && b2) POLOSA.id = 0;
	if (b1 && !b2) POLOSA.id = 1;
	if (!b1 && !b2) POLOSA.id = 3;

	RES.max_ACH = atof1(sympy_eva(ACH_H, "w", ftos(x_max)));

	double y = 0.707 * atof1(sympy_eva(ACH_H, "w", ftos(x_max)));
	switch (POLOSA.id)
	{
	case 0:
		POLOSA.w0 = 0;
		for (shag = 0.01, POLOSA.w1 = x_max; atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w1))) > y; shag *= 1.1)
			POLOSA.w1 += shag;
		for (; shag > 0.0001; shag /= 2)
			if (atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w1))) > y)
				POLOSA.w1 += shag;
			else
				POLOSA.w1 -= shag;
		break;

	case 1:
		POLOSA.w1 = -1;
		for (shag = 0.01, POLOSA.w0 = 0.001; atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w0))) < y; shag *= 1.1)
			POLOSA.w0 += shag;
		for (; shag > 0.0001; shag /= 2)
			if (atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w0))) < y)
				POLOSA.w0 += shag;
			else
				POLOSA.w0 -= shag;
		break;

	case 2:
		for (shag = 0.01, POLOSA.w0 = x_max; atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w0))) > y; shag *= 1.1)
			POLOSA.w0 -= shag;
		for (; shag > 0.0001; shag /= 2)
			if (atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w0))) < y)
				POLOSA.w0 += shag;
			else
				POLOSA.w0 -= shag;

		for (shag = 0.01, POLOSA.w1 = x_max; atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w1))) > y; shag *= 1.1)
			POLOSA.w1 += shag;
		for (; shag > 0.0001; shag /= 2)
			if (atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w1))) > y)
				POLOSA.w1 += shag;
			else
				POLOSA.w1 -= shag;
		break;

	case 3:
		for (shag = 0.01, POLOSA.w0 = 0.001; atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w0))) < y; shag *= 1.1)
			POLOSA.w0 += shag;
		for (; shag > 0.0001; shag /= 2)
			if (atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w0))) > y)
				POLOSA.w0 += shag;
			else
				POLOSA.w0 -= shag;

		for (shag = 0.01, POLOSA.w1 = POLOSA.w0; atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w1))) < y; shag *= 1.1)
			POLOSA.w1 += shag;
		for (; shag > 0.0001; shag /= 2)
			if (atof1(sympy_eva(ACH_H, "w", ftos(POLOSA.w1))) < y)
				POLOSA.w1 += shag;
			else
				POLOSA.w1 -= shag;
		break;
	}

	//cout << " first = " << first << " last = " << last << " min = " << min << " max = " << max << endl;
	POLOSA.show();

	RES.polosa = POLOSA;


	vector<double> vx, vy;
	EASY_TEX tex;
	vx.resize(200);
	vy.resize(200);
	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = (i) / 200.0*(2 * max(POLOSA.w0, POLOSA.w1));
		if (i == 0) vx[i] = 0.0001*max(POLOSA.w0, POLOSA.w1);

		string s = sympy_eva(ACH_H, "w", ftos(vx[i]));
		//cout << "st:" << s;
		myreplace(s, ".", ",");
		vy[i] = atof(s.c_str());
	}


	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "ачх");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();


	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = (i) / 200.0*(2 * max(POLOSA.w0, POLOSA.w1));
		if (i == 0) vx[i] = 0.0001*max(POLOSA.w0, POLOSA.w1);

		string s = sympy_eva(FCH_H, "w", ftos(vx[i]));
		//cout << "st:" << s;
		myreplace(s, ".", ",");
		vy[i] = atof(s.c_str());
		while (vy[i] > M_PI)
			vy[i] -= M_PI * 2;
		while (vy[i] < -M_PI)
			vy[i] += M_PI * 2;
	}


	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "фчх");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();



	//3.3. Найти и построить амплитудный и фазовый спектры входного одиночного импульса.
	//Найти ширину амплитудного спектра по уровню ( )1 max0,1 F j? или критерию, предложенному преподавателем.
	string F1_JW, F1_RE, F1_IM;
	F1_JW = sympy_eva(F1_S, "s", "I*w");
	F1_RE = sympy_re(F1_JW);
	F1_IM = sympy_im(F1_JW);

	string ACH_F1 = "sqrt((" + F1_IM + ")**2+(" + F1_RE + ")**2)";
	string FCH_F1 = "atan((" + F1_IM + ")/(" + F1_RE + "))-Heaviside(-(" + F1_RE + "))*pi";

	x_max = 0; y = 0;

	
	if (sgn_t == SIGN_V || sgn_t == SIGN_D)
	{
		if (sgn_t == SIGN_V)
			RES.dw_AS_F1 = 18.548649279553697*0.5 / sgn_time;
		if (sgn_t == SIGN_D)
			RES.dw_AS_F1 = 456.94248951839495*0.1 / sgn_time;
	}
	else
	{
		int i_max = 0;
		y = atof1(sympy_eva(ACH_F1, "w", ftos(x[0])));
		x_max = x[0];
		for (int i = 0; i < 193; i++)
		{
			temp = atof1(sympy_eva(ACH_F1, "w", ftos(x[i])));
			if (y < temp) { y = temp; x_max = x[i]; i_max = i; }
		}

		cout << "y = " << y << " x_max = " << x_max << " i_max = " << i_max << "\n";
		shag = x[i_max + 1] - x[i_max];
		for (; shag > 0.00001; shag /= 2)
			if (atof1(sympy_eva(ACH_F1, "w", ftos(x_max))) < atof1(sympy_eva(ACH_F1, "w", ftos(x_max + shag))))
				x_max += shag;
			else
				x_max -= shag;
		y = atof1(sympy_lim(ACH_F1, "w", ftos(x_max), "+"));
		//cout << "y = " << y << " x_max = " << x_max << "\n";

		cout << "\n\nAS_F1 = " << ACH_F1;
		cout << "\n\nFS_F1 = " << FCH_F1 << "\n\n";

		double dw_ACH_F1;
		y = 0.1*atof1(sympy_eva(ACH_F1, "w", ftos(x_max)));
		//cout << "0.1A = " << y << "\n";
		dw_ACH_F1 = x_max;
		do
		{
			for (shag = 0.1; atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) > atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1 + shag))); /*shag *= 1.1*/)
			{
				dw_ACH_F1 += shag;
			}

			//cout << "----------\nACH_F1_1 = " << atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) << "\tACH_F1_2 = " << atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1 + 0.01))) << "\n";

			for (shag = 0.1; atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) < atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1 + shag))); /*shag *= 1.1*/)
			{
				dw_ACH_F1 += shag;
			}
			for (; shag > 0.00001; shag /= 2)
				if (atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) < atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1 + shag))))
					dw_ACH_F1 += shag;
				else
					dw_ACH_F1 -= shag;
			if (atof1(sympy_eva(ACH_F1, "w", ftos(x_max))) < atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1)))) x_max = dw_ACH_F1;
			//cout << atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) << " | " << y << "\n";
		} while (atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) > y);

		dw_ACH_F1 = x_max;
		//cout << "ACH_F1_1 = " << atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) << "\t dw_ACH_F1 = " << dw_ACH_F1 << "\n";
		for (shag = 0.01; atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) > y; shag *= 1.1)
		{
			dw_ACH_F1 += shag;
		}
		for (; shag > 0.000001; shag /= 2)
			if (atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) < y)
				dw_ACH_F1 -= shag;
			else
				dw_ACH_F1 += shag;
		cout << "ACH_max = " << atof1(sympy_eva(ACH_F1, "w", ftos(x_max))) << "\t0.1A = " << y << "\nШирина спектра dw_ACH_F1 = " << dw_ACH_F1 << "\tACH_F1 = " << atof1(sympy_eva(ACH_F1, "w", ftos(dw_ACH_F1))) << "\n";

		RES.dw_AS_F1 = dw_ACH_F1;

	}
	
	for (int i = 0; i < vx.size(); i++)
	{
		//vx[i] = i / 60.0;
		vx[i] = (i) / (float)vx.size()*RES.dw_AS_F1 * 2;
		if (i == 0) vx[i] = 0.0001*RES.dw_AS_F1 * 2;
		string s = sympy_eva(ACH_F1, "w", ftos(vx[i]));
		//cout << "\nvx[i]: " << vx[i] << "\t| st: " << s;
		myreplace(s, ".", ",");
		vy[i] = atof(s.c_str());
	}

	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "ас ф1");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	for (int i = 0; i < vx.size(); i++)
	{
		//vx[i] = i / 60.0;
		//if(i == 0) vx[i] += 0.001;
		vx[i] = (i) / (float)vx.size()*RES.dw_AS_F1 * 2;
		if (i == 0) vx[i] = 0.0001*RES.dw_AS_F1 * 2;
		string s = sympy_eva(FCH_F1, "w", ftos(vx[i]));
		//cout << "\nvx[i]: " << vx[i] << "\t| st: " << s;
		myreplace(s, ".", ",");
		vy[i] = atof(s.c_str());
		while (vy[i] > M_PI)
			vy[i] -= M_PI * 2;
		while (vy[i] < -M_PI)
			vy[i] += M_PI * 2;
	}

	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "фс ф1");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();



	//3.4. Сопоставить спектры входного импульса с частотными характеристиками цепи. Дать предварительное заключение
	//об ожидаемых искажениях сигнала на выходе цепи. Сравнить эти качественные оценки с сигналом на выходе, полученным в п. 2.5 задания.

	RES.p34 = "";

	//Приближённая оценка реакции по значению АЧХ цепи на нулевой частоте
	if (-0.0001 < ACH_H_0 && ACH_H_0 < 0.0001)
	{
		RES.p34 += "Суммарная площадь реакции равна нулю\n";
		cout << "Суммарная площадь реакции равна нулю\n";
	}
	else
	{
		RES.p34 += "Площадь реакции в " + ftos(ACH_H_inf) + " раз отличается от площади воздействия\n";
		cout << "Площадь реакции в " << ACH_H_0 << " раз отличается от площади воздействия\n";
	}

	//Приближённая оценка реакции по значению АЧХ цепи на бесконечной частоте
	if (-0.0001 < ACH_H_inf && ACH_H_inf < 0.0001)
	{
		RES.p34 += "Скачок воздействия не пройдёт на выход\n";
		cout << "Скачок воздействия не пройдёт на выход\n";
	}
	else
	{
		RES.p34 += "Скачок реакции в " + ftos(ACH_H_inf) + " раз отличается от скачка воздействия\n";
		cout << "Скачок реакции в " << ACH_H_inf << " раз отличается от скачка воздействия\n";
	}

	//Оценка по АЧХ цепи и АС воздействия

	double w[200], q = RES.dw_AS_F1 / 100, mid_persent = 0, Ew = 0;
	cout << "q = " << q << "\n";
	w[0] = 0;
	/*mid_persent = atof1(sympy_eva(ACH_F1, "w", "0.00001")) / (atof1(sympy_eva(ACH_H, "w", "0.00001")) * atof1(sympy_eva(ACH_F1, "w", "0.00001")));
	cout << "w[0] = " << w[0] << "\n";
	cout << "mid_persent = " << atof1(sympy_eva(ACH_F1, "w", "0.00001")) << " / ( " << atof1(sympy_eva(ACH_H, "w", "0.00001")) << " * " << atof1(sympy_eva(ACH_F1, "w", "0.00001")) << " )\n";*/
	for (int i = 1; i < 200; i++)
	{
		w[i] = w[i - 1] + q;
		/*
		cout << "w[" << i << "] = " << w[i] << "\n";
		cout << "mid_persent = " << mid_persent << " + " << atof1(sympy_eva(ACH_F1, "w", ftos(w[i]))) << " / ( " << atof1(sympy_eva(ACH_H, "w", ftos(w[i]))) << " * " << atof1(sympy_eva(ACH_F1, "w", ftos(w[i]))) << " )\n";
		*/
		mid_persent = mid_persent + 1 / w[i] * (atof1(sympy_eva(ACH_H, "w", ftos(w[i]))) * atof1(sympy_eva(ACH_F1, "w", ftos(w[i])))) / atof1(sympy_eva(ACH_F1, "w", ftos(w[i])));
		Ew += 1/w[i];
	}
	mid_persent = mid_persent * 100 / Ew;
	cout << "mid_persent[AS_F2/AS_F1] = " << mid_persent << " %\n";
	RES.mid_persent = mid_persent;
	RES.p34 += "Примерный процент искажения по амплитуде:" + ftos(100-mid_persent) + "\n";

	q = RES.dw_AS_F1 / 100;
	mid_persent = 0;
	Ew = 0;

	cout << "q = " << q << "\n";
	w[0] = 0;
	double angdiff;
	for (int i = 1; i < 200; i++)
	{
		w[i] = w[i - 1] + q;
		//min(abs(1 - A1 / A2) * 2 + angdiff(ф1 - ф2) / 90, 1) 
		angdiff = (-1)*atof1(sympy_eva(FCH_H, "w", ftos(w[i])));
		while (angdiff > M_PI)
			angdiff -= M_PI * 2;
		while (angdiff < -M_PI)
			angdiff += M_PI * 2;
		mid_persent = mid_persent + (atof1(sympy_eva(ACH_F1, "w", ftos(w[i]))) / w[i]) * min(abs(1 - atof1(sympy_eva(ACH_H, "w", ftos(w[i]))) / (RES.max_ACH)) * 2 + angdiff / 90, 1);
		Ew += atof1(sympy_eva(ACH_F1, "w", ftos(w[i]))) / w[i];
	}
	mid_persent = mid_persent * 100 / Ew;
	cout << "mid_persent = " << mid_persent << " %\n";
	
	RES.p34 += "Примерное искажение :" + ftos(mid_persent) + "\n";

	RES.isk_persent = mid_persent;


	//3.5. Построить графики амплитудного и фазового спектров выходного сигнала, используя графики пп. 3.1, 3.3 задания.
	//Проконтролировать площадь реакции по значению ее спектра при w = 0.
	for (int i = 0; i < vx.size(); i++)
	{
		//vx[i] = i / 60.0;
		//if(i == 0) vx[i] += 0.001;
		vx[i] = (i) / (float)vx.size()*RES.dw_AS_F1 * 2;
		if (i == 0) vx[i] = 0.0001*RES.dw_AS_F1 * 2;

		string s = ftos(atof1(sympy_eva(ACH_H, "w", ftos(vx[i]))) * atof1(sympy_eva(ACH_F1, "w", ftos(vx[i]))));
		//cout << "st:" << s;
		myreplace(s, ".", ",");
		vy[i] = atof(s.c_str());
	}

	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "ас ф2");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();


	for (int i = 0; i < vx.size(); i++)
	{
		//vx[i] = i / 60.0;
		vx[i] = (i) / (float)vx.size()*RES.dw_AS_F1 * 2;
		if (i == 0) vx[i] = 0.0001*RES.dw_AS_F1 * 2;

		string s = ftos(atof1(sympy_eva(FCH_H, "w", ftos(vx[i]))) + atof1(sympy_eva(FCH_F1, "w", ftos(vx[i]))));
		//cout << "st:" << s;
		myreplace(s, ".", ",");
		vy[i] = atof(s.c_str());
		while (vy[i] > M_PI)
			vy[i] -= M_PI * 2;
		while (vy[i] < -M_PI)
			vy[i] += M_PI * 2;
	}

	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "фс ф2");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	cout << "\nПостроены AS_F2, FS_F2.\n";
	double w0_AS_F2 = atof1(sympy_eva(ACH_H, "w", "0.0001")) * atof1(sympy_eva(ACH_F1, "w", "0.0001"));
	cout << "AS_F2[w = 0] = " << w0_AS_F2 << "\n";

	RES.AS_F2_0 = w0_AS_F2;
	//system("pause");
	return RES;
}







class COMP_4_RES
{
public:
	double Ak;
	int lev;
};

COMP_4_RES comp_4(L_S H_S, string F1_T, string F1_S, double T)
{
	vector<double> vx;
	vector<double> vy;
	vector<double> vz;

	EASY_TEX tex;
	COMP_4_RES res;

	cout << F1_T << endl;
	cout << F1_S << endl;
	//H_S это переходная функция тока в лапласе
	string H_JW, H_IM, H_RE;

	H_JW = sympy_eva(H_S.s, "s", "I*w");
	H_IM = sympy_im(H_JW);
	H_RE = sympy_re(H_JW);

	string ACH = "sqrt((" + H_IM + ")**2+(" + H_RE + ")**2)";
	string FCH = "atan((" + H_IM + ")/(" + H_RE + "))-Heaviside(-(" + H_RE + "))*pi";

	string F1_JW = sympy_eva(F1_S, "s", "I*w*k");
	string F1_K = sympy_eva("("+ftos(2/T)+")*("+F1_JW+")", "w", ftos(2 * M_PI / T));

	string A1_IM = sympy_im(F1_K);
	string A1_RE = sympy_re(F1_K);


	

	//4.1.Разложить в ряд Фурье заданный входной периодический сигнал.Построить его амплитудный и фазовый дискретные спектры. 
	//res.f1_A = "sqrt((" + A1_IM + ")**2+(" + A1_RE + ")**2)";
	//res.f1_F = "atan((" + A1_IM + ")/(" + A1_RE + "))-Heaviside(-(" + A1_RE + "))*pi";

	vector<double> arr_IM(50);
	vector<double> arr_RE(50);
	
	for (int i = 0; i < arr_IM.size(); i++)
	{
		string temp;
		if (i != 0)
			temp = sympy_eva(F1_K, "k", ftos(i));
		else
			temp = sympy_lim(F1_K, "k", ftos(0.001),"+");//иногда степени ехр различаются на 0.0000001 и из-за этого вылезает оо

		arr_RE[i] = atof1(sympy_re(temp));
		arr_IM[i] = atof1(sympy_im(temp));

	}

	vector<double> arr_A1_A(50);
	vector<double> arr_A1_F(50);

	for (int i = 0; i < arr_IM.size(); i++)
	{
		arr_A1_A[i] = sqrt(arr_IM[i] * arr_IM[i] + arr_RE[i] * arr_RE[i]);
		arr_A1_F[i] = atan(arr_IM[i] / arr_RE[i]);
		if (arr_RE[i] < 0)
			arr_A1_F[i] += M_PI;

		float f = arr_A1_F[i];


		if (f < 0)
		{
			int ss = 1;
			while (f < 0)
			{
				f = f + M_PI * 2 * ss;
				ss = ss << 1;
			}
		}

		int shag = 1;
		while (f > shag * 2 * M_PI)
		{
			shag <<= 1;
		}

		while (shag >= 1)
		{
			shag >>= 1;

			if (f - M_PI * 2 * shag > 0)
				f = f - M_PI * 2 * shag;
		}

		arr_A1_F[i] = f;
	}


	double Amax = 0;

	for (int i = 1; i < arr_IM.size(); i++)
	{
		if (arr_A1_A[i]>Amax)
			Amax = arr_A1_A[i];
	}

	int lev = 1;

	for (int i = 0; i < arr_IM.size(); i++)
	{
		if (arr_A1_A[i] * 10 > Amax)
			lev = i+1;
		if (arr_A1_A[i] * 10000 < Amax)
		{
			arr_A1_A[i] = 0;
			arr_A1_F[i] = 0;
		}
	}

	res.Ak = Amax;
	res.lev = lev;
	vx.resize(lev);
	for (int i = 0; i < lev; i++)
		vx[i] = i;
	tex = create_diskr_plot(1024, 512, vx, arr_A1_A, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "амплитудная характеристика сигнала");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();


	vx.resize(lev);
	for (int i = 0; i < lev; i++)
		vx[i] = i;
	tex = create_diskr_plot(1024, 512, vx, arr_A1_F, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "частотная характеристика сигнала");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();
	//4.2.Построить на одном графике заданный входной периодический сигнал и его аппроксимацию отрезком ряда Фурье.Число гармоник 
	//отрезка ряда Фурье определяется по уровню 0, 1 km A, где km A – максимальная составляющая амплитудного спектра, 
	//или по другому критерию, предложенному преподавателем. 


	vx.resize(200);

	vz.resize(vx.size());
	vy.resize(vx.size());

	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*T / vx.size();
		string s;
		s = sympy_eva(F1_T, "t", ftos(vx[i]));
		vy[i] = atof1(s);

		double fr = arr_A1_A[0]/2;
		for (int r = 1; r < lev; r++)//есть способ лучше но и так сойдёт ещё можно протейлорить косинус
		{
			fr = fr + arr_A1_A[r] * cos(1 / T * 2 * r *M_PI* vx[i] + arr_A1_F[r]);
		}
		vz[i] = fr;
	}

	tex = create_double_plot(1024, 512, vx, vy, vz);
	tex.numbers(100, tex.gy() - 20, "исходный сигнал");
	tex.numbers(100, tex.gy() - 20 - 18, "и его приближение");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	//4.3.Используя рассчитанные в п. 3.1 задания АЧХ и ФЧХ, найти реакцию цепи в виде отрезка ряда Фурье с 
	//числом гармоник, определенным для входного сигнала. 

	vector<double> arr_A2_A(50);
	vector<double> arr_A2_F(50);

	for (int r = 0; r < lev; r++)
	{
		double HA, HF;
		HA = 0;
		HF = 0;
		if (r > 0)
		{
			HA = atof1(sympy_eva(ACH, "w", ftos(r*M_PI * 2 / T)));
			HF = atof1(sympy_eva(FCH, "w", ftos(r*M_PI * 2 / T)));
		}
		else
		{
			HA = atof1(sympy_eva(ACH, "w", ftos(0.001*M_PI * 2 / T)));
			HF = atof1(sympy_eva(FCH, "w", ftos(0.001*M_PI * 2 / T)));
		}

		float f = HF;


		if (f < 0)
		{
			int ss = 1;
			while (f < 0)
			{
				f = f + M_PI * 2 * ss;
				ss = ss << 1;
			}
		}

		int shag = 1;
		while (f > shag * 2 * M_PI)
		{
			shag <<= 1;
		}

		while (shag >= 1)
		{
			shag >>= 1;

			if (f - M_PI * 2 * shag > 0)
				f = f - M_PI * 2 * shag;
		}

		HF = f;

		//if (r > 0)
		//{
		//	arr_A2_A[r] = arr_A1_A[r] * HA;
		//	arr_A2_F[r] = arr_A1_F[r] + HF;
		//}
		//else
		{
			arr_A2_A[r] = arr_A1_A[r] * HA;
			arr_A2_F[r] = arr_A1_F[r] + HF;

			float f = arr_A2_F[r];


			if (f < 0)
			{
				int ss = 1;
				while (f < 0)
				{
					f = f + M_PI * 2 * ss;
					ss = ss << 1;
				}
			}

			int shag = 1;
			while (f > shag * 2 * M_PI)
			{
				shag <<= 1;
			}

			while (shag >= 1)
			{
				shag >>= 1;

				if (f - M_PI * 2 * shag > 0)
					f = f - M_PI * 2 * shag;
			}

			arr_A2_F[r] = f;
			/*
			double val = 0;
			if (HF > M_PI / 2 && HF < 3 * M_PI / 2)
				val = val - HA;
			else
				val = val + HA;

			if (arr_A2_F[r] > M_PI / 2 && arr_A2_F[r] < 3 * M_PI / 2)
				val = val - arr_A2_A[r];
			else
				val = val + arr_A2_A[r];

			arr_A2_A[r] = abs(val);
			if (val>=0)
				arr_A2_F[r] = 0;
			else
				arr_A2_F[r] = M_PI;
				*/

		}

	}


	vx.resize(lev);
	for (int i = 0; i < lev; i++)
		vx[i] = i;
	tex = create_diskr_plot(1024, 512, vx, arr_A2_A, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "амплитудная характеристика выходного сигнала");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();


	vx.resize(lev);
	for (int i = 0; i < lev; i++)
		vx[i] = i;
	tex = create_diskr_plot(1024, 512, vx, arr_A2_F, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "фазовая характеристика выходного сигнала");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	vx.resize(300);
	vy.resize(300);
	vz.resize(300);

	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*T / vx.size();
		//string s;
		//s = sympy_eva(F2_T, "t", ftos(vx[i]));
		//vy[i] = atof1(s);

		double fr;
		fr = arr_A2_A[0] / 2;
		for (int r = 1; r < lev; r++)
		{
			fr = fr + arr_A2_A[r] * cos(1 / T * 2 * r *M_PI* vx[i] + arr_A2_F[r]);
		}
		//s = fr;
		vz[i] = fr;
	}

	tex = create_plot(1024, 512, vx, vz,0,0,0,0);
	tex.numbers(100, tex.gy() - 20, "выходной сигнал");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();

	//4.4.Построить амплитудный и фазовый дискретные спектры выходного сигнала.Построить график 
	//выходного сигнала, найденного в п. 4.3 задания, в одном масштабе рядом с графиком аппроксимированного входного сигнала. 
	//4.5.Дать заключение об искажении сигнала на выходе цепи.


	return res;
}

class KU_INP
{
public:
	string cha_str;
	int el_id;
	int t_sign;
	double ts, as, T;
};

void print_ku(const EL_CHAIN &cha, const COMP_2_RES &res2, const COMP_3_RES &res3, const COMP_4_RES &res4);




class GR_EL
{
public:
	double x1, x2, y1, y2;
	int t, n, p1, p2;////0-connect 1-no connect 2-R 3- ->; 4- +-; 5-L 6-C
};

void draw_el(EASY_TEX &tex, GR_EL &el)
{
	double c, v, xx, yy, x1, y1, x2, y2;
	switch (el.t)
	{
	case EL_K:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);

		tex.line(el.x1, el.y1, el.x2, el.y2);

		dd(xx, yy, c + M_PI / 4, v*0.15);

		tex.numbers((el.x1 + el.x2) / 2 + xx, (el.y1 + el.y2) / 2 + yy, ftos(el.n));
		break;
	case EL_H:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v*0.35);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);

		dd(x1, y1, c+M_PI*3/4, v*0.2);
		tex.line(el.x1+xx+x1, el.y1+yy+y1, el.x1 + xx, el.y1 + yy);

		dd(xx, yy, c, -v*0.35);
		tex.line(el.x2, el.y2, el.x2 + xx, el.y2 + yy);

		dd(x1, y1, c + M_PI * 1 / 4, v*0.2);
		tex.line(el.x2+xx+x1, el.y2+yy+y1, el.x2 + xx, el.y2 + yy);

		tex.numbers((el.x1 + el.x2) / 2, (el.y1 + el.y2) / 2, ftos(el.n));

		break;
	case EL_U:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v*0.3);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(xx, yy, c, -v*0.3);
		tex.line(el.x2, el.y2, el.x2 + xx, el.y2 + yy);
		dd(x1, y1, c, v*0.5);
		x1 += el.x1;
		y1 += el.y1;

		for (int i = 0; i < 24; i++)
		{
			dd(xx, yy, c + (i + 0) / 24.0 * 2 * M_PI, v*0.2);
			dd(x2, y2, c + (i + 1) / 24.0 * 2 * M_PI, v*0.2);
			tex.line(x1 + x2, y1 + y2, x1 + xx, y1 + yy);
		}

		dd(x1, y1, c, v*0.45);
		dd(x2, y2, c + M_PI / 2, v*0.1);
		tex.line(el.x1 + x1 + x2, el.y1 + y1 + y2, el.x1 + x1 - x2, el.y1 + y1 - y2);
		dd(x1, y1, c, v*0.61);
		dd(x2, y2, c + M_PI / 2, v*0.1);
		tex.line(el.x1 + x1 + x2, el.y1 + y1 + y2, el.x1 + x1 - x2, el.y1 + y1 - y2);
		dd(x1, y1, c, v*0.45);
		dd(x2, y2, c , v*0.1);
		tex.line(el.x1 + x1 + x2, el.y1 + y1 + y2, el.x1 + x1 - x2, el.y1 + y1 - y2);


		break;
	case EL_I:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(x1, y1, c, v*0.5);
		x1 += el.x1;
		y1 += el.y1;

		for (int i = 0; i < 24; i++)
		{
			dd(xx, yy, c + (i + 0) / 24.0 * 2 * M_PI, v*0.2);
			dd(x2, y2, c + (i + 1) / 24.0 * 2 * M_PI, v*0.2);
			tex.line(x1 + x2, y1 + y2, x1 + xx, y1 + yy);
		}

		break;
	case EL_R:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v*0.2);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(xx, yy, c, -v*0.2);
		tex.line(el.x2, el.y2, el.x2 + xx, el.y2 + yy);

		dd(x1, y1, c, v*0.2);
		dd(x2, y2, c, v*0.8);
		dd(xx, yy, c + M_PI / 2, v*0.1);

		tex.line(el.x1 + x1 + xx, el.y1 + y1 + yy, el.x1 + x1 - xx, el.y1 + y1 - yy);
		tex.line(el.x1 + x2 + xx, el.y1 + y2 + yy, el.x1 + x2 - xx, el.y1 + y2 - yy);
		tex.line(el.x1 + x1 - xx, el.y1 + y1 - yy, el.x1 + x2 - xx, el.y1 + y2 - yy);
		tex.line(el.x1 + x2 + xx, el.y1 + y2 + yy, el.x1 + x1 + xx, el.y1 + y1 + yy);

		tex.numbers((el.x1 + el.x2) / 2, (el.y1 + el.y2) / 2, ftos(el.n));
		break;
	case EL_C:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v*0.45);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(xx, yy, c, -v*0.45);
		tex.line(el.x2, el.y2, el.x2 + xx, el.y2 + yy);

		dd(x1, y1, c, v*0.45);
		dd(x2, y2, c, v*0.55);
		dd(xx, yy, c + M_PI / 2, v*0.15);

		tex.line(el.x1 + x1 + xx, el.y1 + y1 + yy, el.x1 + x1 - xx, el.y1 + y1 - yy);
		tex.line(el.x1 + x2 + xx, el.y1 + y2 + yy, el.x1 + x2 - xx, el.y1 + y2 - yy);

		dd(xx, yy, c + M_PI / 4, v*0.2);

		tex.numbers((el.x1 + el.x2) / 2+xx, (el.y1 + el.y2) / 2+yy, ftos(el.n));

		break;
	case EL_L:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v*0.2);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(xx, yy, c, -v*0.2);
		tex.line(el.x2, el.y2, el.x2 + xx, el.y2 + yy);
		dd(x1, y1, c, v*0.5);
		x1 += el.x1;
		y1 += el.y1;

		for (int i = 0; i < 12; i++)
		{
			dd(xx, yy, c + (i + 0) / 12.0 * M_PI, v*0.1);
			dd(x2, y2, c + (i + 1) / 12.0 * M_PI, v*0.1);
			tex.line(x1 + x2, y1 + y2, x1 + xx, y1 + yy);
		}
		dd(x1, y1, c, v*0.7);
		x1 += el.x1;
		y1 += el.y1;

		for (int i = 0; i < 12; i++)
		{
			dd(xx, yy, c + (i + 0) / 12.0 * M_PI, v*0.1);
			dd(x2, y2, c + (i + 1) / 12.0 * M_PI, v*0.1);
			tex.line(x1 + x2, y1 + y2, x1 + xx, y1 + yy);
		}
		dd(x1, y1, c, v*0.3);
		x1 += el.x1;
		y1 += el.y1;

		for (int i = 0; i < 12; i++)
		{
			dd(xx, yy, c + (i + 0) / 12.0 * M_PI, v*0.1);
			dd(x2, y2, c + (i + 1) / 12.0 * M_PI, v*0.1);
			tex.line(x1 + x2, y1 + y2, x1 + xx, y1 + yy);
		}


		tex.numbers((el.x1 + el.x2) / 2, (el.y1 + el.y2) / 2, ftos(el.n));

		break;
	default:
		c = ss(el.x1, el.y1, el.x2, el.y2);
		v = ras(el.x1 - el.x2, el.y1 - el.y2);
		dd(xx, yy, c, v);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(xx, yy, c+0.2, v);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);
		dd(xx, yy, c - 0.2, v);
		tex.line(el.x1, el.y1, el.x1 + xx, el.y1 + yy);

		break;
	}
}

void print_chem(EASY_TEX &tex, const EL_CHAIN &cha)
{
	if (tex.gx() == 0)
	{
		tex.resize(512, 512);
	}

	for (int r = 0; r < tex.gy(); r++)
		for (int i = 0; i < tex.gx(); i++)
		{
			tex.setpixel(i, r, 0, 0, 0);
		}


	vector<GR_EL> el;

	for (int i = 0; i < cha.el.size(); i++)
		if (cha.el[i].p1 != cha.el[i].p2)
		{
			GR_EL elem;

			elem.n = i;
			elem.t = cha.el[i].t;
			elem.p1 = cha.el[i].p1;
			elem.p2 = cha.el[i].p2;

			el.push_back(elem);
		}


	int ssx=5, ssy=2;

	MATRIX<int> m(ssx * 2 - 1, ssy * 2 - 1);

	vector<int> n(el.size());
	vector<vector<vector<int>>> ppos(el.size());
	for (int i = 0; i < el.size(); i++)
	{
		ppos[i].resize(2);
		ppos[i][0].resize(2);
		ppos[i][1].resize(2);
	}

	class H_POS
	{
	public:
		int x, y, t;

		H_POS()
		{
			x = 0;
			y = 0;
			t = 0;
		}
		void next_pos(int sx, int sy)
		{
			if (t == 0)
			{
				x++;
				if (x >= sx)
				{
					x = 0;
					y++;
				}
				if (y >= sy-1)
				{
					x = 0;
					y = 0;
					t = 1;
					return;
				}
			}
			if (t == 1)
			{
				x++;
				if (x >= sx - 1)
				{
					x = 0;
					y++;
				}
				if (y >= sy )
				{
					x = 0;
					y = 1;
					t = 2;
					return;
				}
			}
			if (t == 2)
			{
				x++;
				if (x >= sx)
				{
					x = 0;
					y++;
				}
				if (y >= sy)
				{
					x = 1;
					y = 0;
					t = 3;
					return;
				}
			}
			if (t == 3)
			{
				x++;
				if (x >= sx)
				{
					x = 0;
					y++;
				}
				if (y >= sy)
				{
					x = 0;
					y = 0;
					t = 4;
				}
			}
		}
	};

	vector<H_POS> pos(el.size());
	vector<int> is_see_pos(el.size());

	int exi = 0;
	int curr_t = 0;
	pos[0].x--;

	do
	{
		if (curr_t >= 0 && curr_t < el.size())
		{
			//cout<<m;
			//cout << endl;
			pos[curr_t].next_pos(ssx, ssy);
			bool is_ok_sh = 1;
			
			vector<int> test_p1, test_p2;
			vector<int> test_c_p;
			vector<vector<int>> test_vecs;

			test_vecs.resize(4);
			test_vecs[0] = { 0, 1 };
			test_vecs[1] = { 0, -1 };
			test_vecs[2] = { 1, 0};
			test_vecs[3] = { -1, 0};

			if (pos[curr_t].t == 4)
				is_ok_sh=0;
			if (pos[curr_t].t == 0)
			{
				test_c_p = { 0, 1 };
				test_p1 = { 0, 0 };
				test_p2 = { 0, 1 };
			}
			if (pos[curr_t].t == 1)
			{
				test_c_p = { 1, 0 };
				test_p1 = { 0, 0 };
				test_p2 = { 1, 0 };
			}
			if (pos[curr_t].t == 2)
			{
				test_c_p = { 0, -1 };
				test_p1 = { 0, 0 };
				test_p2 = { 0, -1 };
			}
			if (pos[curr_t].t == 3)
			{
				test_c_p = { -1, 0 };
				test_p1 = { 0, 0 };
				test_p2 = { -1, 0 };
			}

			if (is_ok_sh)
			if (m[pos[curr_t].x * 2 + test_c_p[0]][pos[curr_t].y * 2 + test_c_p[1]] != 0)
				is_ok_sh = 0;

			int b = 0;

			if (is_ok_sh)
			for (int i = 0; i < 4 && b != 1; i++)
			{
				int r = 0;
				b = 0;
				while (b == 0 && r < 2)
				{
					if (pos[curr_t].x * 2 + test_vecs[i][0] * 2 * r >= 0 && pos[curr_t].x * 2 + test_vecs[i][0] * 2 * r < ssx * 2 - 1)
						if (pos[curr_t].y * 2 + test_vecs[i][1] * 2 * r >= 0 && pos[curr_t].y * 2 + test_vecs[i][1] * 2 * r < ssy * 2 - 1)
						{
							int va = m[pos[curr_t].x * 2 + test_vecs[i][0] * 2 * r][pos[curr_t].y * 2 + test_vecs[i][1] * 2 * r];
							if (va == el[curr_t].p1)
								b = 1;
							if (va != el[curr_t].p1)
								if (va != 0)
									b = -1;
						}
					r++;
				}
			}
			if (b != 1)
				if (is_see_pos[el[curr_t].p1])
					is_ok_sh = 0;
			
			
			b = 0;

			if (is_ok_sh)
			for (int i = 0; i < 4 && b != 1; i++)
			{
				int r = 0;
				b = 0;
				while (b == 0 && r < 2)
				{
					if (pos[curr_t].x * 2 + test_p2[0] * 2 + test_vecs[i][0] * 2 * r >= 0 && pos[curr_t].x * 2 + test_p2[0] * 2 + test_vecs[i][0] * 2 * r < ssx * 2 - 1)
						if (pos[curr_t].y * 2 + test_p2[1] * 2 + test_vecs[i][1] * 2 * r >= 0 && pos[curr_t].y * 2 + test_p2[1] * 2 + test_vecs[i][1] * 2 * r < ssy * 2 - 1)
						{
							int va = m[pos[curr_t].x * 2 + test_p2[0] * 2 + test_vecs[i][0] * 2 * r][pos[curr_t].y * 2 + test_p2[1] * 2 + test_vecs[i][1] * 2 * r];
							if (va == el[curr_t].p2)
								b = 1;
							if (va != el[curr_t].p2)
								if (va != 0)
									b = -1;
						}
					r++;
				}
			}

			if (is_ok_sh)
			{
				int va1 = m[pos[curr_t].x * 2 + test_p2[0] * 2][pos[curr_t].y * 2 + test_p2[1] * 2];
				if (va1 != el[curr_t].p2)
					if (va1 != 0)
						is_ok_sh = 0;

				va1 = m[pos[curr_t].x * 2][pos[curr_t].y * 2];
				if (va1 != el[curr_t].p1)
					if (va1 != 0)
						is_ok_sh = 0;
			}


			if (b != 1)
				if (is_see_pos[el[curr_t].p2])
					is_ok_sh = 0;

			if (is_ok_sh)
			{
				is_see_pos[el[curr_t].p1]++;
				is_see_pos[el[curr_t].p2]++;

				m[pos[curr_t].x * 2][pos[curr_t].y * 2] = el[curr_t].p1;
				m[pos[curr_t].x * 2 + test_c_p[0]][pos[curr_t].y * 2 + test_c_p[1]] = 1;
				m[pos[curr_t].x * 2 + test_p2[0] * 2][pos[curr_t].y * 2 + test_p2[1] * 2] = el[curr_t].p2;
				ppos[curr_t][0] = { pos[curr_t].x * 2, pos[curr_t].y * 2 };
				ppos[curr_t][1] = { pos[curr_t].x * 2 + test_p2[0] * 2, pos[curr_t].y * 2 + test_p2[1] * 2 };

				curr_t++;
				if (curr_t < el.size())
				{
					pos[curr_t].t = 0;
					pos[curr_t].x = -1;
					pos[curr_t].y = 0;
				}
			}
			else
			{
				//pos[curr_t].next_pos(ssx,ssy);

				if (pos[curr_t].t == 4)
				{
					if (curr_t <= 0)
						return;
					curr_t--;
					is_see_pos[el[curr_t].p1]--;
					is_see_pos[el[curr_t].p2]--;


					if (pos[curr_t].t == 0)
					{
						test_c_p = { 0, 1 };
						test_p1 = { 0, 0 };
						test_p2 = { 0, 1 };
					}
					if (pos[curr_t].t == 1)
					{
						test_c_p = { 1, 0 };
						test_p1 = { 0, 0 };
						test_p2 = { 1, 0 };
					}
					if (pos[curr_t].t == 2)
					{
						test_c_p = { 0, -1 };
						test_p1 = { 0, 0 };
						test_p2 = { 0, -1 };
					}
					if (pos[curr_t].t == 3)
					{
						test_c_p = { -1, 0 };
						test_p1 = { 0, 0 };
						test_p2 = { -1, 0 };
					}


					int px, py;
					int kol;


					m[pos[curr_t].x * 2 + test_c_p[0]][pos[curr_t].y * 2 + test_c_p[1]] = 0;

					kol = 0;
					px = pos[curr_t].x * 2;
					py = pos[curr_t].y * 2;
					for (int i = 0; i < 4; i++)
						if (test_vecs[i][0] + px >= 0 && test_vecs[i][0] + px < ssx * 2 - 1)
							if (test_vecs[i][1] + py >= 0 && test_vecs[i][1] + py < ssy * 2 - 1)
							{
								if (m[test_vecs[i][0] + px][test_vecs[i][1] + py])
								{
									kol++;
								}
							}
					if (kol == 0)
						m[pos[curr_t].x * 2][pos[curr_t].y * 2] = 0;


					kol = 0;
					px = pos[curr_t].x * 2 + test_p2[0]*2;
					py = pos[curr_t].y * 2 + test_p2[1]*2;
					for (int i = 0; i < 4; i++)
						if (test_vecs[i][0] + px >= 0 && test_vecs[i][0] + px < ssx * 2 - 1)
							if (test_vecs[i][1] + py >= 0 && test_vecs[i][1] + py < ssy * 2 - 1)
							{
								if (m[test_vecs[i][0] + px][test_vecs[i][1] + py])
								{
									kol++;
								}
							}
					if (kol == 0)
						m[pos[curr_t].x * 2 + test_p2[0] * 2][pos[curr_t].y * 2 + test_p2[1] * 2] = 0;


				}


			}


		}
		else
			exi = 1;

	} while (exi==0);

	for (int i = 0; i < el.size(); i++)
	{
		//dd(el[i].x1, el[i].y1, i / (double)el.size()*2*M_PI, 200);
		//dd(el[i].x2, el[i].y2, (i + 0.5) / (double)el.size() * 2 * M_PI, 200);

		vector<int> test_p1, test_p2;
		vector<int> test_c_p;


		if (pos[i].t == 0)
		{
			test_c_p = { 0, 1 };
			test_p1 = { 0, 0 };
			test_p2 = { 0, 1 };
		}
		if (pos[i].t == 1)
		{
			test_c_p = { 1, 0 };
			test_p1 = { 0, 0 };
			test_p2 = { 1, 0 };
		}
		if (pos[i].t == 2)
		{
			test_c_p = { 0, -1 };
			test_p1 = { 0, 0 };
			test_p2 = { 0, -1 };
		}
		if (pos[i].t == 3)
		{
			test_c_p = { -1, 0 };
			test_p1 = { 0, 0 };
			test_p2 = { -1, 0 };
		}

		el[i].x1 = 100 + 100 * (pos[i].x);
		el[i].y1 = tex.gy()-100 - 100 * (pos[i].y);
		el[i].x2 = 100 + 100 * (pos[i].x + test_c_p[0]);
		el[i].y2 = tex.gy() - 100 - 100 * (pos[i].y + test_c_p[1]);
	}

	for (int i = 0; i < el.size(); i++)
	{
		draw_el(tex, el[i]);
	}

	for (int i = 0; i < el.size(); i++)
	{
		for (int r = 0; r < el.size(); r++)
		{
			if (el[i].p1 == el[r].p1)
			{
				tex.line(el[i].x1, el[i].y1, el[r].x1, el[r].y1);
			}
			if (el[i].p1 == el[r].p2)
			{
				tex.line(el[i].x1, el[i].y1, el[r].x2, el[r].y2);
			}
			if (el[i].p2 == el[r].p1)
			{
				tex.line(el[i].x2, el[i].y2, el[r].x1, el[r].y1);
			}
			if (el[i].p2 == el[r].p2)
			{
				tex.line(el[i].x2, el[i].y2, el[r].x2, el[r].y2);
			}
		}
	}

}

/*
class person
{
public:
	person(){};
	explicit person(const string& name):name(name){}
	void set_name(const string& n){ name = n; }
	string get_name()const{ return name; };
	void info()const
	{
		cout << "person name:" << name << endl;
	}
private:
	string name;
};

class student :public person
{
public:
	student(const string &name, const string &passed) :person(name), passed(passed){}
	void info()const
	{
		cout << "student name:" << get_name() << endl;
		cout << "passet things:" << passed << endl;
	}
private:
	string passed;
};

class teacher :public person
{
public:
	void addstudent(person &person);



};

template<class T>
class mylist
{
public:
	bool insert(const T&, size_t index);
	T access(size_t index)const;
private:
	T*buf;
	size_t bufsize_;
};

template<class T>
class myset:
{
public:
	bool insert(const T&, size_t index);
	T access(size_t index)const;
private:
	T*buf;
	size_t bufsize_;
};
*/

int main()
{
	setlocale(0, "RU");


	cout << "Грузись питон"<<endl;
	sympy_init();

	{
		int inn;
		cout << "Сосчитать какуюнить резестивную цепь с одним источником?(0-нет 1-да):";
		cin >> inn;
		if (inn)
		{
			cout << "EX: 1 4 1 I 0  2 3 4 R 0.5  3 1 4 R 2  4 1 2 L 2  5 2 3 R 1  6 3 4 C 4" << endl;

			cout << "вводите:";
			char buf[1000];
			gets_s(buf, 1000);
			gets_s(buf,1000);

			string s = buf;

			//EL_CHAIN cha("1 1 5 U 0  2 1 2 R 0.0625  3 2 5 R 0.25  4 2 3 L 0.025  5 3 5 C 0.4  6 3 4 R 0.25  7 4 5 R 1");//from MU

			EL_CHAIN cha(s);//from MU

			cha.comp_par_1_iu_uns();
			getch();

		}

	}


	thread th(task_graphix);
	th.detach();

	//comp_4((string)"10/s*(1-2*exp(-10*s)+exp(-20*s))");
	//comp_4((string)"2.0/(4.0*s**2 + 8.0*s + 3.5)");
	//th.join();
	//system("pause");
	//sympy_sim("1+3");
	//cout << endl << "Грузись питон";
	//call_python_sympy("sin(t)+1");
	//puts("L: sin(t)");
	//puts(sympy_lap("sin(t)").data());
	//puts("S: (100*s-s**2)*(100+30*s-10)");
	//puts(sympy_sim("(100*s-s**2)*(100+30*s-10)").data());

	KU_INP input;


	//comp_4(
	//	(string)"2.0/(4.0*s**2 + 8.0*s + 3.5)",
	//	"0+(Heaviside((t-(0)))*(10))+(-Heaviside((t-(10))) * 2 * (10))+(Heaviside((t-(20)))*(10))",
	//	"0+(10/s)*exp(-s*(0))+(-20/s)*exp(-s*(10))+(10/s)*exp(-s*(20))",
	//	40);

	/** /   //что это? тестовый вариант наверно
	input.cha_str = "1 4 1 I 0  2 3 4 R 0.5  3 1 4 R 2  4 1 2 L 2  5 2 3 R 1  6 3 4 C 4";
	input.el_id = 2;
	input.t_sign = SIGN_T;
	input.ts = 20;
	input.as = 10;
	input.T = 40;
	/**/

	/** /
	input.cha_str = "1 1 4 U 0  2 2 4 R 1  3 1 3 R 1  4 3 4 R 1  5 1 3 L 1  6 3 2 L 0.25";//14 var
	input.el_id = 2;
	input.t_sign = SIGN_V;
	input.ts = 5;
	input.as = 3;
	input.T = 10;
	/**/

	int va = 0;

	cout << "Доступны варианты 0 2 7 13 14 " << endl;
	cout << "Введите номер варианта:" << endl;
	cin >> va;
	
	
	if (va == 0)
	{
		input.cha_str = "1 4 1 I 0  2 3 4 R 0.5  3 1 4 R 2  4 1 2 L 2  5 2 3 R 1  6 3 4 C 4";
		input.el_id = 2;
		input.t_sign = SIGN_T;
		input.ts = 20;
		input.as = 10;
		input.T = 40;
	}
	if (va == 2)
	{
		input.cha_str = "1 2 1 U 0  2 2 3 C 4  3 3 1 R 1  4 3 4 C 1  5 4 1 R 0.5  6 1 4 R 1";//2 v
		input.el_id = 6;
		input.t_sign = SIGN_D;
		input.ts = 0.1;
		input.as = 20;
		input.T = input.ts * 2;
	}
	if (va == 7)
	{
		input.cha_str = "1 1 2 I 0  2 2 1 R 2  3 2 3 L 1  4 3 1 L 4  5 3 4 R 2  6 1 4 R 1";//13 v
		input.el_id = 6;
		input.t_sign = SIGN_D;
		input.ts = 2;
		input.as = 5;
		input.T = input.ts * 2;
	}
	if (va == 13)
	{
		input.cha_str = "1 1 2 I 0  2 2 1 R 1  3 2 3 C 4  4 3 1 R 2  5 3 4 C 1  6 1 4 R 1";//13 v
		input.el_id = 6;
		input.t_sign = SIGN_D;
		input.ts = 1;
		input.as = 5;
		input.T = input.ts * 2;
	}
	if (va == 14)
	{
		input.cha_str = "1 2 1 U 0  2 2 3 R 1  3 3 1 L 4  4 3 4 R 0.5  5 4 1 L 1  6 4 1 R 1";//14 v
		input.el_id = 6;
		input.t_sign = SIGN_V;
		input.ts = 0.5;
		input.as = 20;
		input.T = input.ts * 2;
	}

	int inn;

	cout << "U - ит" << endl;
	cout << "I - ин" << endl;
	cout << "R - резистор" << endl;
	cout << "C - C-элемент" << endl;
	cout << "L - L-элемент" << endl;
	cout << "K - кз" << endl;
	cout << "H - хх" << endl;
	cout << "сейчас строка цепи:" << endl;
	cout << input.cha_str << endl;
	cout << "Сменить?(0-нет 1-да):";
	cin >> inn;
	if (inn)
	{
		cout << "вводите:";
		char buf[1000];
		gets_s(buf);
		gets_s(buf);
		input.cha_str = buf;

	}
	cout << "сейчас номер элемента реакции:" << endl;
	cout << input.el_id << endl;
	cout << "Сменить?(0-нет 1-да):";
	cin >> inn;
	if (inn)
	{
		cout << "вводите:";
		cin >> input.el_id;
	}
	cout << "A 1-недоступно"<<endl;
	cout << "Б 2-недоступно" << endl;
	cout << "В 3 пирамида" << endl;
	cout << "Г 4-недоступно" << endl;
	cout << "Д 5 анти ступени" << endl;

	cout << "сейчас тип сигнала:" << endl;
	cout << input.t_sign << endl;
	cout << "Сменить?(0-нет 1-да):";
	cin >> inn;
	if (inn)
	{
		cout << "вводите:";
		cin >> input.t_sign;
	}
	cout << "t сигнала:" << endl;
	cout << input.ts << endl;
	cout << "Сменить?(0-нет 1-да):";
	cin >> inn;
	if (inn)
	{
		cout << "вводите:";
		cin >> input.ts;
	}
	cout << "А сигнала:" << endl;
	cout << input.as << endl;
	cout << "Сменить?(0-нет 1-да):";
	cin >> inn;
	if (inn)
	{
		cout << "вводите:";
		cin >> input.as;
	}
	cout << "период Т:" << endl;
	cout << input.T << endl;
	cout << "Сменить?(0-нет 1-да):";
	cin >> inn;
	if (inn)
	{
		cout << "вводите:";
		cin >> input.T;
	}


	/** /
	input.cha_str = "1 1 4 U 0  2 0 0 R 1  3 1 2 R 1  4 2 3 R 1  5 2 4 C 1  6 3 4 C 0.25  7 3 4 R 1";//
	input.el_id = 7;
	input.t_sign = SIGN_Z;
	input.ts = 1;
	input.as = 1;
	input.T = input.ts * 2;
	/**/

	//

	/*
	for (auto i = A.begin_row(); i != A.end_row(); i++)
	{
	for (auto r = (*i).; r != res.end(); r++)
	auto sch = *i;
	printf("r:%lf i:%lf \n", sch.real(), sch.imag());
	}
	}
	*/





	//EL_CHAIN cha_t("1 1 3 R 4  2 1 2 R 2  3 3 2 U 8  4 1 4 R 2  5 2 5 R 4  6 4 5 R 2  7 3 5 K 0");I(U)=4 Rэ=2
	//EL_CHAIN cha_t("1 1 5 R 4000  2 5 6 R 10000  3 2 6 R 11000  4 3 7 U 6  5 4 8 R 2000  6 2 1 K 0  7 3 2 K 0  8 3 4 K 0  9 6 7 K 0  10 8 7 K 0");
	//EL_CHAIN cha_t("1 1 3 I 0.016 2 1 2 R 15000  3 2 3 R 9000  4 2 3 R 6000  5 2 3 R 17000 ");//idz3
	//cha_t.comp_par_1_iu_uns();
	//EL_CHAIN cha("1 1 4 U 0  2 2 4 R 1  3 1 3 R 1  4 3 2 R 1  6 3 4 C 1  5 4 2 C 0.25"); // v1 lin(62)
	//EL_CHAIN cha("1 1 3 I 0.018  2 1 3 R 20000  3 3 4 R 8000  4 1 4 R 14000  5 4 5 R 7000  6 1 5 R 4000  7 1 2 L 0.004  8 5 2 C 0.000014 ");

	//Сначала ввод входных данных

	//POLY p1, p2;
	//	f = "(s**2 + 33*s + 238)/((1*s**3 +30*s**2 + 300*s + 1000))";
	//p1 = get_pol_coef("s**2 + 33*s + 238");
	//p2 = get_pol_coef("1*s**3 + 30*s**2 + 300*s + 1000");
	//NEOPR_COEF d(p1,p2);


	//EL_CHAIN cha("1 1 5 U 0  2 1 2 R 0.0625  3 2 5 R 0.25  4 2 3 L 0.025  5 3 5 C 0.4  6 3 4 R 0.25  7 4 5 R 1");//from MU

	EL_CHAIN cha(input.cha_str);//from MU
	


	EASY_TEX tex;


	tex.resize(512, 256);
	print_chem(tex, cha);
	mtx.lock();
	img.add(tex);
	mtx.unlock();
	//EL_CHAIN cha("1 1 4 U 0  2 2 4 R 1  3 1 3 R 1  4 3 4 R 1  5 1 3 L 1  6 3 2 L 0.25");//from 14V

	//подсчёт уравнений состояния 
	cha.comp_ur_so();

	//h1 аналитически
	cha.comp_h1(input.el_id);

	/*
	double t_3;
	t_3 = -1 / cha.per_s[0].lamd[0].real();
	for (int i = 0; i < cha.per_s[0].lamd.size(); i++)
		t_3 = max(-1 / t_3, -1 / cha.per_s[1].lamd[0].real());

	t_3 *= 3;



	//h1 по лаплассовски
	
	//(238 + 33 s + s^2)/(1000 + 300 s + 30 s^2 + s^3)
	//cha.h1_2.str_l = ((L_S)("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))") );
	//cha.h1_2.sym_str = ((L_S)("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))")).s;
	//sympy_lap("Heaviside(t)*(10)");
	*/
	//H1  и какойто сигнал
	auto res2 = comp_2(cha, input.el_id, input.t_sign, input.as, input.ts);
	///вход2

	/*
	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*t_3 / vx.size();
		vy[i] = atof1(sympy_eva(cha.h1_1.sym_str, "t", ftos(vx[i])));
	}

	mtx.lock();
	tex = create_plot(512, 256, vx, vy, 0, 0, 0, 0);
	img.add(tex);
	mtx.unlock();



	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*t_3 / vx.size();
		vy[i] = atof1(sympy_eva(res2.h1_t, "t", ftos(vx[i])));
		//cout <<vx[i]<<"  "<< sympy_eva(cha.f1_t, "t", ftos(vx[i]))<<endl;
	}

	mtx.lock();
	tex = create_plot(512, 256, vx, vy, 0, 0, 0, 0);
	img.add(tex);
	mtx.unlock();
*/

	auto res3 = comp_3(res2.h1_s*(L_S)"s", res2.f1_t, res2.f1_s, input.t_sign, input.ts);

	auto res4 = comp_4(res2.h1_s*(L_S)"s", res2.f1_t, res2.f1_s, input.T);



	print_ku(cha, res2, res3, res4);


	//EL_CHAIN_L cha1("1 4 1 I 1  2 3 4 R 0.5  3 1 4 R 2  4 1 2 R s*2  5 2 3 R 1  6 3 4 R 1/(s*4)");//from MU

	//cha1.comp_par_1_iu_uns();
	//cha1.comp_h1(2);
	//cha.comp_ur_so();
	//cha.comp_h1(0);

	cout << "press Enter";
	
	char ch = 0;
	do
	{
		if (kbhit())
			ch = getch();
	} while (ch != 13);

	return 0;
}



void print_ku(const EL_CHAIN &cha, const COMP_2_RES &res2, const COMP_3_RES &res3, const COMP_4_RES &res4)
{
	system("cls");
	cout << "1. Анализ цепи во временной области." << endl;
	cout << "1.1. Составить уравнения состояния цепи для t >= 0. " << endl;
	cout << "Для составления уравнений состояния следует заменить L->ИТ C->ИН и выразить через них U/I для заменённых элементов" << endl;
	cout << "Получившиеся I(t)/U(t) надо разделить на C/L чтобы получить U'(t)/I'(t)" << endl;
	cout << "Результат можно записать в матричном виде:" << endl;
	cout << "[f'пс(t)]=[A][fпс(t)]+[B][f1(t)]" << endl;
	cout << "[A]:" << endl;
	cout << cha.ur_so.A << endl;
	cout << "[B]:" << endl;
	cout << cha.ur_so.B << endl;


	M_S f_str;


	f_str.resize(cha.ur_so.id_cl.size(), 1);
	for (int i = 0; i < cha.ur_so.id_cl.size(); i++)
	{
		if (cha.el[cha.ur_so.id_cl[i]].t == EL_C)
		{
			f_str[i][0] = "C";
		}
		if (cha.el[cha.ur_so.id_cl[i]].t == EL_L)
		{
			f_str[i][0] = "L";
		}
		f_str[i][0] = f_str[i][0] + ftos(cha.ur_so.id_cl[i]) + "'(t)";
	}



	cout << "[f'пс(t)]:" << endl;
	cout << f_str << endl;

	f_str.resize(cha.ur_so.id_cl.size(), 1);
	for (int i = 0; i < cha.ur_so.id_cl.size(); i++)
	{
		if (cha.el[cha.ur_so.id_cl[i]].t == EL_C)
		{
			f_str[i][0] = "C";
		}
		if (cha.el[cha.ur_so.id_cl[i]].t == EL_L)
		{
			f_str[i][0] = "L";
		}
		f_str[i][0] = f_str[i][0] + ftos(cha.ur_so.id_cl[i]) + "(t)";
	}

	cout << "[fпс(t)]:" << endl;
	cout << f_str << endl;


	f_str.resize(cha.ur_so.id_ui.size(), 1);
	for (int i = 0; i < cha.ur_so.id_ui.size(); i++)
	{
		if (cha.el[cha.ur_so.id_ui[i]].t == EL_U)
		{
			f_str[i][0] = "U";
		}
		if (cha.el[cha.ur_so.id_ui[i]].t == EL_I)
		{
			f_str[i][0] = "I";
		}
		f_str[i][0] = f_str[i][0] + ftos(cha.ur_so.id_ui[i]) + "(t)";
	}

	cout << "[f1(t)]:" << endl;
	cout << f_str << endl;


	cout << "1.2. По уравнениям состояния аналитическим расчетом во временной области найти переходную характеристику h1(t) для реакции и построить ее график." << endl;

	cout << "Для начала находим собственные числа матрицы A" << endl;
	for (int i = 0; i < cha.per_s[0].lamd.size(); i++)
		cout << "l" + ftos(i + 1) << "=" << comp_to_s(cha.per_s[0].lamd[i]) << endl;

	cout << "Теперь мы знаем что реакция на б(t) имеет такой вид:" << endl;
	cout << cha.obs_vid << endl;

	cout << "Соответственно зная что переходные процессы затухают можно из уравнения состояния можно вычислить fпс(inf):" << endl;
	cout << "ННУ дают возможность также вычислить fпс(0) f'пс(0)" << endl;
	cout << "и составить такую систему линейных уравнений:" << endl;
	cout << "A*fl1(0)+B*fl2(0)=fпс(0)-fпс(inf)" << endl;
	cout << "A*fl1'(0)+B*fl2'(0)=fпс'(0)" << endl;

	cout << "Соответственно решив данную систему для всех переменных состояния получаем:" << endl;
	for (int i = 0; i < cha.per_s.size(); i++)
		cout << cha.per_s[i].cool_str << endl;

	cout << "Выразим требуемое I через переменные состояния и источник" << endl;
	cout << ("h1(t)=I(t)=");


	for (int i = 0; i < cha.ur_so.id_ui.size(); i++)
	{
		string temp;
		if (cha.el[cha.ur_so.id_ui[i]].t == EL_U)
		{
			temp = "U";
		}
		if (cha.el[cha.ur_so.id_ui[i]].t == EL_I)
		{
			temp = "I";
		}
		temp = temp + ftos(cha.ur_so.id_ui[i]) + "(t)*(" + ftos(cha.h1_1.k_u) + ")";
		cout << temp;
	}
	for (int i = 0; i < cha.h1_1.k_per_s.size(); i++)
	{
		string temp;
		if (cha.el[cha.ur_so.id_cl[i]].t == EL_C)
		{
			temp = "U";
		}
		if (cha.el[cha.ur_so.id_cl[i]].t == EL_L)
		{
			temp = "I";
		}
		temp = "+" + temp + ftos(cha.ur_so.id_cl[i]) + "(t)*(" + ftos(cha.h1_1.k_per_s[i]) + ")";
		cout << temp;
	}
	cout << endl;
	cout << ("h1(t)=") << "(" << cha.h1_1.sym_str << ")*b1(t)";
	cout << endl;


	cout << "2. Анализ цепи операторным методом при действии одиночного импульса на входе." << endl;
	cout << "2.1. В соответствии с номером выполняемого варианта определить функцию передачи напряжений или токов. Осуществить проверку функции передачи при s = 0 и ; s -> inf представить соответствующие этим значениям схемы замещения цепи. " << endl;
	cout << "H1I(S):" << endl;
	cout << res2.h1_s << endl;

	cout << "H1I(0):" << endl;
	cout << res2.s0 << endl;

	cout << "H1I(0) по схеме:" << endl;
	cout << res2.s0 << endl;

	cout << "H1I(inf):" << endl;
	cout << res2.s9s << endl;

	cout << "H1I(inf) по схеме:" << endl;
	cout << res2.s9 << endl;


	cout << "2.2. Найти нули и полюсы функции передачи и показать их расположение на плоскости комплексной частоты. По значениям полюсов функции передачи дать заключение о характере и практической длительности переходного процесса. " << endl;

	cout << "POL:";
	for (int i = 0; i < res2.v_polus.size(); i++)
		cout << comp_to_s(res2.v_polus[i]) << "  ";
	cout << endl;

	cout << "ZER:";
	for (int i = 0; i < res2.v_zero.size(); i++)
		cout << comp_to_s(res2.v_zero[i]);
	cout << endl;

	cout << "2.3. Определить переходную h1(t) характеристику цепи, сравнить с найденной в п. 1.2 задания. Проверить h1(0) и h1(inf) по аналитическому выражению h1(t) и непосредственно по схеме цепи." << endl;
	cout << res2.h1_t << endl;

	cout << "h1I(0):" << endl;
	cout << res2.h1_0s << endl;

	cout << "h1I(0) по схеме:" << endl;
	cout << res2.h1_0 << endl;

	cout << "h1I(inf):" << endl;
	cout << res2.h1_9s << endl;

	cout << "h1I(inf) по схеме:" << endl;
	cout << res2.h1_9 << endl;

	cout << "2.4. Определить изображение по Лапласу входного одиночного импульса." << endl;
	cout << "F1(S):" << endl;
	cout << res2.f1_s << endl;
	cout << "f1(t):" << endl;
	cout << res2.f1_t << endl;

	cout << "2.5. Определить изображение выходного сигнала и далее найти реакцию I2(t) или U2(t) во временной области. Построить графики входного и выходного сигналов на одном рисунке." << endl;
	cout << "F2(S):" << endl;
	cout << res2.f2_s << endl;
	cout << "f2(t):" << endl;
	cout << res2.f2_t << endl;



	double t_3;

	t_3 = -1 / cha.per_s[0].lamd[0].real();
	for (int i = 0; i < cha.per_s[0].lamd.size(); i++)
		t_3 = max(-1 / t_3, -1 / cha.per_s[1].lamd[0].real());

	t_3 *= 3;



	//h1 по лаплассовски

	//(238 + 33 s + s^2)/(1000 + 300 s + 30 s^2 + s^3)
	//cha.h1_2.str_l = ((L_S)("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))") );
	//cha.h1_2.sym_str = ((L_S)("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))")).s;
	//sympy_lap("Heaviside(t)*(10)");
	//H1  и какойто сигнал
	///вход2

	vector<double> vx, vy, vz;
	EASY_TEX tex;


	vx.resize(200);
	vy.resize(vx.size());
	vz.resize(vx.size());

	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*t_3 / vx.size();
		vy[i] = atof1(sympy_eva(cha.h1_1.sym_str, "t", ftos(vx[i])));
	}


	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "аш1 от тэ через дифур");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();



	for (int i = 0; i < vx.size(); i++)
	{
		vx[i] = i*t_3 / vx.size();
		vy[i] = atof1(sympy_eva(res2.h1_t, "t", ftos(vx[i])));
		//cout <<vx[i]<<"  "<< sympy_eva(cha.f1_t, "t", ftos(vx[i]))<<endl;
	}

	tex = create_plot(1024, 512, vx, vy, 0, 0, 0, 0);
	tex.numbers(100, tex.gy() - 20, "аш1 от тэ через лапласса");
	tex.numbers(100, tex.gy() - 20 - 18, "");
	mtx.lock();
	img.add(tex);
	mtx.unlock();
 





	cout << "3.1. Используя найденное в 2.1 выражение H(s)," << endl;
	cout << "вычислить и построить графики АЧХ и ФЧХ функций передачи цепи H(jw)." << endl;
	cout << "Произвести проверку АЧХ при w = 0 и w -> inf." << endl;

	cout << "АЧХ:" << res3.ACH_H;
	cout << "ФЧХ:" << res3.FCH_H;
	cout << "АЧХ(0):" << res3.ACH_H_0;
	cout << "АЧХ(0) по схеме:" << res2.h1_0;
	cout << "АЧХ(inf):" << res3.ACH_H_inf;
	cout << "АЧХ(inf) по схеме:" << res2.h1_9;

	cout << "3.2. Определить полосу пропускания цепи по уровню 0,707 max|H(jw)|." << endl;
	cout << "max(A(w))" << res3.max_ACH << endl;
	res3.polosa.show();

	cout << "3.3. Найти и построить амплитудный и фазовый спектры входного" << endl;
	cout << "одиночного импульса.Найти ширину амплитудного спектра по 0,1 max(F(jw))." << endl;
	cout << "A(w)=" << res3.AS_F1<<endl;
	cout << "F(w)=" << res3.FS_F1 << endl;
	cout << "dw=" << res3.dw_AS_F1 << endl;

	cout << "3.4. Сопоставить спектры входного импульса с частотными характеристиками цепи." << endl;
	cout << "Дать предварительное заключение об ожидаемых искажениях сигнала на выходе цепи." << endl;
	cout << "Сравнить эти качественные оценки с сигналом на выходе, полученным в п. 2.5 задания." << endl;
	cout << res3.p34;

	cout << "3.5. Построить графики амплитудного и фазового спектров выходного сигнала, используя" << endl;
	cout << "графики пп. 3.1, 3.3 задания. Проконтролировать площадь реакции по значению ее спектра при w = 0." << endl;
	cout << "Значение спектра при w = 0:" << res3.AS_F2_0 << endl;





	cout << "4. Анализ цепи частотным методом при периодическом воздействии. На вход цепи подан" << endl;
	cout << "сигнал в виде периодической последовательности импульсов напряжения или тока.";

	cout << "4.1.Разложить в ряд Фурье заданный входной периодический сигнал.Построить его" << endl;
	cout << "амплитудный и фазовый дискретные спектры." << endl;
	cout << "4.2. Построить на одном графике заданный входной периодический сигнал и его аппроксимацию" << endl;
	cout << "отрезком ряда Фурье. Число гармоник отрезка ряда Фурье определяется по уровню 0,1 Akm" << endl;
	cout << ", где Akm – максимальная составляющая амплитудного спектра.";
	cout << "Akm:" << res4.Ak << endl;
	cout << "Выбрано гармоник:" << res4.lev << endl;
	cout << "4.3. Используя рассчитанные в п. 3.1 задания АЧХ и ФЧХ, найти реакцию цепи в виде отрезка" << endl;
	cout << "ряда Фурье с числом гармоник, определенным для входного сигнала.";
	cout << "4.4. Построить амплитудный и фазовый дискретные спектры выходного сигнала. Построить" << endl;
	cout << "график выходного сигнала, найденного в п. 4.3 задания, в одном масштабе рядом" << endl;
	cout << "с графиком аппроксимированного входного сигнала." << endl;
	cout << "4.5. Дать заключение об искажении сигнала на выходе цепи. " << endl;
	cout << "";
}





void outche(EL_CHAIN& cha)
{
	EASY_TEX tex;


	tex.resize(512, 256);
	print_chem(tex, cha);
	mtx.lock();
	img.add(tex);
	mtx.unlock();
}
