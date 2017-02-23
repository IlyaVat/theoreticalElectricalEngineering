#undef _DEBUG
#include <python.h>
#pragma comment( lib, "python36.lib" )

#include <stdlib.h>

#include <cmath>

#include "el_cha.h"
#include "py_sympy.h"
#include "usefull.h"


#define _DEBUG

#define EL_R 2
#define EL_I 3
#define EL_U 4
#define EL_L 5
#define EL_C 6

#define SIGN_T 0
#define SIGN_V 0

//using namespace itl;
using namespace std;

string sympy_sim(string funct);
string sympy_lap(string funct);
string sympy_eva(string f1, string f2, string f3);
string sympy_dif(string f1, string f2);




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

class IMAGE_PNG
{
public:
	vector<unsigned char> buff;
	int x, y;
	IMAGE_PNG()
	{
		x = 0;
		y = 0;

	}
};

class EL_CHAIN
{
public:
	vector<POI> po;
	vector<ELEM> el;
	
	IMAGE_PNG img;

	UR_SOST ur_so;
	vector<PER_S> per_s;
	string obs_vid;
	UR_H1_1 h1_1;
	UR_H1_2 h1_2;
	string f2_s, f2_t;

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

	void comp_signl(int t_s, double im, double ti)
	{
		/*
		string b_s,b_l;
		if (t_s == SIGN_T)
		{
			b_s.resize(3);
			b_l.resize(3);
			b_s = "Heaviside(t)*(" + to_string(im) + ")-Heaviside(t - (" + to_string(ti) + ") / 2) * 2 * (" + to_string(im) + ")+Heaviside(t - (" + to_string(ti) + "))*(" + to_string(im) + ")";
			//cout << "\n" << b_s;

			b_l = sympy_lap(b_s);
			//cout << "\n" << b_l;
		}

		//ну теперь перемножим...

		L_S res_l = sympy_sim((b_l * h1_2.str_l*(L_S)"s").s);
		string res_s;
		

		res_s = sympy_sim(sympy_alap((b_l * h1_2.str_l*(L_S)"s").s));



		cout << "\n" << res_l;
		cout << "\n" << res_s;


		//ну результат то получился но он не очень

		*/

		vector<string> a_t,a_l,a_res;
		vector<double> t_sdv;

		string res;

		if (t_s == SIGN_T)
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
			a_t[2] = "+Heaviside(t)*(" + ftos(im) + ")";
			t_sdv[2] = ti;
			//cout << "\n" << b_s;

			res += "Исходный сигнал:\n";
			res += "b1(t)*(" + ftos(im) + ")";
			res += "-b1(t-" + ftos(ti / 2) + ") * 2 * (" + ftos(im) + ")";
			res += "+b1(t-" + ftos(ti) + ")*(" + ftos(im) + ")";


			//b_l = sympy_lap(b_s);
			//cout << "\n" << b_l;
		}

		if (t_s == SIGN_V)
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
			a_t[2] = "+Heaviside(t)*t*(" + ftos(im / ti * 2) + ")";
			t_sdv[2] = ti;
			//cout << "\n" << b_s;

			res += "Исходный сигнал:\n";
			res += "b1(t)*(" + ftos(im) + ")";
			res += "-b1(t-" + ftos(ti / 2) + ") * 2 * (" + ftos(im) + ")";
			res += "+b1(t-" + ftos(ti) + ")*(" + ftos(im) + ")";


			//b_l = sympy_lap(b_s);
			//cout << "\n" << b_l;
		}

		res += "\nСигнал в операторной области:\n";

		for (int i = 0; i < a_t.size(); i++)
		{
			a_l[i] = sympy_lap(a_t[i]);
			string ts;
			if (t_sdv[i] != 0)
				ts = "(" + a_l[i] + ")*exp(-s*(" + ftos(t_sdv[i]) + "))";
			else
				ts = "(" + a_l[i] + ")";

			if (i > 0)
				res += "+";
			res += ts;

		}

		L_S H;
		res += "\nH(S)=H1(S)*S:\n";
		H = h1_2.str_l*(L_S)"s";
		res += H.s;


		res += "\nF2(S)=H(S)*F1(S):\n";
		f2_s = "";
		for (int i = 0; i < a_t.size(); i++)
		{
			L_S temp;
			temp = (a_l[i] * H).s;

			//Формула сдвига лапласса
			//f(t-a)=exp(-a*s)*F(s)

			if (i > 0)
				f2_s += "+";


			if (t_sdv[i]!=0)
				f2_s += "(" + temp.s + ")*exp(-s*(" + ftos(t_sdv[i]) + "))";
			else
				f2_s += "(" + temp.s + ")";
			//
			string rr = crazy_alap(temp.s);
			if (rr == "")
				rr = sympy_alap(temp.s);
			a_res[i] = rr;
		}

		res += f2_s;

		res += "\nf2(t):\n";

		f2_t = "";

		for (int i = 0; i < a_t.size(); i++)
		{
			string temp=a_res[i];
			myreplace(temp, "exp", "hidden_exp");//этот симпи не должен изгаживать экспоненты
			temp = sympy_sim(temp);
			myreplace(temp, "hidden_exp", "exp");

			//тэшку только в конце переписываем а то скобки раскроются
			if (t_sdv[i] != 0)
				myreplace(temp, "t", "(t-" + ftos(t_sdv[i]) + ")");

			myreplace(temp, "((t-" + ftos(t_sdv[i]) + "))", "(t-" + ftos(t_sdv[i]) + ")");

			myreplace(temp, "Heaviside", "b1");

			if (i > 0)
				f2_t += "+";
			f2_t += "(" + temp + ")";


			myreplace(a_res[i], "t", "(t-" + ftos(t_sdv[i]) + ")"); 
		}

		res+=f2_t;


		cout << res;


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
				int id = ur_so.id_cl[i];//unused

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

	void comp_h1_l(int id_res)
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

		h1_2.str_l = cha1.h1_l;
		h1_2.sym_str = cha1.h1_sympy;
		h1_2.v_polus = cha1.v_polus;
		h1_2.v_zero = cha1.v_zero;

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
				double i_u = 0;
				int unk_el = -1;
				int count = 0;
				for (int r = 0; r < el.size(); r++)
					if (el[r].p1 != el[r].p2)
					{
						if (el[r].p1 == i)
						{

							if (el[r].t == 0 && comp[r] == 0 || el[r].t == 4)
								if (comp[r] == 0)
								{
									unk_el = r;
									count++;
								}

							if (el[r].t == 2 || el[r].t == 0 && comp[r] == 1)
							{
								i_u += el[r].I;
							}

						}
						if (el[r].p2 == i)
						{

							if (el[r].t == 0 && comp[r] == 0 || el[r].t == 4)
							{
								unk_el = r;
								count++;
							}

							if (el[r].t == 2 || el[r].t == 3 || el[r].t == 0 && comp[r] == 1)
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

static struct PyModuleDef spammodule = {
	PyModuleDef_HEAD_INIT,
	"spam",   /* name of module */
	NULL, /* module documentation, may be NULL */
	-1,       /* size of per-interpreter state of the module,
			  or -1 if the module keeps state in global variables. */
			  NULL
};
static PyObject *SpamError;

static PyObject *
spam_system(PyObject *self, PyObject *args)
{
	const char *command;
	int sts;

	if (!PyArg_ParseTuple(args, "s", &command))
		return NULL;
	sts = system(command);
	if (sts < 0) {
		PyErr_SetString(SpamError, "System command failed");
		return NULL;
	}
	return PyLong_FromLong(sts);
}

PyMODINIT_FUNC
PyInit_spam(void)
{
	PyObject *m;

	m = PyModule_Create(&spammodule);
	if (m == NULL)
		return NULL;

	SpamError = PyErr_NewException("spam.error", NULL, NULL);
	Py_INCREF(SpamError);
	PyModule_AddObject(m, "error", SpamError);
	return m;
}

static PyObject *my_callback = NULL;

static PyObject *
my_set_callback(PyObject *dummy, PyObject *args)
{
	PyObject *result = NULL;
	PyObject *temp;

	if (PyArg_ParseTuple(args, "O:set_callback", &temp)) {
		if (!PyCallable_Check(temp)) {
			PyErr_SetString(PyExc_TypeError, "parameter must be callable");
			return NULL;
		}
		Py_XINCREF(temp);         /* Add a reference to new callback */
		Py_XDECREF(my_callback);  /* Dispose of previous callback */
		my_callback = temp;       /* Remember new callback */
		/* Boilerplate to return "None" */
		Py_INCREF(Py_None);
		result = Py_None;
	}
	return result;
}

void print_ku(const EL_CHAIN &cha);

int main()
{
	setlocale(0, "RU");
	sympy_init();


	//call_python_sympy("sin(t)+1");
	//puts("L: sin(t)");
	//puts(sympy_lap("sin(t)").data());
	//puts("S: (100*s-s**2)*(100+30*s-10)");
	//puts(sympy_sim("(100*s-s**2)*(100+30*s-10)").data());

	MATRIX<double> mm(4, 3);
	mm[1][0] = 1;
	cout << mm;
	system("pause");



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

	get_pol_coef("-s**2 - 33*s - 238");
	//POLY p1, p2;
	//	f = "(s**2 + 33*s + 238)/((1*s**3 +30*s**2 + 300*s + 1000))";
	//p1 = get_pol_coef("s**2 + 33*s + 238");
	//p2 = get_pol_coef("1*s**3 + 30*s**2 + 300*s + 1000");
	//NEOPR_COEF d(p1,p2);


	//EL_CHAIN cha("1 1 5 U 0  2 1 2 R 0.0625  3 2 5 R 0.25  4 2 3 L 0.025  5 3 5 C 0.4  6 3 4 R 0.25  7 4 5 R 1");//from MU

	//EL_CHAIN cha("1 4 1 I 0  2 3 4 R 0.5  3 1 4 R 2  4 1 2 L 2  5 2 3 R 1  6 3 4 C 4");//from MU
	EL_CHAIN cha("1 1 4 U 0  2 2 4 R 1  3 1 3 R 1  4 3 4 R 1  5 1 3 L 1  6 3 2 L 0.25");//from 14V

	//подсчёт уравнений состояния 
	cha.comp_ur_so();

	//h1 аналитически
	cha.comp_h1(2);

	//h1 по лаплассовски
	cha.comp_h1_l(2);

	//(238 + 33 s + s^2)/(1000 + 300 s + 30 s^2 + s^3)
	//cha.h1_2.str_l = ((L_S)("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))") );
	//cha.h1_2.sym_str = ((L_S)("(409600.0*s**2 + 4915200.0*s + 49152000.0)/(s*(128000.0*s**2 + 512000.0*s + 13312000.0))")).s;

	//H1  и какойто сигнал
	cha.comp_signl(SIGN_V,3,5);





	print_ku(cha);


	//EL_CHAIN_L cha1("1 4 1 I 1  2 3 4 R 0.5  3 1 4 R 2  4 1 2 R s*2  5 2 3 R 1  6 3 4 R 1/(s*4)");//from MU

	//cha1.comp_par_1_iu_uns();
	//cha1.comp_h1(2);
	//cha.comp_ur_so();
	//cha.comp_h1(0);

	system("pause");
	return 0;
}


void print_ku(const EL_CHAIN &cha)
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


	f_str.resize(cha.ur_so.id_cl.size(),1);
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

	cout << "2. Анализ цепи операторным методом при действии одиночного импульса на входе." << endl;
	cout << "2.1. В соответствии с номером выполняемого варианта определить функцию передачи напряжений или токов." << endl;
	cout << "H1I(S):" << endl;
	cout << cha.h1_2.str_l << endl;
	cout << "2.2. Найти нули и полюсы функции передачи и показать их расположение на плоскости комплексной частоты. По значениям полюсов функции передачи дать заключение о характере и практической длительности переходного процесса. " << endl;

	cout << "POL:";
	for (int i = 0; i < cha.h1_2.v_polus.size(); i++)
		cout << comp_to_s(cha.h1_2.v_polus[i])<<"  ";
	cout << endl;

	cout << "ZER:";
	for (int i = 0; i < cha.h1_2.v_zero.size(); i++)
		cout << comp_to_s(cha.h1_2.v_zero[i]);
	cout << endl;

	cout << "2.3. Определить переходную h1(t) характеристику цепи, сравнить с найденной в п. 1.2 задания. Проверить h1(0) и h1(inf) по аналитическому выражению h1(t) и непосредственно по схеме цепи." << endl;
	cout << cha.h1_2.sym_str << endl;
	cout << "2.4. Определить изображение по Лапласу входного одиночного импульса." << endl;
	//
	cout << "2.5. Определить изображение выходного сигнала и далее найти реакцию I2(t) или U2(t) во временной области. Построить графики входного и выходного сигналов на одном рисунке." << endl;
	cout << "F2(S):" << endl;
	cout << cha.f2_s << endl;
	cout << "f2(t):" << endl;
	cout << cha.f2_t << endl;





}



