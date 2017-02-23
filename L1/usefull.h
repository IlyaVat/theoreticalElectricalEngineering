#pragma once
#undef _DEBUG

#include <stdlib.h>
#include <algorithm>
#include <string>
#include <complex>
#include <iostream>
#include <vector>
//#include <boost/numeric/mtl/mtl.hpp>
//#include <boost/numeric/itl/itl.hpp>

#include <cmath>
#include "el_cha.h"

//using namespace mtl;
//using namespace itl;
using namespace std;


#define _DEBUG


template <typename T>
class MATRIX
{
public:
	vector <T> el;
	int x, y;
	MATRIX();
	MATRIX(int _x, int _y);
	T* operator [](int _x);
	friend ostream& operator<<(ostream& os, const MATRIX& dt)
	{
		stringstream ss;
		vector<string> mm;
		string res;

		int maxl = 0;
		mm.resize(dt.el.size());
		for (int i = 0; i < dt.el.size(); i++)
		{
			ss.str("");
			ss.clear();
			//cout << dt.el[i];
			ss << dt.el[i];
			ss >> mm[i];
			maxl = max(maxl, (int)mm[i].length());
		}

		for (int i = 0; i < dt.x; i++)
		{
			res += "[";
			for (int r = 0; r < dt.y; r++)
			{
				res += mm[r + i*dt.y];
				for (int o = 0; o < maxl + 1 - mm[r + i*dt.y].length(); o++)
					res += " ";
			}
			res.resize(res.size() - 1);
			res += "]\n";
		}

		os << res;
		return os;

	}
	int size_x();
	int size_y();
	void resize(int x, int y);
};

template <typename T>
int MATRIX<T>::size_x()
{
	return x;
}

template <typename T>
void MATRIX<T>::resize(int _x, int _y)
{
	if (_x < 1) _x = 1;
	if (_y < 1) _y = 1;

	x = _x;
	y = _y;
	el.resize(x*y);

}

template <typename T>
int MATRIX<T>::size_y()
{
	return y;
}

template <typename T>
MATRIX<T>::MATRIX()
{
	x = 0;
	y = 0;
}
template <typename T>
MATRIX<T>::MATRIX(int _x, int _y)
{
	if (_x < 1) _x = 1;
	if (_y < 1) _y = 1;
	x = _x;
	y = _y;
	el.resize(x*y);
}
template <typename T>
T* MATRIX<T>::operator[](int _x)
{
	return &el[_x*y];
}


typedef MATRIX<double> Matrix;
typedef MATRIX<complex<double>> MatrixC;

typedef MATRIX<string> M_S;

typedef vector<double> Vec;
typedef vector<complex<double>> VecC;




string ftos(float number);
double atof1(string s);
string comp_to_s(complex<double> v);


template <class T>
void putvec(vector<T> V)
{
	for (int i = 0; i < V.size();i++)
	cout << " , "<<V[i];
}

template <class T>
vector<T> solveAXB(MATRIX<T> A, vector<T> B)
{
	vector<T> res(A.size_x(), 0.0);


	if (A.size_x() == 0)
		return res;

	cout << "\nA::\n" << A << endl;
	//cout << "\nB::\n" << B << endl;
	//cout << "\nx::\n" << res << endl;
	cout << "\nB::\n";
	putvec(B);
	cout << endl;

	//1) треугольный вид
	//надо пройтись по всем клеткам матрицы под главной диагональю
	//гаусс вобщето предпологает что возможно бесконечное количество решений
	//но думаю и одного хватит
	//(A[y][x])

	int n = A.size_x();

	for (int i = 0; i < n; i++)
	{
		if (A[i][i]==0)
		{
			int j = n-1;
			for (int r = i + 1; r < n; r++)
			{
				if (!(A[r][i] == 0))
				{
					j = r;
				}
			}

			for (int r = 0; r < n; r++)
			{
				swap(A[j][r], A[i][r]);	
			}
			swap(B[j], B[i]);

			//cout << "\nswap "<<j<<"  "<<i<<"\n";
			//cout << A;
			//cout << "\n========================\n";
		}

		for (int r = i + 1; r < n; r++)
		{//cnk=ank/akk
			T c = A[r][i] / A[i][i];
			for (int o = 0; o < n; o++)
			{
				A[r][o] = A[r][o] - A[i][o] * c;
			}
			B[r] = B[r] - B[i] * c;
		}
		//cout << A;
		//cout << "\n========================\n";
	}

	for (int i = 0; i < n; i++)
		for (int r = 0; r < i; r++)
		{//cnk=ank/akk
			T c = A[r][i] / A[i][i];
			for (int o = 0; o < n; o++)
			{
				A[r][o] = A[r][o] - A[i][o] * c;
			}
			B[r] = B[r] - B[i] * c;
		}

	cout << A;
	putvec(B);
	cout<< endl;

	for (int i = 0; i < n; i++)
		res[i] = B[i] / A[i][i];

	putvec(res);
	cout << endl;
	//system("pause");
	//bicgstab(A, res, B, P, iter);


	return res;
}

