#undef _DEBUG
#include <math.h>			// Header File For Math Library Routines
#include <stdio.h>			// Header File For Standard I/O Routines
#include <stdlib.h>			// Header File For Standard Library
#include <list>

#define _DEBUG

using namespace std;

struct pos
{
public:
	float x, y;
	
	pos()
	{
		x = 0; 
		y = 0;
	}

	pos(float _x,float _y)
	{
		x = _x;
		y = _y;
	}
	
	friend pos operator + (pos _a, pos _b)
	{
		return pos(_a.x + _b.x, _a.y + _b.y);
	}

	friend pos operator - (pos _a, pos _b)
	{
		return pos(_a.x - _b.x, _a.y - _b.y);
	}

	friend pos operator / (pos _a, int _b)
	{
		return pos(_a.x / _b, _a.y / _b);
	}

	friend pos operator * (pos _a, int _b)
	{
		return pos(_a.x * _b, _a.y * _b);
	}

	void norm()
	{
		float r = sqrt(x*x+y*y);
		if (r > 0)
		{
			x /= r;
			y /= r;
		}
		else
		{
			x = 1;
			y = 0;
		}
	}
};

struct M_ACTION
{
	int t;//1-lmd 2-lmu 3-rmd 4-rmu 5-mbd 6-mbu 7-scroll(i[0]) 8-move
	int i[2];
};

class MOUSEINFO
{
public:
	list <M_ACTION> l;
	
	int x, y, bl, br, bm, scroll;
	int xz, yz, blz, brz, bmz;

	MOUSEINFO()
	{
		x = 0; 
		y = 0;
		bl = 0;
		br = 0;
		bm = 0;
		scroll = 0;

		xz = 0;
		yz = 0;
		blz = 0;
		brz = 0;
		bmz = 0;
	}

	void addaction(int _t, int _i0,int _i1)//1-lmd 2-lmu 3-rmd 4-rmu 5-mbd 6-mbu 7-scroll(i[0]) 8-move
	{
		M_ACTION a;
		a.t = _t;
		a.i[0] = _i0;
		a.i[1] = _i1;
		l.push_back(a);
	}

	bool use_action()
	{
		if (l.empty())
		{
			return false;
		}

		M_ACTION a = *l.begin();
		
		scroll = 0;
		blz = bl;
		brz = br;
		bmz = bm;
		xz = x;
		yz = y;

		switch (a.t)
		{
		case 1:
			bl = 1;
			break;
		case 2:
			bl = 2;
			break;
		case 3:
			br = 1;
			break;
		case 4:
			br = 2;
			break;
		case 5:
			bm = 1;
			break;
		case 6:
			bm = 2;
			break;
		case 7:
			scroll = a.i[0];
			break;
		case 8:
			x = a.i[0];
			y = a.i[1];
			break;

		}

		return true;
	}
	//max_size
};


void dd(double& x, double& y, double c, double d);
double ss(double x, double y, double x1, double y1);
double ras(double x, double y);
double ras3(double x, double y, double z);


