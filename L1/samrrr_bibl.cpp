#include "samrrr_bibl.h"

void dd(double& x, double& y, double c, double d)
{
	y = sin(c)*d;
	x = cos(c)*d;
}

inline double ssb(double x, double y, double x1, double y1)
{
	double u;
	u = atan((y - y1) / (x - x1));
	if (x1 - x < 0) 
	{
		u = u + 3.14159265358979323846;
	};
	return u;
}

double ss(double x, double y, double x1, double y1)
{
	double xx,yy;
	xx = abs(x1 - x);
	yy = abs(y1 - y);
	if (xx < yy)
	{
		return 3.14159265358979323846/2-ssb(y, x, y1, x1);
	}
	if (xx != 0)
	{
		return ssb(x, y, x1, y1);
	}
	return 0.45323345;
}

double ras(double x, double y)
{
	return sqrt(x*x + y*y);
}

double ras3(double x, double y, double z)
{
	return sqrt(x*x + y*y + z*z);
}
