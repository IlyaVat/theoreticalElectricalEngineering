

#include "usefull.h"

//using namespace mtl;
//using namespace itl;
using namespace std;


std::string ftos(float number)
{
	std::stringstream buff;
	
	buff << number;
	string s;
	buff >> s;
	return s;
}

double atof1(string s)
{
	myreplace(s, ".", ",");
	for (int i = 0; i < s.length(); i++)
	{
		if ((s[i]<'0' || s[i]>'9') && s[i] != ','&& s[i] != 'e'&& s[i] != 'E'&& s[i] != '-')
		{
			cout << "convert broken  " << s << " isnt number" << endl;
			system("pause");
		}

	}
	return std::atof(s.data());
}
/*
double atof1(string s)
{
	myreplace(s,".",",");
	return std::atof(s.data());
}
*/

string comp_to_s(complex<double> v)
{
	string res;

	if (v._Val[1] == 0)
	{
		return ftos(v._Val[0]);
	}
	else
	{
		if (v._Val[1] > 0)
			return ftos(v._Val[0]) + "+" + ftos(v._Val[1]) + "*I";
		else
			return ftos(v._Val[0]) + ftos(v._Val[1]) + "*I";
	}
}







bool SaveArrFile(const TCHAR* filename, const __int32* arr,
	int width, int height, int bpp)
{

	if ((bpp < 24) || (bpp > 32)) // только 24/32 бит
		return FALSE;

	DWORD p_row = (DWORD)((width * bpp + 31) & ~31) / 8uL;
	DWORD size = (DWORD)(height * p_row);

	// формируем файловый заголовок
	BITMAPFILEHEADER  hdr;
	ZeroMemory(&hdr, sizeof(BITMAPFILEHEADER));
	hdr.bfType = 0x4D42;
	hdr.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	hdr.bfSize = hdr.bfOffBits + size;

	// заголовок описателя растра
	BITMAPINFO dib;
	ZeroMemory(&dib, sizeof(BITMAPINFO));
	dib.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	dib.bmiHeader.biBitCount = bpp;
	dib.bmiHeader.biCompression = BI_RGB;
	dib.bmiHeader.biPlanes = 1u;
	dib.bmiHeader.biWidth = (long)width;
	dib.bmiHeader.biHeight = (long)-height;
	dib.bmiHeader.biSizeImage = size;
	dib.bmiHeader.biXPelsPerMeter = 11811L;
	dib.bmiHeader.biYPelsPerMeter = 11811L;
	dib.bmiHeader.biClrImportant = 0uL;
	dib.bmiHeader.biClrUsed = 0uL;

	// далее запись в файл
	HANDLE fp = CreateFile(filename, GENERIC_WRITE, FILE_SHARE_WRITE, NULL,
		CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (fp == INVALID_HANDLE_VALUE)
		return FALSE;

	// записываем заголовки...
	DWORD  dwr = 0uL;
	WriteFile(fp, (LPCVOID)&hdr, sizeof(BITMAPFILEHEADER), &dwr, NULL);
	WriteFile(fp, (LPCVOID)&dib.bmiHeader, sizeof(BITMAPINFOHEADER), &dwr, NULL);

	// запись массива пикселей
	//if (bpp == 32) // 32-бит

	vector<byte> b;
	b.resize(width*height*4);

	for (int i = 0; i < width*height;i++)
	{
		b[i * 3+0] = ((byte*)arr)[i * 3];
		b[i * 3 + 1] = ((byte*)arr)[i * 3 + 1];
		b[i * 3 + 2] = ((byte*)arr)[i * 3 + 2];
	}

		WriteFile(fp, (LPCVOID)b.data(), size, &dwr, NULL);
		/*
	else if (bpp == 24) { // 24-бит с дополнением до 32-разрядной границы

		BYTE   nil = 0u;
		int   cb = sizeof(RGBQUAD);
		int  align = ((cb - ((width*bpp + 7) / 8) % cb) % cb);

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++)
				WriteFile(fp, (LPCVOID)&arr[y*width + x], sizeof(RGBTRIPLE), &dwr, NULL);

			for (int i = 0; i < align; i++) // до границы DWORD
				WriteFile(fp, (LPCVOID)&nil, sizeof(BYTE), &dwr, NULL);
		}
	}
	*/
	FlushFileBuffers(fp);
	CloseHandle(fp);
	return TRUE;
}

