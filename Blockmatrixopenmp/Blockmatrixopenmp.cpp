//330074 Yunus Emre Çakır 1. öğretim
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <chrono>
#include <ctime>
#include <cmath>
#include <omp.h>
using namespace std;
void matrix_carpim_float(float** X, float** Y, float** Z, int size)
{
	int i, j, k;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				Z[i][j] += X[i][k] * Y[k][j];
			}
		}
	}
}
void matrix_carpim_double(double** X, double** Y, double** Z, int size)
{
	int i, j, k;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				Z[i][j] += X[i][k] * Y[k][j];
			}
		}
	}
}
void matrix_carpim_transpose_float(float** X, float** YT, float** Z, int size)
{
	int i, j, k;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				Z[i][j] += X[i][k] * YT[j][k];
			}
		}
	}
}
void matrix_carpim_transpose_double(double** X, double** YT, double** Z, int size)
{
	int i, j, k;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				Z[i][j] += X[i][k] * YT[j][k];
			}
		}
	}
}
void matrix_carpim_paralel_float(float** X, float** Y, float** Z, int size)
{
	int i, j, k;
#pragma omp parallel shared(X, Y, Z, size) private(i, j, k)
{
#pragma omp for schedule (static) 
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				for (k = 0; k < size; k++)
				{
					Z[i][j] += X[i][k] * Y[k][j];
				}
			}
		}
	}
}
void matrix_carpim_paralel_double(double** X, double** Y, double** Z, int size)
{
	int i, j, k;
#pragma omp parallel shared(X, Y, Z, size) private(i, j, k)
	{
#pragma omp for schedule (static) 
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				for (k = 0; k < size; k++)
				{
					Z[i][j] += X[i][k] * Y[k][j];
				}
			}
		}
	}
}
void matrix_carpim_paralel_transpose_float(float** X, float** YT, float** Z, int size)
{
	int i, j, k;
#pragma omp parallel shared(X, YT, Z, size) private(i, j, k)
	{
#pragma omp for schedule (static) 
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				for (k = 0; k < size; k++)
				{
					Z[i][j] += X[i][k] * YT[j][k];
				}
			}
		}
	}
}
void matrix_carpim_paralel_transpose_double(double** X, double** YT, double** Z, int size)
{
	int i, j, k;
#pragma omp parallel shared(X, YT, Z, size) private(i, j, k)
	{
#pragma omp for schedule (static) 
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				for (k = 0; k < size; k++)
				{
					Z[i][j] += X[i][k] * YT[j][k];
				}
			}
		}
	}
}
void block_matrix_carpim_float(float** X, float** Y, float** Z, int size, int block_size)
{
	int i, j, k;
	for (int j0 = 0; j0 < size; j0 += block_size)
	{
		for (int k0 = 0; k0< size; k0 += block_size)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = j0; j < ((j0 + block_size) > size ? size : (j0 + block_size)); j++)
				{
					for (int k = k0; k < ((k0 + block_size) > size ? size : (k0+ block_size)); k++)
					{
						Z[i][j] += X[i][k] * Y[k][j];
					}
				}
			}
		}
	}
}
void block_matrix_carpim_double(double** X, double** Y, double** Z, int size, int block_size)
{
	int i, j, k;
	for (int j0 = 0; j0< size; j0 += block_size)
	{
		for (int k0 = 0; k0< size; k0 += block_size)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = j0; j < ((j0 + block_size) > size ? size : (j0+ block_size)); j++)
				{
					for (int k = k0; k < ((k0 + block_size) > size ? size : (k0 + block_size)); k++)
					{
						Z[i][j] += X[i][k] * Y[k][j];
					}
				}
			}
		}
	}
}
void block_matrix_carpim_paralel_float(float** X, float** Y, float** Z, int size, int block_size)
{
	int i = 0, j = 0, k = 0, j0 = 0, k0 = 0;
#pragma omp parallel shared(X, Y, Z, size) private(i, j, k, j0, k0)
	{	
#pragma omp for schedule (static) 
		for (j0 = 0; j0 < size; j0 += block_size)
		{
			for (k0 = 0; k0 < size; k0 += block_size)
			{
				for (i = 0; i < size; i++)
				{
					for (j = j0; j < ((j0 + block_size) > size ? size : (j0 + block_size)); j++)
					{
						for (k = k0; k < ((k0 + block_size) > size ? size : (k0 + block_size)); k++)
						{
							Z[i][j] += X[i][k] * Y[k][j];
						}
					}
				}
			}
		}
	}
}
void block_matrix_carpim_paralel_double(double** X, double** Y, double** Z, int size, int block_size)
{
	int i = 0, j = 0, k = 0, j0 = 0, k0 = 0;
#pragma omp parallel shared(X, Y, Z, size) private(i, j, k, j0, k0)
	{
#pragma omp for schedule (static) 
		for (j0 = 0; j0 < size; j0 += block_size)
		{
			for (k0 = 0; k0 < size; k0 += block_size)
			{
				for (i = 0; i < size; i++)
				{
					for (j = j0; j < ((j0 + block_size) > size ? size : (j0 + block_size)); j++)
					{
						for (k = k0; k < ((k0 + block_size) > size ? size : (k0 + block_size)); k++)
						{
							Z[i][j] += X[i][k] * Y[k][j];
						}
					}
				}
			}
		}
	}
}
void block_matrix_carpim_transpose_float(float** X, float** YT, float** Z, int size, int block_size)
{
	int i, j, k;
	for (int j0 = 0; j0 < size; j0 += block_size)
	{
		for (int k0 = 0; k0 < size; k0 += block_size)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = j0; j < ((j0 + block_size) > size ? size : (j0 + block_size)); j++)
				{
					for (int k = k0; k < ((k0 + block_size) > size ? size : (k0 + block_size)); k++)
					{
						Z[i][j] += X[i][k] * YT[j][k];
					}
				}
			}
		}
	}
}
void block_matrix_carpim_transpose_double(double** X, double** YT, double** Z, int size, int block_size)
{
	int i, j, k;
	for (int j0 = 0; j0 < size; j0 += block_size)
	{
		for (int k0 = 0; k0 < size; k0 += block_size)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = j0; j < ((j0 + block_size) > size ? size : (j0 + block_size)); j++)
				{
					for (int k = k0; k < ((k0 + block_size) > size ? size : (k0 + block_size)); k++)
					{
						Z[i][j] += X[i][k] * YT[j][k];
					}
				}
			}
		}
	}
}
void sonuc_matrix_float(float** matrix, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}
void sonuc_matrix_double(double** matrix, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}
void sure_yazdir(double sure) {
	if (sure < 60)
		cout << "Gecen sure " << sure << " saniye" << endl;
	else {
		cout << "Gecen sure "<< (int)sure/60 << " dakika " <<(int)fmod(sure,60)<<" saniye "<< endl;
	}
}
void sifirla_double(double** matrix, int size) {
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = 0.0;
		}
	}
}
void sifirla_float(float** matrix, int size) {
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = 0.0;
		}
	}
}


int main()
{
	int size, block_size;
	cout << "Matris boyutunu gir" << endl;
	cin >> size;
	cout << "Block boyutunu gir" << endl;
	cin >> block_size;
	float** XF = new float* [size];
	float** YF = new float* [size];
	float** YTF = new float* [size];
	float** ZF = new float* [size];

	double** XD = new double* [size];
	double** YD = new double* [size];
	double** YTD = new double* [size];
	double** ZD = new double* [size];
	
	cout << "Bismillahirrahmanirrahim" << endl;
	for (int i = 0; i < size; i++)
	{
		XF[i] = new float[size];
		YF[i] = new float[size];
		YTF[i] = new float[size];
		ZF[i] = new float[size];
	
		XD[i] = new double[size];
		YD[i] = new double[size];
		YTD[i] = new double[size];
		ZD[i] = new double[size];
		
		for (int j = 0; j < size; j++)
		{
			XF[i][j] = 1.0f;
			YF[i][j] = 1.0f;
			ZF[i][j] = 0.0f;
			
			XD[i][j] = 1.0;//
			YD[i][j] = 1.0;
			ZD[i][j] = 0.0;
			
		}
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			YTF[i][j] = YF[j][i];
			YTD[i][j] = YD[j][i];
		}
	}
	auto baslangic = std::chrono::system_clock::now();
	auto bitis = std::chrono::system_clock::now();
	std::chrono::duration<double> gecen_sure = bitis - baslangic;

	cout << "Seri carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_float(XF, YF, ZF, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Seri carpim basladi(double)" << endl;
	 baslangic = std::chrono::system_clock::now();
	matrix_carpim_double(XD, YD, ZD, size);
	 bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);

	cout << "Seri transpose carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_transpose_float(XF, YTF, ZF, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Seri transpose carpim basladi(double)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_transpose_double(XD, YTD, ZD, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());;
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);
	

	cout << "Paralel carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_paralel_float(XF, YF, ZF, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Paralel carpim basladi(double)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_paralel_double(XD, YD, ZD, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);
	

	cout << "Paralel transpose carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_paralel_transpose_float(XF, YTF, ZF, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Paralel transpose carpim basladi(double)" << endl;
	baslangic = std::chrono::system_clock::now();
	matrix_carpim_paralel_transpose_double(XD, YTD, ZD, size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);

	cout << "Seri blok carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	block_matrix_carpim_float(XF, YF, ZF, size,block_size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Seri blok carpim basladi(double)" << endl;
	baslangic = std::chrono::system_clock::now();
	block_matrix_carpim_double(XD, YD, ZD, size, block_size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);

	cout << "Seri blok transpose carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	block_matrix_carpim_transpose_float(XF, YTF, ZF, size, block_size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Seri blok transpose carpim basladi(double)" << endl;
	baslangic = std::chrono::system_clock::now();
	block_matrix_carpim_transpose_double(XD, YTD, ZD, size, block_size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);
	

	cout << "Paralel blok carpim basladi(float)" << endl;
	baslangic = std::chrono::system_clock::now();
	block_matrix_carpim_paralel_float(XF, YF, ZF, size, block_size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_float(ZF, size);
	sifirla_float(ZF, size);

	cout << "Paralel blok carpim basladi(double)" << endl;
	baslangic = std::chrono::system_clock::now();
	block_matrix_carpim_paralel_double(XD, YD, ZD, size, block_size);
	bitis = std::chrono::system_clock::now();
	gecen_sure = bitis - baslangic;
	sure_yazdir(gecen_sure.count());
	sonuc_matrix_double(ZD, size);
	sifirla_double(ZD, size);
	


}

