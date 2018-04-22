#include <iostream>
#include <fstream>
#include <cmath>
#include <Windows.h>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;

#define MAX_ITER 1000
double const EPS = 1e-15;
const double D = 1.0;
const double LAMBDA = 1.0;
const double PI = 3.141592653589793238462643;

double norma(double *p, int size)
{
	double max = fabs(p[0]);

	for (int i = 0; i < size; i++)
	{
		if (fabs(p[i]) > max)
			max = abs(p[i]);
	}

	return max;
}

/*double norma(double *p)                         //norma max z wektora
{
	if (fabs(p[0])>fabs(p[1]) && fabs(p[0])>fabs(p[2]) && fabs(p[0])>fabs(p[3])) return fabs(p[0]);
	else if (fabs(p[1])>fabs(p[2]) && fabs(p[1])>fabs(p[2])) return fabs(p[1]);
	else if (fabs(p[2])>fabs(p[3])) return fabs(p[2]);
	else return fabs(p[3]);
}*/

double* res(double **a, double *x, double *b, int n) //obliczanie wartosci wektora residualnego
{
	double suma = 0.0;
	double *c = new double[n];
	int i, j;

	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
			suma = suma + a[i][j] * x[j];

		c[i] = suma - b[i];
		suma = 0.0;
	}

	return c;
}

double *r(double *a, double *b, int n) {

	double *c = new double[n];
	for (int i = 0; i<n; i++)
		c[i] = a[i] - b[i];

	return c;
}

void SOR(double **A, double *b, double *wyn, double omega, int n, int m) //funkcja stosujaca metode SOR
{
	int i, j, l = 1;
	double *c = new double[m];                                //wektor do przechowywania poprzednich wartosci pierwiastka
	double t = 0.0;                                             //zmienna pomocnicza

	while (1)
	{
		for (i = 0; i<m; i++)
		{
			c[i] = wyn[i];                                     //przechowanie poprzedniej wartosci pierwiastka

			for (j = 0; j <= i - 1; j++)
				t = t + A[i][j] * wyn[j];                   //suma elementow dolnej macierzy pomnozona przez wektor

			for (j = i; j<m; j++)
			{
				if (j == i) t = t + (1 - (1.0 / omega))*A[i][i] * wyn[j];    //elementy przekatnej pomnozone przez parametr
				else t = t + A[i][j] * wyn[j];                       //suma elementow gornej macierzy pomnoz przez wektor
			}

			wyn[i] = (1.0 / A[i][i])*(b[i] - t)*omega;                  //rozwiazanie ukladu z dolna macierza
			t = 0.0;
		}

		if (fabs((norma(res(A, wyn, b, n), m)))<EPS && fabs((norma(r(wyn, c, n), m)))<EPS) break;  //warunek zakonczenia iteracji wartosci wektora residualnego do zera
		if (l>MAX_ITER) break;                                     //warunek zakonczenia iteracji skonczona ilosc iteracji

		l++;
	}
}

void Thomas(double *l, double *d, double *u, double *b, double *x, int m)
{
	for (int i = 2;i<m;i++)
	{
		d[i] = d[i] - (l[i - 1] / d[i - 1])*u[i - 1];
		b[i] = b[i] - (l[i - 1] / d[i - 1])*b[i - 1];
	}
	x[m - 1] = b[m - 1] / d[m - 1];
	for (int i = m - 2;i >= 0;i--)
		x[i] = (b[i] - u[i] * x[i + 1]) / d[i];

}

double *newVector(int n) {
	double *x = new double[n];
	return x;
}

void crank_nicolson_Thomas(double **U, int N, int M, double h, double dt, fstream &file2)
{
	double *a = newVector(M);
	double *b = newVector(M);
	double *c = newVector(M);
	double *d = newVector(M);
	double *u = newVector(M);
	double xh = 0;

	// Brzegowe
	for (int i = 0; i < N; i++)
		U[i][0] = 0.0;
	for (int i = 0; i < N; i++)
		U[i][M-1] = 0.0;

	// Poczatkowe
	for (int i = 0; i < M-1; i++)
	{
		U[0][i] = sin(PI*xh);
		xh = xh + h;
	}

	for (int k = 1; k < N; k++) {
		a[0] = 0.;
		b[0] = U[k - 1][0];
		c[0] = 0.;
		d[0] = 1.;

		for (int i = 1; i < M - 1; i++) {
			a[i] = 1.0 / 2.0;
			b[i] = -(1.0 / 2.0*U[k - 1][i - 1] + (1.0 - 1.0)*U[k - 1][i] + (1.0 / 2.0)*U[k - 1][i + 1]);
			c[i] = 1.0 / 2.0;
			d[i] = -(1 + 1.0);
		}

		a[M - 1] = 0.0;
		b[M - 1] = U[k - 1][M - 1];
		c[M - 2] = 1.0 / 2.0;
		d[M - 1] = 1.0;

		Thomas(a, d, c, b, u, M);
		for (int i = 1; i < M - 1; i++)
			U[k][i] = u[i];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			file2 << U[i][j] << ";";
		}
		file2 << "\n";
	}
}

void crank_nicolson_SOR(double **U, int N, int M, double h, double dt, fstream &file3)
{
	double *a = newVector(M);
	double *b = newVector(M);
	double *c = newVector(M);
	double *d = newVector(M);
	double *u = newVector(M);
	double xh = 0;

	double **C;
	C = new double*[N];
	for (int i = 0;i<N;i++)
		C[i] = new double[M];

	// Brzegowe
	for (int i = 0; i < N; i++)
		U[i][0] = 0.0;
	for (int i = 0; i < N; i++)
		U[i][M - 1] = 0.0;

	// Poczatkowe
	for (int i = 0; i < M - 1; i++)
	{
		U[0][i] = sin(PI*xh);
		xh = xh + h;
	}

	for (int k = 1; k < N; k++) {
		a[0] = 0.;
		b[0] = U[k - 1][0];
		c[0] = 0.;
		d[0] = 1.;

		for (int i = 1; i < M - 1; i++) {
			a[i] = 1.0 / 2.0;
			b[i] = -(1.0 / 2.0*U[k - 1][i - 1] + (1.0 - 1.0)*U[k - 1][i] + (1.0 / 2.0)*U[k - 1][i + 1]);
			c[i] = 1.0 / 2.0;
			d[i] = -(1 + 1.0);
		}

		a[M - 1] = 0.0;
		b[M - 1] = U[k - 1][M - 1];
		c[M - 2] = 1.0 / 2.0;
		d[M - 1] = 1.0;

		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < M; j++)
			{
				C[i][j] = 0.0;
			}
		}
		for (int i = 0; i < M - 1; i++)
		{
			C[i][i] = d[i];
			C[i + 1][i] = a[i];
			C[i][i + 1] = c[i];
		}
		C[M - 1][M - 1] = d[M - 1];

		C[1][0] = a[0];

		SOR(C, b, u, 0.5, M, M);
		for (int i = 1; i < M - 1; i++)
			U[k][i] = u[i];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			file3 << U[i][j] << ";";
		}
		file3 << "\n";
	}
}

void zmienne(double **Zmienne, int n, int m, double h, double dt, fstream &file6)
{
	double x = h;
	double t = 0.0;
	for (int i = 1; i < n - 1; i++)
	{

		Zmienne[i + 1][0] = x;
		file6 << Zmienne[i + 1][0] << ";";
		Zmienne[i + 1][1] = log10(x);
		file6 << Zmienne[i + 1][1] << ";";
		x = h + x;

		Zmienne[i + 1][2] = t;
		file6 << Zmienne[i + 1][2] << ";";
		t = t + dt;
		file6 << "\n";
	}

}

void analitycznie(double **Analit, int n, int m, double h, double dt, fstream &file)
{
	double t = 0.0;
	double x = 0.0;

	for (int i = 0;i<n;i++)
	{
		for (int j = 0;j <= m;j++)
		{
			Analit[i][j] = exp(-PI*PI*D*t)*sin(PI*x);
			x = x + h;
		}
		x = 0.0;
		t = t + dt;
	}

	//ZAPISYWANIE DO PLIKU
	for (int i = 0;i<n;i++)
	{
		for (int j = 0;j <= m;j++)
		{
			file << Analit[i][j] << ";";
		}
		file << "\n";
	}
}

// LICZYMY BLEDY DLA KAZDEGO t z krokiem dt
void blad(double **Analit, double **A, double **Bledy, int n, int m, fstream &file4)
{
	for (int i = 0;i<n;i++)
	{
		for (int j = 0;j<m;j++)
		{
			Bledy[i][j] = log(fabs(A[i][j] - Analit[i][j]));
			//Bledy[i][j] = fabs(A[i][j] - Analit[i][j]);
			file4 << Bledy[i][j] << ";";
		}
		file4 << "\n";
	}
}

//LICZYMY MAXYMALNE BLEDY W WIERSZACH, Do WYKRESU  kroku dt od kroku dx
void maxBlad(double **blad, int N, int M, fstream &file5)
{
	// DLA Bledow_T
	double *maxblad;
	maxblad = new double[N];
	for (int k = 0;k<N;k++)
		maxblad[k] = 1;

	double max;
	for (int i = 0; i < N; i++)
	{
		max = blad[i][0];
		for (int j = 1; j < M; j++)
		{
			if (max < blad[i][j])
				max = blad[i][j];
		}
		maxblad[i] = max;
	}
	for (int i = 0; i < N; i++)
		file5 << maxblad[i] << endl;
}

//MAIN
int main()
{
	fstream file, file2, file3, file4, file5, file6, file7;
	file.open("Analit.csv", ios::out);
	file2.open("CN_T.csv", ios::out);
	file3.open("CN_SOR.csv", ios::out);
	file4.open("Bledy.csv", ios::out);
	file5.open("Bledy_Max.csv", ios::out);
	file6.open("Zmienne.csv", ios::out);
	file7.open("Max_Blad_(1).csv", ios::out);

	double lambda = 1;
	double xmin = 0.0, xmax = 1.0, tmin = 0.0;
	double tmax = 0.5;

	double h = 0.025;
	cout << "h = " << log10(h) << endl;
	double dt = lambda*h*h / D;
	//double dt = 0.001;
	int n = ((tmax - tmin) / dt) + 1;
	int m = ((xmax - xmin) / h);

	cout << "n=" << n << " m=" << m << endl;
	cout << endl << "LICZY!!" << endl;

	double **Analit;
	double **CN_T;
	double **CN_SOR;
	double **Bledy;
	double **Bledy_Max;
	double **Zmienne;

	// SZYKOWANIE PAMIÊCI
	Analit = new double*[n + 1];
	for (int i = 0; i<n + 1; i++)
		Analit[i] = new double[m + 1];

	CN_T = new double*[n+1];
	for (int i = 0;i<n+1;i++)
		CN_T[i] = new double[m+1];

	CN_SOR = new double*[n + 1];
	for (int i = 0;i<n + 1;i++)
		CN_SOR[i] = new double[m + 1];

	Bledy = new double*[n + 1];
	for (int i = 0;i<n + 1;i++)
		Bledy[i] = new double[m + 1];

	Bledy_Max = new double*[n + 1];
	for (int i = 0;i<n + 1;i++)
		Bledy_Max[i] = new double[m + 1];

	Zmienne = new double*[n + 1];
	for (int i = 0; i<n + 1; i++)
		Zmienne[i] = new double[m + 1];

	double **kroczek;
	kroczek = new double*[n + 1];
	for (int i = 0; i<n + 1; i++)
		kroczek[i] = new double[m + 1];

	// WYWOLYWANIE FUNKCJI
	analitycznie(Analit, n, m, h, dt, file);
	//crank_nicolson_Thomas(CN_T, n, m+1, h, dt, file2);
	crank_nicolson_SOR(CN_SOR, n, m+1, h, dt, file3);
	//blad(Analit, CN_T, Bledy, n, m, file4);
	blad(Analit, CN_SOR, Bledy, n, m, file4);
	maxBlad(Bledy, n, m, file5);
	zmienne(Zmienne, n, m, h, dt, file6);

	// Do liczenia wykresow z zadania pierwszego.
	/*for (int i = 0; i < n; i++)
	{
		Analit = new double*[n + 1];
		for (int i = 0; i<n + 1; i++)
			Analit[i] = new double[m + 1];

		CN_T = new double*[n + 1];
		for (int i = 0;i<n + 1;i++)
			CN_T[i] = new double[m + 1];


		Bledy = new double*[n + 1];
		for (int i = 0;i<n + 1;i++)
			Bledy[i] = new double[m + 1];

		analitycznie(Analit, n, m, h, dt, file);
		crank_nicolson_Thomas(CN_T, n, m, h, dt, file2);
		blad(Analit, CN_T, Bledy, n, m, file4);

		for (int j = 0; j < m; j++)
		{
			kroczek[i][j] = Bledy[n - 1][j];
		}
		cout << h << endl << dt << endl << n << endl << m << endl << endl;
		if (i == 17)
			break;
		h = h + 0.001;
		dt = lambda*h*h;
		n = ((tmax - tmin) / dt) + 1;
		m = ((xmax - xmin) / h);
	}

	maxBlad(kroczek, n, m, file7);*/


	file.close();
	file2.close();
	file3.close();
	file4.close();
	file5.close();
	file6.close();
	file7.close();
	for (int i = 0;i<n;i++)
	{
		delete[] Analit[i];
		delete[] CN_T[i];
		delete[] CN_SOR[i];
		delete[] Bledy[i];
		delete[] Bledy_Max[i];
		delete[] Zmienne[i];
		delete[] kroczek[i];
	}
	delete[] Analit;
	delete[] CN_T;
	delete[] CN_SOR;
	delete[] Bledy;
	delete[] Bledy_Max;
	delete[] Zmienne;
	delete[] kroczek;
	Analit = CN_T = CN_SOR = Bledy = Bledy_Max = Zmienne = kroczek = NULL;

	cout << endl << "GOTOWE" << endl;
	system("pause");
	return 0;
}