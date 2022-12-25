//#include "Task I.cpp"
#include "Task II.h"

Vector LU(Matrix& A, Vector& Y)
{
	if (A.getcolumns() != A.getrows()) //работаем только с квадратными матрицами
	{
		throw "\n Your matrix is not quadratic. \n";
	}
	if (Y.getrows() != A.getcolumns()) //проверяем согласованность матрицы и вектора
	{
		throw "\n Wrong values vector. \n";
	}

	int N = A.getcolumns();
	Matrix B = A;
	Vector D = Y;
	cout << "Erfolg!\n";

	Matrix L(N, N, 1);
	Matrix U(N, N);

	for (int i = 1; i <= N; i++) //заполняем матрицы L и U
	{
		for (int j = 1; j < i; j++)
		{
			double s = 0;
			for (int k = 1; k <= j; k++)
			{
				s += L.getelement(i, k) * U.getelement(k, j);
			}
			L.setelement(i, j, 1 / U.getelement(j, j) * (B.getelement(i, j) - s));
		}
		for (int j = i; j <= N; j++)
		{
			double s = 0;
			for (int k = 1; k <= i; k++)
			{
				s += L.getelement(i, k) * U.getelement(k, j);
			}
			U.setelement(i, j, B.getelement(i, j) - s);
		}
	}
	L.print();
	cout << "\n";
	U.print();
	cout << "\n";
	
	Vector Lsol(N);
	Lsol.setelement(1, D.getelement(1) / L.getelement(1, 1));
	for (int j = 2; j <= N; j++)
	{
		double s = 0;
		for (int k = 1; k <= j - 1; k++)
		{
			s += L.getelement(j, k) * Lsol.getelement(k);
		}
		Lsol.setelement(j, (1 / L.getelement(j, j) * (D.getelement(j) - s)));
	}

	Vector Usol(N);
	Usol.setelement(N, Lsol.getelement(N));
	for (int k = N - 1; k >= 1; k--)
	{
		double s = 0;
		for (int j = k + 1; j <= N; j++)
		{
			s += U.getelement(k, j) * Usol.getelement(j);
		}
		Usol.setelement(k, (1/(U.getelement(k, k)) * (Lsol.getelement(k) - s)));
	}
	Usol.print();

	return Usol;
}



vector<Vector> RatInterpolLU(vector<double> Grid, vector<double> GridValues, int p)
{
	int N = Grid.size();
	Matrix SLAU(N, N);
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= p; j++)
		{
			SLAU.setelement(i, j, pow(Grid[i - 1], j - 1));
		}
		for (int j = p + 1; j <= N; j++)
		{
			SLAU.setelement(i, j, -GridValues[i - 1] * pow(Grid[i - 1], j - p - 1));
		}
	}
	cout << "SLAU created!" << endl;
	SLAU.print();
	cout << endl;
	Vector Y(N);
	for (int i = 1; i <= N; i++)
	{
		Y.setelement(i, -pow(Grid[i - 1], p));
	}
	cout << "Y created!\n";

	Vector X = LU(SLAU, Y);
	Vector Anew(p);
	for (int i = 1; i <= p; i++)
	{
		Anew.setelement(i, X.getelement(i));
	}
	cout << "Anew created!\n";
	Vector Bnew(N - p);
	for (int i = 1; i <= N - p; i++)
	{
		Bnew.setelement(i, X.getelement(p + i));
	}
	cout << "Bnew created!\n";

	vector<Vector> Res = { Anew, Bnew };
	return(Res);
}

double RatInterpolFuncLU(vector<Vector> Coeffs, int p, double x)
{
	Vector A = Coeffs[0];
	Vector B = Coeffs[1];
	double P, Q;
	P = pow(x, p);
	Q = 0;
	for (int i = 1; i <= A.getrows(); i++)
	{
		P += A.getelement(i) * pow(x, i - 1);
	}
	for (int i = 1; i <= B.getrows(); i++)
	{
		Q += B.getelement(i) * pow(x, i - 1);
	}
	return (P / Q);
}

int RatInterpolOptLU(vector<double> Grid, vector<double> GridValues, int p)
{
	cout << "Start RIO\n";
	int G = 1000;
	string FileName = "task2/Rational_interpolation_with_p_";
	FileName += to_string(p);
	FileName += ".txt";
	cout << "FileName created!\n";
	ofstream fout;
	fout.open(FileName, ios::out);
	double h = 4.0 / G;
	Vector A(p);
	cout << "A created as nulls!\n";
	Vector B(Grid.size() - p);
	cout << "B created as nulls!\n";
	vector<Vector> Coeffs = RatInterpol(Grid, GridValues, p);
	cout << "RatInterpol!\n";
	for (int i = 0; i <= G; i++)
	{
		double x = -2.0 + i * h;
		fout << x << "\t" << RatInterpolFunc(Coeffs, p, x) << endl;
	}
	fout.close();
	return 0;
}