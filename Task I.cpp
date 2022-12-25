#include "Task I.h"
#pragma once

Vector Gauss(Matrix& A, Vector& Y)
{
	if (A.getcolumns() != A.getrows()) //работаем только с квадратными матрицами
	{
		throw "\n Your matrix is not quadratic. \n";
	}
	if (Y.getrows() != A.getcolumns()) //проверяем согласованность матрицы и вектора
	{
		throw "\n Wrong values vector. \n";
	}
	Matrix B = A; //не трогаем матрицу A
	Vector D = Y; //не трогаем вектор D
	int N = B.getrows();
	for (int j = 1; j < N; j++)
	{
		double dux = abs(B.getelement(j, j)); //поиск ведущего элемента
		int dux_row = j;
		for (int i = j + 1; i <= N; i++)
		{
			if (abs(B.getelement(i, j)) > dux)
			{
				dux = abs(B.getelement(i, j));
				dux_row = i;
			}
		}
		if (dux == 0) throw "\n Your matrix is degenerate. \n"; //если ведущий элемент равен 0, то столбец нулевой => матрица вырождена
		if (dux_row > j) //если ведущий элемент стоит не на главной диагонали, то переставляем текущую строку со строкой ведущего элемента
		{
			B.permutation(j, dux_row);
			D.permutation(j, dux_row);
		}
		B.scalar(j, 1 / B.getelement(j, j)); //нормируем текущую строку
		D.scalar(j, 1 / B.getelement(j, j)); //нормируем вектор решений
		for (int i = j + 1; i <= N; i++) //вычитаем текущую строку из всех последующих с нормировкой
		{
			double b = B.getelement(i, j);
			B.substract(i, j, b);
			D.substract(i, j, b);
		}
	}
	D.setelement(N, D.getelement(N) / B.getelement(N, N)); //на последней итерации нормируем последнюю строку -- таковая состоит из единственного ненулевого элемента B(n;n)
	B.setelement(N, N, 1.0);

	B.print();
	cout << "\n";
	D.print();
	cout << "\n";

	Vector X(N); //строим вектор решений
	for (int i = N; i >= 1; i--)
	{
		double s = 0;
		for (int k = N; k > i; k--)
		{
			s += B.getelement(i, k) * X.getelement(k);
		}
		double x = D.getelement(i) - s;
		X.setelement(i, x);
	}
	X.print();
	cout << "\n";

	return X;
}

vector<Vector> RatInterpol(vector<double> Grid, vector<double> GridValues, int p)
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

	Vector X = Gauss(SLAU, Y);
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

double RatInterpolFunc(vector<Vector> Coeffs, int p, double x)
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

int RatInterpolOpt(vector<double> Grid, vector<double> GridValues, int p)
{
	cout << "Start RIO\n";
	int G = 1000;
	string FileName = "task1/Rational_interpolation_with_p_";
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