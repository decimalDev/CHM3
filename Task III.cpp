#include "Task III.h"

double sgn(double x)
{
	if (x < 0) return -1;
	if (x == 0) return 0;
	else return 1;
}

Vector SRM(Matrix& A, Vector& Y) //дана система AX = Y
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

	Matrix S(N, N); //наша верхнетреугольная матрица S
	Matrix D(N, N); //наша диагональная матрица D

	for (int i = 1; i <= N; i++) //тут заполняем коэффициенты S и D по стандартным формулам
	{
		double s = 0; double a;
		for (int k = 1; k <= i - 1; k++)
		{
			s += D.getelement(k, k) * pow(abs(S.getelement(k, i)), 2.0);
		}
		a = A.getelement(i, i) - s;
		D.setelement(i, i, sgn(a));
		S.setelement(i, i, sqrt(abs(a)));
		double 	b = S.getelement(i, i) * D.getelement(i, i);
		for (int j = i + 1; j <= N; j++)
		{
			double s2 = 0; double a;
			for (int k = 1; k <= i - 1; k++)
			{
				s2 += S.getelement(k, i) * S.getelement(k, j) * D.getelement(k, k);
			}
			a = (A.getelement(i, j) - s2) / b;
			S.setelement(i, j, a);
		}
	}

	Vector Zeta(N); //вектора неизвестных в уравнениях (S^T * Zeta = Y), (D * Eta = Zeta) и (S * Xi = Eta); Xi -- искомый вектор 
	Vector Eta(N);
	Vector Xi(N);

	Matrix ST = S.T();
	for (int i = 1; i <= N; i++) //решаем систему (S^T | Y) как систему с нижнетреугольной матрицей
	{
		double s = 0;
		for (int j = 1; j < i; j++)
		{
			s += ST.getelement(i, j) * Zeta.getelement(j);
		}
		Zeta.setelement(i, (Y.getelement(i) - s) / ST.getelement(i, i));
		Eta.setelement(i, Zeta.getelement(i) * D.getelement(i, i)); //решаем (D | Zeta) -- поскольку D диагональная, причём любой диагональный элемент равен +1 или -1, это делается элементарно
	}

	for (int i = N; i >= 1; i--) //решаем систему (S | Eta) как систему с верхнетреугольной матрицей (стандартный обратный ход)
	{
		double s = 0;
		for (int j = N; j > i; j--)
		{
			s += Xi.getelement(j) * S.getelement(i, j);
		}
		Xi.setelement(i, (Eta.getelement(i) - s) / S.getelement(i, i));
	}

	return Xi;

}