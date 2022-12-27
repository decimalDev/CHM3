#include "Task III.h"

double sgn(double x)
{
	if (x < 0) return -1;
	if (x == 0) return 0;
	else return 1;
}

double f_test(double x)
{
	return cbrt(x + 1) + 1;
}

double expression(double x)
{
	return(exp(x) / (1 + x * x));
}

Vector SRM(Matrix& A, Vector& Y) //дана система AX = Y
{
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

double scalar_product(std::tr1::function<double(double)> f, std::tr1::function<double(double)> g, vector<double> Grid)
{
	double s = 0;
	for (double x : Grid)
	{
		s += f(x) * g(x);
	}
	return s;
}

double scalar_product_polynomes(int n, int m, vector<double> Grid)
{
	double s = 0;
	for (double x : Grid)
	{
		s += pow(x, m) * pow(x, n);
	}
	return s;
}

double scalar_product_func_and_polynome(std::tr1::function<double(double)> f, int m, vector<double> Grid)
{
	double s = 0;
	for (double x : Grid)
	{
		s += f(x) * pow(x, m);
	}
	return s;
}

//строим вектор коэффициентов
Vector SRMApprox_coeffs(vector<double> Grid, int n, std::tr1::function<double(double)> f) //Grid -- сетка, n -- максимальная степень многочленов в базисе
{
	int n_mod = n + 1;
	int N = Grid.size(); //число точек в сетке
	Matrix A(n_mod, n_mod);
	Vector Y(n_mod);
	for (int i = 1; i <= n_mod; i++)
	{
		Y.setelement(i, scalar_product_func_and_polynome(f, i - 1, Grid));
		for (int j = 1; j <= i; j++)
		{
			double a = scalar_product_polynomes(i - 1, j - 1, Grid);
			A.setelement(i, j, a);
			A.setelement(j, i, a);
		}
	}
	Vector X = SRM(A, Y);
	return X;

}

//возвращаем значение аппроксимирующей функции в заданной точке
double SRM_Approx_val(double x, vector<double> Grid, int n, std::tr1::function<double(double)> f)
{
	int N = Grid.size();
	Vector A = SRMApprox_coeffs(Grid, n, f);
	double s = 0;
	for (int i = 1; i <= n + 1; i++)
	{
		s += A.getelement(i) * pow(x, i - 1);
	}
	return s;

}

double SRM_Approx(vector<double> Grid, int n, std::tr1::function<double(double)> f)
{
	string FileName = "C:/Users/student/Desktop/SRM_Approx/SRM_Approx_with_n_";
	FileName += to_string(n);
	FileName += ".txt";
	ofstream fout;
	fout.open(FileName, ios::out);

	double h = (Grid[Grid.size() - 1] - Grid[0]) / 1000;
	double x = 0, appr, exact, div;
	double div_max = -h;
	fout.setf(ios::fixed);
	fout.precision(10);

	for (int i = 0; i <= 1000; i++)
	{
		x += h;
		appr = SRM_Approx_val(x, Grid, n, f);
		exact = f(x);
		div = abs(appr - exact);
		if (div > div_max) div_max = div;
		fout << x << "\t" << appr << "\t" << exact << "\t" << div << endl;
	}

	div = 0;
	for (double x : Grid)
	{
		appr = SRM_Approx_val(x, Grid, n, f);
		exact = f(x);
		div += abs(appr - exact) * abs(appr - exact);
	}

	fout << "div_max: " << div_max;
	fout.close();
	return div;
}

void task_3_2()
{
	vector<double> Grid;
	double h = 0.1;
	double x = -h;
	for (int i = 0; i <= 20; i++)
	{
		x += h;
		Grid.push_back(x);
	}

	double appr_min = 1;
	int n_opt;
	double appr_cur;
	for (int n = 1; n <= 12; n++)
	{
		appr_cur = SRM_Approx(Grid, n, expression);
		if (appr_cur < appr_min) {
			n_opt = n;
			appr_min = appr_cur;
		}
	}
	cout << "optimal n equals " << n_opt;

}