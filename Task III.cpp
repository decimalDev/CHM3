#include "Task III.h"
#include<functional>

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

double expression(double x) { //функция из варианта
	return exp(x) / (1 + x * x);
}

double expression_with_x(double x,int k) { //функция из варианта
	return exp(x)*std::pow(x,k) / (1 + x * x);
}

double Integral_x(int m,int n, double a,double b) {
	return (std::pow(b, m + n + 1) - std::pow(a, m + n + 1)) / (m + n + 1);
}
//template<typename Func>
//double simpson(Func expr, double a, double b, int n) {
double simpson(std::tr1::function<double(double)> expr, double a, double b, int n) {

	if (n % 2 == 1) {
		//std::cout << "error" << std::endl;
		//return 0;
		n--;
	}
	double h = (b - a) / n;
	std::vector<double> func;
	double S = expr(a);
	for (int i = 1; i < n; i += 2) {
		S += 4 * (expr(a + h * i));
	}

	for (int i = 2; i < n; i += 2) {
		S += 2 * (expr(a + h * i));
	}
	S += expr(b);
	return (S * h / 3);
}

int _pow = 0;
double expression_with_xp(double x) { //функция из варианта
	return exp(x) * std::pow(x, _pow) / (1 + x * x);
}

Vector Least_Squares_Method(std::vector<double> nodes,double a,double b,int n) {
	std::vector<std::vector<double>> A;
	std::vector<double> temp;
	for (int i = 0; i <= n; i++) {
		temp.clear();
		for (int j = 0; j <= n; j++)
			temp.push_back(Integral_x(i,j,a,b));
		A.push_back(temp);
	}
	std::vector<double> B;
	for (int i = 0; i <= n; i++) {
		_pow = i;
		B.push_back(simpson(expression_with_xp,0, 2, n));
	}

	Matrix left_side_of_equations(A);
	Vector right_side_of_equations(B);
	return SRM(left_side_of_equations, right_side_of_equations);

}
double Aproximation(double x,Vector coef) {
	double sum = 0;
	for (int i = 0; i < coef.getrows(); i++)
		sum += coef.getelement(i + 1) * std::pow(x, i);
	return sum;
}