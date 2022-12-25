#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Matrix algebra.cpp"
#include "Task I.cpp"
#include "Task II.cpp"
#include "Task III.cpp"
#include "Task IV.cpp"
using namespace std;

void task1() {
	vector<double> Grid = { -2, -1.2, -0.4, 0.4, 1.2, 2 };
	vector<double> GridValues = { 0.30564, 0.79377, 1.5901, -0.50436, -1.2492, -1.3828 };
	RatInterpolOpt(Grid, GridValues, 0);
	RatInterpolOpt(Grid, GridValues, 1);
	RatInterpolOpt(Grid, GridValues, 2);
}

void task2() {
	vector<vector<double>> Av = { {5, 8, 1, -1}, {8, -2, 6, 3}, {1, 6, 1, 2}, {-1, 3, 2, 0 } };//symmetric matrix
	Matrix A(Av);
	//Matrix ASquare = A * A;
	//ASquare.print();
	//cout << "\n";
	//(A.T()).print();
	vector<double> Yv = { 2, -7, -5, 2 };
	Vector Y(Yv);
	Vector X = LU(A, Y);
	//Vector X = SRM(A, Y);
	X.print();

	vector<double> Grid = { -2, -1.2, -0.4, 0.4, 1.2, 2 };
	vector<double> GridValues = { 0.30564, 0.79377, 1.5901, -0.50436, -1.2492, -1.3828 };
	RatInterpolOptLU(Grid, GridValues, 2);
}

void task3() {
	vector<vector<double>> Av = { {5, 8, 1, -1}, {8, -2, 6, 3}, {1, 6, 1, 2}, {-1, 3, 2, 0 } };//symmetric matrix
	Matrix A(Av);
	vector<double> Yv = { 2, -7, -5, 2 };
	Vector Y(Yv);
	Vector X = SRM(A, Y);
	X.print();
}

void task4() {
	vector<vector<double>> Diag5 = { {1, -1, 2}, {2, 2, -2, 3}, {-1, 1, 1, -1, 4}, {2, -2, 2, 2, -5}, {1, -1, 1, 6}, {2, 2, -7} };
	vector<double> Y5 = { 2, 1, 2, 1, 2, 1 };
	vector<double> X5;
	X5 = D5(Diag5, Y5);
	for (double x : X5)
	{
		cout << x << "\n";
	}
}


int main()
{
	cout.setf(ios::fixed);
	cout.precision(5);
	cout.width(10);
	//task1();
	task2();
	//task3();
	//task4();
	return 0;
}