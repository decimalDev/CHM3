#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include<functional>
#include "Matrix algebra.h"

using namespace std;

double sgn(double x);

double expression(double x);

Vector SRM(Matrix& A, Vector& Y);

double scalar_product(std::tr1::function<double(double)> f, std::tr1::function<double(double)> g, vector<double> Grid);

double scalar_product_polynomes(int n, int m, vector<double> Grid);

double scalar_product_func_and_polynome(std::tr1::function<double(double)> f, int m, vector<double> Grid);

Vector SRMApprox_coeffs(vector<double> Grid, int n, std::tr1::function<double(double)> f);

double SRM_Approx_val(double x, vector<double> Grid, int n, std::tr1::function<double(double)> f);

double SRM_Approx(vector<double> Grid, int n, std::tr1::function<double(double)> f);

void task_3_2();
