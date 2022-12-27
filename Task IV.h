#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "Matrix algebra.h"

using namespace std;

vector<double> D5(vector<vector<double>>& M, vector<double>& X);

void fill_randomly(vector<double>& A);

vector<vector<double>> Generate_D5(int N, double q);

