#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "Matrix algebra.h"

using namespace std;

Vector Gauss(Matrix& A, Vector& D);

void RatInterpol(Vector& A, Vector& B, vector<double> Grid, vector<double> GridValues, int p);

vector<Vector> RatInterpolFunc(Vector& A, Vector& B, int p, double x);

int RatInterpolOpt(vector<double> Grid, vector<double> GridValues, int p);
