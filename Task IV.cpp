#include "Task IV.h"
#include<ctime>

//����� �� ���� ������ ���� {
//																				{c_0, d_0, e_0},
//																				{b_1, c_1, d_1, e_1},
//																				{a_2, b_2, c_2, d_2, e_2},
//																				{a_3, b_3, c_3, d_3, e_3},
//																				...
//																				{a_(n-2), b_(n-2), c_(n-2), d_(n-2), e_(n-2)},
//																				{a_(n-1), b_(n-1), c_(n-1), d_(n-1)},
//																				{a_n, b_n, c_n}
//																			}
vector<double> D5(vector<vector<double>>& M, vector<double>& X)
{

	int N = M.size() - 1;
	M[0][1] *= -1;
	M[1][0] *= -1;
	M[1][2] *= -1;
	M[N][1] *= -1;
	for (int i = 2; i <= N - 1; i++)
	{
		M[i][1] *= -1;
		M[i][3] *= -1;
	}

	Matrix ABC(N + 1, 3); //����������� ������������
	ABC.setelement(1, 1, M[0][1] / M[0][0]); // A_1 = d_0/c_0
	ABC.setelement(1, 2, M[0][2] / M[0][0]); // B_1 = e_0/c_0
	ABC.setelement(1, 3, X[0] / M[0][0]); //C_1 = f_0/c_0

	double Delta = M[1][1] - M[1][0] * ABC.getelement(1, 1); //Delta_1 = c_1 - b_1 * A_1
	cout << Delta << endl;
	ABC.setelement(2, 1, (M[1][2] - ABC.getelement(1, 2) * M[1][0]) / Delta); // A_2 = (d_1 - b_1 * B_1)/Delta_1
	ABC.setelement(2, 2, M[1][3] / Delta); // B_2 = e_1 / Delta_1;
	ABC.setelement(2, 3, (X[1] + M[1][0] * ABC.getelement(1, 3)) / Delta); //C_3 = (f_1 + b_1 * C_1) / Delta_1

	for (int i = 2; i <= N - 2; i++)
	{
		Delta = M[i][2] - M[i][0] * ABC.getelement(i - 1, 2) + ABC.getelement(i, 1) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1]); // Delta_i = c_i - a_i * B_(i-1) + A_i * (a_i * A_(i-1) - b_i)
		cout << Delta << endl;
		ABC.setelement(i + 1, 1, ((M[i][3] + ABC.getelement(i, 2) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1])) / Delta)); // A_(i+1) = (1/Delta) * (d_i + B_i * (a_i * A_(i+1) - b_i))
		ABC.setelement(i + 1, 2, (M[i][4] / Delta)); //B_(i+1) = e_i / Delta_i
		ABC.setelement(i + 1, 3, (X[i] - M[i][0] * ABC.getelement(i - 1, 3) - ABC.getelement(i, 3) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1])) / Delta); // C_(i+1) = (1/Delta) * (f_i - a_i * C_(i - 1) - C_i * (a_i * A_(i - 1) - b_i))
	}

	int i = N - 1;
	Delta = M[i][2] - M[i][0] * ABC.getelement(i - 1, 2) + ABC.getelement(i, 1) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1]); // Delta_i = c_i - a_i * B_(i-1) + A_i * (a_i * A_(i-1) - b_i)
	ABC.setelement(i + 1, 1, ((M[i][3] + ABC.getelement(i, 2) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1])) / Delta)); // A_(i+1) = (1/Delta) * (d_i + B_i * (a_i * A_(i+1) - b_i))
	ABC.setelement(i + 1, 3, (X[i] - M[i][0] * ABC.getelement(i - 1, 3) - ABC.getelement(i, 3) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1])) / Delta); // C_(i+1) = (1/Delta) * (f_i - a_i * C_(i - 1) - C_i * (a_i * A_(i - 1) - b_i))

	i = N;
	Delta = M[i][2] - M[i][0] * ABC.getelement(i - 1, 2) + ABC.getelement(i, 1) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1]); // Delta_i = c_i - a_i * B_(i-1) + A_i * (a_i * A_(i-1) - b_i)
	ABC.setelement(i + 1, 3, (X[i] - M[i][0] * ABC.getelement(i - 1, 3) - ABC.getelement(i, 3) * (M[i][0] * ABC.getelement(i - 1, 1) - M[i][1])) / Delta); // C_(i+1) = (1/Delta) * (f_i - a_i * C_(i - 1) - C_i * (a_i * A_(i - 1) - b_i))

	ABC.print();
	cout << "\n";

	vector<double> Y;
	Y.resize(N + 1);
	Y[N] = ABC.getelement(N + 1, 3); // y_N = C_(N+1)
	cout << "y_N = " << Y[N] << endl;
	Y[N - 1] = Y[N] * ABC.getelement(N, 1) + ABC.getelement(N, 3);
	cout << "y_N - 1 = " << Y[N - 1] << endl;
	for (int j = N - 2; j >= 0; j--)
	{
		Y[j] = ABC.getelement(j + 1, 1) * Y[j + 1] - ABC.getelement(j + 1, 2) * Y[j + 2] + ABC.getelement(j + 1, 3); //y_i = A_(i+1) * y_(i+1) - B_(i+1) * y_(i+2) + C_(i+1)
	}

	return Y;
}

void fill_randomly(vector<double>& A)
{
	for (int i = 0; i < A.size(); i++) A[i] = -10 + rand() % 21;
	
}

vector<vector<double>> Generate_D5(int N, double q) // N -- ����� �����, q -- ��������, � ������� �������� ������������ ������������ ������������
{
	srand(time(NULL));
	vector<vector<double>> A(N);

	A[0].resize(3);
	fill_randomly(A[0]);
	double sum = 0;
	sum += abs(A[0][1]) + abs(A[0][2]);
	A[0][0] = q * sum;
	for (double a : A[0]) cout << a << " ";
	cout << endl;


	A[1].resize(4);
	fill_randomly(A[1]);
	sum = 0;
	sum += abs(A[1][0]) + abs(A[1][2]) + abs(A[1][3]);
	A[1][1] = q * sum;
	for (double a : A[1]) cout << a << " ";
	cout << endl;

	for (int i = 2; i < N - 2; i++)
	{
		A[i].resize(5);
		fill_randomly(A[i]);
		A[i][2] = (abs(A[i][0]) + abs(A[i][1]) + (abs(A[i][3]))) * q;
		for (double a : A[i]) cout << a << " ";
		cout << endl;
	}


	A[N - 2].resize(4);
	fill_randomly(A[N - 2]);
	sum = 0;
	sum += abs(A[N - 2][0]) + abs(A[N - 2][1]) + abs(A[N - 2][3]);
	A[N - 2][2] = q * sum;
	for (double a : A[N - 2]) cout << a << " ";
	cout << endl;

	A[N-1].resize(3);
	fill_randomly(A[N - 1]);
	sum = 0;
	sum += abs(A[N - 1][0]) + abs(A[N - 1][1]);
	A[N - 1][2] = q * sum;
	for (double a : A[N - 1]) cout << a << " ";
	cout << endl;
	cout << endl;

	return A;
}

