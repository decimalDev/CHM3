#include "Matrix algebra.h"

// класс ћј“–»÷ј характеризуетс€ трем€ пол€ми: числом строк (rows), числом столбцов (columns) и двумерным массивом элементов (content)

int Matrix::getrows()
{
	return rows;
}

int Matrix::getcolumns()
{
	return columns;
}

void Matrix::setcontent(int m, int n) // по дефолту content -- нулева€ матрица соответствующего размера
{
	rows = m;
	columns = n;
	vector<vector<double>> Theta;
	Theta.resize(m);
	for (int i = 0; i < m; i++)
	{
		Theta[i].resize(n);
		fill(Theta[i].begin(), Theta[i].end(), 0.0);
	}
	content = Theta;
}

//метод, позвол€ющий заполнить выбранную €чейку нужным числом
void Matrix::setelement(int i, int j, double a)
{
	if (i < 1 || j < 1 || i > rows || j > columns)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		content[i-1][j-1] = a;
	}
}

//метод, извлекающий значение из нужной €чейки
double Matrix::getelement(int i, int j)
{
	if (i < 1 || j < 1 || i > rows || j > columns)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		return content[i-1][j-1];
	}
}

//дефолтный конструктор строит матрицу 3x3
Matrix::Matrix()
{
	rows = 3;
	columns = 3;
	setcontent(3,3);
}
	
//конструктор, стро€щий нулевую матрицу m x n
Matrix::Matrix(int m, int n)
{
	rows = m;
	columns = n;
	setcontent(m, n);
}

Matrix::Matrix(int m, int n, double k)
{
	rows = m;
	columns = n;
	vector<vector<double>> E;
	E.resize(m);
	for (int i = 0; i < m; i++)
	{
		E[i].resize(n);
		fill(E[i].begin(), E[i].end(), 0.0);
		if (i < n)
		{
			E[i][i] = k;
		}
	}
	content = E;
}

//конструктор, стро€щий матрицу m x n и заполн€ющий еЄ элементами данного массива C
Matrix::Matrix(vector<vector<double>>& C)
{
	int m = C.size();
	int n = C[0].size();
	setcontent(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			content[i][j] = C[i][j];
		}
	}
	cout << "Matrix created!" << endl;
}

//ѕоложить элемент в вектор
void Vector::setelement(int i, double a)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		content[i-1][0] = a;
	}
}

//получить элемент вектора
double Vector::getelement(int i)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		return content[i-1][0];
	}
}

//дефолтный конструктор вектора -- нулева€ матрица 3 x 1
Vector::Vector()
{
	rows = 3;
	columns = 1;
	setcontent(3, 1);
}

//нулевой вектор на m строках
Vector::Vector(int m)
{
	rows = m;
	columns = 1;
	setcontent(m, 1);
}

//вектор размера m, заполненный элементами массива C
Vector::Vector(vector<double>& C)
{
	int m = C.size();
	columns = 1;
	setcontent(m,1);
	for (int i = 0; i < m; i++)
	{
		content.at(i).at(0) = C.at(i);
	}
}

//оператор присваивани€ дл€ векторов
Vector& Vector::operator= (const Vector& v)
{
	if (rows != v.rows)
	{
		throw "Error.";
	}
	for (int i = 0; i < v.rows; i++)
	{
		content[i-1].at(0) = v.content[i-1].at(0);
	}
	return *this;
}

void Matrix::permutation(int i, int j)
{
	if (i < 1 || i > rows || j  < 1 || j > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		vector <double> buf = content[i - 1];
		content[i - 1] = content[j - 1];
		content[j - 1] = buf;
	}
}

void Matrix::scalar(int i, double a)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		for (int j = 1; j <= columns; j++) content[i - 1][j - 1] = content[i - 1][j - 1] * a;
	}
}

void Matrix::substract(int i, int j, double a)
{
	if (i < 1 || i > rows || j < 1 || j > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		for (int k = 1; k <= columns; k++)
			content[i - 1][k - 1] = content[i - 1][k - 1] - a * content[j - 1][k - 1];
	}
}

void Matrix::print()
{
	for (int i = 1; i <= rows; i++)
	{
		for (int j = 1; j < columns; j++) cout << this->getelement(i, j) << " ";
		cout << this->getelement(i, columns) << endl;
	}
}

Matrix operator* (const Matrix& A, const Matrix& B)
{
	if (A.columns != B.rows)
	{
		throw "\n GIVEN MATRICES ARE NOT MULTIPLIABLE.\n";
	}
	int m = A.rows;
	int n = A.columns;
	int p = B.columns;
	Matrix C(m, p);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			float z = 0;
			for (int k = 0; k < n; k++)
			{
				float x = A.content[i][k];
				float y = B.content[k][j];
				z += (x * y);
			}
			C.content[i][j] = z;
		}
	}
	return C;
}

Vector operator* (const Matrix& A, const Vector& X)
{
	if (A.columns != X.rows)
	{
		throw "\n GIVEN MATRICES ARE NOT MULTIPLIABLE.\n";
	}
	int m = A.rows;
	int n = A.columns;
	Vector Y(m);
	for (int i = 0; i < m; i++)
	{
		float b = 0;
		for (int k = 0; k < n; k++)
		{
			float a = A.content[i][k];
			float x = X.content[0][k];
			b += (a * x);
		}
		Y.content[0][i] = b;
	}
	return Y;
}

Matrix Matrix::T()
{
	Matrix T(columns, rows);
	for (int i = 1; i <= rows; i++)
	{
		for (int j = 1; j <= columns; j++)
		{
			T.setelement(j, i, this->content[i - 1][j - 1]);
		}
	}
	return T;
}