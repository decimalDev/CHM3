#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

#pragma once

class Vector;

// класс ћј“–»÷ј характеризуетс€ трем€ пол€ми: числом строк (rows), числом столбцов (columns) и двумерным массивом элементов (content)
class Matrix
{
protected:
	int rows;
	int columns;
	vector<vector<double>> content;

public:
	int getrows();

	int getcolumns();

	void setcontent(int m, int n); // по дефолту content -- нулева€ матрица соответствующего размера

	//метод, позвол€ющий заполнить выбранную €чейку нужным числом
	void setelement(int i, int j, double a);

	//метод, извлекающий значение из нужной €чейки
	double getelement(int i, int j);

	//метод, осуществл€ющий перестановку строк
	void permutation(int i, int j);

	//метод, осуществл€ющий умножение стоки на скал€р
	void scalar(int i, double a);

	//метод, вычитающий помноженную на a j-ю строку из i-ой
	void substract(int i, int j, double a);

	//транспонирование
	Matrix T();

	//дефолтный конструктор строит матрицу 3x3
	Matrix();

	Matrix(int m, int n);

	Matrix(int m, int n, double k); //строит матрицу m x n с константой k на главной диагонали

	Matrix(vector<vector<double>>& C);

	void print();

	friend Matrix operator* (const Matrix& A, const Matrix& B);

	friend Vector operator* (const Matrix& A, const Vector& X);

};

class Vector : public Matrix
{
public:
	void setelement(int i, double a);

	double getelement(int i);

	Vector();

	Vector(int m);

	Vector(vector<double>& C);

	//friend Vector operator* (const RotationMatrix& A, const Vector& X);

	Vector& operator= (const Vector& v);

	friend Vector operator* (const Matrix& A, const Vector& X);

};

