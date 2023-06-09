#pragma once
#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <iostream>

class Matrix
{
	std::vector<std::vector<float> > m;
	int rows, cols;
public:
	//matrix orf 0s, with specified (or default) row and column
	Matrix(int _rows = 4, int _cols = 4);
	Matrix(Vec3<float> v);	//store vector in matrix object with w
	int num_rows();	//gets rows
	int num_cols();	//get cols
	static Matrix identity(int dimensions); //return the identity
	std::vector<float>& operator[](const int i); //return a row of matrix
	//return matrix multiplication result
	Matrix operator* (const Matrix& a);
	Matrix transpose(); //return the transpose of matrix
	Matrix inverse(); //calculate the inverse of a matrix
	friend std::ostream& operator<<(std::ostream& s, Matrix& m); //print matrix
};

