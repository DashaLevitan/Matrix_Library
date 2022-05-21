//
//  matr_2.hpp
//  matr_2
//
//  Created by Dasha on 21.05.2022.
//

#ifndef matr_2_
#define matr_2_

/* The classes below are exported */
#pragma GCC visibility push(default)

#include <iostream>
#include <vector>
#include <exception>
#include <cmath>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include<cmath>
#pragma once

using namespace std;

class Matrix
{
protected:
    vector<vector<double>> matr;
    int rows;
    int columns;
public:
    Matrix() : rows(0), columns(0) {}
    Matrix(int rows, int columns) : rows(rows), columns(columns)
    {
        matr.resize(rows, vector<double>(columns));
    }
    
    int getRows()
    {
        return this->rows;
    }
    
    int getCols()
    {
        return this->columns;
    }
    
    double &getElem(int i, int j)
    {
        return this->matr[i][j];
    }
    
    friend ostream& operator<<(ostream& out, const Matrix& m);
    friend istream& operator>>(istream& in, Matrix& m);
    friend Matrix operator+(const Matrix& a, const Matrix& b);
    friend Matrix operator-(const Matrix& a, const Matrix& b);
    friend Matrix operator*(double a, const Matrix& b);
    friend Matrix operator*(double a, Matrix& b);
    friend Matrix operator*(const Matrix& a, const Matrix& b);
    Matrix Adamar(const Matrix& a);
    double Trace();
    Matrix Transpose();
    double Scalar(const Matrix& a);
    double MaximumRate();
    double EuclideanNorm();
    double FrobeniusNorm();
    double Angle( Matrix& a);
    Matrix Tril(int,int&);
    double Rank();
    double Det();
    Matrix Inverse();
    Matrix FromFile(string);
    Matrix FromBinFile();
    void ToFile();
    void ToBinFile();
};

class Unit : public Matrix
{
public:
    Unit(int size): Matrix(size,size)
    {
        for(int i = 0; i < size; i++)
        {
            this->matr[i][i]= 1;
        }
    }
};

class Diagonal : public Matrix
{
public:
    Diagonal(int size): Matrix(size,size)
    {
        for(int i = 0; i < size; i++)
        {
            double elem;
            cin >> elem;
            this->matr[i][i]= elem;
        }
    }
};

class UpperTri : public Matrix
{
public:
    UpperTri(int size): Matrix(size,size)
    {
        for(int j = 0; j < size; j++)
        {
            for(int i = 0; i <= j; i++)
            {
                double elem;
                cin >> elem;
                this->matr[i][j] = elem;
            }
            
        }
    }
};

class LowerTri : public Matrix
{
public:
    LowerTri(int size): Matrix(size,size)
    {
        for(int i = 0; i < size; i++)
        {
            for(int j = 0; j <= i; j++)
            {
                double elem;
                cin >> elem;
                this->matr[i][j]= elem;
            }
            
        }
    }
};

class Symmetric : public Matrix
{
public:
    Symmetric(int size): Matrix(size,size)
    {
        for(int i = 0; i < size; i++)
        {
            for(int j = 0; j <= i; j++)
            {
                double elem;
                cin >> elem;
                this->matr[i][j] = this->matr[j][i]= elem;
            }
            
        }
    }
};

class PCA {
private:
    Matrix mat, T, P, E;

public:
    PCA(Matrix m) : mat(m) {}
    
    void centering();
    void print();
    void scaling();
    tuple<Matrix, Matrix,Matrix> nipals(int PC);
    vector<double> scope();
    vector<double> deviation();
    double TRV(double);
    double ERV(double);
    double V0(vector<double>);
};
#pragma GCC visibility pop
#endif
