//
//  matr_2.cpp
//  matr_2
//
//  Created by Dasha on 21.05.2022.
//

#include <iostream>
#include "matr_2.h"


ostream& operator<<(ostream& out, const Matrix& m)
{
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.columns; j++)
            out << m.matr[i][j] << " ";
        out << endl;
    }
    return out;
}

istream& operator>>(istream& in, Matrix& m)
{
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.columns; j++)
        {
            double element;
            in >> element;
            m.matr[i][j] = element;
        }
    }
    return in;
}

Matrix operator+(const Matrix& a, const Matrix& b)
{
    if((a.columns != b.columns) && (a.rows != b.rows))
        throw std::invalid_argument("It is impossible to add matrices of different sizes");
    Matrix res(a.rows, a.columns);
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < a.columns; j++)
        {
            res.matr[i][j] = a.matr[i][j] + b.matr[i][j];
        }
    }
    return res;
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
    if((a.columns != b.columns) && (a.rows != b.rows))
        throw std::invalid_argument("It is impossible to subtract matrices of different sizes");
    Matrix res(a.rows, a.columns);
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < a.columns; j++)
        {
            res.matr[i][j] = a.matr[i][j] - b.matr[i][j];
        }
    }
    return res;
}

Matrix operator*(double a, const Matrix& b)
{
    Matrix res(b.rows, b.columns);
    for(int i = 0; i < b.rows; i++)
    {
        for(int j = 0; j < b.columns; j++)
        {
            res.matr[i][j] = a * b.matr[i][j];
        }
    }
    return res;
}

Matrix operator*(double a, Matrix& b)
{
    Matrix res(b.rows, b.columns);
    for(int i = 0; i < b.rows; i++)
    {
        for(int j = 0; j < b.columns; j++)
        {
            res.matr[i][j] = a * b.matr[i][j];
        }
    }
    return res;
}


Matrix operator*(const Matrix& a, const Matrix& b)
{
    if((a.columns != b.rows) && (a.rows != b.columns))
        throw std::invalid_argument("It is impossible to subtract matrices of different sizes");
    Matrix res(a.rows, b.columns);
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < b.columns; j++)
        {
            res.matr[i][j] = 0;
            for(int k = 0; k < a.columns; k++)
            {
                res.matr[i][j] += a.matr[i][k] * b.matr[k][j];
            }
        }
    }
    return res;
}

Matrix Matrix::Adamar(const Matrix& a)
{
    if((a.columns != this->rows) && (a.rows != this->columns))
        throw std::invalid_argument("It is impossible to calculate the Hadamard product, since the matrices are of different sizes");
    Matrix res(a.rows,a.columns);
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < a.columns; j++)
        {
            res.matr[i][j] = a.matr[i][j] * this->matr[i][j];
        }
    }
    return res;
}


double Matrix::Trace()
{
    if(this->columns != this->rows)
        throw std::invalid_argument("It is impossible to calculate the trace of this matrix, since it is not square");
    double res = 0;
    for(int i = 0; i < this->rows; i++)
    {
        res += this->matr[i][i];
    }
    return res;
}

Matrix Matrix::Transpose()
{
    Matrix res(this->columns, this->rows);
    for(int i = 0; i < this->rows; i++)
    {
        for(int j = 0; j < this->columns; j++)
        {
            res.matr[j][i] = this->matr[i][j];
        }
    }
    return res;
}
    
double Matrix::Scalar(const Matrix& a)
{
    if((this->rows!=1) && (this->columns!=1))
        throw std::invalid_argument("You have not entered a vector");
    if((a.rows!=1) && (a.columns!=1))
        throw std::invalid_argument("You have not entered a vector");
    if((a.rows!=this->columns) && (a.rows!=this->rows))
        throw std::invalid_argument("It is impossible to calculate the scalar product of vectors of different sizes");
    double res = 0.0;
    if (a.rows == this->columns && a.columns ==1) return (*this * a).matr[0][0];
    if (a.columns == this->rows && a.rows ==1) return (a * *this).matr[0][0];
    if (a.rows == 1 && this->rows == 1)
    {
        for(int i = 0; i < a.columns; i++)
        {
            double elem = a.matr[0][i] * this->matr[0][i];
            res += elem;
        }
    }
    else if(a.columns == 1 && this->columns == 1){
        for(int i = 0; i < a.rows; i++)
        {
            double elem = a.matr[i][0] * this->matr[i][0];
            res += elem;
        }
    }
    return res;
}

double Matrix::MaximumRate()
{
    if((this->rows!=1) && (this->columns!=1))
        throw std::invalid_argument("You have not entered a vector");
    double max_elem = abs(this->matr[0][0]);
    if (this->rows == 1)
    {
        for(int i = 1; i < this->columns; i++)
        {
            if(max_elem < abs(this->matr[0][i]))
            {
                max_elem = abs(this->matr[0][i]);
            }
        }
    }
    else{
        for(int i = 1; i < this->rows; i++)
        {
            if(max_elem < abs(this->matr[i][0]))
            {
                max_elem = abs(this->matr[i][0]);
            }
        }
    }
    return max_elem;
}

double Matrix::EuclideanNorm()
{
    if((this->rows!=1) && (this->columns!=1))
        throw std::invalid_argument("You have not entered a vector");
    double res = 0;
    double elem;
    if (this->rows == 1)
    {
        for(int i = 0; i < this->columns; i++)
        {
            elem = pow(this->matr[0][i],2);
            res += elem;
        }
    }
    else{
        for(int i = 0; i < this->rows; i++)
        {
            elem = pow(this->matr[i][0],2);
            res += elem;
        }
    }
    return sqrt(res);
}

double Matrix::FrobeniusNorm()
{
    double res = 0;
    double sum = 0;
    for(int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->columns; j++)
        {
            sum = pow(this->matr[i][j],2);
            res += sum;
        }
    }
    
    return sqrt(res);
}

double Matrix::Angle(Matrix& a)
{
    if((this->rows!=1) && (this->columns!=1))
        throw std::invalid_argument("You have not entered a vector");
    if((a.rows!=1) && (a.columns!=1))
        throw std::invalid_argument("You have not entered a vector");
    return this->Scalar(a)/(a.EuclideanNorm() * this->EuclideanNorm());
}

Matrix Matrix::Tril(int n, int &change)
{
    bool flag = true;
    double max_elem = 0.0;
    int max_index = -1;
    for (int j = n; j < this->rows; j++) {
        if (flag) {
            max_elem = this->matr[j][n];
            max_index = j;
            flag = false;
        }
        else {
            if (this->matr[j][n] > max_elem) {
                max_elem = this->matr[j][n];
                max_index = j;
            }
        }
    }
    if (n != max_index) {
        swap(this->matr[n], this->matr[max_index]);
        change++;
    }
    for (int j = n + 1; j < this->rows; j++) {
        if (this->matr[n][n] == 0)
            break;
        double elem = this->matr[j][n] / this->matr[n][n];
        for (int k = 0; k < this->columns; k++) {
            this->matr[j][k] -= this->matr[n][k] * elem;
            if (abs(this->matr[j][k] - 0) < 0.00001)
                this->matr[j][k] = 0;
        }
    }
    for (int j = n - 1; j >= 0; j--) {
        if (this->matr[n][n] == 0)
            break;
        double elem = this->matr[j][n] / this->matr[n][n];
        for (int k = 0; k < this->columns; k++) {
            this->matr[j][k] -= this->matr[n][k] * elem;
            if (abs(this->matr[j][k] - 0) < 0.00001)
                this->matr[j][k] = 0;
        }
    }
    return *this;
}

double Matrix::Det()
{
    Matrix a = *this;
    if(a.columns != a.rows)
        throw std::invalid_argument("It is impossible to calculate the determinant of this matrix using the Gauss method, since it is not square");
    int change = 0;
    for (int i = 0; i < min(a.rows, a.columns); i++)
        a.Tril(i, change);
    double det = 1;
    for(int i = 0; i < a.rows; i++)
    {
        det *= a.matr[i][i];
    }
    //cout << a << endl;
    return pow(-1,change) * det;
}


double Matrix::Rank()
{
    Matrix a = *this;
    double sum = 0.0;
    int rank = a.rows;
    int change = 0;
    for (int i = 0; i < min(a.rows, a.columns); i++)
        a.Tril(i, change);
    for(int j = a.rows - 1; j >= 0; j-- )
    {
        for(int col =0; col < a.columns;col++)
        {
            sum += a.matr[j][col];
        }
        if(sum == 0)
        {
            rank--;
        }
    }
    cout << a << endl;
    return rank;
}


Matrix Matrix::Inverse()
{
    if(this->columns != this->rows)
        throw std::invalid_argument("It is impossible to calculate the inverse matrix of this matrix, since it is not square");
    if(this->Det() == 0) throw std::invalid_argument("it is impossible to calculate the inverse matrix if its determinant is 0");
    Matrix tmp(this->rows, 2 * this->columns);
    for(int i = 0; i < this->rows; i++)
    {
        for(int j = 0; j < this->columns; j++)
        {
            tmp.matr[i][j] = this->matr[i][j];
        }
    }
    for(int i = 0; i < this->rows; i++)
    {
        for(int j = this->columns; j < 2 * this->columns; j++)
        {
            if(i + this->columns== j)
            {
                tmp.matr[i][j] = 1;
            }
            else{
                tmp.matr[i][j] = 0;
            }
            
        }
    }
    //cout << tmp << endl;
    int change = 0;
    for (int i = 0; i < this->rows; i++)
        tmp.Tril(i, change);
    
    //cout << tmp << endl;
    for (int i = 0; i < this->rows; i++) {
        double elem = 0.0;
        for (int j = 0; j < 2*this->columns; j++) {
            if (i == j)
                elem = tmp.matr[i][j];
            if (j >= this->columns)
                tmp.matr[i][j] /= elem;
        }
    }
    //cout << tmp << endl;
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = columns; j < columns * 2; j++) {
            this->matr[i][j-columns] = tmp.matr[i][j];
        }
    }
    return *this;
}

Matrix Matrix::FromFile(string path)
{
    ifstream f(path);
    if(!f.is_open()) throw std::invalid_argument("Failed to open text file for reading");
    else
    {
        f >> *this;
    }
    return *this;
}

Matrix Matrix::FromBinFile()
{
    ifstream f("matr1.bin", ifstream::binary);
    int row;
    int column;
    if (!f.is_open()) throw std::invalid_argument("Failed to open binary file for reading");
    else
    {
        f.read((char*)&row, sizeof(int));
        f.read((char*)&column, sizeof(int));
        this->matr.resize(row, vector<double>(column));
        for (int i = 0; i < row; i++)
            for (int j = 0; j < column; j++)
                f.read((char*)&this->matr[i][j], sizeof(double));
    }
    cout << "Binary file read successfully" << endl;
    return *this;
}

void Matrix::ToFile()
{
    ofstream f("matr.txt");
    if(!f.is_open()) throw std::invalid_argument("Failed to open text file for writing");
    else
    {
        f << *this;
    }
}

void Matrix::ToBinFile()
{
    ofstream f("matr1.bin");
    if(!f.is_open()) throw std::invalid_argument("Failed to open binary file for writing");
    else
    {
        f.write((char*)&this->rows, sizeof(int));
        f.write((char*)&this->columns, sizeof(int));
        for(int i = 0; i < this->rows; i++)
        {
            for(int j = 0; j < this->columns; j++)
            {
                f.write((char*)&(this->matr[i][j]), sizeof(double));
            }
        }
    }
}


void PCA::centering()
{
    for(int i = 0; i < mat.getCols(); i++)
    {
        double average = 0;
        for(int j = 0; j < mat.getRows(); j++)
        {
            average += mat.getElem(j,i);
        }
        average /= mat.getRows();
        for(int k = 0; k < mat.getRows(); k++)
        {
            mat.getElem(k,i) -= average;
        }
    }
}

void PCA::print()
{
    cout << mat << endl;
}

void PCA::scaling()
{
    for(int i = 0; i < mat.getCols(); i++)
    {
        double average = 0;
        for(int j = 0; j < mat.getRows(); j++)
        {
            average += mat.getElem(j,i);
        }
        average /= mat.getRows();
        double s = 0;
        
        for(int k = 0; k < mat.getRows(); k++)
        {
            s += pow((mat.getElem(k,i) - average),2);
            mat.getElem(k,i) -= average;
        }
        
        s = sqrt(s/(mat.getRows()-1));
        
        for(int p = 0; p < mat.getRows(); p++)
        {
            mat.getElem(p,i) /= s;
        }
    }
}

tuple<Matrix, Matrix,Matrix> PCA::nipals(int PC)
{
    //this->centering();
    //this->scaling();
    E = mat;
    Matrix t(E.getRows(),1), p, t_old;
    vector<Matrix> P_, T_;
    
    for(int h = 0; h < PC; h++)
    {
        for(int j = 0; j < E.getRows(); j++)
        {
            t.getElem(j,0) = E.getElem(j,h);
        }
    
    do{
        p = (t.Transpose()*E).Transpose(); // считаем числитель
        auto q = (t.Transpose()*t).getElem(0,0);
        p = (1/q)*p;
        p = (1/p.EuclideanNorm())*p;
        
        t_old = t;
        t = E*p; // считаем числитель
        auto x = (p.Transpose()*p).getElem(0,0);
        t = (1/x)*t;
    } while((t_old-t).EuclideanNorm() > 0.00000001);
    E = E - t*p.Transpose();
    P_.push_back(p);
    T_.push_back(t);
    }
    
    Matrix T1(E.getRows(), T_.size()), P1(P_[0].getRows(), P_.size());

    for (int i = 0; i < E.getRows(); i++)
    {
        for (int j = 0; j < T_.size(); j++)
        {
            T1.getElem(i, j) = T_[j].getElem(i, 0);
            if (i < P_[0].getRows())
            {
                P1.getElem(i, j) = P_[j].getElem(i, 0);
            }
        }
    }
    T = T1;
    P = P1;
    return {T, P, E};
}

vector<double> PCA::scope()
{
    vector<double> res;
    for(int i = 0; i < T.getRows(); ++i)
    {
        Matrix t(1, T.getCols());
        for(int j = 0; j < T.getCols(); j++)
        {
            t.getElem(0,j) = T.getElem(i,j);
        }
        auto h = (t*((T.Transpose()*T).Inverse()))*t.Transpose();
        res.push_back(h.getElem(0,0));
    }
    
    for(int i = 0; i < res.size(); i++)
    {
        cout << res[i] << endl;
    }
    return res;
}

vector<double> PCA::deviation()
{
    vector<double> v;
    for(int i = 0; i < E.getRows(); i++)
    {
        double square = 0;
        for(int j = 0; j < E.getCols(); j++)
        {
            square += pow(E.getElem(i,j),2);
        }
        v.push_back(square);
    }
    
    for(int i = 0; i < v.size(); i++)
    {
        cout << v[i] << endl;
    }
    
    return v;
}

double PCA::V0(vector<double> v)
{
    double v0 = 0;
    for(int i = 0; i < v.size(); i++)
    {
        v0 += v[i];
    }
    auto res = v0/mat.getRows();
    return res;
}

double PCA::TRV(double res)
{
    res = res/mat.getCols();
    return res;
}

double PCA::ERV(double res)
{
    double denominator = 0;
    for(int i = 0; i < mat.getRows(); i++)
    {
        for(int j = 0; j < mat.getCols(); j++)
        {
            denominator += pow(mat.getElem(i,j),2);
        }
    }
    res = 1 - ((mat.getRows()*res)/denominator);
    return res;
}
