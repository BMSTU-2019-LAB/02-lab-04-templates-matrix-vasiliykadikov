//Copyright 2020 Vasya

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <stdio.h>
#include <math.h>
template <typename T>

class Matrix {
private:
    int n;
    int m;
    T** p;
public:
    Matrix(int n, int m);
    int get_rows() const;
    int get_columns() const;
    Matrix& operator()(Matrix& M);
    T* operator [](int i) const;
    Matrix& operator=(Matrix& M);
    Matrix operator+(Matrix& M);
    Matrix operator-(Matrix& M);
    Matrix operator*(Matrix& M);
    Matrix deletemn(Matrix& M, int row, int column);
    Matrix ad(Matrix& M, int row, int column);
    T det(Matrix& M);
    Matrix Inverse();
};
template<class T>
int Matrix<T>::get_rows() const {
    return n;
}
template<class T>
int Matrix<T>::get_columns() const {
    return m;
}
template<class T>
Matrix<T>::Matrix(int n, int m) {
    this -> p = reinterpret_cast<T**>(malloc(n * sizeof(T*)));
    this->n = n;
    this->m = m;
    for (int i = 0; i < n; i++) {
        p[i] = reinterpret_cast<T*>(malloc(m * sizeof(T)));
        for (int j = 0; j < m; j++) {
            p[i][j] = 0;
        }
    }
}
template<class T>
Matrix<T>& Matrix<T>::operator ()(Matrix& M) {
    for (int i = 0; i < M.get_rows; i++) {
        for (int j = 0; j < T.get_columns; j++) {
            K[i][j] = (*this)[i][j] = M[i][j];
        }
    }
}
template<class T>
T* Matrix<T>::operator [](int i) const {
    return p[i];
}
template<class T>
Matrix<T>& Matrix<T>::operator =(Matrix& M) {
    for (int i = 0; i < M.get_rows; i++) {
        for (int j = 0; j < T.get_columns; j++) {
            K[i][j] = (*this)[i][j] = M[i][j];
        }
    }
}
template<class T>
Matrix<T> Matrix<T>::operator+(Matrix& M) {
    if ((*this).get_columns != M.get_columns ||
     (*this).get_rows != M.get_rows) return 0;
    Matrix<T> K[M.get_rows][M.get_columns];
    for (int i = 0; i < M.get_rows; i++) {
        for (int j = 0; j < (*this).get_columns; j++) {
            K[i][j] = (*this)[i][j] + M[i][j];
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::operator-(Matrix& M) {
    if ((*this).get_columns != M.get_columns ||
    (*this).get_rows != M.get_rows) return 0;
    Matrix<T> K[M.get_rows][M.get_columns];
    for (int i = 0; i < M.get_rows; i++) {
        for (int j = 0; j < (*this).get_columns; j++) {
            K[i][j] = (*this)[i][j] - M[i][j];
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::operator*(Matrix& M) {
    if ((*this).get_rows != M.het_columns) return 0;
    Matrix<T> K[(*this).get_columns][M.get_rows];
    int n = 0;
    for (int i = 0; i < (*this.get_columns); i++) {
        for (int j = 0; j < M.get_rows; j++) {
            n = 0;
            k[i][j] = 0;
            for (int k; k < M.get_columns; k++) {
                K[i][j] = k[i][j](*this)[i][k] + M[k][j];
            }
        }
    }
    return K;
}

template<class T>
Matrix<T> Matrix<T>::deletemn(Matrix& M, int row, int column) {
    K[M.get_rows - 1][M.columns - 1]
        for (int i = 0; i < row - 1; i++) {
            for (int j = 0; j < column - 1; j++) {
                K[i][j] = M[i][j];
            }
            for (int k = column + 1; k < M.get_columns; k++) {
                K[i][k] = M[i][k];
            }
        }
    for (int i = row + 1; i < M.get_rows; i++) {
        for (int j = 0; j < column - 1; j++) {
            K[i][j] = M[i][j];
        }
        for (int k = column + 1; k < M.get_columns; k++) {
            K[i][k] = M[i][k];
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::ad(Matrix& M, int row, int column) {
    return (pow(-1, (row + column)) * M[row][column] * deletemn(M, row, column)
}
template<class T>
T Matrix<T>::det(Matrix& M) {
    T Det;
    if (M.get_rows == 1) {
        Det = M[0][0];
        return Det;
    }
    if (M.get_rows == 2) {
        Det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
        return Det;
    }
    for (int i = 0; i < M.ges_rows; i++) {
        Det += det(ad(M, 0, i));
    }
}
template<class T>
Matrix<T> Matrix<T>::Inverse() {
    if ((*this).get_columns != (*this).get_rows) return 0;
    K[(*this).get_rows][(*this).get_columns];
    for (int i = 0; i < (*this).get_rows; i++) {
        for (int j = 0; j < (*this).get_columns; i++) {
            K[i][j] = ad((*this), i, j);
        }
    }
    M[(*this).get_rows][(*this).get_columns];
    for (int i = 0; i < (*this).get_rows; i++) {
        for (int j = 0; j < (*this).get_columns; i++) {
            M[i][j] = K[j][i];
        }
    }
    T Det = det((*this));
    M = (1 / Det) * M;
    return M;
}

#endif // INCLUDE_MATRIX_HPP_
