//Copyright 2020 Vasya

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <stdio.h>
#include <math.h>
template <typename T>

class Matrix {
    int n;
    int m;
    T** p;
public:
    ~Matrix();
    Matrix(const Matrix& M);
    Matrix(int m, int n);
    int get_rows() const;
    int get_columns() const;
    Matrix& operator()(Matrix<T> &M);
    T* operator [](size_t i) const;
    Matrix& operator=(Matrix<T> &M);
    Matrix operator+(Matrix<T> &M);
    Matrix operator-(Matrix<T> &M);
    Matrix operator*(Matrix<T> &M);
    Matrix deletemn(Matrix<T> &M, int row, int column);
    T det(Matrix<T> &M);
    Matrix Inverse();
    bool Matrix operator ==(const Matrix<T> &M);
    bool Matrix operator !=(const Matrix<T> &M);
};
template<class T>
Matrix<T>::~Matrix(){
    for (int i = 0; i < n; i++){
        free(p[i]);
    }
    free(p);
}
template<class T>
int Matrix<T>::get_rows() const {
    return n;
}
template<class T>
int Matrix<T>::get_columns() const {
    return m;
}
template<class T>
Matrix<T>::Matrix(int m, int n) {
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
Matrix<T>::Matrix(const Matrix<T> &M) {
 n = M.get_columns();
 m = M.get_rows();
 this -> p = reinterpret_cast<T**>(malloc(n * sizeof(T*)));
 for (int i = 0 ; i < M.get_rows(); i++){
  p[i] = reinterpret_cast<T*>(malloc(m * sizeof(T)));
  for (int j = 0 ; j < M.get_columns() ; j++){
   p[i][j] = M[i][j];
  }
 }
}
template<class T>
Matrix<T>& Matrix<T>::operator ()(Matrix<T> &M) {
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            (*this)[i][j] = M[i][j];
        }
    }
}
template<class T>
T* Matrix<T>::operator [](size_t i) const {
    return p[i];
}
template<class T>
Matrix<T>& Matrix<T>::operator =(Matrix<T> &M) {
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            (*this)[i][j] = M[i][j];
        }
    }
    return (*this);
}
template<class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &M) {
    if ((*this).get_columns() != M.get_columns() ||
     (*this).get_rows() != M.get_rows()) {
        Matrix<T> a(0, 0);
        return a;
    }
    Matrix<T> K(M.get_rows(), M.get_columns());
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            K[i][j] = (*this)[i][j] + M[i][j];
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> &M) {
    if ((*this).get_columns() != M.get_columns() ||
    (*this).get_rows() != M.get_rows()) {
        Matrix<T> a(0, 0);
        return a;
    }
    Matrix<T> K(M.get_rows(), M.get_columns());
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            K[i][j] = (*this)[i][j] - M[i][j];
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> &M) {
    if ((*this).get_rows() != M.get_columns()) {
        Matrix<T> a(0, 0);
        return a;
    Matrix<T> K((*this).get_columns(), M.get_rows());
    int n = 0;
    for (int i = 0; i < ((*this).get_columns()); i++) {
        for (int j = 0; j < M.get_rows(); j++) {
            n = 0;
            K[i][j] = 0;
            for (int k; k < M.get_columns(); k++) {
                K[i][j] = K[i][j] + (*this)[i][k] * M[k][j];
            }
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::deletemn(Matrix<T> &M, int row, int column) {
   Matrix<T> K(M.get_rows() - 1, M.get_columns() - 1);
   int a = row-1;
   int b = column-1;
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < b; j++) {
                K[i][j] = M[i][j];
            }
            for (int k = column + 1; k < M.get_columns(); k++) {
                K[i][k] = M[i][k];
            }
        }
    for (int i = row + 1; i < M.get_rows(); i++) {
        for (int j = 0; j < b; j++) {
            K[i][j] = M[i][j];
        }
        for (int k = column + 1; k < M.get_columns(); k++) {
            K[i][k] = M[i][k];
        }
    }
    return K;
}
template<class T>
T Matrix<T>::det(Matrix<T> &M) {
    T Det;
    if (M.get_rows() == 1) {
        Det = M[0][0];
        return Det;
    }
    if (M.get_rows() == 2) {
        Det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
        return Det;
    }
    for (int i = 0; i < M.ges_rows; i++) {
        Det += M[0][i] * pow(-1, i) * det(deletemn(M, 0, i));
   }
}
template<class T>
Matrix<T> Matrix<T>::Inverse() {
    if ((*this).get_columns() != (*this).get_rows()) {
        Matrix<T> a(0, 0);
        return a;
    }
    Matrix<T> K((*this).get_rows(), (*this).get_columns());
    for (int i = 0; i < ((*this).get_rows()); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            K[i][j] = pow(-1, i+j) * (*this)[i][j] * det(
                deletemn((*this), i, j));
        }
    }
    Matrix<T> M((*this).get_rows(), (*this).get_columns());
    for (int i = 0; i < (*this).get_rows(); i++) {
        for (int j = 0; j < (*this).get_columns(); i++) {
            M[i][j] = K[j][i];
        }
    }
    T Det = det((*this));
    double Detrev = 1 / Det;
    for (int i = 0; i < (*this).get_rows; i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            M[i][j] = Detrev * K[i][j];
        }
    }
    return M;
}
template<class T>
bool Matrix<T>::operator ==(const Matrix<T> &M) {
    if ((*this).get_rows != M.get_rows() &&
        (this).get_columns() != M.get_columns())
        return false;
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < M.get_columns(); j++) {
            if (*this)[i][j] != M[i][j] return false;
        }
    }
    return true;
}
template<class T>
bool Matrix<T>::operator !=(const Matrix<T> &M) {
    if ((*this).get_rows != M.get_rows() &&
        (this).get_columns() != M.get_columns())
        return true;
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < M.get_columns(); j++) {
            if (*this)[i][j] != M[i][j] return true;
        }
    }
    return false;
}
#endif // INCLUDE_MATRIX_HPP_
