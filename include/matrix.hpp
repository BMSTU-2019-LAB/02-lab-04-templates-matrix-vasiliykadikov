//Copyright 2020 Vasya

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <stdio.h>
#include <math.h>
#include <limits>
template <class T>

class Matrix {
    int n;
    int m;
    T** p;
public:
    ~Matrix();
    Matrix(const Matrix<T> &M);
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
    double det(Matrix<T> M);
    Matrix Inverse();
    template<class V>
    friend bool operator ==(const Matrix<V> &M, const Matrix<V> &m);
    template<class V>
    friend bool operator !=(const Matrix<V> &M, const Matrix<V> &m);
};
template<class T>
Matrix<T>::~Matrix(){
    for (int i = 0; i < n; i++){
        delete[] p[i];
    }
    delete[] p;
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
Matrix<T>::Matrix(int n, int m) {
    this->n = n;
    this->m = m;
    p = new T *[n];
    for (int i = 0; i < n; i++) {
        p[i] = new T[m];
    }
        for (int i = 0; i < n; i++) {
            for (int j=0; j < m; j++) {
            p[i][j] = 0;
        }
    }
}
template<class T>
Matrix<T>::Matrix(const Matrix<T> &M) {
 this->m = M.get_columns();
 this->n = M.get_rows();
 p = new T *[n];
 for (int i = 0 ; i < M.get_rows(); i++){
  p[i] = new T [m];
 }
  for (int i = 0 ; i < M.get_rows() ; i++){
      for (int j = 0; j < M.get_columns(); j++) {
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
 this->n = M.get_rows();
 this->m = M.get_columns();
 for (int i = 0 ; i < this->n ; i++){
  for (int j = 0 ; j < this->m ; j++){
   this->p[i][j] = M[i][j];
  }
 }
 return *this;
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
        for (int j = 0; j < M.get_columns(); j++) {
            K[i][j] = p[i][j] + M[i][j];
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
    if ((*this).get_columns() != M.get_rows()) {
        Matrix<T> a(0, 0);
        return a;
    }
    Matrix<T> K((*this).get_rows(), M.get_columns());
    for (int i = 0; i < K.get_rows(); i++) {
        for (int j = 0; j < K.get_columns(); j++) {
            for (int z=0; z < (*this).get_columns(); z++) {
                K[i][j] += (*this)[i][z] * M[z][j];
            }
        }
    }
    return K;
}
template<class T>
Matrix<T> Matrix<T>::deletemn(Matrix<T> &M, int row, int column) {
   Matrix<T> K(M.get_rows() - 1, M.get_columns() - 1);
   int a = row;
   int b = column;
        for (int i = 0; i < M.get_rows(); i++) {
            if (i != a) {
                for (int j = 0; j < M.get_columns(); j++) {
                    if (i != b) {
                        K[i][j] = M[i][j];
                    }
                }
            }
        }
    return K;
}
template<class T>
double Matrix<T>::det(Matrix<T> M) {
    double Det = 0;
    if (M.get_rows() == 1) {
        Det = M[0][0];
        return Det;
    }
    if (M.get_rows() == 2) {
        Det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
        return Det;
    }
    for (int i = 0; i < M.get_rows(); i++) {
     Det +=
         M[0][i] * pow(-1, i) * det(deletemn(*this, 0, i));
   }
    return Det;
}
template<class T>
Matrix<T> Matrix<T>::Inverse() {
    double Ad;
    Matrix<T> K((*this).get_rows(), (*this).get_columns());
    for (int i = 0; i < ((*this).get_rows()); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
         Ad = det(deletemn(*this, i, j));
         K[i][j] = pow(-1, i+j) * Ad;
        }
    }
    Matrix<T> M((*this).get_rows(), (*this).get_columns());
    for (int i = 0; i < K.get_rows(); i++) {
        for (int j = 0; j < K.get_columns(); j++) {
            M[i][j] = K[j][i];
        }
    }
    double Det = det(*this);
    double Detrev;
    for (int i = 0; i < (*this).get_rows(); i++) {
        for (int j = 0; j < (*this).get_columns(); j++) {
            Detrev = M[i][j] / Det;
            K[i][j] = Detrev;
        }
    }
    return K;
}
template<class T>
bool operator ==(const Matrix<T> &M, const Matrix<T> &m) {
    if (m.get_rows() != M.get_rows() &&
        m.get_columns() != M.get_columns())
        return false;
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < M.get_columns(); j++) {
            if (M[i][j] != m[i][j]) return false;
        }
    }
    return true;
}
template<>
bool operator ==(const Matrix<double> &M, const Matrix<double> &m) {
    if (m.get_rows() != M.get_rows() &&
        m.get_columns() != M.get_columns())
        return false;
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < M.get_columns(); j++) {
            if (abs(m[i][j] - M[i][j])
            > std::numeric_limits<double>::epsilon()) {
    return false;
            }
        }
    }
    return true;
}
template<>
bool operator ==(const Matrix<float> &M, const Matrix<float> &m) {
    if (m.get_rows() != M.get_rows() &&
        m.get_columns() != M.get_columns())
        return false;
    for (int i = 0; i < M.get_rows(); i++) {
        for (int j = 0; j < M.get_columns(); j++) {
            if (abs(m[i][j] - M[i][j])
            > std::numeric_limits<float>::epsilon()) {
    return false;
            }
        }
    }
    return true;
}
template<class T>
bool operator !=(const Matrix<T> &M, const Matrix<T> &m) {
    return !(M == m);
}
#endif // INCLUDE_MATRIX_HPP_
