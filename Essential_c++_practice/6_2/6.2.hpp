//
// Created by Administrator on 2025/9/14.
//

#ifndef INC_6_2_HPP
#define INC_6_2_HPP

#include <vector>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <random>

template<typename elemType>
class Matrix {
    template<typename T>
    friend Matrix<T> operator+ (const Matrix<T>&, const Matrix<T>&);
    template<typename T>
    friend Matrix<T> operator- (const Matrix<T>&, const Matrix<T>&);
    template<typename T>
    friend Matrix<T> operator* (const Matrix<T>&, const Matrix<T>&);
    template<typename T>
    friend std::ostream& operator<< (std::ostream&, const Matrix<T>&);


public:
    Matrix(int rows, int columns, const elemType& initialValue = elemType());
    Matrix(const Matrix&);
    Matrix(int rows, int columns, elemType min, elemType max);   //随机矩阵
    ~Matrix() = default;

    int rows() const {return _rows;}
    int cols() const {return _cols;}

    std::ostream &print(std::ostream &) const;

    void operator+= (const Matrix&);
    void operator-=(const Matrix&);

    elemType operator() (int row, int column) const {return _matrix[row][column];}
    elemType& operator()(int row, int column) {return _matrix[row][column];}

    static void check_Add_and_Subtraction(const Matrix&, const Matrix&);
    static void checkMultiply(const Matrix&, const Matrix&);
    static Matrix<elemType> sum(const Matrix<elemType>&, elemType c);
    static Matrix<elemType> Multiply(const Matrix<elemType>&, elemType c);
    static Matrix<elemType> transpose(const Matrix<elemType>&);
    static Matrix<elemType> minor(const Matrix<elemType>&, int row, int col);
    static elemType determinant(const Matrix&);
    static Matrix<elemType> inverse(const Matrix<elemType>&);
    static Matrix<elemType> concatenate(const Matrix<elemType>&, const Matrix<elemType>&, int axis);
    static Matrix<elemType> ero_swap(const Matrix<elemType>&, int r1, int r2);
    static Matrix<elemType> ero_multiply(const Matrix<elemType>&, int row, elemType c);
    static Matrix<elemType> ero_sum(const Matrix<elemType>&, int r1, elemType c, int r2);
    static Matrix<elemType> upper_triangular(const Matrix<elemType>&);

private:
    std::vector<std::vector<elemType>> _matrix;
    int _rows;
    int _cols;
};

template<typename elemType>
Matrix<elemType>::Matrix(int rows, int columns, const elemType &initialValue):_rows(rows), _cols(columns), _matrix(rows, std::vector<elemType>(columns, initialValue)){}

template<typename elemType>
Matrix<elemType>::Matrix(const Matrix<elemType> &rhs):_rows(rhs._rows), _cols(rhs._cols), _matrix(rhs._matrix){}

template<typename elemType>
Matrix<elemType>::Matrix(int rows, int columns, elemType min, elemType max):_rows(rows), _cols(columns), _matrix(rows, std::vector<elemType>(columns)) {
    if (min > max) {
        throw std::invalid_argument("max must bigger min");
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min,max);

    for (int ix = 0; ix < _rows; ix++)
        for (int jx = 0; jx < _cols; jx++) {
            _matrix[ix][jx] = dis(gen);
        }
}

template<typename elemType>
void Matrix<elemType>::check_Add_and_Subtraction(const Matrix &m1, const Matrix &m2){
    if(m1._rows != m2._rows || m1._cols != m2._cols) {
        throw std::invalid_argument("Matrix dimensions must be equal");
    }
}

template<typename elemType>
Matrix<elemType> operator+ (const Matrix<elemType> &m1, const Matrix<elemType> &m2) {
    Matrix<elemType>::check_Add_and_Subtraction(m1,m2);
    Matrix<elemType> result(m1);
    result += m2;
    return result;
}

template<typename elemType>
void Matrix<elemType>::operator+= (const Matrix<elemType> &m) {
    Matrix<elemType>::check_Add_and_Subtraction(*this, m);
    for (int ix = 0; ix < m.rows(); ix++) {
        for (int jx = 0; jx <m.cols(); jx++) {
            _matrix[ix][jx] += m(ix,jx);
        }
    }
}

template<typename elemType>
Matrix<elemType> operator- (const Matrix<elemType> &m1, const Matrix<elemType> &m2) {
    Matrix<elemType>::check_Add_and_Subtraction(m1,m2);
    Matrix<elemType> result(m1);
    result -= m2;
    return result;
}

template<typename elemType>
void Matrix<elemType>::operator-=(const Matrix<elemType> &m) {
    check_Add_and_Subtraction(*this, m);
    for (int ix = 0; ix < m.rows(); ix++) {
        for (int jx = 0; jx <m.cols(); jx++) {
            _matrix[ix][jx] -= m(ix,jx);
        }
    }
}

template<typename elemType>
void Matrix<elemType>::checkMultiply(const Matrix &m1, const Matrix &m2) {
    if(m1.cols() != m2.rows()) {
        std::ostringstream oss;
        oss << "Matrix multiplication dimension mismatch: " << "A.cols(" << m1.cols() << "must equal B.rows" << m2.rows() << ")";
        throw std::invalid_argument(oss.str());
    }
}

template<typename elemType>
Matrix<elemType> operator* (const Matrix<elemType> &m1, const Matrix<elemType> &m2) {
    Matrix<elemType>::checkMultiply(m1,m2);
    Matrix<elemType> result(m1.rows(),m2.cols(),elemType());
    for (int ix = 0; ix < m1.rows(); ix++) {
        for (int jx = 0; jx < m2.cols(); jx++) {
            result(ix,jx) = 0;
            for (int kx = 0; kx < m1.cols(); kx++) {
                result(ix,jx) += m1(ix,kx) * m2(kx,jx);
            }
        }
    }
    return result;
}

template<typename elemType>
std::ostream &Matrix<elemType>::print(std::ostream &os) const {
    const double epsilon = 1e-10; // 视为零的阈值

    for (int ix = 0; ix < _rows; ++ix) {
        for (int jx = 0; jx < _cols; ++jx) {
            if (std::abs(_matrix[ix][jx]) < epsilon) {
                os << "0 ";
            } else {
                os << _matrix[ix][jx] << " ";
            }
        }
        os << std::endl;
    }
    return os;
}

template<typename elemType>
std::ostream &operator<< (std::ostream& os, const Matrix<elemType> &m) {
    return m.print(os);
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::sum(const Matrix &m, elemType c) {
    Matrix<elemType> result(m.rows(),m.cols(),0);
    for (int ix = 0; ix < m.rows(); ix++) {
        for (int jx = 0; jx < m.cols(); jx++) {
           result(ix,jx) = m(ix,jx) + c;
        }
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::Multiply(const Matrix &m, elemType c) {
    Matrix<elemType> result(m.rows(),m.cols(),0);
    for (int ix = 0; ix < m.rows(); ix++) {
        for (int jx = 0; jx < m.cols(); jx++) {
            result(ix,jx) = m(ix,jx) * c;
        }
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::transpose(const Matrix<elemType> &m) {
    Matrix<elemType> result(m.cols(),m.rows(),0);
    for (int ix = 0; ix < m.rows(); ++ix) {
        for (int jx = 0; jx < m.cols(); ++jx) {
            result(jx,ix) = m(ix,jx);
        }
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::minor(const Matrix<elemType> &m, int row, int col) {
    if(m.rows() < 2 || m.cols() < 2) {
        throw std::invalid_argument("Cannot compute minor of a matrix smaller than 2x2");
    }
    if(row+1 > m.rows() || row < 0) {
        throw std::invalid_argument("Row index out of range");
    }
    if(col+1 > m.cols() || col < 0){
        throw std::invalid_argument("Col index out of range");
    }

    Matrix<elemType> result(m.rows()-1, m.cols()-1);
    int r_row = 0;
    for (int ix = 0; ix < m.rows(); ++ix) {
        if(ix == row) continue;
        for (int jx = 0, r_col = 0; jx < m.cols(); ++jx) {
           if(jx == col) continue;
            result(r_row, r_col) =  m(ix,jx);
            r_col++;
        }
        r_row++;
    }
    return result;
}

template<typename elemType>
elemType Matrix<elemType>::determinant(const Matrix &m) {
    if(m.rows() != m.cols()) {
        throw std::invalid_argument("Determinant requires a square matrix");
    }

    if(m.rows() == 1 && m.cols() == 1)  return m(0,0);
    if(m.rows() == 2 && m.cols() == 2)  return m(0,0) * m(1,1) - m(0,1) * m(1,0);

    elemType result = 0;

    for(int i = 0; i < m.rows(); ++i) {
        Matrix<elemType> minorMat = minor(m,0,i);

        elemType sign = (i % 2 == 0) ? 1 : -1;

        result += sign * m(0,i) * determinant(minorMat);
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::inverse(const Matrix<elemType> &m) {
    if(m.rows() != m.cols()) {
        throw std::invalid_argument("Determinant requires a square matrix");
    }

    elemType det = Matrix<elemType>::determinant(m);      //m的行列式
    if(det == 0) {
        throw std::runtime_error("Matrix is singular (determinant is zero), no inverse exists");
    }

    if(m.rows() == 1 && m.cols() == 1) {
        Matrix<elemType> inv(1, 1);
        inv(0, 0) = elemType(1) / m(0, 0);
        return inv;
    }

    Matrix<elemType> cofactor(m.rows(), m.cols(),0);
    for(int ix = 0; ix < m.rows(); ++ix) {
        for(int jx = 0; jx < m.cols(); ++jx) {
            Matrix<elemType> minorMat = minor(m,ix,jx);
            int sign = ((ix + jx) % 2 == 0) ? 1 : -1;
            cofactor(ix,jx) = sign * determinant(minorMat);
        }
    }

    Matrix<elemType> adjugate = transpose(cofactor);
    Matrix<elemType> inv = Multiply(adjugate, elemType(1)/det);

    return inv;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::concatenate(const Matrix<elemType> &m1, const Matrix<elemType> &m2, int axis) {
    if(axis != 0 && axis != 1) {
        throw std::invalid_argument("implement this function so that it will concatenate and along the specified axis. (: on top of each other | : alongside each other)");
    }

    int r_rows = 0;
    int r_cols = 0;

    if(axis == 0) {             //垂直拼叠
        if(m1.cols() != m2.cols()) {
            throw std::invalid_argument("For axis=0, number of columns must match");
        }

        Matrix<elemType> result(m1.rows()+m2.rows(), m1.cols());

        for (int ix = 0; ix < m1.rows(); ++ix ) {
            for(int jx = 0; jx < m1.cols(); ++jx) {
                result(r_rows,r_cols) = m1(ix,jx);
                r_cols++;
            }
            r_cols = 0;
            r_rows++;
        }
        for (int ix = 0; ix < m2.rows(); ++ix ) {
            for(int jx = 0; jx < m2.cols(); ++jx) {
                result(r_rows,r_cols) = m2(ix,jx);
                r_cols++;
            }
            r_cols = 0;
            r_rows++;
        }
        return result;
    }
    if(axis == 1) {             //水平拼叠
        if(m1.rows() != m2.rows()) {
            throw std::invalid_argument("For axis=1, number of rows must match");
        }

        Matrix<elemType> result(m1.rows(), m1.cols()+m2.cols());

        for (int ix = 0; ix < m1.rows(); ++ix ) {
            for(int jx = 0; jx < m1.cols(); ++jx) {
                result(r_rows,r_cols) = m1(ix,jx);
                r_cols++;
            }
            r_cols = 0;
            r_rows++;
        }

        r_rows = 0;     //重置r_rows,为result添加下一个矩阵的行索引做准备
        r_cols = m1.cols();     //将r_cols设置为当前储存的列数

        for (int ix = 0; ix < m2.rows(); ++ix ) {
            for(int jx = 0; jx < m2.cols(); ++jx) {
                result(r_rows,r_cols) = m2(ix,jx);
                r_cols++;
            }
            r_cols = m1.cols();
            r_rows++;
        }
        return result;
    }

    throw std::logic_error("Unreachable: invalid axis");
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::ero_swap(const Matrix<elemType> &m, int r1, int r2) {
    if(r1+1 > m.rows() || r2+1 > m.rows()) {
        throw std::invalid_argument("rows out of range");
    }
    if(r1 < 0 || r2 < 0) {
        throw std::invalid_argument("rows must bigger 0");
    }

    Matrix<elemType> result(m);
        for (int jx = 0; jx < m.cols(); ++jx) {
            std::swap(result(r1,jx),result(r2,jx));
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::ero_multiply(const Matrix<elemType> &m, int row, elemType c) {
    if(row + 1 > m.rows()) {
        throw std::invalid_argument("rows out of range");
    }
    if(row < 0 ) {
        throw std::invalid_argument("rows must bigger 0");
    }

    Matrix<elemType> result(m);
        for (int jx = 0; jx < m.cols(); ++jx) {
            result(row,jx) = m(row, jx) * c;
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::ero_sum(const Matrix<elemType> &m, int r1, elemType c, int r2) {
    if(r1+1 > m.rows() || r2+1 > m.rows()) {
        throw std::invalid_argument("rows out of range");
    }
    if(r1 < 0 || r2 < 0) {
        throw std::invalid_argument("rows must bigger 0");
    }

    Matrix<elemType> result(m);
        for (int jx = 0; jx < m.cols(); ++jx) {
            result(r2,jx) += m(r1,jx) * c;
    }
    return result;
}

template<typename elemType>
Matrix<elemType> Matrix<elemType>::upper_triangular(const Matrix<elemType> &m) {
    Matrix<elemType> result(m);

    if (m.rows() == 0 || m.cols() == 0) {
        return result; // 空矩阵直接返回
    }

    size_t k = 0;  //当前处理的行（主元行）

    while(k < m.rows() && k < m.cols()) {
        size_t pivot_row = k;  //在第 k 列中，从第 k 行开始找非零主元

        while (pivot_row < m.rows() && result(pivot_row,k) == elemType(0)) {
            pivot_row++;
        }

        //如果没找到非零主元，说明这一列全为0，跳过
        if(pivot_row == m.rows()) {
            k++;
            continue;
        }

        if (pivot_row != k) {
            result = ero_swap(result, k, pivot_row);
        }

        //用第 k 行消去下面所有行在第 k 列的元素
         for(size_t i = k+1; i < m.rows(); ++i) {
             //如果第 i 行第 k 列已经是 0，跳过
             if(result(i,k) == 0) continue;

             //计算因子：让 result(i,k) 变成 0
             //所以：R_i = R_i - (result(i,k)/result(k,k)) * R_k-------------例如x + factor * p = 0，解的factor = -x / p
             double factor = -(result(i,k)) / (result(k,k));
             result = ero_sum(result,k,factor,i);
         }
        k++;    //处理下一行
    }
    return std::move(result);
}



#endif //INC_6_2_HPP
