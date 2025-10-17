//
// Created by Administrator on 2025/9/14.
//

#ifndef INC_4_5_H
#define INC_4_5_H

#include <iostream>

class Matrix {
    friend Matrix operator+ (const Matrix&, const Matrix&);
    friend Matrix operator* (const Matrix&, const Matrix&);
    friend std::ostream& operator<< (std::ostream& os, const Matrix &m);

public:
     Matrix(const double*);
     Matrix( double = 0., double = 0., double = 0., double = 0.,
             double = 0., double = 0., double = 0., double = 0.,
             double = 0., double = 0., double = 0., double = 0.,
             double = 0., double = 0., double = 0., double = 0.
            );

     int rows() const {return 4;}
     int cols() const {return 4;}

    std::ostream &print(std::ostream &) const;
    void operator+= (const Matrix&);
    double operator() (int row, int column) const {return _matrix[row][column];}
    double& operator()(int row, int column) {return _matrix[row][column];}


private:
     double _matrix[4][4]{};
};




#endif //INC_4_5_H
