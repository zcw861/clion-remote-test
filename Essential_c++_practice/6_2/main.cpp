/*//
// Created by Administrator on 2025/9/17.
//

#include "6.2.hpp"

int main() {
    const Matrix<int> m1(3,4,8);
    std::cout << m1 << std::endl;

    const Matrix<int> m2(m1);
    std::cout << m2 << std::endl;

    std::cout << m1 - m2 << std::endl;
    std::cout << m1 + m2 << std::endl;

    const Matrix<int> m3(4,6,3);
    std::cout << m1 * m3 << std::endl;

    const Matrix<int> m4(6,8,0,100);
    std::cout << m4 << std::endl;

    const Matrix<int> m5(5,4,7);
    std::cout << Matrix<int>::sum(m5,3) << std::endl;

    Matrix<int> m6(2,3,1);
    std::cout << Matrix<int>::Multiply(m6,9) << std::endl;

    const Matrix<int> m7(3,5,0,10);
    std::cout << m7 << "transpose" << std::endl << Matrix<int>::transpose(m7) << std::endl;

    const Matrix<int> m8(4,3,0,10);
    std::cout << m8 << "minor" << std::endl << Matrix<int>::minor(m8,0,1) << std::endl;

    const Matrix<int> m9(3,3,0,10);
    std::cout << m9 << "determinant" << std::endl << Matrix<int>::determinant(m9) << std::endl << std::endl;

    const Matrix<double> m10(5,5,0,10);
    std::cout << m10 << "inverse" << std::endl << Matrix<double>::inverse(m10) << std::endl;

    const Matrix<int> m11(2,4,0,10);
    const Matrix<int> m12(3,4,0,10);
    std::cout << "first Matrix:\n" << m11 << "second Matrix:\n" << m12 << "concatenate" << std::endl << Matrix<int>::concatenate(m11,m12,0) << std::endl;
    const Matrix<int> m13(4,2,0,10);
    const Matrix<int> m14(4,3,0,10);
    std::cout << "first Matrix:\n" << m13 << "second Matrix:\n" << m14 << "concatenate" << std::endl << Matrix<int>::concatenate(m13,m14,1) << std::endl;

    const Matrix<int> m15(5,6,0,10);
    std::cout << m15 << "ero_swap" << std::endl << Matrix<int>::ero_swap(m15,0,2) << std::endl;

    const Matrix<int> m16(4,5,0,10);
    std::cout << m16 << "ero_multiply" << std::endl << Matrix<int>::ero_multiply(m16,1,3) << std::endl;

    const Matrix<int> m17(6,7,0,10);
    std::cout << m17 << "ero_sum" << std::endl << Matrix<int>::ero_sum(m17,1,3,2) << std::endl;

    const Matrix<double> m18(7,6,0,10);
    std::cout << m18 << "upper_triangular" << std::endl << Matrix<double>::upper_triangular(m18) << std::endl;
}*/