//
// Created by Administrator on 2025/10/10.
//

//
// Created by Administrator on 2025/10/10.
//

#include "head.h"
#include <iostream>

int main() {
    int a = 5;
    int b = 3;
    int c = 1;

    std::cout << a << " + " << b << " = " << add(a, b) << std::endl;
    std::cout << a << " - " << b << " = " << sub(a, b) << std::endl;
    std::cout << a << " * " << b << " = " << mul(a, b) << std::endl;
    std::cout << a << " / " << b << " = " << division(a, b) << std::endl;

    return 0;
}