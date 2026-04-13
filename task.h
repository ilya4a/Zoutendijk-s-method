//
// Created by ilya on 4/13/26.
//

#ifndef ZOITENDIJKMETHOD_TASK_H
#define ZOITENDIJKMETHOD_TASK_H
#include <memory>
#include <autodiff/forward/real/eigen.hpp>

#include "diff.h"

using namespace autodiff;

template<typename T>
void print_vector(std::vector<T> const& v, std::string const& m) {
    std::cout << m << " ";
    for (auto& i : v) std::cout << i << " ";
    std::cout << std::endl;
}

class Target : public FuncWrap {
public:
    real f(const VectorXreal& v) const override {
        return 3*v(0) * v(0) + 2*v(1) * v(1) -4*v(0)*v(1) +5*v(0) + 2*v(2)*v(2) + sin(v(2)) + v(3);// + exp(0.1*v(3));
    }
};


class Cond1 : public FuncWrap {
public:
    real f(const VectorXreal& v) const override {
        return v(0) * v(0) + v(1) * v(1) + v(2) * v(2) - 4;
    }
};

class Cond2 : public FuncWrap {
public:
    real f(const VectorXreal& v) const override {
        return v(0) +  2*v(1) - v(2) - 2;
    }
};


class Cond3 : public FuncWrap {
public:
    real f(const VectorXreal& v) const override {
        return -v(1) + 3*v(3) - 5;
    }
};


int start() {
    using Matrix = std::vector<std::vector<double>>;

    AvailDirs solver;

    std::vector<std::unique_ptr<FuncWrap>> functions;
    functions.push_back(std::make_unique<Target>());
    functions.push_back(std::make_unique<Cond1>());
    functions.push_back(std::make_unique<Cond2>());
    functions.push_back(std::make_unique<Cond3>());

    Matrix A = {
        {1.0, 1.0, 1.0, 1.0}
    };
    std::vector<double> b = {3.0};

    if (!solver.load_problem(std::move(functions), A, b)) {
        std::cerr << "Ошибка загрузки задачи\n";
        return 1;
    }

    std::vector<double> solution = solver.solve_problem();

    std::cout << "Найденное решение:\n";
    print_vector(solution, "");

    return 0;
}



#endif //ZOITENDIJKMETHOD_TASK_H