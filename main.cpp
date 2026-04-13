#include <iostream>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include "AvailDirs.h"
// #include "task.h"
#include "task2.h"
#include "task2.h"

using namespace autodiff;


int main() {
    // start2();
    start2();
}


//
// class ObjFunc2 : public FuncWrap {
// public:
//     real f(const VectorXreal& v) const override {
//         return (v(0) - 2) * (v(0) - 2) + v(1) * v(1)*2;
//     }
// };
//
// class ObjFunc3 : public FuncWrap {
// public:
//     real f(const VectorXreal& v) const override {
//         return v(0) * v(1);
//     }
// };
//
// // f0(x) = x1^2 + 4*x2^2
// class ObjFunc : public FuncWrap {
// public:
//     real f(const VectorXreal& v) const override {
//         return v(0) * v(0)*4.0 + v(1) * v(1);
//     }
// };
//
// // f1(x) = -x1 <= 0
// class C1 : public FuncWrap {
// public:
//     real f(const VectorXreal& v) const override {
//         return -v(0);
//     }
// };
//
// // f2(x) = -x2 <= 0
// class C2 : public FuncWrap {
// public:
//     real f(const VectorXreal& v) const override {
//         return -v(1);
//     }
// };

// int main() {
    // using Matrix = std::vector<std::vector<double>>;
    //
    // AvailDirs solver;
    //
    // // ВАЖНО:
    // // functions[0] — целевая функция
    // // functions[1..] — ограничения fi(x) <= 0
    // std::vector<std::unique_ptr<FuncWrap>> functions;
    // functions.push_back(std::make_unique<ObjFunc3>());
    // functions.push_back(std::make_unique<C1>());
    // functions.push_back(std::make_unique<C2>());
    //
    // // Единственное равенство: x1 + x2 = 1
    // Matrix A = {
    //     {1.0, 1.0}
    // };
    // std::vector<double> b = {5.0};
    //
    // if (!solver.load_problem(std::move(functions), A, b)) {
    //     std::cerr << "Ошибка загрузки задачи\n";
    //     return 1;
    // }
    //
    // std::vector<double> solution = solver.solve_problem();
    //
    // std::cout << "Найденное решение:\n";
    // std::cout << "x1 = " << solution[0] << ", x2 = " << solution[1] << "\n";
    //
    // VectorXreal v(2);
    // v << solution[0], solution[1];
    //
    // ObjFunc obj;
    // real val = obj.f(v);
    //
    // std::cout << "Значение целевой функции: " << val << "\n";
    //
    // std::cout << "Проверка ограничений:\n";
    // std::cout << "-x1 = " << -solution[0] << "\n";
    // std::cout << "-x2 = " << -solution[1] << "\n";
    // std::cout << "x1 + x2 = " << solution[0] + solution[1] << "\n";
    //
    // return 0;
// }