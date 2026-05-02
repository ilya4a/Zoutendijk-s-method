#ifndef ZOITENDIJKMETHOD_AVAILDIRSINTERFACE_H
#define ZOITENDIJKMETHOD_AVAILDIRSINTERFACE_H

#include <memory>
#include <stdexcept>
#include <vector>

#include <autodiff/forward/real/eigen.hpp>
using namespace autodiff;

#include "AvailDirs.h"

using Matrix = std::vector<std::vector<double>>;
using Func = std::unique_ptr<FuncWrap>;

#define MAKE_FUNC(ClassName, Body)                              \
class ClassName : public FuncWrap {                          \
public:                                                      \
real f(const VectorXreal& v) const override { Body }     \
};

template<typename Target, typename... Constraints>
std::vector<double> availdirs_solve(
    const Matrix &A_eq,
    const std::vector<double> &b_eq,
    bool print_intermediate_results
) {
    AvailDirs solver;
    std::vector<Func> functions;
    functions.push_back(std::make_unique<Target>());

    auto add = [&](auto func) { functions.push_back(std::move(func)); };
    (add(std::make_unique<Constraints>()), ...);

    if (!solver.load_problem(std::move(functions), A_eq, b_eq)) {
        throw std::runtime_error("AvailDirs: cannot load problem");
    }
    return solver.solve_problem(print_intermediate_results);
}

template<typename T> void print_vector(const std::vector<T> &v, const std::string &m) {
    std::cout << m << " ";
    for (auto &i : v) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

#endif // ZOITENDIJKMETHOD_AVAILDIRSINTERFACE_H
