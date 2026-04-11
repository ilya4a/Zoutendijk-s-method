#include "ClpSimplex.hpp"
#include <iostream>



#include "ClpSimplex.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "diff.h"


class Func1 : public FuncWrap{
    public:
    autodiff::real f(autodiff::VectorXreal const &v) const override {
        return v(1) * v(2);
    };
};

int main() {

    Func1 f1;
    auto grad = f1.gradient_autodiff_lib({1, 2, 3});

    std::cout << "Gradient: ";
    for (double g : grad) std::cout << g << " ";
    std::cout << std::endl;

    // f.gradient_autodiff_lib(std::vector<double> {1, 2, 3});
}