#include "AvailDirsInterface.h"

MAKE_FUNC(Target, 3 * v(0) * v(0) + 2 * v(1) * v(1) - 4 * v(0) * v(1) + 5 * v(0) + 2 * v(2) * v(2) + sin(v(2)) + v(3);)
MAKE_FUNC(Cond1, v(0) * v(0) + v(1) * v(1) + v(2) * v(2) - 4;)
MAKE_FUNC(Cond2, v(0) + 2 * v(1) - v(2) - 2;)
MAKE_FUNC(Cond3, -v(1) + 3 * v(3) - 5;)

int main() {
    Matrix A = {
        { 1.0, 1.0, 1.0, 1.0 }
    };
    std::vector<double> b = { 3.0 };

    std::vector<double> solution = availdirs_solve<Target, Cond1, Cond2, Cond3>(A, b, true);
    print_vector(solution, "Solution:");
}
