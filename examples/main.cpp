#include "AvailDirsInterface.h"

MAKE_FUNC(
    Target, return 3 * v(0) * v(0) + 2 * v(1) * v(1) - 4 * v(0) * v(1) + 5 * v(0) + 2 * v(2) * v(2) + sin(v(2)) + v(3);
)

MAKE_FUNC(Cond1, return v(0) * v(0) + v(1) * v(1) + v(2) * v(2) - 4;)
MAKE_FUNC(Cond2, return v(0) + 2 * v(1) - v(2) - 2;)
MAKE_FUNC(Cond3, return -v(1) + 3 * v(3) - 5;)

int main() {
    Matrix A = {
        { 1.0, 1.0, 1.0, 1.0 }
    };
    std::vector<double> b = { 3.0 };

    try {
        auto solution = availdirs_solve<Target, Cond1, Cond2, Cond3>(A, b, true);
        print_vector(solution, "Solution:");
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
