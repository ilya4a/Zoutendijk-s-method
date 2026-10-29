#include "AvailDirsInterface.h"

// Минимизируем t
MAKE_FUNC(Target, v(3);)

// x1 >= 0, x2 <= 0, x3 >= 0
MAKE_FUNC(Cond1, -v(0);)
MAKE_FUNC(Cond2,  v(1);)
MAKE_FUNC(Cond3, -v(2);)

// q1 = 1.10*x1 + 0.90*x2 + 1.10*x3
MAKE_FUNC(Cond4,
    1.10*v(0) + 0.90*v(1) + 1.10*v(2) - v(3);
)
MAKE_FUNC(Cond5,
   -1.10*v(0) - 0.90*v(1) - 1.10*v(2) - v(3);
)

// q2 = 1.40*x1 + 1.00*x2 + 0.80*x3
MAKE_FUNC(Cond6,
    1.40*v(0) + 1.00*v(1) + 0.80*v(2) - v(3);
)
MAKE_FUNC(Cond7,
   -1.40*v(0) - 1.00*v(1) - 0.80*v(2) - v(3);
)

// q3 = 0.80*x1 + 1.40*x2 + 1.20*x3
MAKE_FUNC(Cond8,
    0.80*v(0) + 1.40*v(1) + 1.20*v(2) - v(3);
)
MAKE_FUNC(Cond9,
   -0.80*v(0) - 1.40*v(1) - 1.20*v(2) - v(3);
)

int main() {
    // x1 - x2 + x3 = 1
    // Коэффициент при t равен нулю.
    Matrix A = {
        { 1.0, -1.0, 1.0, 0.0 }
    };

    std::vector<double> b = { 1.0 };

    std::vector<double> solution =
        availdirs_solve<
            Target,
            Cond1, Cond2, Cond3,
            Cond4, Cond5,
            Cond6, Cond7,
            Cond8, Cond9
        >(A, b, true);

    print_vector(solution, "Solution:");
}


// #include "AvailDirsInterface.h"
//
// MAKE_FUNC(Target, 3 * v(0) * v(0) + 2 * v(1) * v(1) - 4 * v(0) * v(1) + 5 * v(0) + 2 * v(2) * v(2) + sin(v(2)) + v(3);)
//
// MAKE_FUNC(Cond1, v(0) * v(0) + v(1) * v(1) + v(2) * v(2) - 4;)
// MAKE_FUNC(Cond2, v(0) + 2 * v(1) - v(2) - 2;)
// MAKE_FUNC(Cond3, -v(1) + 3 * v(3) - 5;)
//
//
// void f1() {
//     Matrix A = {
//         { 1.0, 1.0, 1.0, 1.0 }
//     };
//     std::vector<double> b = { 3.0 };
//
// }
//
// int main() {
//
//     f1();
//     // Matrix A = {
//     //     { 1.0, 1.0, 1.0, 1.0 }
//     // };
//     // std::vector<double> b = { 3.0 };
//     //
//     // std::vector<double> solution = availdirs_solve<Target, Cond1, Cond2, Cond3>(A, b, true);
//     // print_vector(solution, "Solution:");
// }
