#ifndef ZOITENDIJKMETHOD_AVAILDIRS_H
#define ZOITENDIJKMETHOD_AVAILDIRS_H

#include <memory>
#include <vector>

#include "diff.h"

class AvailDirs {
    using ColSet = std::vector<std::vector<std::pair<int, double>>>;
    using Matrix = std::vector<std::vector<double>>;
    using Functions = std::vector<std::shared_ptr<FuncWrap>>;

    const int MAX_POW = 100;
    const int MAX_ITER = 100;
    const double EPS = 1e-3;

    Functions functions;
    Matrix A;
    std::vector<double> b;

    bool is_problem_exists;
    bool is_first_approx;

    double alpha;
    double lambda;

    double delta;
    double eng;

    void calc_columns_clp_simplex(
        ColSet &columns,
        const std::vector<std::vector<double>> &gradients,
        const std::vector<std::vector<double>> &A_eq,
        const int n,
        const int n_ineq,
        const int n_eq
    );

    bool solve_subproblem(
        std::vector<double> &s_out,
        double &eng_out,
        const std::vector<std::vector<double>> &gradients,
        const std::vector<std::vector<double>> &A_eq
    );

    std::vector<double> solveUnderdeterminedEigen();

    std::vector<int> get_delta_conditions(const std::vector<double> &x);

    double calc_new_alpha(const std::vector<double> &x, const std::vector<double> &s);

    bool check_out_conditions(const std::vector<double> &x);

    std::vector<double> solv_dirs_method(std::vector<double> &x0, bool print_intermediate_results);
    std::vector<double> calc_fist_approx();

    // static std::vector<double> calc_fist_approx(
    //     Functions &functions,
    //     Matrix &A,
    //     std::vector<double> &b,
    //     double alpha,
    //     double lambda
    // );

  public:
    AvailDirs();

    bool load_problem(Functions functions, const Matrix &A, std::vector<double> b);
    std::vector<double> solve_problem(bool print_intermediate_results);
};

#endif // ZOITENDIJKMETHOD_AVAILDIRS_H
