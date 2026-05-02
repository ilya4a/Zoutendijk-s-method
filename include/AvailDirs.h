//
// Created by ilya on 4/9/26.
//

#ifndef ZOITENDIJKMETHOD_AVAILDIRS_H
#define ZOITENDIJKMETHOD_AVAILDIRS_H


#include <iostream>
#include <vector>
#include <memory>

#include "diff.h"

class AvailDirs {
    using ColSet = std::vector<std::vector<std::pair<int, double>>>;
    using Matrix = std::vector<std::vector<double>>;
    using Functions = std::vector<std::unique_ptr<FuncWrap>>;

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
            ColSet& columns,
            const std::vector<std::vector<double>>& gradients,
            const std::vector<std::vector<double>>& A_eq,
            const int n,
            const int n_ineq,
            const int n_eq);

    bool solve_subproblem(
        std::vector<double>& s_out,
        double& eng_out,
        const std::vector<std::vector<double>>& gradients,
        const std::vector<std::vector<double>>& A_eq
        );


    std::vector<double> solveUnderdeterminedEigen();

    std::vector<int> get_delta_conditions(std::vector<double> const& x);

    double calc_new_alpha(std::vector<double> const &x,
        std::vector<double> const &s);

    bool check_out_conditions(std::vector<double> const& x);

    std::vector<double> solv_dirs_method(std::vector<double>& x0);

    double calc_new_alpha_fist_approx(const std::vector<double> &x, const std::vector<double> &s);


    std::vector<double> calc_fist_approx();

public:

    AvailDirs();

    bool load_problem(Functions functions, Matrix const& A, std::vector<double> b);
    std::vector<double> solve_problem();
};


#endif //ZOITENDIJKMETHOD_AVAILDIRS_H