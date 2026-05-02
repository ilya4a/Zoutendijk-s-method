#include "AvailDirs.h"

#include <iomanip>
#include <iostream>
#include <vector>

#include "ClpSimplex.hpp"

void print_iteration_header() {
    std::cout << std::setw(3) << "k"
              << "  " << std::setw(20) << "x"
              << " " << std::setw(10) << "delta"
              << " " << std::setw(10) << "eng"
              << " " << std::setw(10) << "phi_0" << std::endl;
}

void print_iteration(int k, const std::vector<double> &x, double delta, double eng, double phi0) {
    std::cout << std::setw(3) << k << "  ";
    std::cout << "(";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << std::setw(8) << std::fixed << std::setprecision(4) << x[i];
        if (i != x.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")";
    std::cout << " " << std::setw(10) << std::fixed << std::setprecision(4) << delta << " " << std::setw(10)
              << std::fixed << std::setprecision(4) << eng << " " << std::setw(10) << std::fixed << std::setprecision(4)
              << phi0 << std::endl;
}

template<typename T> void print_vector(const std::vector<T> &v, const std::string &m) {
    std::cout << m << " ";
    for (auto &i : v) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

std::vector<double> AvailDirs::solveUnderdeterminedEigen() {
    int m = A.size();
    int n = m > 0 ? A[0].size() : 0;

    Eigen::MatrixXd A_eigen(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            A_eigen(i, j) = A[i][j];
        }
    }

    Eigen::VectorXd b_eigen = Eigen::Map<const Eigen::VectorXd>(b.data(), m);

    Eigen::VectorXd x = A_eigen.completeOrthogonalDecomposition().solve(b_eigen);

    return { x.data(), x.data() + n };
}

std::vector<int> AvailDirs::get_delta_conditions(const std::vector<double> &x) {
    double fistappr = (is_first_approx) ? eng : 0.0;

    std::vector<int> res;

    for (int i = 1; i < functions.size(); i++) {
        double v = (*functions[i])(x) -fistappr;
        if (is_first_approx) {
            if (-delta <= v) {
                res.push_back(i);
            }
        } else {
            if (-delta <= v && v <= 0) {
                res.push_back(i);
            }
        }
    }
    return res;
}

double AvailDirs::calc_new_alpha(const std::vector<double> &x, const std::vector<double> &s) {
    double new_alpha = 1.0;
    for (int k = 1; k < MAX_POW; ++k) {
        std::vector<double> x_new(x.size());
        for (int i = 0; i < x.size(); ++i) {
            x_new[i] = x[i] + new_alpha * s[i];
        }

        bool feasible = true;
        for (int i = 1; i < functions.size(); ++i) {
            if ((*functions[i])(x_new) > EPS) {
                feasible = false;
                break;
            }
        }
        if (!feasible) {
            new_alpha *= lambda;
            continue;
        }

        if ((*functions[0])(x_new) <= (*functions[0])(x) + 0.5 * eng * new_alpha) {
            return new_alpha;
        }
        new_alpha *= lambda;
    }
    return new_alpha;
}

bool AvailDirs::check_out_conditions(const std::vector<double> &x) {
    double maxf = -std::numeric_limits<double>::infinity();
    for (int i = 1; i < functions.size(); ++i) {
        double v = (*functions[i])(x);
        if (v > maxf) {
            maxf = v;
        }
    }

    return fabs(maxf - delta) <= EPS;
}

std::vector<double> AvailDirs::solv_dirs_method(std::vector<double> &x0, bool print_intermediate_results) {
    int num = 0;
    if (print_intermediate_results) {
        print_iteration(num, x0, delta, eng, (*functions[0])(x0));
    }

    while (!check_out_conditions(x0) && num < MAX_ITER) {
        std::vector<int> nearly_to_active_cond_set = get_delta_conditions(x0);

        std::vector<std::vector<double>> gradients(
            nearly_to_active_cond_set.size() + 1,
            std::vector<double>(nearly_to_active_cond_set.size(), 0)
        );

        gradients[0] = (*(functions[0])).gradient_autodiff_lib(x0);
        int j = 1;
        for (const auto &i : nearly_to_active_cond_set) {
            gradients[j] = (*(functions[i])).gradient_autodiff_lib(x0);
            j++;
        }

        std::vector<double> possible_dir(x0.size());
        double eng_out = 0.0;
        solve_subproblem(possible_dir, eng_out, gradients, A);
        eng = eng_out;

        if (eng < -delta) {
            double new_alpha = calc_new_alpha(x0, possible_dir);
            for (int i = 0; i < x0.size(); i++) {
                x0[i] += new_alpha * possible_dir[i];
            }
        } else {
            if (delta < EPS) {
                break;
            }
            delta *= lambda;
        }
        num++;

        if (print_intermediate_results) {
            print_iteration(num, x0, delta, eng, (*functions[0])(x0));
        }
    }

    if (print_intermediate_results) {
        if (num >= MAX_ITER) {
            std::cerr << "limit in AvailDirs::solv_dirs()" << std::endl;
        }
        std::cout << "delta final: " << delta << std::endl;

        int j = 0;
        for (const auto &i : functions) {
            std::cout << "f" << j << " = " << (*i)(x0) << std::endl;
            j++;
        }
    }

    return x0;
}

std::vector<double> AvailDirs::calc_fist_approx() {
    std::vector<double> x0 = solveUnderdeterminedEigen();

    bool feasible = true;
    for (int i = 1; i < functions.size(); ++i) {
        if ((*functions[i])(x0) > EPS) {
            feasible = false;
            break;
        }
    }
    if (feasible) {
        eng = 0.0;
        is_first_approx = false;
        return x0;
    }

    int num = 0;
    double temp_delta = delta;
    is_first_approx = true;

    while (eng >= 0 && num < MAX_ITER) {
        std::vector<int> nearly_to_active_cond_set = get_delta_conditions(x0);

        std::vector<std::vector<double>> gradients(nearly_to_active_cond_set.size(), std::vector<double>(x0.size(), 0));

        int j = 0;
        for (const auto &i : nearly_to_active_cond_set) {
            gradients[j] = (*(functions[i])).gradient_autodiff_lib(x0);
            j++;
        }

        std::vector<double> possible_dir(x0.size());
        double eng_out = 777;
        solve_subproblem(possible_dir, eng, gradients, A);

        if (eng < -delta) {
            double new_alpha = calc_new_alpha(x0, possible_dir);
            for (int i = 0; i < x0.size(); ++i) {
                x0[i] += new_alpha * possible_dir[i];
            }
        } else {
            delta *= lambda;
        }
        num++;

        std::cout << num << " " << delta << " " << eng << " " << (*functions[0])(x0) << std::endl;
    }

    delta = temp_delta;
    is_first_approx = false;

    if (num >= MAX_ITER) {
        std::cerr << "limit in AvailDirs::fist_approx()" << std::endl;
    }
    return x0;
}

std::vector<double> AvailDirs::solve_problem(bool print_intermediate_results) {
    std::vector<double> fist_approx = calc_fist_approx();

    return solv_dirs_method(fist_approx, print_intermediate_results);
}

bool AvailDirs::solve_subproblem(
    std::vector<double> &s_out,
    double &eng_out,
    const std::vector<std::vector<double>> &gradients,
    const std::vector<std::vector<double>> &A_eq
) {
    int n = gradients.empty() ? (A_eq.empty() ? 0 : A_eq[0].size()) : gradients[0].size();
    if (n == 0) {
        return false;
    }

    int n_ineq = gradients.size();
    int n_eq = A_eq.size();
    int num_vars = n + 1;
    int num_rows = n_ineq + n_eq;

    ColSet columns(num_vars);
    calc_columns_clp_simplex(columns, gradients, A_eq, n, n_ineq, n_eq);

    std::vector<int> start(num_vars + 1, 0);
    std::vector<int> index;
    std::vector<double> value;
    std::vector<int> length(num_vars, 0);
    int nnz = 0;
    for (int j = 0; j < num_vars; ++j) {
        start[j] = nnz;
        length[j] = static_cast<int>(columns[j].size());
        for (const auto &p : columns[j]) {
            index.push_back(p.first);
            value.push_back(p.second);
            ++nnz;
        }
    }
    start[num_vars] = nnz;

    std::vector<double> colLower(num_vars, -1);
    std::vector<double> colUpper(num_vars, 1);
    colLower[n] = -COIN_DBL_MAX;
    colUpper[n] = COIN_DBL_MAX;

    std::vector<double> objective(num_vars, 0.0);
    objective[n] = 1.0;

    std::vector<double> rowLower(num_rows, -COIN_DBL_MAX);
    std::vector<double> rowUpper(num_rows, COIN_DBL_MAX);
    for (int i = 0; i < n_ineq; ++i) {
        rowUpper[i] = 0.0;
    }
    for (int k = 0; k < n_eq; ++k) {
        int idx = n_ineq + k;
        rowLower[idx] = 0.0;
        rowUpper[idx] = 0.0;
    }

    ClpSimplex model;
    model.loadProblem(
        num_vars,
        num_rows,
        start.data(),
        index.data(),
        value.data(),
        length.data(),
        colLower.data(),
        colUpper.data(),
        objective.data(),
        rowLower.data(),
        rowUpper.data()
    );

    model.setLogLevel(0);
    model.primal();

    if (model.status() != 0) {
        std::cerr << "CLP error: status = " << model.status() << std::endl;
        return false;
    }

    const double *sol = model.primalColumnSolution();
    s_out.resize(n);
    for (int j = 0; j < n; ++j) {
        s_out[j] = sol[j];
    }
    eng_out = sol[n];

    return true;
}

bool AvailDirs::load_problem(Functions functions, const Matrix &A, std::vector<double> b) {
    this->functions = std::move(functions);
    this->A = A;
    this->b = b;
    is_problem_exists = true;
    return true;
}

AvailDirs::AvailDirs() {
    is_problem_exists = false;
    is_first_approx = false;

    alpha = 0.5;
    lambda = 0.5;
    delta = 0.5;
    eng = 1;
}

void AvailDirs::calc_columns_clp_simplex(
    ColSet &columns,
    const std::vector<std::vector<double>> &gradients,
    const std::vector<std::vector<double>> &A_eq,
    const int n,
    const int n_ineq,
    const int n_eq
) {
    for (int i = 0; i < n_ineq; ++i) {
        const auto &grad = gradients[i];
        for (int j = 0; j < n; ++j) {
            if (std::fabs(grad[j]) > 1e-12) {
                columns[j].emplace_back(i, grad[j]);
            }
        }
        columns[n].emplace_back(i, -1.0);
    }

    for (int k = 0; k < n_eq; ++k) {
        const auto &row = A_eq[k];
        int row_idx = n_ineq + k;
        for (int j = 0; j < n; ++j) {
            if (std::fabs(row[j]) > 1e-12) {
                columns[j].emplace_back(row_idx, row[j]);
            }
        }
    }
}
