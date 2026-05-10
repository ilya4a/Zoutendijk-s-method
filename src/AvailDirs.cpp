#include "AvailDirs.h"

#include <iomanip>
#include <iostream>
#include <vector>

#include "AvailDirsInterface.h"
#include "ClpSimplex.hpp"

void print_iteration_header() {
    std::cout << std::setw(3) << "k"
              << " " << std::setw(10) << "delta"
              << " " << std::setw(10) << "eta"
              << " " << std::setw(10) << "phi_0"
              << "  " << std::setw(20) << "x" << std::endl;
}

void print_iteration(int k, const std::vector<double> &x, double delta, double eng, double phi0) {
    std::cout << std::setw(3) << k << " ";

    std::cout << std::setw(10) << std::fixed << std::setprecision(4) << delta << " " << std::setw(10) << std::fixed
              << std::setprecision(4) << eng << " " << std::setw(10) << std::fixed << std::setprecision(4) << phi0
              << "  ";

    std::cout << "(";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << std::setw(8) << std::fixed << std::setprecision(4) << x[i];
        if (i != x.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;
}

std::vector<double> AvailDirs::solveUnderdeterminedEigen(Matrix &A, std::vector<double> b) {
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
    std::vector<int> res;

    for (int i = 1; i < functions.size(); i++) {
        double v = (*functions[i])(x);

        if (-delta <= v && v <= 0) {
            res.push_back(i);
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
            if ((*functions[i])(x_new) >= EPS) {
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
    return -1;
}

bool AvailDirs::check_out_conditions(const std::vector<double> &x) {
    if (fabs(fabs(eng) - EPS) > EPS) {
        return false;
    }

    double maxf = -std::numeric_limits<double>::infinity();
    for (int i = 1; i < functions.size(); ++i) {
        double v = (*functions[i])(x);

        if (fabs(v) < EPS) {
            continue;
        }

        if (v > maxf) {
            maxf = v;
        }
    }

    return delta < -maxf + EPS;
}

std::vector<double> AvailDirs::solv_dirs_method(
    std::vector<double> x0,
    bool print_intermediate_results,
    bool is_first_approx
) {
    if (print_intermediate_results) {
        print_iteration_header();
    }
    if (print_intermediate_results) {
        print_iteration(0, x0, delta, eng, (*functions[0])(x0));
    }

    int num = 0;

    while (!check_out_conditions(x0) && num < MAX_ITER) {
        if (is_first_approx) {
            if (x0.back() < -EPS) {
                return x0;
            }
        }

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

            if (new_alpha <= 0) {
                delta *= lambda;
            } else {
                for (int i = 0; i < x0.size(); i++) {
                    x0[i] += new_alpha * possible_dir[i];
                }
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
        int j = 0;
        for (const auto &i : functions) {
            std::cout << "f" << j << " = " << (*i)(x0) << std::endl;
            j++;
        }
    }

    return x0;
}

std::vector<double> AvailDirs::calc_fist_approx(
    const Functions &functions,
    Matrix A,
    std::vector<double> b,
    bool print_info
) {
    Functions functions_out;
    functions_out.reserve(functions.size());

    std::vector<double> x0 = solveUnderdeterminedEigen(A, b);

    double eng_start = -std::numeric_limits<double>::infinity();

    for (const auto &i : functions) {
        double v = (*i)(x0);
        if (v > eng_start) {
            eng_start = v;
        }
    }

    for (auto &i : A) {
        i.push_back(0);
    }
    b.push_back(0);
    x0.push_back(eng_start);

    size_t problem_dim = get_problem_dim(A);

    auto target = make_func_wrap([problem_dim](const autodiff::VectorXreal &v) { return v(problem_dim - 1); });

    functions_out.push_back(target);

    bool first = true;
    for (auto &f : functions) {
        if (first) {
            first = false;
            continue;
        }

        std::weak_ptr<FuncWrap> weak_f = f;
        auto cond = make_func_wrap([weak_f, problem_dim](const autodiff::VectorXreal &v) {
            if (auto shared = weak_f.lock()) {
                return shared->f(v) - v(problem_dim - 1);
            }
            throw std::runtime_error("weak_ptr in first_approx expired");
        });
        functions_out.push_back(cond);
    }

    AvailDirs solver;
    solver.load_problem(functions_out, A, b);

    return solver.solv_dirs_method(x0, print_info, true);
}

size_t AvailDirs::get_problem_dim(const Matrix &A) {
    if (A.size() == 0) {
        throw std::runtime_error("matrix A is empty");
    }
    int problem_dim = A[0].size();
    if (problem_dim == 0) {
        throw std::runtime_error("matrix A is empty");
    }

    return problem_dim;
}

std::vector<double> AvailDirs::solve_problem(bool print_intermediate_results) {
    if (print_intermediate_results) {
        std::cout << "initial approximation: " << std::endl;
    }

    std::vector<double> fist_approx = AvailDirs::calc_fist_approx(functions, A, b, print_intermediate_results);

    delta = -fist_approx.back() - EPS9;
    fist_approx.pop_back();

    if (print_intermediate_results) {
        std::cout << "Starting main optimization: " << std::endl;
    }

    return solv_dirs_method(fist_approx, print_intermediate_results, print_intermediate_results);
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
    this->functions = functions;
    this->A = A;
    this->b = b;
    is_problem_exists = true;

    dim = get_problem_dim(A);
    return true;
}

AvailDirs::AvailDirs() {
    is_problem_exists = false;

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
