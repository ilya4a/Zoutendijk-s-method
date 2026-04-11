//
// Created by ilya on 4/9/26.
//

#include "AvailDirs.h"




bool AvailDirs::check_func_condidions(std::vector<double> const &x, double& min_f_x) {
    double min = 1;
    for (int i = 1; i < functions.size(); i++) {
        double min_candidate = (*functions[i])(x);
        if (  min_candidate >= 0 ) return false;
        if (min_candidate < min) min = min_candidate;
    }
    min_f_x = min;
    return true;
}

std::vector<double> AvailDirs::solveUnderdeterminedEigen(const Matrix &A, const std::vector<double> &b) {
    int m = A.size();          // число уравнений
    int n = m > 0 ? A[0].size() : 0; // число переменных

    Eigen::MatrixXd A_eigen(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            A_eigen(i, j) = A[i][j];

    Eigen::VectorXd b_eigen = Eigen::Map<const Eigen::VectorXd>(b.data(), m);

    // Решаем (для недоопределённой системы даёт решение минимальной нормы)
    Eigen::VectorXd x = A_eigen.completeOrthogonalDecomposition().solve(b_eigen);

    return {x.data(), x.data() + n};
}

std::vector<int> AvailDirs::get_delta_conditions(std::vector<double> const &x) {
    std::vector<int> res;
    int j = 0;
    for (auto const& f: functions) {
        double v = (*f)(x);
        if (-delta <= v && v <= 0) {res.push_back(j);}
        j++;
    }
    return res;
}

double AvailDirs::calc_new_alpha(std::vector<double> const &x,
    std::vector<double> const &s,
    std::vector<int> const&nearly_to_active_cond_set) {
    double new_alpha = alpha;
    for (int k = 1; k < MAX_POW; k++) {
        std::vector<double> x_new(x.size());
        for (int i = 0; i < x.size(); i++) {x_new[i] = x[i] + alpha;}

        if ((*functions[0])(x_new) <= 0) {
            for () // CALC APLHA
        }
    }
}

std::vector<double> AvailDirs::solv_dirs_method(std::vector<double>& x0) {
    std::vector<int> nearly_to_active_cond_set = get_delta_conditions(x0);

    std::vector<std::vector<double>> gradients(functions.size(), std::vector<double>(nearly_to_active_cond_set.size(), 0));

    int j = 0;
    for(auto const& i: nearly_to_active_cond_set) {
        gradients[j] = (*(functions[i])).gradient_autodiff_lib(x0);
        j++;
    }

    std::vector<double> possible_dir(x0.size());
    double eng_out = 0.0;
    solve_subproblem(possible_dir, eng_out, gradients, A);
    eng = eng_out;
    if (eng < -delta) {
        // for (auto&p: possible_dir) { p = p*alpha;};


    }

}


bool AvailDirs::solve_subproblem(std::vector<double> &s_out,
                                 double &eng_out,
                                 const std::vector<std::vector<double> > &gradients,
                                 const std::vector<std::vector<double> > &A_eq
)
    {
        int n = gradients.empty() ? (A_eq.empty() ? 0 : A_eq[0].size()) : gradients[0].size();
        if (n == 0) return false;

        int n_ineq = gradients.size();
        int n_eq   = A_eq.size();
        int num_vars = n + 1;          // s1..sn, eng
        int num_rows = n_ineq + n_eq;

    // Сбор разреженной матрицы по столбцам
        ColSet columns(num_vars);
        calc_columns_clp_simplex(columns, gradients, A_eq, n, n_ineq, n_eq);

        // Формирование массивов CLP
        std::vector<int> start(num_vars + 1, 0);
        std::vector<int> index;
        std::vector<double> value;
        std::vector<int> length(num_vars, 0);
        int nnz = 0;
        for (int j = 0; j < num_vars; ++j) {
            start[j] = nnz;
            length[j] = static_cast<int>(columns[j].size());
            for (const auto& p : columns[j]) {
                index.push_back(p.first);
                value.push_back(p.second);
                // ограничение с весом value располагается в index строчке ограничений
                ++nnz;  // первая позиция начала ограничений для другой переменной
            }
        }
        start[num_vars] = nnz;

        // Границы переменных
        std::vector<double> colLower(num_vars, -1);
        std::vector<double> colUpper(num_vars,  1);
        colLower[n] = -COIN_DBL_MAX;
        colUpper[n] = COIN_DBL_MAX;

        // Целевая функция
        std::vector<double> objective(num_vars, 0.0);
        objective[n] = 1.0;   // минимизируем beta

        // Границы строк (ограничений)
        std::vector<double> rowLower(num_rows, -COIN_DBL_MAX);
        std::vector<double> rowUpper(num_rows,  COIN_DBL_MAX);
        for (int i = 0; i < n_ineq; ++i) {
            rowUpper[i] = 0.0;   // grad_i * s - beta <= 0
        }
        for (int k = 0; k < n_eq; ++k) {
            int idx = n_ineq + k;
            rowLower[idx] = 0.0;
            rowUpper[idx] = 0.0;   // A_k * s = 0
        }

        // Создание и решение модели
        ClpSimplex model;
        model.loadProblem(num_vars, num_rows,
                          start.data(), index.data(), value.data(), length.data(),
                          colLower.data(), colUpper.data(),
                          objective.data(),
                          rowLower.data(), rowUpper.data());

        model.primal();

        if (model.status() != 0) {
            std::cerr << "CLP error: status = " << model.status() << std::endl;
            return false;
        }

        const double* sol = model.primalColumnSolution();
        s_out.resize(n);
        for (int j = 0; j < n; ++j) s_out[j] = sol[j];
        eng_out = sol[n];

        return true;
    }

bool AvailDirs::load_problem(Functions functions, Matrix const &A, std::vector<double> b) {
    this->functions = std::move(functions);
    this->A = A;
    this->b = b;
    if (!check_problem()) return false;
    is_problem_exists = true;
    return true;
}

bool AvailDirs::check_problem() {
    return true;
}

AvailDirs::AvailDirs() {
    is_problem_exists = false;
}

void AvailDirs::calc_columns_clp_simplex(ColSet &columns,
                                         const std::vector<std::vector<double> > &gradients,
                                         const std::vector<std::vector<double> > &A_eq,
                                         const int n,
                                         const int n_ineq,
                                         const int n_eq)  {
    // Неравенства (строки 0..n_ineq-1)
    for (int i = 0; i < n_ineq; ++i) {
        const auto& grad = gradients[i];
        for (int j = 0; j < n; ++j) {
            if (std::fabs(grad[j]) > 1e-12) {
                columns[j].emplace_back(i, grad[j]);
                // чтобы не записывать нули в колумс, вставляем индекс где встречается этот не нулевой элемент
            }
        }
        // beta участвует во всех неравенствах с коэффициентом -1
        columns[n].emplace_back(i, -1.0);
    }

    // Равенства (строки n_ineq .. n_ineq+n_eq-1)
    for (int k = 0; k < n_eq; ++k) {
        const auto& row = A_eq[k];
        int row_idx = n_ineq + k;
        for (int j = 0; j < n; ++j) {
            if (std::fabs(row[j]) > 1e-12) {
                columns[j].emplace_back(row_idx, row[j]);
            }
        }
        // eng не входит в равенства
    }
}
