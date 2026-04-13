//
// Created by ilya on 4/9/26.
//

#include "AvailDirs.h"


template<typename T>
void print_vector(std::vector<T> const& v, std::string const& m) {
    std::cout << m << " ";
    for (auto& i : v) std::cout << i << " ";
    std::cout << std::endl;
}

std::vector<double> AvailDirs::solveUnderdeterminedEigen() {
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

    double fistappr = (is_first_approx) ? eng : 0.0;

    std::vector<int> res;

    for (int i = 1; i < functions.size(); i++) {

        double v = (*functions[i])(x) - fistappr;
        if (is_first_approx) {
            if (-delta <= v) {res.push_back(i);}  // ТАК КАК f_i(x) <= 0 не обязаны выполняться при каком-то Ax0 = b
        }else {
            if (-delta <= v && v <= 0) {res.push_back(i);}

        }
    }
    return res;
}



double AvailDirs::calc_new_alpha(std::vector<double> const &x, std::vector<double> const &s) {
    double new_alpha = 1.0;   // начинаем с 1, а не с alpha (0.5)
    for (int k = 1; k < MAX_POW; ++k) {
        std::vector<double> x_new(x.size());
        for (int i = 0; i < x.size(); ++i) x_new[i] = x[i] + new_alpha * s[i];

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
        // Проверка условия Армихо
        if ((*functions[0])(x_new) <= (*functions[0])(x) + 0.5 * eng * new_alpha) {
            return new_alpha;
        }
        new_alpha *= lambda;
    }
    return new_alpha;
}



bool AvailDirs::check_out_conditions(std::vector<double> const &x) {
    double maxf = -std::numeric_limits<double>::infinity();
    for (int i = 1; i < functions.size(); ++i) {
        double v = (*functions[i])(x);
        if (v > maxf) maxf = v;
    }

    return fabs(maxf - delta) <= EPS;
}


std::vector<double> AvailDirs::solv_dirs_method(std::vector<double>& x0) {

    int num = 0;

    while (!check_out_conditions(x0) && num < MAX_ITER ) {
        std::vector<int> nearly_to_active_cond_set = get_delta_conditions(x0);


        std::vector<std::vector<double>> gradients(nearly_to_active_cond_set.size() + 1, std::vector<double>(nearly_to_active_cond_set.size(), 0));

        gradients[0] = (*(functions[0])).gradient_autodiff_lib(x0);
        int j = 1;
        for(auto const& i: nearly_to_active_cond_set) {
            gradients[j] = (*(functions[i])).gradient_autodiff_lib(x0);
            j++;
        }

        std::vector<double> possible_dir(x0.size());
        double eng_out = 0.0;
        solve_subproblem(possible_dir, eng_out, gradients, A);
        eng = eng_out;
        if (eng < -delta) {
            double new_alpha = calc_new_alpha(x0, possible_dir);
            for (int i = 0; i < x0.size(); i++){x0[i] += new_alpha * possible_dir[i];}
        }else {
            if (delta < EPS) {
                break;
            }
            // иногда это условие помогает избежать зацикливания, иногда прерывает слишком рано
            // одна из больших проблем - не правильно работает условие выхода из этого цикла
            delta *= lambda;
        }
        num++;
    }

    std::cout << "NUM ITERATIONS: " << num << std::endl;

    if (num >= MAX_ITER) std::cerr << "limit in AvailDirs::fist_approx()" << std::endl;
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
        std::cout << "Initial point is feasible, no phase I needed.\n";
        return x0;
    }

    int num = 0;
    double temp_delta = delta;
    is_first_approx = true;

    while (eng >= 0 && num < MAX_ITER) {
        std::vector<int> nearly_to_active_cond_set = get_delta_conditions(x0);

        std::vector<std::vector<double>> gradients(nearly_to_active_cond_set.size(), std::vector<double>(x0.size(), 0));

        int j = 0;
        for(auto const& i: nearly_to_active_cond_set) {
            gradients[j] = (*(functions[i])).gradient_autodiff_lib(x0);
            j++;
        }

        std::vector<double> possible_dir(x0.size());
        double eng_out = 777;
        solve_subproblem(possible_dir, eng_out, gradients, A);

        eng = eng_out;

        if (eng < -delta) {
            double new_alpha = calc_new_alpha(x0, possible_dir); // отдельная функция для поиска шага
            for (int i = 0; i < x0.size(); ++i) {
                x0[i] += new_alpha * possible_dir[i];
            }
        } else {
            delta *= lambda;
        }
        num++;
    }

    delta = temp_delta;
    is_first_approx = false;
    // delta = 1e-1; // а оно тут нада?

    if (num >= MAX_ITER) std::cerr << "limit in AvailDirs::fist_approx()" << std::endl;
    return x0;
}


std::vector<double> AvailDirs::solve_problem() {


    std::vector<double> fist_approx = calc_fist_approx();

    std::cout << "first approx: ";
    for (auto& i : fist_approx) std::cout << i << " ";
    std::cout << std::endl;


    std::cout << "eng: " <<  eng <<  std::endl;

    return solv_dirs_method(fist_approx);
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

    //std::cerr << "n_ineq = " << n_ineq << ", n_eq = " << n_eq
    //      << ", n = " << n << std::endl;

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

//         // Создание и решение модели
         ClpSimplex model;
         model.loadProblem(num_vars, num_rows,
                           start.data(), index.data(), value.data(), length.data(),
                           colLower.data(), colUpper.data(),
                           objective.data(),
                           rowLower.data(), rowUpper.data());

        model.setLogLevel(0);
         model.primal();

         // std::cerr << "CLP status = " << model.status()
         //   << ", secondary = " << model.secondaryStatus() << std::endl;


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
    is_first_approx = false;

    alpha = 0.5;
    lambda = 0.5;
    delta = 0.5;
    eng = 1;
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