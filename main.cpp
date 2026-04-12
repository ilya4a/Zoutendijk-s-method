#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Dense>          // для MatrixXd
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

#include "AvailDirs.h"          // ваш заголовочный файл с классами

using namespace autodiff;

// Целевая функция: f(x) = x1 * x2
class MyFunc : public FuncWrap {
public:
    real f(const VectorXreal& v) const override {
        // Индексация с 0: v(0) = x1, v(1) = x2
        return v(0) * v(1);
    }
};

int main() {
    // 1. Создаём решатель
    using Matrix = std::vector<std::vector<double>>;
    AvailDirs solver;

    // 2. Формируем список функций (в данном случае только целевая)
    std::vector<std::unique_ptr<FuncWrap>> functions;
    functions.push_back(std::make_unique<MyFunc>());

    // 3. Линейные ограничения A * x <= b
    //    x1 >= 0  -> -x1 <= 0
    //    x2 >= 0  -> -x2 <= 0
    //    x1 + x2 <= 1
    Matrix A = {
        {-1.0,  0.0},   // первое ограничение: -x1 <= 0
        { 0.0, -1.0},   // второе: -x2 <= 0
        { 1.0,  1.0}    // третье: x1 + x2 <= 1
    };
    std::vector<double> b = {0.0, 0.0, 1.0};

    // 4. Загружаем задачу
    if (!solver.load_problem(std::move(functions), A, b)) {
        std::cerr << "Ошибка загрузки задачи" << std::endl;
        return 1;
    }

    // 5. Решаем
    std::vector<double> solution = solver.solve_problem();

    // 6. Выводим результат
    std::cout << "Найденное решение:\n";
    std::cout << "x1 = " << solution[0] << ", x2 = " << solution[1] << std::endl;

    // 7. Проверка допустимости
    double x1 = solution[0], x2 = solution[1];
    std::cout << "\nПроверка ограничений:\n";
    std::cout << "x1 >= 0 : " << (x1 >= -1e-6) << " (x1 = " << x1 << ")\n";
    std::cout << "x2 >= 0 : " << (x2 >= -1e-6) << " (x2 = " << x2 << ")\n";
    std::cout << "x1 + x2 <= 1 : " << (x1 + x2 <= 1 + 1e-6) << " (сумма = " << x1 + x2 << ")\n";

    // 8. Значение целевой функции
    VectorXreal v(2);
    v << x1, x2;
    real val = MyFunc().f(v);
    std::cout << "\nЗначение целевой функции: " << val << std::endl;

    return 0;
}