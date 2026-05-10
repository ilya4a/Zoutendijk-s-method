#ifndef ZOITENDIJKMETHOD_DIFF_H
#define ZOITENDIJKMETHOD_DIFF_H

#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <vector>

class FuncWrap {
  public:
    double operator()(const std::vector<double> &x) const {
        return f(Eigen::Map<const Eigen::VectorXd>(x.data(), x.size()).cast<autodiff::real>()).val();
    }

    autodiff::real operator()(const autodiff::VectorXreal &v) const { return f(v); }

    virtual autodiff::real f(const autodiff::VectorXreal &v) const = 0;

    std::vector<double> gradient_autodiff_lib(const std::vector<double> &x) const {
        autodiff::VectorXreal x_ad = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size()).cast<autodiff::real>();

        auto func_lambda = [this](const autodiff::VectorXreal &v) -> autodiff::real { return this->f(v); };
        Eigen::VectorXd grad = gradient(func_lambda, wrt(x_ad), at(x_ad));

        return { grad.data(), grad.data() + grad.size() };
    }

    virtual ~FuncWrap() { };
};

template<typename F> class LambdaFuncWrap : public FuncWrap {
    F lambda;

  public:
    explicit LambdaFuncWrap(F f) : lambda(std::move(f)) { }

    autodiff::real f(const autodiff::VectorXreal &v) const override { return lambda(v); }
};

template<typename F> auto make_func_wrap(F &&f) {
    return std::make_shared<LambdaFuncWrap<std::decay_t<F>>>(std::forward<F>(f));
}

#endif // ZOITENDIJKMETHOD_DIFF_H
