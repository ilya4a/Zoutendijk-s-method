#ifndef ZOITENDIJKMETHOD_DIFF_H
#define ZOITENDIJKMETHOD_DIFF_H

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <Eigen/Dense>
#include <vector>


class FuncWrap {

public:

    double operator()(std::vector<double> const& x) const {
        return f(Eigen::Map<const Eigen::VectorXd>(x.data(), x.size()).cast<autodiff::real>()).val();
    }

    autodiff::real operator()(autodiff::VectorXreal const& v) const {
        return f(v);
    }

    virtual autodiff::real f(autodiff::VectorXreal const& v) const = 0;

    std::vector<double> gradient_autodiff_lib(const std::vector<double>& x) const{

        autodiff::VectorXreal x_ad = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size()).cast<autodiff::real>();

        auto func_lambda = [this](autodiff::VectorXreal const& v) -> autodiff::real {
            return this->f(v);
        };
        Eigen::VectorXd grad = gradient(func_lambda, wrt(x_ad), at(x_ad));

        return {grad.data(), grad.data() + grad.size()};
    }

    virtual ~FuncWrap() {};
};


#endif //ZOITENDIJKMETHOD_DIFF_H