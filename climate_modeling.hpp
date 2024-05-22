#include <cmath>
#include <vector>
#include <stddef.h>
#include <fstream>
#include <iostream>
#include <string>

const double kEpsilon{0.62};
const double kSigma{0.567e-7};
const double kS_0{1368.0};
const double kC{0.2e9};
const double kAlpha{0.3};
const double kA{kS_0 / (4 * kC)};
const double kB{kEpsilon * kSigma / kC};

struct Solution
{
    std::vector<double> t;
    std::vector<std::vector<double>> x;
};

/// @brief rhs  of energy-balance model (normalized)
/// @param T
/// @return
double ebm(double T);
/// @brief rhs of energy-balance model with temp-dependent alpha
/// @param T
/// @return
double ebm2(double T);
/// @brief
/// @param rhs
/// @param t_start
/// @param initialValue
/// @param h
/// @return
std::vector<double> predatorPreyEQ(double t, const std::vector<double> &v);

template <class F>
std::vector<double> forwardEulerStep(const F &rhs, const std::vector<double> &initialValue, double initialTime, double h)
{
    size_t dim{initialValue.size()};
    std::vector<double> res{rhs(initialTime, initialValue)};
    for (size_t i{0}; i < dim; ++i)
    {
        res[i] = initialValue[i] + h * res[i];
    }
    return res;
}

template <class F>
std::vector<double> rk4Step(const F &rhs, const std::vector<double> &initialValue, double initialTime, double h)
{
    size_t dim{initialValue.size()};
    std::vector<double> res{initialValue};
    // k1
    std::vector<double> k{rhs(initialTime, initialValue)};
    for (size_t i{0}; i < dim; ++i)
    {
        res[i] += h * k[i] / 6;
        k[i] = initialValue[i] + h * k[i] / 2;
    }
    // k2
    k = rhs(initialTime + h / 2, k);
    for (size_t i{0}; i < dim; ++i)
    {
        res[i] += h * k[i] / 3;
        k[i] = initialValue[i] + h * k[i] / 2;
    }
    // k3
    k = rhs(initialTime + h / 2, k);
    for (size_t i{0}; i < dim; ++i)
    {
        res[i] += h * k[i] / 3;
        k[i] = initialValue[i] + h * k[i];
    }
    // k4
    k = rhs(initialTime + h, k);
    for (size_t i{0}; i < dim; ++i)
    {
        res[i] += h * k[i] / 6;
    }
    return res;
}
/// @brief
/// @param rhs
/// @param t_start
/// @param t_end
/// @param initialValue
/// @param h
/// @return

template <class G>
Solution solveODE(const G &stepMethod, double t_start, double t_end, const std::vector<double> &initialValue, double h)
{
    double t{t_start};
    uint64_t l{static_cast<uint64_t>((t_end - t_start) / h + 1.0)};
    std::vector<double> res_t(l);
    std::vector<std::vector<double>> res_x(l);
    res_t[0] = t_start;
    res_x[0] = initialValue;
    size_t i{1};
    while (t < t_end)
    {
        if (t_end - t < h)
        {
            h = t_end - t;
        }
        res_x[i] = stepMethod(res_x[i - 1], res_t[i - 1], h);
        t += h;
        res_t[i] = t;
        i++;
    }

    return {res_t, res_x};
};
void printVector(const std::vector<double> &v);
void writeSolution(const Solution &s, const std::string &name);
