#include <cmath>
#include <vector>
#include <stddef.h>
#include <fstream>
#include <iostream>
#include <string>

const double k_epsilon{0.62};
const double k_sigma{0.567e-7};
const double k_s_0{1368.0};
const double k_c{0.2e9};
const double k_alpha{0.3};
const double k_a{k_s_0 / (4 * k_c)};
const double k_b{k_epsilon * k_sigma / k_c};

struct Solution
{
    std::vector<double> t;
    std::vector<double> x;
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
double forwardEulerStep(double (*rhs)(double x), double t_start, double initialValue, double h);
/// @brief
/// @param rhs
/// @param t_start
/// @param t_end
/// @param initialValue
/// @param h
/// @return
Solution solveODE(double (*rhs)(double x), double t_start, double t_end, double initialValue, double h);
void printVector(const std::vector<double> &v);
void writeSolution(const Solution &s, const std::string &name);
