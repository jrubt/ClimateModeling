// #include <cstdint>
#include <cmath>
#include <vector>
#include <stddef.h>
#include <fstream>
#include <iostream>
#include <string>

const double EPSILON{0.62};
const double SIGMA{0.567e-7};
const double S_0{1368.0};
const double C{0.2e9};
const double ALPHA{0.3};
const double A{S_0/(4*C)};
const double B{EPSILON*SIGMA/C};

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
void printVector(std::vector<double> v);
void writeSolution(Solution s, std::string name);
