#include "climate_modeling.hpp"

double ebm(double T)
{
    return (1.0 - k_alpha) * k_a - k_b * T * T * T * T;
}

double ebm2(double T)
{
    double alpha{0.5 - 0.2 * tanh((T - 265.0) / 10.0)};
    return (1.0 - alpha) * k_a - k_b * T * T * T * T;
}