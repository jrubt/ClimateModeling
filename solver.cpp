#include "climate_modeling.hpp"

double forwardEulerStep(double (*rhs)(double x), double t_start, double initialValue, double h)
{
    return initialValue + h * rhs(initialValue);
}

Solution solveODE(double (*rhs)(double x), double t_start, double t_end, double initialValue, double h)
{
    double t{t_start};
    uint64_t l{static_cast<uint64_t>((t_end - t_start) / h + 1.0)};
    std::vector<double> res_t(l);
    std::vector<double> res_x(l);
    res_t[0] = t_start;
    res_x[0] = initialValue;
    size_t i{1};
    while (t < t_end)
    {
        if (t_end - t < h)
        {
            h = t_end - t;
        }
        res_x[i] = forwardEulerStep(*rhs, t, res_x[i - 1], h);
        t += h;
        res_t[i] = t;
        i++;
    }

    return {res_t, res_x};
}
