#include "climate_modeling.hpp"

int main()
{
    double t0{0.0};
    double t_end{1e9};
    double h{1e7};
    double T0{320.0};
    Solution s{solveODE(ebm, t0, t_end, T0, h)};//solveODE([](double T)
      //                  { return (1.0 - kAlpha) * kA - kB * T * T * T * T; }, t0, t_end, T0, h)};
    writeSolution(s, "scenario_320K.txt");
    T0 = 280;
    s = solveODE(ebm, t0, t_end, T0, h);
    writeSolution(s, "scenario_280K.txt");
    T0 = 250;
    s = solveODE(ebm2, t0, t_end, T0, h);
    writeSolution(s, "scenario2_250K.txt");
    T0 = 270;
    s = solveODE(ebm2, t0, t_end, T0, h);
    writeSolution(s, "scenario2_270K.txt");
}