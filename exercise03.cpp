#include "climate_modeling.hpp"

int main()
{
    double t0{0.0};
    double t_end{20.0};
    double h{0.001};
    std::vector<double> T0{2.0,1.0};
    Solution s{solveODE(predatorPreyEQ, t0, t_end, T0, h)};//solveODE([](double T)
      //                  { return (1.0 - kAlpha) * kA - kB * T * T * T * T; }, t0, t_end, T0, h)};
    writeSolution(s, "ex03_solution.txt");
    /*T0 = 280;
    s = solveODE(ebm, t0, t_end, T0, h);
    writeSolution(s, "scenario_280K.txt");
    T0 = 250;
    s = solveODE(ebm2, t0, t_end, T0, h);
    writeSolution(s, "scenario2_250K.txt");
    T0 = 270;
    s = solveODE(ebm2, t0, t_end, T0, h);
    writeSolution(s, "scenario2_270K.txt");*/
}