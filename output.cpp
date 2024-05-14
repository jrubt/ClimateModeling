#include "climate_modeling.hpp"

// template <typename T>
void printVector(const std::vector<double> &v)
{
    for (auto i : v)
    {
        std::cout << i << ' ';
    }
    std::cout << '\n';
}

void writeSolution(const Solution &s, const std::string &name)
{
    std::ofstream outf{name};
    uint64_t n{s.t.size()};
    for (size_t i{0}; i < n; ++i)
    {
        outf << s.t[i] << ' ' << s.x[i] << '\n';
    }
}