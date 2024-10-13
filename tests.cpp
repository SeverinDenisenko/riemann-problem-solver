#include <fstream>
#include <iomanip>
#include <ios>
#include <string>
#include <vector>

#include "riemann.hpp"

using namespace riemann;

void test(gas_discontinuity discontinuity, real t, std::string name)
{
    solver solver(discontinuity);

    real x1 = -0.5;
    real x2 = 0.5;
    real n  = 1000;
    real dx = (x2 - x1) / n;

    auto solution = solver.solve();

    std::ofstream out(name);
    out << "x ro v p" << std::endl;

    for (integer i = 0; i < n; ++i) {
        real x          = x1 + dx * i;
        gas_state state = solution(t, x);
        out << std::setw(7) << std::fixed << std::setprecision(5) << x << " " << state.density << " " << state.velocity
            << " " << state.pressure << std::endl;
    }
}

int main()
{
    real gamma = 5.0 / 3.0;

    gas_discontinuity a { .left  = gas_state { .density = 1.0, .velocity = 0.0, .pressure = 3.0 },
                          .right = gas_state { .density = 1.0, .velocity = 0.0, .pressure = 1.0 },
                          .gamma = gamma };
    real a_t = 0.18;
    gas_discontinuity b { .left  = gas_state { .density = 1.0, .velocity = 1.0, .pressure = 3.0 },
                          .right = gas_state { .density = 1.0, .velocity = -1.0, .pressure = 1.0 },
                          .gamma = gamma };
    real b_t = 0.1;
    gas_discontinuity c { .left  = gas_state { .density = 1.0, .velocity = -0.1, .pressure = 1.0 },
                          .right = gas_state { .density = 1.0, .velocity = 0.2, .pressure = 1.0 },
                          .gamma = gamma };
    real c_t = 0.1;
    gas_discontinuity sod { .left  = gas_state { .density = 1.0, .velocity = 0.0, .pressure = 1.0 },
                          .right = gas_state { .density = 0.125, .velocity = 0.0, .pressure = 0.1 },
                          .gamma = 1.4 };
    real sod_t = 0.2;

    test(a, a_t, "first.dat");
    test(b, b_t, "second.dat");
    test(c, c_t, "third.dat");
    test(sod, sod_t, "sod.dat");

    return 0;
}
