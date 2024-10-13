#pragma once

#include <cstddef>
#include <functional>

namespace riemann {

using real    = double;
using integer = size_t;

struct gas_state {
    real density;
    real velocity;
    real pressure;
};

struct gas_discontinuity {
    gas_state left;
    gas_state right;
    real gamma;
};

using gas_state_function = std::function<gas_state(real t, real x)>;

class solver {
public:
    explicit solver(gas_discontinuity discontinuity);

    gas_state_function solve();

private:
    gas_state_function solve_configuarion_a();
    gas_state_function solve_configuarion_b();
    gas_state_function solve_configuarion_c();

    gas_state left_;
    gas_state right_;

    real gamma_;
    real plus_gamma_fraction_;
    real minus_gamma_fraction_;
    real gamma_fraction_;

    real left_sound_velocity_;
    real right_sound_velocity_;
};

}
