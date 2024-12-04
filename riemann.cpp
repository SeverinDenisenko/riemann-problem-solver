#include "riemann.hpp"

#include <cmath>
#include <functional>
#include <utility>

namespace riemann {

namespace newton {
    using function = std::function<real(real)>;
    real derive_function(function f, real point, real step);
    real solve_newton(function f, real guess, integer iterations, real precision, real derive_step);

    real solve_newton(function f, real guess)
    {
        // Not based on any research, taken from my head
        return solve_newton(f, guess, 100, 1e-5, 1e-1);
    }
}

using gas_variable_function = std::function<real(real t, real x)>;

solver::solver(gas_discontinuity discontinuity)
    : left_(discontinuity.left)
    , right_(discontinuity.right)
    , gamma_(discontinuity.gamma)
{
}

gas_state_function solver::solve()
{
    left_sound_velocity_  = sqrt(gamma_ * left_.pressure / left_.density);
    right_sound_velocity_ = sqrt(gamma_ * right_.pressure / right_.density);

    minus_gamma_fraction_ = (gamma_ - 1) / (2 * gamma_);
    plus_gamma_fraction_  = (gamma_ + 1) / (2 * gamma_);
    gamma_fraction_       = (gamma_ - 1) / (gamma_ + 1);

    real velocity_a
        = 2 * left_sound_velocity_ / (gamma_ - 1) * (pow(right_.pressure / left_.pressure, minus_gamma_fraction_) - 1);
    real velocity_b = 1 / (right_.pressure * right_sound_velocity_) * (left_.pressure - right_.pressure)
        / sqrt(plus_gamma_fraction_ * left_.pressure / right_.pressure + minus_gamma_fraction_);

    bool configuration_a
        = velocity_a < (left_.velocity - right_.velocity) && (left_.velocity - right_.velocity) < velocity_b;
    bool configuration_b = (left_.velocity - right_.velocity) > velocity_b;
    bool configuration_c = (left_.velocity - right_.velocity) < velocity_a;

    if (configuration_a) {
        return solve_configuarion_a();
    }

    if (configuration_b) {
        return solve_configuarion_b();
    }

    if (configuration_c) {
        return solve_configuarion_c();
    }

    std::unreachable();
}

gas_state_function solver::solve_configuarion_a()
{
    using namespace newton;

    function velocity_star_left = [this](real pressure_star) -> real {
        return left_.velocity
            - 2 * left_sound_velocity_ / (gamma_ - 1)
            * (pow(pressure_star / left_.pressure, minus_gamma_fraction_) - 1);
    };

    function velocity_star_right = [this](real pressure_star) -> real {
        return right_.velocity
            + 1 / (right_.pressure * right_sound_velocity_) * (pressure_star - right_.pressure)
            / sqrt(plus_gamma_fraction_ * pressure_star / right_.pressure + minus_gamma_fraction_);
    };

    function pressure_solving_function = [velocity_star_left, velocity_star_right](real pressure_star) -> real {
        return velocity_star_left(pressure_star) - velocity_star_right(pressure_star);
    };

    real guess         = (left_.pressure + right_.pressure) / 2;
    real pressure_star = solve_newton(pressure_solving_function, guess);
    real velocity_star = velocity_star_left(pressure_star);

    real density_star_left  = left_.density * pow(pressure_star / left_.pressure, 1 / gamma_);
    real density_star_right = right_.density * ((gamma_ + 1) * pressure_star / right_.pressure + (gamma_ - 1))
        / ((gamma_ + 1) + (gamma_ - 1) * pressure_star / right_.pressure);

    real left_moving_wave_front_speed = left_.velocity - left_sound_velocity_;
    real left_moving_wave_back_speed
        = velocity_star - left_sound_velocity_ * pow(pressure_star / left_.pressure, minus_gamma_fraction_);

    gas_variable_function velocity_inside_left_moving_wave = [this](real t, real x) -> real {
        return left_.velocity * gamma_fraction_ + 2 / (gamma_ + 1) * (x / t + left_sound_velocity_);
    };

    gas_variable_function density_inside_left_moving_wave = [this](real t, real x) -> real {
        return left_.density
            * pow(2 / (gamma_ + 1) + gamma_fraction_ / left_sound_velocity_ * (left_.velocity - x / t),
                  2 / (gamma_ - 1));
    };

    gas_variable_function pressure_inside_left_moving_wave = [this](real t, real x) -> real {
        return left_.pressure
            * pow(2 / (gamma_ + 1) + gamma_fraction_ / left_sound_velocity_ * (left_.velocity - x / t),
                  2 * gamma_ / (gamma_ - 1));
    };

    real right_moving_wave_speed = right_.velocity
        + right_sound_velocity_
            * pow(plus_gamma_fraction_ * pressure_star / right_.pressure + minus_gamma_fraction_, 0.5);

    gas_state_function result_function = [this,
                                          left_moving_wave_front_speed,
                                          left_moving_wave_back_speed,
                                          velocity_star,
                                          right_moving_wave_speed,
                                          velocity_inside_left_moving_wave,
                                          density_inside_left_moving_wave,
                                          pressure_inside_left_moving_wave,
                                          density_star_right,
                                          density_star_left,
                                          pressure_star](real t, real x) -> gas_state {
        real left_moving_wave_front_position = left_moving_wave_front_speed * t;
        real left_moving_wave_back_position  = left_moving_wave_back_speed * t;
        real middle_discontinuity_position   = velocity_star * t;
        real right_moving_wave_position      = right_moving_wave_speed * t;

        if (x < left_moving_wave_front_position) {
            return gas_state {
                .density  = left_.density,
                .velocity = left_.velocity,
                .pressure = left_.pressure,
            };
        } else if (x < left_moving_wave_back_position) {
            return gas_state {
                .density  = density_inside_left_moving_wave(t, x),
                .velocity = velocity_inside_left_moving_wave(t, x),
                .pressure = pressure_inside_left_moving_wave(t, x),
            };
        } else if (x < middle_discontinuity_position) {
            return gas_state {
                .density  = density_star_left,
                .velocity = velocity_star,
                .pressure = pressure_star,
            };
        } else if (x < right_moving_wave_position) {
            return gas_state {
                .density  = density_star_right,
                .velocity = velocity_star,
                .pressure = pressure_star,
            };
        } else {
            return gas_state {
                .density  = right_.density,
                .velocity = right_.velocity,
                .pressure = right_.pressure,
            };
        }
    };

    return result_function;
}

gas_state_function solver::solve_configuarion_b()
{
    using namespace newton;

    function velocity_star_left = [this](real pressure_star) -> real {
        return left_.velocity
            - 1 / (left_.density * left_sound_velocity_) * (pressure_star - left_.pressure)
            / sqrt(plus_gamma_fraction_ * pressure_star / left_.pressure + minus_gamma_fraction_);
    };

    function velocity_star_right = [this](real pressure_star) -> real {
        return right_.velocity
            + 1 / (right_.density * right_sound_velocity_) * (pressure_star - right_.pressure)
            / sqrt(plus_gamma_fraction_ * pressure_star / right_.pressure + minus_gamma_fraction_);
    };

    function pressure_solving_function = [velocity_star_left, velocity_star_right](real pressure_star) -> real {
        return velocity_star_left(pressure_star) - velocity_star_right(pressure_star);
    };

    real guess         = (left_.pressure + right_.pressure) / 2;
    real pressure_star = solve_newton(pressure_solving_function, guess);
    real velocity_star = velocity_star_left(pressure_star);

    real density_star_left = left_.density * ((gamma_ + 1) * pressure_star / left_.pressure + (gamma_ - 1))
        / ((gamma_ + 1) + (gamma_ - 1) * pressure_star / left_.pressure);
    real density_star_right = right_.density * ((gamma_ + 1) * pressure_star / right_.pressure + (gamma_ - 1))
        / ((gamma_ + 1) + (gamma_ - 1) * pressure_star / right_.pressure);

    real left_moving_wave_speed = left_.velocity
        - left_sound_velocity_
            * pow(plus_gamma_fraction_ * pressure_star / left_.pressure + minus_gamma_fraction_, 0.5);

    real right_moving_wave_speed = right_.velocity
        + right_sound_velocity_
            * pow(plus_gamma_fraction_ * pressure_star / right_.pressure + minus_gamma_fraction_, 0.5);

    gas_state_function result_function = [this,
                                          velocity_star,
                                          left_moving_wave_speed,
                                          right_moving_wave_speed,
                                          density_star_right,
                                          density_star_left,
                                          pressure_star](real t, real x) -> gas_state {
        real left_moving_wave_position     = left_moving_wave_speed * t;
        real middle_discontinuity_position = velocity_star * t;
        real right_moving_wave_position    = right_moving_wave_speed * t;

        if (x < left_moving_wave_position) {
            return gas_state {
                .density  = left_.density,
                .velocity = left_.velocity,
                .pressure = left_.pressure,
            };
        } else if (x < middle_discontinuity_position) {
            return gas_state {
                .density  = density_star_left,
                .velocity = velocity_star,
                .pressure = pressure_star,
            };
        } else if (x < right_moving_wave_position) {
            return gas_state {
                .density  = density_star_right,
                .velocity = velocity_star,
                .pressure = pressure_star,
            };
        } else {
            return gas_state {
                .density  = right_.density,
                .velocity = right_.velocity,
                .pressure = right_.pressure,
            };
        }
    };

    return result_function;
}

gas_state_function solver::solve_configuarion_c()
{
    using namespace newton;

    function velocity_star_left = [this](real pressure_star) -> real {
        return left_.velocity
            - 2 * left_sound_velocity_ / (gamma_ - 1)
            * (pow(pressure_star / left_.pressure, minus_gamma_fraction_) - 1);
    };

    function velocity_star_right = [this](real pressure_star) -> real {
        return right_.velocity
            + 2 * right_sound_velocity_ / (gamma_ - 1)
            * (pow(pressure_star / right_.pressure, minus_gamma_fraction_) - 1);
    };

    function pressure_solving_function = [velocity_star_left, velocity_star_right](real pressure_star) -> real {
        return velocity_star_left(pressure_star) - velocity_star_right(pressure_star);
    };

    real guess         = (left_.pressure + right_.pressure) / 2;
    real pressure_star = solve_newton(pressure_solving_function, guess);
    real velocity_star = velocity_star_left(pressure_star);

    real density_star_left  = left_.density * pow(pressure_star / left_.pressure, 1 / gamma_);
    real density_star_right = right_.density * pow(pressure_star / right_.pressure, 1 / gamma_);

    real left_moving_wave_front_speed = left_.velocity - left_sound_velocity_;
    real left_moving_wave_back_speed
        = velocity_star - left_sound_velocity_ * pow(pressure_star / left_.pressure, minus_gamma_fraction_);

    real right_moving_wave_front_speed = right_.velocity + right_sound_velocity_;
    real right_moving_wave_back_speed
        = velocity_star + right_sound_velocity_ * pow(pressure_star / right_.pressure, minus_gamma_fraction_);

    gas_variable_function velocity_inside_left_moving_wave = [this](real t, real x) -> real {
        return left_.velocity * gamma_fraction_ + 2 / (gamma_ + 1) * (x / t + left_sound_velocity_);
    };

    gas_variable_function density_inside_left_moving_wave = [this](real t, real x) -> real {
        return left_.density
            * pow(2 / (gamma_ + 1) + gamma_fraction_ / left_sound_velocity_ * (left_.velocity - x / t),
                  2 / (gamma_ - 1));
    };

    gas_variable_function pressure_inside_left_moving_wave = [this](real t, real x) -> real {
        return left_.pressure
            * pow(2 / (gamma_ + 1) + gamma_fraction_ / left_sound_velocity_ * (left_.velocity - x / t),
                  2 * gamma_ / (gamma_ - 1));
    };

    gas_variable_function velocity_inside_right_moving_wave = [this](real t, real x) -> real {
        return right_.velocity * gamma_fraction_ + 2 / (gamma_ + 1) * (x / t - right_sound_velocity_);
    };

    gas_variable_function density_inside_right_moving_wave = [this](real t, real x) -> real {
        return right_.density
            * pow(2 / (gamma_ + 1) - gamma_fraction_ / right_sound_velocity_ * (right_.velocity - x / t),
                  2 / (gamma_ - 1));
    };

    gas_variable_function pressure_inside_right_moving_wave = [this](real t, real x) -> real {
        return right_.pressure
            * pow(2 / (gamma_ + 1) - gamma_fraction_ / right_sound_velocity_ * (right_.velocity - x / t),
                  2 * gamma_ / (gamma_ - 1));
    };

    gas_state_function result_function = [this,
                                          left_moving_wave_front_speed,
                                          left_moving_wave_back_speed,
                                          right_moving_wave_front_speed,
                                          right_moving_wave_back_speed,
                                          velocity_star,
                                          velocity_inside_left_moving_wave,
                                          density_inside_left_moving_wave,
                                          pressure_inside_left_moving_wave,
                                          velocity_inside_right_moving_wave,
                                          density_inside_right_moving_wave,
                                          pressure_inside_right_moving_wave,
                                          density_star_right,
                                          density_star_left,
                                          pressure_star](real t, real x) -> gas_state {
        real left_moving_wave_front_position  = left_moving_wave_front_speed * t;
        real left_moving_wave_back_position   = left_moving_wave_back_speed * t;
        real middle_discontinuity_position    = velocity_star * t;
        real right_moving_wave_front_position = right_moving_wave_front_speed * t;
        real right_moving_wave_back_position  = right_moving_wave_back_speed * t;

        if (x < left_moving_wave_front_position) {
            return gas_state {
                .density  = left_.density,
                .velocity = left_.velocity,
                .pressure = left_.pressure,
            };
        } else if (x < left_moving_wave_back_position) {
            return gas_state {
                .density  = density_inside_left_moving_wave(t, x),
                .velocity = velocity_inside_left_moving_wave(t, x),
                .pressure = pressure_inside_left_moving_wave(t, x),
            };
        } else if (x < middle_discontinuity_position) {
            return gas_state {
                .density  = density_star_left,
                .velocity = velocity_star,
                .pressure = pressure_star,
            };
        } else if (x < right_moving_wave_back_position) {
            return gas_state {
                .density  = density_star_right,
                .velocity = velocity_star,
                .pressure = pressure_star,
            };
        } else if (x < right_moving_wave_front_position) {
            return gas_state {
                .density  = density_inside_right_moving_wave(t, x),
                .velocity = velocity_inside_right_moving_wave(t, x),
                .pressure = pressure_inside_right_moving_wave(t, x),
            };
        } else {
            return gas_state {
                .density  = right_.density,
                .velocity = right_.velocity,
                .pressure = right_.pressure,
            };
        }
    };

    return result_function;
}

namespace newton {
    real derive_function(function f, real point, real step)
    {
        return (f(point + step) - f(point - step)) / (2 * step);
    }

    real solve_newton(function f, real guess, integer iterations, real precision, real derive_step)
    {
        real res = guess;

        for (integer i = 0; i < iterations; ++i) {
            real derivative = derive_function(f, guess, derive_step);

            res = guess - f(guess) / derivative;
            std::swap(res, guess);

            if (std::abs(res - guess) < precision) {
                break;
            }
        }

        return res;
    }
}

}
