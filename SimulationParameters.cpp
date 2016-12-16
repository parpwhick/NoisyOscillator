/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "SimulationParameters.h"

void Parameters::set_dependent_parameters() {
    time_end = 100000 * dt;
    steps = std::ceil(time_end / dt) - 1;

    time_engine_start = time_end + dt;
    laser_initial_off = 0;
    laser_initial_on = time_end;
}
