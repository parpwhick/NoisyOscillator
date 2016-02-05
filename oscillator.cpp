/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   oscillator.cpp
 * Author: Dawid Crivelli
 *
 * Created on February 3, 2016, 11:01 AM
 */

#include "Simulation.hpp"
#include <vector>
#include <iostream>
#include "Utilities.h"
using namespace std;

void run_average(Simulation &sim, int N, double T) {
    double avg_en = 0.0;

    for (int i = 0; i < N; i++) {
        avg_en += sim.init_state(T);
    }
    avg_en /= N;

    cerr << "Average energy after " << N << " realizations:" << avg_en << endl;
}

int main(int argc, char** argv) {
    using namespace consts;
    Simulation single_ion;

    single_ion.init_state(0.080);
    single_ion.read_state(cerr);
    //    run_average(single_ion, 100000, 1);

    vec freqs;

    auto phys_par = Physical();
    auto local_par = Parameters{phys_par};

    for (double rf_freq = 15; rf_freq <= 100; rf_freq += 5) {
        double rf_omega = 2 * pi * rf_freq * MHz;
        
        single_ion.phys.RF_omega = rf_omega;
        single_ion.sim.dt = 0.05 / rf_omega;
        for (int axis = 0; axis < 3; axis++) {
            freqs[axis] = single_ion.trap_freq(axis);
        }
        cerr << (single_ion.phys.RF_omega / 2 / consts::pi) << " " << freqs[0] << " " <<
                freqs[1] << " " << freqs[2] << " " << endl;
    }
    return 0;
}

