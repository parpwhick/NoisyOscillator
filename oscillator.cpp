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

void scan_frequency(){
    using namespace consts;
    #pragma omp parallel for schedule(dynamic)
    for (int rf_freq = 10; rf_freq <= 100; rf_freq += 5) {
        vec freqs;
        Simulation single_ion;
        double rf_omega = 2 * pi * rf_freq * MHz;
        
        single_ion.phys.RF_omega = rf_omega;
        single_ion.sim.dt = 0.05 / rf_omega;
        for (int axis = 0; axis < 3; axis++) {
            freqs[axis] = single_ion.trap_freq(axis);
        }
        
        #pragma omp critical 
        cout << rf_freq << " " << freqs[0] << " " <<
                freqs[1] << " " << freqs[2] << " " << endl;
    }
}

void scan_frequency_comp(){
    using namespace consts;
    #pragma omp parallel for schedule(dynamic)
    for (int rf_freq = 1; rf_freq <= 100; rf_freq += 1) {
        vec freqs;
        Simulation single_ion;
        double rf_omega = 2 * pi * rf_freq * MHz;
        
        single_ion.phys.RF_omega = rf_omega;
        single_ion.sim.dt = 0.05 / rf_omega / rf_freq;
        single_ion.phys.RF_amplitude = 10*rf_freq; 
        for (int axis = 0; axis < 3; axis++) {
            freqs[axis] = single_ion.trap_freq(axis);
        }
        
        #pragma omp critical 
        cout << rf_freq << " " << freqs[0] << " " <<
                freqs[1] << " " << freqs[2] << " " << endl;
    }
}

void scan_amplitude(){
    using namespace consts;
    #pragma omp parallel for schedule(dynamic)
    for (int rf_ampl = 1; rf_ampl <= 100; rf_ampl += 5) {
        vec freqs;
        Simulation single_ion;
        double rf_omega = 2 * pi * 10 * MHz;
        
        single_ion.phys.RF_omega = rf_omega;
        single_ion.sim.dt = 0.05 / rf_omega;
        single_ion.phys.RF_amplitude = rf_ampl;
        
        for (int axis = 0; axis < 3; axis++) {
            freqs[axis] = single_ion.trap_freq(axis);
        }
        
        #pragma omp critical 
        cout << rf_ampl << " " << freqs[0] << " " <<
                freqs[1] << " " << freqs[2] << " " << endl;
    }
}

/* Perform scan at single frequency, but with different timesteps, 
 * to monitor the converge of the numerical method 
 */
void scan_accuracy(){
    using namespace consts;
    #pragma omp parallel for schedule(dynamic)
    for (int dt = 1; dt <= 100; dt += 1) {
        vec freqs;
        vec pos;
        Simulation single_ion;
        double rf_omega = 2 * pi * 30 * MHz;
        
        single_ion.phys.RF_omega = rf_omega;
        single_ion.sim.dt = 0.01 * dt / rf_omega;
        single_ion.phys.RF_amplitude = 30;
        for (int axis = 0; axis < 3; axis++) {
            freqs[axis] = single_ion.trap_freq(axis);
            pos[axis] = single_ion.x[axis];
        }
        
        #pragma omp critical 
        cout << single_ion.sim.dt << " " << freqs[0] << " " <<
                freqs[1] << " " << freqs[2] << " " 
            << pos[0] << " " << pos[1] << " " << pos[2] 
            << endl;
    }
}

void laser_cool(){
    using namespace consts;
    
    Simulation single_ion;
    single_ion.sim.dt = consts::tau / 20;
    single_ion.init_state(0.2);
    double tottime = 2 * 0.005;
    single_ion.sim.time_end = tottime;
    single_ion.sim.time_engine_start = tottime / 3;
    single_ion.sim.print_every = 5000;

    single_ion.potential = PotentialTypes::Tapered;
    single_ion.calibrateTrapFrequencies();


    single_ion.init_state(0.5);

    // run with lasers
    single_ion.stats = statistics();
    single_ion.phys.saturation = 0.2;
    single_ion.phys.detuning = 13 * MHz;
    single_ion.run();
    single_ion.read_state(std::cerr);

}

void test_taper(){
    using namespace consts;
    
    Simulation sim;
    sim.sim.dt = consts::tau / 100;
    sim.phys.saturation = 0;
    sim.sim.time_end = 0.0001;
    sim.sim.print_every = 100;
    sim.potential = PotentialTypes::Realistic;

    vec freqs;
    for (int axis = 0; axis < 3; axis++) {
        freqs[axis] = sim.trap_freq(axis, 1e-5);
    }
    std::cerr << "Frequencies: " << freqs << std::endl;

    // kick in one direction    
    sim.init_state(0.0);
    sim.x[X] = 100e-6;
    sim.x[Y] = 100e-6;
    sim.run();
    
    sim.read_state(std::cerr);
}

void scan_doppler_temperature(){
    using namespace consts;

    Simulation single_ion;
    single_ion.sim.dt = consts::tau / 20;
    single_ion.init_state(0.2);
    double tottime = 5 * 0.005;
    single_ion.sim.time_end = tottime;
    single_ion.sim.time_engine_start = tottime / 3;
    single_ion.sim.print_every = 50000;

    single_ion.potential = PotentialTypes::Tapered;
    single_ion.calibrateTrapFrequencies();

    for (int detuning = -20; detuning <= -1; detuning++)  {
        single_ion.init_state(0.5);

        // run with lasers
        single_ion.stats = statistics();
        single_ion.phys.saturation = 0.2;
        single_ion.phys.detuning = detuning * MHz;
        single_ion.run();
        single_ion.read_state(std::cerr);
    }
}

int main(int , char** ) {
    using namespace consts;

//    scan_frequency();
//    test_taper();
//    scan_doppler_temperature();

     laser_cool();

    return 0;
}

