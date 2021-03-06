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

void run_average(Simulation &sim, size_t N, double T) {
    double avg_en = 0.0;

    for (size_t i = 0; i < N; i++) {
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
    double tottime = 0.005;
    single_ion.sim.time_end = tottime;
    single_ion.sim.time_engine_start = tottime / 3;
    single_ion.sim.print_every = 500;

    single_ion.potential = PotentialTypes::Tapered;
    single_ion.calibrateTrapFrequencies();

    single_ion.init_state(0.5);

    // run with lasers
    single_ion.stats = statistics();
    single_ion.phys.saturation = {1};
    single_ion.phys.detuning = -10 * MHz;
    single_ion.run();
    single_ion.read_state(std::cerr);
    single_ion.print_history();

}

void test_taper(){
    using namespace consts;
    
    Simulation sim;
    sim.sim.dt = consts::tau / 100;
    sim.phys.saturation = {0};
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
        single_ion.phys.saturation = {0.2};
        single_ion.phys.detuning = detuning * MHz;
        single_ion.run();
        single_ion.read_state(std::cerr);
    }
}


void averaged_runs(){
    using namespace consts;
    const int runs = 4;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        Simulation traj;
        traj.sim.dt = consts::tau / 20;
        double tottime = 0.008;        
        traj.sim.time_end = tottime;
        traj.sim.time_engine_start = tottime / 3;
        traj.sim.print_every = 500;
        traj.potential = PotentialTypes::Tapered;
        traj.calibrateTrapFrequencies(false);

        traj.init_state(0.5);
        // run with lasers
        traj.stats = statistics();
        traj.phys.saturation = {1};
        traj.phys.detuning = -10 * MHz;
        traj.run();

        printf("%3d -- done\n", i);
    }
}


void initial_temperature(){
    using namespace consts;
    const int runs = 100000;
    Simulation traj;


    Eigen::Vector3d v2;
    v2.setZero();

    Eigen::Vector3d x2;
    x2.setZero();

    double ee = 0.0;

    traj.sim.dt = consts::tau / 20;
    double tottime = 0.002;
    traj.sim.time_end = tottime;
    traj.sim.time_engine_start = tottime / 3;
    traj.sim.print_every = 500;
    traj.potential = PotentialTypes::Tapered;
    traj.calibrateTrapFrequencies(false);


    for(int i = 0; i < runs; i++){
        traj.stats = statistics();
        ee += traj.init_state(0.5);
        v2 += traj.v.cwiseAbs2();
        x2 += (traj.x.array() * traj.omega.array()).square().matrix() * 4.8;
        //        std::cerr << "----- RUN " << i << " -------" << std::endl;
        //        std::cerr << "T: " << traj.energy()*4.8 << std::endl;
        //        std::cerr << std::endl;
    }

    std::cerr << v2.transpose() / runs * 4.8 << std::endl;
    std::cerr << (x2/runs).transpose() << std::endl;
    std::cerr << ee / runs * consts::m_over_kb * 1000.0 / 3 << std::endl;
}




void initial_state(){
    using namespace consts;
    const int runs = 10000;
    Simulation traj;

    traj.sim.dt = consts::tau / 20;
    double tottime = 0.002;
    traj.sim.time_end = tottime;
    traj.sim.time_engine_start = tottime / 3;
    traj.sim.print_every = 500;
    traj.potential = PotentialTypes::Tapered;
    traj.calibrateTrapFrequencies(false);


    for(int i = 0; i < runs; i++){
        traj.stats = statistics();
        traj.init_state(0.5);
        std::cout << i << " " << traj.x.transpose() << " " <<
                     traj.v.transpose() << std::endl;
    }
}



void fluorescence(){
    using namespace consts;
    const int runs = 16;

    double tottime = 0.010;
    simpar.dt = 5e-10;
    simpar.time_end = tottime;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 500;
    physical.omega_rad0 = 2 * pi * 0.9 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.200 * MHz;
    physical.saturation = {0.3, 0.3};
    physical.detuning = -20 * MHz;

    physical.lasers = { vec(1,1,0), vec(0,0,1) };
    
    Simulation results;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        Simulation traj;
        traj.potential = PotentialTypes::Tapered;
        traj.calibrateTrapFrequencies(false);
        traj.init_state(0.010);

        // run with lasers
        traj.stats = statistics();
        traj.run();

        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj.read_state(std::cerr);
        std::cerr << std::endl;
        if (i == 0)
            traj.print_history();
        #pragma omp critical 
        {
            results.collect_statistics(traj);
        }
    }    
    results.ensemble_statistics();
}


void test_noise(){
    using namespace consts;
    const int runs = 60;

    double tottime = 0.0004;
    simpar.dt = 5e-10;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 500;
    physical.omega_rad0 = 2 * pi * 0.4 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.080 * MHz;
    physical.saturation = {1};
    physical.detuning = -20 * MHz;
    physical.noise_amp = 0;
    physical.lasers = { vec(1,1,1) };

    double amp = 0.040 /simpar.dt;

    Simulation results;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        Simulation traj;
        traj.potential = PotentialTypes::Tapered;
        traj.calibrateTrapFrequencies(false);
        traj.init_state(0.010);
        traj.stats = statistics();

        for(int tt = 0; tt < 60; tt++){
            traj.phys.noise_amp =  0;
            traj.run(tottime);

            traj.phys.noise_amp =  amp;
            traj.run(tottime);
        }

        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj.read_state(std::cerr);
        std::cerr << std::endl;
        if (i == 0)
            traj.print_history();
        #pragma omp critical 
        {
            results.collect_statistics(traj);
        }
    }    
    results.ensemble_statistics();
}

void noise_ramp(){
    using namespace consts;
    const int runs = 60;

    double tottime = 0.001;
    simpar.dt = 5e-10;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 500;
    physical.omega_rad0 = 2 * pi * 0.4 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.080 * MHz;
    physical.saturation = {2};
    physical.detuning = -20 * MHz;
    physical.noise_amp = 0;
    physical.lasers = { vec(1,1,1) };


    Simulation results;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        Simulation traj;
        traj.potential = PotentialTypes::Tapered;
        traj.calibrateTrapFrequencies(false);
        traj.init_state(0.0001);
        traj.stats = statistics();

        for(int tt = 0; tt < 30; tt++){
            traj.phys.noise_amp =  0;
            traj.run(tottime);
        }


        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj.read_state(std::cerr);
        std::cerr << std::endl;
        if (i == 0)
            traj.print_history();
        #pragma omp critical 
        {
            results.collect_statistics(traj);
        }
    }    
    results.ensemble_statistics();
}


void laser_heating(){
    using namespace consts;
    const int runs = 60;

    double tottime = 0.001;
    simpar.dt = 5e-10;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 500;
    physical.omega_rad0 = 2 * pi * 0.4 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.080 * MHz;
    physical.saturation = {1};
    physical.detuning = -3 * MHz;
    physical.noise_amp = 0;
    physical.lasers = { vec(1,1,0), vec(-1,-1,0), vec(0,0,1), vec(0,0,-1) };

    Simulation results;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        Simulation traj;
        traj.potential = PotentialTypes::Tapered;
        traj.calibrateTrapFrequencies(false);
        traj.init_state(0.010);
        traj.stats = statistics();

        for(int tt = 0; tt < 20; tt++){
            traj.phys.saturation = {1.0, 1.0, 1.0, 1.0};
            traj.run(tottime);

            traj.phys.saturation = {40.0, 40.0, 1.0, 1.0};
            traj.run(tottime);
        }


        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj.read_state(std::cerr);
        std::cerr << std::endl;
        if (i == 0)
            traj.print_history();
        #pragma omp critical 
        {
            results.collect_statistics(traj);
        }
    }
    results.ensemble_statistics();
}

void heat_engine(){
    using namespace consts;
    const int runs = 100;

    double tottime = 0.010;
    simpar.dt = 5e-10;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 200;
    physical.omega_rad0 = 2 * pi * 0.4 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.080 * MHz;
    physical.saturation = {1, 1, 1, 1};
    physical.detuning = -20 * MHz;
    double dutyCycle = 0.1;
    double threshhold = std::sin(pi*(0.5-dutyCycle));
    physical.lasers = { vec(1,1,0), vec(-1,-1,0), vec(0,0,1), vec(0,0,-1) };

    Simulation results;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        Simulation traj;
        traj.potential = PotentialTypes::Tapered;
        traj.calibrateTrapFrequencies(false);
        traj.init_state(0.010);
        traj.stats = statistics();
        double freq = physical.omega_ax0;
    	traj.phys.noise_amp = 0.12 / simpar.dt;
        traj.noiseFun = [threshhold, freq, tottime](double t)->double{
            if(t < tottime / 3)
                return 0;
            else
                return (sin(t * freq) > threshhold);
        };
        traj.run(tottime);


        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj.read_state(std::cerr);
        std::cerr << std::endl;
        if (i == 0)
            traj.print_history();
        #pragma omp critical 
        {
            results.collect_statistics(traj);
        }
    }    
    results.ensemble_statistics();
}

int main(int , char** ) {
    using namespace consts;

//    scan_frequency();
//    test_taper();
//    scan_doppler_temperature();

//     laser_cool();
//    averaged_runs();
//    initial_temperature();
//    initial_state();
//    fluorescence();
   // test_noise();
     laser_heating();
//     heat_engine();
    return 0;
}

