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
    std::vector<Simulation> traj(runs);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        traj[i].sim.dt = consts::tau / 20;
        double tottime = 0.008;        
        traj[i].sim.time_end = tottime;
        traj[i].sim.time_engine_start = tottime / 3;
        traj[i].sim.print_every = 500;
        traj[i].potential = PotentialTypes::Tapered;
        traj[i].calibrateTrapFrequencies(false);

        traj[i].init_state(0.5);
        // run with lasers
        traj[i].stats = statistics();
        traj[i].phys.saturation = {1};
        traj[i].phys.detuning = -10 * MHz;
        traj[i].run();

        printf("%3d -- done\n", i);
    }

    traj[0].print_history();

    for(int i = 0; i < runs; i++){
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj[i].read_state(std::cerr);
        std::cerr << std::endl;
    }

    ensemble_statistics(traj);
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
        //        std::cerr << "T: " << traj[i].energy()*4.8 << std::endl;
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
    std::vector<Simulation> traj(runs);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        traj[i].potential = PotentialTypes::Tapered;
        traj[i].calibrateTrapFrequencies(false);
        traj[i].init_state(0.010);

        // run with lasers
        traj[i].stats = statistics();
        traj[i].run();

        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj[i].read_state(std::cerr);
        std::cerr << std::endl;
    }

    traj[0].print_history();
    ensemble_statistics(traj);
}


void test_noise(){
    using namespace consts;
    const int runs = 20;

    double tottime = 0.0004;
    simpar.dt = 5e-10;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 500;
    physical.omega_rad0 = 2 * pi * 0.9 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.200 * MHz;
    physical.saturation = {1};
    physical.detuning = -20 * MHz;
    physical.noise_amp = 0;
    physical.lasers = { vec(1,1,1) };
    std::vector<Simulation> traj(runs);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        traj[i].potential = PotentialTypes::Tapered;
        traj[i].calibrateTrapFrequencies(false);
        traj[i].init_state(0.010);
        traj[i].stats = statistics();

        for(int tt = 0; tt < 10; tt++){
            traj[i].phys.noise_amp =  0;
            traj[i].run(tottime);

            traj[i].phys.noise_amp =  0.010 /simpar.dt;
            traj[i].run(tottime);
        }


        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj[i].read_state(std::cerr);
        std::cerr << std::endl;
    }

    traj[0].print_history();
    ensemble_statistics(traj);
}

void heat_engine(){
    using namespace consts;
    const int runs = 150;

    double tottime = 0.010;
    simpar.dt = 5e-10;
    simpar.time_engine_start = tottime / 3;
    simpar.print_every = 500;
    physical.omega_rad0 = 2 * pi * 0.9 * MHz;
    physical.omega_ratio = 1.222;
    physical.omega_ax0 = 2 * pi * 0.200 * MHz;
    physical.saturation = {0.5};
    physical.detuning = -20 * MHz;
    physical.noise_amp = 0.050 / simpar.dt;
    double dutyCycle = 0.5;
    double threshhold = std::sin(pi*(0.5-dutyCycle));
    physical.lasers = { vec(1,1,1) };
    std::vector<Simulation> traj(runs);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < runs; i++){
        traj[i].potential = PotentialTypes::Tapered;
        traj[i].calibrateTrapFrequencies(false);
        traj[i].init_state(0.010);
        traj[i].stats = statistics();
        double freq = physical.omega_ax0;
        traj[i].noiseFun = [threshhold, freq, tottime](double t)->double{
            if(t < tottime / 3)
                return 0;
            else
                return (sin(t * freq) > threshhold);
        };
        traj[i].run(tottime);


        printf("%3d -- done\n", i);
        std::cerr << "----- RUN " << i << " -------" << std::endl;
        traj[i].read_state(std::cerr);
        std::cerr << std::endl;
    }

    traj[0].print_history();
    ensemble_statistics(traj);
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
//    test_noise();
    heat_engine();
    return 0;
}

