/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Simulation.cpp
 * Author: dawid
 * 
 * Created on February 3, 2016, 11:06 AM
 */

#include "Simulation.hpp"
#include "Constants.h"
#include <iostream>
#include <ostream>
#include <string>
#include "Utilities.h"


Simulation::Simulation() : phys(), sim(phys) {
//    std::cerr <<  "Default constructor called\n";
//    std::cerr << "dt: " << sim.dt << " s, gamma: " << consts::gamma << " Hz\n";
    init_state(phys.T_init);
}

Simulation::~Simulation() {
}

double window_avg(std::vector<double> values, std::vector<double> w, int from){
    double avg = 0.0; 
    int N1 = values.size();
    int N2 = w.size();
    
    //safeguard
    if (N1 != N2)
        return 0.0; 
    
    for(int i=0; i < N2; i++){        
        avg += w[i] * values[(i + from) % N1];
    }
    
    return avg;
}

double Simulation::trap_freq(int axis) {
    using namespace consts;
    
    v = {0.0, 0.0, 0.0};
    x = {0.0, 0.0, 0.0};
    t = 0.0;
    update_omega(x);
    
    double sigma_x = sqrt(k_B * 1.0 / phys.M / square(omega[axis]));    
    double kick = sigma_x;
    x[axis] = kick;    
    update_state(t);
    
    int zero_crossings = 0; 
    bool inside_threshold = false;
    bool active = true;
    double outer_threshold = 0.05 * kick;
    double inner_threshold = 0.01 * kick;
    double in_time = 0;
    double last_time = t;
    double first_time =sim.time_end;
    
    // filter with window
    int window_len = 10 * std::rint<int>(1 / (sim.dt * phys.RF_omega)); 
    if (window_len < 2){
        std::cerr << "window too short!!!!" << std::endl;
    }
    auto w = lanczos<double>(window_len);
    std::vector<double> buf(window_len);
    
    int steps = 0;
    for (int i = 0; i < window_len; i++){
        buf[i] = x[axis];
        step();
        steps += 1;
    }
    
    int from = 0;
    double t_delta = sim.dt * window_len / 2.0;
    while (t < sim.time_end){
        steps += 1;
        //average with a moving window
        double x_avg = window_avg(buf, w, from);
        from = (from + 1) % window_len; 
//        std::cout << t << " " << x[axis] << " " << t - t_delta << " " << x_avg << std::endl; 
        
        //detection of zero crossing on the smoothened signal
        if (x_avg > 0 && x_avg < inner_threshold && !inside_threshold && active) {
            in_time = t;
            inside_threshold = true;
//            std::cerr << "entered at " << in_time - t_delta << std::endl;
        } else if(x_avg > -inner_threshold && inside_threshold){
            last_time = (t + in_time) / 2 - t_delta;
        } 
        
        if (x_avg < -outer_threshold){
            if (inside_threshold){
                zero_crossings += 1;
                first_time = std::min(last_time, first_time);
//                std::cerr << in_time << " " << last_time << ";" << std::endl;
            }
            inside_threshold = false;
            active = false;
        }
        
        if (x_avg > outer_threshold)
            active = true;
        
        if (zero_crossings >= 41)
            break;
        
        step();
        buf[from] = x[axis];
    }
    
    double freq = (zero_crossings-1) / (last_time-first_time);
    std::cerr << "First time " << first_time << " and last " << last_time << std::endl;
    
    std::cerr << "Axis " << axis << " freq: " << freq << ", with crossings " << zero_crossings << " and steps: " << steps << std::endl;    
    return freq;
}

inline double Simulation::randn() {
    return normal(rng);
}

inline double Simulation::rand() {
    return unif(rng);
}

double Simulation::init_state(double T) {
    using namespace std;
    using namespace consts;
    t = sim.time_start;
    x = {0.0, 0.0, 0.0};
    phys.RF_phi = 2 * pi * rand();
    update_omega(x);
    
    phys.T_init = T;
    
    // https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
    double sigma = std::sqrt(k_B * T / phys.M);
    for (auto & vv : v) {
        vv = sigma * randn();
    }

    for (int i = 0; i < 3; i++) {
        double sigma_i = sqrt(k_B * T / phys.M / square(omega[i]));
        x[i] = sigma_i * randn();
    }
    
    return update_energy();
}

inline const vec & Simulation::update_omega(const vec & pos) {    
    auto omega_x = phys.omega_rad0 / square(1 + pos[Z] * phys.zk);
    auto omega_z = phys.omega_ax0;
    
    // NO MICROMOTION, PONDEROMOTIVE ONLY
    omega = {omega_x, omega_x, omega_z};
    
    return omega;
}

void Simulation::update_state(double tt){
    // radial potential energy / (M/2)
    auto e_rad_pot = (square(omega[X]) * square(x[X]) + square(omega[Y]) * square(x[Y]));
    // total radial energy / (M/2)
    auto e_rad = e_rad_pot + square(v[X]) + square(v[Y]);
    auto e_ax = square(v[Z]) + square(omega[Z]) * square(x[Z]);
    
    double phase = 2 * consts::pi * phys.RF_omega * tt + phys.RF_phi;
    double rfvoltage = cos(phase) * phys.RF_amplitude;
    
    auto z_derivative = -(4/2) * e_rad_pot / (x[Z] + 1 / phys.zk);
    a = { -rfvoltage * square(omega[X]) * x[X], // X
          +rfvoltage * square(omega[Y]) * x[Y], // Y
          -square(omega[Z]) * x[Z] + z_derivative // Z
    };
    
    energies = {0.5 * phys.M * e_rad, 0.5 * phys.M * e_ax, 0.5 * phys.M * (e_rad+e_ax)};
}

double Simulation::update_energy() {
    update_omega(x);
    update_state(t);
    return energies[2];
}

inline double Simulation::speed() {
    return module(v);
}


// Print information about the state in a verbose form
void Simulation::read_state(std::ostream & out){
    using namespace std;
    out << "Time: " << t << ", mass: " << phys.M << endl;
    out << "Position: " << x << endl;
    out << "Velocity: " << v << endl;
    out << "Omegas:   " << omega << endl;
    out << "Energies: " << energies << ", kb T: " << consts::k_B * phys.T_init << endl;
}

// Print a line brief with the information
void Simulation::print_state(std::ostream & out){
//    auto list = ;
    out << t << " ";
    for(auto & k  : {x, v, omega,  energies}){
        for (auto & p : k){
            out << p << " ";
        }
    }
    out << std::endl; 
}

//  Step forward of dt by using the Velocity Verlet algorithm
void Simulation::step(){
    static double dt = sim.dt;
    static double dthalf = 0.5 * dt; 
    
    //update position
    for (int i=0; i<3; i++){
        x[i] += (v[i] + dthalf * a[i]) * dt;
    }
    
    //update acceleration on the new position
    update_omega(x); //it is already x(t+dt)
    vec a_old = a;
    update_state(t+dt);
    
    //update speed
    for (int i=0; i<3; i++){
        v[i] += dthalf * (a[i] + a_old[i]);
    }
    
    //apply stochastic forces
    // ....
    
    //update time
    t += dt;
}

void Simulation::run(){
    if (speed() < 1e-5){
        init_state(phys.T_init);
    }
    
    while (t < sim.time_end){
        print_state();
        step();
    }
}