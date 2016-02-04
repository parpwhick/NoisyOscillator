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
#include <cstdio>
#include <iostream>
#include <ostream>
#include <string>


const int X = 0;
const int Y = 1;
const int Z = 2;

template<typename T> inline T square(const T & x){
    return x*x;
}

Simulation::Simulation() : phys(), sim(phys) {
    printf("Default constructor called\n");
    printf("dt: %.3g s, gamma: %.3g Hz\n", sim.dt, consts::gamma);
    init_state(phys.T_init);
}

Simulation::~Simulation() {
    printf("Deconstructed\n");
}

void Simulation::trap_freq(double kick, int axis) {

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
    t = 0;
    x = {0.0, 0.0, 0.0};
    update_omega(t, x);
    
    phys.T_init = T;
    
    // https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
    double sigma = std::sqrt(k_B * T / phys.M);
    for (auto & vv : v) {
        vv = sigma * randn();
    }

    for (int i = 0; i < 3; i++) {
        // remember that the omega_squared are the micromotion-RMS adjusted frequency,
        // so they need to be divided by sqrt(2))
        double sigma_i = sqrt(k_B * T / phys.M / abs(omega_squared[i] / sqrt(2.0)));
        x[i] = sigma_i * randn();
    }
    
    return update_energy();
}

inline const vec & Simulation::update_omega(double time, vec pos) {
    double t = 0.0;
    double phase = 2 * consts::pi * phys.RF_omega * time + phys.RF_phi;
    double rfvoltage = phys.RF_amplitude * cos(phase);
    
    auto omega_x = phys.omega_rad0 / square(1 + pos[Z] * phys.zk);
    auto omega_z = phys.omega_ax0;
    auto omega_x_t = rfvoltage * omega_x;
    omega_squared = {square(omega_x_t), -square(omega_x_t), square(omega_z)};
    
    return omega_squared;
}

void Simulation::update_state(double t){
    // radial potential energy / (M/2)
    auto e_rad_pot = (abs(omega_squared[X]) * square(x[X]) + abs(omega_squared[Y]) * square(x[Y]));
    // total radial energy / (M/2)
    auto e_rad = e_rad_pot + square(v[X]) + square(v[Y]);
    auto e_ax = square(v[Z]) + omega_squared[Z] * square(x[Z]);
    
    auto z_derivative = -(4/2) * e_rad_pot / (x[Z] + 1 / phys.zk);
    a = { -omega_squared[X] * x[X], // X
          -omega_squared[Y] * x[Y], // Y
          -omega_squared[Z] * x[Z] + z_derivative // Z
    };
    
    energies = {0.5 * phys.M * e_rad, 0.5 * phys.M * e_ax, 0.5 * phys.M * (e_rad+e_ax)};
}

double Simulation::update_energy() {
    update_omega(t, x);
    update_state(t);
    return energies[2];
}

inline double module(vec vector) {
    double m;
    for (auto & x : vector) {
        m += x * x;
    }
    return std::sqrt(m);
}

inline double Simulation::speed() {
    return module(v);
}

template <typename T, std::size_t N>
std::ostream& operator<< (std::ostream& out, const std::array<T, N>& v) {
    out << '[';
    for (int i = 0; i < N; i++)
        out << v[i] << ", "; 
    out << "\b\b]";
  return out;
}

void Simulation::read_state(){
    using namespace std;
    cout << "Time: " << t << ", mass: " << phys.M << endl;
    cout << "Position: " << x << endl;
    cout << "Velocity: " << v << endl;
    cout << "Omegas:   ["; 
    for(auto om : omega_squared) 
        cout << sqrt(abs(om)/2) << ", ";
    cout << "\b\b]" << endl;
    cout << "Energies: " << energies << ", kb T: " << consts::k_B * phys.T_init << endl;
    //printf("Position: (%.3g, %.3g, %.3g)\n", x);
}
