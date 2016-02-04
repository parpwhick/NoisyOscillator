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
#include <cstdio>


const int X = 0;
const int Y = 1;
const int Z = 2;

Simulation::Simulation() : sim(phys) {
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

void Simulation::init_state(double T) {
    double sigma = std::sqrt(0.5 * consts::k_B * T);
    for (auto & vv : v) {
        vv = sigma * randn();
    }

    for (int i = 0; i < 3; i++) {
        double sigma_i = std::sqrt(0.5 * consts::k_B * T * omega[i]);
        x[i] = sigma_i * randn();
    }
}

vec Simulation::get_omega() {
    auto omega_xy = 2 * consts::pi * (x[Z]*(-2241983470.) + 2969444.);
    auto omega_z = 2 * consts::pi * 36003;
    omega[X] = omega[Y] = omega_xy;
    omega[Z] = omega_z;

    return omega;
}

vec force(double t, vec pos, vec vel, vec w, Physical par) {
    double & tan_angle = par.tan_angle;
    vec f = {0.0, 0.0, 0.0};


    double phase = 2 * consts::pi * par.RF_omega * t + par.RF_phi;
    double rfvoltage = par.RF_amplitude * cos(phase);

    double taper_z = par.x0 + pos[Z] * tan_angle;
    double accPot = rfvoltage / (taper_z * taper_z);

    // 
    double accAx = par.Uz / (par.z0 * par.z0);


    f[X] = -2. * par.qDivM * accPot * pos[X];
    f[Y] = 2. * par.qDivM * accPot * pos[Y];
    f[Z] = 2. * par.qDivM * (accPot / taper_z * tan_angle * (pos[X] * pos[X] - pos[Y] * pos[Y])
            - accAx * (pos[Z] - par.zOffset));
    return f;
}

double Simulation::energy(vec pos, vec vel, vec w, Physical par) {
    double erad = 0.0;
    for (int i = 0; i < 2; i++) {
        erad += pos[i] * pos[i] * w[i] * w[i] + vel[i] * vel[i];
    }
    auto ez = vel[Z] * vel[Z] + pos[Z] * pos[Z] * w[Z] * w[Z];

    return 0.5 * par.M * (erad + ez);
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

inline double Simulation::energy() {
    get_omega();
    return energy(x, v, omega, phys);
}