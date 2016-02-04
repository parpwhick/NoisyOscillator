/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Simulation.hpp
 * Author: dawid
 *
 * Created on February 3, 2016, 11:06 AM
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "SimulationParameters.h"
#include <random>
#include <ctime>
#include <array>


typedef std::array<double, 3> vec;

class Simulation {    
private:   
    
    std::random_device rd;
    std::mt19937 rng{rd()};
    std::normal_distribution<> normal{0,1};
    std::uniform_real_distribution<> unif{0,1};
    
public:    
    Simulation();
    ~Simulation();
    
    Physical  phys;
    Parameters sim;
    
    double init_state(double T);
    void trap_freq(double kick, int axis);
    const vec & update_omega(double time, vec pos);
    void update_state(double time);
    double speed();
    double update_energy();
    double randn();
    double rand();
    void read_state();

    
    // Statistics of the simulation
    struct {
        int N;
        int decays;
    } stats;
    
    vec x;
    vec v;
    vec omega_squared;
    vec a;
    double t;
    vec energies;
    
    
};



#endif /* SIMULATION_HPP */

