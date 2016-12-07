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
#include <iostream>

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
    
    // Init state with temperature T on all axes
    double init_state(double T);
    // Init state in elongation kick in all direction
    double init_kick(double kick);
    // Step forward one default time step
    void step();
    // Perform one laser interaction step
    void laserXYZ();
    // Calculate trap frequency on one axis
    double trap_freq(int axis, double kick = 1e-5);
    // Calculate ponderomotive frequencies by measuring them
    void calibrateTrapFrequencies();

    // Update the position-dependent trap frequency
    //const vec & update_omega(const vec & pos);
    
    // Calculate forces and energies
    vec acceleration(double time);
    // Return speed
    double speed();
    // Calculate energy
    double update_energy();
    
    // Random normal variable
    double randn();
    // Random double in [0,1]
    double rand();
    // Read state
    void read_state(std::ostream & out = std::cout);
    // Print state
    void print_state(std::ostream & out = std::cout);
    // Run simulation from t = 0 to phys.end_time
    void run();

    // Statistics of the simulation
    struct {
        int N;
        int decays;
        int printed;
    } stats;
    
    // current time
    double t;
    // current position
    vec x;
    // current speed
    vec v;
    // current acceleration
    vec a;
    vec energies;
    vec omega;
    
    // acceleration at previous step
    vec a_t;
    // acceleration two steps behind
    vec a_tm;
    
};



#endif /* SIMULATION_HPP */

