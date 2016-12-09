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
#include "Utilities.h"
#include <random>
#include <ctime>
#include <array>
#include <iostream>

enum class PotentialTypes {
    Harmonic, 
    Tapered,
    Realistic
};

struct statistics {
    int N;
    int decays;
    int printed;
    int points;

    // averages
    vec avg_x = zeros, avg_x2 = zeros, avg_v = zeros, avg_v2 = zeros;

    statistics(){
        points = N = decays = printed = 0;
    }

    void do_stats(vec &omegas, std::ostream& out);
};

class Simulation {    
private:   
    
    std::random_device rd;
    std::mt19937 rng{rd()};
    std::normal_distribution<> normal{0,1};
    std::uniform_real_distribution<> unif{0,1};
    
public:    
    Simulation();
    ~Simulation();
    
    Physical &  phys;
    Parameters & sim;
    
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
    // Perform statistics
    void do_statistics();

    // define the type of potential, switching below the accelerations below:
    PotentialTypes potential;    
    // Calculate forces and energies
    vec acceleration(double time);
    vec acceleration_taper(double time);
    vec acceleration_microtaper(double time);
    vec acceleration_harmonic(double time);
    // Return speed
    double speed();
    // Calculate energy
    double energy();
    
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
    struct statistics stats;
    
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

