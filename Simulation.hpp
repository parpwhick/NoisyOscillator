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
    Physical phys;
    Parameters sim;
    
    
    std::random_device rd;
    std::mt19937 rng{rd()};
    std::normal_distribution<> normal{0,1};
    std::uniform_real_distribution<> unif{0,1};
    
public:    
    Simulation();
    ~Simulation();
    
    void init_state(double T);
    void trap_freq(double kick, int axis);
    vec get_omega();
    double speed();
    double energy();
    double randn();
    double rand();
    
    // Statistics of the simulation
    struct {
        int N;
        int decays;
    } stats;
    
    vec x;
    vec v;
    vec omega;
    
    static double energy(vec pos, vec vel, vec w, Physical par);
    
    

};



#endif /* SIMULATION_HPP */

