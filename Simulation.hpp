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

#include "Eigen/Dense"
#include "SimulationParameters.h"
#include "Utilities.h"
#include <random>
#include <ctime>
#include <array>
#include <fstream>
#include <iostream>

typedef Eigen::Vector3d vec;

enum class PotentialTypes {
    Harmonic, 
    Tapered,
    Realistic
};

struct statistics {
    size_t points = 0;
    size_t allocated_size = 0;
    size_t total_decays = 0;
    Eigen::MatrixXd table;
	Eigen::MatrixXd table2;
	size_t runs = 0;

    // averages
    Eigen::RowVector3d avg_x, avg_x2, avg_v, avg_v2;

    statistics(){
        points = 0;
        allocated_size = 0;
        avg_v2.setZero();
        avg_x = avg_x2 = avg_v = avg_v2;
    }

    void do_stats(vec &omegas, std::ostream& out);
};

class Simulation {    
private:
    
    std::random_device rd;
    std::mt19937 rng{rd()};
    std::normal_distribution<> normal{0,1};
    std::uniform_real_distribution<> unif{0,1};
    
    std::string fileName;
    std::fstream outFile;

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
    void calibrateTrapFrequencies(bool verbose = true);
    // Perform statistics
    void do_statistics();
    void initializeMatrices();
	void collect_statistics(const Simulation & traj);
	void ensemble_statistics();

    // define the type of potential, switching below the accelerations below:
    PotentialTypes potential;
    // Calculate forces and energies
    vec acceleration(double);
    vec acceleration_taper(double time);
    vec acceleration_microtaper(double time);
    vec acceleration_harmonic(double time);

    // Calculate energy
    double energy();
    
    // Random normal variable
    double randn();
    // Random double in [0,1]
    double rand();
    // Read state
    void read_state(std::ostream & out = std::cout);
    // Print state
    void print_history();
    // Run simulation from t = 0 to phys.end_time
    void run();
    void run(double time);

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
    // current timestep
    size_t N;
    // data written for later averaging
    size_t printed;
    // decays
    size_t decays;
    
    // acceleration at previous step
    vec a_t;
    // acceleration two steps behind
    vec a_tm;
    
    // scattering probabilities per laser beam
    std::vector<double> probs;
    // set up all laser parameters before the simulation
    void setupLaserBeam();

    // Noise giving function
    std::function<double(double)> noiseFun = [](double ) -> double{return 1.0;};
};

// declare global parameters
extern Physical physical;
// global simulation parameters
extern Parameters  simpar;

#endif /* SIMULATION_HPP */

