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


#include "Eigen/Dense"
#include "Simulation.hpp"
#include "Constants.h"
#include <iostream>
#include <ostream>
#include <string>

#include "Utilities.h"

Physical physical;
Parameters  simpar;

Simulation::Simulation() : phys(physical), sim(simpar), potential(PotentialTypes::Tapered) {
//    std::cerr <<  "Default constructor called\n";
//    std::cerr << "dt: " << dt << " s, gamma: " << consts::gamma << " Hz\n";
    init_state(phys.T_init);

    decays = printed = N = 0;
    t = 0;

    omega << phys.omega_rad0, phys.omega_rad0, phys.omega_ax0;
}

Simulation::~Simulation() {
}

double Simulation::trap_freq(int axis, double kick) {
    using namespace consts;

    v.setZero();
    x.setZero();
    t = 0.0;
    dt = sim.dt;


    //double sigma_x = std::sqrt(k_B * 1.0 / phys.M / square(omega[axis]));
    // kick = randn() * sigma_x
    x(axis) = kick;
    a = acceleration(t);
    a_t = a_tm = a;

    if(potential != PotentialTypes::Realistic)
        return omega[axis] / 2 / pi;

    auto old_saturation = phys.saturation;
    phys.saturation = 0;
    
    int zero_crossings = 0;
    bool inside_threshold = false;
    bool active = true;
    double outer_threshold = 0.05 * kick;
    double inner_threshold = 0.01 * kick;
    double in_time = 0;
    double last_time = t;
    double first_time = sim.time_end;

    // filter with exponential 1st order causal filter
    int window_len = 10 * std::rint<int>(1 / (dt * phys.RF_omega));
    if (window_len < 2){
        std::cerr << "window too short!!!!" << std::endl;
    }
    double alpha = 1 - 1/std::exp(1.0/window_len);
    double t_delta = dt * window_len / 2.0;
    double x_avg = x[axis];
    
    while (t < sim.time_end){
//        std::cout << t << " " << x[axis] << " " << t-t_delta << " " << x_avg << std::endl;

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

        if (zero_crossings >= 42)
            break;

        step();
        if(std::abs(x[axis])> 2 * kick){
            throw std::runtime_error("The oscillator has run-away behavior, reduce time steps.");
        }
        x_avg += alpha * (x[axis] - x_avg);
    }

    double freq = (zero_crossings-1) / (last_time-first_time);
//    std::cerr << "First time " << first_time << " and last " << last_time << std::endl;

//    std::cerr << "Axis " << axis << " freq: " << freq << ", with crossings " << zero_crossings << std::endl;
    phys.saturation = old_saturation;
    return freq;
}

void Simulation::calibrateTrapFrequencies(bool verbose) {
	using namespace consts;
	vec freqs;
    for (int axis = 0; axis < 3; axis++) {
        freqs[axis] = trap_freq(axis, 10.0e-6);
    }
    if (verbose)
        std::cerr << "Frequencies: " << freqs.transpose() << std::endl;
    phys.omega_rad0 = 2 * pi * freqs[X];
    phys.omega_ax0 = 2 * pi * freqs[Z];
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
    x.setZero();
    v.setZero();
    a = acceleration(t);
	//  phys.RF_phi = 2 * pi * rand();

    phys.T_init = T;

    // https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_velocity_vector
    double sigma = std::sqrt(k_B * T / phys.M);
    for (int i = 0; i < 3; i++) {
        v[i] = sigma * randn();
    }


    omega = {phys.omega_rad0, phys.omega_rad0 * phys.omega_ratio, phys.omega_ax0};
    for (int i = 0; i < 3; i++) {
        double sigma_i = sqrt(k_B * T / phys.M / square(omega[i]));
        x[i] = sigma_i * randn();
    }

    a = acceleration(t);
    a_t = a;
    a_tm = a;
    return energy();
}

double Simulation::init_kick(double kick) {
    using namespace std;
    using namespace consts;
    t = 0;
    x.setConstant(kick);
    v.setZero();
    a = v;
	//    phys.RF_phi = 2 * pi * rand();

    phys.T_init = 0;

    a_t = acceleration(t);
    a_tm = a;
    return energy();
}

double Simulation::energy() {
    // radial potential energy / (M/2)
    auto e_rad_pot = (square(omega[X]) * square(x[X]) + square(omega[Y]) * square(x[Y]));

    // total radial energy / (M/2)
    auto e_rad = e_rad_pot + square(v[X]) + square(v[Y]);
    auto e_ax = square(v[Z]) + square(omega[Z]) * square(x[Z]-phys.z0);
    auto e_tot = e_rad + e_ax;
    energies << 0.5 * e_rad, 0.5 * e_ax, 0.5 * e_tot;
    return 0.5 * e_tot;
}

vec Simulation::acceleration(double tt) {
	switch(potential) {
		case PotentialTypes::Tapered:
			return acceleration_taper(tt);
		case PotentialTypes::Harmonic:
			return acceleration_harmonic(tt);
		default:
			return acceleration_microtaper(tt);
	}
}

// Pseudopotential with a taper
vec Simulation::acceleration_taper(double){
    auto omega_x = phys.omega_rad0 / square(1 + (x[Z]-phys.z0)/phys.zk);
    auto omega_z = phys.omega_ax0;

    // NO MICROMOTION, PONDEROMOTIVE ONLY
    omega << omega_x, omega_x * phys.omega_ratio, omega_z;

    auto z_derivative = -(square(omega[X] * x[X]) + square(omega[Y] * x[Y])) / ((x[Z]-phys.z0) + phys.zk);
    vec acc {   -square(omega[X]) * x[X],
				-square(omega[Y]) * x[Y],
                -square(omega[Z]) * (x[Z]-phys.z0) - z_derivative };
	return acc;
}


vec Simulation::acceleration_harmonic(double){
    // NO MICROMOTION, PONDEROMOTIVE ONLY
    omega = {phys.omega_rad0, phys.omega_rad0 * phys.omega_ratio, phys.omega_ax0};

    vec acc { -square(omega[X]) * x[X],
				-square(omega[Y]) * x[Y],
				-square(omega[Z]) * (x[Z] - phys.z0) };
	return acc;
}

vec Simulation::acceleration_microtaper(double tt){
    double phase = phys.RF_omega * tt + phys.RF_phi;
    double rad_pot_X = 2 * std::cos(phase) * phys.RF_amplitude *
            square(phys.omega_rad0) / square(1 + (x[Z]-phys.z0)/phys.zk);
    double rad_pot_Y = rad_pot_X * phys.omega_ratio;

    auto z_derivative = -(square(x[X]) * rad_pot_X - square(x[Y]) * rad_pot_Y) / ((x[Z]-phys.z0) + phys.zk);

    vec acc   { -rad_pot_X * x[X], 		// X
       			+rad_pot_Y * x[Y],       // Y
       			-square(phys.omega_ax0) * (x[Z]-phys.z0) - z_derivative     // Z
    };

    return acc;
}


void statistics::do_stats(vec& omegas, std::ostream & out) {
using namespace std;
    double millit = 1000 * consts::m_over_kb;
    out << "Avg pos:  " << avg_x / points << endl;
    out << "Avg vel:  " << avg_v / points << endl;
    vec varx = avg_x2/points - (avg_x/points).cwiseAbs2();
    vec varv = avg_v2/points - (avg_v/points).cwiseAbs2();
    out << "Var x:    " << varx.transpose() << endl;
    out << "Var v:    " << varv.transpose() << endl;
    out << "Temps v:  " << varv.transpose() * millit << endl;
    out << "Temps x:  " << (varx.array() * omegas.array().abs2() * millit).transpose() << endl;
}

// Print information about the state in a verbose form
void Simulation::read_state(std::ostream & out){
    using namespace std;
    acceleration(t);
    energy();
    out << "Time: " << t << ", mass: " << phys.M << endl;
    out << "Freqs:    " << omega.transpose()/(2*consts::pi) << endl;
    out << "Lasers:   ";
    for(auto &l : phys.lasers){
        out << "[" << l.transpose() << "], ";
    }
    out << std::endl;
    out << "Detuning: " << phys.detuning / consts::MHz << endl;
    out << "Steps:    " << N << ", decays: " << stats.total_decays << ", printouts: " << printed << endl;
    stats.do_stats(omega, out);
}

// Print a line brief with the information
void Simulation::do_statistics(){
    if(N % sim.print_every == 0){
        Eigen::Matrix<double, 8, 1> data;
        data << t , x , v , (double) decays;
        stats.table.col(printed) = data;
        printed++;
        // reset the decay counter, per printed value
        stats.total_decays += decays;
        decays = 0;
    }
    if(t < sim.time_engine_start)
        return;
    stats.avg_x += x;
    stats.avg_x2 += x.cwiseAbs2();
    stats.avg_v += v;
    stats.avg_v2 += v.cwiseAbs2();
    stats.points++;
}

void Simulation::print_history() {
    fileName = autoFileName();
    outFile.open(fileName, std::fstream::out);
    if(!outFile.good())
        throw std::runtime_error("Could not open file " + fileName + " for output. Aborting");

    for(int i = 0; i < printed; i++){
        outFile << stats.table.col(i).transpose() << std::endl;
    }
}

//  Step forward of dt by using the Velocity Verlet algorithm
void Simulation::step(){
    double dthalf = 0.5 * dt;

    //update position
    x += (v + dthalf * a) * dt;

    //update acceleration on the new position
    vec a_old = a; //it is already x(t+dt)
    a = acceleration(t+dt);

    //update speed
    v += dthalf * (a + a_old);

    //apply stochastic forces
    if (phys.saturation > 0)
        laserXYZ();

    //update time
    t += dt;
    
    //update counter
    N++;
}

// Step forward of dt by using Beeman algorithm
//void Simulation::step(){
//    static double dt = dt;
//
//    //update position
//    for (int i=0; i<3; i++){
//        x[i] += (v[i] + dt / 6.0 * (4* a[i] - a_tm[i])) * dt;
//    }
//
//    // update acceleration on the new position
//    vec a_old = a; //it is already x(t+dt)
//    a = acceleration(t+dt);
//
//    //update speed
//    for (int i=0; i<3; i++){
//        v[i] += dt / 6.0 * (2 * a[i] + 5 * a_t[i] - a_tm[i]);
//    }
//
//    //apply stochastic forces
//    if (phys.saturation > 0)
//        laserXYZ();
//
//    //update time
//    t += dt;
//    
//    //update counter
//    N++;
//}

void Simulation::laserXYZ() {//Formeln aus Apl. Phys. B 45, 175
    // http://info.phys.unm.edu/~ideutsch/Classes/Phys500S09/Downloads/handpubl.pdf
    using namespace consts;

    // Wavevector of a photon at the Ca+ resonance
    static const double ksp = 2 * pi  / wavelength;
    // total rate per laser, per dt
    static const double rate_per_dt = consts::gamma / 2 * dt / phys.lasers.size();

    // Here, the actual frequency of the laser could be used,
    // but the error is of the order of (delta/actual laser freq) ~ 10^-9
    double klaser = ksp + 2 * pi * phys.detuning / C;
    // freq detuning = vparallel / C
    // delta eff = 2 * pi * (delta - vparallel/C)

    double s = phys.lasers.size() * phys.saturation;
    double tot_prob = 0.0;

    for(unsigned int beam = 0; beam < phys.lasers.size(); beam++){
        double delta_eff = 2 * pi * phys.detuning - klaser * phys.lasers[beam].dot(v);
        double delta_norm = delta_eff / consts::gamma;
        // the probability per beam is ~1/N_beams * total_scattering_probability
        probs[beam] = rate_per_dt *  s / (1 + s + 4 * square(delta_norm) );
        tot_prob += probs[beam];
    }

//    if (dt * scatteringrate > 1)
//        std::cout << "timestep too big for laserinteraction" << std::endl;

    // tot_prob is the probability with the correct saturation for multiple laser beams
    // which is then split between the lasers according to their relative absorbtion likelyhood
    if (rand() < tot_prob) {//photonabsorbed
        // count scattering events here
        decays++;

        // determine laser beam
        int beam = 0;
        double threshold = rand() * tot_prob;
        double psum = probs[beam];
        while(threshold > psum){
            beam++;
            psum += probs[beam];
        }

        auto rand_dir = vec{ randn(), randn(), randn() }.normalized();
        v +=  hbar / MCa * ksp * (phys.lasers[beam] + rand_dir);
    }
}


void Simulation::initializeMatrices(){
    int est_time_steps = (int) std::ceil((sim.time_end - t)/dt);
    int est_prints = est_time_steps / sim.print_every + 1;
    stats.allocated_size += est_prints;

    stats.table.conservativeResize(8, stats.allocated_size);

    probs.resize(phys.lasers.size());
}

void Simulation::run(double time){
    dt = sim.dt;
    if (v.norm() < 1e-5 && x.norm() < 1e-8){
        init_state(phys.T_init);
    }
    for(size_t i = 0; i < phys.lasers.size(); i++)
        phys.lasers[i].normalize();

    initializeMatrices();

    if(sim.time_end - t < time)
        sim.time_end = t + time;
    while (t < sim.time_end){
        do_statistics();
        step();
    }
}

void Simulation::run(){
    run(sim.time_end);
}



void ensemble_statistics(std::vector<Simulation> &traj){
    if (traj.empty())
        return;
    Eigen::MatrixXd avg = traj[0].stats.table;
    Eigen::MatrixXd avg2 = traj[0].stats.table.cwiseAbs2();
    int runs = traj.size();
    for(int i = 1; i < runs; i++){
        avg  += traj[i].stats.table;
        avg2 += traj[i].stats.table.cwiseAbs2();

        traj[0].stats.avg_v += traj[i].stats.avg_v;
        traj[0].stats.avg_v2 += traj[i].stats.avg_v2;
        traj[0].stats.avg_x += traj[i].stats.avg_x;
        traj[0].stats.avg_x2 += traj[i].stats.avg_x2;
        traj[0].stats.points += traj[i].stats.points;
    }
    avg /= runs;
    avg2 /= runs;
    traj[0].read_state();

    print_table("avg", avg);
    print_table("avg2", avg2);
}
