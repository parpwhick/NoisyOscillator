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

    omega = {phys.omega_rad0, phys.omega_rad0, phys.omega_ax0};
    stats = {};
}

Simulation::~Simulation() {
}

double Simulation::trap_freq(int axis, double kick) {
    using namespace consts;

    v = {0.0, 0.0, 0.0};
    x = {0.0, 0.0, 0.0};
    t = 0.0;

    //double sigma_x = std::sqrt(k_B * 1.0 / phys.M / square(omega[axis]));
    // kick = randn() * sigma_x
    x[axis] = kick;
    a = acceleration(t);
    a_t = a_tm = a;
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
    int window_len = 10 * std::rint<int>(1 / (sim.dt * phys.RF_omega));
    if (window_len < 2){
        std::cerr << "window too short!!!!" << std::endl;
    }
    double alpha = 1 - 1/std::exp(1.0/window_len);
    double t_delta = sim.dt * window_len / 2.0;
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

        if (zero_crossings >= 41)
            break;

        step();
        x_avg += alpha * (x[axis] - x_avg);
    }

    double freq = (zero_crossings-1) / (last_time-first_time);
//    std::cerr << "First time " << first_time << " and last " << last_time << std::endl;

//    std::cerr << "Axis " << axis << " freq: " << freq << ", with crossings " << zero_crossings << std::endl;
    phys.saturation = old_saturation;
    return freq;
}

void Simulation::calibrateTrapFrequencies() {
	using namespace consts;
	vec freqs;
    for (int axis = 0; axis < 3; axis++) {
        freqs[axis] = trap_freq(axis, 1e-3);
    }
    std::cerr << "Frequencies: " << freqs << std::endl;
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
    t = sim.time_start;
    x = {0.0, 0.0, 0.0};
	//  phys.RF_phi = 2 * pi * rand();

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

    a = acceleration(t);
    a_t = a;
    a_tm = a;
    return update_energy();
}

double Simulation::init_kick(double kick) {
    using namespace std;
    using namespace consts;
    t = sim.time_start;
    x = {kick, kick, kick};
    v = {0, 0, 0};
    a = v;
	//    phys.RF_phi = 2 * pi * rand();

    phys.T_init = 0;

    a_t = acceleration(t);
    a_tm = a;
    return update_energy();
}

vec Simulation::acceleration(double tt){
	auto omega_x = phys.omega_rad0; // * (phys.real_potential) ? 1. / square(1 + x[Z] * phys.zk) : 1.0;
    auto omega_z = phys.omega_ax0;

    // NO MICROMOTION, PONDEROMOTIVE ONLY
    omega = {omega_x, omega_x, omega_z};

    // radial potential energy / (M/2)
    auto e_rad_pot = (square(omega[X]) * square(x[X]) + square(omega[Y]) * square(x[Y]));

    // total radial energy / (M/2)
    auto e_rad = e_rad_pot + square(v[X]) + square(v[Y]);
    auto e_ax = square(v[Z]) + square(omega[Z]) * square(x[Z]);
    energies = {0.5 * e_rad, 0.5 * e_ax, 0.5 * (e_rad+e_ax)};

    if (phys.real_potential){
    	double phase = phys.RF_omega * tt + phys.RF_phi;
    	double rfvoltage = std::cos(phase) * phys.RF_amplitude;

    	auto z_derivative = -(4/2) * e_rad_pot / (x[Z] + 1 / phys.zk);

    	vec acc = { -rfvoltage * square(omega[X]) * x[X], // X
        	  +rfvoltage * square(omega[Y]) * x[Y],       // Y
          	-square(omega[Z]) * x[Z] + z_derivative     // Z
    	};

    	return acc;
	} else {
		vec acc = { -square(omega[X]) * x[X], 
					-square(omega[Y]) * x[Y],
					-square(omega[Z]) * x[Z] };
		return acc;
	}
}


double Simulation::update_energy() {
    acceleration(t);
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
    out << "Steps:    " << stats.N << ", decays: " << stats.decays << ", printouts: " << stats.printed << endl;
}

// Print a line brief with the information
void Simulation::print_state(std::ostream & out){
//    auto list = ;
    out << t << " ";
    for(auto & k  : {x, v, energies}){
        for (auto & p : k){
            out << p << " ";
        }
    }
    out << std::endl;
    stats.printed++;
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
    vec a_old = a; //it is already x(t+dt)
    a = acceleration(t+dt);

    //update speed
    for (int i=0; i<3; i++){
        v[i] += dthalf * (a[i] + a_old[i]);
    }

    //apply stochastic forces
    if (phys.saturation > 0)
        laserXYZ();

    //update time
    t += dt;
    
    //update counter
    stats.N++;
}

// Step forward of dt by using Beeman algorithm
//void Simulation::step(){
//    static double dt = sim.dt;
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
//    stats.N++;
//}

static const vec directionXYZ = {1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3)};
static const vec directionZ = {0, 0, -1};

void Simulation::laserXYZ() {//Formeln aus Apl. Phys. B 45, 175
    using namespace consts;
    // Momentum of a photon at the Ca+ resonance, about 1.669e-27 [J s / m]
    static const double ksp = 2 * pi  / wavelength;
    static const vec direction = directionZ;
    static const double gamma2 = square(consts::gamma);

    double klaser = ksp + 2 * pi * phys.detuning / C;
    vec laser_k = {klaser * direction[X], klaser * direction[Y], klaser * direction[Z]};

    double term = 2 * pi * phys.detuning - scalar(laser_k, v);
    double saturation = phys.saturation * gamma2 / (square(term) + gamma2);
    // http://info.phys.unm.edu/~ideutsch/Classes/Phys500S09/Downloads/handpubl.pdf
    double scatteringrate = consts::gamma * saturation / (1 + 2 * saturation);

    if (sim.dt * scatteringrate > 1) 
        std::cout << "timestep too big for laserinteraction" << std::endl;

    if (rand() <= sim.dt * scatteringrate) {//photonabsorbed
        // step time after the scattering
        // ...

        // count scattering events here
        stats.decays++;
        
        vec rand_dir;
        for (auto &d : rand_dir)
            d = randn();
        normalize(rand_dir);

        for (int r = 0; r < 3; r++) {
            v[r] += laser_k[r] * hbar / MCa  + ksp * hbar / MCa * rand_dir[r];
        }
    }
}


void Simulation::run(){
    if (speed() < 1e-5 && module(x) < 1e-8){
        init_state(phys.T_init);
    }

    while (t < sim.time_end){
        if(stats.N % sim.print_every == 0)
            print_state();
        step();
    }
}