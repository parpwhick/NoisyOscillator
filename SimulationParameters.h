#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H

#include "Constants.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "Eigen/Dense"
//using namespace consts;

// Physical parameters

typedef Eigen::Vector3d vec;

struct Physical { 
    // Frequency of the RF electrodes
    double RF_omega;
    // Phase of the RF electrodes
    double RF_phi;
    // Amplitude of the RF electrodes
    double RF_amplitude;

    // Laser detuning
    double detuning;
    // Laser intensity normalized [0 .. 10]
    std::vector<double> saturation;

    // Initial temperature
    double T_init;

    // Mass
    double M;
    // Charge over Mass
    double qDivM;

    // Taper angle, in degrees
    double angle;

    // Tangent of the angle
    double tan_angle;

    //Distance rods to trap center
    double R0 = 1.0 * 0.001;
    
    // x0 / tan(alpha) parameter, used for the determination of omega_rad
    double zk;
    
    //Voltage on the endcaps
    //double Uz = 0.16;

    //Center of trap along z axis
    double z0 = 0.0;
    
    // Frequency radial (average, in Hz)
    double omega_rad0;
    
    // Frequency axial (average, in Hz)
    double omega_ax0;
    
    // Detuning ratio for the second radial frequency, 1+eps
    double omega_ratio;

    // Vector of laser directions, automatically normalized
    std::vector<vec> lasers;

    // Noise intensity on the electrodes
    double noise_amp;
    
    void set_dependent_parameters(){        
        using namespace consts;
        double rad_angle = angle / 180. * pi;
        tan_angle = std::tan(rad_angle);
        qDivM = consts::electronCharge / M;
        zk = R0 / tan_angle;  
        //std::cerr << "Setting zk to: " << zk << std::endl;
    }
    
    Physical() {        
        using namespace consts;
        RF_omega = 2 * pi * 22 * MHz;
        RF_amplitude = 38;
        angle = 10.0;
        
        // Pulse (Frequency) radial (average, in Hz)
        omega_rad0 = 2 * pi * 0.440 * MHz;
    
        // Pulse (Frequency) axial (average, in Hz)
        omega_ax0  = 2 * pi * 0.080 * MHz;

        omega_ratio = 1.05;

        // detuning of the laser in Hz
        detuning = -consts::gamma/(2 * pi)/2;

        // Mass of Calcium
        M = MCa;

        // saturation
        saturation = { 0.5 };

        // noise = 0
        noise_amp = 0.0;

        // laser addressing all three directions
        lasers = { vec(1.0, 1.0, 1.0) };

        // initial temperature
        T_init = 0.005; // 5 mK
        
        set_dependent_parameters();
    }
};


// Simulation parameters

class Parameters {
public:

    // Time step for the simulation
    double dt;
    double time_end;
    int steps;
    int print_every = 500;

    double time_engine_start;
    double laser_initial_off;
    double laser_initial_on;

    Parameters() {
        using namespace consts;

        // default 25 steps per micromotion oscillation
        dt = 5e-10;
        set_dependent_parameters();
    }

    void set_dependent_parameters();

};

#endif /* SIMULATIONPARAMETERS_H */

