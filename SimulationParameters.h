#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H

#include "Constants.h"
#include <cmath>
//using namespace consts;

// Physical parameters

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
    double saturation;

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
    double x0 = 1.0 * 0.001;

    //Distance endcaps to trap center
    double z0 = 4.0 * 0.001;
    
    //Voltage on the endcaps
    double Uz = 0.16;

    //Center of trap along z axis
    double zOffset = 0.0;

    
    void set_dependent_parameters(){        
        using namespace consts;
        double rad_angle = angle / 180. * pi;
        tan_angle = std::tan(rad_angle);
        qDivM = consts::electronCharge / M;
    }
    
    Physical() {        
        using namespace consts;
        RF_omega = 22.7 * MHz;

        // detuning of the laser in Hz
        detuning = -consts::gamma;

        // Mass of Calcium
        M = MCa;

        // saturation
        saturation = 5.0;

        // initial temperature
        T_init = 0.080; // 80 mK
        
        set_dependent_parameters();
    }
};


// Simulation parameters

class Parameters {
public:

    // Time step for the simulation
    double dt;
    double time_start;
    double time_end;
    int steps;

    double time_engine_start;
    double laser_initial_off;
    double laser_initial_on;

    Parameters(Physical phys) {
        using namespace consts;

        // default 10 steps per micromotion oscillation
        dt = 0.1 / phys.RF_omega;

        set_dependent_parameters();
    }

    void set_dependent_parameters();

};

#endif /* SIMULATIONPARAMETERS_H */

