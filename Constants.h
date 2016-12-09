/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Constants.h
 * Author: dawid
 *
 * Created on February 3, 2016, 11:12 AM
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

// Constants in SI units
namespace consts {
    typedef const double type;
    // Electron charge [C]
    type electronCharge = 1.602176565e-19;
    //Mass of Calcium [kg]
    type MCa = (40.078 * 1.660538921e-27);
    // m/q ratio
    type qDivM = electronCharge / MCa;
    // Pi
    type pi = 3.141592653589793;
    // Dielectric constant of vacuum
    type epsilon0 = 8.854187818E-12; /* epsilon0 +- .000000071E-12 F/m */
    
    // Force on one ion
    type forceconstant = (electronCharge*electronCharge / (4. * pi*epsilon0)) / MCa;
    
    // speed of light, [m/s]
    type C = 299792458.;
    // millimiter [m]
    type mm = 1. / 1000.;
    // Megahertz SI units [s]
    type MHz = 1000000.;    
    // nano SI prefix
    type nano = 1.0e-9;
    
    // Planck constant, in [J s]
    type hbar = 1.054571726e-34;
    
    // Planck / C, 2 * pi * hbar / C
    type h_div_c = 2 * pi * hbar / C;
    
    // Boltzmann constant, in [J/K]
    type k_B = 1.3806488e-23;
    
    // Resonance lifetime, in [s]
    type tau = 6.904 * nano;
    
    // Resonance linewidth, in [Hz]
    type gamma = 1 / (tau);
    
    // Resonance laser wavelength
    type wavelength = 396.95925 * nano;

    // Resonance laser frequency, in [Hz]
    type resonance_f = C / wavelength;

    type m_over_kb = 0.00480933;
}

#endif /* CONSTANTS_H */

