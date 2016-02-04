/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   oscillator.cpp
 * Author: Dawid Crivelli
 *
 * Created on February 3, 2016, 11:01 AM
 */

#include <cstdlib>
#include <cstdio>
#include "Simulation.hpp"
#include <vector>
#include <iostream>
using namespace std;

void run_average(Simulation &sim, int N, double T){
    double avg_en = 0.0; 
    
    for(int i = 0; i < N; i++){
        avg_en += sim.init_state(T);
    }
    avg_en /= N;
    
    cout << "Average energy after " << N << " realizations:" << avg_en << endl; 
}

int main(int argc, char** argv) {
    Simulation single_ion;
    
    single_ion.init_state(1);
    single_ion.read_state();
    
    run_average(single_ion, 100000, 1);
    return 0;    
}

