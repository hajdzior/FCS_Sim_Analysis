// simFCS.cpp
// Project started: 20.12.2016
// Last modified: 20.12.2016
// Units are micrometers for length, and seconds for time  
// Change Log: 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <numeric>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "FCS.hpp"
#include "config.hpp"

int main ( ) {
    // Read the wall clock time    
    clock_t begin = clock();
    
    // Read the configuration file    
    ConfigFile cfg("config.cfg");

    // Initialize random number generator--------------------------------------
    gsl_rng * rng;
    rng = gsl_rng_alloc ( gsl_rng_mt19937 );
    const std::string set_seed = cfg.getValueOfKey<std::string>("set_seed");
    // Random seed
    long seed = time(0);
    if (set_seed == "\"fixed\"") {
    // fixed seed
        seed = 0;
    } 
    std::cout << "set_seed: " << set_seed << std::endl;
    // Set seed
    gsl_rng_set(rng, seed);
    // ------------------------------------------------------------------------    
     
    // Time duration of the simulation
    const double t_max = cfg.getValueOfKey<double>("t_max"); // seconds
    // Simulation time step
    const double t_step = cfg.getValueOfKey<double>("t_step"); // seconds
    // Diffusion coefficient
    const double D = cfg.getValueOfKey<double>("D"); // um^2/s
    // Set the maximum excitation rate     
    const double k_ex_max = cfg.getValueOfKey<double>("k_ex_max");    
    std::cout << "k_ex_max= " << k_ex_max << std::endl; 
    // Detector efficiency: ??? to be removed ??? 
    const double phi_det = cfg.getValueOfKey<double>("phi_det");   

    const unsigned int n_particles = cfg.getValueOfKey<unsigned int>("n_particles"); // number of particles
    
    const unsigned int n_steps = t_max/t_step;
    std::cout << "n_steps= " << n_steps << std::endl;
    std::vector<unsigned int> counts(n_steps);   
    
    // Set blinking: triplet or powerlaw--------------------------------------- 
    const std::string blinking = cfg.getValueOfKey<std::string>("blinking");
    std::cout << "blinking= " << blinking << std::endl;
    if ( blinking != "\"triplet\"") {
        std::cout << "simFCS.cpp: Powerlaw blinking is not implemented yet. Exiting ..." << std::endl; 
        exit(1); 
    }
    const double phi_f = cfg.getValueOfKey<double>("phi_f"); 
    const double phi_st = cfg.getValueOfKey<double>("phi_st"); 
    const double phi_sb = cfg.getValueOfKey<double>("phi_sb"); 
    const double k_tb = cfg.getValueOfKey<double>("k_tb");
    const double k_ts = cfg.getValueOfKey<double>("k_ts"); 
    //--------------------------------------------------------------------------

    // Set the molecular detection function: Gaussian or semigeometric---------
    const std::string confocal_geometry = cfg.getValueOfKey<std::string>("confocal_geometry");    

    const double w_xy = cfg.getValueOfKey<double>("w_xy");      
    const double w_z = cfg.getValueOfKey<double>("w_z");

    const double NA = cfg.getValueOfKey<double>("NA");                  // numerical aperture
    const double M = cfg.getValueOfKey<double>("M");                    // magnification
    const double n = cfg.getValueOfKey<double>("n");                    // index of refraction
    const double r = cfg.getValueOfKey<double>("r");                    // radius
    const double lambda_ex = cfg.getValueOfKey<double>("lambda_ex");    //excitation wavelength
    const double lambda_em = cfg.getValueOfKey<double>("lambda_em");    //emission wave length
    const double w_0 = cfg.getValueOfKey<double>("w_0");                // ???

    MDF mdf = MDF(confocal_geometry, w_xy, w_z, NA, M, n, r, lambda_ex, lambda_em, w_0, phi_f); 
    //-------------------------------------------------------------------------

    // Set the simulation box--------------------------------------------------
    const double box_size = cfg.getValueOfKey<double>("box_size");;
    Box box = Box(-box_size/2.0, +box_size/2.0, -box_size/2.0, +box_size/2.0, -box_size/2.0, +box_size/2.0);
    //-------------------------------------------------------------------------

    const double AvogadroNum =  6.022140857e+23;  
    const double Molarity = (n_particles/AvogadroNum)/box.volume_L();  
    std::cout << "c= " << Molarity*1e9 << " nM" << std::endl;    

    // Set the simulation system ----------------------------------------------    
    Molecules molecules = Molecules(D, n_particles, box, rng, blinking);
    //-------------------------------------------------------------------------

    std::stringstream dir_name; 
    dir_name << "./data/D=" << D
                            << "-phi_st=" << phi_st
                            << "-phi_sb=" << phi_sb
                            << "-k_ts=" << k_ts/1e3 << "kHz" 
                            << "-k_tb=" << k_tb/1e3 << "kHz"
                            << "-n_particles=" << n_particles
                            << "-c=" << Molarity*1e9 << "nM"
                            << "-t_max=" << t_max
                            << "-t_step=" << t_step/1e-6 << "us"
                            << "/k_ex_max=" << k_ex_max;

    std::stringstream command; // stringstream
    command << "mkdir -p " << dir_name.str(); 

    std::string command_str(command.str()); // convert to string

    const char* command_cstr = command_str.c_str(); // convert to C type string
    std::system(command_cstr); 
 
    const unsigned int num_runs = cfg.getValueOfKey<unsigned int>("num_runs");   
    // Run loop    
    for (unsigned int run=0; run<num_runs; run++) 
    {
        std::cout << "run= " << run << std::endl;        
        std::stringstream counts_file_name_ss;    
        counts_file_name_ss << dir_name.str() << "/" << run << ".dat"; 
        std::cout << "data file: " << counts_file_name_ss.str() << std::endl;

        std::string counts_file_name_str(counts_file_name_ss.str()); // string
        const char* counts_file_name_cstr = counts_file_name_str.c_str(); // C type string

        std::ifstream counts_ifs(counts_file_name_cstr);    
        if (counts_ifs.good()) {
            std::cout << "Count file exists." << std::endl;
            continue;    
        } 
        else {
            counts_ifs.close();
        }
      
        counts[0] = molecules.arrivals(n_particles, k_ex_max, phi_f, phi_det, t_step, mdf, rng);
        // Main simulation loop
        for (size_t step=1; step < n_steps; step++) {
            molecules.update_positions_periodic (n_particles, box, rng, t_step);
            molecules.update_states (n_particles, rng, k_ex_max, phi_st, phi_sb, k_ts, k_tb, t_step, mdf);
            counts[step] = molecules.arrivals(n_particles, k_ex_max, phi_f, phi_det, t_step, mdf, rng);                
        } // End of main loop

        unsigned int counts_sum = std::accumulate(counts.begin(), counts.end(), 0.0);
        std::cout << "avg count rate= " << (counts_sum/t_max)/1e3 << " kHz" << std::endl; 
        std::cout << "num photons= " << counts_sum << std::endl; 
  
        std::ofstream counts_ofs (counts_file_name_cstr);
        if (counts_ofs.is_open()) {
            for (unsigned int i=0; i<n_steps; i++) {
                counts_ofs << counts[i] << std::endl;        
            }        
            counts_ofs.close();
        }
        else std::cout << "Unable to open file";
    } // End of run loop

    // Get the elapsed wall clock time in seconds
    clock_t end = clock();
    double elapsed_mins = double(end - begin) / CLOCKS_PER_SEC / 60.0;
    std::cout <<  "elapsed_mins= " << elapsed_mins << std::endl;

    return 0;
}
