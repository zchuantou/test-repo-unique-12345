/*--------------------------------*- C++ -*----------------------------------*\
|     ____ _                    ____  _           _  ___                      |
|    / ___| |__   ___ _ __ ___ |  _ \| | __ _ ___| |/ (_)_ __                 |
|   | |   | '_ \ / _ \ '_ ` _ \| |_) | |/ _` / __| ' /| | '_ \                |
|   | |___| | | |  __/ | | | | |  __/| | (_| \__ \ . \| | | | |               |
|    \____|_| |_|\___|_| |_| |_|_|   |_|\__,_|___/_|\_\_|_| |_|               |
|                                                                             |
|   A Freeware for Unified Gas-Plasma Kinetics Simulation                     |
|   Version:      1.0.0 (July 2024)                                           |
|   License:      GNU LESSER GENERAL PUBLIC LICENSE, Version 2.1              |
|   Author:       Xiao Shao                                                   |
|   Organization: King Abdullah University of Science and Technology (KAUST)  |
|   Contact:      xiao.shao@kaust.edu.sa                                      |
\*---------------------------------------------------------------------------*/

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/ext/bolos/Logger.h"
#include "cantera/kinetics/Boltzmann.h"
#include <iostream>
#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"
#include <cmath>
#include <fstream>
#include <string>
#include "plasmaReactor.h"
#include "readParameters.h"
#include "cantera/kinetics/ReactionPath.h"
#include "cantera/zerodim.h"
#include "cantera/kinetics.h"
#include <iomanip>
#include "ROPIntegrator.h"
#include "SensitivityAnalysis.h"
#include "ROPIntegratorAll.h"

using namespace Cantera;

/* ------------------------ PREPARE FEW USEFUL FUNCTIONS ------------------------ */
// Get number density of species [#/cm^3]
double getNumberDens(const shared_ptr<ThermoPhase>& gas, const size_t i ){
    return 1e-6 * gas->moleFraction(i) * Avogadro * gas->molarDensity();
}

// Logarithmically changed timestep
// K: grow/decay order | t_N: cover period of 0-t_N | Nsteps: assigned steps | n: step index
double logDynamicTimestep(const double& K, const double& t_N, const int& Nsteps, const int& n ) {
    return (t_N / K) * log10( (Nsteps + (n+1)*(pow(10, K) - 1)) /
                              (Nsteps + n*(pow(10, K) - 1)) );
}

// Triangle wave function with gap for power density (range: [0, A])
double triangle_wave_with_gap(double t, double T_wave, double T_gap, double A) {
    double T_period = T_wave + T_gap;
    double t_mod = fmod(t, T_period);  // Map to one period

    if (t_mod < T_wave) {
        // In effective wave interval, calculate triangle wave
        if (t_mod < T_wave / 2) {
            return (2 * A / T_wave) * t_mod;                  // Rising edge
        } else {
            return (2 * A / T_wave) * (T_wave - t_mod);       // Falling edge
        }
    } else {
        // In gap region, output 0
        return 0.0;
    }
}

// Calculate reduced electric field from power density
double calculate_EN_from_power(double power_density, double electron_density, double gas_density) {
    const double e = ElectronCharge;
    double mu_e = BoltzmannRate::bsolver.mobility();
    return sqrt(power_density / (e * electron_density * (mu_e/(gas_density)))) / (gas_density) * 1e18; // Return EN in Td (V·cm²)
}

// Reaction path analysis function
void reactionPathAnalysis(Cantera::Solution& soln, const std::string& element) {
    auto kin = soln.kinetics();
    Cantera::ReactionPathBuilder builder;
    Cantera::ReactionPathDiagram diagram;
    diagram.title = "Reaction path analysis for element: " + element;
    diagram.threshold = 1e-4;
    std::ostringstream null_log;
    builder.init(null_log, *kin);
    builder.build(*kin, element, null_log, diagram, true);
    std::ofstream file("../examples/NH3H2AIR/" + element + ".dot");
    diagram.exportToDot(file);
    file.close();
    std::cout << "DOT file for " << element << " path exported." << std::endl;
}

int main(int argc, char *argv[]) {
    printHeader();
    
    // Parse command line arguments
    std::string casePath = "../examples/H2O2He";
    std::string logLevel = "INFO";
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-case" && i + 1 < argc) {
            casePath = argv[++i];
        } else if (arg == "-log" && i + 1 < argc) {
            logLevel = argv[++i];
        }
    }
    
    // Set logging level
    if (logLevel == "NONE") {
        Logger::disable();
    } else if (logLevel == "DEBUG") {
        Logger::enable();
        Logger::setLogLevel(Logger::DEBUG);
    } else if (logLevel == "WARNING") {
        Logger::enable();
        Logger::setLogLevel(Logger::WARNING);
    } else {
        Logger::enable();
        Logger::setLogLevel(Logger::INFO);
    }
    
    std::cout << "Case path: " << casePath << std::endl;
    std::cout << "Log level: " << logLevel << std::endl;
    
    try {
        // Load gas solution from mechanism file
        std::string mechFile = casePath + "/chemPlasProperties.yaml";
        auto gas = newSolution(mechFile);
        
        // Set initial conditions
        gas->temperature().set(300.0);  // Initial temperature [K]
        gas->pressure().set(OneAtm);    // Initial pressure [Pa]
        
        // Print initial state
        std::cout << "Initial temperature: " << gas->temperature() << " K" << std::endl;
        std::cout << "Initial pressure: " << gas->pressure() / OneAtm << " atm" << std::endl;
        
        // Create reactor and integrator
        PlasmaReactor reactor;
        reactor.setSolution(*gas);
        
        auto integrator = newIntegrator("CVODE");
        integrator->setMaxStepSize(1e-3);
        
        // Time integration parameters
        double t_end = 1e-3;  // End time [s]
        double dt = 1e-6;      // Initial timestep [s]
        double t = 0.0;        // Current time
        
        std::cout << "Starting simulation..." << std::endl;
        
        // Main time loop
        while (t < t_end) {
            // Integrate one timestep
            integrator->integrate(t + dt);
            t += dt;
            
            // Output progress every 100 steps
            if (static_cast<int>(t / dt) % 100 == 0) {
                std::cout << "Time: " << t * 1e6 << " us, T: " << gas->temperature() << " K" << std::endl;
            }
            
            // Adaptive timestep (simple implementation)
            if (gas->temperature() > 1000.0) {
                dt *= 0.9;  // Reduce timestep for high temperatures
            }
        }
        
        std::cout << "Simulation completed successfully!" << std::endl;
        std::cout << "Final temperature: " << gas->temperature() << " K" << std::endl;
        
    } catch (CanteraError& err) {
        std::cerr << "Cantera error: " << err.what() << std::endl;
        return 1;
    } catch (std::exception& err) {
        std::cerr << "Error: " << err.what() << std::endl;
        return 1;
    }
    
    return 0;
}