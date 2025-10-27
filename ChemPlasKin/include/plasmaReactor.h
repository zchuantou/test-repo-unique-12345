// ChemPlasReactor class. Generate a unified gas-plasma ODE system

#ifndef CHEMPLASKIN_PLASMAREACTOR_H
#define CHEMPLASKIN_PLASMAREACTOR_H

#include "cantera/numerics/FuncEval.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Boltzmann.h"

namespace Cantera
{
class ChemPlasReactor : public FuncEval {
public:
    /*
     * Constructor
     * @param[in] sol Solution object specifying initial system state.
     */
    ChemPlasReactor
    (
            shared_ptr<Solution> sol,
            std::map<std::string, double>& BoltzmannSpecies,
            std::string plasmaHeatModel = " "
    )
    : m_BoltzmannSpecies(BoltzmannSpecies)
    {
        /* ---------------------- INITIALIZE MEMBER VARS ---------------------- */
        // Pointer to system's ThermoPhase object, kinetics manager and Transport object.
        // Updated by the solver during simulation to provide iteration-specific information
        m_gas = sol->thermo();
        m_kinetics = sol->kinetics();
        m_transport = sol->transport();

        // Record pressure and density at initial state, which can be reset using setConstPD()
        m_pressure = m_gas->pressure(); // For constant pressure case (CP)
        m_density = m_gas->density(); // For constant volume case (CV)

        // Number of chemical species in the system.
        m_nSpecies = m_gas->nSpecies();
        std::cout <<"Creating ChemPlasReactor: " << m_nSpecies  << " species and "
        << m_kinetics->nReactions() << " reactions." << std::endl;

        // Resize vector<double> storage containers for
        // species partial molar enthalpies/internal energy
        // and net production rates. Updated and used by the solver per iteration.
        m_hbar.resize(m_nSpecies);
        m_ubar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);
        m_q_net.resize(m_kinetics->nReactions());
        m_deltaH.resize(m_kinetics->nReactions());

        // Number of equations in the ODE system.
        m_nEqs = m_nSpecies + 4;

        electronIndex = findSpeciesIndex({"E", "electron", "Electron", "e"});
        N2Index = findSpeciesIndex("N2");
        NIndex = findSpeciesIndex("N");
        NH3Index = findSpeciesIndex("NH3");

        // Specify plasma heating model.
        // Detailed: calculating from each reaction steps.
        // Flitti & Pancheshnyi 2009 model: global modelling
        detailPlasmaHeatModel = (plasmaHeatModel != "Flitti");

        // Register energy transfer data for plasma reactions
        setPlasmaEnergyTransfer(sol->reaction_section);

        // Create a map for all species
        for (size_t i = 0; i < m_gas->nSpecies(); i++) {
            allSpeciesIndices[m_gas->speciesName(i)] = i;
        }

        // Register species configured by CppBOLOS
        for (const auto& [sp, _] : m_BoltzmannSpecies) {
            bSpeciesIndices[sp] = findSpeciesIndex(sp);
            allSpeciesIndices.erase(sp);  // Remove the species from all-species map
        }
    }

    /**
     * Evaluate the ODE right-hand-side function, ydot = f(t,y).
     *   - overridden from FuncEval, called by the integrator during simulation.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] p sensitivity parameter vector, length nparams()
     *   - note: sensitivity analysis isn't implemented in this example
     */
    void eval(double t, double* y, double* ydot, double* p) override {
        // Time derivative of solution vector [T, he, Ep, e_vib(N2), Y_1, Y_2, ... Y_K],
        // *ydot*, are defined for clear and convenient access to these vectors:
        nSteps ++;
        double temperature = y[0];
        m_he = y[1];
        m_Ep = y[2];
        eVibN2 = y[3];
        double *massFracs = &y[4];

        // Apply external electron number density profile
        if (imposedNe) {
            massFracs[electronIndex] = m_N_e*1e6/Avogadro*m_gas->molecularWeight(electronIndex) / m_density;
        }

        double *dTdt = &ydot[0];
        double *HE_dot = &ydot[1];
        double *Ep_dot = &ydot[2];
        double *eVibN2_dot = &ydot[3];
        double *dYdt = &ydot[4];

        /* ------------------------- UPDATE GAS STATE ------------------------- */
        // The state of the ThermoPhase is updated to reflect the current solution
        // vector, which was calculated by the integrator.
        m_gas->setMassFractions(massFracs);
        if (constPressure) {
            m_gas->setState_TP(temperature, m_pressure); // Constant pressure
        } else {
            m_gas->setState_TD(temperature, m_density); // Constant volume
        }

        /* ------------------------- GET THERMO PROPS ------------------------- */
        m_gas->getPartialMolarEnthalpies(m_hbar.data());
        m_gas->getPartialMolarIntEnergies(m_ubar.data());

        /* ------------------------ GET KINETIC PROPS ------------------------ */
        m_kinetics->getNetProductionRates(m_wdot.data());    // [kmol/m^3/s]
        m_kinetics->getNetRatesOfProgress(m_q_net.data());  // [kmol/m^3/s]
        m_kinetics->getDeltaEnthalpy(m_deltaH.data());      // [J/kmol]

        /* ---------------------- ENERGY EQUATION ---------------------- */
        double sum_cpYdot = 0.0;
        double sum_H_wdot = 0.0;
        for (size_t i = 0; i < m_nSpecies; i++) {
            sum_cpYdot += m_gas->cp_mass(i) * dYdt[i];
            sum_H_wdot += m_hbar[i] * m_wdot[i];
        }

        // Temperature equation
        *dTdt = (sum_H_wdot - sum_cpYdot * temperature) / (m_gas->density() * m_gas->cp_mass());

        // Electron energy equation
        *HE_dot = calculateElectronEnergyBalance(m_wdot, m_q_net);

        // Deposited plasma energy equation
        *Ep_dot = calculatePlasmaEnergyDeposition(m_q_net);

        // N2 vibrational energy equation
        *eVibN2_dot = calculateN2VibrationalEnergy(m_q_net);

        /* ----------------------- SPECIES EQUATIONS ----------------------- */
        for (size_t i = 0; i < m_nSpecies; i++) {
            if (i != electronIndex) {
                dYdt[i] = m_wdot[i] * m_gas->molecularWeight(i) / m_gas->density();
            }
        }

        // Electron density equation
        dYdt[electronIndex] = calculateElectronDensity(m_wdot, m_q_net);
    }

    /**
     * Number of equations in the ODE system.
     */
    size_t neq() const override {
        return m_nEqs;
    }

    /**
     * Set the initial conditions for the ODE system.
     */
    void getInitialConditions(double t0, size_t leny, double* y) override {
        // Get initial mass fractions
        m_gas->getMassFractions(&y[4]);
        
        // Set initial conditions for additional variables
        y[0] = m_gas->temperature();  // Temperature
        y[1] = 0.0;                   // Electron energy
        y[2] = 0.0;                   // Deposited plasma energy
        y[3] = 0.0;                   // N2 vibrational energy
    }

    /**
     * Set whether the reactor operates at constant pressure.
     */
    void setConstPressure(bool const_p = true) {
        constPressure = const_p;
    }

    /**
     * Set constant pressure and density.
     */
    void setConstPD(double P, double rho) {
        m_pressure = P;
        m_density = rho;
    }

    /**
     * Get deposited plasma energy.
     */
    double depositedPlasmaEnergy() const {
        return m_Ep;
    }

    /**
     * Reset deposited plasma energy.
     */
    void resetDepositedPlasmaEnergy() {
        m_Ep = 0.0;
    }

    /**
     * Set external electron number density profile.
     */
    void setExternalElectronDensity(double n_e) {
        m_N_e = n_e;
        imposedNe = true;
    }

    /**
     * Get the electron density.
     */
    double getElectronDensity() const {
        return m_N_e;
    }

private:
    // Helper functions
    double calculateElectronEnergyBalance(const vector<double>& wdot, const vector<double>& q_net) {
        double he_dot = 0.0;
        for (size_t i = 0; i < m_kinetics->nReactions(); i++) {
            if (e_ext_map.find(i) != e_ext_map.end()) {
                he_dot += q_net[i] * e_ext_map[i];
            }
        }
        return he_dot;
    }

    double calculatePlasmaEnergyDeposition(const vector<double>& q_net) {
        double Ep_dot = 0.0;
        for (size_t i = 0; i < m_kinetics->nReactions(); i++) {
            if (e_ext_map.find(i) != e_ext_map.end()) {
                Ep_dot += q_net[i] * e_ext_map[i];
            }
        }
        return Ep_dot;
    }

    double calculateN2VibrationalEnergy(const vector<double>& q_net) {
        double eVibN2_dot = 0.0;
        for (size_t i = 0; i < m_kinetics->nReactions(); i++) {
            if (evib_N2_map.find(i) != evib_N2_map.end()) {
                eVibN2_dot += q_net[i] * evib_N2_map[i];
            }
        }
        return eVibN2_dot;
    }

    double calculateElectronDensity(const vector<double>& wdot, const vector<double>& q_net) {
        // Calculate electron density change
        if (electronIndex != npos) {
            return wdot[electronIndex] * m_gas->molecularWeight(electronIndex) / m_gas->density();
        }
        return 0.0;
    }

    // Find species index with multiple possible names
    size_t findSpeciesIndex(const std::vector<std::string>& names) {
        for (const auto& name : names) {
            size_t index = m_gas->speciesIndex(name);
            if (index != npos) {
                return index;
            }
        }
        return npos;
    }

    size_t findSpeciesIndex(const std::string& name) {
        size_t index = m_gas->speciesIndex(name);
        if (index != npos) {
            return index;
        }
        std::cout << "ChemPlasReactor::findSpeciesIndex ERROR: Register species index of: " << name << std::endl;
        throw std::runtime_error("No valid species name found! Consider adding as a dummy species.");
    }

private:
    // Private member variables, to be used internally.
    shared_ptr<ThermoPhase> m_gas;
    shared_ptr<Kinetics> m_kinetics;
    shared_ptr<Transport> m_transport;
    vector<double> m_hbar;
    vector<double> m_ubar; // internal energy
    vector<double> m_wdot; // net production rate of species, kmol/m^3/s, Length: nSpecies
    vector<double> m_q_net; // net production rate of individual reactions, kmol/m^3/s, Length: nReactions()
    vector<double> m_deltaH;
    double m_pressure;
    size_t m_nSpecies;
    size_t m_nEqs;
    double m_density; // fixed density for constant volume case
    double m_Ep = 0.0; // Deposited energy. (J/m^3)
    double m_he = 0.0;
    size_t N2Index, O2Index, H2Index, OIndex;
    size_t electronIndex, NIndex, NH3Index;

    std::map<std::string, double>& m_BoltzmannSpecies;
    std::map<std::string, size_t> allSpeciesIndices;
    std::map<std::string, size_t> bSpeciesIndices; // indices of species configured by CppBOLOS

    double m_N_e = 0.0; // Number density of e-
    bool imposedNe = false;
    bool constPressure = true;

    bool detailPlasmaHeatModel = true;

    // Maps to store reaction number and corresponding energy_transfer values
    std::map<size_t, double> e_th_map; // energy transfer in reaction steps, can be negative
    std::map<size_t, double> e_ext_map; // Externally deposited plasma energy for Boltzmann reactions, must be positive
    // Plasma energy stored in vibration states of species: N2, O2, NO, NH3, etc.
    std::map<size_t, double> evib_N2_map;
    std::map<size_t, double> evib_O2_map;
    std::map<size_t, double> evib_NO_map;
    std::map<size_t, double> evib_NH3_map;

    double eVibN2 = 0.0;
    size_t nSteps = 0;

    // Register plasma energy transfer data
    void setPlasmaEnergyTransfer(const vector<AnyMap>& reactions) {
        for(size_t i = 0; i < reactions.size(); i++) {
            auto& reaction = reactions[i];
            if(reaction.hasKey("energy_transfer")) {
                auto& energy_transfer = reaction["energy_transfer"].as<AnyMap>();

                // Check for e_th
                if(energy_transfer.hasKey("e_th")) {
                    double e_th = parseEnergyString(energy_transfer["e_th"].as<std::string>());
                    e_th_map[i] = e_th;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = e_th;
                    }
                }

                // Check for evib_N2
                if(energy_transfer.hasKey("evib_N2")) {
                    double evib_N2 = parseEnergyString(energy_transfer["evib_N2"].as<std::string>());
                    evib_N2_map[i] = evib_N2;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = evib_N2;
                    }
                }

                // Check for evib_O2
                if(energy_transfer.hasKey("evib_O2")) {
                    double evib_O2 = parseEnergyString(energy_transfer["evib_O2"].as<std::string>());
                    evib_O2_map[i] = evib_O2;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = evib_O2;
                    }
                }

                // Check for evib_NO
                if(energy_transfer.hasKey("evib_NO")) {
                    double evib_NO = parseEnergyString(energy_transfer["evib_NO"].as<std::string>());
                    evib_NO_map[i] = evib_NO;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = evib_NO;
                    }
                }

                // Check for evib_NH3
                if(energy_transfer.hasKey("evib_NH3")) {
                    double evib_NH3 = parseEnergyString(energy_transfer["evib_NH3"].as<std::string>());
                    evib_NH3_map[i] = evib_NH3;
                    if (reaction.hasKey("type") && reaction["type"] == "Boltzmann") {
                        e_ext_map[i] = evib_NH3;
                    }
                }
            }
        }
    }

    // Utility function to parse custom energy format
    static double parseEnergyString(const std::string& energyString) {
        // Find the position of "_eV"
        size_t pos = energyString.find("_eV");
        if(pos != std::string::npos) {
            // Convert part before "_eV" to double
            return std::stod(energyString.substr(0, pos));
        } else {
            // Throw an exception or return a sentinel value, based on your preference
            throw std::runtime_error("Invalid energy format: " + energyString);
        }
    }

};

}

#endif //CHEMPLASKIN_PLASMAREACTOR_H