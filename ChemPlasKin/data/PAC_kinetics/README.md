# Plasma-Assisted Combustion Kinetics

This directory contains unified gas-plasma kinetic mechanisms in YAML format.

## Directory Structure

```
PAC_kinetics/
├── ZDPlasKin_kinetics/
│   ├── Mao-H2O2He/
│   │   └── H2O2HE_Mao.yaml
│   └── kineticsParser/
│       └── parsePlasKin.py
└── ...
```

## Mechanism Files

The kinetic mechanisms are provided in YAML format, compatible with Cantera:
- Unified gas-phase and plasma chemistry
- Electron impact reactions
- Energy transfer processes
- Vibrational excitation

## Parser Tool

A Python parser is available to convert ZDPlasKin input mechanisms to YAML format:

```bash
cd ChemPlasKin/data/PAC_kinetics/ZDPlasKin_kinetics/kineticsParser
python parsePlasKin.py --input "plasmaH2O2.inp" --output "parsedPlasKin.yaml"
```

## Usage

The mechanism files are referenced in the `chemPlasProperties` files under the `mechFile` parameter:

```bash
mechFile "<case>/../../data/PAC_kinetics/ZDPlasKin_kinetics/Mao-H2O2He/H2O2HE_Mao.yaml";
```