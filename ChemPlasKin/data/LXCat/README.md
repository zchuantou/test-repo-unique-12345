# LXCat Cross-Section Data

This directory contains electron collision cross-section data from LXCat database in the format required by CppBOLOS.

## File Format

The cross-section files should follow the LXCat format:
- Energy-dependent collision cross sections
- Elastic and inelastic processes
- Excitation, ionization, and attachment cross sections

## Available Data

- `bolsigdb_H2O2HE.dat`: Cross-section data for H2-O2-He system

## Data Source

LXCat database: https://nl.lxcat.net/home/

## Usage

The cross-section files are referenced in the `chemPlasProperties` files under the `csDataFile` parameter.

```bash
csDataFile "<case>/../../data/LXCat/bolsigdb_H2O2HE.dat";
```