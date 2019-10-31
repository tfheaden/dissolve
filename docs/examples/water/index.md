---
title: Liquid Water
parent: Examples
has_children: true
nav_order: 3
---
# Liquid Water

30 minute read
{: .label .label-blue }

## Summary

- Run, refine, and analyse the structure of liquid water at 298 K
- Three neutron datasets, measured on the SANDALS diffractometer at the ISIS Pulsed Neutron and Muon Source in 2001
- Reference input file: [water.txt](https://github.com/trisyoungs/dissolve/tree/develop/examples/water/water.txt)

### Data Files
- H<sub>2</sub>O: [SLS18498-H2O.mdcs01](https://github.com/trisyoungs/dissolve/tree/develop/examples/water/data/SLS18498-H2O.mdcs01)
- HDO: [SLS18500-HDO5050.mdcs01](https://github.com/trisyoungs/dissolve/tree/develop/examples/water/data/SLS18500-HDO5050.mdcs01)
- D<sub>2</sub>O: [SLS18502-D2O.mdcs01](https://github.com/trisyoungs/dissolve/tree/develop/examples/water/data/SLS18502-D2O.mdcs01)

## Summary

TODO

Steps:
- Set up a suitable species and configuration
- Equilibrate the system
- Adjust intramolecular geometry to match experiment
- Refine the interatomic potentials
- Analyse the liquid structure of the refined simulation

## Preparation

Create a new directory in a suitable location and download or copy the data files listed above into it.

> Open the Dissolve GUI
{: .action .action_dissolve}
> Choose **Create a New Simulation from Scratch**
{: .action}
> Save As...
{: .action .action_menu}
> Save your input file under a sensible name in the directory containing the example data file
{: .action}

[Let's begin!](step1.md){: .btn }
