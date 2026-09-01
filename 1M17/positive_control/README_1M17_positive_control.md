# 1M17 Positive-Control DOCK6 Workflow

## Overview

This repository contains the workflow used for positive-control
molecular docking of the **1M17 protein--ligand system** using
**DOCK6**. The native ligand used for the positive control is **AQ4**.

The validated docking workflow evaluates sphere radii from **4 Å to 28 Å
in 2 Å increments**. For each radius, the workflow performs sphere
selection, SHOWBOX generation, GRID generation, DOCK6 docking, and
extraction of the final docking metrics.

The final metrics reported are:

-   Grid Score (kcal/mol)
-   Heavy-atom RMSD (Å)
-   Centroid-to-centroid distance (Å)

## Repository Files

``` text
1M17_positive_control/
├── 1M17.pdb
├── AQ4_prepared.mol2
├── 1M17_receptor_prepared.mol2
├── 1M17_master_spheres.sph
├── run_1M17_complete_workflow.sh
├── run_1M17_positive_control_4_to_28.sh
├── 1M17_positive_control_results.tsv
└── README.md
```

### File descriptions

**`1M17.pdb`**\
Raw 1M17 structure used as the starting structure for the complete
workflow.

**`AQ4_prepared.mol2`**\
Prepared native AQ4 ligand used in the validated docking workflow.

**`1M17_receptor_prepared.mol2`**\
Prepared 1M17 receptor used for GRID generation and docking.

**`1M17_master_spheres.sph`**\
Master receptor spheres used for radius-based sphere selection.

**`run_1M17_positive_control_4_to_28.sh`**\
Validated docking script used to perform the 1M17 positive-control
radius series.

**`run_1M17_complete_workflow.sh`**\
Start-to-end automation intended to extend the workflow from raw
structure preparation through sphere generation and docking. See the
validation note below.

**`1M17_positive_control_results.tsv`**\
Final results table generated from the validated positive-control
docking workflow.

## Validated Docking Workflow

The validated docking script starts from the prepared receptor, prepared
AQ4 ligand, and master sphere file.

``` text
AQ4_prepared.mol2
        +
1M17_receptor_prepared.mol2
        +
1M17_master_spheres.sph
        |
        v
Sphere selection
        |
        v
SHOWBOX
        |
        v
GRID
        |
        v
DOCK6
        |
        v
Grid Score + Heavy-atom RMSD + Centroid Distance
```

The sphere radii evaluated are:

``` text
4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28 Å
```

For each radius, the SHOWBOX margin is:

``` text
margin = sphere radius + 6 Å
```

## Running the Validated Docking Workflow

Make the script executable:

``` bash
chmod +x run_1M17_positive_control_4_to_28.sh
```

Run the workflow:

``` bash
./run_1M17_positive_control_4_to_28.sh
```

The script creates a separate directory for each radius and writes the
final summary to:

``` text
1M17_positive_control_results.tsv
```

## Preparation Information

The original receptor and ligand preparation was performed manually
using **UCSF Chimera**.

For AQ4:

-   Charge method: **AM1-BCC**
-   Formal/net charge: **0**

The receptor molecular surface was used for SPHGEN sphere generation.

The SPHGEN input settings used for 1M17 were:

``` text
1M17_receptor.dms
R
X
0.0
4.0
1.4
1M17_receptor.sph
```

The resulting receptor spheres were used to generate the master sphere
file required by the docking workflow.

## Complete Automated Workflow

The repository also contains:

``` text
run_1M17_complete_workflow.sh
```

This script is intended to automate the workflow from the raw `1M17.pdb`
structure through preparation, molecular-surface/sphere generation, and
the positive-control docking series.

Conceptually:

``` text
1M17.pdb
    |
    v
Receptor and AQ4 preparation
    |
    v
Prepared receptor + prepared AQ4
    |
    v
Receptor molecular surface
    |
    v
SPHGEN
    |
    v
Master spheres
    |
    v
Sphere selection (4–28 Å)
    |
    v
SHOWBOX
    |
    v
GRID
    |
    v
DOCK6
    |
    v
Final results
```

### Important validation note

The **docking-only script has been run and validated** using the
prepared input files supplied in this repository.

The **complete automated preparation-to-docking script has not been
validated end-to-end on the original server**, because UCSF Chimera was
not available in that server environment. It is therefore provided as an
intended reproducible automation of the complete workflow and should be
validated in an environment containing all required dependencies before
its outputs are treated as equivalent to the manually prepared workflow.

## Software Requirements

The validated docking stage requires the relevant DOCK6 programs and
parameter files, including:

-   DOCK6
-   `sphere_selector`
-   `showbox`
-   `grid`
-   SPHGEN where sphere regeneration is required
-   Python 3
-   NumPy

The complete preparation-to-docking workflow additionally requires the
molecular-preparation and surface-generation software referenced by the
complete workflow script, including UCSF Chimera where applicable.

Users should check the executable paths and DOCK6 parameter-directory
path in the scripts before running them on another system.

## Output

The main results file is:

``` text
1M17_positive_control_results.tsv
```

It contains:

``` text
Radius
Margin
GridScore
HeavyAtomRMSD
CentroidDistance
Status
```

Heavy-atom RMSD and centroid distance are calculated relative to the
prepared native AQ4 ligand. Hydrogens are excluded from these
calculations.

## Reproducibility

Two routes are provided:

1.  **Validated docking reproduction** --- use the supplied prepared
    receptor, prepared AQ4 ligand, master spheres, and
    `run_1M17_positive_control_4_to_28.sh`.
2.  **Complete workflow automation** --- use `1M17.pdb` with
    `run_1M17_complete_workflow.sh` in an environment containing all
    required dependencies.

The first route corresponds to the docking workflow that was actually
executed and validated. The second route is included to support future
end-to-end automation of the preparation and docking workflow.
