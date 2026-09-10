# 1J3K Negative Control – DOCK6

## Overview

This folder contains the **negative-control molecular docking workflow for the 1J3K–WRA protein–ligand system using DOCK6**.

The purpose of the negative control is to test WRA docking at a deliberately selected site away from the native crystallographic WRA-binding region. This provides a comparison with the positive-control redocking at the known binding site.

## Negative-Control Site

The negative-control reference centre used in this workflow is:

```text
X = 39.053 Å
Y = 17.088 Å
Z = 2.043 Å
```

or:

```text
[39.053, 17.088, 2.043] Å
```

The file `1J3K_negative_site_reference.mol2` represents this fixed location. It is used by `sphere_selector` to select spheres around the negative-control site rather than around the native WRA ligand.

## Workflow

```text
1J3K.pdb
     ↓
Select Chain A
     ↓
Separate receptor and native WRA
     ↓
Prepare receptor and WRA
     ↓
WRA: AM1-BCC charges, net charge 0
     ↓
Generate receptor molecular surface
     ↓
SPHGEN
     ↓
1J3K_master_spheres.sph
       +
1J3K_negative_site_reference.mol2
     ↓
Select spheres around negative-control site
     ↓
SHOWBOX
     ↓
GRID
     ↓
DOCK6
     ↓
Negative-control docked WRA pose
```

## Main Input Files

- `1J3K.pdb` – original 1J3K structure
- `WRA_only.mol2` – prepared native WRA ligand
- `1J3K_chainA_receptor_WRAremoved_charged.mol2` – prepared Chain A receptor with WRA removed
- `1J3K_master_spheres.sph` – master sphere set generated using SPHGEN
- `1J3K_negative_site_reference.mol2` – reference defining the negative-control site
- `run_1J3K_complete_negative_control.sh` – complete preparation-to-negative-control docking workflow

The native ligand **WRA** was prepared using **AM1-BCC partial charges** with a **net charge of 0**.

## SPHGEN Setup

The receptor molecular surface used for sphere generation is:

```text
1J3K_chainA_receptor_WRAremoved.dms
```

The recovered SPHGEN input is:

```text
1J3K_chainA_receptor_WRAremoved.dms
R
X
0.0
4.0
1.4
1J3K_master_spheres.sph
```

This generates `1J3K_master_spheres.sph`. For the negative control, spheres from this master set are selected around `1J3K_negative_site_reference.mol2`.

## Sphere-Radius Testing

The negative-control workflow tests the following sphere-selection radii:

```text
4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28 Å
```

The SHOWBOX margin is defined as:

```text
SHOWBOX margin = sphere radius + 6 Å
```

Examples:

| Sphere radius | SHOWBOX margin |
|---:|---:|
| 4 Å | 10 Å |
| 10 Å | 16 Å |
| 20 Å | 26 Å |
| 28 Å | 34 Å |

GRID calculations use a spacing of **0.5 Å**.

## Running the Workflow

### Executable Command Card

Make the complete script executable:

```bash
chmod +x run_1J3K_complete_negative_control.sh
```

Run the complete negative-control workflow:

```bash
./run_1J3K_complete_negative_control.sh 1J3K.pdb
```

Here, `1J3K.pdb` is the original input structure.

The complete script performs the preparation stages and then runs the negative-control radius series from **4 Å to 28 Å**.

## Output

Each radius creates a separate output directory, for example:

```text
radius10_negative_full/
```

Typical output files include:

```text
selected_spheres_radius<RADIUS>.sph
box_radius<RADIUS>.pdb
grid_radius<RADIUS>.out
dock_radius<RADIUS>.out
WRA_negative_docked_radius<RADIUS>_scored.mol2
```

The complete workflow also produces:

```text
1J3K_negative_control_results.tsv
```

This results table records the radius, SHOWBOX margin, number of selected spheres, GRID score, heavy-atom RMSD, centroid distance and run status.

## Docking Evaluation

The negative-control docking results are evaluated using:

- **GRID score** – the DOCK6 docking score.
- **Heavy-atom RMSD** – difference between the docked WRA pose and the native crystallographic WRA pose.
- **Centroid distance** – distance between the centres of the docked and native WRA ligands.

These measurements allow the negative-control results to be compared with the positive-control redocking results.

## Reproducibility

The complete script combines the major stages from **1J3K structure preparation through negative-control DOCK6 docking**.

The file `1J3K_negative_site_reference.mol2` is intentionally retained as an input because it defines the exact negative-control location used in this project:

```text
[39.053, 17.088, 2.043] Å
```

Keeping the prepared inputs, SPHGEN settings, negative-site reference and executable workflow together allows the negative-control procedure and its parameters to be reviewed and reproduced.
