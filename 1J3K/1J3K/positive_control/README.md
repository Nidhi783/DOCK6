**1J3K Positive Control – DOCK6**

**Overview**

This folder contains the positive-control molecular docking workflow for the 1J3K–WRA protein–ligand system using DOCK6.

The purpose of the positive control is to test whether DOCK6 can reproduce the binding of the native ligand WRA within its known crystallographic binding region.

The workflow starts from the original 1J3K.pdb structure, prepares the receptor and ligand, generates receptor spheres, selects spheres around the native WRA-binding site, creates the docking box and GRID, and finally performs flexible ligand docking using DOCK6.

**Workflow**

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
     ↓
Select spheres around native WRA
     ↓
SHOWBOX
     ↓
GRID
     ↓
DOCK6
     ↓
Docked WRA pose

**Main Input Files**

1J3K.pdb – original 1J3K structure

WRA_only.mol2 – prepared native WRA ligand

1J3K_chainA_receptor_WRAremoved_charged.mol2 – prepared Chain A receptor with WRA removed

1J3K_master_spheres.sph – master sphere set generated using SPHGEN

run_1J3K_complete_positive_control.sh – complete preparation-to-docking workflow

The native ligand WRA was prepared using AM1-BCC partial charges with a net charge of 0.

**SPHGEN Setup**

The receptor molecular surface used for sphere generation is:

1J3K_chainA_receptor_WRAremoved.dms

**The recovered SPHGEN input is:**
1J3K_chainA_receptor_WRAremoved.dms
R
X
0.0
4.0
1.4
1J3K_master_spheres.sph

This produces 1J3K_master_spheres.sph. For the positive control, spheres are selected around the native WRA ligand, so docking is performed within the known ligand-binding region.

**Sphere-Radius Testing**

The positive-control workflow uses even sphere-selection radii from 4 Å to 28 Å:
4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28 Å

**The SHOWBOX margin is defined as:**
SHOWBOX margin = sphere radius + 6 Å

**Examples:**

Sphere radius

SHOWBOX margin

4 Å

10 Å

10 Å

16 Å

20 Å

26 Å

28 Å

34 Å

GRID calculations use a spacing of 0.5 Å.

**Running the Workflow**

Executable Command Card
**Make the script executable:**

chmod +x run_1J3K_complete_positive_control.sh
Run the workflow for a selected radius:

./run_1J3K_complete_positive_control.sh 1J3K.pdb 10

**In this example:**

1J3K.pdb is the original input structure.

10 is the sphere-selection radius in Å.

To use another radius, replace 10 with the required value. For example:

./run_1J3K_complete_positive_control.sh 1J3K.pdb 20

This performs the positive-control workflow using a 20 Å sphere-selection radius.

**Output**

Each radius creates a separate output directory containing the main files generated during sphere selection, SHOWBOX, GRID and DOCK6 docking.

Typical output files include:

selected_spheres_radius<RADIUS>.sph
box_radius<RADIUS>.pdb
grid_radius<RADIUS>.out
dock_radius<RADIUS>.out
WRA_docked_radius<RADIUS>_scored.mol2

The scored MOL2 file contains the final docked WRA pose and its associated docking score.

**Docking Evaluation**

The positive-control docking results can be evaluated using:

GRID score – the DOCK6 docking score.

Heavy-atom RMSD – difference between the docked and crystallographic WRA poses.

Centroid distance – distance between the centres of the docked and crystallographic WRA ligands.

Lower RMSD and centroid-distance values indicate that the docked WRA pose is closer to the native crystallographic binding position.

**Reproducibility**

The complete script combines the major stages of the workflow from 1J3K structure preparation through DOCK6 docking.
The prepared receptor, native WRA ligand, molecular-surface/SPHGEN settings and master spheres are retained so that the positive-control docking workflow and its parameters can be reviewed and reproduced.
