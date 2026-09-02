**3ERT Positive-Control DOCK6 Workflow**

**Overview**

This folder contains the 3ERT positive-control molecular docking workflow using DOCK6. The native ligand used for redocking is OHT.

The repository includes the prepared inputs used in the validated docking work, the original radius-based docking workflow, and a merged start-to-end script that incorporates the initial receptor/ligand preparation and sphere-generation stages before the original docking procedure.

Recommended Repository Structure

3ERT/
└── positive_control/
    ├── README.md
    ├── 3ERT.pdb
    ├── OHT_prepared.mol2
    ├── 3ERT_receptor_prepared.mol2
    ├── 3ERT_master_spheres.sph
    ├── run_3ERT_positive_control.sh
    └── run_3ERT_complete_positive_control.sh

Keep the original filename of the validated docking script if it differs from run_3ERT_positive_control.sh.

**System**

Protein structure: 3ERT

Native/reference ligand: OHT

Docking software: DOCK6

Positive-control approach: redocking of the native OHT ligand

Sphere radii supported by the original script: 4–40 Å in 2 Å increments

SHOWBOX margin: sphere radius + 6 Å

Original Preparation

The original receptor and ligand preparation was performed manually using UCSF Chimera.

OHT preparation

The prepared OHT ligand uses:

Charge method: AM1-BCC
Formal/net charge: +1

The +1 net charge was confirmed from the existing OHT_prepared.mol2 file, for which the summed partial atomic charges were approximately +0.9999.

**Receptor preparation**

The receptor preparation includes hydrogen addition and Amber-based charge/type assignment before use in GRID and DOCK6.
The prepared receptor used by the docking workflow is: 3ERT_receptor_prepared.mol2

Molecular Surface and SPHGEN
The receptor molecular-surface file recovered from the original 3ERT work is:  3ERT_receptor.ms

**The original INSPH settings were:**

3ERT_receptor.ms
R
X
0.0
4.0
1.4
3ERT_master_spheres.sph

SPHGEN therefore generates: 3ERT_master_spheres.sph. Which is subsequently used for radius-based sphere selection.

**Validated Docking Workflow**

The original docking workflow starts from:
3ERT_master_spheres.sph
OHT_prepared.mol2
3ERT_receptor_prepared.mol2

and performs:

Master spheres
      |
      v
sphere_selector
      |
      v
selected_spheres_radius<RADIUS>.sph
      |
      v
SHOWBOX
      |
      v
box_radius<RADIUS>.pdb
      |
      v
GRID
      |
      v
DOCK6
      |
      v
OHT_docked_radius<RADIUS>_scored.mol2

**Radius and Box Settings**

The original script accepts even-numbered sphere radii from:

4, 6, 8, 10, ..., 40 Å

**For each radius:**

MARGIN = RADIUS + 6 Å

**For example:**

Sphere radius = 10 Å
SHOWBOX margin = 16 Å

**GRID Settings**

The original docking script uses the following key GRID settings:

grid_spacing                   0.5
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     all
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
allow_non_integral_charges     yes
bump_filter                    no
vdw_definition_file            vdw_AMBER_parm99.defn

The prepared 3ERT receptor is used as the receptor input.

**DOCK6 Settings**

The original positive-control script uses flexible ligand docking.
Key settings include:

conformer_search_type          flex
calculate_rmsd                 yes
use_rmsd_reference_mol         yes
orient_ligand                  yes
automated_matching             yes
max_orientations               5000
flexible_ligand                yes
grid_score_primary             yes
minimize_ligand                yes
num_scored_conformers          1

The prepared native OHT ligand is also used as the RMSD reference molecule.

**Running the Original Docking Workflow**
The original radius-based script is run with one radius at a time.

Example:
chmod +x run_3ERT_positive_control.sh
./run_3ERT_positive_control.sh 10

The script creates:
radius10_full/   and generates the sphere-selection, SHOWBOX, GRID, docking, and scored ligand files for that radius.

Complete Start-to-End Workflow

**The merged script is:**

run_3ERT_complete_positive_control.sh

It is designed to extend the workflow back to the raw PDB structure:

3ERT.pdb
    |
    v
Extract receptor + native OHT
    |
    v
UCSF Chimera preparation
    |
    +-- receptor preparation
    |
    +-- OHT: AM1-BCC, net charge +1
    |
    v
Prepared receptor + prepared OHT
    |
    v
Molecular-surface generation
    |
    v
3ERT_receptor.ms
    |
    v
SPHGEN
    |
    v
3ERT_master_spheres.sph
    |
    v
Original positive-control docking workflow
    |
    v
Sphere selection -> SHOWBOX -> GRID -> DOCK6

**Running the complete workflow**

The complete script requires both the raw PDB filename and the desired sphere radius.

**For example:**

chmod +x run_3ERT_complete_positive_control.sh
./run_3ERT_complete_positive_control.sh 3ERT.pdb 10

**Software Requirements**

The validated docking stage requires the relevant DOCK6 tools and parameter files, including:

sphere_selector

showbox

grid

dock6

DOCK6 parameter files

The complete start-to-end script additionally requires:

UCSF Chimera

dms

SPHGEN

**Important Reproducibility Note**

The prepared receptor, prepared OHT ligand, master sphere file, and original docking script correspond to the workflow that was actually used for the 3ERT positive-control docking calculations.

The merged run_3ERT_complete_positive_control.sh script incorporates an automated preparation stage so that the workflow can, in principle, begin from the raw 3ERT.pdb structure.

However, the automated preparation portion has not yet been validated end-to-end against the original manually prepared UCSF Chimera files on the current server environment. Therefore, the supplied prepared input files and original docking script should be retained alongside the complete script.

**This provides two reproducibility routes:**

Validated docking reproduction — begin from the supplied prepared receptor, prepared OHT ligand, and master spheres.
Complete workflow reproduction — begin from 3ERT.pdb using the merged preparation-to-docking script in an environment containing all required software.

**Generated Files**
Radius-specific output directories such as:

radius4_full/
radius6_full/
radius8_full/
...
contain generated intermediate and docking output files. These do not need to be stored in GitHub if they can be regenerated using the supplied workflow.

**Typical generated files include:**

selected_spheres_radius<RADIUS>.sph
showbox_radius<RADIUS>.in
showbox_radius<RADIUS>.out
box_radius<RADIUS>.pdb
grid_radius<RADIUS>.in
grid_radius<RADIUS>.out
grid_radius<RADIUS>.nrg
dock_radius<RADIUS>.in
dock_radius<RADIUS>.out
OHT_docked_radius<RADIUS>_scored.mol2
