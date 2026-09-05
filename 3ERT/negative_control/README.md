**3ERT Negative Control – DOCK6**

**Overview**

This folder contains the negative-control docking workflow for the 3ERT–OHT protein–ligand system using DOCK6. The negative control tests OHT docking at a site deliberately positioned away from the crystallographic OHT-binding site.

**Negative-Control Site**
**Parameter                      Value**

OHT heavy atoms                 29

Crystal OHT centroid           [31.571, -1.596, 25.599] Å

Negative-control centre        [13.035, 13.011, 18.731] Å

Crystal-to-negative distance   24.579 Å

The negative-control centre is approximately 24.6 Å from the crystallographic OHT centroid. 3ERT_negative_site_reference.mol2 represents this location and is used by sphere_selector to select spheres around the negative-control site rather than the native OHT-binding site.

**Main Input Files**

3ERT.pdb – original 3ERT structure
OHT_prepared.mol2 – prepared native OHT ligand
3ERT_receptor_prepared.mol2 – prepared receptor
3ERT_master_spheres.sph – master sphere set
3ERT_negative_site_reference.mol2 – negative-site reference
run_3ERT_negative_control.sh – original negative-control docking script
run_3ERT_complete_negative_control.sh – merged preparation-to-docking workflow

OHT was prepared using AM1-BCC charges with a net charge of +1.

**Workflow**

3ERT receptor
      ↓
Receptor preparation
      ↓
Molecular surface → SPHGEN
      ↓
3ERT_master_spheres.sph
      +
3ERT_negative_site_reference.mol2
      ↓
sphere_selector
      ↓
Negative-site spheres
      ↓
SHOWBOX → GRID → DOCK6
      ↓
Docked OHT pose

**The recovered SPHGEN input was**:

3ERT_receptor.ms
R
X
0.0
4.0
1.4
3ERT_master_spheres.sph

**Docking Parameters**

The original negative-control script supports even sphere-selection radii from 4 Å to 30 Å.

Sphere radii: 4, 6, 8, ... 30 Å

SHOWBOX margin: sphere radius + 6 Å

GRID spacing: 0.5 Å

Flexible ligand docking: enabled

Maximum orientations: 5000

Primary scoring: GRID score

Native OHT: RMSD reference

Each radius is run separately and creates its own output directory, such as radius10_full/.

Running the Workflow

**Original negative-control docking script:**

chmod +x run_3ERT_negative_control.sh
./run_3ERT_negative_control.sh 10

**Merged workflow beginning from the original PDB:**

chmod +x run_3ERT_complete_negative_control.sh
./run_3ERT_complete_negative_control.sh 3ERT.pdb 10

Replace 10 with the required even sphere radius between 4 Å and 30 Å.

**Output**

Each radius-specific directory contains the selected spheres, docking box, GRID files, DOCK6 output and scored OHT pose. Principal files include:

selected_spheres_radius<RADIUS>.sph
box_radius<RADIUS>.pdb
grid_radius<RADIUS>.out
dock_radius<RADIUS>.out
OHT_negative_docked_radius<RADIUS>_scored.mol2

**Reproducibility**

The prepared receptor, prepared OHT ligand, master spheres and 3ERT_negative_site_reference.mol2 are retained to reproduce the negative-control docking workflow. The negative-site reference is preserved because it defines the exact negative-control location used in this project.

The merged script additionally includes receptor/ligand preparation, molecular-surface generation and SPHGEN before continuing with the original negative-control docking procedure.
