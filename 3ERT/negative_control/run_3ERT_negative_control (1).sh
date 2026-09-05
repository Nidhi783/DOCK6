#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 3ERT NEGATIVE CONTROL
# Increasing sphere-selection radius
#
# Sphere radius: even values from 4 to 30 Å
# SHOWBOX margin = sphere radius + 6 Å
#
# Examples:
#   ./run_3ERT_negative_control.sh 22
#   ./run_3ERT_negative_control.sh 24
#   ./run_3ERT_negative_control.sh 26
#   ./run_3ERT_negative_control.sh 28
#   ./run_3ERT_negative_control.sh 30
# ============================================================


# ============================================================
# INPUT CHECK
# ============================================================

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 RADIUS"
    echo "Example: $0 22"
    exit 1
fi

RADIUS="$1"

if ! [[ "$RADIUS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Radius must be an integer."
    exit 1
fi

if (( RADIUS < 4 || RADIUS > 30 || RADIUS % 2 != 0 )); then
    echo "ERROR: Radius must be an even integer from 4 to 30."
    exit 1
fi


# ============================================================
# PATHS
# ============================================================

MARGIN=$((RADIUS + 6))
BASE_DIR="$(pwd)"

MASTER_SPHERES="$BASE_DIR/3ERT_master_spheres.sph"
NEGATIVE_REF="$BASE_DIR/3ERT_negative_site_reference.mol2"
LIGAND="$BASE_DIR/OHT_prepared.mol2"
RECEPTOR="$BASE_DIR/3ERT_receptor_prepared.mol2"

PARAM_DIR="/home/mdz/app/dock6/parameters"

VDW_FILE="$PARAM_DIR/vdw_AMBER_parm99.defn"
FLEX_FILE="$PARAM_DIR/flex.defn"
FLEX_DRIVE_FILE="$PARAM_DIR/flex_drive.tbl"

RUN_DIR="$BASE_DIR/radius${RADIUS}_full"


# ============================================================
# CHECK REQUIRED FILES
# ============================================================

for FILE in \
    "$MASTER_SPHERES" \
    "$NEGATIVE_REF" \
    "$LIGAND" \
    "$RECEPTOR" \
    "$VDW_FILE" \
    "$FLEX_FILE" \
    "$FLEX_DRIVE_FILE"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty file:"
        echo "$FILE"
        exit 1
    fi
done


# ============================================================
# CREATE RUN DIRECTORY
# ============================================================

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

cp "$MASTER_SPHERES" master_spheres.sph
cp "$NEGATIVE_REF" negative_site_reference.mol2
cp "$LIGAND" ligand.mol2
cp "$RECEPTOR" receptor.mol2


echo
echo "============================================"
echo "3ERT NEGATIVE CONTROL"
echo "Sphere radius  : ${RADIUS} A"
echo "SHOWBOX margin : ${MARGIN} A"
echo "============================================"


# ============================================================
# STEP 1 — SELECT NEGATIVE-SITE SPHERES
# ============================================================

echo
echo "[1/4] Selecting negative-site spheres..."

rm -f selected_spheres.sph

sphere_selector \
    master_spheres.sph \
    negative_site_reference.mol2 \
    "$RADIUS"

if [[ ! -s selected_spheres.sph ]]; then
    echo "ERROR: sphere_selector failed."
    exit 1
fi

mv selected_spheres.sph \
   "selected_spheres_radius${RADIUS}.sph"


NUM_SPHERES=$(
    awk '/number of spheres in cluster/ {print $NF; exit}' \
    "selected_spheres_radius${RADIUS}.sph"
)

if [[ -z "${NUM_SPHERES:-}" || "${NUM_SPHERES}" == "0" ]]; then
    echo "ERROR: No spheres were selected."
    exit 1
fi

echo "Selected spheres: ${NUM_SPHERES}"


# ============================================================
# STEP 2 — SHOWBOX
# ============================================================

echo
echo "[2/4] Creating SHOWBOX..."

cat > "showbox_radius${RADIUS}.in" <<EOF
Y
${MARGIN}
selected_spheres_radius${RADIUS}.sph
1
box_radius${RADIUS}.pdb
EOF

showbox \
    < "showbox_radius${RADIUS}.in" \
    > "showbox_radius${RADIUS}.out" 2>&1


if [[ ! -s "box_radius${RADIUS}.pdb" ]]; then
    echo "ERROR: SHOWBOX failed."
    cat "showbox_radius${RADIUS}.out"
    exit 1
fi

echo "SHOWBOX complete."


# ============================================================
# STEP 3 — GRID
# ============================================================

echo
echo "[3/4] Running GRID..."

cat > "grid_radius${RADIUS}.in" <<EOF
compute_grids                  yes
grid_spacing                   0.5
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     all
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
allow_non_integral_charges     yes
bump_filter                    no
receptor_file                  receptor.mol2
box_file                       box_radius${RADIUS}.pdb
vdw_definition_file            ${VDW_FILE}
score_grid_prefix              grid_radius${RADIUS}
EOF


grid \
    -i "grid_radius${RADIUS}.in" \
    -o "grid_radius${RADIUS}.out"


if [[ ! -s "grid_radius${RADIUS}.nrg" ]]; then
    echo
    echo "ERROR: GRID failed."
    echo "Last 40 lines of GRID output:"
    tail -40 "grid_radius${RADIUS}.out" || true
    exit 1
fi

echo "GRID complete."


# ============================================================
# STEP 4 — DOCK6
# ============================================================

echo
echo "[4/4] Running DOCK6..."

cat > "dock_radius${RADIUS}.in" <<EOF
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
use_clash_overlap                                            no
write_growth_trees                                           no

ligand_atom_file                                             ligand.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no

calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      ligand.mol2

use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           selected_spheres_radius${RADIUS}.sph

max_orientations                                             5000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no

use_internal_energy                                          yes
internal_energy_rep_exp                                      12
flexible_ligand                                              yes
bump_filter                                                  no

score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no

grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       grid_radius${RADIUS}

multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no

minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1
simplex_trans_step                                           1
simplex_rot_step                                             0.1
simplex_tors_step                                            10
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  250
simplex_grow_tors_premin_iterations                          0
simplex_final_min                                            no
simplex_random_seed                                          0
simplex_restraint_min                                        no

atom_model                                                   all
vdw_defn_file                                                ${VDW_FILE}
flex_defn_file                                               ${FLEX_FILE}
flex_drive_file                                              ${FLEX_DRIVE_FILE}

ligand_outfile_prefix                                        OHT_negative_docked_radius${RADIUS}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF


dock6 \
    -i "dock_radius${RADIUS}.in" \
    -o "dock_radius${RADIUS}.out"


# ============================================================
# CHECK DOCKING OUTPUT
# ============================================================

SCORED_FILE="OHT_negative_docked_radius${RADIUS}_scored.mol2"

if [[ ! -s "$SCORED_FILE" ]]; then
    echo
    echo "ERROR: DOCK6 did not create $SCORED_FILE"
    echo "Last 50 lines of DOCK6 output:"
    tail -50 "dock_radius${RADIUS}.out" || true
    exit 1
fi


# ============================================================
# FINAL RESULTS
# ============================================================

echo
echo "============================================"
echo "RADIUS ${RADIUS} A COMPLETE"
echo "============================================"

echo "Sphere radius    : ${RADIUS} A"
echo "SHOWBOX margin   : ${MARGIN} A"
echo "Selected spheres : ${NUM_SPHERES}"

echo
echo "Grid score:"
grep -i "Grid_Score" "$SCORED_FILE" | head -1 || true

echo
echo "Docked ligand:"
echo "$RUN_DIR/$SCORED_FILE"

echo
echo "Run directory:"
echo "$RUN_DIR"

echo "============================================"