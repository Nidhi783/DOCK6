#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# 1J3K POSITIVE CONTROL — COMPLETE MERGED WORKFLOW
#
# Raw 1J3K.pdb
#   -> extract chain A receptor and native WRA (residue 609)
#   -> UCSF Chimera preparation
#   -> receptor molecular surface generation
#   -> SPHGEN master spheres
#   -> original 1J3K positive-control docking workflow
#
# WRA preparation:
#   AM1-BCC charges
#   formal/net charge = 0
#
# Usage:
#   ./run_1J3K_complete_positive_control.sh 1J3K.pdb RADIUS
#
# Example:
#   ./run_1J3K_complete_positive_control.sh 1J3K.pdb 10
#
# Radius: even values from 4 to 28 Å
# SHOWBOX margin = radius + 6 Å
# ============================================================

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 RAW_PDB RADIUS"
    echo "Example: $0 1J3K.pdb 10"
    exit 1
fi

RAW_PDB="$1"
RADIUS="$2"
BASE_DIR="$(pwd)"

if [[ ! -s "$RAW_PDB" ]]; then
    echo "ERROR: Raw PDB not found or empty: $RAW_PDB"
    exit 1
fi

if ! [[ "$RADIUS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Radius must be an integer."
    exit 1
fi

if (( RADIUS < 4 || RADIUS > 28 || RADIUS % 2 != 0 )); then
    echo "ERROR: Radius must be an even integer from 4 to 28."
    exit 1
fi

CHIMERA_BIN="${CHIMERA_BIN:-$(command -v chimera 2>/dev/null || true)}"
DMS_BIN="${DMS_BIN:-$(command -v dms 2>/dev/null || true)}"

for CMD in sphgen sphere_selector showbox grid dock6; do
    if ! command -v "$CMD" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $CMD"
        exit 1
    fi
done

if [[ -z "$CHIMERA_BIN" || ! -x "$CHIMERA_BIN" ]]; then
    echo "ERROR: UCSF Chimera was not found."
    echo "Set CHIMERA_BIN=/full/path/to/chimera if needed."
    exit 1
fi

if [[ -z "$DMS_BIN" || ! -x "$DMS_BIN" ]]; then
    echo "ERROR: dms was not found."
    echo "Set DMS_BIN=/full/path/to/dms if needed."
    exit 1
fi

RECEPTOR_RAW="$BASE_DIR/1J3K_chainA_receptor_WRAremoved.pdb"
WRA_NATIVE_PDB="$BASE_DIR/WRA_only.pdb"
PREPARED_RECEPTOR="$BASE_DIR/1J3K_chainA_receptor_WRAremoved_charged.mol2"
PREPARED_RECEPTOR_PDB="$BASE_DIR/1J3K_chainA_receptor_WRAremoved_charged.pdb"
PREPARED_LIGAND="$BASE_DIR/WRA_only.mol2"
RECEPTOR_DMS="$BASE_DIR/1J3K_chainA_receptor_WRAremoved.dms"
MASTER_SPHERES_PREP="$BASE_DIR/1J3K_master_spheres.sph"
CHIMERA_SCRIPT="$BASE_DIR/prepare_1J3K_chimera.cmd"

# ============================================================
# STAGE 0A — EXTRACT CHAIN A RECEPTOR AND NATIVE WRA
# ============================================================

echo
echo "============================================"
echo "STAGE 0A — EXTRACTING 1J3K CHAIN A + WRA"
echo "============================================"

# Chain A protein ATOM records only; excludes WRA and other HETATM records.
awk '
/^ATOM  / && substr($0,22,1)=="A" {print}
END {print "END"}
' "$RAW_PDB" > "$RECEPTOR_RAW"

# Native WRA from chain A, residue 609.
awk '
/^HETATM/ && substr($0,18,3)=="WRA" && substr($0,22,1)=="A" && (substr($0,23,4)+0)==609 {print}
END {print "END"}
' "$RAW_PDB" > "$WRA_NATIVE_PDB"

if [[ $(grep -c '^ATOM  ' "$RECEPTOR_RAW" || true) -eq 0 ]]; then
    echo "ERROR: No chain-A receptor atoms were extracted."
    exit 1
fi

if [[ $(grep -c '^HETATM' "$WRA_NATIVE_PDB" || true) -eq 0 ]]; then
    echo "ERROR: WRA chain A residue 609 was not found in $RAW_PDB."
    exit 1
fi

# ============================================================
# STAGE 0B — UCSF CHIMERA PREPARATION
# ============================================================

echo
echo "============================================"
echo "STAGE 0B — PREPARING RECEPTOR AND WRA"
echo "============================================"

cat > "$CHIMERA_SCRIPT" <<EOF
open noprefs $RECEPTOR_RAW
addh
addcharge std chargeModel ff14SB
write format mol2 atomTypes amber #0 $PREPARED_RECEPTOR
write format pdb #0 $PREPARED_RECEPTOR_PDB
close all

open noprefs $WRA_NATIVE_PDB
addh
addcharge nonstd :WRA 0 method am1
write format mol2 atomTypes amber #0 $PREPARED_LIGAND
close all
stop
EOF

"$CHIMERA_BIN" --nogui "$CHIMERA_SCRIPT"

for FILE in "$PREPARED_RECEPTOR" "$PREPARED_RECEPTOR_PDB" "$PREPARED_LIGAND"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Preparation failed; missing or empty file: $FILE"
        exit 1
    fi
done

# ============================================================
# STAGE 0C — GENERATE RECEPTOR DMS SURFACE
# ============================================================

echo
echo "============================================"
echo "STAGE 0C — GENERATING RECEPTOR DMS"
echo "============================================"

"$DMS_BIN" "$RECEPTOR_RAW" -n -w 1.4 -v -o "$RECEPTOR_DMS"

if [[ ! -s "$RECEPTOR_DMS" ]]; then
    echo "ERROR: dms did not create $RECEPTOR_DMS"
    exit 1
fi

# ============================================================
# STAGE 0D — SPHGEN MASTER SPHERES
# Exact recovered INSPH:
#
# 1J3K_chainA_receptor_WRAremoved.dms
# R
# X
# 0.0
# 4.0
# 1.4
# 1J3K_master_spheres.sph
# ============================================================

echo
echo "============================================"
echo "STAGE 0D — GENERATING MASTER SPHERES"
echo "============================================"

cd "$BASE_DIR"

cat > INSPH <<EOF
1J3K_chainA_receptor_WRAremoved.dms
R
X
0.0
4.0
1.4
1J3K_master_spheres.sph
EOF

rm -f "$MASTER_SPHERES_PREP"
sphgen

if [[ ! -s "$MASTER_SPHERES_PREP" ]]; then
    echo "ERROR: SPHGEN did not create $MASTER_SPHERES_PREP"
    exit 1
fi

# ============================================================
# ORIGINAL 1J3K POSITIVE-CONTROL DOCKING WORKFLOW
# The original sphere-selection, SHOWBOX, GRID and DOCK6 logic
# is preserved below.
# ============================================================

# ============================================================
# MARGIN
# ============================================================

MARGIN=$((RADIUS + 6))


# ============================================================
# PATHS
# ============================================================

MASTER_SPHERES="$BASE_DIR/1J3K_master_spheres.sph"
LIGAND="$BASE_DIR/WRA_only.mol2"
RECEPTOR="$BASE_DIR/1J3K_chainA_receptor_WRAremoved_charged.mol2"

PARAM_DIR="${DOCK6_PARAM_DIR:-/home/mdz/app/dock6/parameters}"

VDW_FILE="$PARAM_DIR/vdw_AMBER_parm99.defn"
FLEX_FILE="$PARAM_DIR/flex.defn"
FLEX_DRIVE_FILE="$PARAM_DIR/flex_drive.tbl"

RUN_DIR="$BASE_DIR/radius${RADIUS}_full"


# ============================================================
# CHECK REQUIRED FILES
# ============================================================

for FILE in \
    "$MASTER_SPHERES" \
    "$LIGAND" \
    "$RECEPTOR" \
    "$VDW_FILE" \
    "$FLEX_FILE" \
    "$FLEX_DRIVE_FILE"
do
    if [[ ! -s "$FILE" ]]; then
        echo
        echo "ERROR: Missing or empty file:"
        echo "$FILE"
        exit 1
    fi
done


# ============================================================
# CREATE CLEAN RUN DIRECTORY
# ============================================================

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"

cd "$RUN_DIR"

cp "$RECEPTOR" receptor.mol2
cp "$LIGAND" ligand.mol2
cp "$MASTER_SPHERES" master_spheres.sph


echo
echo "============================================"
echo "1J3K POSITIVE CONTROL"
echo "Sphere radius  : ${RADIUS} A"
echo "SHOWBOX margin : ${MARGIN} A"
echo "============================================"


# ============================================================
# STEP 1 — SPHERE SELECTION
# ============================================================

echo
echo "[1/4] Selecting spheres..."

rm -f selected_spheres.sph

sphere_selector \
    master_spheres.sph \
    ligand.mol2 \
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
    echo
    echo "ERROR: SHOWBOX failed."
    echo
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
    echo
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

ligand_outfile_prefix                                        WRA_docked_radius${RADIUS}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

dock6 \
    -i "dock_radius${RADIUS}.in" \
    -o "dock_radius${RADIUS}.out"


# ============================================================
# VERIFY DOCKING OUTPUT
# ============================================================

SCORED_FILE="WRA_docked_radius${RADIUS}_scored.mol2"

if [[ ! -s "$SCORED_FILE" ]]; then
    echo
    echo "ERROR: DOCK6 did not create:"
    echo "$SCORED_FILE"
    echo
    echo "Last 50 lines of DOCK6 output:"
    tail -50 "dock_radius${RADIUS}.out" || true
    exit 1
fi


# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================"
echo "RADIUS ${RADIUS} A COMPLETE"
echo "============================================"

echo
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

echo
echo "============================================"