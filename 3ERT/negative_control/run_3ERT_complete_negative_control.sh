#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# 3ERT NEGATIVE CONTROL — COMPLETE MERGED WORKFLOW
#
# Usage:
#   ./run_3ERT_complete_negative_control.sh 3ERT.pdb RADIUS
#
# Example:
#   ./run_3ERT_complete_negative_control.sh 3ERT.pdb 10
#
# Radius: even values from 4 to 30 Å
# SHOWBOX margin = radius + 6 Å
#
# NOTE:
# 3ERT_negative_site_reference.mol2 is retained as the validated
# negative-site definition and is NOT regenerated automatically.
# ============================================================

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 RAW_PDB RADIUS"
    echo "Example: $0 3ERT.pdb 10"
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

if (( RADIUS < 4 || RADIUS > 30 || RADIUS % 2 != 0 )); then
    echo "ERROR: Radius must be an even integer from 4 to 30."
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

RECEPTOR_RAW="$BASE_DIR/3ERT_receptor_raw.pdb"
OHT_NATIVE="$BASE_DIR/OHT_native.pdb"

PREPARED_RECEPTOR="$BASE_DIR/3ERT_receptor_prepared.mol2"
PREPARED_RECEPTOR_PDB="$BASE_DIR/3ERT_receptor_prepared.pdb"
PREPARED_LIGAND="$BASE_DIR/OHT_prepared.mol2"

RECEPTOR_MS="$BASE_DIR/3ERT_receptor.ms"
MASTER_SPHERES_PREP="$BASE_DIR/3ERT_master_spheres.sph"
NEGATIVE_REF_PREP="$BASE_DIR/3ERT_negative_site_reference.mol2"
CHIMERA_SCRIPT="$BASE_DIR/prepare_3ERT_negative_control_chimera.cmd"

echo
echo "============================================"
echo "STAGE 0A — EXTRACTING 3ERT RECEPTOR AND OHT"
echo "============================================"

awk '
/^ATOM  / {print}
END {print "END"}
' "$RAW_PDB" > "$RECEPTOR_RAW"

awk '
/^HETATM/ && substr($0,18,3)=="OHT" {print}
END {print "END"}
' "$RAW_PDB" > "$OHT_NATIVE"

if [[ $(grep -c '^ATOM  ' "$RECEPTOR_RAW" || true) -eq 0 ]]; then
    echo "ERROR: No receptor ATOM records were extracted."
    exit 1
fi

if [[ $(grep -c '^HETATM' "$OHT_NATIVE" || true) -eq 0 ]]; then
    echo "ERROR: Native ligand OHT was not found in $RAW_PDB."
    exit 1
fi

echo
echo "============================================"
echo "STAGE 0B — PREPARING RECEPTOR AND OHT"
echo "============================================"

cat > "$CHIMERA_SCRIPT" <<EOF
open noprefs $RECEPTOR_RAW
addh
addcharge std chargeModel ff14SB
write format mol2 atomTypes amber #0 $PREPARED_RECEPTOR
write format pdb #0 $PREPARED_RECEPTOR_PDB
close all

open noprefs $OHT_NATIVE
addh
addcharge nonstd :OHT 1 method am1
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

echo
echo "============================================"
echo "STAGE 0C — GENERATING MOLECULAR SURFACE"
echo "============================================"

"$DMS_BIN" "$RECEPTOR_RAW" -n -w 1.4 -v -o "$RECEPTOR_MS"

if [[ ! -s "$RECEPTOR_MS" ]]; then
    echo "ERROR: dms did not create $RECEPTOR_MS"
    exit 1
fi

echo
echo "============================================"
echo "STAGE 0D — GENERATING MASTER SPHERES"
echo "============================================"

cd "$BASE_DIR"

cat > INSPH <<EOF
3ERT_receptor.ms
R
X
0.0
4.0
1.4
3ERT_master_spheres.sph
EOF

rm -f "$MASTER_SPHERES_PREP"
sphgen

if [[ ! -s "$MASTER_SPHERES_PREP" ]]; then
    echo "ERROR: SPHGEN did not create $MASTER_SPHERES_PREP"
    exit 1
fi

echo
echo "============================================"
echo "STAGE 0E — CHECKING NEGATIVE-SITE REFERENCE"
echo "============================================"

if [[ ! -s "$NEGATIVE_REF_PREP" ]]; then
    echo "ERROR: Missing $NEGATIVE_REF_PREP"
    echo "Keep the validated negative-site reference in this directory."
    exit 1
fi

echo "Negative-site reference: $NEGATIVE_REF_PREP"

# ============================================================
# ORIGINAL 3ERT NEGATIVE-CONTROL DOCKING WORKFLOW
# The original docking logic is preserved below.
# ============================================================

# ============================================================
# PATHS
# ============================================================

MARGIN=$((RADIUS + 6))

MASTER_SPHERES="$BASE_DIR/3ERT_master_spheres.sph"
NEGATIVE_REF="$BASE_DIR/3ERT_negative_site_reference.mol2"
LIGAND="$BASE_DIR/OHT_prepared.mol2"
RECEPTOR="$BASE_DIR/3ERT_receptor_prepared.mol2"

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