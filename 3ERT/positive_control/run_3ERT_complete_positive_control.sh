#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================================
# 3ERT / OHT COMPLETE POSITIVE-CONTROL WORKFLOW
#
# Raw 3ERT.pdb
#   -> receptor + OHT extraction
#   -> UCSF Chimera preparation
#   -> molecular surface generation
#   -> SPHGEN master spheres
#   -> original 3ERT positive-control DOCK6 workflow for ONE requested radius
#
# Usage:
#   ./run_3ERT_complete_positive_control.sh 3ERT.pdb 6
#
# Required:
#   UCSF Chimera
#   dms
#   sphgen
#   sphere_selector
#   showbox
#   grid
#   dock6
# ============================================================================

BASE_DIR="$(pwd)"

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 RAW_PDB RADIUS"
    echo "Example: $0 3ERT.pdb 6"
    exit 1
fi

RAW_PDB="$1"
RADIUS="$2"

if [[ ! -s "$RAW_PDB" ]]; then
    echo "ERROR: Raw PDB not found or empty: $RAW_PDB"
    exit 1
fi

if ! [[ "$RADIUS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Radius must be an integer."
    exit 1
fi

if (( RADIUS < 4 || RADIUS > 40 || RADIUS % 2 != 0 )); then
    echo "ERROR: Radius must be one of: 4, 6, 8, ..., 40."
    exit 1
fi

MARGIN=$((RADIUS + 6))
RUN_DIR="$BASE_DIR/radius${RADIUS}_full"

# --------------------------------------------------------------------------
# Software / parameter paths
# --------------------------------------------------------------------------

CHIMERA_BIN="${CHIMERA_BIN:-$(command -v chimera 2>/dev/null || true)}"
DMS_BIN="${DMS_BIN:-$(command -v dms 2>/dev/null || true)}"
PARAM_DIR="${DOCK6_PARAM_DIR:-/home/mdz/app/dock6/parameters}"

for cmd in sphgen sphere_selector showbox grid dock6; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $cmd"
        exit 1
    fi
done

if [[ -z "$CHIMERA_BIN" || ! -x "$CHIMERA_BIN" ]]; then
    echo "ERROR: UCSF Chimera was not found."
    echo "Install Chimera or set:"
    echo "CHIMERA_BIN=/full/path/to/chimera"
    exit 1
fi

if [[ -z "$DMS_BIN" || ! -x "$DMS_BIN" ]]; then
    echo "ERROR: dms was not found."
    echo "Install dms or set:"
    echo "DMS_BIN=/full/path/to/dms"
    exit 1
fi

# --------------------------------------------------------------------------
# File names
# --------------------------------------------------------------------------

RECEPTOR_RAW="$BASE_DIR/3ERT_receptor_raw.pdb"
OHT_NATIVE="$BASE_DIR/OHT_native.pdb"

PREPARED_RECEPTOR="$BASE_DIR/3ERT_receptor_prepared.mol2"
PREPARED_RECEPTOR_PDB="$BASE_DIR/3ERT_receptor_prepared.pdb"
PREPARED_LIGAND="$BASE_DIR/OHT_prepared.mol2"

RECEPTOR_MS="$BASE_DIR/3ERT_receptor.ms"
MASTER_SPHERES="$BASE_DIR/3ERT_master_spheres.sph"

CHIMERA_SCRIPT="$BASE_DIR/prepare_3ERT_chimera.cmd"

# ============================================================================
# STAGE 0A — EXTRACT RECEPTOR AND NATIVE OHT
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0A — EXTRACTING 3ERT RECEPTOR AND OHT"
echo "============================================================"

# Protein receptor: keep ATOM records only.
# This excludes OHT, waters, ions, and other HETATM records.
awk '
/^ATOM  / {print}
END {print "END"}
' "$RAW_PDB" > "$RECEPTOR_RAW"

# Native ligand OHT.
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

echo "Receptor extracted: $RECEPTOR_RAW"
echo "OHT extracted     : $OHT_NATIVE"

# ============================================================================
# STAGE 0B — UCSF CHIMERA PREPARATION
#
# Receptor:
#   - add hydrogens
#   - assign Amber ff14SB charges/types
#
# OHT:
#   - add hydrogens
#   - AM1-BCC charges
#   - formal/net charge = +1
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0B — PREPARING RECEPTOR AND OHT IN UCSF CHIMERA"
echo "============================================================"

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

if [[ ! -s "$PREPARED_RECEPTOR" ]]; then
    echo "ERROR: Chimera did not create $PREPARED_RECEPTOR"
    exit 1
fi

if [[ ! -s "$PREPARED_RECEPTOR_PDB" ]]; then
    echo "ERROR: Chimera did not create $PREPARED_RECEPTOR_PDB"
    exit 1
fi

if [[ ! -s "$PREPARED_LIGAND" ]]; then
    echo "ERROR: Chimera did not create $PREPARED_LIGAND"
    exit 1
fi

echo "Prepared receptor: $PREPARED_RECEPTOR"
echo "Prepared OHT     : $PREPARED_LIGAND"

# ============================================================================
# STAGE 0C — GENERATE MOLECULAR SURFACE
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0C — GENERATING 3ERT RECEPTOR MOLECULAR SURFACE"
echo "============================================================"

"$DMS_BIN" "$RECEPTOR_RAW" -n -w 1.4 -v -o "$RECEPTOR_MS"

if [[ ! -s "$RECEPTOR_MS" ]]; then
    echo "ERROR: dms did not create $RECEPTOR_MS"
    exit 1
fi

echo "Molecular surface created: $RECEPTOR_MS"

# ============================================================================
# STAGE 0D — SPHGEN
#
# These are the exact INSPH settings recovered from the original 3ERT work:
#
# 3ERT_receptor.ms
# R
# X
# 0.0
# 4.0
# 1.4
# 3ERT_master_spheres.sph
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0D — GENERATING 3ERT MASTER SPHERES"
echo "============================================================"

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

rm -f "$MASTER_SPHERES"
sphgen

if [[ ! -s "$MASTER_SPHERES" ]]; then
    echo "ERROR: SPHGEN did not create $MASTER_SPHERES"
    exit 1
fi

echo "Master spheres created: $MASTER_SPHERES"

# ============================================================================
# ORIGINAL 3ERT POSITIVE-CONTROL DOCKING WORKFLOW
# The docking settings below are preserved from the user's original script.
# ============================================================================

for f in "$MASTER_SPHERES" "$PREPARED_LIGAND" "$PREPARED_RECEPTOR"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: Missing or empty file: $f"
        exit 1
    fi
done

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

echo
echo "============================================================"
echo "3ERT POSITIVE CONTROL — RADIUS ${RADIUS} A"
echo "============================================================"

echo "[1/4] Selecting spheres within ${RADIUS} A"

sphere_selector \
    "$MASTER_SPHERES" \
    "$PREPARED_LIGAND" \
    "$RADIUS"

if [[ ! -s selected_spheres.sph ]]; then
    echo "ERROR: sphere_selector failed"
    exit 1
fi

mv selected_spheres.sph "selected_spheres_radius${RADIUS}.sph"

echo "[2/4] Creating SHOWBOX box with ${MARGIN} A extra margin"

cat > "showbox_radius${RADIUS}.in" <<EOF
Y
${MARGIN}
selected_spheres_radius${RADIUS}.sph
1
box_radius${RADIUS}.pdb
EOF

showbox < "showbox_radius${RADIUS}.in" > "showbox_radius${RADIUS}.out" 2>&1

if [[ ! -s box_radius${RADIUS}.pdb ]]; then
    echo "ERROR: SHOWBOX failed"
    cat showbox_radius${RADIUS}.out
    exit 1
fi

echo "[3/4] Running GRID"

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
receptor_file                  $PREPARED_RECEPTOR
box_file                       box_radius${RADIUS}.pdb
vdw_definition_file            $PARAM_DIR/vdw_AMBER_parm99.defn
score_grid_prefix              grid_radius${RADIUS}
EOF

grid -i "grid_radius${RADIUS}.in" -o "grid_radius${RADIUS}.out"

if [[ ! -s "grid_radius${RADIUS}.nrg" ]]; then
    echo "ERROR: GRID failed"
    tail -40 "grid_radius${RADIUS}.out" || true
    exit 1
fi

echo "[4/4] Running DOCK6"

cat > "dock_radius${RADIUS}.in" <<EOF
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                             no
min_anchor_size                                               5
pruning_use_clustering                                        yes
pruning_max_orients                                           1000
pruning_clustering_cutoff                                     100
use_clash_overlap                                             no
write_growth_trees                                            no
ligand_atom_file                                              $PREPARED_LIGAND
limit_max_ligands                                             no
skip_molecule                                                 no
read_mol_solvation                                            no
calculate_rmsd                                                yes
use_rmsd_reference_mol                                        yes
rmsd_reference_filename                                       $PREPARED_LIGAND
use_database_filter                                           no
orient_ligand                                                 yes
automated_matching                                            yes
receptor_site_file                                            selected_spheres_radius${RADIUS}.sph
max_orientations                                              5000
critical_points                                               no
chemical_matching                                             no
use_ligand_spheres                                            no
use_internal_energy                                           yes
internal_energy_rep_exp                                       12
flexible_ligand                                               yes
bump_filter                                                   no
score_molecules                                               yes
contact_score_primary                                         no
contact_score_secondary                                       no
grid_score_primary                                            yes
grid_score_secondary                                          no
grid_score_rep_rad_scale                                      1
grid_score_vdw_scale                                          1
grid_score_es_scale                                           1
grid_score_grid_prefix                                        grid_radius${RADIUS}
multigrid_score_secondary                                     no
dock3.5_score_secondary                                       no
continuous_score_secondary                                    no
footprint_similarity_score_secondary                          no
pharmacophore_score_secondary                                 no
descriptor_score_secondary                                    no
gbsa_zou_score_secondary                                      no
gbsa_hawkins_score_secondary                                  no
SASA_score_secondary                                          no
amber_score_secondary                                         no
minimize_ligand                                               yes
simplex_max_iterations                                        1000
simplex_tors_premin_iterations                                0
simplex_max_cycles                                            1
simplex_score_converge                                        0.1
simplex_cycle_converge                                        1
simplex_trans_step                                            1
simplex_rot_step                                              0.1
simplex_tors_step                                             10
simplex_anchor_max_iterations                                 500
simplex_grow_max_iterations                                   250
simplex_grow_tors_premin_iterations                           0
simplex_final_min                                             no
simplex_random_seed                                           0
simplex_restraint_min                                         no
atom_model                                                    all
vdw_defn_file                                                 $PARAM_DIR/vdw_AMBER_parm99.defn
flex_defn_file                                                $PARAM_DIR/flex.defn
flex_drive_file                                               $PARAM_DIR/flex_drive.tbl
ligand_outfile_prefix                                         OHT_docked_radius${RADIUS}
write_orientations                                            no
num_scored_conformers                                         1
rank_ligands                                                  no
EOF

dock6 -i "dock_radius${RADIUS}.in" -o "dock_radius${RADIUS}.out"

if [[ ! -s "dock_radius${RADIUS}.out" ]]; then
    echo "ERROR: DOCK6 produced no output"
    exit 1
fi

echo
echo "============================================================"
echo "RADIUS ${RADIUS} FULL WORKFLOW COMPLETE"
echo "============================================================"
echo "Folder: $RUN_DIR"
echo "Sphere file: $RUN_DIR/selected_spheres_radius${RADIUS}.sph"
echo "Box file: $RUN_DIR/box_radius${RADIUS}.pdb"
echo "Grid output: $RUN_DIR/grid_radius${RADIUS}.out"
echo "Dock output: $RUN_DIR/dock_radius${RADIUS}.out"
echo "Docked ligand: $RUN_DIR/OHT_docked_radius${RADIUS}_scored.mol2"
echo "============================================================"
