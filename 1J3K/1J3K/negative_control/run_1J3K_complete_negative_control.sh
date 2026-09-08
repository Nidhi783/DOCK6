#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# 1J3K NEGATIVE CONTROL — COMPLETE START-TO-END WORKFLOW
#
# Workflow:
#   1J3K.pdb
#     -> extract chain A receptor + native WRA
#     -> UCSF Chimera preparation
#     -> receptor DMS surface
#     -> SPHGEN master spheres
#     -> use validated negative-site reference
#     -> sphere selection
#     -> SHOWBOX
#     -> GRID
#     -> DOCK6
#     -> Grid Score + Heavy-atom RMSD + Centroid Distance
#
# Fixed negative-control centre:
#   X = 39.053
#   Y = 17.088
#   Z = 2.043
#
# Radii:
#   4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28 Å
#
# SHOWBOX margin = radius + 6 Å
#
# Usage:
#   ./run_1J3K_complete_negative_control.sh 1J3K.pdb
#
# Notes:
# - WRA is prepared with AM1-BCC charges and net charge 0.
# - 1J3K_negative_site_reference.mol2 is preserved as the
#   validated negative-site definition and is not regenerated.
# ============================================================

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 RAW_PDB"
    echo "Example: $0 1J3K.pdb"
    exit 1
fi

RAW_PDB="$1"
BASE_DIR="$(pwd)"

if [[ ! -s "$RAW_PDB" ]]; then
    echo "ERROR: Raw PDB not found or empty: $RAW_PDB"
    exit 1
fi

RADII=(4 6 8 10 12 14 16 18 20 22 24 26 28)

CHIMERA_BIN="${CHIMERA_BIN:-$(command -v chimera 2>/dev/null || true)}"
DMS_BIN="${DMS_BIN:-$(command -v dms 2>/dev/null || true)}"
PARAM_DIR="${DOCK6_PARAM_DIR:-/home/mdz/app/dock6/parameters}"

for CMD in sphgen sphere_selector showbox grid dock6 python3; do
    if ! command -v "$CMD" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $CMD"
        exit 1
    fi
done

if [[ -z "$CHIMERA_BIN" || ! -x "$CHIMERA_BIN" ]]; then
    echo "ERROR: UCSF Chimera not found."
    echo "Set CHIMERA_BIN=/full/path/to/chimera if required."
    exit 1
fi

if [[ -z "$DMS_BIN" || ! -x "$DMS_BIN" ]]; then
    echo "ERROR: dms not found."
    echo "Set DMS_BIN=/full/path/to/dms if required."
    exit 1
fi

VDW_FILE="$PARAM_DIR/vdw_AMBER_parm99.defn"
FLEX_FILE="$PARAM_DIR/flex.defn"
FLEX_DRIVE_FILE="$PARAM_DIR/flex_drive.tbl"

for FILE in "$VDW_FILE" "$FLEX_FILE" "$FLEX_DRIVE_FILE"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing DOCK6 parameter file: $FILE"
        exit 1
    fi
done

# ============================================================
# PREPARATION FILES
# ============================================================

RECEPTOR_RAW="$BASE_DIR/1J3K_chainA_receptor_WRAremoved.pdb"
WRA_NATIVE="$BASE_DIR/WRA_native.pdb"

RECEPTOR="$BASE_DIR/1J3K_chainA_receptor_WRAremoved_charged.mol2"
RECEPTOR_PREP_PDB="$BASE_DIR/1J3K_chainA_receptor_WRAremoved_prepared.pdb"
LIGAND="$BASE_DIR/WRA_only.mol2"

RECEPTOR_DMS="$BASE_DIR/1J3K_chainA_receptor_WRAremoved.dms"
MASTER_SPHERES="$BASE_DIR/1J3K_master_spheres.sph"
NEGATIVE_REF="$BASE_DIR/1J3K_negative_site_reference.mol2"

CHIMERA_SCRIPT="$BASE_DIR/prepare_1J3K_negative_control_chimera.cmd"
RESULTS="$BASE_DIR/1J3K_negative_control_results.tsv"

# ============================================================
# STAGE 0A — EXTRACT CHAIN A RECEPTOR AND WRA
# ============================================================

echo
echo "===================================================="
echo "STAGE 0A — EXTRACTING CHAIN A RECEPTOR AND WRA"
echo "===================================================="

awk '
/^ATOM  / && substr($0,22,1)=="A" {print}
END {print "END"}
' "$RAW_PDB" > "$RECEPTOR_RAW"

awk '
/^HETATM/ && substr($0,18,3)=="WRA" && substr($0,22,1)=="A" {print}
END {print "END"}
' "$RAW_PDB" > "$WRA_NATIVE"

if [[ $(grep -c '^ATOM  ' "$RECEPTOR_RAW" || true) -eq 0 ]]; then
    echo "ERROR: No chain A receptor ATOM records were extracted."
    exit 1
fi

if [[ $(grep -c '^HETATM' "$WRA_NATIVE" || true) -eq 0 ]]; then
    echo "ERROR: WRA was not found in chain A of $RAW_PDB."
    exit 1
fi

# ============================================================
# STAGE 0B — UCSF CHIMERA PREPARATION
# ============================================================

echo
echo "===================================================="
echo "STAGE 0B — PREPARING RECEPTOR AND WRA"
echo "===================================================="

cat > "$CHIMERA_SCRIPT" <<EOF
open noprefs $RECEPTOR_RAW
addh
addcharge std chargeModel ff14SB
write format mol2 atomTypes amber #0 $RECEPTOR
write format pdb #0 $RECEPTOR_PREP_PDB
close all

open noprefs $WRA_NATIVE
addh
addcharge nonstd :WRA 0 method am1
write format mol2 atomTypes amber #0 $LIGAND
close all
stop
EOF

"$CHIMERA_BIN" --nogui "$CHIMERA_SCRIPT"

for FILE in "$RECEPTOR" "$RECEPTOR_PREP_PDB" "$LIGAND"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Preparation failed; missing or empty file: $FILE"
        exit 1
    fi
done

# ============================================================
# STAGE 0C — GENERATE DMS SURFACE
# ============================================================

echo
echo "===================================================="
echo "STAGE 0C — GENERATING DMS SURFACE"
echo "===================================================="

"$DMS_BIN" "$RECEPTOR_RAW" -n -w 1.4 -v -o "$RECEPTOR_DMS"

if [[ ! -s "$RECEPTOR_DMS" ]]; then
    echo "ERROR: dms did not create $RECEPTOR_DMS"
    exit 1
fi

# ============================================================
# STAGE 0D — SPHGEN MASTER SPHERES
# ============================================================

echo
echo "===================================================="
echo "STAGE 0D — GENERATING MASTER SPHERES"
echo "===================================================="

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

rm -f "$MASTER_SPHERES"
sphgen

if [[ ! -s "$MASTER_SPHERES" ]]; then
    echo "ERROR: SPHGEN did not create $MASTER_SPHERES"
    exit 1
fi

# ============================================================
# STAGE 0E — CHECK NEGATIVE-SITE REFERENCE
# ============================================================

echo
echo "===================================================="
echo "STAGE 0E — CHECKING NEGATIVE-SITE REFERENCE"
echo "===================================================="

if [[ ! -s "$NEGATIVE_REF" ]]; then
    echo "ERROR: Missing validated negative-site reference:"
    echo "$NEGATIVE_REF"
    echo
    echo "Expected centre:"
    echo "X = 39.053"
    echo "Y = 17.088"
    echo "Z = 2.043"
    exit 1
fi

echo "Negative-site centre: [39.053, 17.088, 2.043] Å"

# ============================================================
# RESULTS TABLE
# ============================================================

printf "Radius\tMargin\tSelectedSpheres\tGridScore\tHeavyAtomRMSD\tCentroidDistance\tStatus\n" > "$RESULTS"

# ============================================================
# RMSD + CENTROID FUNCTION
# ============================================================

calculate_metrics () {
python3 - "$1" "$2" <<'PY'
import sys
import numpy as np

reference_file = sys.argv[1]
docked_file = sys.argv[2]

def read_heavy_atoms(filename):
    atoms = {}
    reading = False
    with open(filename) as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                reading = True
                continue
            if reading and line.startswith("@<TRIPOS>"):
                break
            if reading:
                p = line.split()
                if len(p) >= 6:
                    atom_name = p[1]
                    atom_type = p[5]
                    element = atom_type.split(".")[0].upper()
                    if element == "H":
                        continue
                    try:
                        atoms[atom_name] = np.array([
                            float(p[2]),
                            float(p[3]),
                            float(p[4])
                        ])
                    except ValueError:
                        pass
    return atoms

reference = read_heavy_atoms(reference_file)
docked = read_heavy_atoms(docked_file)

common_atoms = [a for a in reference if a in docked]

if len(common_atoms) == 0:
    print("NA\tNA")
    sys.exit()

reference_xyz = np.array([reference[a] for a in common_atoms])
docked_xyz = np.array([docked[a] for a in common_atoms])

rmsd = np.sqrt(
    np.mean(
        np.sum(
            (reference_xyz - docked_xyz) ** 2,
            axis=1
        )
    )
)

reference_centroid = reference_xyz.mean(axis=0)
docked_centroid = docked_xyz.mean(axis=0)
centroid_distance = np.linalg.norm(reference_centroid - docked_centroid)

print(f"{rmsd:.3f}\t{centroid_distance:.3f}")
PY
}

# ============================================================
# START ALL RADII
# ============================================================

for RADIUS in "${RADII[@]}"; do

    MARGIN=$((RADIUS + 6))
    RUN_DIR="$BASE_DIR/radius${RADIUS}_negative_full"

    echo
    echo
    echo "===================================================="
    echo "1J3K NEGATIVE CONTROL"
    echo "STARTING RADIUS ${RADIUS} Å"
    echo "SHOWBOX MARGIN ${MARGIN} Å"
    echo "===================================================="

    rm -rf "$RUN_DIR"
    mkdir -p "$RUN_DIR"
    cd "$RUN_DIR"

    cp "$MASTER_SPHERES" master_spheres.sph
    cp "$NEGATIVE_REF" negative_site_reference.mol2
    cp "$LIGAND" ligand.mol2
    cp "$RECEPTOR" receptor.mol2

    # ========================================================
    # STEP 1 — SPHERE SELECTION
    # ========================================================

    echo
    echo "[1/4] Selecting negative-site spheres..."

    rm -f selected_spheres.sph

    if ! sphere_selector \
        master_spheres.sph \
        negative_site_reference.mol2 \
        "$RADIUS"
    then
        echo "Sphere selection FAILED."
        printf "%s\t%s\tNA\tNA\tNA\tNA\tFAILED_SPHERES\n" \
            "$RADIUS" "$MARGIN" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    if [[ ! -s selected_spheres.sph ]]; then
        echo "No selected sphere file produced."
        printf "%s\t%s\t0\tNA\tNA\tNA\tFAILED_SPHERES\n" \
            "$RADIUS" "$MARGIN" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    mv selected_spheres.sph "selected_spheres_radius${RADIUS}.sph"

    NUM_SPHERES=$(
        awk '/number of spheres in cluster/ {print $NF; exit}' \
        "selected_spheres_radius${RADIUS}.sph"
    )

    if [[ -z "${NUM_SPHERES:-}" ]]; then
        NUM_SPHERES="0"
    fi

    if [[ "$NUM_SPHERES" == "0" ]]; then
        printf "%s\t%s\t0\tNA\tNA\tNA\tFAILED_SPHERES\n" \
            "$RADIUS" "$MARGIN" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    echo "Selected spheres: ${NUM_SPHERES}"

    # ========================================================
    # STEP 2 — SHOWBOX
    # ========================================================

    echo
    echo "[2/4] Creating SHOWBOX..."

    cat > "showbox_radius${RADIUS}.in" <<EOF
Y
${MARGIN}
selected_spheres_radius${RADIUS}.sph
1
box_radius${RADIUS}.pdb
EOF

    if ! showbox \
        < "showbox_radius${RADIUS}.in" \
        > "showbox_radius${RADIUS}.out" 2>&1
    then
        echo "SHOWBOX FAILED."
        printf "%s\t%s\t%s\tNA\tNA\tNA\tFAILED_SHOWBOX\n" \
            "$RADIUS" "$MARGIN" "$NUM_SPHERES" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    if [[ ! -s "box_radius${RADIUS}.pdb" ]]; then
        echo "SHOWBOX did not create box file."
        printf "%s\t%s\t%s\tNA\tNA\tNA\tFAILED_SHOWBOX\n" \
            "$RADIUS" "$MARGIN" "$NUM_SPHERES" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    # ========================================================
    # STEP 3 — GRID
    # ========================================================

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

    if ! grid \
        -i "grid_radius${RADIUS}.in" \
        -o "grid_radius${RADIUS}.out"
    then
        echo "GRID FAILED."
        printf "%s\t%s\t%s\tNA\tNA\tNA\tFAILED_GRID\n" \
            "$RADIUS" "$MARGIN" "$NUM_SPHERES" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    if [[ ! -s "grid_radius${RADIUS}.nrg" ]]; then
        echo "GRID did not create energy grid."
        printf "%s\t%s\t%s\tNA\tNA\tNA\tFAILED_GRID\n" \
            "$RADIUS" "$MARGIN" "$NUM_SPHERES" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    # ========================================================
    # STEP 4 — DOCK6
    # ========================================================

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
simplex_trans_step                                            1
simplex_rot_step                                              0.1
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
ligand_outfile_prefix                                        WRA_negative_docked_radius${RADIUS}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

    dock6 \
        -i "dock_radius${RADIUS}.in" \
        -o "dock_radius${RADIUS}.out"

    DOCK_EXIT=$?
    SCORED_FILE="WRA_negative_docked_radius${RADIUS}_scored.mol2"

    if [[ $DOCK_EXIT -ne 0 || ! -s "$SCORED_FILE" ]]; then
        echo "DOCKING FAILED for radius ${RADIUS} Å."
        printf "%s\t%s\t%s\tNA\tNA\tNA\tFAILED_DOCKING\n" \
            "$RADIUS" "$MARGIN" "$NUM_SPHERES" >> "$RESULTS"
        cd "$BASE_DIR"
        continue
    fi

    GRID_SCORE=$(
        grep -i "Grid_Score" "$SCORED_FILE" |
        head -1 |
        awk '{print $NF}'
    )

    if [[ -z "${GRID_SCORE:-}" ]]; then
        GRID_SCORE="NA"
    fi

    METRICS=$(calculate_metrics "$LIGAND" "$SCORED_FILE")
    RMSD=$(echo "$METRICS" | awk '{print $1}')
    CENTROID=$(echo "$METRICS" | awk '{print $2}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS\n" \
        "$RADIUS" \
        "$MARGIN" \
        "$NUM_SPHERES" \
        "$GRID_SCORE" \
        "$RMSD" \
        "$CENTROID" \
        >> "$RESULTS"

    echo
    echo "----------------------------------------------------"
    echo "RADIUS ${RADIUS} Å COMPLETE"
    echo "----------------------------------------------------"
    echo "Sphere radius     : ${RADIUS} Å"
    echo "SHOWBOX margin    : ${MARGIN} Å"
    echo "Selected spheres  : ${NUM_SPHERES}"
    echo "Grid Score        : ${GRID_SCORE} kcal/mol"
    echo "Heavy-atom RMSD   : ${RMSD} Å"
    echo "Centroid distance : ${CENTROID} Å"
    echo "----------------------------------------------------"

    cd "$BASE_DIR"

done

# ============================================================
# FINAL SUMMARY
# ============================================================

cd "$BASE_DIR"

echo
echo "===================================================="
echo "ALL 1J3K NEGATIVE-CONTROL RUNS FINISHED"
echo "===================================================="
echo

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$RESULTS"
else
    cat "$RESULTS"
fi

echo
echo "Results saved as:"
echo "$RESULTS"
echo
