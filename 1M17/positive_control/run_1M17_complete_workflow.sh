#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================================
# 1M17 / AQ4 COMPLETE DOCK6 POSITIVE-CONTROL WORKFLOW
# Raw PDB -> receptor/ligand preparation -> DMS -> SPHGEN -> radii 4..28
# -> SHOWBOX -> GRID -> DOCK6 -> Grid Score / Heavy-atom RMSD / Centroid distance
# ============================================================================

BASE_DIR="$(pwd)"
RAW_PDB="${1:-$BASE_DIR/1M17.pdb}"
RADII=(4 6 8 10 12 14 16 18 20 22 24 26 28)

# User/site-specific paths may be overridden from the shell, e.g.
#   CHIMERA_BIN=/opt/UCSF-Chimera64-1.17.3/bin/chimera \
#   DMS_BIN=/usr/local/bin/dms \
#   DOCK6_PARAM_DIR=/path/to/dock6/parameters \
#   ./run_1M17_complete_workflow.sh 1M17.pdb
CHIMERA_BIN="${CHIMERA_BIN:-$(command -v chimera 2>/dev/null || true)}"
DMS_BIN="${DMS_BIN:-$(command -v dms 2>/dev/null || true)}"
PARAM_DIR="${DOCK6_PARAM_DIR:-/home/mdz/app/dock6/parameters}"

RECEPTOR_RAW="$BASE_DIR/1M17_receptor_noH.pdb"
LIGAND_RAW="$BASE_DIR/AQ4_native.pdb"
PREPARED_RECEPTOR="$BASE_DIR/1M17_receptor_prepared.mol2"
PREPARED_RECEPTOR_PDB="$BASE_DIR/1M17_receptor_prepared.pdb"
PREPARED_LIGAND="$BASE_DIR/AQ4_prepared.mol2"
PREPARED_LIGAND_PDB="$BASE_DIR/AQ4_prepared.pdb"
DMS_FILE="$BASE_DIR/1M17_receptor.dms"
SPH_FILE="$BASE_DIR/1M17_receptor.sph"
MASTER_SPHERES="$BASE_DIR/1M17_master_spheres.sph"
RESULTS="$BASE_DIR/1M17_positive_control_results.tsv"
CHIMERA_CMD="$BASE_DIR/chimera_prepare_1M17.cmd"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || fail "Required command '$1' was not found in PATH."
}

[[ -s "$RAW_PDB" ]] || fail "Raw PDB not found or empty: $RAW_PDB"

# DOCK6-side dependencies
need_cmd sphgen
need_cmd sphere_selector
need_cmd showbox
need_cmd grid
need_cmd dock6
need_cmd python3

[[ -n "$CHIMERA_BIN" && -x "$CHIMERA_BIN" ]] || fail \
  "UCSF Chimera was not found. Install Chimera or set CHIMERA_BIN=/full/path/to/chimera."
[[ -n "$DMS_BIN" && -x "$DMS_BIN" ]] || fail \
  "The legacy dms program was not found. Install dms or set DMS_BIN=/full/path/to/dms."

for f in \
  "$PARAM_DIR/vdw_AMBER_parm99.defn" \
  "$PARAM_DIR/flex.defn" \
  "$PARAM_DIR/flex_drive.tbl"
do
  [[ -s "$f" ]] || fail "Missing DOCK6 parameter file: $f"
done

# Protect an existing complete run from accidental overwrite.
if [[ -e "$RESULTS" ]]; then
  fail "Results file already exists: $RESULTS. Move/rename old results before a new complete run."
fi
for r in "${RADII[@]}"; do
  [[ ! -e "$BASE_DIR/radius${r}_full" ]] || fail \
    "Run directory already exists: $BASE_DIR/radius${r}_full. Move/rename it before rerunning."
done

# ============================================================================
# STAGE 0A - EXTRACT 1M17 CHAIN A RECEPTOR AND NATIVE AQ4 FROM RAW PDB
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0A - EXTRACTING RECEPTOR AND AQ4"
echo "============================================================"

# Receptor: protein ATOM records from chain A only; waters/ligand/other HETATM
# records are intentionally excluded.
awk '
  /^ATOM  / && substr($0,22,1)=="A" {print}
  END {print "END"}
' "$RAW_PDB" > "$RECEPTOR_RAW"

# Native ligand: AQ4, chain A.
awk '
  /^HETATM/ && substr($0,18,3)=="AQ4" && substr($0,22,1)=="A" {print}
  END {print "END"}
' "$RAW_PDB" > "$LIGAND_RAW"

[[ $(grep -c '^ATOM  ' "$RECEPTOR_RAW" || true) -gt 0 ]] || fail \
  "No chain-A protein ATOM records were extracted from $RAW_PDB."
[[ $(grep -c '^HETATM' "$LIGAND_RAW" || true) -gt 0 ]] || fail \
  "AQ4 in chain A was not found in $RAW_PDB."

echo "Receptor : $RECEPTOR_RAW"
echo "AQ4      : $LIGAND_RAW"

# ============================================================================
# STAGE 0B - CHIMERA PREPARATION
# Receptor: hydrogens + ff14SB charges/types
# AQ4: hydrogens + AM1-BCC charges, formal charge 0
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0B - UCSF CHIMERA PREPARATION"
echo "============================================================"

cat > "$CHIMERA_CMD" <<EOF_CHIMERA
open noprefs $RECEPTOR_RAW
addh
addcharge std chargeModel ff14SB
write format mol2 atomTypes amber #0 $PREPARED_RECEPTOR
write format pdb #0 $PREPARED_RECEPTOR_PDB
close all
open noprefs $LIGAND_RAW
addh
addcharge nonstd :AQ4 0 method am1
write format mol2 atomTypes amber #0 $PREPARED_LIGAND
write format pdb #0 $PREPARED_LIGAND_PDB
stop
EOF_CHIMERA

"$CHIMERA_BIN" --nogui "$CHIMERA_CMD"

[[ -s "$PREPARED_RECEPTOR" ]] || fail "Chimera did not create $PREPARED_RECEPTOR"
[[ -s "$PREPARED_LIGAND" ]] || fail "Chimera did not create $PREPARED_LIGAND"

echo "Prepared receptor : $PREPARED_RECEPTOR"
echo "Prepared AQ4      : $PREPARED_LIGAND"

# ============================================================================
# STAGE 0C - MOLECULAR SURFACE WITH NORMALS
# Standard DOCK6/dms route uses the no-H receptor PDB.
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0C - GENERATING RECEPTOR MOLECULAR SURFACE"
echo "============================================================"

"$DMS_BIN" "$RECEPTOR_RAW" -n -w 1.4 -v -o "$DMS_FILE"
[[ -s "$DMS_FILE" ]] || fail "DMS did not create $DMS_FILE"

echo "DMS surface : $DMS_FILE"

# ============================================================================
# STAGE 0D - SPHGEN
# Exact INSPH values reconstructed from the original 1M17 workflow.
# ============================================================================

echo
echo "============================================================"
echo "STAGE 0D - GENERATING MASTER SPHERES WITH SPHGEN"
echo "============================================================"

cd "$BASE_DIR"
rm -f OUTSPH "$SPH_FILE" "$MASTER_SPHERES"

cat > INSPH <<EOF_INSPH
1M17_receptor.dms
R
X
0.0
4.0
1.4
1M17_receptor.sph
EOF_INSPH

sphgen
[[ -s "$SPH_FILE" ]] || fail "SPHGEN did not create $SPH_FILE"
cp "$SPH_FILE" "$MASTER_SPHERES"
[[ -s "$MASTER_SPHERES" ]] || fail "Could not create $MASTER_SPHERES"

echo "Master spheres : $MASTER_SPHERES"

# ============================================================================
# METRIC FUNCTION
# Heavy atoms matched by atom name; no superposition/alignment.
# ============================================================================

calc_metrics() {
  python3 - "$1" "$2" <<'PY'
import sys
import math

ref_file, dock_file = sys.argv[1], sys.argv[2]

def read_heavy_atoms(path):
    atoms = {}
    reading = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>ATOM"):
                reading = True
                continue
            if reading and line.startswith("@<TRIPOS>"):
                break
            if not reading:
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            name = fields[1]
            atom_type = fields[5]
            element = atom_type.split(".")[0].upper()
            if element == "H":
                continue
            try:
                xyz = tuple(float(fields[i]) for i in (2, 3, 4))
            except ValueError:
                continue
            atoms[name] = xyz
    return atoms

ref = read_heavy_atoms(ref_file)
dock = read_heavy_atoms(dock_file)
common = [name for name in ref if name in dock]

if not common:
    print("NA\tNA")
    raise SystemExit

sq = 0.0
for name in common:
    sq += sum((ref[name][i] - dock[name][i]) ** 2 for i in range(3))
rmsd = math.sqrt(sq / len(common))

ref_centroid = [sum(ref[n][i] for n in common) / len(common) for i in range(3)]
dock_centroid = [sum(dock[n][i] for n in common) / len(common) for i in range(3)]
centroid = math.sqrt(sum((ref_centroid[i] - dock_centroid[i]) ** 2 for i in range(3)))

print(f"{rmsd:.3f}\t{centroid:.3f}")
PY
}

# ============================================================================
# RESULTS FILE
# ============================================================================

printf "Radius\tMargin\tGridScore\tHeavyAtomRMSD\tCentroidDistance\tStatus\n" > "$RESULTS"

# ============================================================================
# STAGE 1 - POSITIVE-CONTROL DOCKING FOR RADII 4..28 A
# ============================================================================

for RADIUS in "${RADII[@]}"; do
  MARGIN=$((RADIUS + 6))
  RUN_DIR="$BASE_DIR/radius${RADIUS}_full"

  echo
  echo "============================================================"
  echo "1M17 POSITIVE CONTROL - RADIUS ${RADIUS} A"
  echo "SHOWBOX margin: ${MARGIN} A"
  echo "============================================================"

  mkdir -p "$RUN_DIR"
  cd "$RUN_DIR"

  # --------------------------------------------------------------------------
  # 1/4 SPHERE SELECTION
  # --------------------------------------------------------------------------
  echo "[1/4] Selecting spheres"
  if ! sphere_selector "$MASTER_SPHERES" "$PREPARED_LIGAND" "$RADIUS"; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_SPHERES\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi
  if [[ ! -s selected_spheres.sph ]]; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_SPHERES\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi
  mv selected_spheres.sph "selected_spheres_radius${RADIUS}.sph"

  # --------------------------------------------------------------------------
  # 2/4 SHOWBOX
  # --------------------------------------------------------------------------
  echo "[2/4] Creating SHOWBOX box"
  cat > "showbox_radius${RADIUS}.in" <<EOF_SHOWBOX
Y
${MARGIN}
selected_spheres_radius${RADIUS}.sph
1
box_radius${RADIUS}.pdb
EOF_SHOWBOX

  if ! showbox < "showbox_radius${RADIUS}.in" > "showbox_radius${RADIUS}.out" 2>&1; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_SHOWBOX\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi
  if [[ ! -s "box_radius${RADIUS}.pdb" ]]; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_SHOWBOX\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi

  # --------------------------------------------------------------------------
  # 3/4 GRID
  # --------------------------------------------------------------------------
  echo "[3/4] Running GRID"
  cat > "grid_radius${RADIUS}.in" <<EOF_GRID
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
EOF_GRID

  if ! grid -i "grid_radius${RADIUS}.in" -o "grid_radius${RADIUS}.out"; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_GRID\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi
  if [[ ! -s "grid_radius${RADIUS}.nrg" ]]; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_GRID\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi

  # --------------------------------------------------------------------------
  # 4/4 DOCK6
  # --------------------------------------------------------------------------
  echo "[4/4] Running DOCK6"
  cat > "dock_radius${RADIUS}.in" <<EOF_DOCK
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
ligand_outfile_prefix                                         AQ4_docked_radius${RADIUS}
write_orientations                                            no
num_scored_conformers                                         1
rank_ligands                                                  no
EOF_DOCK

  if ! dock6 -i "dock_radius${RADIUS}.in" -o "dock_radius${RADIUS}.out"; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_DOCKING\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi

  SCORED_FILE="AQ4_docked_radius${RADIUS}_scored.mol2"
  if [[ ! -s "$SCORED_FILE" ]]; then
    printf "%s\t%s\tNA\tNA\tNA\tFAILED_DOCKING\n" "$RADIUS" "$MARGIN" >> "$RESULTS"
    cd "$BASE_DIR"
    continue
  fi

  GRID_SCORE=$(grep -i "Grid_Score" "$SCORED_FILE" | head -1 | awk '{print $NF}' || true)
  GRID_SCORE="${GRID_SCORE:-NA}"

  METRICS=$(calc_metrics "$PREPARED_LIGAND" "$SCORED_FILE")
  RMSD=$(printf "%s\n" "$METRICS" | awk '{print $1}')
  CENTROID=$(printf "%s\n" "$METRICS" | awk '{print $2}')
  RMSD="${RMSD:-NA}"
  CENTROID="${CENTROID:-NA}"

  printf "%s\t%s\t%s\t%s\t%s\tSUCCESS\n" \
    "$RADIUS" "$MARGIN" "$GRID_SCORE" "$RMSD" "$CENTROID" >> "$RESULTS"

  echo "Radius             : ${RADIUS} A"
  echo "Grid score         : ${GRID_SCORE} kcal/mol"
  echo "Heavy-atom RMSD    : ${RMSD} A"
  echo "Centroid distance  : ${CENTROID} A"

  cd "$BASE_DIR"
done

# ============================================================================
# FINAL RESULTS
# ============================================================================

cd "$BASE_DIR"
echo
echo "============================================================"
echo "1M17 POSITIVE CONTROL - FINAL RESULTS"
echo "============================================================"
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$RESULTS"
else
  cat "$RESULTS"
fi

echo
echo "Results saved to: $RESULTS"
echo "Workflow complete."
echo "============================================================"
