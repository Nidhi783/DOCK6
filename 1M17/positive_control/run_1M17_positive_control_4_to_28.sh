#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# 1M17 POSITIVE CONTROL
#
# Sphere radii:
# 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28 A
#
# SHOWBOX margin = radius + 6 A
#
# Final outputs:
# - Grid Score
# - Heavy-atom RMSD
# - Centroid-to-centroid distance
#
# Summary file:
# 1M17_positive_control_results.tsv
# ============================================================


BASE_DIR="$(pwd)"

RADII=(4 6 8 10 12 14 16 18 20 22 24 26 28)

MASTER_SPHERES="$BASE_DIR/1M17_master_spheres.sph"
PREPARED_LIGAND="$BASE_DIR/AQ4_prepared.mol2"
PREPARED_RECEPTOR="$BASE_DIR/1M17_receptor_prepared.mol2"

PARAM_DIR="/home/mdz/app/dock6/parameters"

RESULTS="$BASE_DIR/1M17_positive_control_results.tsv"


# ============================================================
# CHECK REQUIRED FILES
# ============================================================

for f in \
"$MASTER_SPHERES" \
"$PREPARED_LIGAND" \
"$PREPARED_RECEPTOR"
do
if [[ ! -s "$f" ]]; then
echo "ERROR: Missing or empty file: $f"
exit 1
fi
done


# ============================================================
# CALCULATE HEAVY-ATOM RMSD + CENTROID DISTANCE
#
# Heavy atoms are matched by atom name.
# Hydrogens are excluded.
# No structural superposition/alignment is performed.
# ============================================================

calc_metrics() {
  
  python3 - "$1" "$2" <<'PY'
  
  import sys
  import numpy as np
  
  ref_file = sys.argv[1]
  dock_file = sys.argv[2]
  
  
  def read_heavy_atoms(path):
    
    atoms = {}
  reading = False
  
  with open(path) as f:
    
    for line in f:
    
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
  
  atom_name = fields[1]
  atom_type = fields[5]
  
  element = atom_type.split(".")[0].upper()
  
  # Exclude hydrogens
  if element == "H":
    continue
  
  try:
    
    xyz = np.array(
      [
        float(fields[2]),
        float(fields[3]),
        float(fields[4])
      ],
      dtype=float
    )
  
  except ValueError:
    continue
  
  atoms[atom_name] = xyz
  
  return atoms
  
  
  ref = read_heavy_atoms(ref_file)
  dock = read_heavy_atoms(dock_file)
  
  
  common = [
    name
    for name in ref
    if name in dock
  ]
  
  
  if not common:
    
    print("NA\tNA")
  raise SystemExit
  
  
  ref_xyz = np.array(
    [ref[name] for name in common]
  )
  
  dock_xyz = np.array(
    [dock[name] for name in common]
  )
  
  
  # ============================================================
  # HEAVY-ATOM RMSD
  # ============================================================
  
  rmsd = np.sqrt(
    np.mean(
      np.sum(
        (ref_xyz - dock_xyz) ** 2,
        axis=1
      )
    )
  )
  
  
  # ============================================================
  # CENTROID-TO-CENTROID DISTANCE
  # ============================================================
  
  ref_centroid = ref_xyz.mean(axis=0)
  dock_centroid = dock_xyz.mean(axis=0)
  
  centroid_distance = np.linalg.norm(
    ref_centroid - dock_centroid
  )
  
  
  print(
    f"{rmsd:.3f}\t{centroid_distance:.3f}"
  )
  
  PY
}


# ============================================================
# RESULTS FILE
# ============================================================

if [[ -e "$RESULTS" ]]; then

echo
echo "ERROR: Results file already exists:"
echo "$RESULTS"
echo
echo "Rename or move it before starting a new complete run."
exit 1

fi


printf \
"Radius\tMargin\tGridScore\tHeavyAtomRMSD\tCentroidDistance\tStatus\n" \
> "$RESULTS"


# ============================================================
# RUN ALL RADII
# ============================================================

for RADIUS in "${RADII[@]}"; do

MARGIN=$((RADIUS + 6))
RUN_DIR="$BASE_DIR/radius${RADIUS}_full"


echo
echo "========================================"
echo "1M17 POSITIVE CONTROL"
echo "Sphere radius: ${RADIUS} A"
echo "SHOWBOX margin: ${MARGIN} A"
echo "========================================"


# ========================================================
# SAFETY CHECK
# ========================================================

if [[ -e "$RUN_DIR" ]]; then

echo
echo "ERROR: Run directory already exists:"
echo "$RUN_DIR"
echo
echo "Existing results will NOT be overwritten."
exit 1

fi


mkdir -p "$RUN_DIR"
cd "$RUN_DIR"


# ========================================================
# STEP 1 — SPHERE SELECTION
# ========================================================

echo
echo "[1/4] Selecting spheres"


if ! sphere_selector \
"$MASTER_SPHERES" \
"$PREPARED_LIGAND" \
"$RADIUS"
then

echo "ERROR: sphere_selector failed"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_SPHERES\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


if [[ ! -s selected_spheres.sph ]]; then

echo "ERROR: sphere_selector did not create selected_spheres.sph"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_SPHERES\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


mv \
selected_spheres.sph \
"selected_spheres_radius${RADIUS}.sph"


# ========================================================
# STEP 2 — SHOWBOX
# ========================================================

echo
echo "[2/4] Creating SHOWBOX box"


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

echo "ERROR: SHOWBOX failed"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_SHOWBOX\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


if [[ ! -s "box_radius${RADIUS}.pdb" ]]; then

echo "ERROR: SHOWBOX failed"

cat "showbox_radius${RADIUS}.out"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_SHOWBOX\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


# ========================================================
# STEP 3 — GRID
# ========================================================

echo
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


if ! grid \
-i "grid_radius${RADIUS}.in" \
-o "grid_radius${RADIUS}.out"
then

echo "ERROR: GRID failed"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_GRID\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


if [[ ! -s "grid_radius${RADIUS}.nrg" ]]; then

echo "ERROR: GRID failed"

tail -40 \
"grid_radius${RADIUS}.out" \
|| true

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_GRID\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


# ========================================================
# STEP 4 — DOCK6
# ========================================================

echo
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
ligand_outfile_prefix                                         AQ4_docked_radius${RADIUS}
write_orientations                                            no
num_scored_conformers                                         1
rank_ligands                                                  no
EOF


if ! dock6 \
-i "dock_radius${RADIUS}.in" \
-o "dock_radius${RADIUS}.out"
then

echo "ERROR: DOCK6 failed"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_DOCKING\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


SCORED_FILE="AQ4_docked_radius${RADIUS}_scored.mol2"


if [[ ! -s "$SCORED_FILE" ]]; then

echo "ERROR: DOCK6 did not create:"
echo "$SCORED_FILE"

printf \
"%s\t%s\tNA\tNA\tNA\tFAILED_DOCKING\n" \
"$RADIUS" \
"$MARGIN" \
>> "$RESULTS"

cd "$BASE_DIR"
continue
fi


# ========================================================
# GRID SCORE
# ========================================================

GRID_SCORE=$(
  grep -i "Grid_Score" \
  "$SCORED_FILE" \
  | head -1 \
  | awk '{print $NF}'
)

GRID_SCORE="${GRID_SCORE:-NA}"


# ========================================================
# RMSD + CENTROID DISTANCE
# ========================================================

METRICS=$(
  calc_metrics \
  "$PREPARED_LIGAND" \
  "$SCORED_FILE"
)


RMSD=$(
  printf "%s\n" "$METRICS" \
  | awk '{print $1}'
)


CENTROID=$(
  printf "%s\n" "$METRICS" \
  | awk '{print $2}'
)


RMSD="${RMSD:-NA}"
CENTROID="${CENTROID:-NA}"


# ========================================================
# SAVE RESULTS
# ========================================================

printf \
"%s\t%s\t%s\t%s\t%s\tSUCCESS\n" \
"$RADIUS" \
"$MARGIN" \
"$GRID_SCORE" \
"$RMSD" \
"$CENTROID" \
>> "$RESULTS"


echo
echo "========================================"
echo "RADIUS ${RADIUS} COMPLETE"
echo "========================================"
echo "Margin              : ${MARGIN} A"
echo "Grid score          : ${GRID_SCORE} kcal/mol"
echo "Heavy-atom RMSD     : ${RMSD} A"
echo "Centroid distance   : ${CENTROID} A"
echo "Folder              : $RUN_DIR"
echo "Dock output         : $RUN_DIR/dock_radius${RADIUS}.out"
echo "Docked ligand       : $RUN_DIR/$SCORED_FILE"
echo "========================================"


cd "$BASE_DIR"

done


# ============================================================
# FINAL RESULTS
# ============================================================

cd "$BASE_DIR"

echo
echo
echo "============================================================"
echo "1M17 POSITIVE CONTROL — FINAL RESULTS"
echo "============================================================"
echo


if command -v column >/dev/null 2>&1; then

column \
-t \
-s $'\t' \
"$RESULTS"

else
  
  cat "$RESULTS"

fi


echo
echo "============================================================"
echo "Results saved to:"
echo "$RESULTS"
echo "============================================================"