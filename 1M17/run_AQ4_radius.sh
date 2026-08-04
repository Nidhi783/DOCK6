#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./run_AQ4_radius.sh 5

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <sphere-radius>"
    echo "Example: $0 5"
    exit 1
fi

RADIUS="$1"

if ! [[ "$RADIUS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Radius must be a whole number."
    exit 1
fi

RUN_DIR="radius${RADIUS}"
SPHERE_FILE="selected_spheres_AQ4_radius${RADIUS}.sph"
DOCK_INPUT="dock_AQ4_radius${RADIUS}.in"
DOCK_OUTPUT="dock_AQ4_radius${RADIUS}.out"
DOCK_PREFIX="AQ4_docked_radius${RADIUS}"
DOCKED_FILE="${DOCK_PREFIX}_scored.mol2"
RESULT_FILE="radius${RADIUS}_results.txt"

MASTER_SPHERES="1M17_receptor.sph"
DOCK_TEMPLATE="dock_AQ4_radius4.in"

echo
echo "============================================================"
echo "1M17–AQ4 sphere-radius test"
echo "Sphere-selection radius: ${RADIUS}.0 Å"
echo "============================================================"

# Check required files
for FILE in \
    "$MASTER_SPHERES" \
    "$DOCK_TEMPLATE" \
    "AQ4_prepared.mol2" \
    "AQ4_native.pdb" \
    "1M17_receptor_prepared.mol2" \
    "1M17_receptor_prepared.pdb" \
    "1M17_AQ4_box_m15.pdb" \
    "1M17_grid_m15.bmp" \
    "1M17_grid_m15.nrg"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty file: $FILE"
        exit 1
    fi
done

# Prevent accidental overwriting
if [[ -d "$RUN_DIR" ]]; then
    echo "ERROR: Folder $RUN_DIR already exists."
    echo "Rename it before rerunning this radius."
    exit 1
fi

mkdir "$RUN_DIR"

# Copy fixed inputs
cp "$MASTER_SPHERES" "$RUN_DIR/"
cp AQ4_prepared.mol2 "$RUN_DIR/"
cp AQ4_native.pdb "$RUN_DIR/"
cp 1M17_receptor_prepared.mol2 "$RUN_DIR/"
cp 1M17_receptor_prepared.pdb "$RUN_DIR/"
cp 1M17_AQ4_box_m15.pdb "$RUN_DIR/"
cp 1M17_grid_m15.bmp "$RUN_DIR/"
cp 1M17_grid_m15.nrg "$RUN_DIR/"
cp "$DOCK_TEMPLATE" "$RUN_DIR/$DOCK_INPUT"

cd "$RUN_DIR"

echo
echo "Selecting spheres within ${RADIUS}.0 Å of AQ4..."

rm -f selected_spheres.sph

sphere_selector \
    1M17_receptor.sph \
    AQ4_prepared.mol2 \
    "${RADIUS}.0"

if [[ ! -s selected_spheres.sph ]]; then
    echo "ERROR: Sphere selection failed."
    exit 1
fi

mv selected_spheres.sph "$SPHERE_FILE"

SPHERE_COUNT=$(
    awk '
    /number of spheres in cluster/ {
        print $NF
        exit
    }' "$SPHERE_FILE"
)

echo "Selected sphere count: $SPHERE_COUNT"

# Modify only radius-specific docking settings
sed -i \
    "s|^receptor_site_file.*|receptor_site_file                                           ${SPHERE_FILE}|" \
    "$DOCK_INPUT"

sed -i \
    "s|^ligand_outfile_prefix.*|ligand_outfile_prefix                                        ${DOCK_PREFIX}|" \
    "$DOCK_INPUT"

# Ensure the same fixed grid is used
sed -i \
    "s|^bump_grid_prefix.*|bump_grid_prefix                                             1M17_grid_m15|" \
    "$DOCK_INPUT"

sed -i \
    "s|^grid_score_grid_prefix.*|grid_score_grid_prefix                                       1M17_grid_m15|" \
    "$DOCK_INPUT"

echo
echo "Running DOCK6..."

dock6 -i "$DOCK_INPUT" -o "$DOCK_OUTPUT"

if [[ ! -s "$DOCKED_FILE" ]]; then
    echo "ERROR: Docked MOL2 output was not generated."
    tail -40 "$DOCK_OUTPUT" || true
    exit 1
fi

echo "Docking completed."
echo "Calculating final parameters..."

python3 - "$RADIUS" "$SPHERE_COUNT" <<'PYTHON'
import math
import re
import sys
from pathlib import Path

import numpy as np

radius = int(sys.argv[1])
sphere_count = sys.argv[2]

box_file = Path("1M17_AQ4_box_m15.pdb")
native_file = Path("AQ4_native.pdb")
docked_file = Path(f"AQ4_docked_radius{radius}_scored.mol2")
result_file = Path(f"radius{radius}_results.txt")


def pdb_element(line):
    element = line[76:78].strip() if len(line) >= 78 else ""

    if element:
        return element.upper()

    atom_name = re.sub(r"^[0-9]+", "", line[12:16].strip())
    return atom_name[0].upper() if atom_name else ""


def read_box(path):
    coordinates = []

    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                coordinates.append([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            except ValueError:
                continue

    if not coordinates:
        raise RuntimeError("No coordinates found in box file.")

    return np.asarray(coordinates, dtype=float)


def read_native(path):
    atoms = {}

    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue

        if pdb_element(line) == "H":
            continue

        atom_name = line[12:16].strip()

        try:
            coordinate = np.asarray([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ])
        except ValueError:
            continue

        if atom_name not in atoms:
            atoms[atom_name] = coordinate

    if not atoms:
        raise RuntimeError("No native heavy atoms found.")

    return atoms


def read_docked(path):
    atoms = {}
    properties = {}
    in_atom_section = False
    first_pose_started = False

    for line in path.read_text().splitlines():
        if line.startswith("##########"):
            match = re.match(r"##########\s+([^:]+):\s*(.*)", line)

            if match:
                properties[match.group(1).strip()] = match.group(2).strip()

            continue

        if line.startswith("@<TRIPOS>MOLECULE"):
            if first_pose_started and atoms:
                break

            first_pose_started = True
            in_atom_section = False
            continue

        if line.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            continue

        if line.startswith("@<TRIPOS>"):
            in_atom_section = False
            continue

        if not in_atom_section:
            continue

        fields = line.split()

        if len(fields) < 6:
            continue

        atom_name = fields[1]
        atom_type = fields[5]
        element = atom_type.split(".")[0].upper()

        if element == "H":
            continue

        atoms[atom_name] = np.asarray([
            float(fields[2]),
            float(fields[3]),
            float(fields[4]),
        ])

    if not atoms:
        raise RuntimeError("No docked heavy atoms found.")

    return atoms, properties


def property_number(properties, names):
    for name in names:
        if name not in properties:
            continue

        match = re.search(
            r"[-+]?(?:\d+\.\d+|\d+|\.\d+)(?:[Ee][-+]?\d+)?",
            properties[name],
        )

        if match:
            return float(match.group())

    return None


box_coordinates = read_box(box_file)
box_min = box_coordinates.min(axis=0)
box_max = box_coordinates.max(axis=0)

box_centre = (box_min + box_max) / 2.0
box_dimensions = box_max - box_min

native_atoms = read_native(native_file)
docked_atoms, properties = read_docked(docked_file)

common_atoms = [
    atom_name
    for atom_name in native_atoms
    if atom_name in docked_atoms
]

if len(common_atoms) != len(native_atoms):
    missing = sorted(set(native_atoms) - set(docked_atoms))

    raise RuntimeError(
        f"Matched only {len(common_atoms)} of "
        f"{len(native_atoms)} native heavy atoms.\n"
        f"Missing atoms: {missing}"
    )

native_xyz = np.asarray([
    native_atoms[name] for name in common_atoms
])

docked_xyz = np.asarray([
    docked_atoms[name] for name in common_atoms
])

difference = docked_xyz - native_xyz

heavy_atom_rmsd = math.sqrt(
    np.mean(np.sum(difference ** 2, axis=1))
)

native_centroid = native_xyz.mean(axis=0)
docked_centroid = docked_xyz.mean(axis=0)

centroid_distance = float(
    np.linalg.norm(docked_centroid - native_centroid)
)

grid_score = property_number(
    properties,
    ["Grid_Score", "Grid Score"],
)

grid_vdw = property_number(
    properties,
    ["Grid_vdw_energy", "Grid_vdw"],
)

grid_es = property_number(
    properties,
    ["Grid_es_energy", "Grid_es"],
)

internal_repulsive = property_number(
    properties,
    ["Internal_energy_repulsive", "Internal_energy"],
)


def display(value):
    return "Not reported" if value is None else f"{value:.6f}"


report = f"""============================================================
1M17–AQ4 positive-control sphere-radius test
Sphere-selection radius:    {radius}.0 Å
Selected sphere count:      {sphere_count}
============================================================

Box centre:                  {box_centre[0]:.3f}, {box_centre[1]:.3f}, {box_centre[2]:.3f}
Box dimensions:              {box_dimensions[0]:.3f} × {box_dimensions[1]:.3f} × {box_dimensions[2]:.3f} Å
Grid Score:                  {display(grid_score)} kcal/mol
Grid van der Waals energy:   {display(grid_vdw)}
Grid electrostatic energy:   {display(grid_es)}
Internal repulsive energy:   {display(internal_repulsive)}
Heavy-atom RMSD:             {heavy_atom_rmsd:.3f} Å
Centroid distance:           {centroid_distance:.3f} Å
Heavy-atom count:            {len(common_atoms)}

============================================================"""

print()
print(report)
print()

result_file.write_text(report + "\n")
PYTHON

echo "Radius ${RADIUS} completed successfully."
echo "Results saved in:"
echo "  $(pwd)/$RESULT_FILE"
