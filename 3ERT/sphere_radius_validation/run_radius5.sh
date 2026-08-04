#!/usr/bin/env bash
set -euo pipefail

RADIUS="5.0"
RUN_DIR="radius5"

MASTER_SPHERES="master_sphere_generation/3ERT_master_spheres.sph"
DOCK_TEMPLATE="../dock_OHT_m15.in"

echo "Creating ${RUN_DIR}..."

if [[ -d "$RUN_DIR" ]]; then
    echo "ERROR: $RUN_DIR already exists."
    echo "Rename or remove it before rerunning."
    exit 1
fi

mkdir "$RUN_DIR"

cp "$MASTER_SPHERES" "$RUN_DIR/"
cp OHT_prepared.mol2 "$RUN_DIR/"
cp OHT_native.pdb "$RUN_DIR/"
cp 3ERT_receptor_prepared.mol2 "$RUN_DIR/"
cp 3ERT_OHT_box_m15.pdb "$RUN_DIR/"
cp 3ERT_grid_m15.bmp "$RUN_DIR/"
cp 3ERT_grid_m15.nrg "$RUN_DIR/"
cp "$DOCK_TEMPLATE" "$RUN_DIR/dock_radius5.in"

cd "$RUN_DIR"

echo
echo "Selecting spheres within ${RADIUS} Å of OHT..."

sphere_selector \
    3ERT_master_spheres.sph \
    OHT_prepared.mol2 \
    "$RADIUS"

if [[ ! -s selected_spheres.sph ]]; then
    echo "ERROR: selected_spheres.sph was not generated."
    exit 1
fi

mv selected_spheres.sph selected_spheres_OHT_radius5.sph

SPHERE_COUNT=$(
    awk '/number of spheres in cluster/ {print $NF; exit}' \
    selected_spheres_OHT_radius5.sph
)

echo "Selected sphere count: $SPHERE_COUNT"

sed -i \
's|^receptor_site_file.*|receptor_site_file                                           selected_spheres_OHT_radius5.sph|' \
dock_radius5.in

sed -i \
's|^ligand_outfile_prefix.*|ligand_outfile_prefix                                        OHT_docked_radius5|' \
dock_radius5.in

sed -i \
's|^bump_grid_prefix.*|bump_grid_prefix                                             3ERT_grid_m15|' \
dock_radius5.in

sed -i \
's|^grid_score_grid_prefix.*|grid_score_grid_prefix                                       3ERT_grid_m15|' \
dock_radius5.in

echo
echo "Running DOCK6..."

dock6 -i dock_radius5.in -o dock_radius5.out

DOCKED_FILE="OHT_docked_radius5_scored.mol2"

if [[ ! -s "$DOCKED_FILE" ]]; then
    echo "ERROR: Docked output was not generated."
    tail -40 dock_radius5.out
    exit 1
fi

echo
echo "Calculating final metrics..."

python3 - <<'PYTHON'
import math
import re
from pathlib import Path

import numpy as np

BOX_FILE = Path("3ERT_OHT_box_m15.pdb")
NATIVE_FILE = Path("OHT_native.pdb")
DOCKED_FILE = Path("OHT_docked_radius5_scored.mol2")
RESULT_FILE = Path("radius5_results.txt")


def pdb_element(line):
    element = line[76:78].strip() if len(line) >= 78 else ""

    if element:
        return element.upper()

    name = re.sub(r"^[0-9]+", "", line[12:16].strip())
    return name[0].upper() if name else ""


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
                pass

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
            atoms.setdefault(
                atom_name,
                np.asarray([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            )
        except ValueError:
            continue

    return atoms


def read_docked(path):
    atoms = {}
    properties = {}
    in_atoms = False
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
            in_atoms = False
            continue

        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue

        if line.startswith("@<TRIPOS>"):
            in_atoms = False
            continue

        if not in_atoms:
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

    return atoms, properties


def get_number(properties, names):
    for name in names:
        if name not in properties:
            continue

        match = re.search(
            r"[-+]?(?:\d+\.\d+|\d+|\.\d+)(?:[Ee][-+]?\d+)?",
            properties[name]
        )

        if match:
            return float(match.group())

    return None


box_coordinates = read_box(BOX_FILE)
box_min = box_coordinates.min(axis=0)
box_max = box_coordinates.max(axis=0)
box_centre = (box_min + box_max) / 2
box_dimensions = box_max - box_min

native_atoms = read_native(NATIVE_FILE)
docked_atoms, properties = read_docked(DOCKED_FILE)

common_atoms = [
    atom_name
    for atom_name in native_atoms
    if atom_name in docked_atoms
]

if len(common_atoms) != len(native_atoms):
    missing = sorted(set(native_atoms) - set(docked_atoms))

    raise RuntimeError(
        f"Matched {len(common_atoms)} of {len(native_atoms)} heavy atoms. "
        f"Missing: {missing}"
    )

native_xyz = np.asarray([native_atoms[name] for name in common_atoms])
docked_xyz = np.asarray([docked_atoms[name] for name in common_atoms])

difference = docked_xyz - native_xyz
rmsd = math.sqrt(np.mean(np.sum(difference ** 2, axis=1)))

native_centroid = native_xyz.mean(axis=0)
docked_centroid = docked_xyz.mean(axis=0)
centroid_distance = np.linalg.norm(docked_centroid - native_centroid)

grid_score = get_number(properties, ["Grid_Score"])
grid_vdw = get_number(properties, ["Grid_vdw_energy", "Grid_vdw"])
grid_es = get_number(properties, ["Grid_es_energy", "Grid_es"])
internal_repulsive = get_number(
    properties,
    ["Internal_energy_repulsive", "Internal_energy"]
)


def display(value):
    return "Not reported" if value is None else f"{value:.6f}"


report = f"""============================================================
3ERT positive-control sphere-radius test
Sphere-selection radius: 5.0 Å
============================================================

Box centre:                  {box_centre[0]:.3f}, {box_centre[1]:.3f}, {box_centre[2]:.3f}
Box dimensions:              {box_dimensions[0]:.3f} × {box_dimensions[1]:.3f} × {box_dimensions[2]:.3f} Å
Grid Score:                  {display(grid_score)} kcal/mol
Grid van der Waals energy:   {display(grid_vdw)}
Grid electrostatic energy:   {display(grid_es)}
Internal repulsive energy:   {display(internal_repulsive)}
Heavy-atom RMSD:             {rmsd:.3f} Å
Centroid distance:           {centroid_distance:.3f} Å
Heavy-atom count:            {len(common_atoms)}

============================================================"""

print()
print(report)
print()

RESULT_FILE.write_text(report + "\n")
PYTHON

echo "Radius 5 completed successfully."
echo "Results saved in: $(pwd)/radius5_results.txt"
