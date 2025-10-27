"""
    MDShift.py
    Purpose: Calculate changes in CoM and MoI upon isotopic labeling for explaining/predicting isotopic shifts in high-resolution ion mobility separations
    Usage: MDShift.py <labels.txt> <solvent radius>
    Authors: Noah D. Roberts, Gabe Nagy
    Contact: noah.roberts@utah.edu, gabe.nagy@utah.edu
    See readme.md for more information
"""

import os, sys, glob, re, datetime, csv, math
import pymol2  # Embedded PyMOL API
import numpy as np
import pandas as pd

def main():
    print(f"Running MDShift.py")

    with open("out.txt", "w", encoding="utf-8") as out:
        sys.stdout = out
        sys.stderr = out
        start_time = datetime.datetime.now()
        print(f"Job begins at: {start_time}\n")
        print(sys.argv)

        labels_file = None
        solvent_radius = 1.82  # default
        arg_count = len(sys.argv)

        # --- Argument parsing ---
        if arg_count == 1:
            print("No labels.txt or solvent radius specified.")
            print("Using N₂ as default buffer gas (solvent radius = 1.82 Å)\n")

        elif arg_count == 2:
            # Could be either a labels file or a radius
            arg = sys.argv[1]
            try:
                solvent_radius = float(arg)
                print(f"No labels.txt file specified. Using solvent radius = {solvent_radius} Å\n")
            except ValueError:
                labels_file = arg
                print(f"Labels file specified: {labels_file}")
                print("Using N₂ as default buffer gas (solvent radius = 1.82 Å)\n")

        elif arg_count == 3:
            labels_file = sys.argv[1]
            try:
                solvent_radius = float(sys.argv[2])
                print(f"Labels file specified: {labels_file}")
                print(f"Using custom solvent radius = {solvent_radius} Å\n")
            except ValueError:
                print("Error: Solvent radius must be a decimal number.")
                sys.exit(1)

        else:
            print("Usage: MDShift.py <labels.txt> <solvent radius>")
            sys.exit(1)

        # --- Process all .xyz files ---
        xyz_files = glob.glob("*.xyz")

        if not xyz_files:
            print("No .xyz files found in current directory.")
        else:
            print(f"Found {len(xyz_files)} .xyz files to process:")
            if labels_file:
                for xyz in xyz_files:
                    label(xyz, labels_file)
            else:
                print("No labels.txt file provided — skipping labeling step.\n")

            xyz_files = glob.glob("*.xyz") # Adds in labeled files if created
              
        for xyz in [f for f in xyz_files if not re.search(r"_label+ed", f)]:
            print(f"=== Processing {xyz} ===")
            extract_points(xyz, solvent_radius)
            print(f"=== Finished {xyz} ===\n")

        MoI_CoM_Calc()
        avg_distance()
        avg_spheres()
        avg_distance_each_dim()
        max_dist()
        dMoI_calc()

        end_time = datetime.datetime.now()
        print(f"Job ends at: {end_time}")
        print(f"Total run time: {end_time - start_time}")

############################################
######### PYMOL CONNOLLY CALCULATOR ########
############################################

def extract_points(xyz_file, solvent_radius):
    """
    Generate a Connolly surface for a given .xyz file in PyMOL and
    extract all surface coordinates to <basename>_points.txt.

    Parameters:
        xyz_file (str): input .xyz file. Base molecule with no isotopic labels
        solvent_radius (float): Probe radius for solvent surface generation.
    """
    print("--- PyMOL Surface Extraction Log ---\n")

    dir_path = os.path.dirname(os.path.realpath(__file__))
    basename = os.path.splitext(os.path.basename(xyz_file))[0]
    wrl_file = os.path.join(dir_path, f"{basename}.wrl")
    txt_file = os.path.join(dir_path, f"{basename}_points.txt")

    print(f"Processing {xyz_file} with solvent radius = {solvent_radius} Å")

    try:
        with pymol2.PyMOL() as pymol:
            cmd = pymol.cmd

            # Reset PyMOL
            cmd.reinitialize()

            # Load structure
            cmd.load(xyz_file, "mol")

            # Configure solvent surface parameters
            cmd.set("surface_solvent", 1)
            cmd.set("surface_quality", 1)
            cmd.set("solvent_radius", solvent_radius)
            cmd.set("dot_solvent", 1)

            # Generate and save solvent-accessible surface
            cmd.show("surface", "mol")
            cmd.save(wrl_file, "mol")
            print(f"  Surface saved as {wrl_file}")

            # --- Extract points from .wrl ---
            coords = []
            with open(wrl_file, "r") as f:
                data = f.read()

            matches = re.findall(r'point\s*\[([^\]]+)\]', data, re.S)

            for block in matches:
                for line in block.strip().splitlines():
                    line = line.strip().strip(',')
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            x, y, z = map(float, parts[:3])
                            coords.append((x, y, z))
                        except ValueError:
                            pass

            # Write surface points
            with open(txt_file, "w") as points_out:
                try:
                    if coords:
                        for x, y, z in coords:
                            points_out.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
                        print(f"  Wrote {len(coords)} surface points → {txt_file}")
                    else:
                        points_out.write("No points found\n")
                        print(f"  No points found in {wrl_file}")

                    #Removes the wrl file
                    if os.path.exists(wrl_file):
                        os.remove(wrl_file)

                except Exception as e:
                    print(f"Error in extract_points(): Error processing {xyz_file}: {e}")

    except Exception as e:
        print(f"\nError in extract_points(): PyMOL session error: {e}")

    print(f"\n--- Finished Processing Connolly Points for {basename} and saved as {basename}_points.txt ---")

############################################
##### END PYMOL CONNOLLY CALCULATOR ########
############################################

#############################################
########### COM MOI CALCULATOR ##############
#############################################

def read_xyz(filename):
    atoms = []
    x_vals = []
    y_vals = []
    z_vals = []

    with open(filename, 'r') as f:
        lines = f.readlines()

        # Skip atom count and comment line
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                atoms.append(parts[0])
                x_vals.append(float(parts[1]))
                y_vals.append(float(parts[2]))
                z_vals.append(float(parts[3]))

    return atoms, x_vals, y_vals, z_vals


def calc_CoM(atoms, x_vals, y_vals, z_vals):
    """
    Calculate the center of mass (COM).
    Returns (x_com, y_com, z_com).
    """
    total_mass = 0.0
    x_com = y_com = z_com = 0.0

    for atom, x, y, z in zip(atoms, x_vals, y_vals, z_vals):
        if atom not in mass:
            raise ValueError(f"Atom not found: {atom}")
        m = mass[atom]
        total_mass += m
        x_com += m * x
        y_com += m * y
        z_com += m * z

    x_com /= total_mass
    y_com /= total_mass
    z_com /= total_mass

    return x_com, y_com, z_com


def calc_inertia_tensor(atoms, x_vals, y_vals, z_vals):
    """
    Build and return the 3x3 inertia tensor as a list of lists (matrix).
    You can fill in the math for each element later.
    """

    # Initialize 3x3 matrix with zeros
    I = [[0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0]]

    for atom, x, y, z in zip(atoms, x_vals, y_vals, z_vals):
        if atom not in mass:
            raise ValueError(f"Atom not found: {atom}")
        m = mass[atom]

        I[0][0] += (y**2+z**2)*m  # placeholder for Ixx
        I[1][1] += (x**2+z**2)*m  # placeholder for Iyy
        I[2][2] += (x**2+y**2)*m  # placeholder for Izz
        I[0][1] += x*y*m  # placeholder for Ixy
        I[0][2] += y*z*m  # placeholder for Ixz
        I[1][2] += x*z*m  # placeholder for Iyz

    I[0][1] *= -1
    I[0][2] *= -1
    I[1][2] *= -1

    # Mirror symmetric elements
    I[1][0] = I[0][1]
    I[2][0] = I[0][2]
    I[2][1] = I[1][2]

    return I


# Atomic masses dictionary
mass = {"H": 1.007825, "D": 2.014102, "He": 4.002603, "He3": 3.016029, "Li": 7.016004, "Li6": 6.015122, "Be": 9.012182, "B": 10.012937, "B11": 11.009305, "C": 12.000000, "C13": 13.003355, "N": 14.003074, "N15": 15.000109, "O": 15.994915, "O18": 17.999160, "F": 18.998403, "Ne": 19.992440, "Ne21": 20.993847, "Na": 22.989769, "Mg": 23.985042, "Mg25": 24.985837, "Al": 26.981538, "Si": 27.976926, "Si29": 28.976495, "P": 30.973762, "S": 31.972071, "S34": 33.967867, "Cl": 34.968853, "Cl37": 36.965903, "Ar": 39.962383, "Ar36": 35.967545, "K": 38.963707, "K41": 40.961826, "Ca": 39.962591, "Ca44": 43.955482, "Sc": 44.955910, "Ti": 47.947947, "Ti46": 45.952629, "V": 50.943959, "V50": 49.947163, "Cr": 51.940506, "Cr53": 52.940648, "Mn": 54.938045, "Fe": 55.934938, "Fe57": 56.935394, "Co": 58.933194, "Ni": 57.935342, "Ni62": 61.928345, "Cu": 62.929598, "Cu65": 64.927790, "Zn": 63.929147, "Zn66": 65.926037,
"Ga": 68.925574, "Ga71": 70.924702, "Ge": 73.921178, "Ge72": 71.922076, "As": 74.921596, "Se": 79.916522, "Se77": 76.919914, "Br": 78.918338, "Br81": 80.916290, "Kr": 83.911507, "Kr84": 83.913425, "Rb": 84.911789, "Rb87": 86.909183, "Sr": 87.905613, "Sr86": 85.909260, "Y": 88.905848, "Zr": 89.904704, "Zr91": 90.905645, "Nb": 92.906378, "Mo": 97.905408, "Mo95": 94.905839, "Tc": 98.906254, "Ru": 101.904344, "Ru99": 98.905939, "Rh": 102.905504, "Pd": 105.903483, "Pd104": 103.904036, "Ag": 106.905093, "Ag107": 106.905097, "Cd": 113.903358, "Cd111": 110.904182, "In": 114.903878, "In113": 112.904061, "Sn": 119.902196, "Sn117": 116.902954, "Sb": 120.903816, "Sb121": 120.903812, "Te": 129.906223, "Te128": 127.904461, "I": 126.904468, "Xe": 131.904154, "Xe129": 128.904779, "Cs": 132.905451, "Ba": 137.905247, "Ba134": 133.904503, "La": 138.906348, "La138": 137.907107,
"Ce": 139.905439, "Ce136": 135.907172, "Pr": 140.907652, "Nd": 141.907723, "Nd142": 141.907723, "Pm": 144.912760, "Sm": 151.919739, "Sm149": 148.917183, "Eu": 152.921230, "Eu151": 150.919850, "Gd": 157.924103, "Gd155": 154.922630, "Tb": 158.925346, "Dy": 163.929175, "Dy161": 160.926933, "Ho": 164.930319, "Er": 165.930290, "Er167": 166.932048, "Tm": 168.934217, "Yb": 173.938866, "Yb171": 170.936332, "Lu": 174.940775, "Lu176": 175.942694, "Hf": 179.946549, "Hf176": 175.941409,
"Ta": 180.947992, "Ta180": 179.947464, "W": 183.950930, "W182": 181.948204, "Re": 186.955753, "Re185": 184.952955, "Os": 191.961479, "Os187": 186.955751, "Ir": 192.962924, "Ir191": 190.960591, "Pt": 194.964791, "Pt196": 195.964952, "Au": 196.966569, "Hg": 201.970643, "Hg199": 198.968262, "Tl": 204.974428, "Tl203": 202.972345, "Pb": 207.976652, "Pb204": 203.973043, "Bi": 208.980384, "Th": 232.038055, "Th230": 230.033133, "Pa": 231.035884, "Pa233": 233.040247,
"U": 238.050783, "U235": 235.043929}

def MoI_CoM_Calc():
    # Find and sort all .xyz files
    xyz_files = sorted(glob.glob("*.xyz"))

    if not xyz_files:
        raise FileNotFoundError("No .xyz files found in this directory.")

    # Prepare 5 rows that we'll fill horizontally in 4-column blocks (3 cols data + 1 blank)
    csv_rows = [[""] * (len(xyz_files) * 4) for _ in range(5)]

    for idx, filename in enumerate(xyz_files):
        print("=" * 60)
        print(f"Results for file: {filename}\n")
        
        atoms, x_vals, y_vals, z_vals = read_xyz(filename)

        try:
            com = calc_CoM(atoms, x_vals, y_vals, z_vals)
            I = calc_inertia_tensor(atoms, x_vals, y_vals, z_vals)

            print(f"Center of Mass (Å): ({com[0]:.6f}, {com[1]:.6f}, {com[2]:.6f})\n")

            # raw tensor
            print("Inertia Tensor:")
            for row in I:
                print(str(row))

            # diagonalization
            I_np = np.array(I)
            eigenvalues, eigenvectors = np.linalg.eigh(I_np)
            print("\nDiagonalized inertia tensor/principal moments (amu * \u00C5\u00B2):")
            print(str(np.diag(eigenvalues)) + "\n")
            print("Eigenvectors (principal axes):")
            print(str(eigenvectors) + "\n")

            # ===== CSV placement =====
            # First dataset uses columns A–C; next starts in E (i.e., add a blank spacer column).
            start_col = idx * 5  # 3 data cols + 1 label col in col A + 1 spacer -> next block begins at E
            # Ensure each CSV row is wide enough
            needed_width = start_col + 3  # we will write up to start_col+2 in some rows
            for r in range(len(csv_rows)):
                if len(csv_rows[r]) <= needed_width:
                    csv_rows[r].extend([""] * (needed_width + 1 - len(csv_rows[r])))

            # Row indices (0-based): 0->A1, 1->A2, 2->A3, 3->A4, 4->A5
            csv_rows[0][start_col] = filename                                      # A1: filename
            csv_rows[1][start_col] = "CoM (x,y,z)"                                 # A2: label
            csv_rows[2][start_col] = f"{com[0]:.6f}"                               # A3: x
            csv_rows[2][start_col + 1] = f"{com[1]:.6f}"                           # B3: y
            csv_rows[2][start_col + 2] = f"{com[2]:.6f}"                           # C3: z
            csv_rows[3][start_col] = "MoI (amu * \u00C5\u00B2)"                    # A4: label
            csv_rows[4][start_col] = f"{eigenvalues[0]:.6f}"                       # A5: Ix
            csv_rows[4][start_col + 1] = f"{eigenvalues[1]:.6f}"                   # B5: Iy
            csv_rows[4][start_col + 2] = f"{eigenvalues[2]:.6f}"                   # C5: Iz

        except ValueError as e:
            print(f"Error in file {filename}: {e}\n")

        print("=" * 60 + "\n")

    # Write the CSV/Excel-friendly file
    with open("out.csv", "w", newline="", encoding="utf-8") as fcsv:
        writer = csv.writer(fcsv)
        for row in csv_rows:
            writer.writerow(row)

##################################
######## END COM MOI CALC ########
##################################

def avg_distance():
    """
    Reads out.csv for all molecule CoM coordinates and calculates
    the average distance from each CoM to its corresponding *_points.txt file.
    Writes results into out.csv in cells D3, I3, N3, etc. with header
    'CoMtoSurf' in D2, I2, N2, etc.

    Filenames MUST contain base molecule name (i.e., molecule_points.txt, molecule.xyz, molecule_labeled_1.xyz, etc.)
    - molecule2.xyz and molecule_labeled_1.xyz will calculate relative to the surface of molecule2_points.txt
    """

    if not os.path.exists("out.csv"):
        raise FileNotFoundError("out.csv not found. Run MoI_CoM_Calc() first.")

    df = pd.read_csv("out.csv", header=None)
    print("\n=== Calculating mean distance from CoM to Connolly surface for isotopologues ===")
    # Loop through columns in 5-column blocks (A–E, F–J, etc.)
    for col in range(0, df.shape[1], 5):
        name = str(df.iat[0, col])
        if not isinstance(name, str) or not name.endswith(".xyz"):
            continue

        # Identify base molecule and points file
        base = re.sub(r"_labeled_\d+", "", os.path.splitext(name)[0])
        points_file = f"{base}_points.txt"

        if not os.path.exists(points_file):
            print(f"Warning: {points_file} not found for {name}, skipping.")
            continue

        try:
            # Extract CoM coordinates (row 2 → index 2)
            x = float(df.iat[2, col])
            y = float(df.iat[2, col + 1])
            z = float(df.iat[2, col + 2])
            com = np.array([x, y, z])

            # Load surface points
            coords = np.loadtxt(points_file)
            if coords.ndim == 1:
                coords = coords.reshape(1, -1)

            # Compute mean distance
            distances = np.linalg.norm(coords - com, axis=1)
            avg_dist = np.mean(distances)

            # Place header (row 1) and value (row 2)
            header_col = col + 3  # D, I, N, etc.
            if df.shape[1] <= header_col:
                # Extend if necessary
                new_cols = header_col - df.shape[1] + 1
                df = pd.concat([df, pd.DataFrame([[""] * new_cols for _ in range(df.shape[0])])], axis=1)

            df[header_col] = df[header_col].astype("object")
            df.iat[1, header_col] = "CoMtoSurf"
            df.iat[2, header_col] = f"{avg_dist:.6f}"

            print(f"{name:<25} | Mean distance to {points_file:<25} = {avg_dist:.6f} Å")

        except Exception as e:
            print(f"Error processing {name}: {e}")

    # Overwrite the original out.csv
    df.to_csv("out.csv", index=False, header=False)
    print("=== Updated out.csv with average distance from CoM to Connolly surface. ===\n")

def avg_distance_each_dim():
    """
    Calculates the average distance from the CoM to the Connolly surface
    in each individual dimension (x, y, z).
    Uses read_xyz() and calc_CoM().
    Prints average |Δx|, |Δy|, |Δz| for each isotopologue.
    """
    print("\n=== Calculating average |Δx|, |Δy|, |Δz| from CoM to Connolly surface ===")

    xyz_files = sorted(glob.glob("*.xyz"))
    if not xyz_files:
        print("No .xyz files found in current directory.")
        return

    for filename in xyz_files:
        try:
            atoms, x_vals, y_vals, z_vals = read_xyz(filename)
            com = np.array(calc_CoM(atoms, x_vals, y_vals, z_vals))
            base = re.sub(r"_labeled_\d+", "", os.path.splitext(filename)[0])
            points_file = f"{base}_points.txt"

            if not os.path.exists(points_file):
                print(f"Warning: {points_file} not found for {filename}, skipping.")
                continue

            coords = np.loadtxt(points_file)
            if coords.ndim == 1:
                coords = coords.reshape(1, -1)

            # Compute mean absolute difference per axis
            avg_dx = np.mean(np.abs(coords[:, 0] - com[0]))
            avg_dy = np.mean(np.abs(coords[:, 1] - com[1]))
            avg_dz = np.mean(np.abs(coords[:, 2] - com[2]))

            print(f"{filename:<25} | ⟨|Δx|⟩ = {avg_dx:.6f} Å   ⟨|Δy|⟩ = {avg_dy:.6f} Å   ⟨|Δz|⟩ = {avg_dz:.6f} Å")

        except Exception as e:
            print(f"Error processing {filename}: {e}")

    print("=== Completed average per-axis CoM-to-surface distance calculations. ===")

def avg_spheres():
    """
    Reads out.csv for all molecule CoM coordinates and calculates
    the average of (distance^2 * 4π) from each CoM to its corresponding *_points.txt file.
    Writes results only to out.txt (no CSV modification).

    Filenames MUST contain base molecule name (i.e., molecule_points.txt, molecule.xyz, molecule_labeled_1.xyz, etc.)
    """

    if not os.path.exists("out.csv"):
        raise FileNotFoundError("out.csv not found. Run MoI_CoM_Calc() first.")

    df = pd.read_csv("out.csv", header=None)
    print("\n=== Calculating ⟨distance² × 4π⟩ from CoM to Connolly surface for isotopologues ===")

    # Loop through columns in 5-column blocks (A–E, F–J, etc.)
    for col in range(0, df.shape[1], 5):
        name = str(df.iat[0, col])
        if not isinstance(name, str) or not name.endswith(".xyz"):
            continue

        # Identify base molecule and points file
        base = re.sub(r"_labeled_\d+", "", os.path.splitext(name)[0])
        points_file = f"{base}_points.txt"

        if not os.path.exists(points_file):
            print(f"Warning: {points_file} not found for {name}, skipping.")
            continue

        try:
            # Extract CoM coordinates (row 2 → index 2)
            x = float(df.iat[2, col])
            y = float(df.iat[2, col + 1])
            z = float(df.iat[2, col + 2])
            com = np.array([x, y, z])

            # Load surface points
            coords = np.loadtxt(points_file)
            if coords.ndim == 1:
                coords = coords.reshape(1, -1)

            # Compute ⟨r² × 4π⟩
            distances = np.linalg.norm(coords - com, axis=1)
            sphere_vals = (distances ** 2) * (4 * math.pi)
            avg_sphere = np.mean(sphere_vals)

            print(f"{name:<25} | ⟨r² × 4π⟩ to {points_file:<25} = {avg_sphere:.6f} Å²")

        except Exception as e:
            print(f"Error processing {name}: {e}")

    print("=== Completed ⟨r² × 4π⟩ calculations for all isotopologues. ===")

def max_dist():
    """
    Reads all xyz files (including isotopologues) and calculates
    the distance from the CoM to the furthest atom.
    Prints data to out.txt.
    """
    print("\n=== Calculating maximum CoM → atom distance for all molecules ===")

    xyz_files = sorted(glob.glob("*.xyz"))
    if not xyz_files:
        print("No .xyz files found in current directory.")
        return

    for filename in xyz_files:
        try:
            atoms, x_vals, y_vals, z_vals = read_xyz(filename)
            com = calc_CoM(atoms, x_vals, y_vals, z_vals)
            com = np.array(com)

            coords = np.column_stack((x_vals, y_vals, z_vals))
            distances = np.linalg.norm(coords - com, axis=1)
            max_distance = np.max(distances)
            max_atom_index = np.argmax(distances) + 1  # 1-based indexing

            print(f"{filename:<25} | Max CoM→atom distance = {max_distance:.6f} Å (atom {max_atom_index})")

        except Exception as e:
            print(f"Error processing {filename}: {e}")

    print("=== Completed maximum CoM → atom distance calculations. ===")


def dMoI_calc():
    """
    Calculates average relative change in MoI (dMoI) for each isotopologue
    relative to its base molecule (no '_labeled_' in name).
    Writes results into out.csv and prints to out.txt/console.

    Example:
        Base: Palmitate2.xyz → cols A–C
        Isotope: Palmitate2_labeled_1.xyz → cols F–H
        dMoI = (1/3)*[(F5/A5)+(G5/B5)+(H5/C5)]
        Stored in I4="dMoI", I5=<value>
    """
    import pandas as pd

    if not os.path.exists("out.csv"):
        raise FileNotFoundError("out.csv not found. Run MoI_CoM_Calc() first.")

    df = pd.read_csv("out.csv", header=None)
    ncols = df.shape[1]

    # Helper to get base molecule name
    def base_name(name):
        return re.sub(r"_labeled_\d+", "", name)

    # Map base molecules to their column indices
    molecule_blocks = {}
    for col in range(0, ncols, 5):
        name = str(df.iat[0, col])
        if isinstance(name, str) and name.endswith(".xyz"):
            molecule_blocks[name] = col

    # Group by base names
    grouped = {}
    for name in molecule_blocks:
        base = base_name(name)
        grouped.setdefault(base, []).append(name)

    print("\n=== Calculating dMoI for isotopologues ===")

    for base, group in grouped.items():
        base_col = None
        for g in group:
            if not re.search(r"_labeled_", g):
                base_col = molecule_blocks[g]
                break
        if base_col is None:
            print(f"No base molecule found for {base}, skipping group.")
            continue

        Ix_base = float(df.iat[4, base_col])
        Iy_base = float(df.iat[4, base_col + 1])
        Iz_base = float(df.iat[4, base_col + 2])

        for g in group:
            if g == base:
                continue
            if not re.search(r"_labeled_", g):
                continue

            iso_col = molecule_blocks[g]
            try:
                Ix = float(df.iat[4, iso_col])
                Iy = float(df.iat[4, iso_col + 1])
                Iz = float(df.iat[4, iso_col + 2])

                dMoI = (1/3) * ((Ix/Ix_base) + (Iy/Iy_base) + (Iz/Iz_base))

                header_col = iso_col + 3  # like I column in example
                df[header_col] = df[header_col].astype("object")
                df.iat[3, header_col] = "dMoI"
                df.iat[4, header_col] = f"{dMoI:.6f}"

                print(f"{g:<30} | dMoI vs {base:<25} = {dMoI:.6f}")

            except Exception as e:
                print(f"Error computing dMoI for {g}: {e}")

    # Write updated CSV
    df.to_csv("out.csv", index=False, header=False)
    print("=== Updated out.csv with dMoI values. ===")

def label(xyz_file, label_file):
    """
    Apply isotopic or atom labels from a label.txt file to an .xyz structure.

    Each pair of lines in label.txt defines a separate labeling operation that 
    generates a new <basename>_labeled_<n>.xyz file.

    Example label.txt:
        1 2 3
        D D D
        5
        13C
    → produces molecule_labeled_1.xyz (atoms 1–3 → D)
      and molecule_labeled_2.xyz (atom 5 → 13C)

    If only 1 element label, that atom will be applied to all atom numbers in above row.
    Example label.txt:
        1 2 3
        D
    → produces molecule_labeled_1.xyz (atoms 1—3 → D)
    """

    # Skip if no label file provided
    if label_file is None:
        print("No label file provided. Skipping labeling step.")
        return

    # Read XYZ file
    with open(xyz_file, "r") as f:
        lines = f.readlines()

    header1 = lines[0].strip()
    header2 = lines[1].strip()
    atoms_original = [line.strip().split() for line in lines[2:]]

    # Read label file
    with open(label_file, "r") as f:
        label_lines = [line.strip() for line in f if line.strip()]

    # Process each label pair separately
    for i in range(0, len(label_lines), 2):
        atom_indices = [int(x) - 1 for x in label_lines[i].split()]
        labels = label_lines[i + 1].split()

        # Allow single-label assignment to multiple atoms
        if len(labels) == 1 and len(atom_indices) > 1:
            labels = labels * len(atom_indices)
        elif len(atom_indices) != len(labels):
            print(f"Warning: mismatch in number of atoms and labels on pair {i//2 + 1}")
            continue

        # Create a copy of the atom list for this labeling variant
        atoms = [atom[:] for atom in atoms_original]

        # Apply the labeling
        for idx, new_label in zip(atom_indices, labels):
            if 0 <= idx < len(atoms):
                atoms[idx][0] = new_label
            else:
                print(f"Warning: atom index {idx+1} out of range for {xyz_file}")

        # Write a uniquely numbered labeled file
        basename = os.path.splitext(os.path.basename(xyz_file))[0]
        labeled_file = f"{basename}_labeled_{(i // 2) + 1}.xyz"

        with open(labeled_file, "w") as f:
            f.write(f"{header1}\n{header2}\n")
            for atom in atoms:
                f.write(" ".join(atom) + "\n")

        print(f"Labeled structure saved as {labeled_file}")
    print("\n")

if __name__ == "__main__":
    main()
