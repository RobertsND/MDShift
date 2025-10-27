# MDShift.py
**Purpose:** The purpose of this program is to calculate the changes in mass distribution upon isotopic labeling in order to explain/predict isotopic shifts in high-resolution ion mobility separations. For relevant literature visit nagylab.com.  
**Authors:** Noah D. Roberts and Gabe Nagy  
**Contact:** noah.roberts@utah.edu, gabe.nagy@utah.edu  
**Citation**: Roberts, N.D., Ramos, K., Jijieshvili, G., Armentrout, P.B., Nagy, G. Isotopic dimethylation coupled with high-resolution ion mobility separations enables site-specific isotopic shifts in peptide isomers. 2026, *Anal. Chem.*  

Notes: While it can be run on any operating system, this readme file contains written instructions for running MDShift on Windows. This program has not been written for parallelization, but it requires very little computational power and can be run locally; runtime may be on the order of hundreds of milliseconds per molecule.


# Table of Contents
1. [Requirements](#req)
2. [Implementation](#imp)
3. [Troubleshooting](#trou)
4. [Background](#background)


# 1. [Requirements](#req)
**PyMOL and Python**
PyMOL is required for surface generation. A free version for evaluation can be downloaded at python.org. The easiest way to implement MDShift is to add the PyMOL embedded Python API to your system's environment variables. This way, Python will be run within PyMOL and have access to all of PyMOL's libraries. See [Troubleshooting](#trou) for more information.

You will need to install numpy and pandas into your PyMOL Python environment.

` pip install numpy pandas `

# 2. [Implementation](#imp)

Run the script from the directory containing your `.xyz` and optional `labels.txt` files:

` python MDShift_v2.py <labels.txt> <solvent_radius> `

This will automatically apply `labels.txt` to all `.xyz` files in that directory and execute calcuations on each resulting isotopologue. Both arguments are optional and can be provided in any order. For example, if only a single argument is given and it is a decimal, no labels file will be used and the decimal will be stored in `solvent_radius`.

**Arguments**
* ``<labels.txt>`` (optional)
Enter the name of your `.txt` file which specifies the isotopic substitutions to be made. The file should be setup so that each pair of lines corresponds to a unique isotopologue to be generated. For example, in the input file below, three isotopologues would be generated:
	- (1) substitution with deuterium at atoms 1 through 6
	- (2) substitution with with carbon-13 at atoms 7 through 9
	- (3) substitution with deuterium at atoms 1 through 6 *and* carbon-14 at atoms 7 through 9
		`labels.txt`
		```
		1 2 3 4 5 6
		D
		7 8 9
		C13 C13 C13
		1 2 3 4 5 6 7 8 9
		D D D D D D C13 C13 C13
		```
		IMPORTANT NOTES:
		- If only one atom label is specified, that atom will be applied to every atom number above. Alternatively, each atom number can have a corresponding atom label. See `Examples` for more.
		- If you elect not to use a a labels file, each `.xyz` file will get a surface calculation and all rotational information will be calculated. However, in order to use the capability of the program to automatically calculate changes relative to the unlabeled isotopologue, it is recommended to use a labels file or manually create them. If you elect to manually create isotopically labeled `.xyz` files, make sure to save them as `molecule_labeled_1.xyz`, `molecule_labeled_2.xyz`, etc. Please see `Examples` for more.
* ``<solvent_radius>`` (optional)
	This specifies the solvent/probe radius in Å used for surface generation. If no value is entered, the default collision gas for surface generation is N<sub>2</sub> which uses a solvent radius of 1.82 Å. Any unitless decimal can be entered here for alternate collision gases.

**Output**
| File | Contents|
|:---:|:---:|
| `out.txt` | Descriptive output of all computed values.
| `out.csv` | Outputs center of mass cooridnates and principal moments of inertia.
| `*_labeled_n.xyz` | Each molecule (*) that has been isotopically labeled gets a unique filename. **DO NOT** rename these files as they are automatically linked to base molecules.
| `*_points.txt` | Connolly surface points of each unlabeled molecule (*). **DO NOT** rename these files as they are automatically linked to base molecules.|

A temporary `.wrl` file is generated during surface generation and deleted when the program is terminated. If you would like to save this file, delete or comment the following line in the source code:

`os.remove(wrl_file)`

**Optional implementation**
Included is an optional script `Run_MDShift.ps1` for Windows PowerShell which allows MDShift to be run from any directory. Replace `$ScriptPath` with the path to MDShift.py and implement with the path to where your `.xyz` and optional `.txt` files are:

`.\Run_MDShift.ps1 "C:\Users\Henry Eyring\data" <labels.txt> <solvent_radius>`

To use the current working directory, use `pwd` `cwd` `.` or `./` in place of the command line path.

# 3. [Troubleshooting](#trou)
**ModuleNotFoundError: No module named 'pymol'**
- First, make sure PyMOL is installed and added to your PATH.
- In order for MDShift to use the PyMOL API, the script needs to find PyMOL's version of Python first before your regular version of Python.
- To test this, run the following command in PowerShell: `python -c "import pymol"`
- If this returns `ModuleNotFoundError`, your system is likely using python outside of PyMOL. To fix this, simply navigate to your PATH environment variables and move all PyMOL paths to the top.
Check that command again and if no errors are returned, then you are good to go.
- If you want to use Python in the future outside of PyMOL, simply use `py` instead of `python` when running .py scripts from the command line.

# 4. [Background](#background)
Isotopic shifts in ion mobility separations have been difficult to model due to very subtle differences in arrival time of isotopically labeled species. This work attempts to explain isotopic shifts through calculating the changes in center of mass (CoM) and moments of inertia (MoI) upon isotopic substitution.

**CoM**  
It is insufficient to simply describe the change in CoM as only a distance. This is because if the CoM shifts more towards the molecule's geometric center, we would anticipate a reduction in cross section. Thus, it becomes necessary to describe the change in CoM in terms of both its magnitude *and* direction. Quantifying the CoM shift in a particular direction is trivial, but it becomes difficult to decide on a singular quantity to describe the shift in irregular 3D objects. Thus, we elected to describe the shift as a change in the the average distance from the CoM to the ion-neutral collision surface (i.e., the Connolly surface or the solvent accessible surface, SAS). We begin with the base, unlabeled isotopologue by calculating the average distance from its CoM to the SAS and call this value r̄<sub>p<sub>light</sub></sub>. Then, upon isotopic substitution, the CoM is recalculated and r̄<sub>p<sub>heavy</sub></sub> is calculated relative to the original SAS. Each r̄<sub>p</sub> value is reported in `out.csv` as `CoMtoSurf`. `Examples\Palmitate` contains `Palmitate.xyz` and `label.txt` which labels the Palmitate anion with two variants of isotopic substitution.

Another variation of this calculation can be done where each distance r<sub>p</sub> is represented as a sphere before averaging by squaring each value of r<sub>p</sub> and multiplying by 4π. This done with `avg_spheres()` and the output is shown in `out.txt` under `Calculating ⟨distance² × 4π⟩ from CoM to Connolly surface for isotopologues`.

Another function, `avg_distance_each_dim()`, calculates the average distance from the center of mass to the solvent-accessible/Connolly surface along each Cartesian axis. These calculations are performed in the molecule’s native coordinate frame and do not apply any transformation to align with the principal moments of inertia. Therefore, the resulting values should be interpreted with caution.

Finally, `max_dist()` calculates the distance from the CoM to the furthest atom. This strategy attempts to explain the shift in CoM assuming that the collision surface is a perfect sphere made by the outermost point of the molecule. This strategy does not account for the shape of the molecule (e.g., no difference between oblate and prolate spheroids) which may be an oversight and should also be interpreted with caution.

**MoI**  
In order to obtain a single quantity which describes the change in the rotational inertia, we elected to calculate dMoI which is an average change in the each moment of inertia relative to the starting isotopologue.
$$
 dMoI =\frac{1}{3}\left (\frac{I_{x, heavy}}{I_{x, light}}+\frac{I_{y, heavy}}{I_{y, light}}+\frac{I_{z, heavy}}{I_{z, light}} \right )
$$
 Moments of inertia are calculated by building the inertia tensor in the molecule's coordinate system and diagonalizing using numpy.linalg.eigh(). Each `.xyz` file is processed and `out.txt` contains the initial inertia tensor, the diagonalized inertia tensor, and the eigenvectors. The columns of the eigenvectors correspond to the principal axes in the original coordinate system. The diagonal elements of the diagonalized inertia tensor correspond to the principal moments in each of the three principal axes.




