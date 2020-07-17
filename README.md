# Following particle exchange during Molecular Dynamics

This program works on a **.xyz trajectory file** with an **orthorhombic cell** to track proton transfers in water by counting for each oxygen center the number of hydrogens for which it is the closest oxygen, at each MD step. It can be easily modified to track exchange of different particles and/or different centers etc. 

The program outputs an .xyz file with the hydrogen coordination of each atom (0 for particles other than oxygen atoms) as an additional column. 

This additional column can be read by the Visual Molecular Dynamics (VMD) software so that the proton exchange can be dynamically visualized (see for instance [this video of a OH- moving through water](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09708-7/MediaObjects/41467_2019_9708_MOESM4_ESM.mov)). 



## Requirements 

1. .xyz trajectory format

2. gfortran (can be installed on Linux using `apt install gfortran` or ` brew install gfortran` on MacOS)

3. Download [XYZ_exchange_tracker.f90]() and [XYZ_exchange_tracker.sh]() and place both in the same directory

4. Execute `./XYZ_exchange_tracker.sh path/to/your_trajectory.xyz` 
   The output trajectory is named ` your_trajectory_exchange_track.xyz`
   Unless there is a CP2K input file named `md.in` in the directory of execution, the program will ask you for the orthorhombic cell parameters.

5. Visualize with VMD (see below)

   

## Visualization

1. Install the [readXYZUser](https://github.com/tonigi/vmd_extensions/blob/master/VMDextensions.tcl) plugin for [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) (see [here how to install VMD plugins](https://gist.github.com/tonigi/a9cfaf7642a7fbc13293))
2. Load your XYZ trajectory file in VMD
3. In the VMD console type `readXYZUser file_with_user_column.xyz `
4. Visualize with a new *representation*, with a *selection* such as for instance `within 1.8 of user 2 3` (selects H2O and H3O+ in the present case, as it would select atoms at 1.8 from oxygen centers with hydrogen coordinations of 2 or 3). Make sure you ticked *Update Selection Every Frame* and *Update Color Every Frame* under the *Trajectory* tab of the *Graphical Representations* window (under *Graphics > Representations*), otherwise VMD won't follow the exchange of particles. 

