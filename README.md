# waveaxe

##How to Run the Program?
-	Run: FCHFX_WB2_2Crack_dynamic_MultiplePulse.m

##How to model inputs?
   * Model Inputs are prescribed in ConfigFile_Twocrack.m
   * The config file prescribes the input/output path (IOPATH) to the folder containing the model mesh and model output
   * The config file prescribes the name of the file containing the mesh (MeshFileName) 
		% mesh file name and location.
   * Code reads a gmsh file mesh files, e.g., 1WB9in_Q4S_v2.msh
   * This file must be in the folder given by the IOPATH.

##Before you Run the Program
1.	Create a folder for IO.
2.	Copy the mesh file (1WB9in_Q4S_v2.msh) into this folder
3.	Modify the IOPATH in the config file.