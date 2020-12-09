clear all
clc
format long

global SMesh CMesh Domain Material Control ConfigFileName

disp(' X-FEM HYDRAULIC FRACTURE SIMULATOR ')
tic

disp([num2str(toc),': Loading Mesh and Config file...']);
ConfigFileName = 'ConfigFile_Singlecrack';
BuildMesh_v13