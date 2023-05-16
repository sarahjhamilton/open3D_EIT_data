-----------------------------------------
1. TITLE: open3D_EIT_data
-----------------------------------------

-----------------------------------------
2. DESCRIPTION: 
-----------------------------------------
This is an open source collection of experimental data from electrical impedance tomography (EIT) collected on the ACT 5 EIT system on 5/28/2021 at University of Albany, Troy, NY, USA.

ACT 5 SYSTEM: 
O. Rajabi Shishvan, A. Abdelwahab, and G.J. Saulnier. ACT5 EIT system. In Proceedings of the 21st International Conference on Biomedical Applications of Electrical Impedance Tomography. Zenodo, June 2021. doi: 10.5281/zenodo.4635480. URL https://doi.org/10.5281/zenodo.4635480.


This is the data used in the following publication:

Sarah J Hamilton, P A Muller, D Isaacson, V Kolehmainen, J Newell, O Rajabi Shishvan, G Saulnier,
and J Toivanen, Fast absolute 3d CGO-based electrical impedance tomography on experimental tank data, Physiological Measurement 43 (2022), no. 12, 124001.

https://iopscience.iop.org/article/10.1088/1361-6579/aca26b/meta
DOI 10.1088/1361-6579/aca26b


A complete description of the data files included are given in the file Data-description.txt

Briefly, 3 experimental setups are considered for a 3D plexiglass box with 32 electrodes.  The first is saline only (no targets), the second is one target only, and the third has two targets.  

The purpose of this data is to test the robustness of reconstruction methods for EIT algorithms.  The data provided allows for testing of "absolute" EIT methods as well as "difference" EIT methods.  Making the data freely available allows more direct comparison of reconstruction algorithms on the same data, and for targets of known, measured, conductivity.

-----------------------------------------
3. How to cite this data
-----------------------------------------
When using this data, please give proper citation to the following two references:

** DATA + files **
Sarah J Hamilton, P A Muller, D Isaacson, V Kolehmainen, J Newell, O Rajabi Shishvan, G Saulnier,
and J Toivanen, Fast absolute 3d CGO-based electrical impedance tomography on experimental tank data, Physiological Measurement 43 (2022), no. 12, 124001.


** DEVICE **
O. Rajabi Shishvan, A. Abdelwahab, and G.J. Saulnier. ACT5 EIT system. In Proceedings of the 21st International Conference on Biomedical Applications of Electrical Impedance Tomography. Zenodo, June 2021. doi: 10.5281/zenodo.4635480. URL https://doi.org/10.5281/zenodo.4635480.


-----------------------------------------
4. How to run/use this data
-----------------------------------------

The EIT data provided are real-valued voltage and real-valued current files from data using the ACT 5 device.  They are in Matlab .mat files.  Details are given in the data-description.txt file as to the contents of each .mat file.  

Additionally, a finite element mesh for the measurement domain (3D box), compatible with EIDORS, is given.  The locations of the electrodes are given in the mesh as well as the excel file and locations.m Matlab file which plots the locations of the electrodes.

Furthermore, a sample EIDORS script (in Matlab) for simulating voltage data on the measurement domain is given for the specified mesh.

The data is provided to test reconstruction algorithms for electrical impedance tomography (EIT) on real, experimental, data with targets of known conductivity values (agar targets).

-----------------------------------------
5. Collaborators for data collection
-----------------------------------------
The following individuals were involved with the project:

Sarah Hamilton	 	Marquette University
Peter Muller 		Villanova University
David Isaacson		Rensselaer Polytechnic Institute
Jon Newell		Rensselaer Polytechnic Institute
Gary Saulnier		University of Albany
Omid Rajabi Shishvan	University of Albany


This research was supported by the National Institute of Biomedical Imaging and Bioengineering of the National Institutes of Health under award numbers R21EB028064 (SH and JN) 1R01EB026710-01A1 (GS, DI, JN, ORS, and the development of the ACT5 device). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. 
