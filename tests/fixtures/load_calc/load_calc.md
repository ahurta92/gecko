# Information about each fixture

01  contains an mra hyperpolarizability calcualtion in one of the orginal formats containing BH2Cl. 
The output file is in outputs.json and for this type of directory, we also expect an input molecule file because
at this time we did not include the molecule in the outputs.json file.

01_mra-d04_n2 contains an mra hyperpolarizability calculation from the NLO project which the n2 molecule.  n2 is 
just a name we used to identify the molecule, it is not the actual n2 molecule.  

02_aug-cc-pVDZ_n2 contains an aug-cc-pVDZ basis set hyperpolarizability calculation from the NLO project which the n2 molecule.  

These three examples are a good example of what I want to be able to do with the load_calc function.  Ultimately
for basis set calculations, and mra calculations I want to identify the molecule or the calculation, the basis set of
the calculation, the method, hartree-fock or dft, and the hyperpolarizability results, ground-state energy, dipole moment, and polarizability and any other properties of interest.  If we can put this together in a common format that would be great.


03_mra-raman_h2o contains an mra hyperpolarizability calculation from the NLO project which the h2o molecule.  This calculation is using the latest version of madness.  It contains a geomtery optimization result followed by the raman results all in a single mad.h2o_gop.calc_info.json file.  This is a good example of how we want to store multiple calculation steps in a single calculation json.  The file also contains the input file and under the optimization block there is the final optimized geometry.  

05_dalton_raman_h2o contains a dalton raman calculation for h2o.  In this example we have the first mol file which is H2O_d-aug-cc-pV6Z.mol, this is the initial geometry used for the calculation.  Then we have optimization ouput in optimize_H2O_d-aug-cc-pV6Z.out and the parsed output opt_H2O_d-aug-cc-pV6Z.mol followed by the raman result in raman_opt_H2O_d-aug-cc-pV6Z.out


This dalton example is an example of a multi-step workflow.  

