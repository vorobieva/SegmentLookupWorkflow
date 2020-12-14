# SegmentLookupWorkflow
Rosetta and python scripts describing the workflow to graft loops on a protein with DirectSegmentLookup
These scripts were developed by Anastassia Vorobieva. The Rosetta Movers are available free of charge for the academic users from the Rosetta Molecular Modelling Suite. 

# Step1: Prepare the scaffolds. 
If your complex is symmetric, generate the symmetry difinition file. 

For a dimeric protein with C2 symmetry:
make_symmdef_file.pl -m NCS -p ../$file -a A -i B -r 10 > C2.symm

If your complex contains a small-molecule, generate the appropriate .param files (fa and cen) or use the default parameters in Rosetta. 

Relax the crystal structure while applying coordinate constraints (relax_coordinate_cst.xml).

If the goal is to diversify an existing loop, delete a few residues of that loop in the PDB. The DirectSegmentLookup mover will look for "TER" lines in the PDB to close the chain. Several segments can be connected at the same time. The search for segments matching the N- and C- terminii of the gap is exhaustive. Since the diversity of the returned solutions will strictly depend on the N- and C-terminii, it is a good idea to generate several scaffolds with different deletion end points as starting PDBs for the DirectSegmentLookup search. 

# Step2: Run DirectSegmentLookup.


