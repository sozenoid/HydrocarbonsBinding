# Binding affinity of rigid hydrocarbons with CB[7]

This repository is meant to hold the codes used to process hydrocarbon molecules while their binding affinity with CB[7] is assessed.
The code mostly serves as a binding between AutoDock Vina, the RDKIT, Open Babel and native python libraries. 
Hugues Lambert

# Convert from SDF to PDB then PDBQT
Using Open Babel convert the SDFs to PDBQT

```for f in *GUEST*sdf; do obabel $f -O $f.pdbqt; echo $f; done```


Then using Autodock Vina do the docking 
```for f in ALL_RIGID/*.pdbqt; do ./bin/vina --config noligconf.txt --ligand=$f --out=$f-ALL.pdbqt > $f-report-vina; done```
