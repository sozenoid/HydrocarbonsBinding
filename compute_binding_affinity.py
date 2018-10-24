################
################
################
################
####
#### This piece of code is meant to compute the binding affinity of a single of a molecule with CB[7] 
#### PRE: A pdbqt file containing the best n pose of the molecule with CB[7] is given as generated by autodock vina (fed a pdbqt file produced using prepare_ligand4.py from a rdkit pdb), this pdbqt file is used to reconstruct the conformations of the guest molecules
####	  In their free state as well as docked
####      The binding affinity computed AMBER from AmberTools. The complex are generated again by combining the CB_candidate.pdbqt coordinates with the coordinates of the molecule poses
####	  
####
#### POST: The binding affinity for all the poses are computed and the added to an output log file
################
################
################
################


class guestMolecule:
    def __init__(self, pdbqtblock, sdfref):
        """
		PRE : pdbqtblock is a valid pdbqt file containing n conformations for the same guest molecule (usually 9 or less).
			  The rdkit module should be available. The order of the atoms in sdfref and pdbqt block should be the same.
			  SDFREF contains the information about hydrogen placements and bonding, this is needed as the passage through pdbqt format breaks connectivity records and hydrogen positions
		POST : an object guestMolecule is created that contains 
		"""
        self.workPath = "/home/macenrola/Desktop/"
        self.pathtocbpdbfile = "/home/macenrola/Documents/Thesis/ScreeningManuscriptFinalData/HydrocarbonsBindingPython/CB_candidate_data/CB_candidate.pdb"
        self.pdbstring = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:>6s}"
        self.molName = sdfref.split('/')[-1][:-4]
        # self.baseConf = [Chem.AddHs(x, addCoords=True) for x in Chem.SDMolSupplier(sdfref)]
        print 'about to get the base comps'
        self.baseConf = self.generate_mols_withH_withCO_from_ALL_pdbqt(pdbqtblock, sdfref)
        print 'ending to get base comps'
        self.showmollist(self.baseConf)
        self.dockedConf = self.get_docked_conformations(self.baseConf)
        self.produce_formatted_PDBs(self.baseConf, self.dockedConf)

    def showmollist(self, listOfMols):
        """
		:param listOfMols: a list of mols in rdkit format
		:return: prints the mols as sdf blocks
		"""
        if len(listOfMols) == []:
            print 'the list is empty'
        elif len(listOfMols) == 1:
            print Chem.MolToMolBlock(listOfMols)
        else:
            for mol in listOfMols:
                print "smth should show"
                print Chem.MolToMolBlock(mol)

    def get_docked_conformations(self, listOfMols):
        """
		PRE: Takes in a list of valid rdkit mol formats
		POST: Returns another list of those molecules merged with the pdb coordinates of CB[7]
		"""
        dockedConfs = []
        CB_host = Chem.MolFromPDBFile(self.pathtocbpdbfile, removeHs=False)
        for mol in listOfMols:
            dockedConfs.append(Chem.CombineMols(mol, CB_host))
        self.write_mol_to_pdb(CB_host, fout="{}CB_candidate_formatted.pdb".format(self.workPath),
                              residueName=('CB7', 'GST'), residueNumber=(2, 1))
        return dockedConfs

    def get_sphere_circle(self, listOfDockedMols):
        """
		PRE: Takes in a list of complexes
		POST: Returns their enclosing sphere diameter
		"""
        pass

    def produce_formatted_PDBs(self, listOfMols, listOfDockedMols):
        """
		PRE: Takes in rdkit mols of guests and their docked versions
		PST: Produces one valid pdb file for each molecule and docked complex 
		"""
        for lst in zip(["guests", "complex"], [listOfMols, listOfDockedMols]):
            for i, mol in enumerate(lst[1]):
                self.write_mol_to_pdb(mol, fout="{}{}_{}Pose{}.pdb".format(self.workPath, self.molName, lst[0], i))

    def remove_CONNECT_LINES(self, fname):
        """
		PRE: fname contains a valid PDB file of a molecule with 3D coordinates, MOLECULES NEED TO BE NAMED
		POST: The lines containing CONNECT and MASTER are removed from the original file, the original file IS MODIFIED
		"""
        with open(fname, 'rb') as r:
            lines = r.readlines()
        with open(fname, 'wb') as w:
            w.writelines([x for x in lines if 'CONECT' not in x if 'MASTER' not in x][1:])

    def write_mol_to_pdb(self, molin, fout, converge=False, residueName=('GST', 'CB7'), residueNumber=(1, 2)):
        """
		PRE : Takes in the absolute path to a SDF file, if a single molecule is converted only the first element of both residueName and residueNumber is used, both in the case of a complex
		POST: Produces a PDB file by converting the original molecule from the SDF file, DOES NOT optimize the conformation by default
		"""
        mol = molin
        flavour = 28
        if converge:
            # converge_molecule(mol)
            print "Can't converge molecule, this function isn't defined"
        isComplex = (len(Chem.GetMolFrags(mol)) != 1)

        print '{} contains {} Fragments'.format("current mol", len(Chem.GetMolFrags(mol)))

        if not isComplex:
            atm_dic = {}
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in atm_dic:
                    atm_dic[atom.GetSymbol()] = 1
                else:
                    atm_dic[atom.GetSymbol()] += 1

                atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} {}'.format(
                    atom.GetSymbol() + str(atm_dic[atom.GetSymbol()]) + ' ' * int(atm_dic[atom.GetSymbol()] < 10),
                    residueName[0])))
                atom.GetMonomerInfo().SetResidueNumber(residueNumber[0])
            Chem.MolToPDBFile(mol, fout, flavor=flavour)
            self.fix_PDB_spacing(fout, hostname=residueName[1])
            return
        else:
            comp = Chem.GetMolFrags(mol, asMols=True)
            if comp[0].GetNumAtoms() < comp[1].GetNumAtoms():
                guest, host = comp[0], comp[1]
            else:
                guest, host = comp[1], comp[0]
            for res in zip((guest, host), [0, 1]):
                atm_dic = {}
                for atom in res[0].GetAtoms():
                    if atom.GetSymbol() not in atm_dic:
                        atm_dic[atom.GetSymbol()] = 1
                    else:
                        atm_dic[atom.GetSymbol()] += 1

                    atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} {}'.format(
                        atom.GetSymbol() + str(atm_dic[atom.GetSymbol()]) + ' ' * int(atm_dic[atom.GetSymbol()] < 10),
                        residueName[res[1]])))
                    atom.GetMonomerInfo().SetResidueNumber(residueNumber[res[1]])

            comp_annotated = Chem.CombineMols(guest, host)
            comp_annotated.SetProp('_Name', mol.GetProp('_Name'))
            Chem.MolToPDBFile(comp_annotated, fout, flavor=flavour)
            self.remove_CONNECT_LINES(fout)
            self.fix_PDB_spacing(fout, hostname=residueName[1])

    def fix_PDB_spacing(self, fname, hostname='CB7', strandname='A'):
        """
		PRE: The PDB files is formated by rdkit MolToPDBFile flavour 28 without any MASTER, TER, CONECT or Charge flags, 
		POST: The file IS MODIFIED to be formated with the adequate amount of blank space to be read by AMBER, This method destroys the original file
		"""
        raw_spc = [7, 2, 4, 4, 4, 12, 8, 8, 6, 6, 12, 0]
        new_lines = []
        cb_yet = False
        with open(fname, 'rb') as r:
            for line in r:
                if hostname in line and not cb_yet:
                    cb_yet = True
                    new_lines.append('TER')
                    strandname = 'B'
                if 'ATOM' in line:
                    line_content = line.split()
                    line_content.insert(4, strandname)
                    # print line_content
                    ls = [len(x) for x in line_content]
                    actual_spacing = [raw_spc[0] - ls[1],  # after ATOM
                                      raw_spc[1],  # after ATOM# (11)
                                      raw_spc[2] - ls[2],  # after ATOM NAME (H41)
                                      raw_spc[3] - ls[3],  # after RESIDUE NAME (CUC)
                                      raw_spc[4] - ls[4],  # after CHAIN ID (A)
                                      raw_spc[5] - ls[6],  # after RESIDUE NUMBER (1)
                                      raw_spc[6] - ls[7],  # after X CART COORDINATE (6.171)
                                      raw_spc[7] - ls[8],  # after Y CART COORDINATE (3.377)
                                      raw_spc[8] - ls[9],  # after Z CART COORDINATE (21.096)
                                      raw_spc[9] - ls[10],  # after enigmatic number (1.00)
                                      raw_spc[10] - ls[11],  # after partial charge (0.00)
                                      raw_spc[11],  # after ATOMIC SYMBOL (N)
                                      ]
                    new_lines.append(''.join([x[0] + ' ' * x[1] for x in zip(line_content, actual_spacing)]))
                # print new_lines[-1]

                else:
                    new_lines.append(line.strip())
                # print new_lines[-1]
        with open(fname, 'wb') as w:
            w.writelines('\n'.join(new_lines) + '\n')

    def generate_mols_withH_withCO_from_ALL_pdbqt(self, ALL_pdbqt_file, original_sdf_file):
        """
		PRE: Takes in one sdf file and a pdbqt file, one sdf file with correct bonding information and arbitrary coordinates and another one containing all ligand docked conformations in pdbqt format
			 where the coordinates are desireable but there is no hydrogen or connect records
		POST: Produces a list of mols with connect records taken from the pdb file that has both hydrogens and connect records and the coordinates from the
			  pdbqt file. The mols are under the rdkit format as python objects

		comment: this method is necessary as autock vina breaks down the connect records by reverting to raw cartesian coordinates

		"""

        def get_coordinates_block_from_pdbqt_file():
            """
			PRE: Takes in a pdbqt file as a piece of text
			POST: returns a list of list where each sublist contains the coordinate lines from a pdbqt conformation (in practice a fraction of the ALL file)
				  each conformation is encoded as a list of strings. All lists should have the same number of elements
			"""
            list_of_conf_as_string = []
            print ALL_pdbqt_file
            with open(ALL_pdbqt_file, 'rb') as r:
                for line in r:
                    temp = []
                    if 'MODEL' in line:
                        for i in range(100):
                            if 'HETATM' in line:
                                temp.append(line[0:77])  # (77) Necessary to trim the final A in several pdbqt files
                                line = r.next().strip()
                            elif "ENDMDL" in line:
                                list_of_conf_as_string.append(temp)
                                print temp
                                temp = []
                                break
                            else:
                                line = r.next()

            return list_of_conf_as_string

        def sort_coordinate_list_by_carbon_number(string_of_coordinates_line):
            """
			:pre: takes a list of strings that represent the atomic coordinates in a pdb format
			'HETATM    5  C5  UNL     1       2.952   0.220  -0.069  1.00  0.00           C  '
			:return: return those very same lines aligned by carbon number, in this case C5
			"""
            sortingarray = []
            list_sorted_strings = []
            for coordstring in string_of_coordinates_line:
                sortingarray.append((int(coordstring.split()[2][1:]), coordstring))
            for i, els in enumerate(sorted(sortingarray)):
                temp = els[1].split()
                temp.insert(3, "")
                temp.insert(5, "")
                temp.insert(7, "")
                temp.insert(-1, "")
                if i > 8:
                    temp[2] = " " + temp[2]
                # print self.pdbstring.format(temp[0], int(temp[1]), temp[2], temp[3], temp[4],
                # 							temp[5], int(temp[6]), temp[7], float(temp[8]), float(temp[9]),
                # 							float(temp[10]), float(temp[11]), float(temp[12]), (temp[13]), (temp[14]))
                list_sorted_strings.append(self.pdbstring.format(temp[0], i + 1, temp[2], temp[3], temp[4],
                                                                 temp[5], int(temp[6]), temp[7], float(temp[8]),
                                                                 float(temp[9]),
                                                                 float(temp[10]), float(temp[11]), float(temp[12]),
                                                                 (temp[13]), (temp[14])))
            return list_sorted_strings

        # return [self.pdbstring.format(x[1].split().insert(5,'').insert(14,'')) for x in sorted(sortingarray)]

        def sort_list_of_list_of_coordinate_string(list_of_list_coords):
            """
			:parms: takes in a list of list of strings where each sublist corresponds to a conformation of a given molecule while the strings are the atomic coordiantes in pdb format
			:return: the very same list of list where the individual strings are sorted within each sublist according to the carbon number (sort_coordinate_list_by_carbon_number)
			"""
            temp = []
            for confs in list_of_list_coords:
                temp.append(sort_coordinate_list_by_carbon_number(confs))
            return temp

        def generate_pdb_noH_withConnect():
            """
			PRE: Takes in a sdf file with hydrogens and bonds connect
			POST: Generates a pdb block as a list of strings without the hydrogens yet still the hydrogens for the remaining carbons
			"""
            sdf_H_CO = Chem.SDMolSupplier(original_sdf_file, removeHs=True)[0]
            return [x.strip() for x in Chem.MolToPDBBlock(sdf_H_CO).split("\n")]

        conf_list = get_coordinates_block_from_pdbqt_file()
        skeleton = generate_pdb_noH_withConnect()
        conf_list_sorted = sort_list_of_list_of_coordinate_string(conf_list)
        # print '\n'.join(skeleton)
        list_of_mols = []
        for conf in conf_list_sorted:
            pdb_conf = [skeleton[0]] + conf + skeleton[len(conf) + 1:-1]
            pdb_conf = '\n'.join(pdb_conf)
            print pdb_conf
            list_of_mols.append(Chem.AddHs(Chem.MolFromPDBBlock(pdb_conf), addCoords=True))
        return list_of_mols


def convert_sdf_to_pdb(sdf_h_and_connect, pdb_noh_and_connect=None):
    """
	PRE: Takes in a pdb file with hydrogens and bonds connect
	POST: Generates a pdb file without the hydrogens yet still the hydrogens for the remaining carbons
	"""
    if pdb_noh_and_connect == None:
        pdb_noh_and_connect = sdf_h_and_connect[:-4] + "_noH.pdb"

    sdf_H_CO = Chem.MolFromMolFile(sdf_h_and_connect, removeHs=True)
    Chem.MolToPDBFile(sdf_H_CO, pdb_noh_and_connect, flavor=28)


if __name__ == "__main__":
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import glob

    guestMolecule('/home/macenrola/Desktop/7734xcw_OUT_GUEST28_noH.pdbqt-ALL.pdbqt',
                  '/home/macenrola/Desktop/7734xcw_OUT_GUEST28.sdf')
# flist = glob.glob('/home/macenrola/Thesis/ScreeningManuscriptFinalData/Generation3Ddata/pdb_2_pdbqt_w_prepare_ligand4/*.sdf')
# for f in flist:
# 	convert_sdf_to_pdb(f)
# testGuest = guestMolecule("/home/macenrola/Thesis/ScreeningManuscriptFinalData/HydrocarbonsBindingPython/tests/docking_w_connect_records/8473xeb_OUT_GUEST166_adtools-ALL.sdf")
# generate_pdb_noH_withConnect("/home/macenrola/Thesis/ScreeningManuscriptFinalData/HydrocarbonsBindingPython/tests/8473xeb_OUT_GUEST166.sdf")
# listofmols = generate_mols_withH_withCO_from_ALL_pdbqt('/home/macenrola/Thesis/ScreeningManuscriptFinalData/HydrocarbonsBindingPython/tests/dock_w_prepare_ligand4/931xae_OUT_GUEST15.sdf',
# 	'/home/macenrola/Thesis/ScreeningManuscriptFinalData/HydrocarbonsBindingPython/tests/dock_w_prepare_ligand4/931xae_OUT_GUEST15_noH_out.pdbqt')
# for i, m in enumerate(listofmols):
# 	print Chem.MolToMolFile(m, "/home/macenrola/Desktop/testPose{}.sdf".format(i))
# convert_sdf_to_pdb("/home/macenrola/Desktop/602809xziop_OUT_GUEST22.sdf")
