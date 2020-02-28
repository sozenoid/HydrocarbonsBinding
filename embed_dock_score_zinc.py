def split_main_file_in_block_of_1000(fname):
	"""
	PRE: Will split the main file in blocks of 1000 easier to process
	"""
	with open(fname, "r") as r:
		reader = csv.reader(r, delimiter=',')
		for i,row in enumerate(reader):
			if i%1000==0:
				aname = "/home/macenrola/HydrocarbonsBinding/zinc_pieces/{}-{}-zinc-pieces".format(i, i+1000)
			with open(aname, "a") as w:
				writer = csv.writer(w, delimiter=',')
				writer.writerow(row)
			
# =============================================================================
# 			if i%3000==0 and i>1:
# 				break
# =============================================================================

def embed_mol(a,b,c,d,e):
	"""
	PRE: Takes a smiles
	POST: Returns a rdkit mol
	"""
	smiles = a
	name = e
	mol = Chem.MolFromSmiles(smiles)
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	AllChem.MMFFOptimizeMolecule(mol)
	mol.SetProp("_Name", name)
	Chem.MolToMolFile(mol, "/home/macenrola/Desktop/MACHINE_LEARNING/Zinc_Dock/{}.sdf".format(name))
	return mol

def dock_with_ledock(mol, i, ledock_path):
	"""
	PRE: Takes in a rdkit mol
	POST: Will dock it and add the result to a list
	"""
	
	return
	
def process_a_batch(fname, ledock_path):
	"""
	PRE: Takes in a batch
	POST: Will go through it and process it
	
	"""
	recap_dic = {}
	recap_list = []
	with open(fname, "r") as r:
		reader = csv.reader(r, delimiter=',')
		for i,row in enumerate(reader):
			namestring = "{0:07d}".format(i)
			print(namestring, row)
			recap_dic[namestring] = row
			row.append(namestring)
			recap_list.append(row)
	return recap_list

def multiproc_embed_mol():
	ledock_path = "/home/macenrola/Desktop/MACHINE_LEARNING/ledock/lepro_linux_x86"


	#split_main_file_in_block_of_1000("/home/macenrola/HydrocarbonsBinding/250k_rndm_zinc_drugs_clean_3.csv")
	recap_list = process_a_batch("/home/macenrola/Desktop/MACHINE_LEARNING/Zinc_Dock/250k_rndm_zinc_drugs_clean_3.csv", ledock_path)
	print(recap_list)
	pickle.dump( recap_list, open( "/home/macenrola/Desktop/MACHINE_LEARNING/Zinc_Dock/recap_list_250k", "wb" ) )
	with multiprocessing.Pool(processes=16) as pool:
		results = pool.starmap(embed_mol, recap_list)


def associate_a_dok_w_pdb(receptor_pdb, dokfile):
	"""
	PRE  : Takes in a dok file
	POST : Associates it with the original receptor pdb file
	"""
	receptor = Chem.MolFromPDBFile(receptor_pdb, removeHs=False)
	print(Chem.MolToMolBlock(receptor))
	with open(dokfile, "rb") as r:
		dokcontent = r.readlines()
	indices = [i for i, string in enumerate(dokcontent) if b"REMARK" in string ]# finds occurences of REMARK
	for j, k in enumerate(indices[1:-1]): # discard the first occurence as it is only the processing time
		PDBBlock = (b''.join(dokcontent[indices[1:][j]:indices[1:][j+1]])) # Gets each of the poses one by one based on the location of the REMARK token 
		temp = tempfile.NamedTemporaryFile() # creates a temporary file for each pose because the MolFromPDBFile works better than MolFromPDBBlock method
		temp.write(PDBBlock) # write it to the named temporary file 
		temp.seek(0) # needs to rewind to be able to read it from without
		#print(temp.read())
		mol = Chem.MolFromPDBFile(temp.name, removeHs=False) # Build the molecule from the PDB file and creates the connectivity records
		#print(Chem.MolToMolBlock(mol))
		temp.close()

		### Now rebuilds the complex based on the original receptor PDB and the reconstruced PDB obtained from the dok file
		comp = Chem.CombineMols(receptor, mol)
		
		print(Chem.MolToMolBlock(comp))
		Chem.MolToMolFile(comp, dokfile+"LOL{}.sdf".format(j))
	
		
if __name__ == "__main__":
	import csv
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import multiprocessing
	import pickle
	import tempfile

	associate_a_dok_w_pdb("/home/macenrola/Desktop/MACHINE_LEARNING/Zinc_Dock/CB7_pure_het.pdb", "/home/macenrola/Desktop/MACHINE_LEARNING/Zinc_Dock/0026211.sdf.dok")