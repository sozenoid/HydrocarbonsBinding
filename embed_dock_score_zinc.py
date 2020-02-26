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

def embed_mol(smiles):
	"""
	PRE: Takes a smiles
	POST: Returns a rdkit mol
	"""
	mol = Chem.MolFromSmiles(smiles)
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	AllChem.MMFFOptimizeMolecule(mol)
	return mol

def dock_with_ledock(mol, ledock_path):
	"""
	PRE: Takes in a rdkit mol
	POST: Will dock it and add the result to a list
	"""
	
	return
	
def process_a_batch(fname):
	"""
	PRE: Takes in a batch
	POST: Will go through it and process it
	
	"""
	with open(fname, "r") as r:
		reader = csv.reader(r, delimiter=',')
		for i,row in enumerate(reader):
			print i, row
			mol = embed_mol(row[0])
			
		
		
if __name__ == "__main__":
	import csv
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	#split_main_file_in_block_of_1000("/home/macenrola/HydrocarbonsBinding/250k_rndm_zinc_drugs_clean_3.csv")
	process_a_batch("/home/macenrola/HydrocarbonsBinding/zinc_pieces/1000-2000-zinc-pieces")