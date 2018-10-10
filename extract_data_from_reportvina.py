def get_list_vina_report(path="/home/macenrola/Thesis/autodock_vina/autodock_vina_1_1_2_linux_x86/ALL_RIGID/", suffix="pdbqt-report-vina"):
	"""
	PRE: -
	POST: returns the list of vina reports
	"""

	import glob
	flist = glob.glob('{}*{}'.format(path, suffix))
	# print '\n'.join(flist)
	return flist

def get_best_pose_val_from_vina_report(vinareport):
	with open(vinareport, 'rb') as r:
		bestpose = None
		for lines in r:
			parts = lines.split()
			# print parts
			# bestpose = None
			if len(parts) == 0: continue
			if parts[0] == '1':
				bestpose = parts[1]
				break
			else: pass
		if bestpose == None:
			bestpose = 1
		return float(bestpose)

def produce_picture_grid_from_SDFFILE_list(LISTOFSDF):
	"""
	Pre: takes a list of SDF Files 
	Post: produces png images containing a fixed amount of guest molecules (no complexes) their pubchem names and their best binding energies as legend
	"""
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem, Draw
	block = []
	legends = []
	SUMFILE = '/'.join(LISTOFSDF[0][1].split('/')[:-1])+'/'+'BESTSVINA'
	print SUMFILE
	for i, line in enumerate(LISTOFSDF):

		pb_number = line[1].strip().split('/')[-1].split('x')[0]
		print pb_number
		binding = line[0]
		sdfFileName = '.'.join(line[1].split('.')[:-1])
		print sdfFileName
		# binding = float(parts[7])
		# e_guest = float(parts[2])
		# guest_nbr = parts[3]

		tmp_mol = Chem.MolFromMolFile(sdfFileName)
		try:
			AllChem.Compute2DCoords(tmp_mol)
		except: continue
		block.append(tmp_mol)
		# legends.append('PB_NBR:{0}/E_BIND:{1:4.4f}({2})/E_G:{3:4.4f}({4})'.format(pb_number, binding, complex_nbr, e_guest, guest_nbr ))
		legends.append('PB_NBR:{0}/BE:{1:4.4f}'.format(pb_number, float(binding)))
		if i%100==99: 
			print SUMFILE+'{}-{}.png'.format(i-100+1, i+1)
			img = Chem.Draw._MolsToGridImage(block, legends=legends, molsPerRow=10, subImgSize=(400,400))
			img.save(SUMFILE+'-{}-{}.png'.format(i-100+1, i+1))
			block = []
			legends = []
	print 'block length is {}'.format(len(block))
	if len(block) != 0:		
		img = Chem.Draw._MolsToGridImage(block, legends=legends, molsPerRow=10, subImgSize=(400,400))
		img.save(SUMFILE+'-{}-{}.png'.format(i-100+1, i+1))
		block = []
		legends = []


if __name__ == "__main__":
	flist = get_list_vina_report()
	poses = [get_best_pose_val_from_vina_report(f) for f in flist]
	couples = zip(poses, flist)
	slist = sorted(couples)[:100]
	produce_picture_grid_from_SDFFILE_list(slist)
	for s in slist:
		print s 
	# bestpose = min([float(x) for x in poses if x != None])
	# for i, f in enumerate(flist):
		# print f, poses[i]

	 # poses.index(min(poses))
