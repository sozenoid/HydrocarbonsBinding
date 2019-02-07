
def dummy_sort(fname):
	# fname = '/home/macenrola/Thesis/hydrocarbons/only_c_errors/pictures_errors_with_Hs/ERRORS_SDF_with_HS.sdf_SUM_SORTED'
	fname = fname
	fout = fname + "-rigid"

	# fname = '/home/macenrola/Thesis/GEOM_OBABEL_ON_THE_FLY_FOR_RDKIT/first_sorted.smi'
	# fout = fname+'_ONLY_RELEVANT'
	with open(fname, 'rb') as r:
		with open(fout, 'wb') as a:
			acc=0
			for i, lines in enumerate(r):
				line = lines.strip().split()
				if i%1000==0:
					print zip(line, range(len(line)))
				if line[-1] == '1':
					a.write('\t'.join(line)+'\n')
				# if float(line[10])>float(line[11]):
				# 	line[10], line[11] = line[11], line[10]
				# a.write('\t'.join(line)+'\n')
				# if min(float(line[11]), float(line[10])) <=-24.0:
				# 	print i, lines, min(float(line[11]), float(line[10]))
				# # 	# a.write(lines)
				# 	a.write('\t'.join(line[:2])+'\n')
				# if i==10:
				# 	break
				
				# if line[0].count('+')==0 and line[0].count('-')==0:
				# 	a.write(lines)
				# a.write(' '.join([line[0], line[1], 'E{} S{}'.format(float(line[10]), float(line[16]))])+'\n')
				# 	a.write(lines)
				# a.write('\t'.join(line[:2])+'\n')
					# a.write(lines)

def assign_ZINCS(canfile, zincfile):
	oufile = canfile+'ZINC15OK'
	with open(canfile, 'rb') as r:	
		canlines = r.readlines()
		smilines = [x.strip().split()[0] for x in canlines]	
	with open(oufile, 'wb'):pass
	with open(zincfile, 'rb') as r:
		for i,line in enumerate(r):
			current_can = line.strip().split()[0]
			if i%100000==0: print i
			if current_can in smilines:
				print 'FOUND {} {}'.format(line, canlines[smilines.index(current_can)])
				with open(oufile, 'ab') as a:
					a.write(canlines[smilines.index(current_can)].strip()+'\t'+line.strip().split()[1]+'\n')

def sort_ZINCs15(ZINC_match_file):
	oufile = ZINC_match_file+'_sorted'
	with open(oufile, 'wb'): pass
	with open(ZINC_match_file, 'rb') as r:
		for lines in r:
			l = lines.strip().split(',')
			with open(oufile, 'ab') as a:
				a.write(l[1]+'\t'+l[0]+'\n') 


############################### FORCE FIELD COMPARE
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
N_STEPS = 1000
tol = 1e-12
import urllib
from subprocess import call
import cPickle
import glob
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem, TorsionFingerprints
from mol_ops import align_mol, create_pdb_and_topology_from_sdf, get_frequency_report, align_mol
# from get_solvation import get_solvation_for_pdb

########### 1) Get the structure from pubchem and print its energy
def get_sdf_from_pbsite(pubchem_id):
	"""Retrieves the pubchem number, saves it to the current file with the number as name
	Takes the pubchem number as string
	"""
	try:
		open('{}.sdf'.format(pubchem_id))
		print('NO NEED TO DOWNLOAD')
	except:
		main_address = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{0}/record/SDF/?record_type=3d&response_type=save&response_basename=Structure3D_CID_{0}'
		urllib.urlcleanup()
		urllib.urlretrieve(main_address.format(pubchem_id), pubchem_id+'.sdf')
	supp = Chem.SDMolSupplier('{}.sdf'.format(pubchem_id))
	# print '{} PUBCHEM ENERGY is {}'.format(pubchem_id, supp[0].GetProp('PUBCHEM_MMFF94_ENERGY'))
	return supp[0].GetProp('PUBCHEM_MMFF94_ENERGY')




def get_ref_energies(reffile='/home/macenrola/Thesis/Calibration_force_field/MMFF94.energies', sdfilevalid='/home/macenrola/Thesis/Calibration_force_field/MMFF94_dative.sdf'):
	"""
	PRE:  takes in the file from mmff94 energies and the corresponding rdkit energies file, THE TWO FILES NEED TO BE SORTED! :)
	POST: checks the energies sdf file and compares the two
	"""
	N_STEPS = 100000
	tol = 1e-12
	with open(reffile, 'rb') as r:
		reflines = r.readlines()[16:]
	supp = Chem.SDMolSupplier(sdfilevalid, removeHs=False)
	i = 0
	for mols in supp:
		ff = AllChem.MMFFGetMoleculeForceField(mols, AllChem.MMFFGetMoleculeProperties(mols,mmffVariant='mmff94'))
		ff.Initialize()
		ff.Minimize(N_STEPS, tol,tol)
		print('Energies {3}vs{4} are ref:{0:>12.4f}\tand rdkit:{1:>12.4f}\tdiff is \t{2:>12.6f}'.format(float(reflines[i].split()[1]), ff.CalcEnergy(), float(reflines[i].split()[1])-ff.CalcEnergy(), reflines[i].split()[0], mols.GetProp('_Name') ) )
		i+=1



############ 2) Optimize the energy using the RDKIT and print its energy
N_STEPS = 10000
tol = 1e-12
def get_rdkit_energy(sdfile):
	def converge(yo):
		AllChem.MMFFSanitizeMolecule(yo)
		ff = AllChem.MMFFGetMoleculeForceField(yo, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(yo, mmffVariant='MMFF94'), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
		converged = ff.Minimize(N_STEPS, tol, tol)
		# converged = 0
		return converged, ff.CalcEnergy()

	mol = Chem.SDMolSupplier(sdfile, removeHs=False)[0]
	flags = converge(mol)
	return flags

#### RDKIT ERRORS
### AMHTAR01 BAOXLM01 BBSPRT10 CAGREH10 CAMALD03





########### Optimize the energy using openbabel and print its energy
def get_OB_ENERGY(sdfile):
	outxyz = 'outfile_obabel.xyz'
	outenergy = 'outfile_energy'
	with open(outxyz, 'wb') as w:
		call('obminimize -ff MMFF94 -cg -n 20000 -c {} {}'.format(tol, sdfile).split(' '), stdout=w)
	with open(outxyz, 'rb') as r:
		lines = r.readlines()
		print(lines[0])
		if b'WARNING' in lines[0]:
			with open(outxyz+'_no_warning', 'wb') as w:
				w.writelines(lines[1:])
			outxyz = outxyz+'_no_warning'
	# outxyz = sdfile

	with open(outenergy, 'wb') as w:
		call('obenergy -ff MMFF94 {}'.format(outxyz).split(' '), stdout=w)
	with open(outenergy, 'rb') as r:
		lines = r.readlines()
		return lines[-1]

def get_energies(smilist):
	with open(smilist, 'rb') as r:
		acc = 0
		for lines in r:
			acc+=1
			# if acc==5:break
			print('#'*30+'\n'+'#'*30)
			pb_number = lines.strip().split('\t')[-1]
			try:
				pb_energy = get_sdf_from_pbsite(pb_number)
			except IndexError:
				print('{} does not have a sdf file on PUBCHEM'.format(pb_number))
				continue
			rdkit_flags = get_rdkit_energy('{}.sdf'.format(pb_number))
			# ob_energy = get_OB_ENERGY('{}.sdf'.format(pb_number))
			print(lines.strip())
			print('PB_ENERGY is {}\nRDKIT_ENERGY is {} kcal/mol\nOBABEL_ENERGY {}'.format(pb_energy, rdkit_flags[1], ob_energy))

def get_energies_calibration(calibration_suite_file):
	supp = Chem.SDMolSupplier(calibration_suite_file, removeHs=False)
	for mol in supp:
		if mol is not None:
			# fname = 'xxx{}.sdf'.format(mol.GetProp('_Name'))
			# w = Chem.SDWriter(fname)
			# w.write(mol)
			# w.close()

			print(fname)
			try:
				rdkit_flags = get_rdkit_energy(fname)
			except: print('RDKIT Failure at {}'.format(fname))
			# ob_energy = get_OB_ENERGY(fname)
			print('RDKIT_ENERGY is {} kcal/mol\nOBABEL_ENERGY {}'.format(rdkit_flags[1], ob_energy))
		else:
			print('{} is NNNNNOOOOOONE'.format(mol))

def keep_only_MMFF94_compatible(smifile):
	with open(smifile, 'rb') as r:
		outfile = smifile + '_only_MMFF94'
		with open(outfile, 'wb') as w:
			for i, line in enumerate(r):
				if i%2000==0: print i, line
				# if i==10: break
				smi = line.split('\t')[0]
				mol = Chem.MolFromSmiles(smi)
				invalid = False
				if mol is not None:
					for atom in mol.GetAtoms():
						if atom.GetAtomicNum() not in [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]:
							invalid = True
							break
						if '+' in smi or '-' in smi:
							invalid = True
							break

					if not invalid:
						w.write(line)

def keep_only_hydrocarbons(smifile):
	with open(smifile, 'rb') as r:
		outfile = smifile + '_only_C_atoms'
		with open(outfile, 'wb') as w:
			for i, line in enumerate(r):
				if i%2000==0: print i, line
				# if i==10: break
				smi = line.split('\t')[0]
				mol = Chem.MolFromSmiles(smi)
				invalid = False
				if mol is not None:
					if mol.GetNumHeavyAtoms() > 20:
						invalid=True

					# for atom in mol.GetAtoms():
					# 	if atom.GetAtomicNum() not in [1,6]:#[1, 6, 7, 8, 9, 15, 16, 17, 35, 53]:
					# 		invalid = True
					# 		break
					# 	if '+' in smi or '-' in smi:
					# 		invalid = True
					# 		break

					if not invalid:
						w.write(line)

def make_sdf_from_dic(pubchem_id, dic):
	block = dic[pubchem_id]
	w = Chem.SDWriter('./{}.sdf'.format(pubchem_id))
	w.write(Chem.MolFromMolBlock(block))
	w.close()

def make_sdf_from_OBABEL(SMILES, pubchem_id):
	call('obabel -:{0} -O {1}.sdf --gen3d'.format(SMILES, pubchem_id).split(' '))
	
def make_conformers_from_sdf(pubchem_id):
	call('obabel {0}.sdf -O {0}_confs.sdf --conformer --nconf 1000 --writeconformers'.format(pubchem_id).split(' '))

def converge_conformers(pubchem_id):
	call('obabel {0}_confs.sdf -O {0}_min_confs.sdf --minimize --ff MMFF94 --steps 50000 --crit 1e-8 --cg --rvdw 20.0 --rele 20.0'.format(pubchem_id).split(' '))

def merge_dic():
	import glob
	dic = {}
	se = set()
	dic_list = glob.glob('/media/macenrola/FAT/tests_dic_failures/failure_dics_part_3/*.sdf_failures_dictionary')
	for its in dic_list:
		print its
		with open(its, 'rb') as r:
			temp_dic = cPickle.load(r)
		for key in temp_dic:
			se.add(key)
	dic = dict.fromkeys(se)
	for its in dic_list:
		print its
		with open(its, 'rb') as r:
			temp_dic = cPickle.load(r)
		for key in temp_dic:
			dic[key] = temp_dic[key]
	with open('/'.join(dic_list[0].split('/')[:-2])+'part3_merged', 'wb') as w:
		cPickle.dump(dic, w)

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def make_error_dic():
	import glob
	with open('/home/macenrola/Thesis/Calibration_force_field/EXTREMES_RETREAT/xaa_xab_failures.can_only_mmff94_atoms.can', 'rb') as r:
		lines = r.readlines()
		lines = [int(x.strip().split('\t')[1]) for x in lines]
		print len(lines)
	return
	dic_list = glob.glob('/media/macenrola/FAT/tests_dic_failures/failure_dics_part_*/*.sdf_failures_dictionary')
	print len(dic_list)
	error_dic = {}
	for dic in dic_list:
		with open(dic, 'rb') as r:
			print 'loading {}'.format(dic)
			temp_dic = cPickle.load(r)
			print 'loaded {}'.format(dic)
			for els in lines:
				if els in temp_dic:
					print els
					error_dic[els] = temp_dic[els]
	with open('/home/macenrola/Thesis/Calibration_force_field/EXTREMES_RETREAT/error_dic', 'wb') as w:
		cPickle.dump(error_dic, w)

def keep_only_best_SDFs(SDFILE, SUM_file):
	"""IN: TAKES A SINGLE SDF FILE WITH COMPLEXES AND HOSTS AND WRITES TO A NEW SDF FILES ONLY BEST E COMPLEXES AND BEST COMPLEXES (2 SDFs total per molecule, less if invalid complexes)"""
	list_of_best_complexes = get_list_of_best_complexes(SUM_file)
	list_of_best_guests = get_list_of_best_guests(SUM_file)
	# for l in list_of_best_complexes:
	# 	print l
	# print list_of_best_guests
	list_of_best_complexes_truncated = [x for x in list_of_best_complexes]
	supp = Chem.SDMolSupplier(SDFILE, removeHs=False)
	writer = Chem.SDWriter(SDFILE+'_ONLY_BEST.sdf')
	a = 0
	for mol in supp:
		a+=1
		# if a==100: break
		if mol is None:
			continue
		else:
			if mol.GetProp('_Name') in list_of_best_guests:
				# print '{} was found'.format(mol.GetProp('_Name'))
				writer.write(mol)
			elif mol.GetProp('_Name') in list_of_best_complexes_truncated:
				# if mol.GetProp('COMPLEX_NUMBER') == list_of_best_complexes[list_of_best_complexes_truncated.index(mol.GetProp('_Name'))][-1]:
				print '{} was found'.format(mol.GetProp('_Name'))
				writer.write(mol)
			writer.flush()
	writer.close()

def get_list_of_best_complexes(SUM_file):
	""" 
		IN a sum file with only two complexes per molecule (no)
	"""
	best_complexes_names = []
	with open(SUM_file, 'rb') as r:
		for i, line in enumerate(r):
			l = line.split()
			sphere_radii = [x for x in zip(range(2),l[16:18]) if float(x[1])<=10.0]
			# if i%100==0: print sphere_radii
			if sphere_radii == []:
				continue
			energies = [(int(k[0]), l[10+int(k[0])]) for k in sphere_radii]
			complex_number_best = min(energies, key=lambda x: float(x[1]))
			# i_valid_radii = find(sphere_radii<=8.0)
			# energies = l[10:13]
			# best_E = min(energies[i_valid_radii])

			best_complexes_names.append('{}-inclusion-complex-attempt-{}-{}'.format(l[1], complex_number_best[0]<2, complex_number_best[0]))
	return best_complexes_names

def get_list_of_best_guests(SUM_file):
	best_guests_names = []
	with open(SUM_file, 'rb') as r:
		for i, line in enumerate(r):
			l = line.split()
			best_guests_names.append('{}-guest-best-E'.format(l[1]))
	return best_guests_names

def get_SDFs_from_canlist(pubchem_list, tar_list, destination_sdf):
	"""T
	Pre: Takes in a list of pubchem numbers and a list of tar files containing each a sdf file with molecules.
		 The molecules have a name spelled as pubchem_number-*
	Post: Writes all the molecules' sdf that match a condition to the destination_sdf file
	"""
	def open_tar(tarfile):
		call('tar -xf {} -C {}'.format(tarfile, '/'.join(tarfile.split('/')[:-1])), shell=True)
	def rm_sdf(sdfile):
		call('rm {}'.format(sdfile), shell=True)
	def write_canSDFs_form_sdfile(pubchem_list, sdwriter, sdfile):
		supp = Chem.SDMolSupplier(sdfile, removeHs =False)
		for mol in supp:
			number = mol.GetProp('PUBCHEM_NUMBER')
			# print number
			if number in pubchem_list:
				print '{} was found'.format(number)
				sdwriter.write(mol)

	def treat_all_tars(pubchem_list, tar_list):
		w = Chem.SDWriter(destination_sdf)
		tarlen = len(tar_list)
		for i, tars in enumerate(tar_list):
			print 'steps {}/{}:\tTreating {}'.format(i, tarlen, tars)
			sdfname = tars[:-4]+'_OUT.sdf'
			open_tar(tars)
			write_canSDFs_form_sdfile(pubchem_list, w, sdfname)
			extract_guests_from_complex(sdfname, sdfname+'_ONLY_GUESTS.sdf')
			rm_sdf(sdfname)
		w.close()

	treat_all_tars(pubchem_list, tar_list)

def get_pubchem_list(sumfile):
	"""
	Pre: From a SUMFILE, extracts a list of pubchem numbers
	POST: returns the list of pubchem numbers
	"""
	with open(sumfile, 'rb') as r:
		lines = r.readlines()
		numbers = [x.split()[1] for x in lines]
		return numbers

def extract_guests_from_complex(SDFfileIn, SDFfileOut):
	"""
	PRE: Takes in a SDFfile with complexes
	POST: Extracts the guests from it, converges it, computes its energy, and associates as molecule properties the old guest energy, the old complex energy along with the recent energy of the guest 
	"""
	n_steps = 100000
	tol = 1e-9
	supp = Chem.SDMolSupplier(SDFfileIn, removeHs=False)
	wr = Chem.SDWriter(SDFfileOut)
	mols = []
	for mol in supp:
		Guest = Chem.GetMolFrags(mol, asMols=True)[1]
		AllChem.MMFFSanitizeMolecule(Guest)  
		ff = AllChem.MMFFGetMoleculeForceField(Guest, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(Guest, mmffVariant='MMFF94'), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
		ff.Initialize()
		cf = ff.Minimize(n_steps, tol, tol)
		Guest.SetProp('GUEST_ENERGY_WHEN_EXTRACTED', '{0:4.4f}'.format(ff.CalcEnergy()))
		prop_dict = mol.GetPropsAsDict(includePrivate=True)
		for prop in prop_dict:
			Guest.SetProp(prop, str(prop_dict[prop]))
		wr.write(Guest)
	wr.close()


def keep_guests(SDFfile):
	"""
	PRE:  Takes in a SDF with diverse molecules
	POST: Produces a SDF file SDFfile_GUESTS with only the molecule with a name containing 'best-guest-E'
	"""
	supp = Chem.SDMolSupplier(SDFfile, removeHs=False)
	w = Chem.SDWriter(SDFfile+'_GUESTS')
	for mol in supp:
		if mol is None: continue
		if '-guest-best-E' in mol.GetProp('_Name'):
			w.write(mol)
		else:
			pass
	w.close()


def split_sdf_file_individually(sdfile):
	supp = Chem.SDMolSupplier(sdfile, removeHs = False)
	acc = 0
	for mol in supp:
		acc+=1
		# if acc==10: break
		if mol==None:
			continue
		else:
			name = mol.GetProp('_Name')
			print name
			try:
				c_nbr = '-'+mol.GetProp('COMPLEX_NUMBER')
			except:
				c_nbr = ''
			w = Chem.SDWriter('/'.join(sdfile.split('/')[:-1])+'/'+name+c_nbr+'.sdf')
			w.write(mol)
			w.close()

def make_pdb_complex_with_named_residues(sdf_guest, sdf_CB):
	supp_guest = Chem.SDMolSupplier(sdf_guest, removeHs=False)
	guest = supp_guest[0]
	# guest.SetProp('_Name', 'GUE')
	# guest =Chem.MolFromPDBFile('/home/macenrola/Thesis/AMBER/mege_test/cb7_b2.pdb')
	# print Chem.MolToPDBBlock(guest)
	atm_dic = {}
	for atom in guest.GetAtoms():
		# print atom.GetSymbol()
		if atom.GetSymbol() not in atm_dic:
			atm_dic[atom.GetSymbol()] = 1
		else: atm_dic[atom.GetSymbol()] += 1

		atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} GST'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10))))
		# atom.SetMonomerInfo(Chem.rdchem.AtomMonomerInfo().SetName('{} GST '.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10))))
		# atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo(' {}   GST '.format(atom.GetSymbol())))
		atom.GetMonomerInfo().SetResidueNumber(1)
		# atom.SetProp(MonomerInfo, 42)
	flavour = 28
	print Chem.MolToPDBBlock(guest, flavor = flavour)
	# writer = Chem.PDBWriter(sdf_guest.rsplit('.',1)[0]+'.pdb')
	# writer.write(guest)
	# writer.close()
	Chem.MolToPDBFile(guest, sdf_guest.rsplit('.',1)[0]+'.pdb', flavor=flavour)
	fix_PDB_spacing(sdf_guest.rsplit('.',1)[0]+'.pdb')

	supp_CB = Chem.SDMolSupplier(sdf_CB, removeHs=False)
	CB = supp_CB[0]
	# CB.SetProp('_Name', 'CUC')
	atm_dic = {}
	for atom in CB.GetAtoms():
		if atom.GetSymbol() not in atm_dic:
			atm_dic[atom.GetSymbol()] = 1
		else: atm_dic[atom.GetSymbol()] += 1
		# atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('CB'))
		atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} CB7'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10))))
		# atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo(' {}   CB7 '.format(atom.GetSymbol())))
		atom.GetMonomerInfo().SetResidueNumber(2)
	print Chem.MolToPDBBlock(CB, flavor = flavour)
	Chem.MolToPDBFile(CB, sdf_CB.rsplit('/',1)[0]+'/CB.pdb', flavor=flavour)
	fix_PDB_spacing(sdf_CB.rsplit('/',1)[0]+'/CB.pdb')
	# writer = Chem.PDBWriter(sdf_CB.rsplit('/',1)[0]+'/CB.pdb')
	# writer.write(CB)
	# writer.close()

	complex_CB_host = Chem.CombineMols(guest, CB)
	yo = Chem.GetMolFrags(complex_CB_host)
	for els in yo:
		print els
	complex_CB_host.SetProp('_Name', 'CMP')
	print Chem.MolToPDBBlock(complex_CB_host, flavor=flavour)
	Chem.MolToPDBFile(complex_CB_host, sdf_CB.rsplit('/',1)[0]+'/COMPLEX.pdb', flavor=flavour)
	remove_CONNECT_LINES(sdf_CB.rsplit('/',1)[0]+'/COMPLEX.pdb')
	fix_PDB_spacing(sdf_CB.rsplit('/',1)[0]+'/COMPLEX.pdb')
	# writer = Chem.PDBWriter(sdf_CB.rsplit('/',1)[0]+'/COMPLEX.pdb')
	# writer.write(complex_CB_host)
	# writer.close()


def remove_CONNECT_LINES(fname):
	with open(fname, 'rb') as r:
		lines = r.readlines()
	with open(fname, 'wb') as w:
		w.writelines([x for x in lines if 'CONECT' not in x if 'MASTER' not in x][1:])

def fix_PDB_spacing(fname):
	"""The PDB files is formated by rdkit MolToPDBFile flavour 28 without any MASTER, TER, CONECT or Charge flags, 
	it gets out formated with the adequate amount of blank space to be read by AMBER, 
	This method destroys the original file"""
	raw_spc = [7, 2, 4, 4, 4, 12, 8, 8, 6, 6, 12, 0]
	new_lines = []
	cb_yet = False
	with open(fname, 'rb') as r:
		for line in r:
			if 'CB7' in line and not cb_yet:
				cb_yet = True
				new_lines.append('TER')
			if 'ATOM' in line:
				line_content = line.split()
				line_content.insert(4, 'A')
				# print line_content
				ls = [len(x) for x in line_content]
				actual_spacing = [raw_spc[0]-ls[1], # after ATOM 
				raw_spc[1],  # after ATOM# (11)
				raw_spc[2]-ls[2], # after ATOM NAME (H41)
				raw_spc[3]-ls[3], # after RESIDUE NAME (CUC)
				raw_spc[4]-ls[4], # after CHAIN ID (A)
				raw_spc[5]-ls[6], # after RESIDUE NUMBER (1)
				raw_spc[6]-ls[7], # after X CART COORDINATE (6.171)
				raw_spc[7]-ls[8], # after Y CART COORDINATE (3.377)
				raw_spc[8]-ls[9], # after Z CART COORDINATE (21.096)
				raw_spc[9]-ls[10], # after enigmatic number (1.00)
				raw_spc[10]-ls[11], # after partial charge (0.00)
				raw_spc[11], # after ATOMIC SYMBOL (N)
				]
				new_lines.append(''.join([x[0]+' '*x[1] for x in zip(line_content, actual_spacing)]))
				print new_lines[-1]

			else:
				new_lines.append(line.strip())
				print new_lines[-1]
	with open(fname, 'wb') as w:
		w.writelines('\n'.join(new_lines)+'\n')


def get_infos_from_sdffilelist(sdffilelist):
	"""
	PRE: Takes in a list of sdffiles with possible multiple molecules per file
	POST: Extrats the relevant informations, pubchem number and energy
	"""
	molIds = {}
	# print sdffilelist
	for fname in sdffilelist:
		with open(fname, 'rb') as r:
			print len(r.readlines()), fname
			if r.readlines() == ['']: continue
		supp = Chem.SDMolSupplier(fname, removeHs=False)
		for mol in supp:
			if mol is None:
				continue
			try:
				pubchemNumber = mol.GetProp('PUBCHEM_NUMBER')
			except: 
				print 'error here with: {}'.format(fname)
				continue
			complexNumber = int(mol.GetProp('COMPLEX_NUMBER'))
			bindingEnergy = mol.GetProp('BINDING_ENERGY')
			if complexNumber in [0,1]:
				if pubchemNumber in molIds:
					if molIds[pubchemNumber][0] > float(bindingEnergy):
						molIds[pubchemNumber] = (float(bindingEnergy), complexNumber)
				else:
					molIds[pubchemNumber] = (float(bindingEnergy), complexNumber)

	cPickle.dump(molIds, open('/'+'/'.join(fname.split('/')[:-1])+'/'+'molecules_to_keep_from_extracted', 'wb'))

def keep_mols_that_match_pbnb_bindingE(sdfflistin, infodic, sdfout):
	"""
	PRE: Takes in a list of sdf files, a dictionary containing the pubchem numbers as keys and the binding energies as values of the molecules to be kept, and an outfile name for the output sdf with the molecules kept.
	POST: Writes out all the molecules that match the informations contained in infodic to the sdf file sdfout
	"""
	infodic = cPickle.load(open(infodic, 'rb'))
	print len(infodic.keys())
	w = Chem.SDWriter(sdfout)
	for fname in sdfflistin:
		supp = Chem.SDMolSupplier(fname, removeHs = False)
		for mol in supp:
			if mol is None:
				continue
			try:
				pubchemNumber = mol.GetProp('PUBCHEM_NUMBER')
				complexNumber = mol.GetProp('COMPLEX_NUMBER')
			except: 
				print 'error here with: {}'.format(fname)
				continue
			bindingEnergy = float(mol.GetProp('BINDING_ENERGY'))
			if pubchemNumber in infodic:
				# print infodic[pubchemNumber]
				# print  (bindingEnergy, complexNumber)
				if infodic[pubchemNumber] == (bindingEnergy, int(complexNumber)):
					w.write(mol)
	w.close()
def swap_id_smiles(fname):
	"""
	PRE: Takes in a file where pubchem id and smiles are in the wrong order 
	POST: And put them back in the right order
	"""
	with open(fname , 'rb') as r:
		lines = r.readlines()

	for l in lines:
		print l.split()
	with open(fname+'hahhha', 'wb') as w:
		w.writelines([x.strip().split()[1] +'\t'+x.strip().split()[0]+'\n' for x in lines])

def converge_AMBER_GAFF(pdbfile):
	call('{}/one_minimalist.sh {}'.format('/'.join(pdbfile.split('/')[:-1]), pdbfile.split('/')[-1]).split(' '))

def produce_picture_grid(sdfile):
	"""
	PRE: Takes in a SDF File with binding energy information
	POST: Produce a png grid image, with legends as names and binding energy and inclusion
	"""
	supp =  Chem.SDMolSupplier(sdfile, removeHs=False)
	legends = []
	for mol in supp:
		AllChem.Compute2DCoords(mol)
		if 'inclusion' in mol.GetProp('_Name'): binding = ' ' + mol.GetProp('BINDING_ENERGY')
		else: binding = ''
		legends.append('-'.join(mol.GetProp('_Name').split('-')[:2])+binding)
	img = Chem.Draw._MolsToGridImage(supp, legends=legends, molsPerRow=10, subImgSize=(400,400))
	img.save(sdfile+'img.png')

def produce_picture_grid_from_SUM(SUMFILE):
	"""
	Pre: takes a SUM file formatted as follows:
			[OH+]1CC[OH+]CC[OH+]CC[OH+]CC[OH+]CC[OH+]CC1	5150062	1.0000 877.3256 5.8075 11.7531	3.0000 889.0905 5.5935 12.2175	-481.9217 -482.5407	6	True True True	 7.6088  7.7197
	Post: produces png images containing a fixed amount of guest molecules (no complexes) their pubchem names and their best binding energies as legend
	"""
	with open(SUMFILE, 'rb') as r:
		block = []
		legends = []
		for i, line in enumerate(r):
			parts = line.strip().split()

			pb_number = parts[1]
			smiles = parts[0]
			print smiles
			# binding = float(parts[7])
			binding = parts[-1]
			# e_guest = float(parts[2])
			# guest_nbr = parts[3]

			tmp_mol = Chem.MolFromSmiles(smiles)
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



def add_new_binding_energy(MOLDIC1, MOLDIC2, SUMFILE):
	"""
	PRE : Takes in a molecule dictionary with pubchem numbers as keys and [energy of extracted guest, energy of best guest as values] plus a sum file
	POST : Generates a new sumfile with a new binding energy appened together with the former energy of the best guest (as a sanity check)
	E_bind = E_complex - E_guest - E_CB
	E_bind_new = E_bind + E_guest_old - E_guest_new
	"""
	moldic = cPickle.load(open(MOLDIC1, 'rb'))
	moldic.update(cPickle.load(open(MOLDIC2, 'rb')))
	with open(SUMFILE+'_ERRORS', 'wb'): pass 
	print 'moldic has {} kkeys'.format(len(moldic.keys()))
	with open(SUMFILE+'_new_binding', 'wb') as w:
		with open(SUMFILE, 'rb') as r:
			for i, line in enumerate(r):
				parts = line.strip().split()
				try:
					new_E = float(moldic[parts[1]][0])
				except KeyError as e:
					with open(SUMFILE+'_ERRORS', 'ab') as a:
						a.write(line)
				w.write('\t'.join(parts) + '\t{0:.4f}\t{1:.4f}\n'.format(new_E, float(parts[10]) + float(parts[3]) - new_E))
				if i%10000==0:
					print '{}/{}\t{}'.format(i, 4844587, parts[1])

def get_the_guests(sumfile, sdflist):
	"""
	PRE : Takes in a sum file list and a list of sdf files
	POST : Retains only the sdfs that occur in the pubchem_list
	"""
	pubchemlist = get_pubchem_list(sumfile)
	for fname in sdflist:
		print fname
		supp = Chem.SDMolSupplier(fname, removeHs=False)
		w = Chem.SDWriter(fname+'_EXTRACTED.sdf')
		for mol in supp:
			if mol.GetProp('PUBCHEM_NUMBER') in pubchemlist:
				w.write(mol)


def show_conf_vs_opti(fname):
	"""
	PRE : Takes in a SUM file obtained for two standard complexes and assumes it has been treated to yield 10 complexes so that another file SUM exists
	POST : Shows the difference in best conformer energy and binding energy for the two approaches
	"""
	import ast
	with open(fname, 'rb') as r:
		with open(fname+'_SUM') as rc:
			block = []
			legends = []
			IMG_SIZE = 500
			for i, line in enumerate(r):
				parts = line.strip().split()
				l = rc.readline().strip().split('\t')
				pb_number = parts[1]
				best_E200 = parts[2]
				best_E10 = float(l[2])
				smiles = parts[0]
				binding = min(float(parts[10]), float(parts[11]))
				try:
					binding10 = ast.literal_eval(l[-1].split(',')[1].strip())
				except:
					binding10 = 0.0
								
				# print pb_number, smiles, binding
				tmp_mol = Chem.MolFromSmiles(smiles)
				AllChem.Compute2DCoords(tmp_mol)
				block.append(tmp_mol)
				legends.append('PB:{0}#EB:{1}/{2}#EG:{3}/{4:.4f}'.format(pb_number, binding, binding10, best_E200, best_E10))
				if i%100==99: 
					print fname+'{}-{}.png'.format(i-100+1, i+1)
					img = Chem.Draw._MolsToGridImage(block, legends=legends, molsPerRow=10, subImgSize=(IMG_SIZE,IMG_SIZE))
					img.save(fname+'-{}-{}.png'.format(i-100+1, i+1))
					block = []
					legends = []
			print 'block length is {}'.format(len(block))
			if len(block) != 0:		
				img = Chem.Draw._MolsToGridImage(block, legends=legends, molsPerRow=10, subImgSize=(IMG_SIZE,IMG_SIZE))
				img.save(fname+'-{}-{}.png'.format(i-100+1, i+1))
				block = []
				legends = []

def read_text_as_repr(fname):
	"""
	PRE: Reads a file
	POST: prints the content of the file as a string for hardcoding purposes
	"""
	with open(fname, 'rb') as r:
		print ''.join(r.readlines()).__repr__()

def get_rotable_bond_number(SMILES):
	"""
	PRE : Takes in a smiles and 
	POST :	returns the number of rotatable bonds of that molecule as defined by the rdkit
	"""
	mol = Chem.MolFromSmiles(SMILES)
	mol = Chem.AddHs(mol)
	print Chem.rdMolDescriptors.CalcNumRotatableBonds(mol, Chem.rdMolDescriptors.NumRotatableBondsOptions.Strict)
	AllChem.EmbedMolecule(mol, AllChem.ETKDG())
	Chem.MolToMolFile(mol, '/home/macenrola/Desktop/ROTABLE.sdf')
	return

def extract_best_host_guest_sdf(fnameSUM):
	"""
	PRE  : Takes in a SUM file that tells the number of the best guest and best host, it is assumed that the name of the associated SDF is fnameSUM[:-4]+'_OUT.sdf' 
		   This sum file contains ONE line
	POST : Return the best complex and best guest as defined from the SDF file; both as isolated molecules, namely fnameSUM + '_GUEST{}'.format(bestguestnumber) and fnameSUM + '_COMPLEX{}'.format(besthostnumber)
	"""
	with open(fnameSUM, 'rb') as r:
		lines = r.readlines()[0].strip().split()
	guestnumber = lines[3]
	complexnumber =lines[9]
	
	fnameSDF = fnameSUM[:-4]+'_OUT.sdf'
	extract_host_guest_from_SDF(fnameSDF, [guestnumber], [complexnumber])
	return

def extract_host_guest_from_SDF(fnameSDF, guestnumber, hostnumber, newdestination):
	"""
	PRE  : Takes in a SDF file and the number of a guest and of a hostnumber, guests need to be named 91540321-conf-142 (pubchemnumber-conf-guestnumber) and complexes 91540321-inclusion-complex-attempt-193
		   guestnumber and hostnumber are STRINGS.

	POST : Return the best complex and best guest as defined from the SDF file; both as isolated molecules, namely fnameSDF[:-4]+'_GUEST{}.sdf'.format(guestnumber) and fnameSDF[:-4]+'_COMPLEX{}.sdf'.format(hostnumber)
	"""
	name_sdf = fnameSDF.split('/')[-1]
	print name_sdf
	supp = Chem.SDMolSupplier(fnameSDF, removeHs=False)
	for mol in supp:
		name = mol.GetProp('_Name').split('-')
		if 'conf' == name[1] and name[2] == guestnumber:
			Chem.MolToMolFile(mol, newdestination + '/' + name_sdf[:-4]+'_GUEST{}.sdf'.format(name[2]))
		if 'inclusion' == name[1] and name[-1] == hostnumber:
				Chem.MolToMolFile(mol, newdestination + '/' + name_sdf[:-4]+'_COMPLEX{}.sdf'.format(name[-1])) 

def reorganizes_sum_file_list(fname):
	"""
	PRE  : Takes a file and rewrites it as follows
	POST : Take a line and if the next line doesn't contain the SUM writes it along, writes nothing if the next line contains a SUM
	"""
	with open(fname, 'rb') as r:
		for k,_ in enumerate(r):
			pass
	print k
	with open(fname, 'rb') as r:
		with open(fname+'_rearranged', 'wb') as w:
			i = 0
			previousline = ''
			while i<k+1:
				line = r.next().strip()
				if 'SUM' in previousline:
					if 'SUM' in line:
						w.write(previousline + '\n')
					else:
						w.write(previousline + '\t' + line + '\n')
				previousline = line
				i+=1

def return_carbon_filtration_string():
	"""
	RETURNS THE CARBON AND HYDROGEN ONLY string for open below_20_pubchem_no_undesirable_only_mmff94_atoms_no_salts
	"""
	"""
	obabel sample.can -osmi -v'[!#6;!#1]'obabel sample.can -osmi -v'[!#6;!#1]'
	"""
	pass 

def keep_only_the_carbons(fname):
	"""
	PRE: Takes in a can file
	POST : returns only the all carbons ones
		   DOES NOT WORK TOO WELL..........................
		   PLEASE FAVOUR OPENBABEL WITH FILTRATION STRING ABOVE
	"""
	okaylist = ['C', 'c', '(', ')', '[', ']','#','=','\\','/','H']
	for num in range(10):
		okaylist.append(str(num))
	with open(fname, 'rb') as r:
		with open(fname+'_only_carbs_lol.can', 'wb') as w:
			for line in r:
				invalid = False
				for els in line.strip().split()[0]:
					# print els
					if els not in okaylist:
						invalid = True
						break
				if not invalid:
					w.write(line)


def get_the_hydrocarbons_not_treated(fname_treated_list, fname_all_of_them):
	"""
	PRE: Takes in a sum file with all the treated hydrocarbons (without ERRs then), takes as well a list of all hydrocarbons
	POST : Creates a list of all those not cited in the treated list
	"""
	not_treated_set = set()
	with open(fname_all_of_them, 'rb') as r:
		for line in r:
			current = line.strip().split()[1]
			not_treated_set.add(current)

	with open(fname_treated_list, 'rb') as r:
		for line in r:
			current = line.strip().split()[1]
			if current in not_treated_set:
				not_treated_set.remove(current)

	for i, number in enumerate(not_treated_set):
		print i, number
		if number == '9238':
			print 'ADAMANTANE NOT BLOODY TREATED'+'\n #####'*30
			break

	with open(fname_all_of_them, 'rb') as r:
		with open(fname_all_of_them+'_remain_to_be_treated.can', 'wb') as w:
			for line in r:
				current = line.strip().split()[1] 
				if current in not_treated_set:
					w.write(line)

def return_uncorrelated_conformations(fnameSDF, correlationThresh = 0.01):
	"""
	PRE  : Takes in a sdf file where there are multiple guest conformations with name as 22022682-conf-53, and complexes named as 22022682-inclusion-complex-attempt-154
		   The TFD is performed using the RDKIT, himself inspired from 
		   Scharfer, Christin, et al. "CONFECT: conformations from an expert collection of torsion patterns." ChemMedChem 8.10 (2013): 1690-1700.
	POST : Builds two molecules with various conformers, cluster them using BUTANA RMS clustering and returns the various conformations uncorrelated to a certain rms threshold as a SDF file: fname uncorrelated
	"""
	#########
	#   FIRST BUILDS THE MOLECULES
	#########
	make_sample(fnameSDF)
	supp = Chem.SDMolSupplier(fnameSDF, removeHs=False)
	isComplex = (len(Chem.GetMolFrags(supp[0])) != 1)  
	if isComplex:
		return []
	guests = None
	guestsE = {}
	complexes = None
	complexesE = {}
	for mol in supp:
		name = mol.GetProp('_Name')
		if 'conf' in name:
			if guests == None:
				guests = mol
				guestsE[0] =  mol.GetProp('Guest:Total')
			else:
				guestsE[guests.AddConformer(mol.GetConformer(), assignId=True)] =  mol.GetProp('Guest:Total')
		if 'inclusion' in name:
			if complexes == None:
				complexes = mol
				complexesE[0] = mol.GetProp('Complex:Total')
			else:
				complexesE[complexes.AddConformer(mol.GetConformer(), assignId=True)] = mol.GetProp('Complex:Total')
	
	########
	# Then aligns the molecules PCA style
	# Not needed for a tfd approach
	########	
	# guests_aligned = align_all_conformers_of_molecule(guests)
	guests_aligned = guests

	########
	# Then creates the correlation matrix
	########

	num = guests_aligned.GetNumConformers()
	# RMS APPROACH
	# rmsmat = AllChem.GetConformerRMSMatrix(guests_aligned, prealigned=False)
	# rms_clusters = Butina.ClusterData(rmsmat, num, correlationThresh, isDistData=True, reordering=True)

	# TFD APPROACH
	tfdmat = TorsionFingerprints.GetTFDMatrix(guests_aligned)
	tfd_clusters = Butina.ClusterData(tfdmat, num, correlationThresh, isDistData=True, reordering=True)

	res = [['{}:{}'.format(i, guestsE[i]) for i in x] for x in tfd_clusters]
	# for x in res:
	# 	print x
	return [x[0].split(':')[0] for x in res]


def mv_guest_from_a_to_b(fname, folderA='/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/splitby1hydrocarbon-19', folderB='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/butina'):
	"""
	PRE: Takes a file fname to copy from foldernameA to foldernameB
	POST: idem pre
	"""
	try:
		open('{}/{}'.format(folderB, fname), 'rb')
	except:
		'{} not there'.format(fname)
		call('cp {}/{} {}/'.format(folderA, fname, folderB), shell=True)
		call('tar -xvf {0}/{1} -C {0}'.format(folderB, fname), shell=True)
	return

def untar(fname, folderA='/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS'):
	"""
	PRE : Takes in a file fname and its location folder A
	POST: untars it on location to yield the sum file error file and sdf file
	"""
	if fname[:3]=='not':	
		folderA = folderA+'/'+fname.split('/')[0]
		fname = 'not_treated_'+fname.split('/')[1]
		print folderA, fname
		cmd = 'tar -xvf {0}/{1} -C {0}'.format(folderA, fname)
		call(cmd, shell=True)
		print cmd
	else:
		folderA = folderA+'/'+fname.split('/')[0]
		fname = fname.split('/')[1]
		print folderA, fname
		cmd = 'tar -xvf {0}/{1} -C {0}'.format(folderA, fname)
		call(cmd, shell=True)
		print cmd		

def align_all_conformers_of_molecule(rdkitmol):
	"""
	PRE : Takes in a rdkit molecule with many conformers and 
	POST : aligns them using pca and 0 centering the output is another rdkit mol 
	"""
	conflist = []
	aligned_conflist = []
	for i in range(rdkitmol.GetNumConformers()):
		conflist.append(Chem.MolToMolBlock(rdkitmol, confId=i))
		aligned_conflist.append(align_mol(conflist[-1]))
 
	aligned_rdkit_mol = None

	for i in range(rdkitmol.GetNumConformers()):
		if i == 0:
			aligned_rdkit_mol = Chem.MolFromMolBlock(aligned_conflist[i], removeHs=False)
		else:
			aligned_rdkit_mol.AddConformer(Chem.MolFromMolBlock(aligned_conflist[i], removeHs=False).GetConformer(), assignId=True)

	return aligned_rdkit_mol

def make_sample(sdfname):
	"""
	PRE : Takes in a molecule with plenty of conformers 
	POST : Writes out a sample for visual inspection as an sdf file
	"""
	supp = Chem.SDMolSupplier(sdfname, removeHs=False)
	Chem.MolToMolFile(supp[0], '{}_SAMPLE.sdf'.format(sdfname))


def getindex(index_moleculenumber_filenames='/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_NAMED_HYDROCARBONS_rearranged'):
	"""
	PRE  : 
	PAST :
	"""
	index = {}
	with open(index_moleculenumber_filenames, 'rb') as i:
		for line in i:
			parts = line.strip().split()
			tarname = parts[0][:-4]
			index[parts[-1]] = tarname
	return index

def extract_best_guest_and_complex_from_tars(sumlist, location_of_tars='/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS', 
	newlocation_of_tars='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid'):
	"""
	PRE  : Takes in a sumlist, a file with the tars corresponding to that sumlist and a new folder where to copy the best guests and complexes as sdf files
	POST : Will extract from each line of the sum file the number of the best complexes and best guest and write them to the newlocation_of_tars destination as two separate sdf files
	"""
	index = getindex()
	with open(sumlist, 'rb') as r:
		for i, line in enumerate(r):
			# if i==3:break
			try:
				parts = line.strip().split()
				fname = index[parts[1]]
				print 'processing {} {}'.format(i,fname)

				untar('{}.tar'.format(fname, folderA=location_of_tars))

				extract_host_guest_from_SDF('{}/{}_OUT.sdf'.format(location_of_tars, fname), parts[3], parts[9], newlocation_of_tars)

				call('rm {}/{}_OUT.sdf'.format(location_of_tars, fname), shell=True)
			except Exception:
				with open(newlocation_of_tars+'/'+'ERRFILE', 'ab') as a:
					a.write('COPYSDF\t'+line) 



def get_no_rotatable_bonds(SUMlist, location_of_tars='/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS', 
	index_moleculenumber_filenames='/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_NAMED_HYDROCARBONS_rearranged',
	newlocation_of_tars='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid'):
	"""
	PRE  : Reads in a SUMlist formatted as follows: 
		  c1(c2ccccccccc2)ccccccccc1	123229859	347.935127228	97	True	(0, -1184.6390863392066)	-81.1591	 7.7442	196	

	POST : Will create a new file with for each line the number of significantly different conformers found (1 if no rotatable bonds)
		   This will be done by:
		   - finding the filename in the index,
		   - copying the relevant matching tar file and unpacking it,
		   - extracting the number of uncorrleated conformers
		   - writing out the numer to a new sum list 
	"""
	fout = SUMlist + '_with_uncorrelated_conformers'
	index = getindex(index_moleculenumber_filenames)


	with open(SUMlist, 'rb') as r:
		with open(fout, 'ab') as w:
			for k,line in enumerate(r):
				# if k<45640: continue
				parts = line.strip().split()
				fname = index[parts[1]]
				print k, fname
				# mv_guest_from_a_to_b('{}.tar'.format(fname, folderA=location_of_tars, folderB=newlocation_of_tars)) # FOR MOVING AROUND THE FILES
				print '{}.tar'.format(fname, folderA=location_of_tars)
				try:
					untar('{}.tar'.format(fname, folderA=location_of_tars))
					# uncorr_confs_num = len(return_uncorrelated_conformations('{}/{}_OUT.sdf'.format(newlocation_of_tars, fname))) # FOR MOVING AROUND THE FILES
					uncorr_confs_num = len(return_uncorrelated_conformations('{}/{}_OUT.sdf'.format(location_of_tars, fname)))

					w.write('\t'.join(parts) + '\t{}\n'.format(uncorr_confs_num))
					rmcmd = 'rm {}/{}_OUT.sdf'.format(location_of_tars, fname)
					print rmcmd
					call('rm {}/{}_OUT.sdf'.format(location_of_tars, fname), shell=True)
				except Exception:
					w.write('\t'.join(parts) + '\t{}\n'.format('ERROR'))

def add_solvation_to_sumlist(SUMlist, location_of_sdf='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid'):
	"""
	PRE : The sdf in question are the "best SDFs" obtained from a big SDF file containing 200 conformations of the guests and complexes obtained using "extract_best_host_guest_sdf" method hence formatted as '{}/{}_OUT_GUEST{}.sdf'.format(location_of_sdf, fname, parts[3])
		  The sumfile in question is expected to contain the index number of the best complexes and SDF (to have been produced along those SDF in the same process of docking basically)
		  get_solvation.py must be imported and in path
	POST : Produces a new sum file with the solvation energy of the molecules added as the tail of the sum file as apolar guest , polar guest, apolar complex, polar complex, total apolar, total elec, and total apolar + total elec
		   The last number basically represent the energy penalty upon binding
	"""
	fout = SUMlist + '_with_solvation'
	index=getindex()

	with open(SUMlist, 'rb') as r:
		with open(fout, 'wb') as w:
			for k,line in enumerate(r):
				try:
				# if k<680: continue
					parts = line.strip().split()
					fname = index[parts[1]].split('/')[-1]
					print k, fname
					# if k == 5: break
					guest_file = '{}/{}_OUT_GUEST{}.sdf'.format(location_of_sdf, fname, parts[3])
					complex_file = '{}/{}_OUT_COMPLEX{}.sdf'.format(location_of_sdf, fname, parts[9])
					create_pdb_and_topology_from_sdf(guest_file)
					create_pdb_and_topology_from_sdf(complex_file)
					guest_solv = get_solvation_for_pdb(guest_file[:-4]+'.pdb')
					complex_solv = get_solvation_for_pdb(complex_file[:-4]+'.pdb')
					total_solv = (float(complex_solv[0])-float(guest_solv[0])-10.6300, float(complex_solv[1])-float(guest_solv[1])-(-33.0814))
					solv_penalty = float(total_solv[0]) + float(total_solv[1])
					parts.extend([guest_solv[0], guest_solv[1], complex_solv[0], complex_solv[1], total_solv[0], total_solv[1], solv_penalty])

					w.write('\t'.join([str(x) for x in parts]) + '\n')
				except Exception:
					with open(location_of_sdf+'/'+'ERRFILE', 'ab') as a:
						a.write('SOLVATION\t'+line) 

def add_entropy_to_sumlist(SUMlist, location_of_sdf='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid'):
	"""
	PRE : The sdf in question are the "best SDFs" obtained from a big SDF file containing 200 conformations of the guests and complexes obtained using "extract_best_host_guest_sdf" method hence formatted as '{}/{}_OUT_GUEST{}.sdf'.format(location_of_sdf, fname, parts[3])
		  The sumfile in question is expected to contain the index number of the best complexes and SDF (to have been produced along those SDF in the same process of docking basically)
		  THOSE SDFs ARE ASSUMED TO HAVE BEEN PROCESSED TO YIELD PDB AND PRMTOP FILES such as '{}/{}_OUT_GUEST{}.pdb'.format(location_of_sdf, fname, parts[3]) '{}/{}_OUT_GUEST{}.prmtop'.format(location_of_sdf, fname, parts[3])
		  mol_ops.py must be imported and in path
	POST : Produces a new sum file with the solvation energy of the molecules added as the tail of the sum file as apolar guest , polar guest, apolar complex, polar complex, total apolar, total elec, and total apolar + total elec
		   The last number basically represent the energy penalty upon binding
	"""
	fout = SUMlist + '_with_entropy'
	index = getindex()

	with open(SUMlist, 'rb') as r:
		with open(fout, 'ab') as w:
			for k,line in enumerate(r):
				if k<7410: continue
				parts = line.strip().split()
				fname = index[parts[1]].split('/')[-1]
				print k, fname
				# if k == 20: break
				guest_file_pdb = '{}/{}_OUT_GUEST{}.pdb'.format(location_of_sdf, fname, parts[3])
				guest_file_prmtop = '{}/{}_OUT_GUEST{}.prmtop'.format(location_of_sdf, fname, parts[3])
				complex_file_pdb = '{}/{}_OUT_COMPLEX{}.pdb'.format(location_of_sdf, fname, parts[9])
				complex_file_prmtop = '{}/{}_OUT_COMPLEX{}.prmtop'.format(location_of_sdf, fname, parts[9])
				guest_entropy, guest_energy = get_frequency_report(guest_file_pdb, guest_file_prmtop)
				complex_entropy, complex_energy = get_frequency_report(complex_file_pdb, complex_file_prmtop)
				if 'IMPROPER_FREQS' in [guest_entropy, complex_entropy]:
					w.write('\t'.join(parts)+'\tEntropyFailed\n')
					continue
				entropy_penalty = float(complex_entropy) - float(guest_entropy) - 319.095#230.549  cal/mol.K for just vibrational or 319.095 for total
				energy_binding_amber = float(complex_energy) - float(guest_energy) - 2.5739720226e+00 # 2.5739720226e+00 for lone CB
				entropy_penalty300 = -entropy_penalty * 300/1000.0 # kcal/mol@300K == -TdS
				parts.extend([guest_entropy, complex_entropy, entropy_penalty, entropy_penalty300, energy_binding_amber])
				w.write('\t'.join([str(x) for x in parts]) + '\n')


def reformat_after_solvation_entropy(SUMFILE):
	"""
	PRE : Takes in a sum file that has gone through solvation and energy appending
	POST : Returns another SUM file with the total energy and removes the irreleavant individual contributions
	"""
	SUMout = SUMFILE + '_reformatted'
	with open(SUMFILE, 'rb') as r:
		with open(SUMout, 'wb') as w:
			for i, line in enumerate(r):
				parts = line.strip().split()
				print len(parts)
				if len(parts)!=23:
					continue
				# print parts
				relevantpars = parts[:11]
				additionnal_parts = [float(x) for x in [parts[17], parts[21], parts[22]]]
				relevantpars.extend(additionnal_parts)
				relevantpars.append(sum(additionnal_parts))
				# print relevantpars
				w.write('\t'.join([str(x) for x in relevantpars])+'\n')
				# if i==5: break


def remove_too_large_to_fit(SUMFILE):
	"""
	PRE   : Takes in a sumfile
	POST  : Will write a sumfile+nottoolarge and and a sumfiletoolarge with those that are small enough to fit inside of the sphere of 8 of radius or not
	"""
	sout = SUMFILE+'nottolarge'
	eout = SUMFILE+'toolarge'
	with open(SUMFILE, 'rb') as r:
		with open(sout, 'wb') as w:
			with open(eout, 'wb') as e:
				for i, line in enumerate(r):
					# if i==3: break
					parts = line.strip().split()
					sphere = parts[8]
					if float(sphere)<7:
						w.write(line)
					else:
						e.write(line)

def remove_3_membered_rings(SUMFILE):
	"""
	PRE   : Takes in a sumfile
	POST  : Will write a sumfile+nothreerings and and a threerings with those that are small enough to fit inside of the sphere of 8 of radius or not
	"""

	sout = SUMFILE+'notthreerings'
	eout = SUMFILE+'threerings'
	with open(SUMFILE, 'rb') as r:
		with open(sout, 'wb') as w:
			with open(eout, 'wb') as e:
				for i, line in enumerate(r):
					# if i==20: break
					parts = line.strip().split()

					smi = parts[0]
					# print line
					try:
						mol = Chem.MolFromSmiles(smi)
						rings = Chem.GetSymmSSSR(mol)
						# print rings
						hasring3 = False
						for els in rings:
							rsize = len(list(els))
							if rsize==3:
								hasring3=True
								e.write(line)
								break
						
						if not hasring3:
							print i, line
							w.write(line)
					except:
						e.write(line)

def remove_allenes(SUMFILE):
	"""
	PRE   : Takes in a sumfile
	POST  : Will write a sumfile+nothreerings and and a threerings with those that are small enough to fit inside of the sphere of 8 of radius or not
	"""

	sout = SUMFILE+'noallenes'
	eout = SUMFILE+'yesallenes'
	with open(SUMFILE, 'rb') as r:
		with open(sout, 'wb') as w:
			with open(eout, 'wb') as e:
				for i, line in enumerate(r):
					# if i==20: break
					parts = line.strip().split()

					smi = parts[0]
					patt = Chem.MolFromSmarts('[#6]=[#6]=[#6]')
					try:
						mol = Chem.MolFromSmiles(smi)

						hasmatch = mol.HasSubstructMatch(patt)
						if hasmatch:
							print i, line
							e.write(line)
						else:
							w.write(line)
					except:
						e.write(line)


def remove_duplicate_from_sum(SUMFILE):
	"""
	PRE : Takes in a sum file
	POST: Writes a new sum file with only the first occurence of each pubchem number in a line
	"""
	sumfilenoduplicate = SUMFILE + '_NODUPLICATE'
	alreadyin = {}
	with open(SUMFILE, 'rb') as r:
		with open(sumfilenoduplicate, 'wb') as w:
			for line in r:
				parts = line.split()
				if parts[1] in alreadyin:
					continue
				else:
					alreadyin[parts[1]] = 0
					w.write(line)

def keep_only_not_rotatable(sumfile):
	"""
	PRE   : takes a sumfile
	POST  : writes out only the molecules that have a single cluster of TFD 
	"""
	sumnotrotatable = sumfile+'_NOTROTATABLE'
	with open(sumfile, 'rb') as r:
		with open(sumnotrotatable, 'wb') as w:
			for line in r:
				parts = line.strip().split()
				if parts[-1] == '1':
					w.write(line)

	pass


def copy_sdfs_for_analysis(SUM_FILE, location_of_sdf='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid/', 
	newlocationofsdf='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid_no_allenes_no3rings_noexclusion_notover200kcal_allsdfallpictures/'):
	"""
	PRE : Takes in a sum file and will take all the files ending in sdf that match the index_pbnumber_file_name in the location_of_sdf to the newlocationofsdf
	POST: SEE up
	"""
	index=getindex()

	with open(SUM_FILE, 'rb') as r:
		for line in r:
			parts = line.strip().split()
			fname = index[parts[1]].split('/')[-1]
			sdflist = glob.glob(location_of_sdf+fname+'*frequency_report')
			for s in sdflist:
				# mol = Chem.MolFromMolFile(s, removeHs=False)
				fname = s.split('/')[-1]
				# Chem.MolToMolFile(mol, newlocationofsdf+parts[1]+fname)
				call('cp {} {}'.format(s, newlocationofsdf+parts[1]+fname), shell=True)

def split_molecule_into_fragments(fname):
	"""
	PRE: Takes in a sdf or mol2 file 
	POST: Writes out as fname-frag1 fname-frag2 etc...
	"""
	if fname[-4:] == '.sdf':
		mol = Chem.MolFromMolFile(fname, removeHs=False)
	elif fname[-4:] == 'mol2':
		mol = Chem.MolFromMol2File(fname, removeHs=False)

	frags = Chem.GetMolFrags(mol, asMols=True)
	print frags
	for i, m in enumerate(frags):
		print i
		if  fname[-4:] == '.sdf':
			Chem.MolToMolFile(m, fname+'-frag{}'.format(i))
		elif  fname[-4:] == 'mol2':
			AllChem.ComputeGasteigerCharges(m)
			Chem.MolToPDBFile(m, fname+'-frag{}.pdb'.format(i))

def associate_complex_charges_to_mol2_files(complexmol2, mol2fraglist):
	"""
	PRE   : Takes in a complex mol2 file with the electrostatic charges and a LIST of mol2 fragments, the mol2 must be given in the same order as the mols are present in the unique complex mol2 file 
	POST  : Will write the charges of the complex to each atom concerned 
	"""
	def make_charge_list(complexmol2):
		"""
		PRE : Takes in a complexmol2 file 
		POST: Will return a dictionary where the keys are the xyz values formatted as a string as '{0:4.4f}{1:4.4f}{2:4.4f}'.format(xflat,yfloat, zfloat) and the value the charge as a float
		"""
		chargelist = []
		with open(complexmol2, 'rb') as r:
			for line in r:
				parts = line.strip().split()
				if len(parts)==9:
					x = float(parts[2])
					y = float(parts[3])
					z = float(parts[4])
					charge = line[-11:-1]
					# chargedic['{0:.0f}{1:.0f}{2:.0f}'.format(x,y,z)] = charge
					chargelist.append(charge)
		return chargelist

	# chargedic = make_dictinary_charge(complexmol2)
	chargelist = make_charge_list(complexmol2)
	for key in chargelist:
		print key
	
	alreadyseen = 0
	for mols in mol2fraglist:
		outname = mols[:-5]+'_withcharges.mol2'
		print outname
		j = 0
		with open(mols, 'rb') as r:
			with open(outname, 'wb') as w:
				for i, line in enumerate(r):
					
					parts = line.strip().split()
					# print len(parts)
					
					if len(parts)==9 and 'TEMP' not in line:
						print j+alreadyseen
						x = float(parts[2])
						y = float(parts[3])
						z = float(parts[4])
						# chargekey = '{0:.0f}{1:.0f}{2:.0f}'.format(x,y,z)
						charge = chargelist[j+alreadyseen]
						newlineend = list(line)
						newlineend[-11:-1] = charge
						line = ''.join(newlineend) 
						j+=1
					w.write(line)
		alreadyseen += j

def compare_RMS(mol1, mol2, includeCB):
	"""
	PRE : Takes in two mols
	POST: Returns the rms variation
	"""
	print mol1
	mol1 = Chem.MolFromMol2File(mol1, removeHs=False)
	mol2 = Chem.MolFromMol2File(mol2, removeHs=False)

	if includeCB:
		nocbfragmol1 = []
		mol1frags = Chem.GetMolFrags(mol1, asMols=True)
		print len(mol1frags)
		for frag in mol1frags:
			if frag.GetNumHeavyAtoms() != 84:
				nocbfragmol1.append(frag)
				print frag.GetNumHeavyAtoms()
		if len(nocbfragmol1) == 2:
			mol1 = Chem.CombineMols(nocbfragmol1[0], nocbfragmol1[1])
		else:
			mol1 = nocbfragmol1[0]

		nocbfragmol2 = []
		mol2frags = Chem.GetMolFrags(mol2, asMols=True)
		for frag in mol2frags:
			if frag.GetNumHeavyAtoms() != 84:
				nocbfragmol2.append(frag)
				print frag.GetNumHeavyAtoms()
		if len(nocbfragmol2) == 2:
			mol2 = Chem.CombineMols(nocbfragmol2[0], nocbfragmol2[1])
		else:
			mol2 = nocbfragmol2[0]






	mol1.AddConformer(mol2.GetConformer(), assignId=True)

	print 'The RMS value is {}'.format(AllChem.GetConformerRMS(mol1, 0, 1, prealigned=False))

	# print mol1.GetConformers()



def pick_random_linesfromfile(fname, nlines):
	"""
	PRE : Takes in a file and writes with lines
	POST: Writes out another file fname+'randomlines{}'.format(nlines)
	"""
	import random
	with open(fname, 'rb') as r:
		lines = r.readlines()

	with open(fname+'_randomelines{}'.format(nlines), 'wb') as w:
		w.writelines(random.sample(lines, nlines))

	return

def add_numbers_to_sumfile(fname):
	"""
	PRE : Takes in a sumfile
	POST: Appends at the end of each line the number of the line starting at 1
	"""
	with open(fname, 'rb') as r:
		with open(fname+'_line_number', 'wb') as w:
			for i, line in enumerate(r):
				line = line.strip()
				w.write(line+'\t{}\n'.format(i+1))


def get_volume(rdkitmol):
	"""
	PRE : Takes in a rdkitmol, the molecule needs to have the hydrogens up 
	POST: Returns its volume based on grid encoding as computed by the rdkit 
	"""
	AllChem.MMFFOptimizeMolecule(rdkitmol)
	return AllChem.ComputeMolVolume(rdkitmol)

def expand_sumfile_with_physical_features(sumfile, location_of_sdf='/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid/'):
	"""
	PRE :Takes in a sumfile formatted as c12c3c(c1cc1c2cccc1)cc1c3cccc1	127262304	131.628796223	119	True	(0,	-1346.148456048429)	-26.3622	7.7353	0	1	1.4179070322	15.246	-125.259267393	-108.595360361
	POST: Will append to that line the xyz of the PCA bounding box, the radius of the projected molecule and the radius of the sphere and the volume and the number of heavy atoms, in that order 
	"""
	index=getindex()
	outf = sumfile+'_extendedphysicalfeatures'
	with open(outf, 'wb') as w:
		with open(sumfile, 'rb') as r: 
			for i, line in enumerate(r):
				print i, line 
				parts = line.strip().split()

				fname = index[parts[1]].split('/')[-1]
				sdflist = glob.glob(location_of_sdf+fname+'*GUEST*.sdf')

				mol = Chem.MolFromMolFile(sdflist[0], removeHs=False)

				MOLBLOCK, xthick, ythick, zthick, rad, sphere_radius = align_mol(Chem.MolToMolBlock(mol))
				nheavy = mol.GetNumHeavyAtoms()
				volume = get_volume(mol)

				parts.extend([xthick, ythick, zthick, rad, sphere_radius, volume, nheavy])

				w.write('\t'.join([str(x) for x in parts])+'\n')


def make_heat_map_input_from_expandedsumfile(sumfile):
	"""
	PRE  : Takes in a expanded sumfile as in expand_sumfile_with_physical_features
	POST : Will output a file to be treated as a input for heat map with just ones and zeros
	"""
	outputheat = sumfile+'_inputforheat'
	allthebools = []
	allconditions = []
	firsttime=True
	with open(sumfile, 'rb') as r:
		with open(outputheat, 'wb') as w:
			for i, line in enumerate(r):
				parts = line.strip().split()
				# print i,line
				# print zip(parts, range(len(parts)))
				# if i==5:break
				datadic = {}
				datadic['volume'] = float(parts[22]) # RDKIT volume
				datadic['xthick'] = float(parts[15]) # Largest  one
				datadic['ythick'] = float(parts[16]) # Middle   one
				datadic['zthick'] = float(parts[17]) # Smallest one
				datadic['heavyatoms'] = float(parts[23])  # Heavy atoms 
				datadic['guestsphere'] = float(parts[21]) # Sphere radius
				datadic['radiusguest'] = float(parts[20][:-1]) # Radius of the guest

				bracketdic = {}
				bracketdic['volume'] = range(72,359, 50) #(72.12, 358.72, 50)
				bracketdic['zthick'] = range(3, 13, 1) #(3.3375, 10.4544, 1)
				bracketdic['ythick'] = range(3, 13, 1)#(3.4,12.82,1)
				bracketdic['heavyatoms'] = range(4,20, 1)
				bracketdic['guestsphere'] = range(1,8,1)#(1.39, 8.1,1)
				bracketdic['radiusguest'] = range(1,8,1)#(1.7,6.4126,1)

				sortedkeys = (sorted(bracketdic.keys()))
				# print sortedkeys
				booleanlist = []
				for key in sortedkeys:
					booleanlist.append(datadic[key]<bracketdic[key][0])
					if firsttime: allconditions.append('{}<{}'.format(key, bracketdic[key][0]))
					for i in range(len(bracketdic[key][0:-1])):
						# booleanlist.append((datadic[key]>=bracketdic[key][i] and datadic[key]<bracketdic[key][i+1], datadic[key], bracketdic[key][i]))
						booleanlist.append((datadic[key]>=bracketdic[key][i] and datadic[key]<bracketdic[key][i+1]))
						if firsttime: allconditions.append('{}<={}<{}'.format(bracketdic[key][i], key, bracketdic[key][i+1]))
					booleanlist.append(datadic[key]>=bracketdic[key][-1])
					if firsttime: allconditions.append('{}<={}'.format(key, bracketdic[key][-1]))
				# print booleanlist
				allthebools.append([str(int(x)) for x in booleanlist])
				firsttime=False
				w.write('\t'.join(allthebools[-1])+'\n')

	# print allthebools
	return allconditions

def plot_heat_map_from_binary_file(binaryfile, conditions):
	"""
	PRE : Takes in a binary file with only ones and zeros per line 
	POST: Plots it as a heatmap
	"""
	with open(binaryfile, 'rb') as r: data = [[int(y) for y in x.strip().split()] for x in r.readlines()]
	data.extend(conditions)
	print data
	import seaborn as sns; sns.set()
	import matplotlib.pyplot as plt
	ax = sns.heatmap(data)
	plt.show()

def test_charge_representation_in_pdb(smi='C1CNC(=[NH+]1)C12CC3CC(C1)CC(C2)(C3)C1=[NH+]CCN1'):
	"""
	:param smi: A SMI test
	:return: produces a pdb file containing the 3D structure of the molecule, just to see how it is reprensented
	"""
	mol = Chem.MolFromSmiles(smi)
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol)
	AllChem.MMFFOptimizeMolecule(mol)
	Chem.MolToPDBFile(mol, '/home/macenrola/Desktop/test_charge.pdb')
	for at in mol.GetAtoms():
		print at.GetExplicitValence(), at.GetFormalCharge(), at.GetSymbol(), at.GetImplicitValence()

def test_read_invalid_obabel_mol(f='/home/macenrola/Desktop/391-orig_guestsPose1.pdb'):
	"""
	:param f: the invalid mol pdbfile
	:return: None
	"""
	mol = Chem.MolFromPDBFile(f, removeHs=False)
	Chem.MolToMolFile(mol, f[:-4]+".sdf")
	return

def find_untreated(flist):
	"""
	:param flist: takes in a list of prmtops such as 276-orig_guestsPose1.prmtop 248-orig_complexPose1-run_tops.sh
	:return: a list of files that correspond to the .sh file that should have run to produce them such as 276-orig_guestsPose1-run_tops.sh
		276-orig_guestsPose1-run_tops.sh
	"""
	with open('/home/macenrola/Desktop/unprocessed', 'wb'): pass
	for f in flist:
		print f[:-12]
		try:
			with open(f[:-12]+'-freq.out', 'rb') as r: pass
		except:
			with open('/home/macenrola/Desktop/unprocessed', 'ab') as a:
				a.write('{}\n'.format(f))

def write_sphere_radii_for_best_complexes(best_complexes_file, format_string_pdbfile='/home/macenrola/Documents/amberconvergedmols/all_pdbs_and_prmtops/{}*_OUT_GUEST*_complexPose{}.pdb'):
	"""
	:param best_complexes_file: a file formatted as follows
		guest best pose; complex best pose; free energy difference; bad;      vdW;     elect;   nonpolar;   genBorn; entropy; guest number ; pubchem number
		7	1	-39.5088353	-42.77	-0.06	-45.34	0.64	-1.38	3.37	3.261	G1	101803327
		2	1	-39.0214055	-42.31	-0.04	-44.9	0.59	-1.34	3.38	3.289	G2	101803326
		1	7	-38.77519885	-42.24	-0.05	-44.8	0.69	-1.34	3.27	3.465	G3	144631
		4	2	-38.50895125	-50.83	2.28	-55.3	3.94	-1.99	0.25	12.321	G4	15196674
	:param format_string_pdbfile: a file location such as /home/macenrola/Documents/amberconvergedmols/all_pdbs_and_prmtops/241xaa_OUT_GUEST158_complexPose0.pdb
	:return: a copy of the initial file with the radii of the sphere that fits the smallest sphere where the best complex would fit
	"""
	from miniball_example_containers import doit
	from mol_ops import get_atoms_coords
	with open(best_complexes_file+'_WITH_RADII', 'wb') as w:
		with open(best_complexes_file, 'rb') as r:
			for i, line in enumerate(r):
				try:
					els = line.strip().split()
					complex_num = els[1]
					pubchem_num = els[11]
					fname = format_string_pdbfile.format(pubchem_num, complex_num)
					truefname = glob.glob(fname)[0]
					mol = Chem.MolFromPDBFile(truefname, removeHs=False)
					# res = align_mol(Chem.MolToMolBlock(mol))
					# print res[-1]
					atm_list, atm_coords = get_atoms_coords(Chem.MolToMolBlock(mol))
					miniball_data = doit(atm_coords)
					radius = miniball_data[-1]**.5
					els.append(radius)
					els.append(truefname)
					w.write('{}\n'.format('\t'.join([str(x) for x in els])))
				except:
					print 'error at step {} at line {}'.format(i, line)

def plot_distributions_for_endo_exo(fin='/home/macenrola/Documents/amberconvergedmols/datamanuscript/sumdic_with_apolar_breakdown-processedfreeenergy-sorted_with_smi_nohighbad.txt_WITH_RADII_with_centroid_diff_pca_atominside',
									threshold=[16.0, 6, 90]):
	"""
	:param fin: gives a file formatted as by  write_sphere_radii_for_best_complexes with binding affinities and complex radii
	:param threshold: the radius threshold to mark a complex as exo, above is exo
	:return: prints a graph with the density of exo and endo according to the binding affinity
	"""
	from scipy.stats import binned_statistic
	import numpy as np
	from scipy.interpolate import interp1d

	import matplotlib.pyplot as plt
	crits = []
	with open(fin, 'rb') as r:
		for line in r:
			els = line.strip().split()
			# print zip(els, range(len(els)))
			crits.append((float(els[2]),float(els[14]), int(els[15]), float(els[16]), float(els[17]), int(els[11])))
			# print crits[-1]
	endo = []
	exo = []
	# return
	for els in crits:
		valid =  [els[1]<threshold[0], els[2]>=threshold[1], els[4]<threshold[2]]
		print els, valid
		if all(valid):
			endo.append(els)
		else:
			exo.append(els)

	bin_val_endo, bin_edges_endo, bin_d_endo = binned_statistic([x[0] for x in endo], [1 for x in endo], 'count', 15)
	bin_width_endo = bin_edges_endo[1]-bin_edges_endo[0]
	bin_centers_endo = bin_edges_endo[:-1]-bin_width_endo/2.0

	bin_val_exo, bin_edges_exo, bin_d_exo =  binned_statistic([x[0] for x in exo], [1 for x in exo], 'count', 15)
	bin_width_exo = bin_edges_exo[1]-bin_edges_exo[0]
	bin_centers_exo = bin_edges_exo[:-1] - bin_width_exo / 2.0

	# print '{}\t{}\t{}\t{}'.format('cendo', 'vendo', 'cexo', 'vexo')
	# for i in range(len(bin_centers_endo)):
	# 	print '{}\t{}\t{}\t{}'.format(bin_centers_endo[i], bin_val_endo[i], bin_centers_exo[i], bin_val_exo[i])

	plt.figure('Endo and Exo distribution')
	fcub = interp1d(bin_centers_endo,bin_val_endo, kind="cubic")
	xnew = np.linspace(bin_centers_endo.min(),bin_centers_endo.max(), 200)
	plt.plot(bin_centers_endo,bin_val_endo, 'o', xnew, fcub(xnew), '-')

	fcub = interp1d(bin_centers_exo,bin_val_exo, kind="cubic")
	xnew = np.linspace(bin_centers_exo.min(),bin_centers_exo.max(), 200)
	plt.plot(bin_centers_exo, bin_val_exo, 'x', xnew, fcub(xnew), '--')
	plt.legend(['Endo points','Endo interpolation','Exo points','Exo interpolation'])
	plt.ylabel('Count (#)')
	plt.xlabel('Binding affinity (kcal/mol)')
	plt.show()
	print '{}\t{}\t{}\t{}'.format('cendo', 'vendo', 'cexo', 'vexo')
	for i in range(len(bin_centers_endo)):
		print '{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}'.format(bin_centers_endo[i], bin_val_endo[i], bin_centers_exo[i], bin_val_exo[i])

def plot_volume_distrib(fin='/home/macenrola/Documents/amberconvergedmols/Pubchem_C/pdbs/PUBCHEM_ONLY_C_only_19_atomswith_volume'):
	"""
	:param fin: takes in a file formatted as c12c(c(cc3c1cccc3)C)cc1c(c2)cccc1	9422	275.90
	:return: plots a binned version of the volume distribution
	"""
	from scipy.stats import binned_statistic
	import numpy as np
	from scipy.interpolate import interp1d
	import matplotlib.pyplot as plt
	vols = []
	with open(fin, 'rb') as r:
		for i,line in enumerate(r):
			els = line.strip().split()
			print i, els
			vols.append(float(els[-1]))

	print vols
	bin_val, bin_edges, bin_d = binned_statistic(vols, [1 for x in vols], 'count', 15)
	bin_width = bin_edges[1]-bin_edges[0]
	bin_centers = bin_edges[:-1]-bin_width/2.0

	plt.figure('Vols distribution')
	fcub = interp1d(bin_centers,bin_val, kind="cubic")
	xnew = np.linspace(bin_centers.min(),bin_centers.max(), 200)
	plt.plot(bin_centers,bin_val, 'o', xnew, fcub(xnew), '-')
	plt.show()


def copy_files_for_first_mols(fin='/home/macenrola/Documents/amberconvergedmols/datamanuscript/sumdic_with_apolar_breakdown-processedfreeenergy-sorted_with_smi_nohighbad.txt_WITH_RADII_with_centroid_diff_pca_atominside'):
	"""
	:param fin: a file that contain the number of best guests and complexes 
	:return: will copy the files of the best guests and complexes to another file
	"""
	with open(fin, 'rb') as r:
		for i,line in enumerate(r):
			if i ==50:
				break
			gnum, cnum = line.split()[:2]
			pbnum = line.split()[11]
			rank = line.split()[10]
			cfile = line.split()[13]
			gfile = cfile[:-16]+'guestsPose{}.pdb'.format(gnum)
			guest = Chem.MolFromPDBFile(gfile, removeHs=False)
			complex = Chem.MolFromPDBFile(cfile, removeHs=False)
			print guest, complex
			Chem.MolToMolFile(guest, '/home/macenrola/Documents/amberconvergedmols/top50/RANK_{}_PBNUM_{}_GUEST.sdf'.format(rank, pbnum))
			Chem.MolToMolFile(complex, '/home/macenrola/Documents/amberconvergedmols/top50/RANK_{}_PBNUM_{}_COMPLEX.sdf'.format(rank, pbnum))
if __name__ == '__main__':       
	import glob
	# plot_distributions_for_endo_exo()

	copy_files_for_first_mols()

	# exes('/home/macenrola/Documents/amberconvergedmols/datamanuscript/sumdic_with_apolar_breakdown-processedfreeenergy-sorted_with_smi_nohighbad.txt')
	# plot_volume_distrib()
	# test_charge_representation_in_pdb()
	# test_read_invalid_obabel_mol()
	# flist = glob.glob('/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/prmtops_freqs/*sh')
	# find_untreated(flist)
	# conditions = make_heat_map_input_from_expandedsumfile('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal_extendedphysicalfeatures')
	# plot_heat_map_from_binary_file('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal_extendedphysicalfeatures_inputforheat', conditions)
	# add_numbers_to_sumfile('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/details-amber-solv-entropy')
	# expand_sumfile_with_physical_features('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal')
	# make_heat_map_input_from_expandedsumfile('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal_extendedphysicalfeatures')
	# remove_too_large_to_fit('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal')
	# pick_random_linesfromfile('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/DFTrigids/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal', 30)
	# produce_picture_grid_from_SUM('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/sortednottolargenotthreeringsnoallenesbelow200kcalnottolarge/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcalnottolarge')
	# compare_RMS('/home/macenrola/Thesis/ESCAPETIME/CB7_DBOA_1H_in0_TS-C2H4_wB97XD_6-31Gs_1c.mol2', '/home/macenrola/Thesis/ESCAPETIME/CB7_DBOA_1H_in0_PD-C2H4_wB97XD_6-31Gs.mol2', iscb)
	# compare_RMS('/home/macenrola/Thesis/ESCAPETIME/CB7_DBOA_1H_in0_TS-N2_wB97XD_6-31Gs_1b.mol2', '/home/macenrola/Thesis/ESCAPETIME/CB7_DBOA_1H_in0_PD-N2_wB97XD_6-31Gs.mol2', iscb)
	# compare_RMS('/home/macenrola/Thesis/ESCAPETIME/CB7_DBOA_1H_in0_TS-NH3_wB97XD_6-31Gs_1b.sdf', '/home/macenrola/Thesis/ESCAPETIME/CB7_DBOA_1H_in0_PD-NH3_wB97XD_6-31Gs.mol2', iscb)
	# associate_complex_charges_to_mol2_files('/home/macenrola/Thesis/ESCAPETIME/NH3-no-charge/CB7_DBOA_1H_in0_PD-NH3_wB97XD_6-31Gs.mol2', 
	# 	['/home/macenrola/Thesis/ESCAPETIME/NH3-no-charge/CB7_DBOA_1H_in0_PD-NH3_wB97XD_6-31Gs_frag0_SDF2PDB.mol2',
	# 	'/home/macenrola/Thesis/ESCAPETIME/NH3-no-charge/CB7_DBOA_1H_in0_PD-NH3_wB97XD_6-31Gs_frag1_SDF2PDB.mol2',
	# 	'/home/macenrola/Thesis/ESCAPETIME/NH3-no-charge/CB7_DBOA_1H_in0_PD-NH3_wB97XD_6-31Gs_frag2_SDF2PDB.mol2'])
	# remove_allenes('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreerings')
	# remove_3_membered_rings('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolarge')
	# copy_sdfs_for_analysis('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal_extendedphysicalfeatures')
	# extract_best_guest_and_complex_from_tars('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE')
	# add_solvation_to_sumlist('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE')
	# add_entropy_to_sumlist('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation')
	# produce_picture_grid_from_SUM('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE_reformatted_sortednottolargenotthreeringsnoallenesbelow200kcal-20first')
	# reformat_after_solvation_entropy('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE_NOTROTATABLE_with_solvation_with_entropy_NODUPLICATE')
	# extract_host_guest_from_SDF('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/rigid/xzzazpq_OUT.sdf', 2, 3, '/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/test')
	# keep_only_not_rotatable('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS_with_uncorrelated_conformers_NODUPLICATE')
	# remove_duplicate_from_sum('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS_with_uncorrelated_conformers')
	# reorganizes_sum_file_list('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_NAMED_HYDROCARBONS')
	# get_no_rotatable_bonds('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/2000firstSUM_with_uncorrelated_conformers-rigid')
	# flist = glob.glob('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/butina/*SUM')
	# for i,f in enumerate(flist):
	# 	if i<165: continue
	# 	print i, f
	# 	extract_best_host_guest_sdf(f)
		# print f
	# get_no_rotatable_bonds('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/ALL_SUMS')
	# produce_picture_grid_from_SUM('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/pictures/amber_force')
	# reformat_after_solvation_entropy('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/2000firstSUM_with_uncorrelated_conformers-rigid_with_solvation_with_entropy')
	# add_entropy_to_sumlist('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/2000firstSUM_with_uncorrelated_conformers-rigid_with_solvation')
	# dummy_sort('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/2000firstSUM_with_uncorrelated_conformers')
	# mv_guest_from_a_to_b('xzzefzn.tar')
	# fname = '/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/butina/xzzbmyo_OUT.sdf'
	# confs_uncorrelated = return_uncorrelated_conformations(fname)
	# extract_host_guest_from_SDF(fname, guestnumber=confs_uncorrelated, hostnumber=())
	# reorganizes_sum_file_list('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/HYDROCARBONS_ALL_CONFS/NAMED_SUM_ALL_HYDRO')
	##### Test align conformers
	# 
	# mol = Chem.MolFromSmiles('CCCCCCCCC')
	# mol = Chem.AddHs(mol)
	# AllChem.EmbedMultipleConfs(mol, 50, AllChem.ETKDG())
	# AllChem.MMFFOptimizeMoleculeConfs(mol)
	# align_all_conformers_of_molecule(mol)


	# read_text_as_repr('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd_aligned.sdf')
	# get_the_hydrocarbons_not_treated('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/sorted_ALLCONFHYDROSUM', '/media/macenrola/backup_ext4/PUBCHEM/ALLCANSPUBCHEM_unique_no_stereo_below20_only_carbs.can')
	# produce_picture_grid_from_SUM('/home/macenrola/Thesis/VIBRATIONS/hosts_salts_stripped.can')
	# extract_host_guest_from_SDF('/home/macenrola/Thesis/VIBRATIONS/test_vibrations/four_armed_benzene.sdf', return_uncorrelated_conformations('/home/macenrola/Thesis/VIBRATIONS/test_vibrations/four_armed_benzene.sdf'), ())
	# keep_only_the_carbons('/media/macenrola/backup_ext4/PUBCHEM/ALLCANSPUBCHEM_unique_no_stereo_below20')
	# extract_best_host_guest_sdf('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/highbinding/xzzehqn_SUM')
	# print get_rotable_bond_number('c1ccccc1')
	# produce_picture_grid_from_SUM('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/pictures/2000_70kto72k')
	# fnames = [
	# '/home/macenrola/Thesis/XAAXAB_SUM/allconformations/all_conformations/48first2000',
	# '/home/macenrola/Thesis/XAAXAB_SUM/allconformations/all_conformations/48first2000-neutral',
	# '/home/macenrola/Thesis/XAAXAB_SUM/allconformations/all_conformations/48last2000',
	# '/home/macenrola/Thesis/XAAXAB_SUM/allconformations/all_conformations/48last2000-neutral',
	# '/home/macenrola/Thesis/XAAXAB_SUM/allconformations/all_conformations/48middle2000',
	# '/home/macenrola/Thesis/XAAXAB_SUM/allconformations/all_conformations/48middle2000-neutral',
	# ]
	# for fname in fnames:
	# 	show_conf_vs_opti(fname)
	# sdflist = glob.glob('/media/macenrola/fat/ONLY_GUESTS/INCLUSION_EXCLUSION/*_ONLY_GUESTS.sdf') + glob.glob('/media/macenrola/fat/ONLY_GUESTS/XAA_INCLUSION_EXCLUSION/*_ONLY_GUESTS.sdf')
	# get_the_guests('/home/macenrola/Thesis/XAAXAB_SUM/all_newguests', sdflist)

	# add_new_binding_energy('/home/macenrola/Thesis/XAAXAB_SUM/MOLDIC_GUESTS', '/home/macenrola/Thesis/XAAXAB_SUM/MOLDIC_GUESTS_2' ,'/home/macenrola/Thesis/XAAXAB_SUM/xaaxabSUM-bestE-sorted')
	# get_pubchem_list('/home/macenrola/Thesis/XAAXAB_SUM/48first2000-neutral')
	# for i, fname in zip([1,2,3], ['INCLUSION_EXCLUSION', 'INCLUSION_EXCLUSION', 'XAA_INCLUSION_EXCLUSION']):
	# 	get_SDFs_from_canlist(get_pubchem_list('/media/macenrola/fat/all_guests_of_interest'), glob.glob('/media/macenrola/fat/{}/{}/x*.tar'.format(i, fname)), '/media/macenrola/fat/All_of_interest_{}.sdf'.format(i))
	# get_infos_from_sdffilelist(glob.glob('/home/macenrola/Thesis/XAAXAB_SUM/EXTRACTED_2000/*.sdf'))
	# keep_mols_that_match_pbnb_bindingE(glob.glob('/home/macenrola/Thesis/XAAXAB_SUM/EXTRACTED_2000/*.sdf'), '/home/macenrola/Thesis/XAAXAB_SUM/EXTRACTED_2000/molecules_to_keep_from_extracted', '/home/macenrola/Thesis/XAAXAB_SUM/EXTRACTED_2000/best_complexes')
	# extract_guests_from_complex('/media/macenrola/fat/xaaanc/xaaanc_OUT.sdf', '/media/macenrola/fat/xaaanc/xaaanc_OUT_ONLY_GUESTS.sdf')
	# dummy_sort('/home/macenrola/Thesis/XAAXAB_SUM/xaaxabSUM-bestE-sorted_new_binding_sorted')
	# produce_picture_grid('/home/macenrola/Thesis/PHARMACY_BAD_TASTING/bad_taste_medecine_no_salts_other_diclofenac.can (2)/bad_taste_medecine_no_salts_other_diclofenac.can_OUT.sdf_ONLY_BEST.sdf')
	# flist = glob.glob('/media/macenrola/FAT/Pubchem_DB_SMILES/*Compound*')
	# for f in flist:
	# 	keep_only_hydrocarbons(f)
	# keep_only_hydrocarbons('/media/macenrola/FAT/PUBCHEM_ONLY_C')
	# get_energies_calibration('/home/macenrola/Thesis/Calibration_force_field/MMFF94_hypervalent_rdkit.sdf')
	# get_ref_energies()
	# keep_guests('/home/macenrola/Thesis/GENERATE_FINGERPRINTS/xzamc-min_20_only_mmff94_OUT.sdf')'/media/macenrola/fat/3/XAA_INCLUSION_EXCLUSION'
	# swap_id_smiles('/home/macenrola/Thesis/PHARMACY_BAD_TASTING/bad_taste_medecine')
	# split_sdf_file_individually('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/individual_good_sdfs/cyclodextrin_guests_OUT_ONLY_GOOD.sdf')
	# fix_PDB_spacing('/home/macenrola/Thesis/AMBER/converge_with_amber/COMPLEX.pdb')
	# dummy_sort()
	# make_pdb_complex_with_named_residues('/home/macenrola/Thesis/AMBER/converge_with_amber/242.sdf','/home/macenrola/Thesis/AMBER/converge_with_amber/CB_candidate.sdf')
	# keep_only_best_SDFs('/home/macenrola/Thesis/PHARMACY_BAD_TASTING/bad_taste_medecine_no_salts_other_diclofenac.can (2)/bad_taste_medecine_no_salts_other_diclofenac.can_OUT.sdf', '/home/macenrola/Thesis/PHARMACY_BAD_TASTING/bad_taste_medecine_no_salts_other_diclofenac.can (2)/bad_taste_medecine_no_salts_other_diclofenac.can_SUM')
	# converge_AMBER_GAFF('/home/macenrola/Thesis/AMBER/converge_with_amber/test.pdb')
	# keep_only_hydrocarbons('/home/macenrola/Thesis/hydrocarbons/below_20_pubchem_no_undesirable_only_mmff94_atoms_no_salts.can')
	# keep_only_best_SDFs('/home/macenrola/Thesis/XAAXAB_SUM/reran2000/48reran/48middle2000-neutral/48middle2000-neutral_all.sdf', '/home/macenrola/Thesis/XAAXAB_SUM/reran2000/48reran/48middle2000-neutral/48middle2000-neutral_SUM')
	# complexfile = '/home/macenrola/Thesis/drugbank/drugbank_treated/drugbank_SUM_best_complex_names'
	# guestfile = '/home/macenrola/Thesis/drugbank/drugbank_treated/drugbank_SUM_best_guests_names'
	# with open(complexfile, 'rb') as r:
	# 	pubchem_list = [x.strip() for x in r.readlines()]
	# with open(guestfile, 'rb') as r:
	# 	pubchem_list.extend([x.strip() for x in r.readlines()])
	# tarlist = glob.glob('/home/macenrola/Thesis/drugbank/drugbank_treated/*tar')
	# get_SDFs_from_canlist(pubchem_list, tarlist)
	# get_list_of_best_complexes('/home/macenrola/Thesis/drugbank/drugbank_treated/drugbank_SUM')
	# get_list_of_best_guests('/home/macenrola/Thesis/drugbank/drugbank_treated/drugbank_SUM')
	# dummy_sort()
	# sort_ZINCs15('/home/macenrola/Gaussian_stuff/sanity_check_inclusion_rdkit/all_carbon_energy_sort_split/rerun_all_carbon_zinc_compatible/all_carbon_zinc15/all_carbon_zinc15')
	# assign_ZINCS('/home/macenrola/Gaussian_stuff/sanity_check_inclusion_rdkit/all_carbon_energy_sort_split/all_carbon_energy_sort', '/home/macenrola/Gaussian_stuff/sanity_check_inclusion_rdkit/all_carbon_energy_sort_split/rerun_all_carbon_zinc_compatible/all_carbon_zinc15/all_carbon_zinc15_sorted')
	# split_sdf_file_individually('/home/macenrola/Thesis/XAAXAB_SUM/EXTRACTED_2000/Best_complexes/best_complexes')
	# for fname in glob.glob('/home/macenrola/Thesis/XAAXAB_SUM/pictures_of_2000/48*'):
	# 	produce_picture_grid_from_SUM(fname)
	# produce_picture_grid_from_SUM('/home/macenrola/Gaussian_stuff/validation_RDKIT_C_+2/random_hydrocarbons/Force-Field-Check/SUM')
	# keep_only_MMFF94_compatible('/home/macenrola/Thesis/Calibration_force_field/SUM_BEFORE_ERRORS_XAA_XAB_FULLY_NEUTRAL_SORTED')
	# complexesInfo = cPickle.load(open('/home/macenrola/Thesis/XAAXAB_SUM/EXTRACTED_2000/molecules_to_keep_from_extracted', 'rb'))
	# for k in complexesInfo:
	# 	print k, complexesInfo[k]
