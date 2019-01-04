def keepdifference(allparams, treatedparams, fout):
	alreadytreated = set()
	with open(treatedparams, 'rb') as r:
		for line in r:
			alreadytreated.add(line.split('.tar')[0][2:])
	difference = []
	all = set()
	with open(allparams, 'rb') as r:
		for line in r:
			candidate = line.strip().split()[-1]
			all.add(candidate)
			if candidate not in alreadytreated:
				difference.append(candidate)
	print len(difference)
	print len(alreadytreated)
	print len(all)
	with open(fout, 'wb') as w:
		for i, els in enumerate(sorted(difference)):
			w.write('{:04d}\t{}\n'.format(i, els))

def make_parralel_script(listoffile):
	"""
	:param listPattern: a pattern for a list of scripts such as /home/macenrola/Documents/Thesis/ScreeningManuscriptFinalData/Generation3Ddata/PDBs_and_docked/*guestsPose0-run_tops.sh
	matching /home/macenrola/Documents/Thesis/ScreeningManuscriptFinalData/Generation3Ddata/PDBs_and_docked/241xaa_OUT_GUEST158_guestsPose0-run_tops.sh
	:return: a single script to run them all in parallel
	"""
	with open('/home/macenrola/Documents/amberconvergedmols/all_pdbs_and_prmtops/parallel_converge_scriptby1000.sh', 'wb') as w:
		k = 0
		for i,f in enumerate(listoffile):
			converge_script ='echo {0} && nab {0} -o {0}.out && ./{0}.out > {0}-freqreport && rm {0}.out  & \n'
			convert_script = "echo {0}; sed 's/\/home\/macenrola\/Desktop\/all_pdbs_and_prmtops\//\/home\/macenrola\/Documents\/amberconvergedmols\/all_pdbs_and_prmtops\//g' {0} -i; " \
							 "sed 's/gb=0, dielc=1.0, gbsa=0/gb=1, dielc=1.0, gbsa=1/g' {0} -i & \n"
			fname = f.split('/')[-1]
			try:
				with open(f+'-freqreport', 'rb') as r:
					print "this file exists {}, I'm skippint it".format(fname)
					print fname
			except:
				k+=1
				w.write(converge_script.format(fname))
				print fname
				if k % 200 == 0:
					k=0
					w.write('wait\n')



def parse_amber_report(fname, indexofff):
	"""
	:param fname: the file name of the output of a nab minimization with gb=1 pbsa=1 gaff and normal mode analysis
				  the file name is formatted like this
				  /home/macenrola/Desktop/12864396xzzabef_OUT_GUEST5_guestsPose1-freq.nab-freqreport
				  indexofff is the index of the ff line, for some reason, the final one doesn't include the apolar contribution
				  if indexofff is -1, the last ff line is used and the nonpolar contribution is omitted, if -3 is used, the contribution is included
	:return: a tuple of the values Entropy, Sovlation energy, Force field energy etc...
			 the breakdown split is formatted as:
			       iter    Total       bad      vdW     elect   nonpolar   genBorn      frms
			ff:     1     64.90     61.96      4.70     -0.77      0.00     -0.99  1.55e-02

			the total entropy is formatted as:
			            cm**-1       kcal/mol       cal/mol-K     cal/mol-K
			Total:                   197.487         31.212         83.509

			the pubchem number, pose number and type of molecule (guest or complex) is obtained from the file name
	"""
	breakdownline = ""
	total = ""
	type = ""
	posenum = ""
	number = fname.split('/')[-1].split('x')[0]
	if fname.split("_")[-1][0] == "g":
		type = "GUEST"
	else:
		type = "COMPLEX"
	posenum = fname.split("Pose")[-1][0]
	breakdownlines = []
	with open(fname, 'rb') as r:
		for l in r:
			# if l[:10] == "ff:     1 ":
			if l[:3] == "ff:":
				# breakdownline = l
				breakdownlines.append(l)
				# print breakdownline

			elif l[:10] == "Total:    ":
				total = l
			else:
				# print l[:10]
				pass
	return int(number), type, int(posenum), [float(x) for x in breakdownlines[indexofff].split()[2:]], [float(x) for x in total.split()[1:]]

def get_number_from_fname(fname="/home/macenrola/Desktop/12883016xzzabfs_OUT_GUEST17_complexPose0-freq.nab-freqreport"):
	"""
	:param fname: the file name is formatted like this
				  /home/macenrola/Desktop/12864396xzzabef_OUT_GUEST5_guestsPose1-freq.nab-freqreport
	:return: the pubchem number contained in the file name
	"""
	return int(fname.split('/')[-1].split('x')[0])

def make_summary_dic_of_numbers_and_poses(listoffiles, indexofff=-3):
	"""
	:param listoffiles: a list of files formatted as
						/home/macenrola/Desktop/12864396xzzabef_OUT_GUEST5_guestsPose1-freq.nab-freqreport
	:return: a dictionnary that contains one key for each pubchem number, that links to a dictionary containing each pose for GUEST and COMPLEX which and the energy breakdown and total entropy are stored as lists

						       iter    Total       bad      vdW     elect   nonpolar   genBorn      frms
						ff:     1     64.90     61.96      4.70     -0.77      0.00     -0.99  1.55e-02

			the total entropy is formatted as:
			            cm**-1       kcal/mol       cal/mol-K     cal/mol-K
			Total:                   197.487         31.212         83.509
	"""
	import cPickle as pk
	presort = set()
	for i,f in enumerate(listoffiles):
		number = get_number_from_fname(f)
		presort.add(number)
		if i%1000==0:
			print "step {}".format(i)
	sumdic = dict.fromkeys(presort)

	for i,f in enumerate(listoffiles):
		if i%1000==0:
			print "step {}/{}".format(i,f)
		try:
			number, type, posenum, breakdown, total = parse_amber_report(f, indexofff)
		except:
			print "error with step {}/{}".format(i,f)
			continue
		if sumdic[number] == None:
			sumdic[number] = {}
		if type not in sumdic[number]:
			sumdic[number][type] = {}
		sumdic[number][type][posenum] = (breakdown, total)

	with open("/home/macenrola/Desktop/sumdic_with_apolar", "wb") as handle:
		pk.dump(sumdic, handle)

	# with open("/home/macenrola/Desktop/sumdic", "rb") as r:
	# 	print pk.load(r)

def treatdic(sumdic_file = "/home/macenrola/Desktop/sumdic_with_apolar_breakdown"):
	"""

	:param sumdic_file: takes in a dic formatted as in make_summary_dic_of_numbers_and_poses and produces the binding energies
	:return: a file with the binding energies
	"""
	import cPickle
	outfile = sumdic_file + "-processedfreeenergy"
	with open(sumdic_file, 'rb') as r:
		sumdic = cPickle.load(r)
	with open(outfile, "wb") as w:
		for i, k in enumerate(sorted(sumdic)[:]):
			if i%1000==0: print "step {}".format(i)
			minEguest = 1e8
			minEcomplex = 1e8
			gbestpose = -1
			cbestpose = -1
			print k, sumdic[k]
			if sumdic[k] == None: continue
			# if i==10: break
			for p in sumdic[k]["GUEST"]:
				breakdown, total = sumdic[k]["GUEST"][p]
				if breakdown == [] or total == []: continue
				tempE = breakdown[0] - total[-1]*298.15/1000
				# print tempE
				if tempE<minEguest:
					minEguest=tempE
					gbestpose = p
			for p in sumdic[k]["COMPLEX"]:
				breakdown, total = sumdic[k]["COMPLEX"][p]
				if breakdown == [] or total == []: continue
				tempE = breakdown[0] - total[-1]*298.15/1000
				# print tempE
				if tempE<minEcomplex:
					minEcomplex=tempE
					cbestpose = p

			# print minEguest, minEcomplex
			CB_breakdown = [-47.78,    187.19,    -47.83,   -135.56,      3.45,    -55.04, -296.926*298.15/1000]
			if gbestpose == -1 or cbestpose == -1: continue
			G_breakdown, G_total = sumdic[k]["GUEST"][gbestpose]
			G_breakdown[-1] = (-G_total[-1]*298.15/1000)
			C_breakdown, C_total = sumdic[k]["COMPLEX"][cbestpose]
			C_breakdown[-1] = (-C_total[-1]*298.15/1000)
			# print len(CB_breakdown), CB_breakdown
			# print len(G_breakdown), G_breakdown
			# print len(C_breakdown), C_breakdown
			#      iter    Total       bad      vdW     elect   nonpolar   genBorn      frms
			# ff:    60    -47.78    187.19    -47.83   -135.56      3.45    -55.04  1.00e-03
			#			freq		E			Cv				S
			#			cm ** -1	kcal / mol	cal / mol - K	cal / mol - K
			# Total:                585.346		243.982			296.926
			w.write("{}\t{}\t{}\t{}\t{}\r\n".format(k, gbestpose, cbestpose, minEcomplex - minEguest - (-47.78-296.926*298.15/1000), '\t'.join(['{:6.3f}'.format(x[0]-x[1]-x[2]) for x in zip(C_breakdown, G_breakdown, CB_breakdown)])))


def generate_3d(SMIile):
	"""
	:param SDFile: Takes in a SMIFile
	:return: many sdf File with their 3D data
	"""
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	# w = Chem.SDWriter('/home/macenrola/Documents/docked_for_data_analysis/500krandomless25hvatoms.sdf')
	with open('/home/macenrola/Documents/docked_for_data_analysis/500krandomless25hvatoms_ERRORS', 'wb') as errs:
		with open(SMIile, 'rb') as r:
			for i, line in enumerate(r):
				print line
				try:
					str, nb = line.strip().split()
					mol = Chem.MolFromSmiles(str)
					tempmol = Chem.AddHs(mol)
					AllChem.EmbedMolecule(tempmol)
					AllChem.MMFFOptimizeMolecule(tempmol)
					tempmol.SetProp('_Name', nb)
					Chem.MolToPDBFile(tempmol, '/home/macenrola/Documents/docked_for_data_analysis/500k_pdbs/{}.pdb'.format(nb))
					# w.write(tempmol)
				except:
					errs.write(line)
					print 'ERROR LINE:{}'.format(line)
				if i%10==0:
					print '{}th step lowl'.format(i)


def split_bigsdf_to_individual_pdbs(SDFile):
	"""
	:param SDFile: Takes in a sdf file
	:return: produces an individual file as pdb for each entry
	"""
	supp = Chem.SDMolSupplier(SDFile)


def addSmi(basefile, pubsmi):
	"""
	:param basefile: takes in a file with each line starting with a pubchem number and separeted from the rest by tab
	:param pubsmi: takes in a file starting with smis and separated from their pubchem number
	:return: produces a third file which adds the smi associated with the pubchem number on each line from the basefile
	"""
	# step 1 builds a dic
	predic = set()
	with open(basefile, 'rb') as r:
		for i, line in enumerate(r):
			if i % 1000 == 0: print i
			pbnbr = line.strip().split('\t')[0]
			predic.add(int(pbnbr))

	smidic = dict.fromkeys(predic)

	with open(pubsmi, 'rb') as r:
		for line in r:
			smi, pbnbr = line.strip().split('\t')
			if int(pbnbr) in smidic:
				smidic[int(pbnbr)] = smi
				print smi, pbnbr

	with open(basefile, 'rb') as r:
		with open(basefile+'_with_smi', 'wb') as w:
			for line in r:
				pbnbr = int(line.split('\t')[0])
				w.write(line.strip()+'\t{}\n'.format(smidic[pbnbr]))


def remove_high_bad(basefile):
	"""
	:param basefile: takes in a file formatted as
	91413520	0	2	-44.6241473	-51.320	-18.690	-34.560	-0.510	-1.490	 3.950	 6.696	c12c(cc3=CC=CC=CC=c3c2)C=CC1
	:return: returns a copy of the file where the bad energy (here -18.690) is above 10, it means probably something was misfolded
	"""

	with open(basefile, 'rb') as r:
		with open(basefile+'_nohighbad', 'wb') as w:
			for i,line in enumerate(r):
				# if i==10: break
				els =

				print els, els[5]
				if abs(float(els[5])) > 10.0:
					continue
				else:
					w.write(line)

def number_lines(f):
	"""
	:param f: takes in a text file
	:return: another text file "-numbered" with lines formatted as 000001
	"""
	with open(f,'rb') as r:
		with open(f+'-numbered', 'wb') as w:
			for i, line in enumerate(r):
				w.write("{0:05d}\t{1}".format(i, line))

if __name__ == "__main__":

	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem
	import glob
	# make_parralel_script(sorted(glob.glob('/home/macenrola/Documents/amberconvergedmols/all_pdbs_and_prmtops/*nab')))
	# keepdifference('/home/macenrola/Documents/amberconvergedmols/allminus1stround',
	# 			   '/home/macenrola/Documents/amberconvergedmols/2ndroundtar',
	# 			   '/home/macenrola/Documents/amberconvergedmols/allminus1stminus2nd')
	# print parse_amber_report()
	# make_summary_dic_of_numbers_and_poses(glob.glob("/home/macenrola/Documents/amberconvergedmols/all_pdbs_and_prmtops/*report"))
	# print get_number_from_fname()
	# treatdic()
	# generate_3d("/home/macenrola/Documents/docked_for_data_analysis/500krandomless25hvatoms")
	# treatdic()
	# number_lines('/home/macenrola/Documents/docked_for_data_analysis/400k-500klist_pdbqt')
	# addSmi('/home/macenrola/Documents/docked_for_data_analysis/sumdic_with_apolar_breakdown-processedfreeenergy-sorted',
		   # '/home/macenrola/Documents/docked_for_data_analysis/pubchem_smis')
	# remove_high_bad('/home/macenrola/Documents/docked_for_data_analysis/sumdic_with_apolar_breakdown-processedfreeenergy-sorted_with_smi')