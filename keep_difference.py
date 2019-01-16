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
			.
			.
			.
			vibrational:             866.574        310.409        289.623


			the pubchem number, pose number and type of molecule (guest or complex) is obtained from the file name
	"""
	breakdownline = ""
	total = ""
	vibrational = ""
	type = ""
	posenum = ""
	number = get_number_from_fname(fname)
	# if fname.split("_")[-1][0] == "g":
	if 'guests' in fname:
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
			elif l[:12] == "vibrational:":
				vibrational = l
			else:
				# print l[:10]
				pass
	return int(number), type, int(posenum), [float(x) for x in breakdownlines[indexofff].split()[2:]], [float(x) for x in total.split()[1:]], [float(x) for x in vibrational.split()[1:]]

def get_number_from_fname(fname="/home/macenrola/Desktop/12883016xzzabfs_OUT_GUEST17_complexPose0-freq.nab-freqreport"):
	"""
	:param fname: the file name is formatted like this
				  /home/macenrola/Desktop/12864396xzzabef_OUT_GUEST5_guestsPose1-freq.nab-freqreport
	:return: the pubchem number contained in the file name
	"""
	return int(fname.split('/')[-1].split('-')[0])

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
			for example:
			241 {'COMPLEX': {0: ([-63.66, 187.38, -64.69, -134.7, 3.42, -55.07, 0.0635], [639.04, 270.186, 344.779], [704.337, 264.229, 256.047]), 1: ([-64.81, 187.45, -65.85, -134.33, 3.43, -55.51, 0.0169], [637.817, 270.212, 345.442], [704.28, 264.254, 256.718]), 2: ([-64.76, 187.48, -65.83, -134.32, 3.42, -55.51, 0.108], [636.687, 266.243, 331.443], [703.094, 260.285, 242.72]), 3: ([-64.82, 187.45, -65.85, -134.33, 3.42, -55.51, 0.0171], [637.817, 270.212, 345.441], [704.28, 264.255, 256.717]), 4: ([-64.82, 187.45, -65.85, -134.33, 3.42, -55.51, 0.0177], [638.411, 272.197, 355.796], [704.874, 266.239, 267.071]), 5: ([nan, nan, nan, nan, 0.0, nan, nan], [nan, nan, nan], [nan, nan, 0.0]), 6: ([-61.52, 187.31, -63.27, -135.3, 3.56, -53.82, 0.0223], [639.863, 266.155, 328.389], [703.159, 260.197, 239.65]), 7: ([101358087626.72, 101358088121.49, -0.09, -6.41, 0.0, -488.27, 374000.0], [], []), 8: ([-61.57, 187.32, -63.32, -135.32, 3.56, -53.8, 0.044], [639.796, 266.159, 328.855], [703.148, 260.202, 240.116])}, 'GUEST': {0: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 1: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 2: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 3: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 4: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 5: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 6: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 7: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226]), 8: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86], [64.582, 10.542, 4.226])}}


	"""
	import cPickle as pk
	presort = set()
	for i,f in enumerate(listoffiles):
		try:
			number = get_number_from_fname(f)
			presort.add(number)
		except:
			print "error with step {}/{}".format(i,f)
			pass
		if i%1000==0:
			print "step {}".format(i)
	sumdic = dict.fromkeys(presort)

	for i,f in enumerate(listoffiles):
		if i%1000==0:
			print "step {}/{}".format(i,f)
		try:
			number, type, posenum, breakdown, total, vibrational = parse_amber_report(f, indexofff)
		except:
			print "error with step {}/{}".format(i,f)
			continue
		if sumdic[number] == None:
			sumdic[number] = {}
		if type not in sumdic[number]:
			sumdic[number][type] = {}
		sumdic[number][type][posenum] = (breakdown, total, vibrational)

	with open("/home/macenrola/Desktop/sumdic_with_apolar_known_guests_tweaked_params", "wb") as handle:
		pk.dump(sumdic, handle)

	# with open("/home/macenrola/Desktop/sumdic", "rb") as r:
	# 	print pk.load(r)

def treatdic(sumdic_file = "/home/macenrola/Desktop/sumdic_with_apolar_known_guests"):
	"""

	:param sumdic_file: takes in a dic formatted as in make_summary_dic_of_numbers_and_poses and produces the binding energies
	one line looks like this
	key 241: {'COMPLEX': {0: ([-63.66, 187.38, -64.69, -134.7, 3.42, -55.07, 0.0635], [639.04, 270.186, 344.779]), 1: ([-64.81, 187.45, -65.85, -134.33, 3.43, -55.51, 0.0169], [637.817, 270.212, 345.442]), 2: ([-64.76, 187.48, -65.83, -134.32, 3.42, -55.51, 0.108], [636.687, 266.243, 331.443]), 3: ([-64.82, 187.45, -65.85, -134.33, 3.42, -55.51, 0.0171], [637.817, 270.212, 345.441]), 4: ([-64.82, 187.45, -65.85, -134.33, 3.42, -55.51, 0.0177], [638.411, 272.197, 355.796]), 5: ([nan, nan, nan, nan, 0.0, nan, nan], [nan, nan, nan]), 6: ([-61.52, 187.31, -63.27, -135.3, 3.56, -53.82, 0.0223], [639.863, 266.155, 328.389]), 7: ([101358087626.72, 101358088121.49, -0.09, -6.41, 0.0, -488.27, 374000.0], []), 8: ([-61.57, 187.32, -63.32, -135.32, 3.56, -53.8, 0.044], [639.796, 266.159, 328.855])}, 'GUEST': {0: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 1: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 2: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 3: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 4: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 5: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 6: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 7: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86]), 8: ([4.02, 0.13, 3.14, 0.54, 1.05, -0.85, 0.016], [69.324, 16.499, 68.86])}}

	:return: a file with the binding energies
	"""
	import cPickle
	outfile = sumdic_file + "-processedfreeenergy"
	with open(sumdic_file, 'rb') as r:
		sumdic = cPickle.load(r)
	with open(outfile, "wb") as w:
		for i, k in enumerate(sorted(sumdic)[:]):
			print i, k, sumdic[k]
			# print k, sumdic[k]
			# return
			if i%1000==0: print "step {}".format(i)
			minEguest = 1e8
			minEcomplex = 1e8
			gbestpose = -1
			cbestpose = -1
			# print k, sumdic[k]
			if sumdic[k] == None: continue
			# if i==10: break

			for p in sumdic[k]["GUEST"]:
				breakdown, total, vibrational = sumdic[k]["GUEST"][p]
				if breakdown == [] or total == [] or vibrational == []: continue
				tempE = breakdown[0] - total[-1]*298.15/1000
				# tempE = breakdown[0] - total[-1] * 298.15 / 1000
				# print tempE
				if tempE<minEguest:
					minEguest=tempE
					gbestpose = p

			for p in sumdic[k]["COMPLEX"]:
				breakdown, total, vibrational = sumdic[k]["COMPLEX"][p]
				if breakdown == [] or total == [] or vibrational == []: continue
				tempE = breakdown[0] - total[-1] * 298.15/1000
				# tempE = breakdown[0] - total[-1]*298.15/1000
				# print tempE
				if tempE<minEcomplex:
					minEcomplex=tempE
					cbestpose = p

			# print minEguest, minEcomplex
			CB_breakdown = [-47.78,    187.19,    -47.83,   -135.56,      3.45,    -55.04, -296.926*298.15/1000] # Standard using ff energy and overall entropy
			if gbestpose == -1 or cbestpose == -1: continue
			G_breakdown, G_total, G_vibrational = sumdic[k]["GUEST"][gbestpose]
			G_breakdown[-1] = (-G_total[-1]*298.15/1000)
			C_breakdown, C_total, C_vibrational = sumdic[k]["COMPLEX"][cbestpose]
			C_breakdown[-1] = (-C_total[-1]*298.15/1000)
			# print len(CB_breakdown), CB_breakdown
			# print len(G_breakdown), G_breakdown
			# print len(C_breakdown), C_breakdown
			#      iter    Total       bad      vdW     elect   nonpolar   genBorn      frms
			# ff:    60    -47.78    187.19    -47.83   -135.56      3.45    -55.04  1.00e-03
			#			freq		E			Cv				S
			#			cm ** -1	kcal / mol	cal / mol - K	cal / mol - K
			# Total:                585.346		243.982			296.926
			w.write("{}\t{}\t{}\t{}\t{}\r\n".format(k, gbestpose, cbestpose, sum(C_breakdown[1:]) - sum(G_breakdown[1:]) - sum(CB_breakdown[1:]), '\t'.join(['{:6.3f}'.format(x[0]-x[1]-x[2]) for x in zip(C_breakdown, G_breakdown, CB_breakdown)])))


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

def get_vina_affinities_report(listoffiles, fout='/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/RES_VINA_KNOWN_HOSTS'):
	"""
	:param listoffiles: takes in a list of files formatted as: /home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/pdbs/158-orig.pdbqt-LOG
	:return: prints a file fout with the results formatted as number \tab affinity in kcal/mol
	"""
	with open(fout, 'wb') as w:
		for f in listoffiles:
			num = get_number_from_fname(f)
			affinity = parse_vina_log(f)
			w.write('{}\t{}\n'.format(num, affinity))


def parse_vina_log(fname):
	"""
	:param fname: a vina log formatted as: /home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/pdbs/140-orig.pdbqt-LOG
	:return: the affinity of the best pose from the report
	"""
	with open(fname, 'rb') as r:
		for line in r:
			print line[:4]
			if line[:4]=='   1':
				affinity = line.strip().split()[1]
				return affinity
	raise Exception('NO AFFINITY FOUND')

def make_vina_bc_ourmethod(vina_report, bcfile, methodsummary, fout):
	"""
	:param vina_report: /home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/RES_VINA_KNOWN_HOSTS a vina report formatted '{}\t{}\n'.format(number, affinity)
	:param bcfile:/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/res_host_with_BC.can formatted [CH-]1C=CC=C1	140	dES:91.54	dFF:-12.74	BC:6
	:param methodsummary: /home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/sumdic_with_apolar_known_guests-processedfreeenergy 140	0	8	-24.761762	-3.450	 0.000	 0.010	 0.040	-3.450	-0.050	-21.312
	:return: a merge of the three files including a conversion of bc to affinities
	"""
	import math
	vinadic = {}
	bcdic = {}
	methoddic ={}

	with open(vina_report, 'rb') as r:
		for line in r:
			num, affinity = line.strip().split()
			vinadic[int(num)] = float(affinity)

	with open(bcfile, 'rb') as r:
		for line in r:
			els = line.strip().split()
			num = els[1]
			smi = els[0]
			print els[-1][3:]
			bc = -float(els[-1][3:])/1000/4.18*8.3145*298.15*2.30258509299
			bcdic[int(num)] = (bc, smi)

	with open(methodsummary, 'rb') as r:
		for line in r:
			els = line.strip().split()
			num = (els[0])
			affinity = els[3]
			methoddic[int(num)] = affinity

	with open(fout, 'wb') as w:
		w.write('{}\t{}\t{}\t{}\t{}\n'.format('number', 'smi', 'ref', 'vina', 'method'))
		for k in sorted(bcdic.keys()):
			print k
			try:
				w.write('{}\t{}\t{}\t{}\t{}\n'.format(k, bcdic[k][1], bcdic[k][0], vinadic[k], methoddic[k]))
			except: print '{} posed problem'.format(k)
# def remove_high_bad(basefile):
# 	"""
# 	:param basefile: takes in a file formatted as
# 	91413520	0	2	-44.6241473	-51.320	-18.690	-34.560	-0.510	-1.490	 3.950	 6.696	c12c(cc3=CC=CC=CC=c3c2)C=CC1
# 	:return: returns a copy of the file where the bad energy (here -18.690) is above 10, it means probably something was misfolded
# 	"""
#
# 	with open(basefile, 'rb') as r:
# 		with open(basefile+'_nohighbad', 'wb') as w:
# 			for i,line in enumerate(r):
# 				# if i==10: break
# 				els =
#
# 				print els, els[5]
# 				if abs(float(els[5])) > 10.0:
# 					continue
# 				else:
# 					w.write(line)

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
	# treatdic()
	# make_parralel_script(sorted(glob.glob('/home/macenrola/Documents/amberconvergedmols/all_pdbs_and_prmtops/*nab')))
	# keepdifference('/home/macenrola/Documents/amberconvergedmols/allminus1stround',
	# 			   '/home/macenrola/Documents/amberconvergedmols/2ndroundtar',
	# 			   '/home/macenrola/Documents/amberconvergedmols/allminus1stminus2nd')
	# print parse_amber_report()
	# make_summary_dic_of_numbers_and_poses(glob.glob("/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/prmtops_freqs_tweak_parm/*.nab.out-freqreport"))
	# # print get_number_from_fname()
	# flist=glob.glob('/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/pdbs/*-orig.pdbqt-LOG')
	# get_vina_affinities_report(flist)


	make_vina_bc_ourmethod('/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/RES_VINA_KNOWN_HOSTS',
						   '/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/res_host_with_BC.can',
						   '/home/macenrola/Desktop/sumdic_with_apolar_known_guests_tweaked_params-processedfreeenergy',
						   '/home/macenrola/Documents/amberconvergedmols/VinaVsOurMethodVsExp/results/all_together_report_tweakedParams')

	# generate_3d("/home/macenrola/Documents/docked_for_data_analysis/500krandomless25hvatoms")
	# treatdic('/home/macenrola/Desktop/sumdic_with_apolar_known_guests_tweaked_params')
	# number_lines('/home/macenrola/Documents/docked_for_data_analysis/400k-500klist_pdbqt')
	# addSmi('/home/macenrola/Documents/docked_for_data_analysis/sumdic_with_apolar_breakdown-processedfreeenergy-sorted',
		   # '/home/macenrola/Documents/docked_for_data_analysis/pubchem_smis')
	# remove_high_bad('/home/macenrola/Documents/docked_for_data_analysis/sumdic_with_apolar_breakdown-processedfreeenergy-sorted_with_smi')