import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

def keep_only_low_E_from_SDF(SDFname):
	"""From the 4 complexes (2 in - 2 out) keeps the best in and the best out"""
	supp = Chem.SDMolSupplier(SDFname, removeHs=False)
	best_complexes = {}

	for mol in supp:
		cname = mol.GetProp('_Name').strip().split('-')[0]
		n = int(mol.GetProp('COMPLEX_NUMBER'))
		print cname
		bestInKey = cname+'-bestIn'
		bestOutKey = cname+'-bestOut'
		if n <=2:
			if bestInKey not in best_complexes:
				best_complexes[bestInKey] = mol
			elif float(mol.GetProp('BINDING_ENERGY')) < float(best_complexes[bestInKey].GetProp('BINDING_ENERGY')):
				best_complexes[bestInKey] = mol

		else:
			if float(mol.GetProp('BINDING_ENERGY')) > 100: continue
			if bestOutKey not in best_complexes:
				best_complexes[bestOutKey] = mol
			elif float(mol.GetProp('BINDING_ENERGY')) < float(best_complexes[bestOutKey].GetProp('BINDING_ENERGY')):
				best_complexes[bestOutKey] = mol
	w = Chem.SDWriter(SDFname+'out.sdf')
	for key in best_complexes:
		w.write(best_complexes[key])
	w.close()


def create_G09(title, xyz, path):
	"""This code will take each of the molecules in the file and create a gaussian file with it"""
	title  = title.replace(' ', '').replace('\t', '').split('.')[0]
	# title.replace('\t', '')
	g09_backbone = [
	'%Chk={0}.chk',
	'%NProcShared=12',
	'%mem=2gb',
	'#n wB97XD/6-31G* Opt=Tight Int=Ultrafine SCF=QC',
	# '%rwf=DNT_wB97XD_6-31Gd_Raman.rwf',
	# '%NoSave',
	# '%Chk=DNT_wB97XD_6-31Gd_Raman.chk',
	# '%NProcShared=64',
	# '%mem=2gb',
	# '#n wB97XD/6-31G* Opt Freq=Raman',
	# #n HF/6-31G* Opt Freq=Raman
	# '# Opt=Restart',
 	'',
 	'{0}',
 	'',
 	'0 1',
 	'{1}'
 	'\n',
 	'\n'
 	# '{0}.wfn',
 	# '\n'
 	]
 	g09_backbone = '\n'.join(g09_backbone).format(title, '\n'.join(xyz))
 	with open(path+'/{}_wB97XD_631Gd.com'.format(title), 'wb') as w:
		w.write(g09_backbone)

def parse_xyz(fname):
	"""
	PRE: Takes in a xyz file obtained by translating a rdkit sdf file into an obabel xyz
	POST: 
	"""
	path = '/'.join(fname.split('/')[:-1])
	XYZs = {}
	with open(fname, 'rb') as r:
		for line in r:
			try:
				n = int(line.strip())
				print n
				pb_nb = r.next().strip()
				print pb_nb
				XYZs[pb_nb] = []
				for i in range(n):
					l = r.next().strip()
					print l
					XYZs[pb_nb].append([x for x in l.split(' ') if x])
			except Exception as e:
				print e
	for key in XYZs:
		# print key, XYZs[key]
		for i, at in enumerate(XYZs[key]):
			a, b, c, d = at
			bs, cs, ds = '','',''
			if float(b)>=0: bs=' '
			if float(c)>=0: cs=' '
			if float(d)>=0: ds=' '
			XYZs[key][i] = a +' '*9+ bs + b + ' '*7 + cs + c + ' '*7 + ds + d
			# print XYZs[key][i]
		create_G09(key, XYZs[key], path)

def get_errors_or_unprocessed(folderpath):
	import glob
	ferr = glob.glob(folderpath + '*ERR')
	ferr = [e.split('/')[-1] for e in ferr]
	fsdf = glob.glob(folderpath + '*OUT.sdf')
	fsdf = [e.split('/')[-1] for e in fsdf]
	fout = 'xaa_all_failures'
	with open(folderpath+fout, 'wb'): pass
	for f in ferr:
		rad = f[:-4]
		if rad+'_OUT.sdf' in fsdf:
			print rad +'!'*30
			with open(folderpath+fout, 'ab') as a:
				with open(folderpath+rad, 'rb') as r:
					for line in r:
						a.write(line)
		else:
			print f
			with open(folderpath+fout, 'ab') as a:
				with open(folderpath+f, 'rb') as r:
					for line in r:
						a.write(line)			




def make_gaussian_arrayjob(arrayjobfilename):
	"""
	PRE : Takes in an arrayjobfilename and assumes that in the same folder there are com files 
	POST: Creates the paramfilegaussian and arrayjobfile name 
	"""
	import glob
	parse_xyz(arrayjobfilename)
	prefix = '/'.join(arrayjobfilename.split('/')[:-1])+'/'
	paramfile = 'paramfile'
	paramfilelist = glob.glob(prefix+'*com')
	paramfilelist = [x.split('/')[-1] for x in paramfilelist]
	script_name = 'adamantanelike.sh'
	canon_string = [
	"#!/bin/bash -l",
	"#$ -S /bin/bash",
	"#$ -cwd",
	"#$ -l h_rt=24:0:0",
	"#$ -l mem=2G",
	"#$ -l tmpfs=10G",
	"#$ -N {}".format(script_name[:-3]),
	"#$ -pe smp 12"
	]

	gaussian_string = [
	#'module load gaussian/g09-d01/pgi-2015.7',
	#'source $g09root/g09/bsd/g09.profile',
	'module load gaussian/g16-a03/pgi-2016.5',
	''
	'mkdir -p $GAUSS_SCRDIR',
	# Run g09 job
	'echo "GAUSS_SCRDIR = $GAUSS_SCRDIR"',
	'echo "Running: g16 < $g09infile > $g09outfile"',
	'export GAUSS_MEMDEF=76MW',
	'g16 < $g09infile > $g09outfile'
	]



	def generate_array_job():
		with open(prefix+paramfile, 'wb') as w:
			for i, fnames in enumerate(paramfilelist):
				w.write('{0:04}\t{1}\n'.format(i+1,fnames))
		with open(prefix+script_name, 'wb') as w:
			w.write('\n'.join(canon_string) + '\n')
			w.write("#$ -t 1-{0}\n".format(len(paramfilelist)))
			s = ["number=$SGE_TASK_ID",
				"paramfile=$(pwd)/{0}".format(paramfile),
				"index=$(sed -n ${number}p $paramfile | awk '{print $1}')",
				"g09infile=$(sed -n ${number}p $paramfile | awk '{print $2}')",
				"g09outfile=$g09infile'_OUT.out'"]
			print s
			w.write('\n'.join(s) +'\n')
			w.write('\n'.join(gaussian_string))

	generate_array_job()

if __name__ == '__main__':
	# SDFname = '/home/macenrola/Thesis/openbabel_for_failed_rdkit/sample_complex.sdfout.sdf'
	# supp = Chem.SDMolSupplier(SDFname, removeHs=False)
	# for mol in supp:
		# print mol.GetProp('_Name'), mol.GetProp('BINDING_ENERGY')
	# parse_xyz('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/20first/mypart10-20/all_firstamber.xyz')
	# get_errors_or_unprocessed('/home/macenrola/Thesis/XAA_INCLUSION_EXCLUSION/')
	# keep_only_low_E_from_SDF('/home/macenrola/Thesis/openbabel_for_failed_rdkit/sample_complex.sdf')
	make_gaussian_arrayjob('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/adamantane/adamantanelike.xyz')