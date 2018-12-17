import os
import glob
paramfile = 'paramfile'
script_name = 'CYCLO_VIBRATIONS.sh'
base_python_script = "vibrations.py"

canon_string = [
"#!/bin/bash -l",
"#$ -S /bin/bash",
"#$ -cwd",
"#$ -l h_rt=12:00:00",
"#$ -l mem=5G",
"#$ -l tmpfs=20G",
"#$ -N {}".format(script_name),
]
library_string = [
			"export LD_LIBRARY_PATH=/home/uccahcl/bin/lib:$LD_LIBRARY_PATH",
			"export PATH=/home/uccahcl/anaconda2/bin:$PATH",
			"export PATH='/home/uccahcl/bin/bin':$PATH"
			]

extension_change = 'for file in *.html; do mv "$file" "$(basename "$file" .html).txt"; done'
python_path_venv = "/home/uccahcl/python/venv/bin/python"
python_path_conda = "~/anaconda2/envs/my-rdkit-env/bin/python"
canon_string = '\n'.join(canon_string)
library_string = '\n'.join(library_string)
# base_python_script = "rdkit_experimentation.py"
def generate_list_of_sh_scripts_cylinders(list_of_files):
	"""takes in input a list of files to be treated using the find_cylinder.py script, 
	then for each of them produces a .sh file to be submitted as a stand-alone job"""
	for files in list_of_files:
		with open(files+'_job.sh', 'wb') as w:
			w.write(canon_string +'\n')
			w.write(' '.join([python_path_conda, base_python_script, "'%s'"%files, '\n']))

def generate_array_job(parameter_files, name, path):
	with open(prefix+paramfile, 'wb') as w:
		for i, fnames in enumerate(parameter_files):
			w.write('{0:04}\t{1}\n'.format(i+1,fnames))
	with open(prefix+name, 'wb') as w:
		w.write(canon_string + '\n')
		w.write("#$ -t 1-{0}\n".format(len(parameter_files)))
		s = ["number=$SGE_TASK_ID",
			"paramfile=$(pwd)/{0}".format(paramfile),
			"index=$(sed -n ${number}p $paramfile | awk '{print $1}')",
			"variable1=$(sed -n ${number}p $paramfile | awk '{print $2}')"]
		print s
		w.write('\n'.join(s) +'\n')
		w.write(library_string + '\n')
		w.write(' '.join([python_path_conda, base_python_script, '$variable1', '\n']))
		# w.write('obabel -ican {} -O {} --gen3D\n'.format('$variable1', '${variable1}.xyz'))
		# w.write('obabel {} -O {} --filter "atoms < 20" -d\n'.format('$variable1', '$variable1.can' ))
		# w.write('tar -cf {}.tar {}'.format('$variable1.can'))
		# w.write('tar -cf {}.tar {}'.format('$variable1.smi'))


# def submit_list_of_sh_scripts(list_of_scripts):
# 	with open('thrower.sh', 'wb') as w:
# 		w.write("export PATH=$PATH:$(pwd)\n")
# 		w.write("chmod +x *.sh\n")
# 		for scripts in list_of_scripts:
# 			w.write("qsub %s\n"%scripts)
# def converter_job(in_file, out_file, prefix):
# 	conversion_string = [
# 	"#!/bin/bash -l",
# 	"#$ -S /bin/bash",
# 	"#$ -cwd",
# 	"#$ -l h_rt=10:00:00",
# 	"#$ -l mem=5G",
# 	"#$ -l tmpfs=5G",
# 	"export LD_LIBRARY_PATH=/home/uccahcl/bin/lib:$LD_LIBRARY_PATH",
# 	"export PATH='/home/uccahcl/bin/bin':$PATH"
# 	]
# 	with open(prefix+'converter_job.sh', 'wb') as w:
# 		w.write('\n'.join(conversion_string)+'\n')
# 		# w.write('obabel {} -O {} --unique /nostereo\n'.format(in_file, out_file))
# 		w.write('obabel {} -O {} --filter "atoms < 20" -d\n'.format(in_file, out_file))
# 		w.write('tar -cf {0}.tar {0}\n'.format(in_file))
# 		w.write('tar -cf {0}.tar {0}'.format(out_file))

# prefix = './smiles_search/filtration_all/withHetHet_unique_nonNS/'
prefix = '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/individual_good_sdfs/'			
list_of_files_for_scripts = glob.glob(prefix+'*sdf')
list_of_files_for_scripts = [elems.split("/")[-1] for elems in list_of_files_for_scripts]
# converter_job('zinc15_merged.smi', 'zinc15_merged.can', prefix)
# print len(list_of_files_for_scripts), list_of_files_for_scripts
generate_array_job(list_of_files_for_scripts, script_name, prefix)

# generate_list_of_sh_scripts_cylinders(list_of_files_for_scripts)
# list_of_files_to_throw = glob.glob('*job.sh')
# print list_of_files_to_throw
# list_of_files_to_throw = [elems.split("\\")[1] for elems in list_of_files_to_throw]
# print list_of_files_to_throw
# submit_list_of_sh_scripts(list_of_files_to_throw)
