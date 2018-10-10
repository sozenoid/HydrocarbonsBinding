import rdkit 
from rdkit import Chem
from rdkit.Chem import AllChem
import scipy
from sklearn.decomposition import PCA
import numpy as np
from subprocess import call
from SmallestEnclosingCircle_CASTING import returnCircleAsTuple
from miniball_example_containers import doit

def get_atoms_coords(RDKIT_BLOCK):
	"""Takes as input an RDKIT BLOCK and returns a list of atoms with a numpy array containing the coordinates"""
	RDKIT_BLOCK = RDKIT_BLOCK.split('\n')
	atm_number = int(RDKIT_BLOCK[3][:3])
	RDKIT_BLOCK = [x.split() for x in RDKIT_BLOCK]
	atm_list = []
	coords_array = np.zeros([atm_number, 3], dtype=float)
	for i, line in enumerate(RDKIT_BLOCK[4:4+atm_number]):
	    coords_atm = line
	    atm_list.append(coords_atm[3])
	    coords_array[i, :] = coords_atm[:3]
	return atm_list, coords_array


def align_mol(RDKIT_BLOCK):
	"""
	PRE: Takes as input a RDKIT_BLOCK with valid 3D coordinates (basically a SDF file in text format)
	POST: Returns a RDKIT_BLOCK with the 3D coordinates rotated so that the main axis of the molecule (PCA wise) is aligned with the z axis and all coordinates are 0 centered
	"""
	def align_xyz(list_atoms, coord_matrix):
		"""This method uses a PCA on the atoms of the molecules as represented by the nuclear coordinates as points 
		and their electron cloud as a sphere of points around those nuclear coordinates
		The function returnCircleAsTuple is imported via a pybind11 generated module from c++ code"""
		sphere_points = get_sphere_points(list_atoms, coord_matrix)
		total_matrix = np.concatenate((coord_matrix, sphere_points),axis = 0)
		pca = PCA(n_components = 3)

		transform = pca.fit_transform(total_matrix)
		transform_coord = pca.transform(coord_matrix) 

		point_cloud = zip(transform.T[1][:], transform.T[2][:])
		xthick = np.max(transform.T[0][:]) - np.min(transform.T[0][:])
		ythick = np.max(transform.T[1][:]) - np.min(transform.T[1][:])
		zthick = np.max(transform.T[2][:]) - np.min(transform.T[2][:])
		# rad = make_circle(point_cloud)
		rad = returnCircleAsTuple(point_cloud)

		transform_coord_centered = transform_coord.copy()
		transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
		transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
		transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])
		return transform_coord_centered, xthick, ythick, zthick, rad

	def get_sphere_points(list_atoms, coord_matrix):
		"""Find the thinnest cylinder dimesions where the molecule could fit add clouds of points around the atoms before the PCA (pi/3 in every direction with atomic radius)
		x = r*sin(theta)*cos(theta)
		y = r*sin(theta)*sin(theta)        theta is inclination from top (0 to pi) and phi is azimuth (0 to 2pi)
		z = r*cos(theta)
		This method uses list comprehension and return all the points representing the spheres around the atoms
		"""
		index_of_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06} 
		angle_nbs = 6 # for angle_nbs = 3 there will be points spaced by pi/3 so three sections in inclinations and six in azimutal 
		sphere_points = []
		for i in xrange(len(list_atoms)):
			radius = index_of_radii[list_atoms[i]]
			top_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(0)]
			bottom_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(np.pi)]
			sphere_points.append(top_sphere_point)
			sphere_points.append(bottom_sphere_point)


		new_sphere_points = [[coord_matrix[i][0]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.cos(azymuth_angle), coord_matrix[i][1]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.sin(azymuth_angle), coord_matrix[i][2]+index_of_radii[list_atoms[i]]*np.cos(inclination_angle)] 
							for azymuth_angle in np.linspace(0, np.pi *2, angle_nbs*2) for inclination_angle in np.linspace(0, np.pi, angle_nbs + 1)[1:-1] for i in range(len(list_atoms))]
		sphere_points.extend(new_sphere_points)
		return sphere_points

	def generate_mol_from_MDL(RDKIT_BLOCK_IN, coord_matrix):
		"""Will write the MDL of the mol file then replace the xyz coordinates from the coord_matrix"""
		RDKIT_BLOCK = [x+'\n' for x in RDKIT_BLOCK_IN.split('\n')]
		atm_number = int(RDKIT_BLOCK[3][:3])
		for i in range(0,atm_number):
			j = i+4
			RDKIT_BLOCK[j] = RDKIT_BLOCK[j].split()
			RDKIT_BLOCK[j][:3] = coord_matrix[i, :]
			RDKIT_BLOCK[j] = (' '*(3+int(np.sign(RDKIT_BLOCK[j][0])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][0])+
					' '*(3+int(np.sign(RDKIT_BLOCK[j][1])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][1])+
					' '*(3+int(np.sign(RDKIT_BLOCK[j][2])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][2])+
					' {}   '.format(RDKIT_BLOCK[j][3]) + '  '.join(RDKIT_BLOCK[j][4:]) + '\n'
					)

		RDKIT_BLOCK_OUT = ''.join(RDKIT_BLOCK)

		return RDKIT_BLOCK_OUT  

	def set_miniball_data(atom_list, atom_coords):
		atom_coords_as_tuples = [tuple(x) for x in atom_coords]
		miniball_data = doit(atom_coords_as_tuples)
		miniball_data[-1] = miniball_data[-1]**.5
		return miniball_data[-1]


	atom_coords = get_atoms_coords(RDKIT_BLOCK)
	transformed_coords, xthick, ythick, zthick, rad = align_xyz(atom_coords[0], atom_coords[1])
	sphere_radius = set_miniball_data(atom_coords[0], atom_coords[1])
	return generate_mol_from_MDL(RDKIT_BLOCK, transformed_coords), xthick, ythick, zthick, rad, sphere_radius

def get_CB_BLOCK():
	"""
	PRE: -
	POST: Returns the RDKIT_BLOCK of CB7 as a string
	"""
	return 'CB7\n     RDKit          3D\n\n126147  0  0  0  0  0  0  0  0999 V2000\n    1.1730   -4.4176    2.9511 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0750   -4.6349    3.6688 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0765   -3.4911    4.7539 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1723   -2.7865    4.5031 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9191   -3.3583    3.4714 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3124   -4.4064    2.9369 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0522   -3.3353    3.4421 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3130   -2.7704    4.4841 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0691   -3.0468    3.1437 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1939   -3.0104    3.0984 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1854   -0.4434    5.2692 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0591   -0.0113    5.8877 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0441    1.5519    5.6863 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2039    1.7887    4.9742 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9401    0.6201    4.7698 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2994   -0.4233    5.2450 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0305    0.6539    4.7401 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2819    1.8110    4.9656 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0894    0.5504    4.3212 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1741    0.6051    4.2737 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2527   -3.8136   -3.7077 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0124   -3.8818   -4.4671 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0312   -2.5540   -5.3160 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2779   -1.9153   -4.9189 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0091   -2.6783   -4.0054 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2331   -3.7867   -3.7189 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9639   -2.6375   -4.0281 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2071   -1.8908   -4.9333 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1545   -2.4399   -3.6082 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1102   -2.3785   -3.6467 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3096    4.4630   -2.8373 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0731    5.2303   -2.7776 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0598    5.7883   -1.3030 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2961    5.2698   -0.7356 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0425    4.5245   -1.6504 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1763    4.4888   -2.8686 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9303    4.5503   -1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1898    5.2750   -0.7590 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1845    4.0859   -1.4747 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0830    4.1303   -1.5454 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2378    3.8620    3.6411 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0016    4.6239    3.6943 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0176    5.4457    2.3491 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2631    5.0331    1.7185 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9935    4.1349    2.4992 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2476    3.8752    3.6162 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9787    4.1602    2.4607 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2224    5.0509    1.6961 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1395    3.7339    2.2694 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1253    3.7745    2.2088 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8611   -0.7617   -5.5808 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8794   -1.7513    5.3501 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7603   -0.7297   -5.6064 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7412   -1.7795    5.3807 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1982   -5.2981    0.6479 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0498   -5.9266    0.2399 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0334   -5.7852   -1.3302 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2138   -5.0813   -1.5932 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9522   -4.8379   -0.4333 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2882   -5.2657    0.6262 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0214   -4.8068   -0.4705 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2708   -5.0748   -1.6171 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1034   -4.3917   -0.3828 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1644   -4.3384   -0.4391 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7455   -5.3457    1.9916 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8706   -5.3200    1.9558 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8157   -4.8707   -2.9471 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7943   -4.9008   -2.9121 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7863    3.0997    4.7479 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8328    3.1315    4.7180 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8453    5.6668    0.5485 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7668    5.6860    0.5091 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1731    2.6985   -4.5639 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0766    2.5733   -5.2995 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0673    1.0736   -5.7843 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1835    0.5544   -5.2506 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9264    1.5226   -4.5717 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3114    2.6728   -4.5336 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0441    1.4838   -4.5307 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3002    0.5315   -5.2313 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0788    1.3907   -4.1449 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1843    1.3310   -4.0791 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7319    3.9561   -4.0995 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8841    3.9189   -4.0557 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0817   -5.6424    4.1014 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0872   -3.8701    5.7828 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0702   -0.3130    6.9419 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0393    2.1106    6.6299 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0054   -4.7893   -5.0823 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0345   -2.7314   -6.3981 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0901    6.0195   -3.5388 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0617    6.8837   -1.2538 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0074    5.2580    4.5891 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0242    6.5312    2.5046 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7799   -0.8998   -6.6666 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9121   -0.7356   -5.2918 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9301   -1.6557    5.0751 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7985   -2.0876    6.3912 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8160   -0.6826   -5.3367 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6624   -0.8715   -6.6906 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7974   -1.7001    5.1217 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6395   -2.1199    6.4192 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0657   -6.9691    0.5804 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0296   -6.7488   -1.8530 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8037   -5.0954    1.9086 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6356   -6.3657    2.3808 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9239   -5.0599    1.8472 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7805   -6.3449    2.3374 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8721   -4.6322   -2.8198 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7133   -5.8028   -3.5171 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8522   -4.6814   -2.7612 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6872   -5.8367   -3.4748 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8415    2.9390    4.5246 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6891    3.6916    5.6670 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8862    2.9907    4.4754 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7409    3.7289    5.6339 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9003    5.3920    0.5370 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7470    6.7545    0.6523 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6574    6.7728    0.6099 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8247    5.4251    0.4716 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0949    3.2963   -6.1235 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0742    0.9704   -6.8759 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7919    3.7784   -3.9142 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6170    4.7060   -4.8933 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8027    4.6722   -4.8487 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9349    3.7196   -3.8436 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1 65  1  0\n  2  3  1  0\n  2  6  1  0\n  2 85  1  0\n  3  4  1  0\n  3 86  1  0\n  4  5  1  0\n  4 54  1  0\n  5  1  1  0\n  6  7  1  0\n  6 66  1  0\n  7  8  1  0\n  8  3  1  0\n  9  5  2  0\n 10  7  2  0\n 11 12  1  0\n 12 13  1  0\n 12 16  1  0\n 12 87  1  0\n 13 14  1  0\n 13 88  1  0\n 14 15  1  0\n 14 69  1  0\n 15 11  1  0\n 16 17  1  0\n 17 18  1  0\n 18 13  1  0\n 18 70  1  0\n 19 15  2  0\n 20 17  2  0\n 21 22  1  0\n 21 68  1  0\n 22 23  1  0\n 22 26  1  0\n 22 89  1  0\n 23 24  1  0\n 23 90  1  0\n 24 25  1  0\n 25 21  1  0\n 26 27  1  0\n 26 67  1  0\n 27 28  1  0\n 28 23  1  0\n 28 53  1  0\n 29 25  2  0\n 30 27  2  0\n 31 32  1  0\n 31 84  1  0\n 32 33  1  0\n 32 36  1  0\n 32 91  1  0\n 33 34  1  0\n 33 92  1  0\n 34 35  1  0\n 34 71  1  0\n 35 31  1  0\n 36 37  1  0\n 36 83  1  0\n 37 38  1  0\n 38 33  1  0\n 38 72  1  0\n 39 35  2  0\n 40 37  2  0\n 41 42  1  0\n 41 69  1  0\n 42 43  1  0\n 42 46  1  0\n 42 93  1  0\n 43 44  1  0\n 43 94  1  0\n 44 45  1  0\n 45 41  1  0\n 46 47  1  0\n 46 70  1  0\n 47 48  1  0\n 48 43  1  0\n 48 72  1  0\n 49 45  2  0\n 50 47  2  0\n 51 24  1  0\n 51 80  1  0\n 51 95  1  0\n 51 96  1  0\n 52  8  1  0\n 52 16  1  0\n 52 97  1  0\n 52 98  1  0\n 53 99  1  0\n 53100  1  0\n 54 11  1  0\n 54101  1  0\n 54102  1  0\n 55 56  1  0\n 55 65  1  0\n 56 57  1  0\n 56 60  1  0\n 56103  1  0\n 57 58  1  0\n 57104  1  0\n 58 59  1  0\n 58 68  1  0\n 59 55  1  0\n 60 61  1  0\n 60 66  1  0\n 61 62  1  0\n 62 57  1  0\n 62 67  1  0\n 63 59  2  0\n 64 61  2  0\n 65105  1  0\n 65106  1  0\n 66107  1  0\n 66108  1  0\n 67109  1  0\n 67110  1  0\n 68111  1  0\n 68112  1  0\n 69113  1  0\n 69114  1  0\n 70115  1  0\n 70116  1  0\n 71 44  1  0\n 71117  1  0\n 71118  1  0\n 72119  1  0\n 72120  1  0\n 73 74  1  0\n 74 75  1  0\n 74 78  1  0\n 74121  1  0\n 75 76  1  0\n 75122  1  0\n 76 53  1  0\n 76 77  1  0\n 77 73  1  0\n 77 81  2  0\n 78 79  1  0\n 79 80  1  0\n 80 75  1  0\n 82 79  2  0\n 83 73  1  0\n 83123  1  0\n 83124  1  0\n 84 78  1  0\n 84125  1  0\n 84126  1  0\nM  END\n'

def get_CD_BLOCK():
	"""
	PRE: - 
	POST: Returns the RDKIT_BLOCK of CD-beta as a string
	"""
	return 'beta-cd.pdb\n     RDKit          3D\n\n147154  0  0  1  0  0  0  0  0999 V2000\n    2.4563    6.1625    0.1033 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6286   -5.5350   -0.0520 C   0  0  2  0  0  0  0  0  0  0  0  0\n    3.5580   -5.3006   -1.2367 C   0  0  1  0  0  0  0  0  0  0  0  0\n    3.9837   -3.8317   -1.3843 C   0  0  2  0  0  0  0  0  0  0  0  0\n    4.5019   -3.3183   -0.0435 C   0  0  1  0  0  0  0  0  0  0  0  0\n    3.5467   -3.5748    1.1232 C   0  0  1  0  0  0  0  0  0  0  0  0\n    4.1733   -3.2899    2.4960 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9966   -5.7434   -2.4582 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0165   -3.7773   -2.3494 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6552   -1.9193   -0.1817 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1297   -4.9389    1.1294 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9142   -1.9631    2.8937 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6070   -6.9320    0.2034 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.4352   -5.9046   -1.0043 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1286   -3.2497   -1.7329 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.4178   -3.8630    0.1725 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6540   -2.9530    1.0304 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.2483   -3.4688    2.4675 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7307   -3.9551    3.2374 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6285   -5.5433   -3.1567 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.3322   -2.8665   -2.4107 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.3690   -1.3727    2.2840 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6886   -5.5112   -0.0972 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -1.9222   -6.0856   -1.2834 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -0.5071   -5.5014   -1.4207 C   0  0  2  0  0  0  0  0  0  0  0  0\n    0.2120   -5.5934   -0.0765 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -0.5895   -5.0119    1.0884 C   0  0  1  0  0  0  0  0  0  0  0  0\n    0.0178   -5.3292    2.4631 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.6142   -5.9152   -2.5059 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1846   -6.2703   -2.3870 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4020   -4.8403   -0.2073 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9156   -5.5381    1.0869 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8745   -4.2889    2.8770 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.7958   -6.3657    0.1484 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8489   -7.1500   -1.0550 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5838   -4.4684   -1.7645 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3570   -6.6504    0.1345 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6610   -3.9267    1.0005 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5627   -6.2723    2.4319 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7828   -5.4159    3.1991 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0668   -6.2882   -3.2039 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0906   -5.9452   -2.4453 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6230   -4.2591    2.2725 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.9816   -1.3349   -0.1372 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -5.9462   -2.2860   -1.3280 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -4.6083   -3.0299   -1.4640 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -4.2400   -3.6569   -0.1213 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -4.2901   -2.6743    1.0499 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -4.1663   -3.3560    2.4198 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.2354   -1.6320   -2.5497 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7738   -4.0426   -2.4370 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9087   -4.1178   -0.2482 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.5271   -1.9627    1.0455 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8211   -3.3762    2.8433 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -7.3407   -1.0028    0.1029 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.7342   -3.0071   -1.1095 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.8452   -2.3240   -1.7976 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.9774   -4.4298    0.0818 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4841   -1.9432    0.9691 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5610   -4.3711    2.3788 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7399   -2.7908    3.1561 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.1805   -2.2877   -3.2510 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.9578   -4.5558   -2.4899 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.3251   -3.9378    2.2378 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7679    3.8433   -0.1485 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -5.4791    3.2263   -1.3470 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -5.2263    1.7167   -1.4858 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -5.4966    1.0327   -0.1479 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -4.7688    1.6812    1.0309 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -5.2326    1.1538    2.3963 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.1374    3.8626   -2.5641 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.1174    1.2199   -2.4651 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.0270   -0.2947   -0.2760 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.9850    3.0915    1.0297 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.4244    0.0775    2.8136 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.3586    5.1116    0.0927 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.5366    3.3938   -1.1380 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.1975    1.5611   -1.8153 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.5621    1.1286    0.0472 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.6926    1.5091    0.9572 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.2759    0.8418    2.3495 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.1418    1.9479    3.1385 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -5.6043    3.4064   -3.2709 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -6.0086    0.2625   -2.5226 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.5614   -0.6540    2.2022 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0413    6.1038   -0.1180 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -0.8834    6.3045   -1.3173 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -1.9024    5.1642   -1.4682 C   0  0  2  0  0  0  0  0  0  0  0  0\n   -2.6165    4.9431   -0.1373 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -1.6646    4.7757    1.0482 C   0  0  1  0  0  0  0  0  0  0  0  0\n   -2.3783    4.8084    2.4078 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1646    6.4421   -2.5290 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8384    5.5524   -2.4554 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3591    3.7475   -0.2755 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6979    5.8237    1.0599 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7134    3.5039    2.8235 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7357    7.3201    0.1216 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4153    7.2336   -1.1067 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3779    4.2644   -1.7961 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.2078    5.8345    0.0565 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1287    3.8271    0.9781 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.2765    5.4233    2.3528 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7091    5.2378    3.1539 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8091    6.5250   -3.2386 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5143    4.8658   -2.5248 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.3647    3.1515    2.2079 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.8135    3.7985   -0.0619 C   0  0  2  0  0  0  0  0  0  0  0  0\n    4.3879    4.6328   -1.2637 C   0  0  1  0  0  0  0  0  0  0  0  0\n    2.8622    4.7231   -1.4219 C   0  0  2  0  0  0  0  0  0  0  0  0\n    2.2126    5.1187   -0.0956 C   0  0  1  0  0  0  0  0  0  0  0  0\n    2.6930    4.2811    1.0931 C   0  0  1  0  0  0  0  0  0  0  0  0\n    2.2727    4.8573    2.4527 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9478    4.1628   -2.4761 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5881    5.7035   -2.4044 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8215    4.9238   -0.2435 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1155    4.1750    1.1088 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0252    4.3334    2.8487 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.1788    4.0895    0.1970 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7828    5.6262   -1.0446 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4854    3.7543   -1.7586 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2833    3.2711    1.0178 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2221    5.9450    2.4050 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0111    4.5758    3.2047 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6130    4.7207   -3.1843 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6307    5.8063   -2.4788 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3615    4.6392    2.2212 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.9654   -1.3939   -0.0334 C   0  0  2  0  0  0  0  0  0  0  0  0\n    6.3606   -0.5330   -1.2274 C   0  0  1  0  0  0  0  0  0  0  0  0\n    5.4793    0.7165   -1.3851 C   0  0  2  0  0  0  0  0  0  0  0  0\n    5.4039    1.4542   -0.0509 C   0  0  1  0  0  0  0  0  0  0  0  0\n    5.0108    0.5581    1.1251 C   0  0  1  0  0  0  0  0  0  0  0  0\n    5.1846    1.2381    2.4905 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3530   -1.2571   -2.4432 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0789    1.5475   -2.3604 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.4069    2.4453   -0.1971 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.8144   -0.6198    1.1411 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9899    1.8755    2.8829 O   0  0  0  0  0  0  0  0  0  0  0  0\n    7.0428   -2.2799    0.2277 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3810   -0.2224   -0.9987 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.4892    0.4076   -1.7274 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4015    1.8318    0.1594 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9668    0.2500    1.0389 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.9994    1.9616    2.4544 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.4239    0.4838    3.2406 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5869   -0.6430   -3.1467 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5689    2.3652   -2.4229 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8152    2.5892    2.2600 H   0  0  0  0  0  0  0  0  0  0  0  0\n  2  3  1  0\n  2 11  1  0\n  2 31  1  0\n  2 13  1  6\n  3  4  1  0\n  3  8  1  0\n  3 14  1  1\n  4  5  1  0\n  4  9  1  0\n  4 15  1  6\n  5  6  1  0\n  5 10  1  0\n  5 16  1  1\n  6  7  1  0\n  6 11  1  0\n  6 17  1  6\n  7 12  1  0\n  7 18  1  0\n  7 19  1  0\n  8 20  1  0\n  9 21  1  0\n 10127  1  0\n 12 22  1  0\n 23 24  1  0\n 23 32  1  0\n 23 52  1  0\n 23 34  1  6\n 24 25  1  0\n 24 29  1  0\n 24 35  1  1\n 25 26  1  0\n 25 30  1  0\n 25 36  1  6\n 26 27  1  0\n 26 31  1  0\n 26 37  1  1\n 27 28  1  0\n 27 32  1  0\n 27 38  1  6\n 28 33  1  0\n 28 39  1  0\n 28 40  1  0\n 29 41  1  0\n 30 42  1  0\n 33 43  1  0\n 44 45  1  0\n 44 53  1  0\n 44 73  1  0\n 44 55  1  6\n 45 46  1  0\n 45 50  1  0\n 45 56  1  6\n 46 47  1  0\n 46 51  1  0\n 46 57  1  1\n 47 48  1  0\n 47 52  1  0\n 47 58  1  6\n 48 49  1  0\n 48 53  1  0\n 48 59  1  1\n 49 54  1  0\n 49 60  1  0\n 49 61  1  0\n 50 62  1  0\n 51 63  1  0\n 54 64  1  0\n 65 66  1  0\n 65 74  1  0\n 65 94  1  0\n 65 76  1  6\n 66 67  1  0\n 66 71  1  0\n 66 77  1  6\n 67 68  1  0\n 67 72  1  0\n 67 78  1  1\n 68 69  1  0\n 68 73  1  0\n 68 79  1  6\n 69 70  1  0\n 69 74  1  0\n 69 80  1  1\n 70 75  1  0\n 70 81  1  0\n 70 82  1  0\n 71 83  1  0\n 72 84  1  0\n 75 85  1  0\n 86 87  1  0\n 86 95  1  0\n 86115  1  0\n 86 97  1  1\n 87 88  1  0\n 87 92  1  0\n 87 98  1  6\n 88 89  1  0\n 88 93  1  0\n 88 99  1  1\n 89 90  1  0\n 89 94  1  0\n 89100  1  6\n 90 91  1  0\n 90 95  1  0\n 90101  1  1\n 91 96  1  0\n 91102  1  0\n 91103  1  0\n 92104  1  0\n 93105  1  0\n 96106  1  0\n107108  1  0\n107116  1  0\n107135  1  0\n107118  1  1\n108109  1  0\n108113  1  0\n108119  1  1\n109110  1  0\n109114  1  0\n109120  1  6\n110  1  1  6\n110111  1  0\n110115  1  0\n111112  1  0\n111116  1  0\n111121  1  6\n112117  1  0\n112122  1  0\n112123  1  0\n113124  1  0\n114125  1  0\n117126  1  0\n127128  1  0\n127136  1  0\n127138  1  1\n128129  1  0\n128133  1  0\n128139  1  1\n129130  1  0\n129134  1  0\n129140  1  6\n130131  1  0\n130135  1  0\n130141  1  1\n131132  1  0\n131136  1  0\n131142  1  6\n132137  1  0\n132143  1  0\n132144  1  0\n133145  1  0\n134146  1  0\n137147  1  0\nM  END\n\n'

def remove_CONNECT_LINES(fname):
	"""
	PRE: fname contains a valid PDB file of a molecule with 3D coordinates, MOLECULES NEED TO BE NAMED
	POST: The lines containing CONNECT and MASTER are removed from the original file, the original file IS MODIFIED
	"""
	with open(fname, 'rb') as r:
		lines = r.readlines()
	with open(fname, 'wb') as w:
		w.writelines([x for x in lines if 'CONECT' not in x if 'MASTER' not in x][1:])

def fix_PDB_spacing(fname, hostname = 'CB7', strandname = 'A'):
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
				#print new_lines[-1]

			else:
				new_lines.append(line.strip())
				#print new_lines[-1]
	with open(fname, 'wb') as w:
		w.writelines('\n'.join(new_lines)+'\n')

def make_pdb_complex_with_named_residues(RDKIT_BLOCK_GUEST, pdb_file_guest, pdb_file_cb, pdb_file_complex, isCB=True):
	"""
	PRE: Takes a RDKIT_BLOCK with valid 3D coordinates (basically a SDF file in text format), centered around 0 and with principal axis (in the PCA sense) aligned along the z axis
	POST: Creates three PDB files with properly named residues (default GST for the guest, CB7 for CB7) and connection records to be used in AMBER
	"""
	flavour = 28
	guest = Chem.MolFromMolBlock(RDKIT_BLOCK_GUEST, removeHs=False)
	converge_molecule(guest)
	atm_dic = {}
	for atom in guest.GetAtoms():
		if atom.GetSymbol() not in atm_dic:
			atm_dic[atom.GetSymbol()] = 1
		else: atm_dic[atom.GetSymbol()] += 1

		atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} GST'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10))))
		atom.GetMonomerInfo().SetResidueNumber(1)
	Chem.MolToPDBFile(guest, pdb_file_guest, flavor=flavour)
	fix_PDB_spacing(pdb_file_guest)
	print 'GUEST IS DONE'

	if isCB==True:
		CB = Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False)
	else:
		CB = Chem.MolFromMolBlock(get_CD_BLOCK(), removeHs=False)
	converge_molecule(CB)
	atm_dic = {}
	for atom in CB.GetAtoms():
		if atom.GetSymbol() not in atm_dic:
			atm_dic[atom.GetSymbol()] = 1
		else: atm_dic[atom.GetSymbol()] += 1
		atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} CB7'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10))))
		atom.GetMonomerInfo().SetResidueNumber(2)

	Chem.MolToPDBFile(CB, pdb_file_cb, flavor=flavour)
	fix_PDB_spacing(pdb_file_cb)
	print 'CB is DONE'

	complex_CB_host = Chem.CombineMols(guest, CB)
	converge_molecule(complex_CB_host)
	# yo = Chem.GetMolFrags(complex_CB_host)
	# for els in yo:
	# 	print els
	complex_CB_host.SetProp('_Name', 'CMP')
	Chem.MolToPDBFile(complex_CB_host, pdb_file_complex, flavor=flavour)
	remove_CONNECT_LINES(pdb_file_complex)
	fix_PDB_spacing(pdb_file_complex)
	print 'COMP is DONE'

def get_binding_energy_withCB7(sdfile):
	"""
	Takes a valid SDF file in input and returns its binding enregy with CB7
	"""
	def converge(mol, n_steps = 100000, tol=1e-9):
		"""Converges all of the molecule built-in conformations using n_threads"""
		# converged = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=100000, ignoreInterfragInteractions=False, nonBondedThresh=100.0)
		converged = []
		for ids in range(mol.GetNumConformers()):
			AllChem.MMFFSanitizeMolecule(mol)  
			ff = AllChem.MMFFGetMoleculeForceField(mol, confId=ids, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94', mmffVerbosity = 1), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
			ff.Initialize()
			cf = ff.Minimize(n_steps, tol, tol)
			converged.append((cf, ff.CalcEnergy()))
		return sum([x[0] for x in converged]), ff.CalcEnergy()

	mol = Chem.MolFromMolFile(sdfile, removeHs=False)
	mol = Chem.MolFromMolBlock(align_mol(Chem.MolToMolBlock(mol)[0]), removeHs = False)
	converged_mol, mol_E = converge(mol)
	complex_with_CB = Chem.CombineMols(mol, Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False))
	Chem.GetSSSR(complex_with_CB)
	converged_comp, complex_E = converge(complex_with_CB)
	Chem.MolToMolFile(mol, sdfile+'_MOL.sdf')
	Chem.MolToMolFile(complex_with_CB, sdfile+'_COMPLEX.sdf')
	print 'MMFF94 BINDING ENERGY is {0:.4f} (GUEST_ENERGY:{3}, converged=={1}; COMPLEX_ENERGY:{4}, COMPLEX_CONVERGE=={2})'.format(complex_E-mol_E-(-1451.415064), converged_mol==0, converged_comp==0, mol_E, complex_E)

def assign_features(rdkitmol):
	"""
	PRE: Takes in a valid rdkitmol with 3d coordinates
	POST: Returns the following properties as mol properties:
		  - AROMATIC CYCLES...
	"""
	pass

def get_energy_contributions(rdkitmol, label='Guest:', fout='./sample.sdf'):
	"""
	PRE: Takes in a rdkit mol
	POST: Assigns to it energy properties of MMFF94
	"""
	mp = AllChem.MMFFGetMoleculeProperties(rdkitmol)
	for i in range(7):
		termList = [['BondStretch', False], ['AngleBend', False],
		['StretchBend', False], ['OopBend', False], ['Torsion', False],
		['VdW', False], ['Electrostatic', False]]
		termList[i][1] = True
		mp.SetMMFFBondTerm(termList[0][1])
		mp.SetMMFFAngleTerm(termList[1][1])
		mp.SetMMFFStretchBendTerm(termList[2][1])
		mp.SetMMFFOopTerm(termList[3][1])
		mp.SetMMFFTorsionTerm(termList[4][1])
		mp.SetMMFFVdWTerm(termList[5][1])
		mp.SetMMFFEleTerm(termList[6][1])
		ff = AllChem.MMFFGetMoleculeForceField(rdkitmol, mp)
		rdkitmol.SetProp(label+termList[i][0], '{0:12.4f}'.format(ff.CalcEnergy()))
		# print '{0:>16s} energy: {1:12.4f} kcal/mol'.format(termList[i][0],ff.CalcEnergy())
	ff = AllChem.MMFFGetMoleculeForceField(rdkitmol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(rdkitmol, mmffVariant='MMFF94', mmffVerbosity = 1), ignoreInterfragInteractions=False, nonBondedThresh=100.0) 
	rdkitmol.SetProp(label+'Total', '{0:12.4f}'.format(ff.CalcEnergy()))
	w=Chem.SDWriter(fout)
	w.write(rdkitmol)
	w.close()
	return rdkitmol


def convert_sdf_to_pdb(fin, fout='', converge=False, residueName=('GST', 'CB7'), residueNumber=(1,2)):
	"""
	PRE : Takes in the absolute path to a SDF file, if a single molecule is converted only the first element of both residueName and residueNumber is used, both in the case of a complex
	POST: Produces a PDB file by converting the original molecule from the SDF file, DOES NOT optimize the conformation by default
	"""
	if fout == '':
		fout = fin[:-4]+'_SDF2PDB.pdb'
	mol = Chem.MolFromMolFile(fin, removeHs = False)
	flavour = 28
	if converge:
		converge_molecule(mol)
	
	isComplex = (len(Chem.GetMolFrags(mol)) != 1)  

	print '{} contains {} Fragments'.format(fin, len(Chem.GetMolFrags(mol)))

	if not isComplex:
		atm_dic = {}
		for atom in mol.GetAtoms():
			if atom.GetSymbol() not in atm_dic:
				atm_dic[atom.GetSymbol()] = 1
			else: atm_dic[atom.GetSymbol()] += 1

			atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} {}'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10), residueName[0])))
			atom.GetMonomerInfo().SetResidueNumber(residueNumber[0])
		Chem.MolToPDBFile(mol, fout, flavor=flavour)
		fix_PDB_spacing(fout, hostname=residueName[1])
		return
	else:
		comp = Chem.GetMolFrags(mol, asMols=True)
		if comp[0].GetNumAtoms() < comp[1].GetNumAtoms(): guest,host = comp[0], comp[1]
		else: guest,host = comp[1], comp[0]
		for res in zip((guest, host),[0,1]):
			atm_dic = {}
			for atom in res[0].GetAtoms():
				if atom.GetSymbol() not in atm_dic:
					atm_dic[atom.GetSymbol()] = 1
				else: atm_dic[atom.GetSymbol()] += 1

				atom.SetMonomerInfo(Chem.rdchem.AtomPDBResidueInfo('{} {}'.format(atom.GetSymbol()+str(atm_dic[atom.GetSymbol()])+' '*int(atm_dic[atom.GetSymbol()]<10), residueName[res[1]])))
				atom.GetMonomerInfo().SetResidueNumber(residueNumber[res[1]])

		comp_annotated = Chem.CombineMols(guest, host)
		comp_annotated.SetProp('_Name', mol.GetProp('_Name'))
		Chem.MolToPDBFile(comp_annotated, fout, flavor=flavour)
		remove_CONNECT_LINES(fout)
		fix_PDB_spacing(fout, hostname=residueName[1])

def create_pdb_and_topology_from_sdf(fin, resName = ('GST', 'HST'), resNumber = (1,2)):
	"""
	PRE  : Takes in a sdf file, needs to have a name, if the molecule is not a complex only the first element of both resName and resNumber is used to define the residue
	POST : Produces its amber topology files and pdb files. The files ".sh" produced by make_topology need to be run (they in turn call tleap) to create the actual topology files (using chmod +x $file) may help
		   NOTE: The sh file needs to be run with the amber programs in the path (tleap, antechamber and the like, or AmbertTools) 
		   NOTE: The suffixes '-GUEST.sdf', '-HOST.sdf', '-SDF2PDB.pdb', '-COMPLEX-SDF2PDB.pdb' and combinations thereof are all hardcoded. Don't get lost!
		   NOTE: Now the output files are just sdf[:-4]+'.pdb', much easier!
	"""
	print fin
	mol = Chem.MolFromMolFile(fin, removeHs = False)
	isComplex = (len(Chem.GetMolFrags(mol)) != 1)
	print 'isComplex? {}'.format(isComplex)
	if isComplex:
		# Gets the host and guest
		comp = Chem.GetMolFrags(mol, asMols=True)
		if comp[0].GetNumAtoms() < comp[1].GetNumAtoms(): guest,host = comp[0], comp[1]
		else: guest,host = comp[1], comp[0]
		# Writes them to two additional sdf files
		outpdbfiles = []
		for m in zip((guest,host),(fin[:-4]+'-guest.sdf', fin[:-4]+'-host.sdf'), resName, resNumber):
			m[0].SetProp('_Name', m[2])
			Chem.MolToMolFile(m[0], m[1])
			# converts them to proper pdb
			outpdbfiles.append(m[1][:-4]+'-SDF2PDB.pdb')
			convert_sdf_to_pdb(m[1], residueName=(m[2], resName[1]), residueNumber=(m[3],m[3]), fout=outpdbfiles[-1])
		
		outpdbfiles.append(fin[:-4]+'.pdb')
		convert_sdf_to_pdb(fin, residueName=resName, residueNumber=resNumber, fout=outpdbfiles[-1])

		make_topology(outpdbfiles[0], outpdbfiles[1], outpdbfiles[2])

	else:
		fout = fin[:-4]+'.pdb'
		convert_sdf_to_pdb(fin, residueName=resName, residueNumber=resNumber, fout=fout)
		make_topology(fout)

def make_topology(guestPDB, hostPDB = '', compPDB = '', tleapexecpath='/home/macenrola/Thesis/AMBER/amber16/bin/'):  
	"""
	PRE  : If a guest, just takes in a pdb file, if a complex, takes three corresponding to the guest, host and complex. The complex would have no connect records while the guest and host pdb would have some 
		   If this method is used as a chain DON't forget to run the .sh scripts of the host guest subfiles from the complex before running the code to get the final complex prmtop. Alternatively, run everything TWICE or more haha
	POST : Writes an amber script that generates prmtop, frcmod for the compPDB file if complex or the guestPDB file if single molecule
		   NOTE: Currently the script to generate the topologies is named guestPDB without '.pdb' 
	"""
	isComplex = (hostPDB !='')

	def make_topology_one_mol(tleapexecpath, targetPDB):
		"""
		PRE  : Takes in one molecule, that has a NAMED residue number, it is a valid file from which the charge can be extracted. USING GASTEIGER CHARGES BECAUSE BCC TOO LONG
		POST : writes its tleap script for topology generation
			   NOTE: The output file is named as the pdbfile without '.pdb' + -run_tops.sh 
		"""
		script_bin = [
			'# Default atom type is GAFF -at gaff',
			'{0}antechamber -i {2}.pdb -fi pdb -o {2}.mol2 -fo mol2 -c gas -nc {1} -s 2 -j 4',
			'# Unit names in the text are not variables but actual residue/unit names ><',
			'{0}parmchk2 -i {2}.mol2 -f mol2 -o {2}.frcmod -f frcmod',
			'{0}tleap -s -f {2}-tleap.in > {2}-tleap.log',
			]

		script_tleap = [
			'source leaprc.gaff',
			'loadAmberParams {2}.frcmod',
			'{3} = loadMol2 {2}.mol2',
			'check {3}',
			'loadamberparams gaff.dat',
			'saveamberparm {3} {2}.prmtop {2}.rst7',
			'saveoff {3} {2}.lib',
			'quit'
			]
		m = Chem.MolFromPDBFile(targetPDB, removeHs=False)
		target = targetPDB[:-4]
		charge = Chem.GetFormalCharge(m)
		residuename = Chem.SplitMolByPDBResidues(m).keys()[0]

		script_bin = [x.format(tleapexecpath, charge, target, residuename)+'\n' for x in script_bin]
		script_tleap = [x.format(tleapexecpath, charge, target, residuename)+'\n' for x in script_tleap]

		fout1 = target+'-run_tops.sh'
		fout2 = target+'-tleap.in'

		with open(fout1, 'wb') as w:
			w.writelines(script_bin)
		with open(fout2, 'wb') as w:
			w.writelines(script_tleap)


		call('chmod +x {}'.format(fout1), shell=True)
		call('{}'.format(fout1),shell=True)


	def make_topology_mol_combination(guestPDB, hostPDB, compPDB, tleapexecpath):
		"""
		PRE: Takes in the pdb files of the two building blocks AND the pdb file without connect of the complex of the combinations and 
		     supposes that all associated files ie target.frcmod, target.prmtop, target.mol2, target.lib exist for target.pdb
		POST : Generates the rst7 and prmtop files for the complex, they will use the path and prefix/radical of the pdb files
		"""
		script_tleap = [
		'source leaprc.gaff',
		'loadoff {0}.lib',
		'loadoff {1}.lib',
		'loadamberparams {0}.frcmod',
		'loadamberparams {1}.frcmod',
		'{2} = loadmol2 {0}.mol2',
		'{3} = loadmol2 {1}.mol2',
		'COMPLEX = loadPDB {4}.pdb',
		'savemol2 COMPLEX {4}.mol2 1',
		'saveamberparm COMPLEX {4}.prmtop {4}.rst7',
		'quit'
		]
		script_bin = ['{5}tleap -s -f {4}-tleap.in > {4}-tleap.log']

		guest_name = guestPDB[:-4]
		host_name = hostPDB[:-4]
		comp_name = compPDB[:-4]

		guest = Chem.MolFromPDBFile(guestPDB, removeHs=False)
		guest_res = Chem.SplitMolByPDBResidues(guest).keys()[0]
		host = Chem.MolFromPDBFile(hostPDB, removeHs=False)
		host_res = Chem.SplitMolByPDBResidues(host).keys()[0]
		print guest_res, host_res

		script_tleap = [x.format(guest_name, host_name, guest_res, host_res, comp_name, tleapexecpath)+'\n' for x in script_tleap]
		script_bin = [x.format(guest_name, host_name, guest_res, host_res, comp_name, tleapexecpath)+'\n' for x in script_bin]


		fout1 = comp_name+'-run_tops.sh'
		fout2 = comp_name+'-tleap.in'

		with open(fout1, 'wb') as w:
			w.writelines(script_bin)
		with open(fout2, 'wb') as w:
			w.writelines(script_tleap)

		call('chmod +x {}'.format(fout1), shell=True)
		call('{}'.format(fout1),shell=True)

	if not isComplex:
		make_topology_one_mol(tleapexecpath, guestPDB)
	else:
		make_topology_one_mol(tleapexecpath, guestPDB)
		make_topology_one_mol(tleapexecpath, hostPDB)
		make_topology_mol_combination(guestPDB, hostPDB, compPDB, tleapexecpath)

	# if tleapexecpath == '':
	# 	raise Exception('tleap path not set')
	return

def get_multiple_conformations(SMILES='c1cccc(CCCCC(C)C)c1', N=200):
	"""
	PRE: Takes in a SMILES
	POST: Finds its N lowest energy conformers and minimizes them using MMFF94
	"""
	mol = Chem.MolFromSmiles(SMILES)
	mol = Chem.AddHs(mol)
	ids = AllChem.EmbedMultipleConfs(mol, N, AllChem.ETKDG())
	min_e = (-1, 1e9)
	for i in ids:
		ff = AllChem.MMFFGetMoleculeForceField(mol, confId=i, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol))
		ff.Minimize()
		e = ff.CalcEnergy()
		print i, e
		if e<min_e[1]: min_e=(i,e)
		Chem.MolToMolFile(mol, confId=i, filename='/home/macenrola/Desktop/illustration/{}.sdf'.format(i))
	print 'best energy is {}'.format(min_e)

def test_make_PDB_files(SMILES = 'C[N+](C12CC3C4C(C2)C2C(C1)C3CC(C4)(C2)[N+](C)(C)C)(C)C', path = '/home/macenrola/Thesis/AMBER/converge_with_amber/'):
	"""
	PRE: Takes a SMILES in input 
	POST: Writes 3 PDB files to path, one for the guest, one for CB and one for the complex. The PDB files are formatted to be read by AMBER.
	"""
	from rdkit.Chem import AllChem
	mol = Chem.MolFromSmiles(SMILES)
	mol = Chem.AddHs(mol)
	AllChem.EmbedMolecule(mol, AllChem.ETKDG())
	print Chem.MolToMolBlock(mol)
	out_block =  align_mol(Chem.MolToMolBlock(mol))[0]
	guest_name = 'DADA' # Diadamantyl diammonium
	fguest = path+guest_name+'.pdb'
	fCB = path+'CB.pdb'
	fcomplex = path+guest_name+'-CB.pdb'
	make_pdb_complex_with_named_residues(out_block, fguest, fCB, fcomplex)

def test_align_molecule(SMILES='C[N+](C12CC3C4C(C2)C2C(C1)C3CC(C4)(C2)[N+](C)(C)C)(C)C', inSDF='', fname=''):
	"""
	PRE: Takes a SMILES in input 
	POST: Returns an RDKIT_BLOCK with 3D coordinates 0 centered and with principal axis aligned with z
	"""
	from rdkit.Chem import AllChem
	if inSDF == '':
		mol = Chem.MolFromSmiles(SMILES)
		mol = Chem.AddHs(mol)
		AllChem.EmbedMolecule(mol, AllChem.ETKDG())
		AllChem.MMFFOptimizeMolecule(mol)

	else:
		mol = Chem.MolFromMolFile(inSDF, removeHs=False)

	# get_energy_contributions(mol)
	print Chem.MolToMolBlock(mol)
	out_block =  align_mol(Chem.MolToMolBlock(mol))[0]
	if fname !='':
		Chem.MolToMolFile(Chem.MolFromMolBlock(out_block, removeHs=False), fname)
		with open(fname, 'rb') as r:
			print r.readlines()
	return out_block

def minimize_molecule(fname):
	"""
	PRE : Takes in a SDF file
	POST: Prints the energy of the molecule according to MMFF94
	"""
	mol = Chem.MolFromMolFile(fname, removeHs=False)
	ff = AllChem.MMFFGetMoleculeForceField(mol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94', mmffVerbosity = 1), ignoreInterfragInteractions=False, nonBondedThresh=100.0) 	
	converged = ff.Minimize(1000000)
	print 'Energy is {}, converged is {}'.format(ff.CalcEnergy(), converged)

def converge_molecule(molecule, n_steps=100000, tol=1e-8):
	"""
	PRE: Takes in a molecule
	POST: Returns its energy and converges the molecule in situ
	"""
	AllChem.MMFFSanitizeMolecule(molecule)  
	ff = AllChem.MMFFGetMoleculeForceField(molecule, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant='MMFF94'), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
	ff.Initialize()
	cf = ff.Minimize(n_steps, tol, tol)
	return cf, ff.CalcEnergy()	


def get_frequency_report(target_pdb, target_prmtop, nabpath='/home/macenrola/Thesis/AMBER/amber16/bin/nab'):
	"""
	PRE : Takes in a target pdb file and its prmtop file
		  NOTE: Should try to include the frcmod file as well.
	POST : create a nab short code that returns a callable c program (callable as ./a.out by default)the normal mode analysis for that particular molecule within the AMBER GAFF force field
		   To tune the name of the a.out file 
		   Returns the total entropy including the translational and rotational + obviously vibrational
	"""
	script_nab = [
	'molecule m;',
	'float x[4000], fret;',
	'// m = getpdb_prm( "{}.pdb", "leaprc.gaff", "", 0);',

	'm = getpdb( "{0}" );',
	'readparm( m, "{1}");',
	'mm_options( "cut=999., ntpr=50, nsnb=99999, diel=C, gb=0, dielc=1.0" );',
	'mme_init( m, NULL, "::Z", x, NULL);',
	'setxyz_from_mol( m, NULL, x );',

	'// conjugate gradient minimization',
	# 'dgrad = 3*natoms*0.001;',
	'conjgrad(x, 3*m.natoms, fret, mme, 0.001, 0.001, 200000 );',

	'// Newton-Raphson minimization\fP',
	'mm_options( "ntpr=1" );',
	'newton( x, 3*m.natoms, fret, mme, mme2, 0.00000001, 0.0, 60 );',

	'// get the normal modes:',
	'nmode( x, 3*m.natoms, mme2, 0, 0, 0.0, 0.0, 0);'
	]

	script_nab = [x.format(target_pdb, target_prmtop)+'\n' for x in script_nab]

	fout = target_pdb[:-4]+'-freq.nab'
	frequency_report = target_pdb[:-4]+'frequency_report'
	with open(fout, 'wb') as w:
		w.writelines(script_nab)
	cmd = '{} {} -o {}'.format(nabpath, fout, fout[:-4]+'.out')
	call(cmd, shell=True)

	cmd2 = 'cd {} && ./{} > {}'.format('/'.join(fout.split('/')[:-1]), fout.split('/')[-1][:-4]+'.out', frequency_report)
	print cmd2
	call(cmd2, shell=True)

	# return return_entropy_from_frequency_report(frequency_report)
	with open(frequency_report, 'rb') as r:
		for line in r:
			if 'Total:' in line:
				Entropy = line.strip().split()[-1]
				break
			if 'Energy       =' in line:
				Energy = line.strip().split()[-1]
	print Entropy, Energy
	return Entropy, Energy
	
def return_entropy_from_frequency_report(frequency_report, wheretostartifcomplex=12):
	"""
	PRE  : Takes an AMBER frequency report as generated by get_frequency_report
		   wheretostart indicates how many degrees of freedom are lost upon binding, 6 means all degrees of freedom of the guest are lost, 12 means nothing is lost, 8 means only 2 arent lost like for the symmetry of CB
	POST : returns the entropy associated to the given frequencies 
		 : NOTE it is assumed that 2 degrees of freedom are not lost upon binding, rotation according to the axis that goes along CB and translation along that same axis
		   NOTE returns the entrop in cal/mol/K like in AMBER
	"""
	frequencies = []
	import math
	with open(frequency_report, 'rb') as r:
		pastff=False
		for i, line in enumerate(r):
			if i==0: 
				if 'COMPLEX' in line:
					isComplex=True
				else: isComplex=False
			if pastff:
				frequencies.append(line.split()[1])
			if 'ff   energy:' in line:
				pastff = True

	if isComplex:
		whereToStart=wheretostartifcomplex
	else:
		whereToStart=6
	frequencies = [float(x)*1e2*2.99e8 for x in frequencies[whereToStart:]]

	if frequencies[whereToStart]<=0:
		return 'IMPROPER_FREQS'

	k = 1.38E-23
	R = 8.3145
	T = 300
	beta = 1/(k*T)
	h = 6.626E-34
	Srrho = 0
	for f in frequencies:
	    inc = R*(beta*h*f*(math.exp(beta*h*f)-1)**-1 - math.log(1-math.exp(-beta*h*f)))
	    Srrho = inc+Srrho

	print 'Srroh is {} cal/mol.K or {} kcal/mol@300K'.format(Srrho/4.184, Srrho/4.184*300/1000)
	return Srrho/4.184


def make_job_kinetic_barrier():
	"""
	PRE  : Takes an amber line and 
	POST : create a job to run it several times with different initial velocities
	"""
	amberline = "sander -O -i md2.in -p CB7_DBOA_1H_in0_PD-NH3_wB97XD_6-31Gs_SDF2PDB.prmtop -c md_heat.rst7 -o md_eq_{0}.out -r md_eq_{0}.rst7 -x md_eq_{0}.nc"
	script_name = "kinetic_amber"
	canon_string = [
	"#!/bin/bash -l",
	"#$ -S /bin/bash",
	"#$ -cwd",
	"#$ -l h_rt=1:00:00",
	"#$ -l mem=2G",
	"#$ -l tmpfs=5G",
	"#$ -N Kinetic_barrier",
	"#$ -t 1-50",
	# "#$ -pe mpi 48",
	"module load python/2.7.9",
	"module unload amber",
	"module load amber/16/mpi/intel-2015-update2",
	amberline.format("$SGE_TASK_ID")]
	with open("/home/macenrola/Desktop/{}".format(script_name+'.sh'), 'wb') as w:
		w.writelines('\n'.join(canon_string))


def get_transition(onedatfile):
	"""
	PRE  : Provide data rmsd for three trajectories inside CB7  
	POST : Return a file for the transition to 2A rmsd
	"""
	prefix =  '/'.join(onedatfile.split('/')[:-1]) +'/'
	flist = glob.glob(prefix+'*.dat')
	outfile = prefix+'steps_To_2A_RMSD'
	with open(outfile, 'wb'): pass
	for i, f in enumerate(flist):
		with open(f, 'rb') as r:
			for line in r:
				if 'F' in line:
					continue
				parts = line.split()
				if float(parts[1]) > 3:
					print i, f, parts[0]
					with open(outfile, 'ab') as a:
						a.write(line)
					break



if __name__ == "__main__":
	import glob
	from subprocess import call
	# make_job_kinetic_barrier()
	get_transition('/home/macenrola/Desktop/MD_trajectories/DFT/C2H4/rmstraj38.dat')
	# create_pdb_and_topology_from_sdf('/home/macenrola/Thesis/hydrocarbons/CB_REFERENCE_VALUES/CB_candidate.sdf', ('CB7', 'haha'))
	# get_frequency_report('/home/macenrola/Thesis/hydrocarbons/CB_REFERENCE_VALUES/CB_candidate.pdb', '/home/macenrola/Thesis/hydrocarbons/CB_REFERENCE_VALUES/CB_candidate.prmtop')
	# print return_entropy_from_frequency_report('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/test/xaf_OUT_COMPLEX0frequency_report')
	# get_frequency_report('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/test/CB_candidate.pdb', '/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/test/CB_candidate.prmtop')
	# create_pdb_and_topology_from_sdf('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/test/CB_candidate.sdf')
	# get_frequency_report('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/butina/xaf_OUT.sdf_SAMPLE-SDF2PDB.pdb', '/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/guest_of_interest/butina/xaf_OUT.sdf_SAMPLE-SDF2PDB.prmtop')
	# minimize_molecule('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/CYCLODEXTRIN_VIBRATIONS/beta-cd_aligned_e_contrib.sdf')
	# flist = glob.glob('/home/macenrola/Thesis/VIBRATIONS/test_vibrations/FOUR_BENZENE/*.sdf')
	# print flist
	# # ['/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/241-best-complex.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/5054-best-complex.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/4409-best-conf.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/156391-best-complex.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/4409-best-complex.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/241-best-conf.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/5054-best-conf.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/3394-best-complex.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/3394-best-conf.sdf', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/200confs/BEST_CONF_STUDY_CD/156391-best-conf.sdf']

	# for fname in flist:
	# 	create_pdb_and_topology_from_sdf(fname)
	# 	if 'SDF2PDB' in fname:
	# 		continue
	# 	if 'complex' in fname:
	# 		get_frequency_report(fname[:-4]+'-COMPLEX-SDF2PDB.pdb', fname[:-4] + '-COMPLEX-SDF2PDB.prmtop')
	# 	else:
	# 		get_frequency_report(fname[:-4]+'-SDF2PDB.pdb', fname[:-4] + '-SDF2PDB.prmtop')

	# create_pdb_and_topology_from_sdf('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/CB7_VIBRATIONS/CB_candidate.sdf')
	# get_frequency_report('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/CB7_VIBRATIONS/CB_candidate-SDF2PDB.pdb', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/CB7_VIBRATIONS/CB_candidate-SDF2PDB.prmtop')
	# mol = Chem.MolFromMolFile('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd_aligned.sdf', removeHs=False)
	# converge_molecule(mol)
	# get_energy_contributions(mol, fout='/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd_aligned_e_contrib.sdf')
	# aligned_mol = Chem.MolFromMolBlock(align_mol(Chem.MolToMolBlock(mol)), removeHs=False)
	# Chem.MolToMolFile(aligned_mol, '/home/macenrola/Desktop/illustration/91_aligned.sdf')
	# get_binding_energy_withCB7('/home/macenrola/Desktop/illustration/91_aligned.sdf')
	# test_align_molecule('', '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd.sdf','/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/BETACYCLODEXTRINE')
	# get_multiple_conformations()
	# minimize_molecule('/home/macenrola/Desktop/sample.sdf')
	# Chem.MolToMolFile(Chem.MolFromMolBlock(align_mol(Chem.MolToMolBlock(Chem.MolFromMolFile('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd.sdf', removeHs=False))),removeHs=False), '/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd_aligned.sdf')
	# print align_mol(Chem.MolToMolBlock(Chem.MolFromMolFile('/home/macenrola/Dropbox/Vibrations/CYCLODEXTRINE/beta-cd.sdf', removeHs=False))).__repr__()

	pass
