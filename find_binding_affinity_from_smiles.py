# ==============================
# Scoring function

# pre: Will take a smiles in argument via sys.argv[1]
# post: will return an estimate of the binding affinity of the molecule with CB7

# ==============================

def align_mol(mol):
	"""
	PRE: Takes as input a RDKIT_mol with valid 3D coordinates (basically a SDF file in text format)
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

		transform_coord_centered = transform_coord.copy()
		transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
		transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
		transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])
		return transform_coord_centered, xthick, ythick, zthick

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
		for i in range(len(list_atoms)):
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


	# print RDKIT_BLOCK
	# mol=Chem.MolFromSmiles(s)
	# mol=Chem.AddHs(mol)
	# AllChem.EmbedMolecule(mol)
	RDKIT_BLOCK=Chem.MolToMolBlock(mol)
	atom_coords = get_atoms_coords(RDKIT_BLOCK)
	transformed_coords, xthick, ythick, zthick = align_xyz(atom_coords[0], atom_coords[1])
	aligned_mol=Chem.MolFromMolBlock(generate_mol_from_MDL(RDKIT_BLOCK, transformed_coords), removeHs=False)
	return aligned_mol, xthick, ythick, zthick

def align_mol_from_atm_list_and_xyz(xyz_file):
	"""
	PRE: Takes in input an atom list and a list of coordinates,
	POST: Returns an atom list and the transformed coordinates with the 3D coordinates rotated so that the main axis of the molecule (PCA wise) is aligned with the z axis and all coordinates are 0 centered
		This method exists to dock molecules that are not recognized as valid by the RDKIT
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

		transform_coord_centered = transform_coord.copy()
		transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
		transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
		transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])
		return transform_coord_centered, xthick, ythick, zthick

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
		for i in range(len(list_atoms)):
			radius = index_of_radii[list_atoms[i]]
			top_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(0)]
			bottom_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(np.pi)]
			sphere_points.append(top_sphere_point)
			sphere_points.append(bottom_sphere_point)


		new_sphere_points = [[coord_matrix[i][0]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.cos(azymuth_angle), coord_matrix[i][1]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.sin(azymuth_angle), coord_matrix[i][2]+index_of_radii[list_atoms[i]]*np.cos(inclination_angle)]
							for azymuth_angle in np.linspace(0, np.pi *2, angle_nbs*2) for inclination_angle in np.linspace(0, np.pi, angle_nbs + 1)[1:-1] for i in range(len(list_atoms))]
		sphere_points.extend(new_sphere_points)
		return sphere_points



	def get_atoms_coords(xyz_file):
		"""Takes as input an RDKIT BLOCK and returns a list of atoms with a numpy array containing the coordinates"""
		with open(xyz_file, 'rb') as r:
			block=[x.strip() for x in r.readlines()]
		atm_number = int(block[0])
		coords = [x.split() for x in block[2:]]
		atm_list = [x[0] for x in coords]
		coords_array = np.array([np.array(xi[1:]) for xi in coords], dtype=float)
		return atm_list, coords_array


	atom_coords = get_atoms_coords(xyz_file)
	transformed_coords, xthick, ythick, zthick = align_xyz(atom_coords[0], atom_coords[1])
	return atom_coords[0], transformed_coords, xthick, ythick, zthick

def get_CB_BLOCK():
	"""
	PRE: -
	POST: Returns the RDKIT_BLOCK of CB7 as a string
	"""
	return 'CB7\n     RDKit          3D\n\n126147  0  0  0  0  0  0  0  0999 V2000\n    1.1730   -4.4176    2.9511 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0750   -4.6349    3.6688 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0765   -3.4911    4.7539 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1723   -2.7865    4.5031 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9191   -3.3583    3.4714 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3124   -4.4064    2.9369 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0522   -3.3353    3.4421 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3130   -2.7704    4.4841 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0691   -3.0468    3.1437 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1939   -3.0104    3.0984 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1854   -0.4434    5.2692 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0591   -0.0113    5.8877 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0441    1.5519    5.6863 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2039    1.7887    4.9742 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9401    0.6201    4.7698 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2994   -0.4233    5.2450 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0305    0.6539    4.7401 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2819    1.8110    4.9656 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0894    0.5504    4.3212 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1741    0.6051    4.2737 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2527   -3.8136   -3.7077 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0124   -3.8818   -4.4671 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0312   -2.5540   -5.3160 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2779   -1.9153   -4.9189 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0091   -2.6783   -4.0054 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2331   -3.7867   -3.7189 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9639   -2.6375   -4.0281 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2071   -1.8908   -4.9333 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1545   -2.4399   -3.6082 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1102   -2.3785   -3.6467 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3096    4.4630   -2.8373 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0731    5.2303   -2.7776 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0598    5.7883   -1.3030 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2961    5.2698   -0.7356 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0425    4.5245   -1.6504 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1763    4.4888   -2.8686 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9303    4.5503   -1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1898    5.2750   -0.7590 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1845    4.0859   -1.4747 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0830    4.1303   -1.5454 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2378    3.8620    3.6411 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0016    4.6239    3.6943 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0176    5.4457    2.3491 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2631    5.0331    1.7185 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9935    4.1349    2.4992 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2476    3.8752    3.6162 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9787    4.1602    2.4607 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2224    5.0509    1.6961 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1395    3.7339    2.2694 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1253    3.7745    2.2088 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8611   -0.7617   -5.5808 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8794   -1.7513    5.3501 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7603   -0.7297   -5.6064 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7412   -1.7795    5.3807 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1982   -5.2981    0.6479 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0498   -5.9266    0.2399 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0334   -5.7852   -1.3302 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2138   -5.0813   -1.5932 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9522   -4.8379   -0.4333 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2882   -5.2657    0.6262 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0214   -4.8068   -0.4705 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2708   -5.0748   -1.6171 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1034   -4.3917   -0.3828 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1644   -4.3384   -0.4391 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7455   -5.3457    1.9916 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8706   -5.3200    1.9558 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8157   -4.8707   -2.9471 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7943   -4.9008   -2.9121 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7863    3.0997    4.7479 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8328    3.1315    4.7180 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8453    5.6668    0.5485 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7668    5.6860    0.5091 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1731    2.6985   -4.5639 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0766    2.5733   -5.2995 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0673    1.0736   -5.7843 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1835    0.5544   -5.2506 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9264    1.5226   -4.5717 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3114    2.6728   -4.5336 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0441    1.4838   -4.5307 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3002    0.5315   -5.2313 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0788    1.3907   -4.1449 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1843    1.3310   -4.0791 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7319    3.9561   -4.0995 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8841    3.9189   -4.0557 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0817   -5.6424    4.1014 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0872   -3.8701    5.7828 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0702   -0.3130    6.9419 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0393    2.1106    6.6299 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0054   -4.7893   -5.0823 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0345   -2.7314   -6.3981 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0901    6.0195   -3.5388 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0617    6.8837   -1.2538 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0074    5.2580    4.5891 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0242    6.5312    2.5046 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7799   -0.8998   -6.6666 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9121   -0.7356   -5.2918 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9301   -1.6557    5.0751 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7985   -2.0876    6.3912 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8160   -0.6826   -5.3367 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6624   -0.8715   -6.6906 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7974   -1.7001    5.1217 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6395   -2.1199    6.4192 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0657   -6.9691    0.5804 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0296   -6.7488   -1.8530 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8037   -5.0954    1.9086 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6356   -6.3657    2.3808 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9239   -5.0599    1.8472 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7805   -6.3449    2.3374 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8721   -4.6322   -2.8198 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7133   -5.8028   -3.5171 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8522   -4.6814   -2.7612 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6872   -5.8367   -3.4748 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8415    2.9390    4.5246 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6891    3.6916    5.6670 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8862    2.9907    4.4754 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7409    3.7289    5.6339 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9003    5.3920    0.5370 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7470    6.7545    0.6523 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6574    6.7728    0.6099 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8247    5.4251    0.4716 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0949    3.2963   -6.1235 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0742    0.9704   -6.8759 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7919    3.7784   -3.9142 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6170    4.7060   -4.8933 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8027    4.6722   -4.8487 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9349    3.7196   -3.8436 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1 65  1  0\n  2  3  1  0\n  2  6  1  0\n  2 85  1  0\n  3  4  1  0\n  3 86  1  0\n  4  5  1  0\n  4 54  1  0\n  5  1  1  0\n  6  7  1  0\n  6 66  1  0\n  7  8  1  0\n  8  3  1  0\n  9  5  2  0\n 10  7  2  0\n 11 12  1  0\n 12 13  1  0\n 12 16  1  0\n 12 87  1  0\n 13 14  1  0\n 13 88  1  0\n 14 15  1  0\n 14 69  1  0\n 15 11  1  0\n 16 17  1  0\n 17 18  1  0\n 18 13  1  0\n 18 70  1  0\n 19 15  2  0\n 20 17  2  0\n 21 22  1  0\n 21 68  1  0\n 22 23  1  0\n 22 26  1  0\n 22 89  1  0\n 23 24  1  0\n 23 90  1  0\n 24 25  1  0\n 25 21  1  0\n 26 27  1  0\n 26 67  1  0\n 27 28  1  0\n 28 23  1  0\n 28 53  1  0\n 29 25  2  0\n 30 27  2  0\n 31 32  1  0\n 31 84  1  0\n 32 33  1  0\n 32 36  1  0\n 32 91  1  0\n 33 34  1  0\n 33 92  1  0\n 34 35  1  0\n 34 71  1  0\n 35 31  1  0\n 36 37  1  0\n 36 83  1  0\n 37 38  1  0\n 38 33  1  0\n 38 72  1  0\n 39 35  2  0\n 40 37  2  0\n 41 42  1  0\n 41 69  1  0\n 42 43  1  0\n 42 46  1  0\n 42 93  1  0\n 43 44  1  0\n 43 94  1  0\n 44 45  1  0\n 45 41  1  0\n 46 47  1  0\n 46 70  1  0\n 47 48  1  0\n 48 43  1  0\n 48 72  1  0\n 49 45  2  0\n 50 47  2  0\n 51 24  1  0\n 51 80  1  0\n 51 95  1  0\n 51 96  1  0\n 52  8  1  0\n 52 16  1  0\n 52 97  1  0\n 52 98  1  0\n 53 99  1  0\n 53100  1  0\n 54 11  1  0\n 54101  1  0\n 54102  1  0\n 55 56  1  0\n 55 65  1  0\n 56 57  1  0\n 56 60  1  0\n 56103  1  0\n 57 58  1  0\n 57104  1  0\n 58 59  1  0\n 58 68  1  0\n 59 55  1  0\n 60 61  1  0\n 60 66  1  0\n 61 62  1  0\n 62 57  1  0\n 62 67  1  0\n 63 59  2  0\n 64 61  2  0\n 65105  1  0\n 65106  1  0\n 66107  1  0\n 66108  1  0\n 67109  1  0\n 67110  1  0\n 68111  1  0\n 68112  1  0\n 69113  1  0\n 69114  1  0\n 70115  1  0\n 70116  1  0\n 71 44  1  0\n 71117  1  0\n 71118  1  0\n 72119  1  0\n 72120  1  0\n 73 74  1  0\n 74 75  1  0\n 74 78  1  0\n 74121  1  0\n 75 76  1  0\n 75122  1  0\n 76 53  1  0\n 76 77  1  0\n 77 73  1  0\n 77 81  2  0\n 78 79  1  0\n 79 80  1  0\n 80 75  1  0\n 82 79  2  0\n 83 73  1  0\n 83123  1  0\n 83124  1  0\n 84 78  1  0\n 84125  1  0\n 84126  1  0\nM  END\n'


def get_CB_xyz():
	"""
	RETURNS the xyz coords of the CB7 as aligned for docking
	"""
	return 'N          1.17300       -4.41760        2.95110\nC         -0.07500       -4.63490        3.66880\nC         -0.07650       -3.49110        4.75390\nN          1.17230       -2.78650        4.50310\nC          1.91910       -3.35830        3.47140\nN         -1.31240       -4.40640        2.93690\nC         -2.05220       -3.33530        3.44210\nN         -1.31300       -2.77040        4.48410\nO          3.06910       -3.04680        3.14370\nO         -3.19390       -3.01040        3.09840\nN          1.18540       -0.44340        5.26920\nC         -0.05910       -0.01130        5.88770\nC         -0.04410        1.55190        5.68630\nN          1.20390        1.78870        4.97420\nC          1.94010        0.62010        4.76980\nN         -1.29940       -0.42330        5.24500\nC         -2.03050        0.65390        4.74010\nN         -1.28190        1.81100        4.96560\nO          3.08940        0.55040        4.32120\nO         -3.17410        0.60510        4.27370\nN          1.25270       -3.81360       -3.70770\nC          0.01240       -3.88180       -4.46710\nC          0.03120       -2.55400       -5.31600\nN          1.27790       -1.91530       -4.91890\nC          2.00910       -2.67830       -4.00540\nN         -1.23310       -3.78670       -3.71890\nC         -1.96390       -2.63750       -4.02810\nN         -1.20710       -1.89080       -4.93330\nO          3.15450       -2.43990       -3.60820\nO         -3.11020       -2.37850       -3.64670\nN          1.30960        4.46300       -2.83730\nC          0.07310        5.23030       -2.77760\nC          0.05980        5.78830       -1.30300\nN          1.29610        5.26980       -0.73560\nC          2.04250        4.52450       -1.65040\nN         -1.17630        4.48880       -2.86860\nC         -1.93030        4.55030       -1.69480\nN         -1.18980        5.27500       -0.75900\nO          3.18450        4.08590       -1.47470\nO         -3.08300        4.13030       -1.54540\nN          1.23780        3.86200        3.64110\nC         -0.00160        4.62390        3.69430\nC          0.01760        5.44570        2.34910\nN          1.26310        5.03310        1.71850\nC          1.99350        4.13490        2.49920\nN         -1.24760        3.87520        3.61620\nC         -1.97870        4.16020        2.46070\nN         -1.22240        5.05090        1.69610\nO          3.13950        3.73390        2.26940\nO         -3.12530        3.77450        2.20880\nC          1.86110       -0.76170       -5.58080\nC         -1.87940       -1.75130        5.35010\nC         -1.76030       -0.72970       -5.60640\nC          1.74120       -1.77950        5.38070\nN          1.19820       -5.29810        0.64790\nC         -0.04980       -5.92660        0.23990\nC         -0.03340       -5.78520       -1.33020\nN          1.21380       -5.08130       -1.59320\nC          1.95220       -4.83790       -0.43330\nN         -1.28820       -5.26570        0.62620\nC         -2.02140       -4.80680       -0.47050\nN         -1.27080       -5.07480       -1.61710\nO          3.10340       -4.39170       -0.38280\nO         -3.16440       -4.33840       -0.43910\nC          1.74550       -5.34570        1.99160\nC         -1.87060       -5.32000        1.95580\nC         -1.81570       -4.87070       -2.94710\nC          1.79430       -4.90080       -2.91210\nC          1.78630        3.09970        4.74790\nC         -1.83280        3.13150        4.71800\nC          1.84530        5.66680        0.54850\nC         -1.76680        5.68600        0.50910\nN         -1.17310        2.69850       -4.56390\nC          0.07660        2.57330       -5.29950\nC          0.06730        1.07360       -5.78430\nN         -1.18350        0.55440       -5.25060\nC         -1.92640        1.52260       -4.57170\nN          1.31140        2.67280       -4.53360\nC          2.04410        1.48380       -4.53070\nN          1.30020        0.53150       -5.23130\nO         -3.07880        1.39070       -4.14490\nO          3.18430        1.33100       -4.07910\nC         -1.73190        3.95610       -4.09950\nC          1.88410        3.91890       -4.05570\nH         -0.08170       -5.64240        4.10140\nH         -0.08720       -3.87010        5.78280\nH         -0.07020       -0.31300        6.94190\nH         -0.03930        2.11060        6.62990\nH          0.00540       -4.78930       -5.08230\nH          0.03450       -2.73140       -6.39810\nH          0.09010        6.01950       -3.53880\nH          0.06170        6.88370       -1.25380\nH         -0.00740        5.25800        4.58910\nH          0.02420        6.53120        2.50460\nH          1.77990       -0.89980       -6.66660\nH          2.91210       -0.73560       -5.29180\nH         -2.93010       -1.65570        5.07510\nH         -1.79850       -2.08760        6.39120\nH         -2.81600       -0.68260       -5.33670\nH         -1.66240       -0.87150       -6.69060\nH          2.79740       -1.70010        5.12170\nH          1.63950       -2.11990        6.41920\nH         -0.06570       -6.96910        0.58040\nH         -0.02960       -6.74880       -1.85300\nH          2.80370       -5.09540        1.90860\nH          1.63560       -6.36570        2.38080\nH         -2.92390       -5.05990        1.84720\nH         -1.78050       -6.34490        2.33740\nH         -2.87210       -4.63220       -2.81980\nH         -1.71330       -5.80280       -3.51710\nH          2.85220       -4.68140       -2.76120\nH          1.68720       -5.83670       -3.47480\nH          2.84150        2.93900        4.52460\nH          1.68910        3.69160        5.66700\nH         -2.88620        2.99070        4.47540\nH         -1.74090        3.72890        5.63390\nH          2.90030        5.39200        0.53700\nH          1.74700        6.75450        0.65230\nH         -1.65740        6.77280        0.60990\nH         -2.82470        5.42510        0.47160\nH          0.09490        3.29630       -6.12350\nH          0.07420        0.97040       -6.87590\nH         -2.79190        3.77840       -3.91420\nH         -1.61700        4.70600       -4.89330\nH          1.80270        4.67220       -4.84870\nH          2.93490        3.71960       -3.84360\n'

def make_complex(mol):
	"""
	PRE: Takes in a molecule aligned with the macrocycle
	POST: combines it with the macrocycle and returns a complex molecule
	"""
	macrocycle=Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False)
	complex_CB_guest = Chem.CombineMols(mol, macrocycle)
	return complex_CB_guest

def converge_molecule(mol):
	"""
	Pre: takes  in a molecule and
	Post : converges it and returns its energy in MMFF94
	"""
	n_steps=1000000
	tol=1e-8
	AllChem.MMFFSanitizeMolecule(mol)
	ff = AllChem.MMFFGetMoleculeForceField(mol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94', mmffVerbosity = 0), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
	ff.Initialize()
	cf = ff.Minimize(n_steps, tol, tol)
	return cf, ff.CalcEnergy()

def dock_molecule(mol, converge=False):
	"""
	PRE: Takes in a rdkit molecule
	POST: returns a docked complex with CB7
	"""
	cb = Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False)
	sal=align_mol(mol)[0]
	comp = make_complex(sal)
	if converge:
		converge_molecule(comp)
	return comp

def align_xyz_from_file_to_file(f):
	"""
	PRE: Takes in a xyz file
	POST: will align the corresponding molecule, write it to another file called f+'_aligned'.
		The function will also write another file that corresponds to the complex by concatenating the guest and host files.
	"""
	res = align_mol_from_atm_list_and_xyz(f)

	with open(f+'_aligned.xyz', "wb") as w:
		w.write(str(len(res[0]))+"\n")
		w.write("name is {}\n".format(f))
		for i in range(len(res[0])):
			w.write("{0}\t{1:+2.12f}\t{2:+2.12f}\t{3:+2.12f}\n".format(res[0][i], *res[1][i]))

	with open(f+'_docked.xyz', "wb") as w:
		w.write(str(len(res[0])+126)+"\n")
		w.write("name is {}-docked\n".format(f))
		for i in range(len(res[0])):
			w.write("{0}\t{1:+2.12f}\t{2:+2.12f}\t{3:+2.12f}\n".format(res[0][i], *res[1][i]))
		w.write(get_CB_xyz())
	return

def format_gaussian_input_from_xyz(xyz_file):
	"""
	PRE: Takes in a valid xyz file with a +1 charge implicitly
	POST: will produce a gaussian input file
	"""
	name=xyz_file.split("/")[-1]
	route="#N wB97XD/6-31G(d) opt=(ts, noeigentest, modredundant, calcfc, maxcyc=999) freq=noraman"
	freeze=" D       2       3       9      10 F"
	# route="#N wB97XD/6-31G(d) opt=(ReadOptimize)"
	# route="#N PM6D3 opt=(ReadOptimize)"
	# route="#N PM6D3 opt=(ts, noeigentest, modredundant, calcfc, maxcyc=999) freq=noraman"
	checkpoint="%Chk={}.chk".format(name)
	mem="%mem=120gb"
	procs="%NProcShared=24"
	# route="#N PM6D3 freq=noraman"
	if 'diradical' in f or 'N2' in f:
		charge_mult="1 3"
		print "{} is diradical or n2".format(f)
	else:
		charge_mult = "2 1"
	if 'TS' in f:
		route="#n wB97XD/6-31G(d) opt=(ts, noeigentest, modredundant, calcfc, maxcyc=999) maxdisk=100GB freq"
	else:
		route="#N wB97XD/6-31G(d) opt freq=noraman"

	with open(xyz_file, 'rb') as r:
		coords = r.readlines()[2:]

	with open(xyz_file+".com", "wb") as w:
		w.write(procs+"\n")
		w.write(checkpoint+"\n")
		w.write(mem+"\n")
		w.write(route+"\n\n")
		w.write(name+"\n\n")
		w.write(charge_mult+"\n")
		w.writelines(coords)
		w.write("\n")
		# w.write("notatoms=1-{}\n".format(len(coords)-126))
		w.write(freeze+"\n")
		w.write("\n")

def convert_smiles_to_sdf(smifile):
	"""
	PRE   : Takes in a smiles file
	POST  : Returns a sdf file with all the said molecules
	"""
	w = Chem.SDWriter(smifile[:-4]+".sdf")
	with open(smifile, "rb") as r:
		for line in r:
			smile = line.strip().split()[0]
			name = line.strip().split()[1:]
			smile = smile.split('.')[0]
			print smile, name 
			try:
				mol = Chem.MolFromSmiles(smile)
				mol = Chem.AddHs(mol)
				AllChem.EmbedMolecule(mol)
				AllChem.MMFFOptimizeMolecule(mol)
				mol.SetProp("_Name", "-".join(name))
				w.write(mol)
			except:
				print "NOT TREATING {}".format(name)
	w.close()



if __name__ =="__main__":
	import rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem, Draw
	import numpy as np
	import scipy
	from sklearn.decomposition import PCA
	import sys
	import glob

	mol=Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False )
	Chem.MolToPDBFile(mol, "/home/macenrola/Documents/MACHINE_LEARNING/ledock/CB7_pure.pdb")
# =============================================================================
# 	convert_smiles_to_sdf("/home/macenrola/Documents/MACHINE_LEARNING/ledock/todock_jiaqi.can")
# =============================================================================
# =============================================================================
# 	if len(sys.argv)==1:
# 		flist=[]
# 	elif len(sys.argv)==2:
# 		flist=[sys.argv[1]]
# 	else:
# 		flist=sys.argv[1:]
# 
# # =============================================================================
# # 	flist=glob.glob('/home/macenrola/Documents/DBOA/doubly-charged/base_pdbs/*_docked.xyz')
# # =============================================================================
# 	for f in flist:
# 		# print
# 		format_gaussian_input_from_xyz(f)
# 		# align_xyz_from_file_to_file(f)
# =============================================================================
########## HS COMPENSATION
#	with open("/home/macenrola/Documents/H-S-compensation/C19_sample_200.can", "rb") as r:
#		for i, line in enumerate(r):
#			smi, pb = line.strip().split()
#			try:				
#				print i, pb
#				mol = Chem.MolFromSmiles(smi)
#				mol = Chem.AddHs(mol)
#				AllChem.EmbedMolecule(mol)
#				converge_molecule(mol)
#				comp = dock_molecule(mol, True)
#				wm = Chem.SDWriter('/home/macenrola/Documents/H-S-compensation/sdfs/{}-mol.sdf'.format(pb))
#				wc = Chem.SDWriter('/home/macenrola/Documents/H-S-compensation/sdfs/{}-comp.sdf'.format(pb))
#				wm.write(mol)
#				wc.write(comp)
#				wm.close()
#				wc.close()
#			except:
#				print "Number {} failed".format(pb)
############## HS COMPENSATION
	# cb=Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False)
	# Chem.MolToMolFile(cb, "/home/macenrola/Documents/DBOA/DBOA_TO_DOCK/CB7.sdf")
	# with open('/home/macenrola/Documents/DBOA/DBOA_TO_DOCK/CB7.xyz', 'rb') as r:
	# 	lines=r.read()
	# 	print lines.__repr__()
	# s="C1CC2(CCC1(CC2)C(=O)O)C(=O)O"
	#
	#
	#
	# salenergy = converge_molecule(sal)[1]
	# compenergy = converge_molecule(comp)[1]
	# cbenergy = converge_molecule(cb)[1]
	#
	# print("energy binding is {}".format(compenergy-salenergy-cbenergy))
