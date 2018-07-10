from organism import Organism
import networkx as nx

def load(string_id,ppi_infile_path,threshold,ess_infile_path=None,name=None):
	'''
	Parameters
	------------
	string_id - The STRING database ID of the orgaism.

	ppi_infile_path - Path of the file containing the protein protein interaction (PPI) data.
	For this package this file is assumed to be obtained from STRING database.
	The format is-
	Potein1 Protein2 Score1 Score2 ... Overall Score

	ess_file_path - Path of the file containing the list of essential protiens. 
	For this package this file is assumed to be obtained from
	The format is-
	Protein 1
	Protein 2
	.
	.
	.
	Protein N

	threshold - The cut-off to be consedered for making edges between any two proteins.
	For example if the threshold is 700, all the PPIs below this overall score are neglected.
	This would result in what is generally refered to as a high confidence network. 

	name - Taxonomical name of the organism.

	'''

	G=nx.Graph()

	with open(ppi_infile_path) as f:
		line=''
		while True:
			line=f.readline()
			if line=='':
				break

			split_line=line.strip().split()

			if int(split_line[len(split_line)-1])>=threshold:
				G.add_edge(split_line[0],split_line[1],weight=int(split_line[len(split_line)-1]))

	essential_protiens=[]

	if ess_infile_path is not None:
		with open(ess_infile_path) as f:
			line=''
			while True:
				line=f.readline()
				if line=='':
					break

				split_line=line.strip().split()

				# print (split_line)

				essential_protiens.append(split_line[0])

	if name is not None:
		org=Organism(string_id,name=name)
	else:
		org=Organism(string_id)

	org.graph=G

	# print (essential_protiens)

	org.essential_proteins+=essential_protiens

	return org









