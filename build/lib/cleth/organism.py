import networkx as nx
import os, errno
import pickle 


class Organism:

	def __init__(self,string_id,name=None):
		self.string_id=string_id

		if name is not None:
			self.name=name

		self.graph=nx.Graph()
		self.essential_proteins=[]


		
	def node_attributes(self,centralities='primary'):
		'''
		Parameters
		-------------

		centralities - The centralities which need to be calculated or need to be looked up from the cache.

		Returns
		-------------

		A dictionary of node attributes or the centralities as specified in the parameter list.

		'''

		if centralities=='primary':

			primary_cents={}
			if os.path.isfile('cache/'+str(self.string_id)+'/centralities/flag.txt'):
				with open('cache/'+str(self.string_id)+'/centralities/primary.pickle', 'rb') as handle:
					primary_cents = pickle.load(handle)

				return primary_cents




			else:
				try:
					os.makedirs('cache/'+str(self.string_id)+'/centralities')
				except OSError as e:
	  			    if e.errno != errno.EEXIST:
	        			raise
				
				with open('cache/'+str(self.string_id)+'/centralities/flag.txt','w') as f:
					f.write('This is a flag file to indicate that this path is actually present')


	        	

				primary_cents['betweenness']=nx.algorithms.centrality.betweenness_centrality(self.graph)
				primary_cents['degree']=nx.degree_centrality(self.graph)
				primary_cents['eigenvector']=nx.eigenvector_centrality(self.graph)
	 
				with open('cache/'+str(self.string_id)+'/centralities/primary.pickle','wb') as handle:
					pickle.dump(primary_cents, handle, protocol=pickle.HIGHEST_PROTOCOL)

				return primary_cents






