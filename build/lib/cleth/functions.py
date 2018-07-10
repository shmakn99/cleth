import numpy as np
import random as rd
import os, errno
import pickle 
from matplotlib import pyplot as plt

from auxiliary import compare

def predict_essential(organism,cut_off=90,centralities='primary'):
	'''
	Parameters
	-------------

	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	cut_off - The percentile cut-off to be considered while selecting the essential proteins.

	centralities - The set of centralties or node (protein) attributes to be used based on which the proteins 
	will be segregated.

	Returns
	-------------

	A dictionary. Where the keys are the centralities or the node (protein) attributes considered and the values 
	are the list of essential proteins for each centraility.

	'''

	predicted_essential_cent_wise={}

	predicted_essential_concur=[]

	primary_cent=organism.node_attributes(centralities=centralities)

	for cent_name in primary_cent.keys():

		current_cent=primary_cent[cent_name]

		cut_off_value=np.percentile(list(current_cent.values()),cut_off)

		predicted_essential_cent_wise[cent_name] = [protein for protein in current_cent.keys() if current_cent[protein] > cut_off_value]

	return predicted_essential_cent_wise



def percentilewise_distribution(organism,centralities='primary',steps=10):
	'''
	Parameters
	-------------

	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	centralities - The set of centralties or node (protein) attributes to be used based on which the proteins 
	will be segregated.

	steps - The steps in which percentile brackets will increase.
	For example if steps=10, the percentiles consedered will be - 10,20,30 ... 90

	Returns
	-------------

	A dictionary. Where the keys are different centralities and the values are dictionaries themselves.

	Where the keys are different cut-offs and the values are the percentage of nodes (proteins) which are above the 
	cut-off and are essential. The higher the number for higher percentile the better the centrality measure is. 


	'''

	if len(organism.essential_proteins) == 0:
		print ('ERROR: Can Not calculate percentile wise distribution without the actual essentiality data.')
		return 0

	primary_cent = organism.node_attributes(centralities=centralities)

	centralitywise_dict={}

	if os.path.isfile('cache/'+str(organism.string_id)+'/percentilewise_distribution/flag'+centralities+str(steps)+'.txt'):

		with open('cache/'+str(organism.string_id)+'/percentilewise_distribution/'+centralities+str(steps)+'.pickle', 'rb') as handle:
			centralitywise_dict=pickle.load(handle)

		return centralitywise_dict

	else:
		try:
			os.makedirs('cache/'+str(organism.string_id)+'/percentilewise_distribution')
		except OSError as e:
	  		if e.errno != errno.EEXIST:
	  			raise

		with open('cache/'+str(organism.string_id)+'/percentilewise_distribution/flag'+centralities+str(steps)+'.txt','w') as f:
			f.write('This is a flag file to indicate that this path is actually present')

		for cent_name in primary_cent.keys():

			percentilewise_distributed_dict={}

			current_cent=primary_cent[cent_name]

			for cut_off in range(steps,100,steps):

				cut_off_value=np.percentile(list(current_cent.values()),cut_off)

				proteins_above_cut_off = [protein for protein in current_cent.keys() if current_cent[protein]>cut_off_value]
				ess_proteins_above_cut_off = [protein for protein in proteins_above_cut_off if protein in organism.essential_proteins]

				if len(proteins_above_cut_off) == 0:
					percentilewise_distributed_dict[cut_off]=0

				else:
					percentilewise_distributed_dict[cut_off]=len(ess_proteins_above_cut_off)/len(proteins_above_cut_off)

			centralitywise_dict[cent_name]=percentilewise_distributed_dict

		with open('cache/'+str(organism.string_id)+'/percentilewise_distribution/'+centralities+str(steps)+'.pickle','wb') as handle:
			pickle.dump(centralitywise_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


		return centralitywise_dict


def validate(organism,centralities='primary',niter=100,param='mean',zsco_threshold=2.33,get_from_cache=False):
	'''
	Parameters
	------------

	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	centralities - The set of centralties or node (protein) attributes to be used for which the pvals will
	be calculated

	param - The statistical test used for comparing the significance of in silico predictions.

	'mean' - The mean of the node properties of proteins is compared with that of randomly obtained population. 
	'median' - The median of the node properties of proteins is compared with that of randomly obtained population.
	'zsco' - z-scores are obtained for comparision and a threshold value is used to segregate significant distributions.

	get_from_cache - If true a pre calculated test value is supplied. The values for every might differ slightly
	as there are random processed involved.

	niter - The number of times the test is performed.

	Returns
	-------------

	A dictionary where the keys are the are different centralities and the values are the p-values based on whatever
	parameter is supplied for testing. 


	'''

	if len(organism.essential_proteins)==0:
		print ('ERROR: Can not validate as list of essential protiens not provided')

		return 0

	if get_from_cache is True and os.path.isfile('cache/'+str(organism.string_id)+'/pvals/flag'+centralities+param+str(niter)+'.txt'):
		with open('cache/'+str(organism.string_id)+'/pvals/'+param+str(niter)+'.pickle', 'rb') as handle:
			centralitywise_pval = pickle.load(handle)

		return centralitywise_pval


	else:
		try:
			os.makedirs('cache/'+str(organism.string_id)+'/pvals')
		except OSError as e:
	  		if e.errno != errno.EEXIST:
	  			raise

		with open('cache/'+str(organism.string_id)+'/pvals/flag'+centralities+param+str(niter)+'.txt','w') as f:
			f.write('This is a flag file to indicate that this path is actually present')


		primary_cent=organism.node_attributes(centralities=centralities)

		centralitywise_pval={}

		for cent_name in primary_cent.keys():

			centralitywise_pval[cent_name]=0
			current_cent=primary_cent[cent_name]
			comparisions=[]

			for i in range(niter):

				distribution_essential = [current_cent[protein] for protein in current_cent.keys() if protein in organism.essential_proteins]
				randomly_chosen_proteins = rd.sample(list(current_cent.keys()),len(distribution_essential))
				distribution_random = [current_cent[protein] for protein in randomly_chosen_proteins]

				comparisions.append(compare(distribution_essential,distribution_random,param=param))


			if param == 'zsco':
				centralitywise_pval[cent_name]=len([1 for instance in comparisions if instance>zsco_threshold])
			else:
				centralitywise_pval[cent_name]=len([1 for instance in comparisions if instance>0])

		with open('cache/'+str(organism.string_id)+'/pvals/'+centralities+param+str(niter)+'.pickle','wb') as handle:
			pickle.dump(centralitywise_pval, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return centralitywise_pval


def graphit(organism,centralities='primary',steps=10):
	'''
	Parameters
	------------

	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	centralities - The set of centralties or node (protein) attributes to be used based on which the proteins 
	will be segregated.

	steps - The steps in which percentile brackets will increase.
	For example if steps=10, the percentiles consedered will be - 10,20,30 ... 90

	Returns
	------------

	Nothing. Publishes graphs to a directory.
	'''

	if len(organism.essential_proteins) == 0:
		print ('ERROR: Can Not calculate percentile wise distribution without the actual essentiality data.')
		return 0

	percentilewise_dstrbtn=percentilewise_distribution(organism,centralities=centralities,steps=steps)


	if os.path.isfile('cache/'+str(organism.string_id)+'/graphs/flag'+centralities+str(steps)+'.txt'):
		print ('Graphs for this organism with this step value are already present')

		return 0

	else:
		try:
	   		os.makedirs('cache/'+str(organism.string_id)+'/graphs')
		except OSError as e:
	  		if e.errno != errno.EEXIST:
	  			raise

		
		for cent_name in percentilewise_dstrbtn.keys():
			curr_dist=percentilewise_dstrbtn[cent_name]

			fig=plt.figure()
			plt.scatter(list(curr_dist.keys()),list(curr_dist.values()))
			fig.suptitle(cent_name+' '+str(organism.string_id))
			plt.xlabel('Percentile of nodes\' degree')
			plt.ylabel('Fraction of essential nodes with a higher degree')
			plt.xlim(0,100)
			plt.ylim(0,1)
			fig.savefig('cache/'+str(organism.string_id)+'/graphs/steps='+str(steps)+' '+cent_name)

		with open('cache/'+str(organism.string_id)+'/graphs/flag'+centralities+str(steps)+'.txt','w') as f:
			f.write('This is a flag file to indicate that this path is actually present')

		return 0
