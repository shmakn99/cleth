import numpy as np 



def compare(distribution_original,distribution_random,param='mean'):
	if param == 'mean':
		if np.mean(distribution_original) < np.mean(distribution_random):
			return 1
		else:
			return 0

	if param == 'median':
		if np.median(distribution_original) < np.median(distribution_random):
			return 1
		else:
			return 0

	if param == 'zsco':
		x1=np.mean(np.array(distribution_original))
		x2=np.mean(np.array(distribution_random))
		ssig1=(np.std(np.array(distribution_original))/np.sqrt(len(distribution_random)))**2
		ssig2=(np.std(np.array(distribution_original))/np.sqrt(len(distribution_random)))**2
		return (x1-x2)/np.sqrt(np.abs(ssig1-ssig2))


