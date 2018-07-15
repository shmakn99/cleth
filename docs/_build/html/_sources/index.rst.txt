Overview
========
CLeth is a Python package which implements 'Centrality Lethality Hypothesis'. The package renders
organisms as weighted graphs of protein protein interaction scores and tries to elicit essential
proteins using centrality measures as computed using the NetworkX package in Python. 

Installation
============
The recommended way is to install from PyPI (virtual environment recommended):

	>>> pip install cleth

Usage
=====
One can either predict essential proteins or check the significance of proteins predicted by the algorithm
using statistical tests by providing a list of experimentally obtained essential proteins. 

Organism
--------

This is a Python class which contains the String ID, Name, List of essential proteins and a graph of all the interactions as provided in the PPI database. An organism variable can be declared using the load module.

	>>> import cleth.cleth as cl
	>>> from cl import _io
	>>> E_coli=_io.load(57271,ppi_infile_path,threshold,ess_infile_path=None,name=E_coli)

This will load the model in the variable E_coli.
User needs to provide the threshold, interaction above the threshold score will only be considered.

Node Attributes
---------------

	>>> E_coli.node_attributes(centralities='primary')

This will return a dictionary of primary centralities (degree, betweeness and eigen-vector).
This command will also create a pickle file containing this dictionary in the cache directory which will be present in the current directory. 

Predict Essential Proteins
--------------------------

	>>> from cl import functions as f
	>>> predicted = f.predict_essential(E_coli,cut_off=90,centralities='primary')

Parameters
	
	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	cut_off - The percentile cut-off to be considered while selecting the essential proteins.

	centralities - The set of centralties or node (protein) attributes to be used based on which the proteins 
	will be segregated.

Returns

	A dictionary. Where the keys are the centralities or the node (protein) attributes considered and the values 
	are the list of essential proteins for each centraility.

Percentile Wise Distribution
----------------------------

	>>> distribution = f.percentilewise_distribution(E_coli,centralities='primary',steps=10)



Parameters

	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	centralities - The set of centralties or node (protein) attributes to be used based on which the proteins 
	will be segregated.

	steps - The steps in which percentile brackets will increase.
	For example if steps=10, the percentiles consedered will be - 10,20,30 ... 90

Returns

	A dictionary. Where the keys are different centralities and the values are dictionaries themselves.

	Where the keys are different cut-offs and the values are the percentage of nodes (proteins) which are above the 
	cut-off and are essential. The higher the number for higher percentile the better the centrality measure is. 


A pickle file containing the returned dictionary is also saved in the cache directory.

Validate
--------

	>>> p_vals = f.validate(E_coli,centralities='primary',niter=100,param='mean',zsco_threshold=2.33,get_from_cache=False)


Parameters

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

	A dictionary where the keys are the are different centralities and the values are the p-values based on whatever
	parameter is supplied for testing. 

A pickle file containing the returned dictionary is also saved in the cache directory.

Graphing the Results
--------------------

	>>> f.graphit(organism,centralities='primary',steps=10)

Parameters

	organism - object of class organism, it brings with itself the node_attrunites function used in
	the function.

	centralities - The set of centralties or node (protein) attributes to be used based on which the proteins 
	will be segregated.

	steps - The steps in which percentile brackets will increase.
	For example if steps=10, the percentiles consedered will be - 10,20,30 ... 90

Returns

	Nothing. Publishes graphs to cache directory.
