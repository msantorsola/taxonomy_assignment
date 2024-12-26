#!/usr/bin/python

import json
import re
from cStringIO import StringIO
import pandas as pd
from Bio import Phylo
from pandas.io.json import json_normalize
import numpy as np


#### read RAxML tree output
with open('RAxML_portableTree.LTPs128_SkinMicrobiome_7aprile.jplace') as json_file:
	data = json.load(json_file)

pl = pd.io.json.json_normalize(data, 'placements')
pl.rename(columns={pl.columns[0]: 'Query', pl.columns[1]:'Placements'}, inplace = True)

pl2 = pd.DataFrame()
for row in pl.index:
	col = pl.ix[row, 'Placements']
	N = len(col)
	pl2 = pl2.append(pd.DataFrame([[i, pl.ix[row,'Query']] for i in col],
    	index=[row] * N, columns = ['Placements','Query']))

# add 'edge_num' and 'like_weight_ratio' columns
pl2['edge_num'] = [i[0] for i in pl2['Placements']]
pl2['like_weight_ratio'] = [i[2] for i in pl2['Placements']]
pl2['Query'] = [i[0] for i in pl2['Query']]

edge_num_list = pl2['edge_num'].tolist()

del pl2['Placements']
pl2.set_index("edge_num", inplace=True)
del pl2.index.name

# list of edge_num to search in the ref_tree (json)
ref_tree = data['tree']
ref_tree = re.sub(r"[\}]", "]", ref_tree)
ref_tree = re.sub(r"[\{]", "[", ref_tree)

handle = StringIO(ref_tree)
ref_tree = Phylo.read(handle, "newick")

def get_parent(tree, child_clade):
	node_path = tree.get_path(child_clade)
	return node_path[-2]

parents = {}
for i in set(edge_num_list):
	ins_node = ref_tree.find_elements(comment=str(i)).next()
	parent = get_parent(ref_tree, ins_node)
	parents.setdefault(parent,[]).append(i)

### Get Silva full lineages for Bacteria

import xml.etree.ElementTree as ET

tree = ET.parse('All_SILVA_fullLineage_5maggio.xml')
root = tree.getroot()
taxa = root.findall("./ROOT/taxon")

taxonomy={}
#lineageToTree={}
for i in taxa:
	tax_lin=[]
	for j in i.findall("./lineage/taxon"):
    	if "rank" in j.attrib:
        	tax_lin.append((j.attrib["scientificName"],j.attrib["rank"]))
	taxonomy[i.attrib["taxId"]]=tax_lin

import time

### Silva lineages to pandas dataframe

DF=pd.DataFrame()
frames=[]
counter=0
old=A=time.time()
for t,tt in taxonomy.items():
	counter+=1
	if (counter/100)==(counter/100.0):
    	new=time.time()
    	print new-A
    	print new-old
    	old=new
    	print counter
    	print DF.shape
	DFtemp=pd.DataFrame(tt)
	DFtemp.set_index(1,inplace=True)
	DFtemp.columns=[t]
	DFtemp=DFtemp.transpose()
	frames.append(DFtemp)


DF = pd.concat(frames)
DF=DF[["superkingdom", 'phylum', 'subphylum','class', 'subclass','order','suborder','family','subfamily', 'tribe', 'genus','species group', 'species subgroup','species' ,'subspecies']]


silva = pd.read_table('LTPs128_SSU_taxId.csv', sep='\t', header=None)
silva = silva.iloc[:,[0, 1]]
silva.columns = ['Accession', 'TaxID']
silva['TaxID'] = silva['TaxID'].astype(str)

DF2 = DF.merge(silva, left_index=True, right_on='TaxID', copy=True, indicator=False, how='outer')
DF2.set_index('Accession', inplace=True)
del DF2['TaxID']

parentsName = {}
for i in parents:
	refname=DF2.loc[i.get_terminals()[0].name,:].tolist()
	value=[set(DF2.loc[x.name,:].tolist()) for x in i.get_terminals()]
	CommonTaxonName=set.intersection(*value)
	CommonTaxonName.discard(np.nan)
	parentsName[i] = sorted(CommonTaxonName, key=lambda x: refname.index(x))[-1]
	parentsName[i] = parentsName[i], DF2.columns[refname.index(parentsName[i])]


taxon = [[parentsName[x]]*len(e) for x, e in parents.items()]
taxon2, rank = zip(*sum(taxon,[]))
assign_df = pd.DataFrame({"edge_num": sum(parents.values(),[]), "Assignment":taxon2, "Rank":rank})
assign_df.set_index("edge_num", inplace=True)

output = pl2.merge(assign_df, left_index=True, right_index=True, copy=True, indicator=False)

grouped = output.groupby(['Query', "Assignment", "Rank"])['like_weight_ratio'].sum()

assignment = grouped.to_frame()
assignment = assignment.reset_index()
assignment.set_index('Query', inplace=True, drop=False)

ByQueryAssignment=assignment.groupby(["Query"])

phyloT_taxonomy = Phylo.read('phyloT_generated_tree_1495789075_newick.txt', 'newick')

##aggregating likelihood values for parent nodes by NCBI taxonomy tree

Bits=[]
for Query, ByTaxon in ByQueryAssignment:
	if (ByTaxon.shape[0]!=1):
    	TestN=[phyloT_taxonomy.find_elements(name=re.sub("\)","",re.sub("\(","",re.sub("[ ]", "_",i)))).next() for i in ByTaxon.Assignment]
    	AggregatedScore=[sum([ByTaxon.like_weight_ratio[ByTaxon.Assignment==i.name.replace('_', ' ').replace('MAC', "(MAC)")] for i in TestN if x.get_path(i) != None]) for x in TestN]
    	AggregatedValue= [i.values for i in AggregatedScore]
	else:
    	AggregatedValue=ByTaxon.like_weight_ratio
	ByTaxon["AggregLikeWeightRatio"]=AggregatedValue
	Bits.append(ByTaxon)

AssignmentF=pd.concat(Bits)

AssignmentF['Likelihood_Weight_Ratio'] = AssignmentF['AggregLikeWeightRatio'].str.get(0).astype(float)
AssignmentF.Likelihood_Weight_Ratio.fillna(AssignmentF['AggregLikeWeightRatio'], inplace=True)

del AssignmentF['like_weight_ratio']
del AssignmentF['Query']
del AssignmentF.index.name
del AssignmentF['AggregLikeWeightRatio']

report_html = AssignmentF.to_html(justify = 'left')
report_file = open('PAN_taxonomy_assignment.html','w')
report_file.write(report_html)
report_file.close()

