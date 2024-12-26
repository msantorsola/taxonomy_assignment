import pandas as pd
from pandas_ml import ConfusionMatrix


out = pd.read_csv('PAN_output_TRUE_predicted_lineages.csv', sep='\t')
out_blast = pd.read_csv('Blastn_output_TRUE_predicted_lineages.csv', sep='\t')
out_rdp = pd.read_csv('RDP_output_TRUE_predicted_lineages.csv', sep='\t')



###for each rank

###PAN
out_95 = out[(out.Likelihood_Weight_Ratio >= 0.95)]

actual =out_95.loc[:,'species_TRUE'].astype(str).tolist()
predicted = out_95.loc[:,'species'].astype(str).tolist()

cm = ConfusionMatrix(actual, predicted)
A = cm.stats()
A = pd.DataFrame(A.items()[2][1])
del A['nan']

A['mean'] = A.mean(axis=1)


###Blast
out_blast95 = out_blast[(out_blast['bit score'] >= 2946.32)]

actual =out_blast95.loc[:,'family_TRUE'].astype(str).tolist()
predicted = out_blast95.loc[:,'family'].astype(str).tolist()

cm = ConfusionMatrix(actual, predicted)
A = cm.stats()
A = pd.DataFrame(A.items()[2][1])
del A['nan']

A['mean'] = A.mean(axis=1)


####RDP

out_rdp95 = out_rdp[(out_rdp.confidence_phylum >= 0.95)]

actual =out_blast95.loc[:,'phylum_TRUE'].astype(str).tolist()
predicted = out_blast95.loc[:,'phylum'].astype(str).tolist()

cm = ConfusionMatrix(actual, predicted)
A = cm.stats()
A = pd.DataFrame(A.items()[2][1])
del A['nan']

A['mean'] = A.mean(axis=1)

