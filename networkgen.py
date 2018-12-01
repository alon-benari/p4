from chemUtils import Tanimoto, Pvalue
import argparse
import pandas as pd
from itertools import combinations

class networkgen():
  """
  A class to generate a network of protein based on their p-value from tanimoto score
  """
  def __init__(self,**kwargs): 
    """ set parameters for networkgen.py
    Params
    ------
    Input: kwargs from command line,  which are the drugs.csv from DrugBank, Targets.csv detailing the molecules,
            and the file detailing the proteins to be examined for the generation of the network.
    Output: None
    """
    self.net_params = {}
    for k,v in kwargs.items():
      self.net_params[k] = v
    
    self.data = pd.read_csv(self.net_params['protein_nodes_file'])
    self.data.columns = ['uniprot_id','id', 'pathology']
    self.prot_list = combinations(self.data.uniprot_id,2) # make unique combinations of protein pairs

  def get_pv(self):
    """
    A method to interate over a unique set of protein pairs, compute the p-value for them. If p-value<0.5 than add them to an edge list.

    Parametes
    ----------
    Input: None
    Output: list of tuple forming an edge list.
    """
    #set parameters for Pvalue class
    edge_list = []
    pvalue_params = {}
    pvalue_params['seed'] = 214
    pvalue_params['N'] = 500
    pvalue_params['drug_file'] = self.net_params['drug_file']
    pvalue_params['targets_file'] = self.net_params['targets_file']
    for (p1,p2) in self.prot_list:
      if p1 != p2:
        pvalue_params['proteinA'] = p1
        pvalue_params['proteinB'] = p2
        # Instantiate the Pvalue class
        pval = Pvalue(**pvalue_params)
        # print(pval.get_pvalue())
        if pval.get_pvalue() <= 0.05:
          # print(p1,p2)
          edge_list.append((p1,p2))

    return edge_list



if __name__ =="__main__":
  parser = argparse.ArgumentParser(description = 'Params for pvalue.py for p4')
  #
  # number of iterations
  #
  

  parser.add_argument('drug_file',type = str,help = 'specifies the drug file')

  parser.add_argument('targets_file',type = str,help = 'specify the target file')

  parser.add_argument('protein_nodes_file', type = str,help = 'specify the protein_nodes file')
  
  
  res = parser.parse_args()
  #res = parser.parse_args(['drugs.csv','targets.csv','protein_nodes.csv'])
  kwargs = {k:v for k,v in res._get_kwargs()}
  protein_G = networkgen(**kwargs)
  el = protein_G.get_pv()
  elist = [e1+','+e2+'\n' for e1,e2 in el]

  with open('network_edgelist.txt','w') as e:
    e.writelines(elist)

