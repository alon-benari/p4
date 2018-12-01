# from chemUtils import Tanimoto
# import argparse
# import sys
# from random import seed, random
# import numpy as np

# class Pvalue():
#   """
#   Main Method to calculate pvalue for the Tanimoto summary
#   """
#   def __init__(self, **kwargs):

#     self.params = {}
#     for key, value in kwargs.items():
#       self.params[key] = value
#     seed(self.params['seed'])
  
#     self.tanimoto_obj = Tanimoto(self.params['drug_file'],self.params['targets_file'])

#   def get_summary(self):
#     """
#     A method to compute the Tanimoto summary for protein A and protein B.
#     It leverages the proper method from the Tanimoto Object.
    
#     Params:
#     ------
#     Input:None
#     output:none
#     """
#     p1 = self.params['proteinA']
#     p2 = self.params['proteinB']
#     return (self.tanimoto_obj.tanimoto_summary(str(self.params['proteinA']),str(self.params['proteinB'])))


#   def generate_sample(self):
#     """
#     A method to generate a distibution form N random samples of nA ligands for protein A and nB ligands for protein B. It uses the method get_distirib 

#     Params
#     ------
#     Input: None
#     Output: p-value as defined in stem
#     """
#     count = 0
#     T_summary = self.get_summary()
#     #
#     for k in range(self.params['N']):
#       sample = np.array(self.tanimoto_obj.get_distrib(self.params['proteinA'],self.params['proteinB']))
#       count += (sample>T_summary).sum()
#     return count/float(self.params['N'])
      
from chemUtils import Pvalue
import argparse
import argparse
import sys
import time

if __name__ =="__main__":
  parser = argparse.ArgumentParser(description = 'Params for pvalue.py for p4')
  #
  # number of iterations
  #
  parser.add_argument('-n' ,
                      action = 'store',dest = 'N', type = int,
                      required = False,default = 500,
                      help = 'determine sample size of bootstrap distribution')

  parser.add_argument('-r',action = 'store',dest = 'seed', default = 214,type = int,
                      help = 'a parameter to set the seed for the samplig (default 214)')

  parser.add_argument('drug_file',type = str,help = 'specifies the drug file')

  parser.add_argument('targets_file',type = str,help = 'specify the target file')

  parser.add_argument('proteinA', type = str,help = 'specify the uniprot_accession ID for proteinA')
  
  parser.add_argument('proteinB', type = str,help = 'specify the uniprot_accession ID for proteinB')
  
  res = parser.parse_args()
  #res = parser.parse_args(['-n' ,'100','-r', '214','drugs.csv','targets.csv','P21918','P18089'])
  kwargs = {k:v for k,v in res._get_kwargs()}
  pvalue = Pvalue(**kwargs)
  T_summary = pvalue.get_summary()
  p_value = pvalue.get_pvalue()
  print(p_value)
