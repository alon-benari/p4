import sys
import pandas as pd
import numpy as np
from random import sample, choice
import networkx as nx
from itertools import combinations
import csv 
from matplotlib import pyplot as plt


class Tanimoto():
    """
    A class with methods to implement calculation of tanimoto score, and generate a bootstrap distribution
    """
    def __init__(self, drugs = 'drugs.csv',targets = 'targets.csv',output = 'output_file.csv'):
      self.cutoff = 0.5
      self.message = 'hello'
      self.protein_ligand = {} # a dictionary to hold the ligans
      self.protein_ligands_num ={} # dictionary to hold protein and the number of ligands associated with it
      drgs = open(drugs)
    
       
      self.drugs_bitmap = {l.split(',')[0]: l.split(',')[2].strip('\n').split() for l in drgs.readlines()} # A dictionary to hold the 
      
      #  = {for k.v}
      #
      self.targets = targets
      self.map_protein_to_ligand()
      self.output = output

    def make_hist(self,hist_type):
      """
      A method to render the histograms for P4. The method will output a .png file according to the hist_type passed to it a parameter.
      Params
      ------
      Input: hist_type, (integer) 1- All,2-share target 3-not shared
      Output: None.
      """
      SunetId = 'abenari'
      data = pd.read_csv(self.output) # read data
      data.columns = ['ID1','ID2','score','shared']
      if (hist_type == 2):
        data = data[data.shared == 1]
        fname = 'shared__tanimoto.png'
        suffix = 'Shared'
      elif (hist_type == 3):
        data = data[data.shared == 0]
        fname = 'notshared__tanimoto.png'
        suffix = 'Not Shared'
      elif (hist_type ==1):
        suffix = 'All'
        fname = 'all_tanimoto.png'
      binwidth = (data.score.max() - data.score.min())/float(np.sqrt(data.shape[0]))
      bins = np.arange(data.score.min(),data.score.max(),binwidth) # array of bins
      fig, ax = plt.subplots()
      ax.set_title(SunetId+suffix)
      ax.set_xlabel('tanimoto score')
      ax.set_ylabel('count')
      data['score'].hist(bins = bins)
     
      plt.savefig(fname)
      plt.show()
    

    
    def map_protein_to_ligand(self):
      """
      Method to map ligands that interact with a given protein.
      This method will populate a dictionary keys by  protein and value is a list of ligands
      It will also create a directed graph  from protein to ligands. 
      where the source is the protein and the target is the small molecule ligand.

      Parameters:
      -----------
      input: None
      output: None
      """
      targets = pd.read_csv(self.targets)
      uniq_targets = (targets.uniprot_accession).unique().tolist()
      
      #initialize dictionary
      
      for protein in uniq_targets:
        line = targets[targets.uniprot_accession == protein].db_id.tolist()
        self.protein_ligand[protein] =  line #[int(item) for item in line.strip('\n').split()]
        self.protein_ligands_num[protein] = len(line)
      
      # 
      self.ligand_graph = nx.DiGraph()  # graph object to hold protein to ligand relationship ()
      for id in targets.index:
          self.ligand_graph.add_edge(targets.iloc[id].uniprot_accession,targets.iloc[id].db_id)


      

    def get_tanimoto(self, molA, molB):
      """
      A method to compute the Tanimoto score between two molecules 
      Parametes:
      ----------
      Input : molA, molB , molecules between which we wish to compute the said score
      Output: Tanimoto score for molA nd molB

      """
      intersect = set(self.drugs_bitmap[molA]).intersection(set(self.drugs_bitmap[molB]))
      union = set(self.drugs_bitmap[molA]).union(set(self.drugs_bitmap[molB]))
      return  len(intersect)/float(len(union))
    

    def tanimoto_summary(self,protA,protB):
      """
      A method to  compute the summrary of teh the Tanimoto score for two proteins given all the ligands tha interact with them

      Params:
      -------
      input: protA, prot B (strings) names of molecules for who the summary needs be computed.
      output: A summary score of all Tanimoto score greater than a predetrmined cutoff (in this case 0.5)
      """
      tn_sum = []
      #  
      for i in self.ligand_graph.neighbors(protA):
        for j in self.ligand_graph.neighbors(protB):
          tn_sum.append(self.get_tanimoto(i,j))
      return  np.array(tn_sum)[np.array(tn_sum )> self.cutoff].sum()
    
    def output_tanimoto(self):
      """
      A method to output the tanimoto score for any unique pair of ligand molecules
      
      Parameters
      ------------
      Input: None
      Output: None
      """
      rows = []
      uniq_pairs = combinations(self.drugs_bitmap.keys(),2)  # create a generator with all unique pairs
      for (lig1,lig2) in uniq_pairs:
        rows.append([lig1,lig2,format(self.get_tanimoto(lig1,lig2),'.6f'),self.share_target(lig1,lig2)])

      # write to csv file
      with open(self.output,mode = 'w') as output_file:
        tanimoto_writer = csv.writer(output_file,delimiter = ',')
        tanimoto_writer.writerows(rows)
      
      return rows
    def share_target(self,molA, molB):
      """
      A method to check if two ligands share a common protein target

      Parameters
      ----------
      Input: molA , molB (strings)  ligands  to examine if shared target exists
      output: binary, 1 - if at least common target exists 0- otherwise.

      """
      
      try:
         setA = set(self.ligand_graph.predecessors(molA))
      except:
        setA = set()
     
      try:
        setB = set(self.ligand_graph.predecessors(molB))
      except:
        setB = set()

      if (len(setA.intersection(setB))>=1):
        flag = 1
      else:
        flag = 0 
      return flag


    def get_distrib(self,protA,protB):
      """
      A method to create a bootstrap distribution with replacement.
      This will create a sample by sampling a random set of ligands the size of ligands for A, and a similar set for B
      We compute the Tanimoto score for  both and create a distribution.

      Parameters
      ----------
      Input: protA, protA (string), names of proteins for whom we wish to compute the empirical  distribution
      output: list of tanimoto scores for the 
      """
      summary_dist = [] 
      # find how many ligand molecules are associated with protA  and protB, nA and NB respectively.
      nA = self.ligand_graph.out_degree(protA) # how many out edges from node protA
      nB = self.ligand_graph.out_degree(protB) # same for proteinB
      
        

      for ligand1 in [choice(list(self.drugs_bitmap.keys())) for _ in range(nA)]:
        for ligand2 in [choice(list(self.drugs_bitmap.keys())) for _ in range(nB)]:
          summary_dist.append(self.get_tanimoto(ligand1,ligand2)) 
      
      return summary_dist

#################################Class for pvalue.py  #####################################
from chemUtils import Tanimoto

from random import seed, random
import numpy as np

class Pvalue():
  """
  Main Method to calculate pvalue for the Tanimoto summary
  """
  def __init__(self, **kwargs):

    self.params = {}
    for key, value in kwargs.items():
      self.params[key] = value
    seed(self.params['seed'])
  
    self.tanimoto_obj = Tanimoto(self.params['drug_file'],self.params['targets_file'])

  def get_summary(self):
    """
    A method to compute the Tanimoto summary for protein A and protein B.
    It leverages the proper method from the Tanimoto Object.
    
    Params:
    ------
    Input:None
    output:return the T_summary
    """
    p1 = self.params['proteinA']
    p2 = self.params['proteinB']
    return (self.tanimoto_obj.tanimoto_summary(str(self.params['proteinA']),str(self.params['proteinB'])))


  def get_pvalue(self):
    """
    A method to generate a distibution form N random samples of nA ligands for protein A and nB ligands for protein B. It uses the method get_distirib 

    Params
    ------
    Input: None
    Output: p-value as defined in stem
    """
    count = 0
    T_summary = self.get_summary()
    #
    for k in range(self.params['N']):
      sample = np.array(self.tanimoto_obj.get_distrib(self.params['proteinA'],self.params['proteinB']))
      if sample[sample > 0.5].sum() > T_summary:
        count += 1 
      
      
    return count/float(self.params['N'])



###################################PlotGraph####################################################
import networkx as nx
import argparse 
import pandas as pd
import matplotlib.pyplot as plt


class PlotGraph():
  """
  A class with methods to render an un-directed graph where nodes are proteins and egdes exist if the tanimoto score between them is <0.5
  """
  def __init__(self,**kwargs):
    """
    A method to initialize  the class, with CLI parameters and other useful data structures
    """
    self.net_params = {}
    for k,v in kwargs.items():
      self.net_params[k] = v  
    #
    # node color
    self.node_color = {"bp":"red","bp;cholesterol":"green","bp;cholesterol;diabetes":"blue", "bp;diabetes":"purple"}

    self.el = pd.read_csv(self.net_params['edgelist_file'],sep = ',')
    self.el.columns = ['origin','target']
    self.edgelist = [(self.el.ix[i,:]['origin'],self.el.ix[i,:]['target']) for i in self.el.index]
    #
    #
    self.nodes = pd.read_csv(self.net_params['protein_nodes_file']).set_index('uniprot_accession').to_dict()
    #
    self.node2color ={k:self.node_color[v] for k, v in self.nodes['indications'].items()}
    self.node2label = self.nodes['uniprot_id']
    #
    # Build the Graph
    #
    self.G = nx.Graph()
    self.G.add_nodes_from([n for n in self.node2color.keys()])
    for n in self.G.nodes:
      self.G.node[n]['color'] = self.node2color[n]
      self.G.node[n]['label'] = self.node2label[n]
      self.G.nodes[n]['indications'] = self.nodes['indications'][n]
    self.G.add_edges_from(self.edgelist)
  
  def draw_graph(self):
    """
    A method to output the graph image
    """
    pos = nx.kamada_kawai_layout(self.G)
    fig,ax = plt.subplots(figsize = (8,8))
    
  
    nx.draw_networkx_nodes(self.G,pos,ax=ax,
                          nodelist = self.G.nodes(),
                          node_color  = [self.G.node[n]['color'] for n in self.G.nodes()]
                           )
    nx.draw_networkx_edges(self.G,pos,width = 1.0)
    nx.draw_networkx_labels(self.G,pos,self.node2label,font_size=10)
      


