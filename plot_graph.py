import networkx as nx
import argparse 
import pandas as pd
import matplotlib.pyplot as plt
from chemUtils import PlotGraph

if __name__ =='__main__':
  parser = argparse.ArgumentParser(description = 'Params for plot_graph for p4')
  parser.add_argument('edgelist_file',type = str,help = 'specifies the edgelist file')
  parser.add_argument('protein_nodes_file', type = str,help = 'specify the protein_nodes file')
  parser.add_argument('output_graph',default = 'network.png',type = str,help = 'specify the outputfile with a .png suffix')
  res = parser.parse_args()
  #res = parser.parse_args(['network_edgelist.txt','protein_nodes.csv','p4_graph.png'])
  kwargs = {k:v for k,v in res._get_kwargs()}

  plotG = PlotGraph(**kwargs)
  plotG.draw_graph()
  # plt.rcParams["figure.figsize"] = [8,8]

  plt.savefig(kwargs['output_graph'],dpi = 150)
  plt.show()
  #
  #[{'from':s,'to':t} for s,t in plotG.G.edges()]
  #[plotG.G.node[n] for i,n in enumerate(plotG.G.nodes())]

