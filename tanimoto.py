import sys
from chemUtils import Tanimoto 



def main():
  drugs_csv =  sys.argv[1]    #'drugs.csv' 
  targets_csv =  sys.argv[2]  #'targets.csv'
  try:
    outputfile = sys.argv[3]  #'outputfile.csv' 
  except:
    outputfile = 'output_file.csv'

  tanimoto = Tanimoto(drugs_csv,targets_csv,outputfile)

  tanimoto.output_tanimoto()
  return tanimoto

if __name__ == '__main__':
  t = main()  