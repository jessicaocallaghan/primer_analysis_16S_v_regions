import pandas as pd
import sys

#Usage: python parser.py
if len(sys.argv) < 3:    
    # Show the correct usage.
    print('Usage: python parser.py <blast_tabular_output> <output_file>')
else:
  taxa = {}
  #97_otu_taxonomy file is part of Greenenes bundle
  #parse taxonomy file
  with open('97_otu_taxonomy.txt','r') as tax:
      for line in tax:
          data = line.rstrip().split()
          taxa[int(data[0])] = data[1]
  #read in BLAST results
  frame = pd.read_csv(sys.argv[1], sep= '\t',header = None,index_col = 0)
  #generate taxonomy information for each row
  taxonomy = frame[1].apply(lambda x: taxa.get(x))
  #rename columns
  frame.columns = ['Reference_Name','Percent_ID', 'Alignment_Length','Mismatches','Gap_Openings','Query_Start','Query_End',\
        'Ref_Start','Ref_End','Bit_Score','E_Value']
  #add taxonomy info
  frame['Taxonomy'] = taxonomy
  #write to output file
  frame.to_csv(sys.argv[2], sep='\t')
