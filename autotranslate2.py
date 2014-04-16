#! /usr/bin/python
import modules2

def autotranslate2(inputfile, outputfile):
    
    # should be exception handling here.. out of time!
    #also more error checking...
    
    
    if ((modules2.read_fasta_file(inputfile))!= -1):
      (chromoname, chromo) = modules2.read_fasta_file(inputfile)
    else:
      print("Error in fasta file or file not found")
      return -1
      
    orf = modules2.get_ORF(chromo)
    
    gene = modules2.get_gene_by_ORF(chromo, orf)
    
    translation = modules2.translate(gene)
    
    protname = "protein_" + chromoname
    
    protfasta = modules2.get_fasta(protname,translation)
    
    # should be exception handling here too... out of time!
    
    w = open(str(outputfile), "w")
    w.writelines(protfasta)
    w.close()
    
    print("file created " + outputfile)
    return 0	
