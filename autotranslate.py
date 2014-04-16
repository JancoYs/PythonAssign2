#! /usr/bin/python

def autotranslate(inputfile, outputfile):
    import modules
    if modules.read_fasta_file(inputfile)!= -1:
      (chromoname, chromo) = modules.read_fasta_file(inputfile)
    else:
      print("Error in fasta file or file not found")
      return -1
      
    orf = modules.get_ORF(chromo);
    
    gene = modules.get_gene_by_ORF(chrom, orf)
    
    translation = modules.translate(gene)
    
    protname = "protein_" + "chromoname"
    
    protfasta = modules.getfasta(protname,translation)
    
    
    # should be exception handling here.. out of time!
    
    w = open(outputfile, "w")
    w.writelines(protfasta)
    
    print("file created " + outputfile)
    return 0	
