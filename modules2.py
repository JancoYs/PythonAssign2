#! /usr/bin/python


#QUESTION 1
def read_fasta_file(filename):
  try:
   fastaFile = open(filename , "r")
  except IOError:
    print"Error, filename ", filename, "not found"
    return -1
  else:  
    print "File was found."
    
    #I am assuming for now that the fasta files you give us contain a single entry, and start with firstline ">'name'"
    line = fastaFile.readline()
    if line[:1] == ">":
      chromoName = line[1:].rstrip('\n')

    else:  
	print "This file is not in the correct format, or is possibly not a fasta file."
	return -1
	
  chromo = ""
  
  for line in fastaFile: 
    chromo += line.rstrip('\n')
  return(chromoName, chromo)
  
  
  
### END OF QUESTION 1 

  
  
# QUSTION 2  
# I have done the reverse_complement in two ways - one using the String translate functions, and one without...
def reverse_complement(dna):
  from string import maketrans
  from string import translate
  
  #reverse the string
  dna2 = dna[::-1]
  
  
  # for various reasona, using this in the translation table creates errors, so must be done seperately
  dna2 = dna.replace("u", "t")
  dna2 = dna2.replace("U", "T")
  
  #create translation table
  trtable = maketrans("ACGTacgt","TGCAtgca") 
  
  # use trtable to translate
  revc = translate(dna2 , trtable)
  print "The reverse complement of", dna, "is", revc
  return revc  
 
 
#Question 2 Option 2
# this second one doesn't use the string translate functions..
def reverse_complement2(dna):

    #just in case,replace the cases when there are U's or u's with T's or t's which should be unlikely.
    #+ we do this first, because otherwise the T's and t's may not get translated into A's... 
    dna2 = dna.replace("u", "t")
    dna2 = dna2.replace("U", "T")
    
    #now we temporarily replace lower case letters with numbers...
    dna2 = dna2.replace("a", "1")
    dna2 = dna2.replace("c", "2")
    dna2 = dna2.replace("g", "3")
    dna2 = dna2.replace("t", "4")
    
    #make all the letters that remain lowercase (otherwise, we get the situation where we replace A's with T's,
    #+ and the later, we replace all the T's with A's, so all A's or T's end up as A's which is not what we want!
    dna2 = dna2.lower();
    
    # make replacements.
    dna2 = dna2.replace("g","C")
    dna2 = dna2.replace("c","G")
    dna2 = dna2.replace("a","T")
    dna2 = dna2.replace("t","A")
    
    # and now replace the numbers (for lowercase) to their complements...
    dna2 = dna2.replace("1","t")
    dna2 = dna2.replace("2","g")
    dna2 = dna2.replace("3","c")
    dna2 = dna2.replace("4","a")
    
    # and now reverse the string
    dna3 = dna2[::-1]
    
    print "The reverse complement of", dna, "is", dna3
    return dna3

### End of Question 2

#Question 3:

def get_ORF(dna):
  
  #for simplicity sake, I make everything upper case:
  dna2 = dna.upper()
  
  #now we find the start codon
  removed = 0 
  length = 0
  startpos = []
  openframe = False
  while (not openframe):
    if dna2[removed:].find("ATG")!=1:
      pos = dna2[removed:].find("ATG")
      frame = pos%3 + removed%3
    else:
      print("There are either no start codons or no ORFs") #this  should break the loop, so no imfinite loops	
      return -1
    for I in range(pos+3,len(dna2)):
      if I%3 == frame:
	  continue
      if ((dna2[I:I+3] == "TAA") or (dna2[I:I+3] == "TGA") or (dna2[I:I+3] == "TAG")):
	orf = frame;
	return orf
	   
  


    
    #END of question 3


### Question 4

def get_gene_by_ORF(dna,orf):
  
  dna2 = dna.upper()
  start = 0
  for x in range(0,len(dna2)):
    if x%3 == orf:
      if (dna2[x:x+3] == 'ATG'): 
	start = x
	print "start = ", start
	break
      
  for y in range((start+3),len(dna2)):
    if y%3 == orf:
      if ((dna2[y:y+3] == 'TAG') or  (dna2[y:y+3] == 'TGA') or (dna2[y:y+3] == 'TAA')):	
	end = y			
	break
  return dna[x:y+3]
  #I wasn't sure if you wanted to include the start and stop codons - if not, one can just "retunr dna[x+3:y]

  
  
# Question 5 (probably could do this with biopython, but I thought I;d make a dictionary)
  
def translate(dna):

  Gdict={ 
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "Stop",
    "TAG": "Stop",
    "TGT": "C",
    "TGC": "C",
    "TGA": "Stop",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"}

  translation = ""
  for a in range(0,len(dna)-3):
    if a%3 == 0:
      translation = translation + (Gdict[(dna[a:a+3])])
  return translation
																																
	
### Question 6

def get_fasta(name, seq):
  fastastr = ">"+name+"\n" 
  
  n=0
  while n  < len(seq):
    fastastr = fastastr + seq[n:n+59]+"\n"
    n += 60
  return fastastr
	
