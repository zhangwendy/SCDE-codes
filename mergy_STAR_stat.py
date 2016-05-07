#!/usr/bin/env python


import glob
import sys
import re
import os
from collections import defaultdict
from argparse import ArgumentParser

def merge(input_dir,outfile,library): 

    gene_dict = defaultdict(dict)
    samples = defaultdict(dict)

    if library == "F":
       col = 2
    elif library == "R":
       col = 3
    else:
       col = 1 
    
    files = glob.glob(input_dir+"*ReadsPerGene.out.tab")

    for file in files:
        name = re.search(r'([^\/]+)_ReadsPerGene.out.tab',file);
        sample = name.group(1)
        samples[sample] = 1
        with open(file, "r") as f:
            for line in f:
                a = line.strip().split()
                if a[0] == "N_unmapped" or a[0] == "N_multimapping" or a[0] == "N_noFeature" or a[0] == "N_ambiguous": 
                    continue
                gene = a[0]             
                gene_dict[gene][sample] = int(a[col])
             

    fout = open(outfile,"w")
    for sample in sorted(samples):
        fout.write("\t"+ sample)
    fout.write("\n")

    for gene in sorted(gene_dict):
        fout.write(gene)
        for sample in sorted(gene_dict[gene]):
            fout.write("\t"+ str(gene_dict[gene][sample])) 
        fout.write("\n")
    fout.close()


def main(argv):

   description = "A script to merge mappable reads to each gene in each samples, and generate a matrix table for SCDE"

   parser = ArgumentParser(description=description)
   parser.add_argument('-i', '--input_dir', required=True, help="Input: a folder containing STAR alignmet files *ReadsPerGene.out.tab, like ./")
   parser.add_argument('-o', '--outfile', required=True, help="Specify a file for the output (Required)")
   parser.add_argument('-l', '--library', required=True, help="RNAseq library type. Three values F (reads map to sense strand), R (reads map to antisense strand), and U (reads map to both strands).")

   args = parser.parse_args()

   if not os.path.exists(args.input_dir):
       parser.error("You must supply a right -i argument for a folder containing STAR alignmet files *ReadsPerGene.out.tab!")

   if not args.outfile:
       parser.error("You must supply the -o argument for output file!")

   if not args.library:
       parser.error("You must supply the -l argument for RNAseq library! Three values F (reads map to sense strand), R (reads map to antisense strand), and U (reads map to both strands)")

   merge(args.input_dir,args.outfile,args.library)


if __name__ == '__main__':
    main(sys.argv)



