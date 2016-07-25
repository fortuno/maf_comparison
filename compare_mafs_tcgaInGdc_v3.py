#!/usr/bin/python

from pyliftover import LiftOver
import sys
import re
import glob
import pybedtools
from pybedtools import BedTool

# Print USAGE
if len(sys.argv) < 3:
   print """\
   There are some missing arguments.
   Usage: compare_mafs MAF_FILE_GDC MAF_FILE_TCGA
        MAF_FILE_1: Path for GDC maf file   
        MAF_FILE_2: Path for TCGA maf file 
   """
   sys.exit()
else:
   gdc_maf_project = sys.argv[1]
   tcga_maf_file  = sys.argv[2]

# Read files in GDC path
gdc_maf_files = glob.glob('../kossproject/*_maf_files_tcga/TCGA.' + gdc_maf_project + '*.maf')
nfiles_gdc = len(gdc_maf_files)

# Read crossing reference
lo = LiftOver('hg19', 'hg38')
fastaRef = pybedtools.example_filename('/mnt/GDCpaper/Homo_sapiens.GRCh38.dna.primary_assembly.fa')

# Variables for count FP, FN, TP, TN
pair_list = {}
TP=0
FP=0
total=0
noncross=0
diffref=0

# Reading each file separately
gdc_var_files_list = [None] * nfiles_gdc
gdc_pairs = []
file = 0
for maf_file in gdc_maf_files:
   # Reading GDC maf file
   gdc_var_list = []
   print "Retrieving variant keys from " + maf_file + " in GDC ..."
   with open(maf_file) as f:
      for line in f:

          # Read columns for each variant in the MAF file
          columns = line.split('\t')

          # Filter empty rows and headers
          if len(columns)>2 and columns[0] != "Hugo_Symbol":

              # Filtering variants in GDC
              # 1) SNPs
              # 2) Filter "PASS" or 'common_mutation'
              if columns[9] == "SNP" and (columns[108] == "PASS" or columns[108] == "common_variant"):        
 
                 samples_pair = ' '.join([columns[15], columns[16]])
                 # gdc_var_list[' '.join([samples_pair, columns[4], columns[5], columns[6], columns[7]])] = columns[10] 
                 gdc_var_list.append(' '.join([columns[4], columns[5], columns[6], columns[7], samples_pair])) 

                 # position = columns[4].replace('chr', '') + ':' + columns[5] + '-' + columns[6]
                 # refbase = BedTool.seq(position, fastaRef) 
                 # print "{0} {1} {2}".format(position, columns[10], refbase)
                       
                 # Check samples in GDC
                 if samples_pair not in gdc_pairs:                
                     gdc_pairs += [samples_pair]

      # Close GDC MAF file
      f.close()

   print "{0} GDC variants considered in {1}".format(len(gdc_var_list), maf_file)
   gdc_var_files_list[file] = gdc_var_list
   file += 1

# Reading TCGA maf file
print "Checking variants keys in TCGA..."
total_variants = 0
with open(tcga_maf_file) as f:
   for line in f:
      
       # Read columns for each variant in the MAF file
       columns = line.split('\t')
  
       # Filter empty rows and headers
       if len(columns)>2 and columns[0] != "Hugo_Symbol":

          pair_key = columns[15] + ' ' + columns[16]

          # Filtering variants in TCGA
          # 1) SNPs
          # 2) This sample comparison exists in GDC
          if columns[9] == "SNP" and pair_key in gdc_pairs:

             start = lo.convert_coordinate('chr' + columns[4], int(columns[5]))
             end = lo.convert_coordinate('chr' + columns[4], int(columns[6]))
             total_variants += 1

             # Check if reference has been correctly crossed
             if start is not None and end is not None and len(start)==1 and len(end)==1:
         
                 refbase = BedTool.seq(start[0][0].replace('chr','') + ':' + str(start[0][1]) + '-' + str(end[0][1]), fastaRef)

                 # Check if reference in TCGA is the same in hg38 ref
                 if refbase == columns[10]:

                     variant_key = ' '.join([start[0][0], str(start[0][1]), str(end[0][1]), start[0][2], columns[15], columns[16]])
                
                     # Create pair if it is not created
                     if pair_key in pair_list:
                         pair_list[pair_key][4] += 1 
                     else:    
                         pair_list[pair_key] = [0] * nfiles_gdc + [1]
                                              
                     # Check if this is a TP in  all gdc files       
                     for i in range(0,nfiles_gdc):
                         if variant_key in gdc_var_files_list[i]:
                             pair_list[pair_key][i] += 1

                     # if pair_list[pair_key][4] % 100 == 0:
		     #     print "PARTIAL TP={0}  TOTAL={1} Recall={2}".format(pair_list[pair_key][0], pair_list[pair_key][4], float(pair_list[pair_key][0])/float(pair_list[pair_key][4])) 
                 else:
                     diffref +=1    
             else:
                  noncross += 1
  
   # Close GDC MAF file
   f.close()   
    
total = 0
tp = [0] * nfiles_gdc
recall_file = gdc_maf_project + '_recall' + '.txt'
f = open(recall_file, 'w')
f.write("Samples\tMuSE\tMuTect2\tSomaticSniper\tVarScan2\tTotal\n")

for key in pair_list:    
    total += pair_list[key][4]
    tp = [tp[0]+pair_list[key][0], tp[1]+pair_list[key][1], tp[2]+pair_list[key][2], tp[3]+pair_list[key][3]] 
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(key, pair_list[key][0], pair_list[key][1], pair_list[key][2], pair_list[key][3], pair_list[key][4])) 

f.close()

print "TOTAL: {0}\t{1}\t{2}\t{3}".format(float(tp[0])/float(total),float(tp[1])/float(total),float(tp[2])/float(total),float(tp[3])/float(total))
print "Cases = {0}".format(len(pair_list))
print "Non Crossed variants = {0}/{1}".format(noncross,total_variants)
print "Ref seq not matching variants = {0}/{1}".format(diffref,total_variants)
