#!/usr/bin/python

import os, re, random, sys, argparse
import csv
from Bio import SeqIO


import forgi.utilities.stuff as fus
import forgi.graph.bulge_graph as fgb

# Depricated
EXTEND = 50


def main():
	parser = argparse.ArgumentParser(description="Get stem loop prediction from sequence\n",
									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	parser.add_argument('-s', '--stemsize', dest = 'stemsize', action = 'store', default = 25,
						help = "Define the stem size")
						
	parser.add_argument('-m', '--mismatches', dest = 'mismatches', action = 'store', default = 10,
						help = "Determine the max number of mismatches in the stem")
						
	parser.add_argument('-ms', '--mismatche_size', dest = 'mm_size', action = 'store', default = 6,
						help = "Determine the max size of mismatches in the stem")
	
	parser.add_argument('--multiple', dest = 'multiple', action = 'store_true', default = False,
						help = "Allow multiple suboptimal secondary structures within a user defined energy range above the minimum free energy ")
	
	parser.add_argument('--energy', dest = 'energy', action = 'store', default = 1,
						help = "Compute suboptimal structures with energy in a certain range of the optimum (kcal/mol).")
						
	parser.add_argument('-i', '--input', dest = 'input', action = 'store', default = 'predicted.fold',
						help = "Input file name (EX: regions.fa)")
						
	parser.add_argument('-o', '--outDir', dest = 'outDir', action = 'store', default = '../results/',
						help = "Input file name (EX: regions.fa)")
	
	parser.add_argument('-d', '--inDir', dest = 'inDir', action = 'store', default = '../data/',
						help = "Directory containing input files")
	
	options = parser.parse_args()


	filename = options.inDir + options.input

	read_file(filename,options)
	

def read_file(filename,options):

	x = ''	
	
	base=os.path.basename(options.input)

	
	file_name = str(os.path.splitext(base)[0]) + '_stem' + str(options.stemsize) + '_mm' + str(options.mismatches) + "_bs" +  str(options.mm_size) + ".csv"
	
	data = ['name','hairpin_loop_size','stem_length','sum_mm','max_mm_size']
	
	file_name_all_structures =  options.outDir +  file_name
	
	# all structures
	create_file_and_header(file_name_all_structures,data)
	

	with open(filename, 'r') as f:
		for line in f:
			line = line.strip()
			# ~ if line.startswith('>'):
			if line.startswith('>'):
				name = line
				x = 'name'
			elif line.startswith('(') or line.startswith('.') :
				fold = line
				x = 'fold'
			else: 
				seq = line 
				x = 'seq'
			if (x == 'fold'):
				stem(name,fold,seq,options,file_name_all_structures)

def create_file_and_header(file_name,data):
	with open(file_name, 'w') as csvfile:
		spamwriter = csv.writer(csvfile, delimiter='\t')
		if (data):
			spamwriter.writerow(data)
	
	
def stem(name,fold,seq,options,file_name_all_structures):
	

	all_structures = list()
	
	bg = fgb.BulgeGraph.from_dotbracket(fold)

	hloop_positions = list()
	
	
	for hloop in bg.hloop_iterator():
		(i,j) = bg.get_flanking_region(hloop)
		hloop_pos = (i+j)/2
		hloop_positions.append(hloop_pos)
	
	original_name = name.replace('>','').split(" ")
	
			
	for hloop_pos in  hloop_positions:
		found = False
		
		stem_length = 0
		mm_sizes = list()
		mm_positions = list()
		last_element_type = ''
		hairpin_loop_size=0
		nro_loops_ant = 0
		
		for element in bg.iter_elements_along_backbone(hloop_pos):
			
			element_type = element[0]
			# hairpin loop
			if (element_type == 'h'):
				hairpin_loop_size = bg.get_length(element)
				continue
				
			# Interior loop - mismatches
			elif(element_type == 'i'):
				
				i,j = bg.get_bulge_dimensions(element)
				if (i == 0):
					# Ommit this because is already count in iter_elements_along_backbone
					mm_sizes.append((bg.get_length(element)*2))
				else:
					mm_sizes.append(bg.get_length(element))
				if (j == 0):
					continue

			# Stem 
			elif(element_type == 's'):

				stem_length = stem_length + bg.get_length(element)
				
				
				a,b = bg.get_side_nucleotides(element, 0)
				x = a - 1 
				y = b

			elif(element_type == 'm'):
				# "Multi-loops"
				p,q =  bg.flanking_nucleotides(element)

		
				x = min(a,b) - 1
				y = max(a,b)
				break
			else:
				# print "Warning", bg.get_length(element)
				found = False
				break
			last_element_type = element_type
			

			loops,nro_loops_ant = int_loop(fold[x:y], seq[x:y],nro_loops_ant)

			
			if (loops):
				mm_sizes.append(loops)

			all_mm = mm_sizes
			
				
			if (stem_length >= int(options.stemsize) ):
				if (sum(all_mm[:-1]) <=  int(options.mismatches)):
					if (all_mm[:-1] and (max(all_mm[:-1]) <= int(options.mm_size))):
						found = True
					else:
						if (found):
							break
						else:
							continue
				else:
					if (found):
						break
					else:
						continue
			else:
				if (found):
					break
				else:
					continue
		

		if (found):
			if (all_mm[:-1]):
				max_mm_size =  str(max(all_mm[:-1]))
			else:
				max_mm_size =  "0"
			
			data = [original_name[4],hairpin_loop_size,stem_length,sum(all_mm[:-1]),max_mm_size,seq[x:y],fold[x:y]]
			
			print(data)
			all_structures.append(data)
			
			
	for st in all_structures:
		write_to_csv(file_name_all_structures,st)
		
		hairpin = st[1]
		stem = st[2]
		mm = st[3]
		seq = st[5]
		struct = st[6]


def write_to_csv(file_name,data):

	# append
	with open(file_name, 'a') as csvfile:
		spamwriter = csv.writer(csvfile, delimiter='\t')
		spamwriter.writerow(data)

			
def int_loop(fold,seq,nro_loops_ant):

	foldStr = ''.join([str(elem) for elem in fold])

	bg = fgb.BulgeGraph.from_dotbracket(foldStr)
	

	loops = list()
		
	actual_loops = 0
	
	for iloop in bg.iloop_iterator():
		i,j = bg.get_bulge_dimensions(iloop)

		if (i == 0):
			# Ommit this because is already count in iter_elements_along_backbone
			continue
			
		if (j == 0):

			loops.append(float(i))
	
	actual_loops = sum(loops)
	if (nro_loops_ant != actual_loops):
		diff = actual_loops - nro_loops_ant
		return diff,actual_loops
	else:
		return 0,nro_loops_ant

		

if __name__ == "__main__":
	main()
