from Bio import AlignIO
from Bio.Seq import Seq
import pandas as pd
import glob
import seaborn as sns
import numpy as np 
import matplotlib.pyplot as plt
import subprocess
import glob, os
from pathlib import Path
from Bio import SeqIO
import argparse
import shutil
import re


# Records_id
def get_seq_ids(aln):
	seq_ids = list()
	for record in aln:
		seq_ids.append(record.id)
	return seq_ids

# Get only aligned positions
def get_aligned_positions(aln):
	aln_pos = list()
	for i in range(aln.get_alignment_length()):
		# Aligned positions
		if len(set(aln[:, i])) == 1:
			aln_pos.append(i)
	
	# I've return only position of proteins
	return aln_pos

def column_from_residue_number(aln, id, res_no):
	rec = next((r for r in aln if r.id == id), None)
	j = 0
	for i, res in enumerate(rec.seq):
		if res!='-':
			if i==res_no:
				return j
			j+=1

# Get dictionary with aligned positions relative to nextPARS score
def dic_from_aligned_positions(aln):
	aln_pos = get_aligned_positions(aln)
	old_seq_id = ''
		
	# Generate empty dictionary with gene names 
	res = dict.fromkeys(get_seq_ids(aln), "")
	
	for seq_id in get_seq_ids(aln):
		# List of corrected positions
		corrected_pos = [column_from_residue_number(aln,seq_id,aln_p) for aln_p in aln_pos]
		# Convert protein aligned positions to nucleot
		t  = [[x*3,x*3+1,x*3+2] for x in corrected_pos]
		aln_pos_to_nucl = [item for sublist in t for item in sublist]	
		
		if(seq_id.startswith("transcript")):
			old_seq_id = seq_id
			seq_id = seq_id.replace('transcript_','transcript:')
			del res[old_seq_id]
		
		res[seq_id] = (aln_pos_to_nucl)
	return res

	
######################
##### ALIGNMENT ######
######################

def tcoffee(in_fasta_path,out_path_MSA,file_name):
	
	temp_file = in_fasta_path + 'temp.fasta'
	shutil.copyfile(in_fasta_path + file_name + '.fasta',temp_file )

	output = subprocess.run(['t_coffee', temp_file + " -method=slow_pair -output=clustalw_aln " + "-run_name="+out_path_MSA + file_name ], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

	return output.stdout
	
def perc_identity(aln):
    i = 0
    for a in range(0,len(aln[0])):
        s = aln[:,a]
        if s == len(s) * s[0]:
            i += 1
    return 100*i/float(len(aln[0]))

def get_strand(gff_path,name):
	
	os.chdir(gff_path)

	ext_name = "ID=" + name + ";"
	gff_files = []
	for file in glob.glob('*gff'):
		gff_files.append(file)
	
	for gff in gff_files:
		file = open(gff, "r")
		for line in file:
			# ~ if re.search(name, line):
			if re.search(ext_name, line):
				return(line.split("\t")[6])

def get_st_end_positions(gff_path,name,st_end):
	
	os.chdir(gff_path)
	res = list()	
	if(name.startswith("transcript")):
		ext_name = "ID=" + name + ";"
	else:
		ext_name = name
	gff_files = []
	for file in glob.glob('*gff'):
		gff_files.append(file)
	
	for gff in gff_files:
		file = open(gff, "r")
		
		for line in file:
			if re.search(ext_name, line):
			# ~ if re.search(ext_name, line):
				if (st_end):
					# Get end
					res.append(line.split("\t")[4])
				else:
					# Get start
					res.append(line.split("\t")[3])
	return res

######################
####### DF ###########
######################

	
def plot_correlation(df,to_save,csv_save,csv_all_save,csv_all_save_shuffle,perc_id):
	
	# Fix AttributeError: 'PandasArray' object has no attribute '_str_len'
	pd.set_option("display.max_columns", None)
	
	# SHUFFLE DATASET
	# Copy	
	df_shuffle = df
	# Shuffle and reset index
	shuffle_row = df_shuffle[df_shuffle.columns[0]].sample(frac=1)
	shuffle_row.index = (list(range(0,len(shuffle_row))))
	
	df_shuffle = df_shuffle.drop(df_shuffle.columns[0], 1)
	df_shuffle['shuffle'] = shuffle_row
	

	# monotonic relationship
	corr = df.corr(method='spearman')
	corr_shuffle = df_shuffle.corr(method='spearman')
	
	
	corr.to_csv(csv_save,index=False)
	
	sns.set(style="white")
	mask = np.zeros_like(corr, dtype=bool)
	mask[np.triu_indices_from(mask,k=1)] = True
	ax = plt.axes()
	cmap = sns.diverging_palette(10, 220, as_cmap=True)

	
	# Values of correlation
	svm = sns.heatmap(corr, mask=mask, ax = ax, cmap=cmap, square=True, linewidths=.5,vmin=0, vmax=1,xticklabels=True,yticklabels=True)

	svm.set_yticklabels(svm.get_yticklabels(), rotation = 0, fontsize = 12)
	svm.set_xticklabels(svm.get_xticklabels(), rotation = 90, fontsize = 12)

	figure = svm.get_figure()    


	# fix for mpl bug that cuts off top/bottom of seaborn viz
	b, t = plt.ylim() # discover the values for bottom and top
	b += 0.5 # Add 0.5 to the bottom
	t -= 0.5 # Subtract 0.5 from the top
	plt.ylim(b, t) # update the ylim(bottom, top) values
	
	figure = plt.gcf() # get current figure
	figure.set_size_inches(12,12)
	figure.savefig(to_save)

	plt.close()
	
	# Return the values of correlation matrix in the descending order, not diagonal, not NA
	ds = (corr[corr < 1].unstack().transpose().sort_values( ascending=False).drop_duplicates().dropna())
	print("Corr",ds)
	
	# Return the values of correlation matrix in the descending order, not diagonal, not NA for shuffled dataset
	ds_shuffle = (corr_shuffle[corr_shuffle < 1].unstack().transpose().sort_values( ascending=False).drop_duplicates().dropna())

	# Add percentage of identity
	perc_df = pd.DataFrame([perc_id])
	
	if ds.empty:
		print('Warning: DataFrame is empty!')
	else:
		perc_df.to_csv(csv_all_save, mode='a',header=False,index=False,line_terminator=',',)
		ds.to_csv(csv_all_save, mode='a',header=False)
		
		# Save randomized dataset
		perc_df.to_csv(csv_all_save_shuffle, mode='a',header=False,index=False,line_terminator=',',)
		ds_shuffle.to_csv(csv_all_save_shuffle, mode='a',header=False)
	
######################
#### PROCESS #########
######################

def main(in_path_fasta,csv_all_save,csv_all_save_shuffle,gff_path,threshold):
	
	all_corr = list()
	os.chdir(in_path_fasta)
	perc_id = 0
	

	for file in glob.glob("*.fasta"):		



		fasta_proteins = list()
		multi_fasta_name = Path(file).stem
		fasta_len = dict()
		
		# Get fasta lengths
		input_file = in_path_fasta + multi_fasta_name + ".fasta"
		fasta_sequences = SeqIO.parse(open(input_file),'fasta')

		last_fasta_name =''
		last_fasta_seq =''
		
		# Fasta sequnce only have one of the exons, but it is OK, then I will catch it before to transform to protein
		# And I will use that information for the nextPARS
		for fasta in fasta_sequences:
			name, sequence = fasta.id, fasta.seq
			
			strand = get_strand(gff_path,name)
			# if(strand == '-'):
			
			# Get positions to know if there is one or more protein (XXX-P) withe the same name 
			st_pos = get_st_end_positions(gff_path,name,0)
			if (len(st_pos) == 1):
				if len(sequence) %3 ==1:
					sequence = sequence + Seq('NN')
				elif (len(sequence)) %3 == 2:
					sequence = sequence + Seq('N')
			
				fasta_proteins.append(">" + name)
				fasta_proteins.append(str(sequence.translate()))
				
				fasta_len[name] = len(sequence)

			elif(len(st_pos)==2) :
				print("Sequence " + name + " is duplicated, merge the two sequences")

				next_element = next(fasta_sequences)
			
				# Append the next sequence, at the end if strand is +
				if(st_pos[0] < st_pos[1] ):
					if(strand == '+'):
						merged_seq = sequence + next_element.seq
					else:
						merged_seq = next_element.seq + sequence
				
				# Append the sequence, at the begining if strand is +
				else:
					if(strand == '+'):
						merged_seq = next_element.seq + sequence
					else:
						merged_seq = sequence + next_element.seq
				
				if len(merged_seq) %3 ==1:
					fasta_proteins = merged_seq + Seq('NN')
				elif (len(merged_seq)) %3 == 2:
					fasta_proteins= merged_seq + Seq('N')
				
				fasta_proteins.append(">" + name)
				fasta_proteins.append(str(merged_seq.translate()))
			
				fasta_len[name] = len(merged_seq)

			elif(len(st_pos)==3) :
				print("Sequence " + name + "is duplicated, merge the three sequences")

				next_element = next(fasta_sequences)
				last_element = next(fasta_sequences)
			
				# Append the next sequence, at the end if strand is +
				if(st_pos[0] < st_pos[1] ):
					if(strand == '+'):
						merged_seq = sequence + next_element.seq + last_element.seq
					else:
						merged_seq = last_element.seq + next_element.seq + sequence
				
				# Append the sequence, at the begining if strand is +
				else:
					if(strand == '+'):
						merged_seq = last_element.seq + next_element.seq + sequence
					else:
						merged_seq = sequence + next_element.seq + last_element.seq
				
				if len(merged_seq) %3 ==1:
					fasta_proteins = merged_seq + Seq('NN')
				elif (len(merged_seq)) %3 == 2:
					fasta_proteins= merged_seq + Seq('N')
				
				fasta_proteins.append(">" + name)
				fasta_proteins.append(str(merged_seq.translate()))
				
				fasta_len[name] = len(merged_seq)
			
			else:
				print("Warning I have to solve this issue where are more than 3 proteins with the same name")
		
		# Generate proteins alignments		
		in_path_fasta_protein = in_path_fasta + "proteins/"
		protein_fasta = in_path_fasta_protein + multi_fasta_name + ".fasta"
		f = open(protein_fasta, "w")
		f.write('\n'.join(fasta_proteins))
		f.close()
	
		out_path_MSA_proteins = out_path_MSA + "proteins/"
		
		# Generate alignments and output to file
		tcoffee(in_path_fasta_protein,out_path_MSA_proteins,multi_fasta_name)


		# Get alignment in clustal format
		clustal_file_name = out_path_MSA_proteins + multi_fasta_name + ".clustalw_aln"
		try:
			aln = AlignIO.read(clustal_file_name, "clustal")
		except:
			print("File: ",clustal_file_name," is missing")
			
		
		# Calculate percent identity between the already-aligned sequences
		perc_id = round(perc_identity(aln), 2)

		print("Processing " , perc_id, multi_fasta_name)

		# Get dictionary with aligned positions relative to nextPARS score
		result_dic = dic_from_aligned_positions(aln)
		# ~ print(result_dic)

		df_list = list()
		index_names = list()

		# Here I have to get all the nextPARS scores for all the exons
		
		# Only one score file
		for key in result_dic:
			# Tcoffe modified replace : by _
			if(key.startswith("transcript")):
				old_key = key
				key = key.replace("transcript:","CDS:")
				# Avoid to look for YEL009C-A
				key = key.replace("_mRNA",";")
				
			st_pos = get_st_end_positions(gff_path,key,0)
			end_pos = get_st_end_positions(gff_path,key,1)
			strand = get_strand(gff_path,key)
			
		
			if(key.startswith("CDS:")):
				transcript = old_key.replace('transcript_','transcript:')
				make_sure_transcript = old_key.replace('transcript_','transcript:')
				key = transcript
			else:
				transcript = key.replace('-P','-T')
				make_sure_transcript = key.replace('-P','-T-')

			st_pos_transcript = get_st_end_positions(gff_path,make_sure_transcript,0)
			end_pos_transcript = get_st_end_positions(gff_path,make_sure_transcript,1)
			
			
			if (len(st_pos) == 1):
				text_files = glob.glob(nextPARS_path + "/**/RNN/" + transcript + "*.tab", recursive = True)
				if not text_files:
					print ("Warning: file not found, skiping...")
					continue

				specie = text_files[0].split("/")[-3]
				name = specie + "_" + key
				index_names.append(name)

				nextPARS_file = text_files[0]
				df = pd.read_csv(nextPARS_file, header=None, sep=';',index_col=0).dropna(axis='columns', how='all')
				
				# Protein start and transcript start don't match, probably because UTRs
				if st_pos[0] != st_pos_transcript[0]:
					if(strand == '+'):
						diff_beg = int(st_pos[0]) - int(st_pos_transcript[0])
						diff_end = int(end_pos_transcript[0]) - int(end_pos[0])
					# Swap beg and end
					else:
						diff_end = int(st_pos[0]) - int(st_pos_transcript[0])
						diff_beg = int(end_pos_transcript[0]) - int(end_pos[0])

								
					# Drop first N columns of dataframe
					df = df.iloc[: , diff_beg:]
					
					# Drop last N columns of dataframe
					df = df.iloc[: , :-diff_end]

				# Rename and put -1 to the first element to be consistent with the 0-base index
				df.columns = list(range(0,len(df.columns)))
				
				# Check lengths of fasta and nextPARS
				# ~ # If they are not the same is because nextPARS is from an exon instead hole trasncript
				if (len(df.columns) != fasta_len[key]):
					print ("Warning: length are not matching ",key,fasta_len[key],len(df.columns))

				# Filter position
				try:
					filtered_df = df[result_dic[key]]
				except:
					print("An exception occurred, most probably because length of fasta and nextPARS are not equal",len(df.columns),fasta_len[key])
				
				
				# Rename positions from 0 to length of matched positions
				filtered_df.columns = list(range(0,len(filtered_df.columns)))

				df_list.append(filtered_df)
				
			# At least one of the two genes, have more than one exon. Concatanate score files
			elif(len(st_pos)==2) :
				# ~ print("B - two gene")
				print (key,transcript)
				# TODO: improve this
				text_file1 = glob.glob(nextPARS_path + "/**/RNN/" + transcript + "-E1.RNN.tab", recursive = True)
				text_file2 = glob.glob(nextPARS_path + "/**/RNN/" + transcript + "-E2.RNN.tab", recursive = True)
				
				if not text_file1 and not text_file2:
					print ("Warning: file",transcript," E1 and E2 not found, skiping...")
					continue
				
				elif not text_file1:
								
					specie = text_file2[0].split("/")[-3]
					print ("Warning: file",transcript," E1 not found, skiping...")
					nextPARS_file = text_file2[0]
					df = pd.read_csv(nextPARS_file, header=None, sep=';',index_col=0).dropna(axis='columns', how='all')
					
					# I have UTR but one is missing (E1 is the first one)
					if st_pos[0] != st_pos_transcript[0]:
						if(strand == '+'):
							diff_end = int(end_pos_transcript[1]) - int(end_pos[1]) # 
						else:
							diff_end = int(st_pos[0]) - int(st_pos_transcript[0])

						# Drop last N columns of dataframe (3 UTR) - from 2nd T
						df = df.iloc[: , :-diff_end]

					# Add NA at the begining
					diff_length = fasta_len[key] -len(df.columns)
					# Rename so empty diff_length space at the begining 
					df.columns = list(range(diff_length+1,len(df.columns)+diff_length+1))
					df = df.reindex(list(range(1,fasta_len[key]+1)), axis="columns")

				# E2 is missing, fill with Nan at the end
				elif not text_file2:

					# TODO, I have to add if st_pos[0] != st_pos_transcript[0]: as before
					# I didn't do it because for the moment there is no case like this
					specie = text_file1[0].split("/")[-3]
					print ("Warning: file",transcript," E2 not found, skiping...")
					nextPARS_file = text_file1[0]
					df = pd.read_csv(nextPARS_file, header=None, sep=';',index_col=0).dropna(axis='columns', how='all')
					
					# Add NaN at the end until complete fasta length
					df = df.reindex(list(range(1,fasta_len[key]+1)), axis="columns")
					
				# I have to concatanate both files
				else:
					specie = text_file1[0].split("/")[-3]
					nextPARS_file1 = text_file1[0]
					nextPARS_file2 = text_file2[0]
					df1 = pd.read_csv(nextPARS_file1, header=None, sep=';',index_col=0).dropna(axis='columns', how='all')
					df2 = pd.read_csv(nextPARS_file2, header=None, sep=';',index_col=0).dropna(axis='columns', how='all')
					
					# Rename df2 to start at the end
					try:
						df2.columns = list(range(len(df1.columns)+1,fasta_len[key]+1))
					except:
						print("Warning Glabrata and in this case T are distinct than P")
						continue
					
					# Rename index to merge equal -T
					df1.index = [transcript] 
					df2.index = [transcript]
					
					frames = [df1.T, df2.T]
					df = pd.concat(frames,ignore_index=True)
					df = df.T
					
				name = specie + "_" + key
				index_names.append(name)

				# Rename and put -1 to the first element to be consistent with the 0-base index
				df.columns = list(range(0,len(df.columns)))
					

				# Check lengths of fasta and nextPARS
				if (len(df.columns) != fasta_len[key]):
					print ("Warning: length are not matching ",key,fasta_len[key],len(df.columns))

				# Filter position
				try:
					filtered_df = df[result_dic[key]]
				except:
					print("An exception occurred, most probably because length of fasta and nextPARS are not equal",len(df.columns),fasta_len[key])
				
				# Rename positions from 0 to length of matched positions
				filtered_df.columns = list(range(0,len(filtered_df.columns)))

				df_list.append(filtered_df)
				
			elif(len(st_pos)==3) :
				print("Warning they are three")
			elif(len(st_pos)==0) :
				print("Warning empty positions",make_sure_transcript)
			else :
				print("Warning they are many ")
	
		try:
			all_df = pd.concat(df_list,ignore_index=True)
		except:
			print("No object to concatanate, most probably because both files are missing")
			continue
		
		all_df.index = index_names

		save_pdf = save_path + multi_fasta_name + ".pdf"
		save_csv = save_path + multi_fasta_name + ".csv"

		# Drop values using threshold
		all_df = all_df[(all_df[all_df.columns] >= threshold) | (all_df[all_df.columns] <= -threshold)]
		all_df = all_df.dropna(axis=1)


		if all_df.empty:
			print('Warning: Filtered dataFrame is empty!')
		else:
			plot_correlation(all_df.T,save_pdf,save_csv,csv_all_save,csv_all_save_shuffle,perc_id)
		print()

		# REMOVE THIS HACK
		# ~ all_df.to_csv('/home/uchorostecki/lab/uchorostecki/projects/MULTI-FOLDS/nextPARS/temperatures/processing/scripts_candida/good_information/df_temp1.csv',index=False)

		


######################
###### MAIN ##########
######################		


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-t","--type",dest="type",action="store",required=True, default="None",help="msa or pair")
args = parser.parse_args()

# Define main paths
path = "PATH"
nextPARS_path = path + "nextPARS_score/"
gff_path = path + "orthologs_metaphors/gff_files/"

threshold = 0.2
# threshold = 0

if __name__=="__main__":
	if args.type == "msa":
		print("MSA")
		in_path_fasta = path + "orthologs_metaphors/fasta/multi_fasta_orthologs/"
		save_path = path + "orthologs_metaphors/MSA/proteins/plots/"
		out_path_MSA = path + "orthologs_metaphors/MSA/"
		
		csv_all_save = path + "orthologs_metaphors/MSA/proteins/all.csv"
		csv_all_save_shuffle = path + "orthologs_metaphors/MSA/proteins/all_shuffle.csv"
		# Clear file
		open(csv_all_save, 'w').close()
		open(csv_all_save_shuffle, 'w').close()
		main(in_path_fasta,csv_all_save,csv_all_save_shuffle,gff_path,threshold)

	elif args.type == "pair" :
		print("pair")

		in_path_fasta = path + "orthologs_metaphors/fasta/pair_fasta_orthologs/"
		save_path = path + "orthologs_metaphors/PSA/proteins/plots/"
		out_path_MSA = path + "orthologs_metaphors/PSA/"
	
		csv_all_save = path + "orthologs_metaphors/PSA/proteins/all.csv"
		csv_all_save_shuffle = path + "orthologs_metaphors/PSA/proteins/all_shuffle.csv"
		# Clear file
		open(csv_all_save, 'w').close()
		open(csv_all_save_shuffle, 'w').close()
		main(in_path_fasta,csv_all_save,csv_all_save_shuffle,gff_path,threshold)
	else:
		print("msa or pair?")
    
