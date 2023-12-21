from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os.path
import itertools    
import re
import glob, os
import argparse


parser = argparse.ArgumentParser(description='Compare score')
parser.add_argument('-sp', dest='specie', action='store', help='Get specie')
parser.add_argument('-threshold', dest='threshold', action='store', help='Get threshold')
parser.add_argument('-window', dest='window', action='store', help='Get window size')
parser.add_argument('-utr', dest='utr', action='store', help='Get window size')

def map_int(n):
	try:
		return int(n)
	except ValueError:
		pass
		
def get_score(score_dir,nextPARS_name):
	file_open = score_dir + "/"  + nextPARS_name + ".csv"
	# ~ print(file_open)
	score  = pd.read_csv(file_open,sep=';',header=None,index_col = 0) 
	# Remove last element
	return score.iloc[:, :-1]

def get_positions(specie,target,window,utr):
	file_open = variants_dir + "/" + specie + '_uniq_variation_annot_'+ utr + 'UTR.tab'
	var_annot = pd.read_csv(file_open,index_col = 0,sep='\t')
	
	if  specie == 'tropicalis':
		target = target.replace("exon-", "")[0:-2]

	elif specie == 'parapsilosis':
		target = target.replace("exon-", "rna-")[0:-2]

	var_annot = var_annot[var_annot['Gene']==target]
	
	return var_annot['cDNA_position'].tolist()
	
score_path ='/home/ucielp/ChorosteckiLab/GabaldonLab/projects/MULTI-FOLDS/nextPARS/temperatures/processing/' 
variants_dir = '/home/ucielp/ChorosteckiLab/GabaldonLab/projects/MULTI-FOLDS/nextPARS/temperatures/processing/scripts_candida/miki_data/parse_data/'

def get_data_by_specie(specie,threshold,score_dir,pattern_to_search,window,utr):

	os.chdir(score_dir)
	score_list = list()
	cols = ['All', 'SNP_pos', 'no_SNP_pos']
	df_list = list()

	for file in glob.glob(pattern_to_search):
		base = os.path.basename(file)
		
		# get target
		target = os.path.splitext(base)[0]
		nextPARS_name = target.replace('-E1','')	
		
		score_df = get_score(score_dir,target)

		if (window != 0):
			# Rolling version
			score_df = score_df.T.rolling(window).mean().T
		
		# Filter by score Threshold
		gene_name = score_df.index
		score_df = score_df.T[(score_df.T[gene_name] >= threshold) | (score_df.T[gene_name] <= -threshold)]
		score_df = score_df.T


		pos = get_positions(specie,nextPARS_name,window,utr)
		# Map to int
		positions = list(map(map_int, pos)) 
		
		df_pos = score_df[score_df.columns.intersection(positions)].T
		
		all_positions = range(1,score_df.shape[1]+1,1)
		
		# 
		neg_positions = [x for x in all_positions if x not in positions]
		df_neg = score_df[score_df.columns.intersection(neg_positions)].T
		
		# instead of saving the mean i will save all the data
		# score_list.append((score_df.T[target].mean(),df_pos[target].mean(),df_neg[target].mean()))
		
		
		# Ommit all
		df_neg['SNP'] = 0
		df_pos['SNP'] = 1
		
		df_pos.rename(columns={target:'nextPARS'}, inplace=True)
		df_neg.rename(columns={target:'nextPARS'}, inplace=True)

		frames = (df_pos,df_neg)	
		df_result = pd.concat(frames)
		

		df_result['gene'] = target
		
		# ~ print(df_result)
		
		df_list.append(df_result)
		
	all_df = pd.concat(df_list,ignore_index=True)
	out_csv = variants_dir + '../cinta_output/' + specie + '_' + str(window) + '_' + utr + 'UTR.csv'
	all_df.to_csv(out_csv,index=False)	
	
	return all_df

# Este método es de otro lado pero me puede servir
def get_strand(gene_full,suffix):

	import pyranges as pr

	bed_file = '/home/uchorostecki/lab/uchorostecki/projects/MULTI-FOLDS/nextPARS/temperatures/processing/DB/C_glabatra/candida_glabrata' +suffix  + '.bed'
	
	gr = pr.read_bed(bed_file)
	
	df = gr.df
	gene_row = df.loc[df['Name'] == gene_full]
	return gene_row 

def density_plot(df,specie,window,utr):

	if specie == 'all':
		sns.set_style("whitegrid") 
	
	else:
		df = df.drop('gene', axis=1)
		sns.set_style("whitegrid") 
	
		ax = sns.kdeplot(df.nextPARS[df.SNP == 1],label='SNP');
		ax = sns.kdeplot(df.nextPARS[df.SNP == 0],label='no SNP');

		plt.xlabel('value')
		plt.ylabel('density')

		figure = ax.get_figure()
		figure.set_size_inches(15, 12)

		pdf_out = variants_dir + '../cinta_output/' + specie + "_" + str(window) +  '_' + utr +'UTR.pdf'
		figure.savefig(pdf_out, dpi=400)
		figure.clf()
		
def compare_score(specie,threshold,window,utr):
	
	if specie == 'all':
		
		list_of_species = ('glabrata','parapsilosis')
		frames = list()
		for specie in list_of_species:
			print("ES la especie: ",specie)
			if specie == 'albicans':
				score_dir = score_path + '/score_files_' + specie + '_polyA_one_alelle/'
				pattern_to_search = "*-T-E1.csv"
				# IMPORTANT TODO: modify this
				# pattern_to_search = "C1_09*-T-E1.csv"
			elif specie == 'parapsilosis':
				# TODO: modify the temperature
				score_dir = score_path + '/score_files_' + specie + '_candida_mine/23'
				pattern_to_search = "exon*.csv"
			elif specie == 'tropicalis':
				score_dir = score_path + '/score_files_' + specie + '_candida_mine/'
				pattern_to_search = "exon*.csv"
				# ~ pattern_to_search = "exon-XM_0025451*.csv"
			else:
				# ~ score_dir = score_path + '/score_files_' + specie 
				score_dir = score_path + '/score_files_' + specie + '_' + utr + 'UTR'
				pattern_to_search = "*.csv"
				# IMPORTANT TODO: modify this
				# pattern_to_search = "C1_08*-T-E1.csv"
			score_df = get_data_by_specie(specie,threshold,score_dir,pattern_to_search,window,utr)
			
			# Old plot
			# bar_plot(score_df,specie)
			density_plot(score_df,specie,window,utr)
		
			score_df['specie'] = specie
			frames.append(score_df)
			
		# TODO: check this
		all_df = pd.concat(frames)
		out_csv = variants_dir + '../cinta_output/all_' + str(window) + '_' + utr + 'UTR.csv'
		all_df.to_csv(out_csv,index=False)	
		
		# Ommit this plot
		# bar_plot(all_df,'all')
	
	else:
		if specie == 'albicans':
			score_dir = score_path + '/score_files_' + specie + '_polyA_one_alelle/'
			pattern_to_search = "*-T-E1.csv"
			# TODO: modify this
			# pattern_to_search = "C1_09*-T-E1.csv"
		elif specie == 'parapsilosis':
			# TODO:check other temperatures
			score_dir = score_path + '/score_files_' + specie + '_candida_mine/23'
			pattern_to_search = "exon*.csv"
		elif specie == 'tropicalis':
			score_dir = score_path + '/score_files_' + specie + '_candida_mine/'
			pattern_to_search = "exon*.csv"
		else:
			score_dir = score_path + '/score_files_' + specie 
			pattern_to_search = "*-T-E1.csv"
		

	
args = parser.parse_args()
compare_score(args.specie,float(args.threshold),int(args.window),str(args.utr))
