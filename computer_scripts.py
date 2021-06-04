#!/usr/bin/python3
#tested on version: Python 3.7.9


#import of packages
#general packages
import numpy as np
import pandas as pd
import re
import sys
import os
import os.path
import glob
import itertools
from collections import Counter
import abc
import requests
import time
import gzip
from functools import reduce
from igraph import *

#graphics packages
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import venn #(from https://github.com/tctianchi/pyvenn)
import cairo
from mpl_toolkits.axes_grid.axislines import SubplotZero
#from wordcloud import WordCloud

#machine learning and statistics packages
from sklearn.ensemble import RandomForestClassifier #scikit-learn version 0.19.1 used
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import learning_curve
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster import hierarchy
from sklearn import metrics
from sklearn.base import clone
from scipy.spatial.distance import cdist
from scipy.stats import fisher_exact
from scipy.stats import pearsonr
from math import floor, ceil
from random import sample, choices

#GUI package
from tkinter import *
from tkinter.filedialog import askdirectory, askopenfilename
from tkinter import messagebox


#Following there are various classes to handle database communication, extraction of experimental data from files,
#generation of the Training set and its subsets, machine learning predictions and evaluations, network analysis and analysis of data
#all generated classes:
'''
PrepareData					class to prepare folder structure and data for analysis
OverwriteException			exception in try-except expressions
IDTranslations				contains methods to translate protein and metabolite IDs
DBHandler					instance that trims STITCH database to score cutoff
FileHandler					instance for extraction of experimental data
DBFileInteractions			instance constructing intersections between DBHandler and FileHandler
XSet						general class containing common methods for "Set"-classes, as feature engineering etc.
XSet_sampled				class to generate sampled Positive or Negative Set by finding always positive/negative predicted protein-metabolite pairs with random training
NegativeSet					class to construct negative subset from PubChem
PositiveSet					class to construct positive subset from STITCH
PositiveSet_PubChem			class to construct positive subset from PubChem
TrainingSet					class to construct training set from subsets, including creation of gold standard
CompleteSet					class that constructs feature engineered profiles for all pairs in experimental data to make predictions on
MLClassifiers				Machine learning application and evaluation
NetworkAnalysis				Analysis of network from predictions
DataAnalysis				Analysis of experimental data, predictions and main construction of figures for thesis
'''

#Often occuring attributes of classes:
'''
self.experimental			folder with experimental data files
self.databases				folder with experiment-specific databases and subsets
self.analyses				folder with analyses of experiment
self.normalized				string with normalisation method for profiles
self.feature				string describing feature engineering approach
self.proteinwise			boolean: True for prediction of interactions on protein-metabolite pairs, False for OG-metabolite pairs
self.protcol				coloumn name depending on self.proteinwise, either "protein" or "OG"
self.approach				string describing balancing and confidence in training set
self.method					profile pooling method
'''
###############################################################
###############################################################
###############################################################


#prepares folder structure and splits excel sheet with datasets into separate files for analysis
#Folder structure used to handle program in-/outputs
'''
in a workspace folder, construct a subfolder with the program script(s) manually, from which all programmatic action is executed.
A folder with databases (including self generated ones for ID-translations) called "databases",
a folder "overall_analysis" for the overall_analysis over different analysis experiments and
per analysis experiment will be created automatically
	a folder with experimental data files, e.g. ../experimental_data_yeast/
	a folder with databases specific for your experimental data (training subsets), e.g. ../databases_yeast/
	and a folder for analyses, e.g. ../analyses_yeast
	
The folder databases needs to be filled manually with the databases (Stitch, String if necessary and manually curated database for ID-Translations
'''
class PrepareData():
	#create general folder structure
	def create_folders(self,appendix):
		self.experimental="../experimental_data"+appendix
		self.databases="../databases"+appendix
		self.analyses="../analyses"+appendix
		databases="../databases"
		overall_analysis="../overall_analysis"
		if not os.path.isdir(self.analyses):
			os.makedirs(self.analyses)
			print("folder "+self.analyses+" created")
		if not os.path.isdir(self.databases):
			os.makedirs(self.databases)
			print("folder "+self.databases+" created")
		if not os.path.isdir(self.experimental):
			os.makedirs(self.experimental)
			print("folder "+self.experimental+" created")
		if not os.path.isdir(databases):
			os.makedirs(databases)
			print("folder "+databases+" created")
		if not os.path.isdir(overall_analysis):
			os.makedirs(overall_analysis)
			print("folder "+overall_analysis+" created")
		return
	
	#select an excel spread sheet via GUI
	def selectfile(self):
		print("Please select file with experimental data")
		root=Tk()
		root.withdraw()
		inputfile=askopenfilename(initialdir=self.experimental,filetypes=[("Excel Spread Sheets","*.xlsx")])
		return(inputfile)
	
	
	#prepare experimental data from excel file supplemental_datases (extracts only raw data)
	#organism={"S.cerevisiae","A.thaliana"}
	def prepare_expdata(self,organism,supplemental_datasets=None):
		if supplemental_datasets is None:
			supplemental_datasets=self.selectfile()
		excel_file=pd.ExcelFile(supplemental_datasets)
		if organism=="S.cerevisiae":
			ds1_writer=pd.ExcelWriter(self.experimental+"S.cerevisiae_dataset1.xlsx")
			ds2_writer=pd.ExcelWriter(self.experimental+"S.cerevisiae_dataset2.xlsx")
			ds3_writer=pd.ExcelWriter(self.experimental+"S.cerevisiae_dataset3.xlsx")
			ds4_writer=pd.ExcelWriter(self.experimental+"S.cerevisiae_dataset4.xlsx")
			ds1_met=excel_file.parse("DatasetS2A_dataset1")
			ds1_prot=excel_file.parse("DatasetS1A_dataset1")
			ds2_met=excel_file.parse("DatasetS2A_dataset2")
			ds2_prot=excel_file.parse("DatasetS1A_dataset2")
			ds3_met=excel_file.parse("DatasetS2A_dataset3")
			ds3_prot=excel_file.parse("DatasetS1A_dataset3")
			ds4_met=excel_file.parse("DatasetS2A_dataset4")
			ds4_prot=excel_file.parse("DatasetS1A_dataset4")
		elif organism=="A.thaliana":
			ds1_writer=pd.ExcelWriter(self.experimental+"A.thaliana_dataset1.xlsx")
			ds2_writer=pd.ExcelWriter(self.experimental+"A.thaliana_dataset2.xlsx")
			ds3_writer=pd.ExcelWriter(self.experimental+"A.thaliana_dataset3.xlsx")
			ds1_met=excel_file.parse("DatasetS2B_dataset1")
			ds1_prot=excel_file.parse("DatasetS1B_dataset1")
			ds2_met=excel_file.parse("DatasetS2B_dataset2")
			ds2_prot=excel_file.parse("DatasetS1B_dataset2")
			ds3_met=excel_file.parse("DatasetS2B_dataset3")
			ds3_prot=excel_file.parse("DatasetS1B_dataset3")
		ds1_met.to_excel(ds1_writer,"Metabolites_profiles_rep_raw",index=False)
		ds1_prot.to_excel(ds1_writer,"Protein_profiles_rep_raw",index=False)
		ds2_met.to_excel(ds2_writer,"Metabolites_profiles_rep_raw",index=False)
		ds2_prot.to_excel(ds2_writer,"Protein_profiles_rep_raw",index=False)
		ds3_met.to_excel(ds3_writer,"Metabolites_profiles_rep_raw",index=False)
		ds3_prot.to_excel(ds3_writer,"Protein_profiles_rep_raw",index=False)
		ds1_writer.save()
		ds2_writer.save()
		ds3_writer.save()
		ds1_writer.close()
		ds2_writer.close()
		ds3_writer.close()
		if organism=="S.cerevisiae":
			ds4_met.to_excel(ds4_writer,"Metabolites_profiles_rep_raw",index=False)
			ds4_prot.to_excel(ds4_writer,"Protein_profiles_rep_raw",index=False)
			ds4_writer.save()
			ds4_writer.close()
		return
	
	
	#prepare databases
	def prepare_databases(self,supplemental_datasets=None):
		if supplemental_datasets is None:
			supplemental_datasets=self.selectfile()
		excel_file=pd.ExcelFile(supplemental_datasets)
		ID_translations=excel_file.parse("DatasetS10")
		ID_translations.to_csv("../databases/ID_translations.tsv",sep="\t")
		print("please put manually the following files into the folder 'databases':")
		#print("MetaboliteLibrary_ChemForm_CID.txt with manually curated Stitch ID translations")
		print("Sac_4932.protein_chemical.links.v5.0.tsv.gz (Stitch database for S.cerevisiae)")
		print("Ara_3702.protein_chemical.links.v5.0.tsv.gz (Stitch database for A.thaliana)")
		print("metabolite_kegg_classes_pathways.tsv (manually created database with kegg classes and pathways for metabolites)")
		print("COG.mappings.v10.5.txt from the STRING database")
		return

###############################################################
###############################################################
###############################################################


#Exception for try-except expressions if overwriting is requested by user
class OverwriteException(Exception):
	def __init__(self,message):
		super().__init__(message)


'''
in inheriting file

from OverwriteException import OverwriteException

try:
	if overwrite==True:
		raise OverwriteException("overwriting of existing file requested")
	load("try") # not reached if OverwriteException
except (FileNotFoundError,OverwriteException):
	pass

'''

###############################################################
###############################################################
###############################################################

#class providing methods to translate protein and metabolite IDs
class IDTranslations():
	
	#get files with experimental data
	def get_expfiles(self):
		expfiles=glob.glob(self.experimental+"*.xlsx") # if you want to predict data different from training data, adjust this path to prediction data 
		if expfiles==[]:
			root=Tk()
			root.withdraw()
			expfiles=glob.glob(askdirectory(initialdir="../")+"/*.xlsx") 
		expfiles.sort()#sort them in alphabetical order
		return(expfiles)
		
		
	#transforms set (of proteins or metabolites) to one string
	def ids_to_string(self,ids):
		str_ids=""
		for id in ids:
			str_ids=str_ids+" "+str(id)
		return(str_ids)
		
		
	#translates protein IDs of query and prints UniProt ID, STRING ID and gene ID to file, in_id_type = string of input id type, out_id_types = list of strings for output ids
	#see https://www.uniprot.org/help/api_idmapping
	def translate_proteins(self,ids,in_id_type,out_id_types=["","ID","STRING_ID","ENSEMBLGENOME_PRO_ID","GENENAME"],save=False): 
		#out_id_type as list
		#UniProtID = "ID" for input, "" for output, e.g. P50077
		#STRING ID = "STRING_ID", e.g. 3702.AT5G63570.1
		#Ensembl Genome Protein = "ENSEMBLGENOME_PRO_ID", e.g. AT5G63570.1
		#gene name = "GENENAME" e.g. CCH1
		
		df=pd.DataFrame(index=range(len(ids)),columns=out_id_types)
		df[in_id_type]=ids
		df=df.set_index(in_id_type,drop=False)
		url="https://www.uniprot.org/uploadlists/"
		str_ids=self.ids_to_string(ids)
		for out_id_type in out_id_types:
			if out_id_type!=in_id_type:
				params={"from": in_id_type,"to": out_id_type,"format":"tab","query":str_ids}
				while True: ###
					try:
						resp=requests.post(url,params)
						break
					except ConnectionError:
						print("ConnectionError to uniprot.org") 
						time.sleep(5)
				out=resp.text
				out=out[8:].replace("\n","\t").split("\t")
				out.remove("")
				ids_found=out[::2] #every even element corresponds to input id
				out=out[1:][::2] #every odd element corresponds to output id
				for i in range(len(ids_found)):
					df.loc[ids_found[i],out_id_type]=out[i]
		if save!=False:
			df.to_csv(save,header=True,index=False,sep="\t")
		return(df)
		
	
	
	#padd cids with zeros in front so that they match structure from database (8 digits)
	def padd_cids(self,cids):
		cids_padded=list()
		for cid in cids:
			if type(cid)!=str: #check whether given cid is integer or string
				cid=str(cid)
			if cid[0]!="0": #check whether cid is already padded
				cid=str(int(eval(cid))) # convert cids like 6083.0 to 6083
			if len(cid)<8:
				cid_padded="0"*(8-len(cid))+cid
				cids_padded.append(cid_padded)
			else:
				cids_padded.append(cid)
		return(cids_padded)
			
	
	#translates either given list/set of metabolites or metabolites from file (self.metabolites) (translates only distinct metabolites, if translation of all wanted, modify at "continue")
	def translate_metabolites_pubchem(self,metabolites=None):
		if metabolites==None:
			metabolites=self.get_metabolites()
		t1=0
		notfound=0
		not_translated_metabolites=""
		df=pd.DataFrame(columns=["CIDs"])
		for name in metabolites:
			names=name.split("/")
			if "" in names: #sometimes at the end of metabolite name there is a "/" in file
				names.remove("")
			if len(names)>1:
				continue #if more than 1 metabolite given (MS data not distinct), no need for translation
			else:
				metabolite=names[0].split(" (")[0]
				t2=time.time()
				if t2-t1<=60/400 and t1!=0: #on PubChem, no more than 400 requests per minute
					time.sleep(60/400-(t2-t1))
				t1=time.time()
				query="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+metabolite+"/cids/TXT"
				resp=requests.post(query)
				out=resp.text
				if out[:11]=="Status: 400":
					print(metabolite)
				if out[:11]!="Status: 404" and out[:11]!="Status: 400": #...if compound found
					cids=list(out.split("\n"))
					cids.remove("")
					#padd cids with zeros in front so that they match structure from database (8 digits)
					cids_padded=self.padd_cids(cids)
					df.loc[name,"CIDs"]=cids_padded
				else:
					notfound=notfound+1
					not_translated_metabolites=not_translated_metabolites+metabolite+"; "
		'''
		print(str(notfound)+" metabolites out of "+str(len(metabolites))+" not translated")
		if notfound>0:
			print("Not translated metabolites: "+not_translated_metabolites)
		'''
		#df.to_csv(self.databases+self.organism+"_metabolites.tsv",header=True,index=True,sep="\t")
		return(df)
		
		
	#translate PubChem CIDs online on Pubchem
	def translate_cids_pubchem(self,cids):
		t1=0
		notfound=0
		not_translated_cids=""
		df=pd.DataFrame(columns=["Name"])
		for cid in cids:
			t2=time.time()
			if t2-t1<=60/400 and t1!=0: #on PubChem, no more than 400 requests per minute
				time.sleep(60/400-(t2-t1))
			t1=time.time()
			query="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(cid)+"/synonyms/TXT"
			resp=requests.post(query)
			out=resp.text
			if out[:11]!="Status: 404": #...if compound found
				names=list(out.split("\n"))
				df.loc[cid,"Name"]=names[0]
			else:
				notfound=notfound+1
				not_translated_cids=not_translated_cids+"; "+str(cid)
		
		print(str(notfound)+" cids out of "+str(len(cids))+" not translated")
		if notfound>0:
			print("Not translated metabolites: "+not_translated_cids)
		
		#df.to_csv(self.databases+self.organism+"_metabolites.tsv",header=True,index=True,sep="\t")
		return(df)
		
		
		
		
	#translates given list of metabolites
	def translate_metabolites_offline(self,metabolites=None,online_translations=False):
		try:
			if metabolites==None: #if no metabolites given as argument, try to load them with inbuild method, if possible
				metabolites=self.get_metabolites()
		except:
			print("no metabolites given for translation")
		df=pd.read_csv("../databases/ID_translations.tsv",sep="\t").set_index("metabolite_name")
		df=df.rename(columns={"Stitch_CID":"CID"})
		return(df.loc[metabolites,"CID"])
	
	
	'''
	#with extensive manually curated textfile with ID translations
	def translate_metabolites_offline(self,metabolites=None,online_translations=True):
		try:
			if metabolites==None: #if no metabolites given as argument, try to load them with inbuild method, if possible
				metabolites=self.get_metabolites()
		except:
			print("no metabolites given for translation")
		df=pd.DataFrame(columns=["CID"])
		try:
			db=open("../databases/MetaboliteLibrary_ChemForm_CID.txt","rt")
		except FileNotFoundError:
			print("Please select Marcin's database with metabolite name - CID - mapping")
			root=Tk()
			root.withdraw()
			db=askopenfilename(initialdir=self.databases,filetypes=[(".txt","*.txt")])
		for line in db:
			coloumns=line.split("\t")
			if coloumns[0] in metabolites:
				df.loc[coloumns[0],"CID"]=coloumns[2][4:-1]
		
		not_translated_metabolites=set(metabolites)-set(df.index.tolist()) #difference between metabolites and df.index
		#print(str(len(metabolites)-len(df))+" metabolites out of "+str(len(metabolites))+" not translated")
		if len(not_translated_metabolites)>0:
			print("Not translated metabolites: "+";".join(str(m) for m in list(not_translated_metabolites)))
			if online_translations:
				print("trying online translation on pubchem")
				#better to translate manually and add them to offline library
				pubchem_translation=self.translate_metabolites_pubchem(not_translated_metabolites)
				for i in pubchem_translation.index:
					df.loc[i,"CID"]=pubchem_translation.loc[i,"CIDs"][0] #take first CID found on pubchem
				
				print("Finally "+str(len(metabolites)-len(df))+" metabolites out of "+str(len(metabolites))+" not translated")
				not_translated_metabolites=set(metabolites)-set(df.index.tolist()) #difference between metabolites and df.index
				if len(not_translated_metabolites)>0:
					print("Not translated metabolites: "+";".join(str(m) for m in list(not_translated_metabolites)))
		return(df)
		'''
	
	
	#translates list or set of metabolite CID to metabolite name occuring in excel sheet with Marcin's library
	def translate_CIDs_offline(self,cids):
		#metabolites=cids
		df=pd.read_csv("../databases/ID_translations.tsv",sep="\t",dtype={"Stitch_CID":str,"PubChem_CID":str}).set_index("Stitch_CID")
		df=df.rename(columns={"metabolite_name":"Name"})
		return(df.loc[list(map(lambda x: str(x),cids))])
		
		
	
	#translate CIDs from Stitch to PubChem CIDs
	def translate_CIDs_from_Stitch_to_PubChem(self,cids):
		#translate CIDs from Stitch to PubChem CIDs
		metabolite_names=self.translate_CIDs_offline(cids)
		metabolite_cids_pubchem=self.translate_metabolites_pubchem(metabolite_names["Name"].tolist())
		#metabolite_pubchem_high_corr=list(map(lambda x: int(x),list(itertools.chain.from_iterable(metabolite_cids_pubchem["CIDs"].tolist())))) #including all CIDs found per name
		#metabolite_pubchem_high_corr=list(map(lambda x: metabolite_cids_pubchem.loc[x,"CIDs"][0],metabolite_cids_pubchem.index)) #extract always first CID
		metabolite_names.reset_index(inplace=True)
		metabolite_names=metabolite_names.rename(columns={"index": "Stitch_CID"})
		metabolite_names["PubChem_CID"]=1
		for i in metabolite_names.index:
			metabolite_names.loc[i,"Stitch_CID"]=int(metabolite_names.loc[i,"Stitch_CID"])
			if metabolite_names.loc[i,"Stitch_CID"]==92823: #2'3' Cyclic Guanosine monophosphate
				metabolite_names.loc[i,"PubChem_CID"]= 135398728
			elif metabolite_names.loc[i,"Stitch_CID"]==6802: #Guanosine
				metabolite_names.loc[i,"PubChem_CID"]= 135398635
			elif metabolite_names.loc[i,"Stitch_CID"]==101812: #Adenosine-2,3-cyclic monophosphate
				metabolite_names.loc[i,"PubChem_CID"]= 23666344  
			elif metabolite_names.loc[i,"Stitch_CID"]==6804: #Guanosine 5'-Monophosphate
				metabolite_names.loc[i,"PubChem_CID"]= 135398631 
			elif metabolite_names.loc[i,"Stitch_CID"]==7009629: #Leu-Glu
				metabolite_names.loc[i,"PubChem_CID"]= 5259590
			elif metabolite_names.loc[i,"Stitch_CID"]==6992295: #Leu-Ala
				metabolite_names.loc[i,"PubChem_CID"]= 81721
			elif metabolite_names.loc[i,"Stitch_CID"]==92843: #Leu-Gly 
				metabolite_names.loc[i,"PubChem_CID"]= 97364
			elif metabolite_names.loc[i,"Stitch_CID"]==435718: #Leu-Ile
				metabolite_names.loc[i,"PubChem_CID"]= 7010534
			elif metabolite_names.loc[i,"Stitch_CID"]==4682588: #Lys-Leu
				metabolite_names.loc[i,"PubChem_CID"]= 7016103
			elif metabolite_names.loc[i,"Stitch_CID"]==3425788: #Met-Leu
				metabolite_names.loc[i,"PubChem_CID"]= 444206
			elif metabolite_names.loc[i,"Stitch_CID"]==4422358: #Phe-Glu 
				metabolite_names.loc[i,"PubChem_CID"]= 151134
			elif metabolite_names.loc[i,"Stitch_CID"]==4078229: #Phe-Leu 
				metabolite_names.loc[i,"PubChem_CID"]= 76808
			elif metabolite_names.loc[i,"Stitch_CID"]==5246010: #Ile-Ala
				metabolite_names.loc[i,"PubChem_CID"]= 7009577
			elif metabolite_names.loc[i,"Stitch_CID"]==13817313: #Ile-Leu
				metabolite_names.loc[i,"PubChem_CID"]= 7019083
			elif metabolite_names.loc[i,"Stitch_CID"]==435728: #Ile-Phe
				metabolite_names.loc[i,"PubChem_CID"]= 7009596
			elif metabolite_names.loc[i,"Stitch_CID"]==6021: #Inosine
				metabolite_names.loc[i,"PubChem_CID"]= 135398641
			elif metabolite_names.loc[i,"Stitch_CID"]==18218240: #Ser-Ile
				metabolite_names.loc[i,"PubChem_CID"]= 71429009
			elif metabolite_names.loc[i,"Stitch_CID"]==13919048: #Ser-Leu
				metabolite_names.loc[i,"PubChem_CID"]= 7015695
			elif metabolite_names.loc[i,"Stitch_CID"]==416721: #Thr-Val
				metabolite_names.loc[i,"PubChem_CID"]= 7020902
			elif metabolite_names.loc[i,"Stitch_CID"]==6992386: #Gly-Ile
				metabolite_names.loc[i,"PubChem_CID"]= 88079
			elif metabolite_names.loc[i,"Stitch_CID"]==3459978: #Lys-Phe
				metabolite_names.loc[i,"PubChem_CID"]= 151410
			elif metabolite_names.loc[i,"Stitch_CID"]==4441256: #Arg-Phe
				metabolite_names.loc[i,"PubChem_CID"]= 150964
			elif metabolite_names.loc[i,"Stitch_CID"]==417358: #Ala-Ile
				metabolite_names.loc[i,"PubChem_CID"]= 7408079
			elif metabolite_names.loc[i,"Stitch_CID"]==5886: #NADP
				metabolite_names.loc[i,"PubChem_CID"]= 5885
			else:
				try:
					metabolite_names.loc[i,"PubChem_CID"]=int(metabolite_cids_pubchem[metabolite_cids_pubchem.index==metabolite_names.loc[i,"Name"]].iloc[0,0][0])
				except IndexError: #empty Dataframe
					metabolite_names.loc[i,"PubChem_CID"]=metabolite_names.loc[i,"Stitch_CID"]
		return(metabolite_names)
		


###############################################################
###############################################################
###############################################################

#class to handle database files
#an instance of this class is linked to the database, further attributes are metabolite and protein names
#Terminology: 	raw database: 			Stitch database for Protein-Metabolite-interactions with score for each interaction
#				compressed database:	raw database shortened to interactions above a given score (score_cutoff)


class DBHandler(IDTranslations):
	#orgstring corresponds to filename including path, e.g. "organism.xlsx", "../experimental_data/organism.xlsx"
	#self.file is compressed database
	def __init__(self,orgstring,score_cutoff=None,overwrite=False,simulation=False,experimental="../experimental_data/",databases="../databases/",analyses="../analyses/"):
		self.experimental=experimental
		self.databases=databases
		if score_cutoff!=None:
			self.score_cutoff=int(score_cutoff)
		else:
			self.score_cutoff=800 
		self.overwrite=overwrite #overwrite existing compressed database
		self.simulation=simulation
		#self.organisms=list(filter(lambda w: w not in ["[",",","]",""], re.split(",",orgstring))) #excluding parenthesis and comma from input "list" orgstring
		self.organism=orgstring[:-5].split("/")[-1] #self.organism only to filename
		if self.simulation==False:
			self.file="../databases/"+self.organism.split("_")[0]+"_"+str(self.score_cutoff)+"_interactions.txt"
		else:
			self.file="../simulation_data/"+self.organism.split("_")[0]+"_"+str(self.score_cutoff)+"_interactions.txt"
		self.extract_interactions()
		self.metabolites=set() # as CID
		self.proteins=set() #as String ID
		self.extract_meta_and_prots()
		try:
			self.proteins=pd.read_csv("../databases/"+self.organism.split("_")[0]+"_ID_collection.tsv",sep="\t") ### new: .split("_")[0], to make only one per organism
			print("database with different IDs for "+self.organism+" found and loaded")
		except FileNotFoundError: 
			print("loading and translating protein IDs for "+self.organism)
			self.proteins=set(self.translate_proteins(self.proteins,in_id_type="STRING_ID",out_id_types=[""],save=False)[""].values.tolist())
			self.proteins=self.translate_proteins(self.proteins,in_id_type="ID",out_id_types=["ID","STRING_ID","ENSEMBLGENOME_PRO_ID","GENENAME"],save="../databases/"+self.organism.split("_")[0]+"_ID_collection.tsv")### new: .split("_")[0], to make only one per organism
			print("database with different IDs for "+self.organism+" created and loaded")
	
	
	
	#load raw Stitch database file for one organism
	def load_raw_db(self):
		#if self.organism in ["A.thaliana", "S.cerevisiae", "HeLa", "E.coli"] or "A.thaliana" in self.organism:
		if self.simulation==False:
			if "A.thaliana" in self.organism or "S.cerevisiae" in self.organism or "HeLa" in self.organism or "E.coli" in self.organism: 
				try:
					if "A.thaliana" in self.organism:
						print("loading database for A.thaliana")
						db=gzip.open("../databases/Ara_3702.protein_chemical.links.v5.0.tsv.gz","rt")
					if "S.cerevisiae" in self.organism:
						db=gzip.open("../databases/Sac_4932.protein_chemical.links.v5.0.tsv.gz","rt")
					if "E.coli" in self.organism:
						db=gzip.open("../databases/Ecoli_511145.protein_chemical.links.v5.0.tsv.gz","rt")
				except FileNotFoundError:
					print("raw Stitch Database for "+self.organism+" not found. Please select manually")
					root=Tk()
					root.withdraw()
					db=askopenfilename(initialdir=self.databases,filetypes=[(".gz","*.gz")])
			else:
				print("no raw database downloaded for " +self.organism)
				print("check for correct spelling")
				print("available DBs: A.thaliana, S.cerevisiae, HeLa, E.coli") #HeLa not yet available
				print("terminated in function load_raw_db")
				#sys.exit(0)
				return
		else:
			print("loading simulation database for "+self.organism)
			if "A.thaliana" in self.organism or "S.cerevisiae" in self.organism or "HeLa" in self.organism or "E.coli" in self.organism: 
				if "A.thaliana" in self.organism:
					db=open("../simulation_data/Ara_3702_testdatabase.txt","rt")
				if "S.cerevisiae" in self.organism:
					db=open("../simulation_data/Sac_4932_testdatabase.txt","rt")
		print("initial 'raw' database for "+self.organism+" loaded")
		return(db)
					
	
	#compressing database for one organism to interactions with score higher score_cutoff and saving reduced db to outputfile
	def compress_db(self,db):
		out=open(self.file,"w")
		out.write("chemical\tprotein\tcombined_score\n")
		print("Compressing database to empirically significant interactions for "+self.organism+" with score_cutoff= "+str(self.score_cutoff)+"...")
		for line in db:
			score=re.search("[0-9]{1,4}\n$",line)
			if score:
				if int(score.group(0))>=self.score_cutoff:
					out.write(line)
		print("database for "+self.organism+" with score_cutoff= "+str(self.score_cutoff)+" compressed for empirical significant interactions")
		
		
	#extract metabolites and proteins
	def extract_meta_and_prots(self):
		print("reading metabolites and proteins from database...")
		out=open(self.file,"rt")
		for line in out:
			coloumns=line.split()
			self.metabolites.add(coloumns[0])
			#protein=re.search("\.(.+)$",coloumns[1]) ensembl genome protein ID =>self.proteins.add(protein.group(0)[1:])
			self.proteins.add(coloumns[1]) #STRING ID
		self.metabolites.remove("chemical")
		
		
	#create (if necessary) databases for only significant (according to score) interactions
	def extract_interactions(self):
		if not os.path.isfile(self.file) or self.overwrite==True:
			db=self.load_raw_db() #load stitch database
			self.compress_db(db) #select only significant interactions
		else:
			print("compressed database "+self.organism+" with score cutoff = "+str(self.score_cutoff)+" found")
		
		
	#returns metabolites for a given organism
	def get_metabolites(self):
		return(self.metabolites)
	
	
	#returns protein UniProt IDs for a given organism
	def get_proteins(self):
		return(set(self.proteins.ID.values.tolist()))
	
	
	#opens the COG database where clusters of orthologous group for each protein could be find
	def open_cog_db(self):
		try:
			cog_db=open("../databases/COG.mappings.v10.5.txt","rt")
		except FileNotFoundError:
			print("Database for orthologous groups (COG) for "+self.organism+" not found. Please select manually")
			root=Tk()
			root.withdraw()
			cog_db=askopenfilename(initialdir=self.databases,filetypes=[(".gz","*.gz"),("text files","*.txt")])
		return(cog_db)
		
		
	
	
	#find orthologous group from String database, returns result as dataframe (including Uniprot IDs)
	def get_orthologs(self, string_ids=None,save=False):
		try:
			df=pd.read_csv("../databases/"+self.organism.split("_")[0]+"_string_uniprot_cog.tsv",sep="\t") ### new: .split("_")[0], to make only one per organism
			df=df.set_index("STRING_ID",drop=True)
			if string_ids is None:
				return(df)
			else:
				return(df.loc[string_ids])
			'''
			if string_ids is not None:
				print("check and correct code, DBHandler function get_orthologs called with input string_ids")
				
			else:
				try:
					df=pd.read_csv("../databases/"+self.organism+"_string_uniprot_cog.tsv",sep="\t")
					return(df)
				except:
					pass
			'''
		except:
			print("getting string orthologs for "+self.organism)
			if string_ids==None:
				string_ids=self.proteins.STRING_ID.values.tolist()
				save=True
			df=pd.DataFrame(index=range(len(string_ids)),columns=["STRING_ID","COG","CID"])
			df["STRING_ID"]=string_ids
			df=df.set_index("STRING_ID",drop=True)
			cog_db=self.open_cog_db()
			match=0 
			for cog_line in cog_db:
				cog_coloumns=cog_line.split()
				for i in range(len(string_ids)):
					if string_ids[i]==cog_coloumns[0]:
						match=match+1 
						#print(string_ids[i]+"\t"+cog_coloumns[3]) ###
						try:
							if pd.isna(df.loc[string_ids[i],"COG"]):###
								df.loc[string_ids[i],"COG"]=cog_coloumns[3]
							else:
								df.loc[string_ids[i],"COG"]=df.loc[string_ids[i],"COG"]+";"+cog_coloumns[3]
						except ValueError: #if one string id occurs several times (some uniprotIDs point to different UniProt IDs)
							if pd.isna(df.loc[string_ids[i],"COG"].iloc[0]):
								df.loc[string_ids[i],"COG"]=[cog_coloumns[3]]*len(df.loc[string_ids[i]])
							else:
								df.loc[string_ids[i],"COG"]=list(map(lambda x: x+";"+cog_coloumns[3],df.loc[string_ids[i],"COG"]))
			print(str(match)+" annotated clusters of orthologous group for given protein ids found")
			if save==True:
				uniprots=self.translate_proteins(df.index.tolist(),in_id_type="STRING_ID",out_id_types=[""])
				df["UniProt_ID"]=list(map(lambda s: uniprots.loc[s,""],df.index.tolist()))
				df.drop("CID",axis=1).to_csv("../databases/"+self.organism.split("_")[0]+"_string_uniprot_cog.tsv",header=True,index=True,sep="\t") ###
			return(df)
		
		
	
	#gets dataframe with string IDs and COG and initialized Chemical IDs (CID) and adds corresponding CIDs from database
	def find_metabolites(self,df_string_cogs):
		out=open(self.file,"rt")
		for line in out:
			coloumns=line.split()
			if coloumns[1] in df_string_cogs.index:
				if pd.isna(df_string_cogs.loc[coloumns[1],"CID"]):
					df_string_cogs.loc[coloumns[1],"CID"]=coloumns[0]
				else:
					df_string_cogs.loc[coloumns[1],"CID"]=df_string_cogs.loc[coloumns[1],"CID"]+";"+coloumns[0]
		if np.sum(pd.isna(df_string_cogs.CID))!=0:
			print(str(np.sum(pd.isna(df_string_cogs.CID)))+"proteins without metabolites annotated") # to check: if not 0, no corresponding metabolite found
		return(df_string_cogs)
	
	
	
	
	#trim stitch database to intersection with YMDB to get rid of non-native metabolites, e.g. drugs
	def trim_db_to_YMDB_containing_metabolites(self):
		try:
			df=pd.read_csv("../databases/YMDB_ID_translation.tsv",sep="\t")
		except FileNotFoundError:
			#initiate df with id translations
			df=pd.DataFrame(columns=["YMDB_ID","InChI","CID"])
			#pubchem queries
			query="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON"
			header={"Content-Type": "application/x-www-form-urlencoded"}
			#scan through YMDB and extract YMDB ID and InChI
			ymdb_file="../databases/ymdb.sdf"
			ymdb=open(ymdb_file,"rt")
			for line in ymdb:
				if line=="> <DATABASE_ID>\n":
					index=next(ymdb)[:-1]
				elif line=="> <INCHI_IDENTIFIER>\n":
					df.loc[index,"InChI"]=next(ymdb)[:-1]
				elif line=="> <YMDB_ID>\n":
					df.loc[index,"YMDB_ID"]=next(ymdb)[:-1]
			#translate InChI to PubChem CIDs
			for i in df.index:
				t1=0
				if pd.isna(df.loc[i,"CID"]):
					inchi=df.loc[i,"InChI"]
					#inchi="InChI=1S/C5H10N2O3/c6-3(5(9)10)1-2-4(7)8/h3H,1-2,6H2,(H2,7,8)(H,9,10)/t3-/m0/s1"
					t2=time.time()
					if t2-t1<=60/400 and t1!=0: #on PubChem, no more than 400 requests per minute
						time.sleep(60/400-(t2-t1))
					t1=time.time()
					resp=requests.post(query,data={"inchi":inchi},headers=header)
					if resp.status_code==200:
						out=resp.json()["IdentifierList"]["CID"]
						if len(out)==1: #one match
							df.loc[i,"CID"]=out[0]
					else:
						print("error with ID: "+i)
						print("pubchem error code: "+str(resp.status_code))
			#save id translations
			df.dropna().to_csv("../databases/YMDB_ID_translation.tsv",sep="\t",header=True,index=False)
		#load trimmed database and save it with intersecting metabolites
		stitch_db=pd.read_csv(self.file,sep="\t")
		ymdb_metas=list(map(lambda x: "CIDs"+x,self.padd_cids(df["CID"].tolist())))
		for i in stitch_db.index:
			if stitch_db.loc[i,"chemical"] not in ymdb_metas:
				stitch_db.drop(i,axis=0,inplace=True)
		stitch_db.to_csv("../databases/"+self.organism.split("_")[0]+"_"+str(self.score_cutoff)+"_native_interactions.txt",sep="\t",header=True,index=False)
		return
	
	
	#reads and translates STITCH database with selected score threshold and translates IDs to save it like ML predictions
	#under construction
	def make_STITCH_predictions_file(self):
		db=pd.read_csv("../databases/Sac_4932.protein_chemical.links.v5.0.tsv",sep="\t")
		uniprot_df=self.translate_proteins(db["protein"].tolist(),in_id_type="STRING_ID",out_id_types=["","ID"],save=False)
		db["protein"]
		db["metabolite"]=list(map(lambda m: int(m[4:]),db["chemical"]))
		pass
	
	
	def annotate_STITCH_predictions(self):
		pass
	
	
	
	def helpfun():
		print("====Usage:====")
		print("first argument: organism for corresponding database")
		print("second optional argument: interaction score cutoff")
		print("second/third optional argument: set True if you want to overwrite existing compressed database (e.g. with new score_cutoff)")
		print("e.g.: DBHandler.py [A.thaliana,S.cerevisiae] 500 True")
	
	
	




###############################################################
###############################################################
###############################################################


# class to handle experimental data, e.g. loading excel sheets, normalizing and pooling profiles, etc
# an instance of this class is an excel sheet, attributes are file, metabolite and protein names
# filenames start with organism, e.g S.cerevisiae_cell_cultures.xlsx
class FileHandler(IDTranslations):
	#'experimental' corresponds to folder including files with experimental data, 'databases' corresponds to folder containing databases from given organism
	def __init__(self,inputfile=None,experimental="../experimental_data/",databases="../databases/"):
		self.experimental=experimental
		self.databases=databases
		if inputfile==None:
			self.inputfile=self.selectfile()
		else:
			self.inputfile=inputfile
		self.file=self.load_file()
		self.organism=self.inputfile.split("/")[-1][:-5] #filename has to be organism_specification.xlsx, e.g.: "A.thaliana_rosettes.xlsx"
		#self.load_proteins() #UniProtIDs
		#self.load_metabolites() #names of metabolites
		#print(str(len(self.proteins))+" proteins, "+str(len(self.metabolites))+" metabolites in file "+self.organism)
	
	#select an excel spread sheet via GUI
	def selectfile(self):
		print("Please select file with experimental data")
		root=Tk()
		root.withdraw()
		inputfile=askopenfilename(initialdir=self.experimental,filetypes=[("Excel Spread Sheets","*.xlsx")])
		return(inputfile)
		
		
	#load excel file
	def load_file(self):
		print("reading excel file "+self.inputfile)
		return(pd.ExcelFile(self.inputfile))
		
		
	#returns available sheets in file
	def get_sheets(self):
		return(self.file.sheet_names)
		
		
	#load a sheet from excel file
	def load_sheet(self,sheet,**kwargs):
		print("loading excel sheet "+sheet+" for "+self.organism)
		return(self.file.parse(sheet,**kwargs))
		
		
	#loads protein UniProtIDs from sheet and stores them in self.proteins
	def load_proteins(self):
		proteinlist=list()
		protsheet=self.load_sheet("Protein_profiles_rep_raw")
		proteins=protsheet["UniProt_IDs"][1:].dropna().values.tolist() 
		nans=protsheet[protsheet["UniProt_IDs"].isna().tolist()]
		if len(nans)>0:
			print(str(len(nans))+" Proteins with missing UniProt ID in excel sheet")
			print("Following proteins have no UniProtIDs:")
			print(nans)
		for i in proteins:
			for j in re.split(";|/",i): #split i by multiple delimiters ";" and "/"
				proteinlist.append(j)
		self.proteins=set(proteinlist)
		print(str(len(self.proteins))+" proteins loaded for "+self.organism)
		
	
	
	#loads metabolite names from sheet and stores them in self.metabolites
	def load_metabolites(self):
		metsheet=self.load_sheet("Metabolites_profiles_rep_raw",dtype={"Stitch_CID":str})
		self.metabolites=set(metsheet["Stitch_CID"][1:].values.tolist())
		if np.nan in self.metabolites:
			self.metabolites.remove(np.nan)
		print(str(len(self.metabolites))+" metabolites loaded for "+self.organism)
		
	
	#returns list of metabolite names
	def get_metabolites(self):
		if not hasattr(self,"metabolites"):
			self.load_metabolites()
		return(self.metabolites)
		
		
	#returns list of protein UniProt IDs
	def get_proteins(self):
		if not hasattr(self,"proteins"):
			self.load_proteins()
		return(self.proteins)
		
	
	
	#extract profiles per repetition for proteins
	def load_prot_profiles_and_split_into_reps(self,profiles="raw"):
		#select sheets
		if profiles=="raw":
			prot_sheetname="Protein_profiles_rep_raw"
		elif profiles=="deconvoluted":
			prot_sheetname="Protein_profiles_rep_deconv"
		else:
			print("no "+profiles+" available")
			print("available profiles: raw or deconvoluted")
			sys.exit(0)
		#PROTEINS
		#load sheets and set indices
		prot_id_coloumn="UniProt_IDs"
		prot_sheet=self.load_sheet(prot_sheetname)
		#prot_sheet=prot_sheet[1:].dropna().set_index(prot_id_coloumn,drop=True) #set coloumn with ids as index and drop NaN
		#prot_sheet=prot_sheet[1:].set_index(prot_id_coloumn,drop=True)
		#drop first line with headers and sizes of fractions
		prot_sheet=prot_sheet.drop(prot_sheet[prot_sheet[prot_id_coloumn]==prot_id_coloumn].index)
		prot_sheet=prot_sheet.drop(prot_sheet[prot_sheet[prot_id_coloumn]==" "+prot_id_coloumn].index)
		prot_sheet=prot_sheet.set_index(prot_id_coloumn,drop=True)
		#split coloumnnames into different repetitions (=Experiments in Excel sheet)
		#l=list(filter(lambda col: "Experiment" in col,prot_sheet.columns.tolist()))
		#l=list() #coloumns with profile data
		prot_reps=list() #length of this list shows number of experiments
		for prot_col in prot_sheet.columns.tolist():
			if "Experiment" in prot_col:
				#l.append(prot_col)
				if prot_col[:len("Experiment")+1] not in prot_reps: #adjust if more than 9 repetitions done
					prot_reps.append(prot_col[:len("Experiment")+1])
		
		#split data into dataframes for each repetition
		prot_repetitions=list() #list where every element corresponds to one dataframe for one repetition
		for prot_exp in prot_reps:
			double_ids=list() #ids which occur twice or more and are therefore left out for further analysis :(
			nan_uids=0 #missing uniprot ids
			prot_exp_cols=prot_sheet.columns[prot_sheet.columns.str.contains(prot_exp)]
			prot_sheet_ids_separated=pd.DataFrame(columns=prot_exp_cols) #in excel sheet multiple UniProt_IDs per profile, these profiles are doublicated, so that every UniProt_ID has a profile
			for uids in prot_sheet.index.tolist(): #skip column title
				if not pd.isna(uids):
					for uid in uids.split("/"):
						try:
							prot_sheet_ids_separated=prot_sheet_ids_separated.append(prot_sheet.loc[uids,prot_exp_cols].rename(uid),ignore_index=False)###
						except TypeError: #if uid occurs more than once, it is left out from further analaysis :(
							double_ids.append(uid)
				else:
					nan_uids+=1
			prot_repetitions.append(prot_sheet_ids_separated)
		
		print("number of multiple times occuring uniprot IDs: "+str(len(double_ids))) #same for every repetition
		print("number of missing uniprot IDs: "+str(nan_uids))
		return([prot_sheet,prot_repetitions])
		
		
	#extract profiles per repetition for metabolites
	def load_met_profiles_and_split_into_reps(self,profiles="raw"):
		#METABOLITES
		#select sheets
		if profiles=="raw":
			met_sheetname="Metabolites_profiles_rep_raw"
		elif profiles=="deconvoluted":
			met_sheetname="Metabolites_profiles_rep_deconv"
		else:
			print("no "+profiles+" available")
			print("available profiles: raw or deconvoluted")
			sys.exit(0)
		#load sheet and set index
		met_id_coloumn="Stitch_CID"
		met_sheet=self.load_sheet(met_sheetname,dtype={"Stitch_CID":str})
		#drop first line with headers and sizes of fractions
		#met_sheet=met_sheet.drop(met_sheet[met_sheet["Metabolite Names"]=="Metabolite Names"].index)
		#met_sheet=met_sheet.drop(met_sheet[met_sheet["Metabolite Names"]==" Metabolite Names"].index)
		met_sheet=met_sheet.drop(met_sheet[met_sheet["Stitch_CID"]=="Stitch_CID"].index)
		#met_sheet=met_sheet.dropna() ### if dropped and one chromatogram is missing, everything will be missing
		met_sheet=met_sheet.set_index("Stitch_CID",drop=True)
		#met_sheet=met_sheet.drop(np.nan)
		#split coloumnnames into different repetitions (=Experiments in Excel sheet)
		#l=list() #coloumns with profile data
		met_reps=list()
		for met_col in met_sheet.columns.tolist():
			if "Experiment" in met_col:
				#l.append(met_col)
				if met_col[:len("Experiment")+1] not in met_reps: #adjust if more than 9 repetitions done
					met_reps.append(met_col[:len("Experiment")+1])
		met_repetitions=list()
		for met_exp in met_reps:
			met_repetitions.append(met_sheet.iloc[:,met_sheet.columns.str.contains(met_exp)])
		return([met_sheet,met_repetitions])
	
	
	
	#get protein profiles pooled over repetitions with method and pool protein profiles belonging to one OG, as given in xset, normalize every profile
	#xset_profiles is dataframe with protein-metabolite pairs to extract
	#profiles={"raw","deconvoluted"}
	#normalized can be "sum" to normalize profiles to an area of one or "max" to normalize to maximum peak intensity
	def extract_prot_profiles(self,xset,xset_profiles=None,method=np.mean,profiles="raw",normalized=""): ###xset==None?
		prot_sheet,prot_repetitions=self.load_prot_profiles_and_split_into_reps(profiles)
		if normalized=="":
			normalized="sum"
		#getting xset and xset_profiles if all relations in file should be observed
		if xset_profiles is None:
			#create dataframe with all protein metabolite pairs as indices, what is og for training is protein for test data
			prot_index=list()
			metabolite_index=list()
			for i in prot_sheet.index: #i=uniprot ids for every protein profile
				for j in met_sheet.index: #j=metabolite ids
					prot_index.append(i)
					metabolite_index.append(j)
			xset_profiles=pd.DataFrame(columns=["protein","metabolite"])
			xset_profiles["protein"]=prot_index
			xset_profiles["metabolite"]=metabolite_index
			xset_profiles=xset_profiles.set_index(["protein","metabolite"],drop=True)
		for prot in xset.index.unique(): ####for prot in xset.index: #xset_profiles.index.get_level_values(0): #for every prot(=og) in index in xset_profiles
			#prot_profile=list() #combination (method) for all repetitions of all proteins in one og for all fractions in a list
			for coloumn in range(len(list(prot_repetitions[0]))): #all repetitions must have the same size and amount of fractions, coloumn=fraction
				#protein profiles
				combined_prot_fraction=list() #combination (method) for all repetitions of all proteins to one OG in one fraction
				if xset is None: #no ogs given, so only one protein per og_met
					for prot_df in prot_repetitions:
						combined_prot_fraction.append(prot_df.loc[prot,list(prot_df)[coloumn]])
				else:
					for uniprot_id in xset.loc[prot,"UniProt_IDs"].split(";"): #for every uniprot ID from given organism to corresponding OG
						prot_fraction=list() #list of repetitions for one protein uniprot ID
						for prot_df in prot_repetitions:
							try:
								if prot_df.loc[uniprot_id].sum()!=0:
									if not pd.isna(prot_df.loc[uniprot_id,list(prot_df)[coloumn]]) and type(prot_df.loc[uniprot_id,list(prot_df)[coloumn]])!=str: #if no number is given in excel sheet (NaN)
										eval("prot_fraction.append(prot_df.loc[uniprot_id,list(prot_df)[coloumn]]/prot_df.loc[uniprot_id]."+normalized+"())") #normalized fractions intensity
								else:
									prot_fraction.append(0.0)
									#print("sum of profile equals zero: "+uniprot_id)
									continue
							except (KeyError,ValueError): #if only one of IDs is not in 
								try: #if one uniprot id occurs several times, skip it
									uniprot_id2=list(filter(lambda s: uniprot_id in s,prot_df.index.tolist()))
									if len(uniprot_id2)==1:
										if prot_df.loc[uniprot_id2[0]].sum()!=0:
											eval("prot_fraction.append(prot_df.loc[uniprot_id2[0],list(prot_df)[coloumn]]/prot_df.loc[uniprot_id2[0]]."+normalized+"())") #normalized fractions intensity
										else:
											prot_fraction.append(0.0)
											#print("sum of profile equals zero: "+uniprot_id)
											continue
									else:
										break
								except:
									break
						if prot_fraction==list(): #if uniprot ID filtered out or occurs multiple times, do not take it into consideration
							continue
						else:
							combined_prot_fraction.append(method(prot_fraction))
				#prot_profile.append(method(combined_prot_fraction))
				if combined_prot_fraction==list():
					#xset_profiles=xset_profiles.drop(index=og)
					#xset_profiles=xset_profiles.drop(og,axis=0)
					if prot in xset_profiles.index: #may be already excluded from other organism
						for met in xset_profiles.loc[prot].index:
							xset_profiles=xset_profiles.drop((prot,met))
					print("skipped following Orthologous Group due to multiple times occuring id or NaNs in profile: "+str(prot))
					break
				else:
					#xset_profiles.ix[og_met,self.organism+"_fraction_"+str(coloumn+1)]=[method(combined_prot_fraction),method(met_fraction)]
					#xset_profiles.loc[prot,self.organism+"_proteins_"+str(coloumn)]=method(combined_prot_fraction) #proteins only
					xset_profiles.loc[xset_profiles.index.isin([prot],level=0),self.organism+"_proteins_"+str(coloumn)]=method(combined_prot_fraction) ####
		return(xset_profiles)
		
		
	
	#get metabolite profiles pooled over repetitions and normalize them
	#xset_profiles is dataframe with protein-metabolite pairs to extract
	#profiles={"raw","deconvoluted"}
	#normalized can be "sum" to normalize profiles to an area of one or "max" to normalize to maximum peak intensity
	def extract_met_profiles(self,xset_profiles=None,method=np.mean,profiles="raw",normalized=""):
		met_sheet,met_repetitions=self.load_met_profiles_and_split_into_reps(profiles)
		if normalized=="":
			normalized="sum"
		#getting xset and xset_profiles if all relations in file should be observed
		if xset_profiles is None:
			#create dataframe with all protein metabolite pairs as indices, what is og for training is protein for test data
			prot_sheet,prot_repetitions=self.load_prot_profiles_and_split_into_reps(profiles)
			prot_index=list()
			metabolite_index=list()
			for i in prot_sheet.index: #i=uniprot ids for every protein profile
				for j in met_sheet.index: #j=metabolite ids
					prot_index.append(i)
					metabolite_index.append(j)
			xset_profiles=pd.DataFrame(columns=["protein","metabolite"])
			xset_profiles["protein"]=prot_index
			xset_profiles["metabolite"]=metabolite_index
			xset_profiles=xset_profiles.set_index(["protein","metabolite"],drop=True)
		#metabolite profiles
		errorneous_metas=set()
		met_previous="none"
		for met in xset_profiles.index.get_level_values(1).unique(): #for every metabolite occuring in index in xset_profiles
			for coloumn in range(len(list(met_repetitions[0]))): #all repetitions must have the same size and amount of fractions, coloumn=fraction
				#xset_profiles.loc[:,self.organism+"_metabolites_"+str(coloumn)]=np.nan #np.empty((len(xset_profiles),0)).tolist()
				met_fraction=list()
				for met_df in met_repetitions:
					try: ###
						if met_df.loc[met].sum()!=0:
							if not pd.isna(met_df.loc[met,list(met_df)[coloumn]]) and type(met_df.loc[met,list(met_df)[coloumn]])!=str: #if NaN in excelsheet, skip the value (e.g. one chromatogram in one experiment missing, so do not include in pool)
								eval("met_fraction.append(met_df.loc[met,list(met_df)[coloumn]]/met_df.loc[met]."+normalized+"())") #combination/method of all repetitions for one metabolite (normalized)
						else:
							met_fraction.append(0.0)
							if met!=met_previous:
								#print("sum of profile equals zero: "+met)
								met_previous=met
							continue
					except:### possibly metabolite name occuring more than once in excel sheet
						errorneous_metas.add(met)
						#print("check for "+met)###
				if met_fraction==[]:
					errorneous_metas.add(met)
				else:
					xset_profiles.loc[xset_profiles.index.isin([met],level=1),self.organism+"_metabolites_"+str(coloumn)]=method(met_fraction)
		print("check for following metabolites:")
		print(errorneous_metas) ###
		return(xset_profiles)
		
	
	#calculate RV coefficient as similarity measure over matrices (used to get correlation between profiles)
	#repetitionlist is list of matrices
	def rv_coeff(self,repetitionlist):
		reps=list(map(lambda x: "Experiment_"+str(x+1),range(len(repetitionlist))))
		rv=pd.DataFrame(index=reps,columns=reps)
		for i in range(len(repetitionlist)):
			for j in range(i,len(repetitionlist)):
				#covmat_ij=np.cov(repetitionlist[i],repetitionlist[j])
				#covmat_ji=np.cov(repetitionlist[j],repetitionlist[i])
				#covv=np.trace(covmat_ij*covmat_ji) # scalar-valued covariance
				#vav_i=np.trace(covmat_ij*covmat_ij) #scalar-valued variance
				#vav_j=np.trace(covmat_ji*covmat_ji)
				#rv.ix[i,j]=covv/np.sqrt(vav_i*vav_j)
				rv.iloc[i,j]=np.multiply(repetitionlist[i],repetitionlist[j]).sum().sum()/np.sqrt(np.multiply(repetitionlist[i],repetitionlist[i]).sum().sum()*np.multiply(repetitionlist[j],repetitionlist[j]).sum().sum()) #IX
				rv.iloc[j,i]=rv.iloc[i,j] #IX
		return(rv)
		
		
		
	#calculate rv-coefficient for correlation of profiles for repetitions
	def profiles_correlation_reps(self,profiles="raw"):
		#load sheets
		prot_sheet,prot_repetitions=self.load_prot_profiles_and_split_into_reps(profiles)
		met_sheet,met_repetitions=self.load_met_profiles_and_split_into_reps(profiles)
		#if one chromatogramm missing, delete columns of same size in all experiments for proteins and metabolites
		##find columns with NaNs
		drop_cols=list()
		for prep in prot_repetitions:
			drop_cols+=list(map(lambda c: c.split(".")[-1],prep.columns[prep.isna().sum()==len(prep)]))
		for mrep in met_repetitions:
			drop_cols+=list(map(lambda c: c.split(".")[-1],mrep.columns[mrep.isna().sum()==len(mrep)]))
		##delete found columns
		for p in range(len(prot_repetitions)):
			prot_repetitions[p]=prot_repetitions[p].drop(list(filter(lambda c: c.split(".")[-1] in drop_cols,prot_repetitions[p].columns)),axis=1)
		for m in range(len(met_repetitions)):
			met_repetitions[m]=met_repetitions[m].drop(list(filter(lambda c: c.split(".")[-1] in drop_cols,met_repetitions[m].columns)),axis=1)
		#if errorneous output: exclusive nans in profiles are not dropped
		#calculate rv-coefficent for proteins and metabolites
		rv_prot=self.rv_coeff(prot_repetitions)
		rv_met=self.rv_coeff(met_repetitions)
		#if necessary, add calculation of histogram for correlations inside OGs
		
		return([rv_prot,rv_met])
		
		
		
	#calculates the correlation of profiles within one orthologous group per file and returns list of all correlations, which can then be plotted as a histogram
	#also calculates the number of proteins per OG
	def profiles_correlation_ogs(self,xset,profiles="raw",method=np.mean):
		#load sheets
		if not hasattr(self,"prot_repetitions"):
			self.prot_sheet,self.prot_repetitions=self.load_prot_profiles_and_split_into_reps(profiles)
		#met_sheet,met_repetitions=self.load_met_profiles_and_split_into_reps(profiles)
		correlations=pd.DataFrame(index=xset.index,columns=["correlations"])
		correlations["correlations"]=np.empty((len(correlations),0)).tolist() #initialize empty lists
		num_prots=pd.DataFrame(columns=["num_prots_in_OG",self.organism],data=[[1,0]]).set_index("num_prots_in_OG")
		for og in xset.index: #for every orthologous group
			df=pd.DataFrame(columns=xset.loc[og,"UniProt_IDs"].split(";"))
			#counting assigned proteins and add them
			if len(xset.loc[og,"UniProt_IDs"].split(";")) in num_prots.index:
				num_prots.loc[len(xset.loc[og,"UniProt_IDs"].split(";")),self.organism]+=1
			else:
				num_prots.loc[len(xset.loc[og,"UniProt_IDs"].split(";")),self.organism]=1
			if len(xset.loc[og,"UniProt_IDs"].split(";"))>1: #only need to do further calculations if more than one protein in orthologous group
				for prot in xset.loc[og,"UniProt_IDs"].split(";"):
					if prot in self.prot_sheet.index:
						prot_profile=list()
						for fraction in range(len(list(self.prot_repetitions[0]))):
							prot_fraction=list()
							for exp in self.prot_repetitions:
								try: #sometimes one uniprot ID occurs multiple times in data set, skip it :(
									prot_fraction.append(exp.loc[prot,list(self.prot_repetitions[0])[fraction]]) #IX
								except:
									break
							if prot_fraction==list():
								df=df.drop(prot,axis=1)
								break
							else:
								prot_profile.append(method(prot_fraction))
						if prot_profile!=list():
							df[prot]=prot_profile
				correlations.loc[og,"correlations"]=df.corr().where(np.triu(np.ones(df.corr().shape),k=1).astype(np.bool)).stack().values.tolist() #upper triangle without diagonal of correlation matrix
		return([correlations,num_prots])
	
	
	
	
	#for given list of uniprot IDs, calculate correlation
	def profiles_correlation_proteins(self,proteins,profiles="raw",method=np.mean):
		df=pd.DataFrame(columns=proteins)
		if not hasattr(self,"prot_repetitions"):
			self.prot_sheet,self.prot_repetitions=self.load_prot_profiles_and_split_into_reps(profiles)
		for prot in proteins:
			if prot in self.prot_sheet.index:
				prot_profile=list()
				for fraction in range(len(list(self.prot_repetitions[0]))):
					prot_fraction=list()
					for exp in self.prot_repetitions:
						try: #sometimes one uniprot ID occurs multiple times in data set, skip it :(
							prot_fraction.append(exp.loc[prot,list(self.prot_repetitions[0])[fraction]]) #IX
						except:
							break
					if prot_fraction==list():
						df=df.drop(prot,axis=1)
						break
					else:
						prot_profile.append(method(prot_fraction))
				if prot_profile!=list():
					df[prot]=prot_profile
		return(df.corr())
	
	
	
	#averages profiles for a list of metabolites and a dataframe of OGs and proteins and returns list where list[0] is average protein profile and list[1] is average metabolite profile 
	#metabolites need to be padded with zeros in front
	def average_meta_and_prot_profiles(self,ogs,metabolites,profiles="raw",normalized=""):
		if normalized=="":
			normalized="sum"
		#OGs
		prot_sheet,prot_repetitions=self.load_prot_profiles_and_split_into_reps(profiles)
		all_og_profiles=list()
		for og in ogs.index:
			og_rep=list()
			for prep in prot_repetitions:
				for prot in ogs.loc[og,"proteins"].split(";"):
					if prot in prep.index:
						try:
							eval("og_rep.append(prep.loc[prot]/prep.loc[prot]."+normalized+"())")
						except:
							pass #multiple times occuring uniprot id
			if og_rep!=[]:
				all_og_profiles.append(list(np.nanmean(og_rep,axis=0)))
		
		#metabolites
		met_sheet,met_repetitions=self.load_met_profiles_and_split_into_reps(profiles)
		all_metabolite_profiles=list()
		for met in self.padd_cids(metabolites):
			met_rep=list()
			for mrep in met_repetitions:
				if sum(mrep.loc[met])!=0:
					eval("met_rep.append(mrep.loc[met]/mrep.loc[met]."+normalized+"())")
				else:
					met_rep.append([0]*len(mrep.loc[met]))
			all_metabolite_profiles.append(np.mean(met_rep,axis=0))
		
		return([list(np.mean(all_og_profiles,axis=0)),list(np.mean(all_metabolite_profiles,axis=0))])
		
	
	
	#deletes proteins and metabolites with zero-profiles in excelfile and saves changes to excelfile
	#before applying, make manually a copy of file in advance to get back headers and format
	def change_excel_sheets(self):
		writer=pd.ExcelWriter(self.inputfile)
		#read excel file and change sheets as wanted (here: remove rows with all zeros in "ExperimentE columns)
		met_sheet=self.load_sheet("Metabolites_profiles_rep_raw")
		met_expcols=list(filter(lambda c: "Experiment" in c,met_sheet.columns))
		new_met_sheet=met_sheet.loc[met_sheet[met_expcols].sum(axis=1)!=0]
		prot_sheet=self.load_sheet("Protein_profiles_rep_raw")
		prot_expcols=list(filter(lambda c: "Experiment" in c,prot_sheet.columns))
		new_prot_sheet=prot_sheet.loc[prot_sheet[prot_expcols].sum(axis=1)!=0]
		#save dataframes to sheets
		new_met_sheet.to_excel(writer,"Metabolites_profiles_rep_raw",index=False)
		new_prot_sheet.to_excel(writer,"Protein_profiles_rep_raw",index=False)
		writer.save()
		print("Change back column headers and formatting manually in Excel!")
		return
	
	
	
	def helpfun():
		print("Provide as optional argument the filename")




###############################################################
###############################################################
###############################################################

#class for interactions between one file with experimental data and the corresponding data base
#used to generate positive training set
class DBFileInteractions():
	#simulation only to check if script runs well
	#'experimental'=folder with experimental data files
	#'databases'=folder with corresponding databases
	def __init__(self,expfilename=None,simulation=False,experimental="../experimental_data/",databases="../databases/"):
		self.simulation=simulation
		self.experimental=experimental
		self.databases=databases
		if expfilename!=None:
			self.expfile=self.get_expfile(expfilename)
			self.database=self.get_database(expfilename)
			
		
	
	#loads experimental data to self.expfile
	def get_expfile(self,expfilename):
		return(FileHandler(expfilename,experimental=self.experimental,databases=self.databases))
		
		
	#loads database to self.database
	def get_database(self,expfilename):
		#return(dh.DBHandler(expfilename[:-5].split("/")[-1].split("_")[0])) #1.split: extract filename from path, 2.split: take name till first underscore (to load same database for A.thaliana_rosettes and A.thaliana_cell_culture)
		return(DBHandler(expfilename,simulation=self.simulation,experimental=self.experimental,databases=self.databases))
		
	#extracts from database PMIs for proteins which occur in experimental data (e.g. no hydrophobic/TM proteins)
	def intersect(self,dbset,expdataset):
		return(set(dbset & expdataset))
		
	#get union of two sets
	def union(self,dbset,expdataset):
		return(dbset | expdataset)
		
	#gets set for known interactions
	def find_positive_candidates(self,pooling="intersect"):
		positive_candidates=eval("self."+pooling+"(self.database.get_proteins(),self.expfile.get_proteins())") #:=proteins suitable for (+) set, pooling (intersect or union) between database and experimental data
		print(str(len(positive_candidates))+" proteins with interactionscore >= "+str(self.database.score_cutoff)+" found in intersection of database and experimental data")
		return(positive_candidates) #return UniProtIDs
		
		
	#get orthologs and save dataframe to file
	def extract_and_save_orthologs(self,uniprot_set,outfilename,overwrite=False,include_metas=True):
		if overwrite==True or not os.path.isfile(outfilename):
			try:
				string_set=self.database.translate_proteins(uniprot_set,"ID",["ID","STRING_ID","ENSEMBLGENOME_PRO_ID","GENENAME"],save=False) #"../databases/"+self.database.organism+"_intersection.tsv")
			except:
				print("ConnectionError")
				t=5
				while True:
					time.sleep(t)
					try:
						string_set=self.database.translate_proteins(uniprot_set,"ID",["ID","STRING_ID","ENSEMBLGENOME_PRO_ID","GENENAME"],save=False) #save="../databases/"+self.database.organism+"_intersection.tsv")
						break
					except:
						t=t+5
						#time.sleep(10)
			try:
				string_uniprot_cog=pd.read_csv("../databases/"+self.database.organism.split("_")[0]+"_string_uniprot_cog.tsv",sep="\t") ### new: .split("_")[0], to make only one per organism
				string_uniprot_cog=string_uniprot_cog.set_index("STRING_ID")
			except:
				string_uniprot_cog=self.database.get_orthologs(save=True)
			string_orthos=string_uniprot_cog.loc[string_set["STRING_ID"].dropna()]
			string_orthos.to_csv(outfilename,header=True,index=True,sep="\t")
			if include_metas==True:
				string_orthos["CID"]=np.nan
				string_orthos_with_metabolites=self.database.find_metabolites(string_orthos)
				#string_orthos=string_orthos.dropna() #get rid of rows where are NaNs
				string_orthos_with_metabolites.to_csv(outfilename,header=True,index=True,sep="\t")
				#return(string_orthos)
		else:
			print("orthologs for "+self.database.organism+" found and not overwritten")
			#return(pd.read_csv(outfilename,sep="\t"))
		





###############################################################
###############################################################
###############################################################

#class for generating a set with common methods for PositiveSet and NegativeSet
# includes generation of feature engineered profiles, class annotations, etc. 


class XSet(IDTranslations):
	#feature={"","_product","_full_convoluted","_same_convoluted","_same_crosscorr","_same_crosscorr_not_normalized_binsize4", etc.}, describes feature engineering procedure
	#normalized is normalization of every profile by itself, normalized2 normalization of pooled profile (over repetitions or OGs)
	#experimental is folder with experimental data files
	#databases is folder with databases, here are training subsets saved
	#methods is list of pooling methods
	#analyses is folder for results
	#proteinwise = True for analysis of protein-metabolite pairs, False for analysis of OG-metabolite pairs
	def __init__(self,simulation=False,overwrite=False,experimental="../experimental_data/",databases="../databases/",analyses="../analyses/",feature="",methods=[np.mean],normalized="",normalized2=None,proteinwise=False):
		self.simulation=simulation
		self.overwrite=overwrite
		self.experimental=experimental
		self.databases=databases
		self.analyses=analyses
		self.feature=feature
		self.normalized=normalized
		if normalized2 is None:
			self.normalized2=normalized
		else:
			self.normalized2=normalized2
		self.proteinwise=proteinwise #if proteinwise=True, compute XSet for protein-metabolite, else OG-metabolite pair
		if self.proteinwise==True:
			self.protcol="protein"
		else:
			self.protcol="OG"
		if "full" in feature:
			self.mode="full"
		else:
			self.mode="same"
		self.methods=methods
		#if folder not existing, create them:
		if not os.path.isdir(self.analyses):
			os.makedirs(self.analyses)
			print("folder "+self.analyses+" created")
		if not os.path.isdir(self.databases):
			os.makedirs(self.databases)
			print("folder "+self.databases+" created")
		if not os.path.isdir(self.experimental):
			os.makedirs(self.experimental)
			print("folder "+self.experimental+" created")
		self.expfiles=self.get_expfiles()
		
		
	
	
	#reset feature
	def set_feature(self,feature):
		self.feature=feature
		return
	
	#extraction of orthologs?
	#in PositiveSet get_orthologs_between_species, in NegativeSet find_orthologs_between_files
	
	#loops through files with experimental data, finds orthologous group for every protein and creates intersection between proteins of all files, saves OG and corresponding STRING ID for each
	#and returns set of orthologous groups that are present in every file
	def find_orthologs_between_files(self):
		cogs_all_orgs=list()
		for i in range(len(self.expfiles)):
			expfilename=self.expfiles[i]
			outfilename=self.databases+expfilename[:-5].split("/")[-1]+"_string_cog_from_expdata.tsv"
			try: #try to load dataframe
				df_orthologs=pd.read_csv(outfilename,sep="\t")
				print("file with STRING IDs and corresponding orthologous groups for "+expfilename[:-5].split("/")[-1]+" found")
			except FileNotFoundError: 
				expfile=FileHandler(expfilename,experimental=self.experimental,databases=self.databases)
				db=DBHandler(expfilename)
				#extract proteins from experimental data and translate them to STRING IDs
				proteins=expfile.get_proteins()
				string_ids=set(db.translate_proteins(proteins,in_id_type="ID",out_id_types=["STRING_ID"]).STRING_ID.values.tolist())
				string_ids.remove(np.nan)
				#find orthologous group for STRING IDs and save them
				df_orthologs=db.get_orthologs(list(string_ids))
				#del df_orthologs["CID"] #get rid of coloumn CID
				if "CID" in df_orthologs:
					df_orthologs.drop("CID",axis=1,inplace=True)
				df_orthologs=df_orthologs.dropna() #get rid of missing entries (no OG found for STRING ID)
				df_orthologs.to_csv(outfilename,header=True,index=True,sep="\t")
				print("file with STRING IDs and corresponding orthologous groups for "+expfilename[:-5].split("/")[-1]+" created")
			cogs_list=df_orthologs.COG.values.tolist()
			cogs_org=list()
			for og in cogs_list:
				cogs=str(og).split(";")
				for element in cogs:
					cogs_org.append(element)
					cogs_all_orgs.append(element)
			if i==0:
				cogs_intersect=set(cogs_org)
			else:
				cogs_intersect=set(set(cogs_org) & cogs_intersect)
			print("In "+expfilename[:-5].split("/")[-1]+" "+str(len(cogs_org))+" orthologous groups annotated")
		cogs_set=set(cogs_all_orgs)
		print("Over all excel sheets with experimental data "+str(len(cogs_set))+" orthologous groups found")
		print("In intersection over all excel sheets "+str(len(cogs_intersect))+" orthologous groups found")
		return(cogs_set)
	
	
	
	
	#builds intersection between metabolites from experimental data files,returns metabolites as set
	def meta_intersect_files(self):
		#build intersection of metabolites over all files of experimental data
		for i in range(len(self.expfiles)):
			expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
			if i==0:
				metabolites=expdata.get_metabolites() #set(expdata.translate_metabolites_offline()["CID"].values.tolist())
				#metabolites=set(list(itertools.chain.from_iterable(expdata.translate_metabolites_pubchem()["CIDs"].values.tolist())))
			else:
				metabolites=set(metabolites & expdata.get_metabolites()) #set(metabolites & set(expdata.translate_metabolites_offline()["CID"].values.tolist()))
				#metabolites=set(metabolites & set(list(itertools.chain.from_iterable(expdata.translate_metabolites_pubchem()["CIDs"].values.tolist()))))
		print(str(len(metabolites))+" metabolites in intersection between files of experimental data")
		return(metabolites)
		
	
	#builds intersection of protein UniProt IDs among experimental data files, returns proteins as set
	def prot_intersect_files(self):
		for i in range(len(self.expfiles)):
			expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
			if i==0:
				proteins=set(expdata.get_proteins())
			else:
				proteins=set(proteins & set(expdata.get_proteins()))
		print(str(len(proteins))+" proteins in intersection between files of experimental data")
		return(proteins)
		
		
	
	#builds union between metabolites from experimental data,returns metabolites as set
	def meta_union_files(self):
		#build intersection of metabolites over all files of experimental data
		for i in range(len(self.expfiles)):
			expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
			if i==0:
				metabolites=expdata.get_metabolites()#set(expdata.translate_metabolites_offline()["CID"].values.tolist())
				#metabolites=set(list(itertools.chain.from_iterable(expdata.translate_metabolites_pubchem()["CIDs"].values.tolist())))
			else:
				metabolites=metabolites.union(expdata.get_metabolites()) #metabolites.union(set(expdata.translate_metabolites_offline()["CID"].values.tolist()))
				#metabolites=set(metabolites & set(list(itertools.chain.from_iterable(expdata.translate_metabolites_pubchem()["CIDs"].values.tolist()))))
		print(str(len(metabolites))+" metabolites in union between files of experimental data")
		return(metabolites)
		
	
	#abstract class to generate and save xset
	@abc.abstractmethod
	def save_set(self,tabfile=None,df=None):
		"""generate and save set"""
		
	
	
	#tries to load existing set from file, if not found generates it
	def load_set(self,tabfile,method=None):
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			try:
				xset=pd.read_csv(tabfile,sep="\t",dtype={"CIDs":str})
			except: #maybe not neccessary, if column CIDs not existent
				xset=pd.read_csv(tabfile,sep="\t")
			xset=xset.set_index(list(xset)[0],drop=True)
			'''
			try:
				xset=xset.set_index(["OG","metabolite"],drop=True)###
			except KeyError:
				xset=xset.set_index(["OG","CIDs"],drop=True)###
			'''
		except (OverwriteException,FileNotFoundError):
			xset=self.save_set(tabfile,method=method)
		return(xset)
			
	
	
	
	#gets all uniprot_IDs into one textfile, which can be uploaded on Uniprot to query EC numbers and Gene Ontology
	def collect_uniprot_IDs(self,xset=None,predictions=None):
		if predictions is None:
			uniprots=list()
			for column in xset.columns:
				uniprots.append(";".join(xset[column].dropna().tolist()))
			newline_sep_uniprots="\n".join(";".join(uniprots).split(";"))
		elif xset is None:
			if self.proteinwise==True:
				uniprots=predictions.index.get_level_values(0).unique().tolist()
				newline_sep_uniprots="\n".join(uniprots)
		uniprot_file=open(self.analyses+"all_uniprot_ids.txt","w")
		uniprot_file.write(newline_sep_uniprots)
		uniprot_file.close()
		return
	
	
	
	#retrieve GO-Terms and EC-numbers for proteins
	def retrieve_EC_and_GO(self):
		print("1. run self.collect_uniprot_IDs(xset) or self.collect_uniprot_IDs(predictions)")
		print("2. upload "+self.analyses+"all_uniprot_ids.txt to https://www.uniprot.org/uploadlists/")
		print("3. From UniProtKB AC/ID to UniProtKB")
		print("4. select columns EC number, GO... etc")
		print("5. save result to "+self.analyses+"all_protein_annotations.tsv")
	
	
	
	#translates given uniprot IDs to KEGG IDs and classifies them with BRITE
	#org is KEGG organism identifier {"ath","sce","eco"}
	def classify_proteins(self,uniprot_ids=None,org=None):
		try:
			df=pd.read_csv(self.analyses+"uniprot_kegg_briteclasses_keggpathways.tsv",sep="\t").set_index("UniProt ID")
		except FileNotFoundError:
			if uniprot_ids is None:
				print("provide function arguments")
				return
			df=pd.DataFrame(index=uniprot_ids,columns=["KEGG ID","BRITE classes","KEGG pathways"])
			df.index.name="UniProt ID"
			for uniprot_id in uniprot_ids:
				#translate uniprot ID to KEGG ID
				conv_url="http://rest.kegg.jp/conv/"+org+"/uniprot:"+uniprot_id
				conv_resp=requests.post(conv_url)
				conv_out=conv_resp.text
				kegg_id=re.split("\t|\n",conv_out)[1]
				try:
					if kegg_id[3]==":": #check if ID was translated successfully
						df.loc[uniprot_id,"KEGG ID"]=kegg_id
						#query for BRITE classes
						class_url="http://rest.kegg.jp/link/brite/"+kegg_id
						class_resp=requests.post(class_url)
						try:
							class_out=re.split("\t|\n",class_resp.text)
							brite_classes=list(filter(lambda x: "br:" in x,class_out))
							df.loc[uniprot_id,"BRITE classes"]=";".join(brite_classes)
						except: 
							pass
						#query for KEGG pathway
						path_url="http://rest.kegg.jp/link/pathway/"+kegg_id
						path_resp=requests.post(path_url)
						try:
							path_out=re.split("\t|\n",path_resp.text)
							pathways=list(filter(lambda x: "path:" in x,path_out))
							df.loc[uniprot_id,"KEGG pathways"]=";".join(pathways)
						except: 
							pass
				except IndexError:
					pass
			df.to_csv(self.analyses+"uniprot_kegg_briteclasses_keggpathways.tsv",sep="\t",index=True,header=True)
		return(df)
	
	
	#assign classes to OGs in complete set
	def assign_classes_and_pathways_to_OGs(self,complete_set,kegg_df):
		exp_cols=complete_set.columns
		for cog in complete_set.index:
			#get union of UniProt IDs
			complete_set.loc[cog,"union UniProt IDs"]=";".join(set(";".join(complete_set.loc[cog,exp_cols]).split(";")))
			#assign classes to uniprot IDs in union
			union_brite_classes=set() #empty set
			union_pathways=set()
			i=0
			while True: #initialize intersect_brite_classes
				#initialize intersect_brite_classes as Brite classes of first uniprot ID that has Brite classes assigned (not nan)
				if complete_set.loc[cog,"union UniProt IDs"].split(";")[i] in kegg_df.index:
					if not pd.isna(kegg_df.loc[complete_set.loc[cog,"union UniProt IDs"].split(";")[i],"BRITE classes"]):
						intersect_brite_classes=set(kegg_df.loc[complete_set.loc[cog,"union UniProt IDs"].split(";")[i],"BRITE classes"].split(";")) 
						break
				if i<len(complete_set.loc[cog,"union UniProt IDs"].split(";"))-1:
					i+=1
				else: #if no UniProt ID has a brite class annotated, intersection is empty
					intersect_brite_classes=set()
					break
			while True: #initialize intersect_pathways
				#initialize intersect_pathways as pathways of first uniprot ID that has pathways assigned (not nan)
				if complete_set.loc[cog,"union UniProt IDs"].split(";")[i] in kegg_df.index:
					if not pd.isna(kegg_df.loc[complete_set.loc[cog,"union UniProt IDs"].split(";")[i],"KEGG pathways"]):
						intersect_pathways=set(kegg_df.loc[complete_set.loc[cog,"union UniProt IDs"].split(";")[i],"KEGG pathways"].split(";")) 
						break
				if i<len(complete_set.loc[cog,"union UniProt IDs"].split(";"))-1:
					i+=1
				else: #if no UniProt ID has a pathway annotated, intersection is empty
					intersect_pathways=set()
					break
			for uniprot_id in complete_set.loc[cog,"union UniProt IDs"].split(";"):
				if complete_set.loc[cog,"union UniProt IDs"].split(";")[i] in kegg_df.index:
					if uniprot_id in kegg_df.index:
						if not pd.isna(kegg_df.loc[uniprot_id,"BRITE classes"]): #if no class annotated to protein, skip it
							union_brite_classes=union_brite_classes | set(kegg_df.loc[uniprot_id,"BRITE classes"].split(";"))
							intersect_brite_classes=intersect_brite_classes & set(kegg_df.loc[uniprot_id,"BRITE classes"].split(";"))
						if not pd.isna(kegg_df.loc[uniprot_id,"KEGG pathways"]):
							union_pathways=union_pathways | set(kegg_df.loc[uniprot_id,"KEGG pathways"].split(";"))
							intersect_pathways=intersect_pathways & set(kegg_df.loc[uniprot_id,"KEGG pathways"].split(";"))
					else:
						print(uniprot_id+" missing in kegg_df")
			complete_set.loc[cog,"union BRITE classes"]=";".join(union_brite_classes)
			complete_set.loc[cog,"intersect BRITE classes"]=";".join(intersect_brite_classes)
			complete_set.loc[cog,"union KEGG pathways"]=";".join(union_pathways)
			complete_set.loc[cog,"intersect KEGG pathways"]=";".join(intersect_pathways)
		return(complete_set)
	
	
	
	#run after assign_classes_and_pathways_to_OGs to additionally assign EC classes
	def assign_EC_numbers_to_OGs(self,complete_set):
		try:
			all_protein_annotations=pd.read_csv(self.analyses+"all_protein_annotations.tsv",sep="\t").set_index("Entry")
			for og in complete_set.index:
				union_EC_numbers=set(all_protein_annotations.loc[complete_set.loc[og,"union UniProt IDs"].split(";"),"EC number"])
				if np.nan in union_EC_numbers: #if no EC number annotated
					union_EC_numbers.remove(np.nan)
				complete_set.loc[og,"union EC numbers"]=";".join(union_EC_numbers)
				if union_EC_numbers!=set():
					intersect_EC_numbers=set([list(union_EC_numbers)[0].split(".")[0]])
					for ec in union_EC_numbers:
						intersect_EC_numbers=set(ec.split(".")[0]) & intersect_EC_numbers
					complete_set.loc[og,"intersect EC class"]=";".join(intersect_EC_numbers)
				else:
					complete_set.loc[og,"intersect EC class"]=""
			return(complete_set)
		except FileNotFoundError:
			self.retrieve_EC_and_GO()
			return
	
	
	#assign EC numbers to proteins (for proteinwise analysis)
	#df is dataframe with uniprot ids as index, optionally with annotated BRITE classes and KEGG pathways
	def assign_EC_numbers_to_proteins(self,df):
		try:
			all_protein_annotations=pd.read_csv(self.analyses+"all_protein_annotations.tsv",sep="\t").set_index("Entry")
			for uniprot in df.index:
				if uniprot in all_protein_annotations.index:
					if not pd.isna(all_protein_annotations.loc[uniprot,"EC number"]):
						df.loc[uniprot,"EC number"]=all_protein_annotations.loc[uniprot,"EC number"]
						df.loc[uniprot,"EC class"]=all_protein_annotations.loc[uniprot,"EC number"].split(".")[0]
			return(df)
		except FileNotFoundError:
			self.retrieve_EC_and_GO()
			return
	
	
	#finds pathways for metabolites (KEGG compound IDs) on KEGG 
	def find_pathways_metabolites(self,kegg_ids):
		df=pd.DataFrame(index=kegg_ids,columns=["KEGG pathways"])
		df.index.name="KEGG_CPD_ID"
		for kegg_id in kegg_ids: #for cid in cids:
			'''
			#translate PubChem CID to KEGG ID
			conv_url="http://rest.kegg.jp/conv/compound/pubchem:"+str(cid)
			conv_resp=requests.post(conv_url)
			conv_out=conv_resp.text
			kegg_id=re.split("\t|\n",conv_out)[1]#[4:]
			'''
			#query for KEGG pathway
			path_url="http://rest.kegg.jp/link/pathway/cpd:"+kegg_id
			path_resp=requests.post(path_url)
			try:
				path_out=re.split("\t|\n",path_resp.text)
				pathways=list(filter(lambda x: "path:" in x,path_out))
				df.loc[kegg_id,"KEGG pathways"]=";".join(pathways)
			except: 
				pass
		df.to_csv(self.analyses+"metabolite_kegg_pathways.tsv",sep="\t",index=True,header=True)
		return(df)
	
	
	
	#query CIDs on pubchem for metabolite classification
	def retrieve_CID_classes(self,cids=None):
		try: 
			metabolite_df=pd.read_csv(self.analyses+"metabolite_classifications.tsv",sep="\t")
		except FileNotFoundError:
			if cids is None:
				complete_set=pd.read_csv(self.databases+"complete_set_normprofiles_"+cs.feature+"_"+method.__name__+".tsv",sep="\t")
				cids=set(complete_set["metabolite"].tolist())
			metabolite_df=pd.DataFrame(index=cids)
			for cid in cids:
				url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(cid)+"/classification/JSON"
				resp=requests.post(url)
				out=resp.json()
				if "Hierarchies" in out:
					for tree in out["Hierarchies"]["Hierarchy"]:
						nodes=list()
						for node in tree["Node"]:
							try:
								nodes.append(node["Information"]["Name"])
							except KeyError:
								nodes.append(node["Information"]["Description"][0])
						metabolite_df.loc[cid,tree["SourceName"]]=";".join(nodes)
				else:
					print(cid)
			metabolite_df_trimmed=metabolite_df.dropna(axis=1)
			#pubchem tree is no classification tree, but page record
			metabolite_df.to_csv(self.analyses+"metabolite_classifications.tsv",index=True,header=True,sep="\t")
		return(metabolite_df)
	
	
	
	
	
	#stacks elution profiles, calls for every excelsheet (organism) fh.extract_profiles() and stacks/concatenates them over organisms 
	def stack_profiles(self,xset,method=np.mean,profiles="raw"):
		#create index for dataframe with profiles (combined index of OG and metabolite)
		statusfile=open("./status.txt","w") #print progress to file (necessary when no print statements executed, e.g. when working on calc9)
		prot_index=list()
		metabolite_index=list()
		for i in xset.index: #i=OG
			cids_list=str(xset.loc[i,"CIDs"]).split(";")
			for j in cids_list: #j=CID
				prot_index.append(i)
				metabolite_index.append(str(int(j)))
		xset_profiles=pd.DataFrame(columns=[self.protcol,"metabolite"])
		xset_profiles[self.protcol]=prot_index
		xset_profiles["metabolite"]=metabolite_index
		xset_profiles=xset_profiles.set_index([self.protcol,"metabolite"],drop=True)
		
		#initialize concatenated profiles as empty lists
		#xset_profiles["concat protein profile"]=np.empty((len(xset_profiles),0)).tolist() ###
		#xset_profiles["concat metabolite profile"]=np.empty((len(xset_profiles),0)).tolist() ###
		#extract combined profiles for every organism to xset_profiles
		
		expfilenames=self.expfiles
		statusfile.close()
		
		for expfilename in expfilenames: #self.expfiles
			expdata=FileHandler(expfilename,experimental=self.experimental,databases=self.databases)
			organism=expfilename[:-5].split("/")[-1]
			if self.proteinwise==True:
				xset_org=xset
			else:
				xset_org=xset[[organism,"CIDs"]]
			xset_org.columns=["UniProt_IDs","CIDs"]
			#if COG/protein skipped due to missing ID in some experimental data file, take it out of analysis
			for prot in xset_org.index:
				if prot not in xset_profiles.index:
					xset_org=xset_org.drop(prot)
			statusfile=open("./status.txt","a")
			statusfile.write(time.asctime(time.localtime(time.time()))+"\textracting profiles for "+organism+"...\n")
			statusfile.close()
			xset_profiles=expdata.extract_prot_profiles(xset_org,xset_profiles,method,profiles,normalized=self.normalized)
			xset_profiles.to_csv(self.databases+"xset_profiles_beta.tsv",header=True,index=True,sep="\t") #backup 
			statusfile=open("./status.txt","a")
			statusfile.write(time.asctime(time.localtime(time.time()))+"\tprotein profiles for "+organism+" extracted\n")
			statusfile.close()
			xset_profiles=expdata.extract_met_profiles(xset_profiles,method,profiles,normalized=self.normalized)
			xset_profiles.to_csv(self.databases+"/xset_profiles_beta.tsv",header=True,index=True,sep="\t") #backup
			statusfile=open("./status.txt","a")
			statusfile.write(time.asctime(time.localtime(time.time()))+"\tmetabolite profiles for "+organism+" extracted\n")
			statusfile.close()
		return(xset_profiles)
		
		
		
	#normalizes profile for each organism to an area of 1 of the concatenated profiles
	def normalize_concat_profiles_per_org(self,xset_profiles):
		if self.normalized2=="":
			normalized="sum"
		else:
			normalized=self.normalized2
		for expfilename in self.expfiles:
			organism=expfilename[:-5].split("/")[-1]
			columns=list(filter(lambda col: organism in col,xset_profiles.columns.tolist()))
			metabolite_cols=list(filter(lambda col: "metabolites" in col,columns))
			protein_cols=list(filter(lambda col: "proteins" in col,columns))
			for i in xset_profiles.index:
				if xset_profiles.loc[i,metabolite_cols].sum()!=0:
					xset_profiles.loc[i,metabolite_cols]=eval("xset_profiles.loc[i,metabolite_cols]/xset_profiles.loc[i,metabolite_cols]."+normalized+"()")
				if xset_profiles.loc[i,protein_cols].sum()!=0:
					xset_profiles.loc[i,protein_cols]=eval("xset_profiles.loc[i,protein_cols]/xset_profiles.loc[i,protein_cols]."+normalized+"()")
		return(xset_profiles)
	
	
	
	#constructs mean of metabolite and protein profiles for given xset and returns list where list[0] is over organisms concatenated protein profile and list[1] metabolite profile
	def average_meta_and_prot_profiles_per_file(self,xset,profiles="raw"):
		#all metabolites in xset
		metabolites=set(itertools.chain.from_iterable(map(lambda x: x.split(";"),xset["CIDs"].unique().tolist())))
		#all OGs in xset
		ogs=pd.DataFrame(index=xset.index,columns=["proteins"])
		#get and concatenate profiles
		avg_og_profile=list()
		avg_met_profile=list()
		for expfilename in self.expfiles:
			expdata=FileHandler(expfilename,experimental=self.experimental,databases=self.databases)
			organism=expfilename[:-5].split("/")[-1]
			print(organism) #####
			ogs["proteins"]=xset[organism]
			#proteins=set(itertools.chain.from_iterable(map(lambda x: x.split(";"),xset[organism].unique().tolist())))
			avg_og_profile_org,avg_met_profile_org=expdata.average_meta_and_prot_profiles(ogs,metabolites,profiles,normalized=self.normalized)
			avg_og_profile.append(avg_og_profile_org)
			avg_met_profile.append(avg_met_profile_org)
		return([list(itertools.chain.from_iterable(avg_og_profile)),list(itertools.chain.from_iterable(avg_met_profile))])
	
	
	#multiplies corresponding fraction of protein an metabolite and normalizes (optionally) them
	#or convolutes protein and metabolite profiles
	#input is x_set with profiles and ["OG","metabolite"] as index, interaction can be included
	#mode is parameter for convolution: {"same","full"}
	def refeaturize_x_set(self,x_set):
		if len(x_set)==0:
			print("empty DataFrame given to refeaturize_x_set()")
			return(x_set)
		else:
			if self.normalized=="":
				normalized="sum"
			else:
				normalized=self.normalized
			if self.feature=="":
				print("no feature selected, so no feature engineering will be performed")
				return
			organisms=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
			#initialize x_set_fracs dataframe
			x_set_fracs=pd.DataFrame(index=x_set.index) 
			#x_set_fracs["OG"]=x_set.index.get_level_values(0)
			#x_set_fracs["metabolite"]=x_set.index.get_level_values(1)
			#x_set_fracs=x_set_fracs.set_index(["OG","metabolite"],drop=True)
			if "convoluted" in self.feature or "crosscorr" in self.feature:
				for org in organisms:
					#separate into one dataframe for protein and one for metabolite profiles
					metabolite_cols=list(filter(lambda col: org+"_metabolites" in col,x_set.columns))
					protein_cols=list(filter(lambda col: org+"_proteins" in col,x_set.columns))
					metabolite_xset=x_set[metabolite_cols]
					protein_xset=x_set[protein_cols]
					exp_cols=list()
					if self.mode=="same":
						num_fractions=len(list(filter(lambda c: org in c,x_set.columns)))/2
					elif self.mode=="full":
						num_fractions=len(list(filter(lambda c: org in c,x_set.columns)))-1
					for fraction in range(int(num_fractions)):
						exp_cols.append(org+"_"+str(fraction))
					x_set_fracs_org=pd.DataFrame(index=x_set.index,columns=exp_cols)
					for i in range(len(x_set.index)):
						if "convoluted" in self.feature:
							x_set_fracs_org.iloc[i]=np.convolve(protein_xset.iloc[i,:].tolist(),metabolite_xset.iloc[i,:].tolist(),mode=self.mode) #mode "same" or "full", same truncates to length of single profile, full has len of sum of both profiles -1
						elif "crosscorr" in self.feature:
							x_set_fracs_org.iloc[i]=np.correlate(protein_xset.iloc[i,:].tolist(),metabolite_xset.iloc[i,:].tolist(),mode=self.mode) 
					if "derivative" in self.feature: #calculate difference between every fraction and its previous fraction
						x_set_fracs_org=x_set_fracs_org.diff(axis=1)
						x_set_fracs_org=x_set_fracs_org.drop(list(x_set_fracs_org)[0],axis=1)
					if "binsize" in self.feature:
						binsize=int(re.search("[0-9]+",self.feature)[0])
						middlecol=np.median(range(x_set_fracs_org.shape[1]))
						x_set_fracs_org=x_set_fracs_org.iloc[:,int(np.ceil(middlecol-binsize/2)):int(np.floor(middlecol+binsize/2)+1)]
					x_set_fracs=x_set_fracs.join(x_set_fracs_org)
			else:
				#x_set_fracs=pd.DataFrame(index=x_set.index)
				for org in organisms:
					num_fractions=len(list(filter(lambda c: org in c,x_set.columns)))/2
					for fraction in range(int(num_fractions)):
						#multiply fraction pairs
						if "_product" in self.feature: #case: multiplication: "_product"
							x_set_fracs[org+"_"+str(fraction)]=x_set[org+"_metabolites_"+str(fraction)]*x_set[org+"_proteins_"+str(fraction)]
						#if necessary, include here other methods
					
			#normalize x_set_fracs
			if not "not_normalized" in self.feature:
				for org in organisms:
					org_cols=list(filter(lambda o: org in o,x_set_fracs.columns))
					for i in x_set_fracs.index:
						try: #case dataframe
							if x_set_fracs.loc[i,org_cols].sum(axis=1).iloc[0]!=0:
								x_set_fracs.loc[i,org_cols]=eval("x_set_fracs.loc[i,org_cols]/x_set_fracs.loc[i,org_cols]."+normalized+"(axis=1).iloc[0]")
							else:
								x_set_fracs.loc[i,org_cols]=x_set_fracs.loc[i,org_cols]
						except ValueError: #case Series
							if x_set_fracs.loc[i,org_cols].sum()!=0:
								x_set_fracs.loc[i,org_cols]=eval("x_set_fracs.loc[i,org_cols]/x_set_fracs.loc[i,org_cols]."+normalized+"()")
							else:
								x_set_fracs.loc[i,org_cols]=x_set_fracs.loc[i,org_cols]
			if "interaction" in x_set.columns:
				x_set_fracs["interaction"]=x_set["interaction"]
			return(x_set_fracs)
	
	
	

###############################################################
###############################################################
###############################################################

# class to build positive or negative sampled set (subsets of training set)
# with aid of one given set (y_set) predictions are made by random sampling of profiles as x_set
# those pairs always predicted as "x" are included returned as sampled x_set

#x_set is set which needs to be sampled, y_set is the one which is available
class XSet_sampled(XSet):
	
	def __init__(self,x,y_set,balanced=True,simulation=False,overwrite=False,experimental="../experimental_data/",databases="../databases/",analyses="../analyses/",feature="",normalized="",normalized2=None,proteinwise=False):
		super().__init__(simulation=simulation,overwrite=overwrite,experimental=experimental,databases=databases,analyses=analyses,feature=feature,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
		self.x=x #definition of this set: {"positive","negative"}
		self.balanced=balanced
		self.y_set=y_set.dropna() #other set
		if len(self.y_set)==0:
			print("given set is empty")
	
	
	#remaining_set is complete_set excluding y_set, n_samples determines how many xs should be sampled
	def sample_xs(self,remaining_set,n_samples):
		indices=sample(range(len(remaining_set)),n_samples)
		x_set=remaining_set.iloc[indices,:]
		return(x_set)
		
		
	#Random Forest classifier
	def rf_clf(self):
		return(RandomForestClassifier(oob_score=True))
		#return(RandomForestClassifier())
	
	
	#SVM classifier, kernel={"linear",,"poly,"rbf"}
	def svm_clf(self,C=1,kernel="linear",param=None): 
		return(svm.SVC(C=C,probability=False,kernel=kernel))
		
	
	#train classifier several times with sampled negative training data and extract constistently predicted as xs
	#n_trainings determines how many times classifiers are trained to get
	#classifiers=["rf_clf()","svm_clf(kernel='linear')","svm_clf(kernel='rbf')"]
	def train_to_get_xs(self,n_trainings,classifiers=["svm_clf(kernel='linear')"],method=np.mean):
		complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
		complete_set=complete_set.set_index([self.protcol,"metabolite"])
		#construct set to sample from
		n_samples=len(self.y_set)
		intersect_indices=set(set(self.y_set.index.tolist()) & set(complete_set.index.tolist()))
		remaining_set=complete_set.drop(intersect_indices)
		#exclude metabolites not occuring in y_set
		if self.balanced==True:
			metas_in_yset=set(self.y_set.index.get_level_values(1))
			remaining_set=remaining_set.loc[list(filter(lambda i: i[1] in metas_in_yset,remaining_set.index.tolist()))]
		for classifier in classifiers: #if several classifiers are given, take intersecting predictions among them
			for i in range(n_trainings):
				x_set=self.sample_xs(remaining_set,n_samples)
				clf=eval("self."+classifier)
				clf.fit(self.y_set.append(x_set),[True]*len(self.y_set)+[False]*len(x_set)) #x is labelled False, y is labelled True
				predictions=clf.predict(remaining_set)
				false_predictions=list(map(lambda a: bool(abs(1-a)),predictions))
				if i==0:
					always_xs=set(remaining_set.loc[false_predictions].index.tolist())
				else:
					always_xs=set(always_xs & set(remaining_set.loc[false_predictions].index.tolist()))
			if classifier==classifiers[0]:
				always_xs_tot=always_xs
			else:
				always_xs_tot=set(always_xs_tot & always_xs)
		return(always_xs_tot)
	
	
	
	#plots size of negative set over number of trainings n_vec for given x={"negative", "positive"} and given y_set
	def plot_n_trainings(self,n_vec,save=False,method=np.mean):
		f=plt.figure(figsize=(10,10))
		n_x=list()
		if save==True:
			complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
			complete_set=complete_set.set_index([self.protcol,"metabolite"])
		for i in n_vec:
			always_xs=self.train_to_get_xs(n_trainings=i)
			if save==True:
				complete_set.loc[list(always_xs)].to_csv(self.databases+self.x+"_sampled_"+str(i)+"_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			n_x.append(len(always_xs))
		plt.plot(n_vec,n_x)
		plt.ylabel="size of consistent "+self.x+"s"
		plt.xlabel="number of trained classifiers"
		f.savefig(self.analyses+"number_trainings_for_sampling_"+self.x+"s.png")
		
		
	
	#returns x_set, x={"negative", "positive"} and given y_set
	def save_set(self,tabfile,n_trainings=50000,method=np.mean):
		complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
		complete_set=complete_set.set_index([self.protcol,"metabolite"])
		always_xs=self.train_to_get_xs(n_trainings=n_trainings,method=method)
		complete_set.loc[list(always_xs)].to_csv(tabfile,index=True,header=True,sep="\t")
		return(complete_set.loc[list(always_xs)])
	
	
	#gets intersection between x_set from databases and sampled x_set
	def intersect_xs_from_sampling_and_db(self,method=np.mean):
		x_sampled_set=pd.read_csv(self.databases+self.x+"_sampled_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		x_sampled_set=x_sampled_set.set_index([self.protcol,"metabolite"])
		x_db_set=pd.read_csv(self.databases+self.x+"_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		x_db_set=x_db_set.set_index([self.protcol,"metabolite"])
		intersection=set(set(x_sampled_set.index.tolist()) & set(x_db_set.index.tolist()))
		print("length of intersection:\t"+str(len(intersection)))
		print("length of "+self.x+" sampled set:\t"+str(len(x_sampled_set)))
		print("length of "+self.x+" set from database:\t"+str(len(x_db_set)))




###############################################################
###############################################################
###############################################################

# class constructing subset of training set containing experimentally verified absence of interactions (negative set)

class NegativeSet(XSet):
	
	def __init__(self,simulation=False,overwrite=False,experimental="../experimental_data/",databases="../databases/",analyses="../analyses/",feature="",db_organism=("9606","H.sapiens"),normalized="sum",normalized2=None,proteinwise=False):
		super().__init__(feature=feature,simulation=simulation,overwrite=overwrite,experimental=experimental,databases=databases,analyses=analyses,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
		self.db_organism=db_organism
			
			
	#loads database for mapping of orthologous groups
	def open_cog_db(self):
		try:
			cog_db=open("../databases/COG.mappings.v10.5.txt","rt")
		except FileNotFoundError:
			print("Database for orthologous groups (COG) not found. Please select manually")
			root=Tk()
			root.withdraw()
			cog_db=askopenfilename(initialdir=self.databases,filetypes=[(".gz","*.gz"),("text files","*.txt")])
		return(cog_db)
		
		
	#queries database with COG-mappings with given set of orthologous groups for STRING IDs of OGs in human and saves them to dataframe
	def get_orthologs_in_human(self,cogs_set=None):
		if cogs_set is None:
			print("looking for orthologs between experimental data")
			cogs_set=self.find_orthologs_between_files()
		df=pd.DataFrame(columns=[self.db_organism[1],"OG"])
		cog_db=self.open_cog_db()
		#for every line of database, if orthologeous group present in intersection and humane protein found, put it into dataframe
		match=0 #number of proteins found
		for line in cog_db:
			coloumns=line.split()
			#protein_string_id=coloumns[0]
			#orthologous_group=coloumns[3]
			if coloumns[0][:len(self.db_organism[0])]==self.db_organism[0]:
				if coloumns[3] in cogs_set:
					match=match+1
					if coloumns[3] not in df.OG.values.tolist():
						df=df.append(pd.DataFrame(data={self.db_organism[1]:[coloumns[0]],"OG":[coloumns[3]]}),ignore_index=True)
					else: #if cog already found, append new string id to previous
						df.loc[df.loc[df["OG"]==coloumns[3]].index[0],self.db_organism[1]]=df.loc[df["OG"]==coloumns[3]][self.db_organism[1]].values[0]+";"+coloumns[0] #IX
		df=df.set_index("OG",drop=False)
		print(str(match)+" "+self.db_organism[1]+" proteins found")
		return(df)
		
		
	#for given dataframe with (human) STRING IDs to orthologous groups, read files with OGs and STRING IDs from other organisms and concatenate them
	def collect_all_string_ids(self,df=None):
		if df is None:
			df=self.get_orthologs_in_human()
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			df=pd.read_csv(self.databases+"neg_string_og_"+self.db_organism[1]+".tsv",sep="\t")
			print("collection of STRING IDs from all experiments to corresponding intersection of orthologous groups between all experiments found")
		except (OverwriteException,FileNotFoundError):
			for expfilename in self.expfiles:
				df_org_x=pd.read_csv(self.databases+expfilename[:-5].split("/")[-1]+"_string_cog_from_expdata.tsv",sep="\t")
				df_org_x=df_org_x.set_index("COG",drop=False)
				org_cogs=df_org_x.COG.values.tolist()
				#for every OG...
				match=0
				for cog in set(df.index.tolist()):
					for i in range(len(org_cogs)):
						if cog in org_cogs[i]: #if string of cog in ;-separated string of org_cogs
							try: #if cog occurs only once
								df.loc[cog,expfilename[:-5].split("/")[-1]]=df_org_x.loc[org_cogs[i],"STRING_ID"]
								match=match+len(df_org_x.loc[org_cogs[i],"STRING_ID"].split(";"))
							except: #if cog occurs multiple times
								coglist=df_org_x.loc[org_cogs[i],"STRING_ID"].values.tolist()
								df.loc[cog,expfilename[:-5].split("/")[-1]]=";".join(coglist)
								match=match+len(";".join(coglist).split(";"))
				#print("For "+expfilename[:-5].split("/")[-1]+" "+str(match)+" proteins STRING_IDs found in set of OGs")
			df=df.dropna() # delete OG which are not present in all organisms, since dataframe initialized with OGs over all organisms, not which are present in all organisms
			df.to_csv(self.databases+"neg_string_og_"+self.db_organism[1]+".tsv",header=True,index=False,sep="\t")
		organisms=list(df)
		organisms.remove("OG")
		organisms.remove(self.db_organism[1])
		for org in organisms:
			print("For "+org+" "+str(len(";".join(df[org].values.tolist()).split(";")))+" STRING_IDs found in set of OGs")
		print("collection of STRING IDs from all experiments to corresponding intersection of orthologous groups between all experiments created")
		print(str(len(df))+" orthologous groups present in all organisms")
		return(df)
					
			
	#given a DataFrame with Orthologous Groups and corresponding ;-separated STRING IDs, translate them to UniProt IDs
	def translate_string_df_to_uniprot(self,df_string=None):
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			df_uniprot=pd.read_csv(self.databases+"neg_uniprot_og_"+self.db_organism[1]+".tsv",sep="\t")
			df_uniprot=df_uniprot.set_index("OG",drop=True)
			print("DataFrame with Uniprot IDs per organism found and loaded")
		except (OverwriteException,FileNotFoundError):
			if df_string is None:
				df_string=self.collect_all_string_ids()
			df_uniprot=pd.DataFrame(index=df_string.OG,columns=list(df_string))
			del df_uniprot["OG"]
			df_string=df_string.set_index("OG",drop=True)
			#db=dh.DBHandler("A.thaliana.xlsx") # not nice programming, but creating artificial instance of dh to access translate_proteins function, file not even existent
			''' #all in one query => Connection Error arises
			for coloumn in list(df_string):
				df_uniprot[coloumn]=db.translate_proteins(df_string[coloumn].values.tolist(),in_id_type="STRING_ID",out_id_types=[""])[""].values.tolist()
			'''
			for coloumn in list(df_string):
				match=0
				for i in df_string.index:
					try:
						#uniprot_ids=db.translate_proteins(df_string.loc[i,coloumn].split(";"),in_id_type="STRING_ID",out_id_types=[""])[""].values.tolist()
						uniprot_ids=self.translate_proteins(df_string.loc[i,coloumn].split(";"),in_id_type="STRING_ID",out_id_types=[""])[""].values.tolist()
					except:
						print("ConnectionError")
						t=5
						while True:
							time.sleep(t)
							try:
								#uniprot_ids=db.translate_proteins(df_string.loc[i,coloumn].split(";"),in_id_type="STRING_ID",out_id_types=[""])[""].values.tolist()
								uniprot_ids=self.translate_proteins(df_string.loc[i,coloumn].split(";"),in_id_type="STRING_ID",out_id_types=[""])[""].values.tolist()
								break
							except:
								print("double Error at "+df_string.loc[i,coloumn])
								t=t+5
								#time.sleep(10)
								#uniprot_ids=db.translate_proteins(df_string.ix[i,coloumn].split(";"),in_id_type="STRING_ID",out_id_types=[""])[""].values.tolist()
					try:
						df_uniprot.loc[i,coloumn]=";".join(uniprot_ids)
						match=match+len(uniprot_ids)
					except TypeError: #translation not found, therefore nan in uniprot_ids
						uniprot_ids_string=""
						for uniprot_id in uniprot_ids:
							if pd.isna(uniprot_id)==False:
								if uniprot_ids_string!="":
									uniprot_ids_string=uniprot_ids_string+";"+uniprot_id
									match=match+1
								else:
									uniprot_ids_string=uniprot_id
									match=match+1
						df_uniprot.loc[i,coloumn]=uniprot_ids_string
				print("coloumn "+coloumn+" successfully translated")
				print("For "+coloumn+" "+str(match)+" proteins translated")
			df_uniprot.to_csv(self.databases+"neg_uniprot_og_"+self.db_organism[1]+".tsv",header=True,index=True,sep="\t")
			print("DataFrame with Uniprot IDs per organism created with "+str(len(df_uniprot))+" OGs")
		return(df_uniprot)
			
			
	#query given list of proteins uniprot IDs at PubChem for bioassays and extract compounds (metabolites) that are inactive (=>interaction score =0)
	def find_nonint_cids_in_pubchem(self,uniprot_ids):
		df=pd.DataFrame(columns=["UniProtID","CIDs"])
		#query every protein at PubChem for metabolite CIDs where inactivity (no effect on activity of protein) were determined
		t1=0
		for uniprot_id in uniprot_ids:
			t2=time.time()
			if t2-t1>=60/400 or t1==0: #on PubChem, no more than 400 requests per minute
				t1=time.time()
				query_for_cids="https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/"+uniprot_id+"/cids/TXT?cids_type=inactive"
				resp_cids=requests.post(query_for_cids)
				out_cids=resp_cids.text
				#...if protein and corresponding "inactive" CIDs found on PubChem, eliminate doubled CIDs
				if out_cids[:11]!="Status: 404":
					cids=set(out_cids.split("\n"))
					cids.remove("")
					df=df.append(pd.DataFrame(data={"UniProtID":[uniprot_id],"CIDs":[cids]}))
			else:
				time.sleep(60/400-(t2-t1))
		return(df)
		
		
	#fill dataframe (df_uniprot) with CIDs for each OG
	def find_CIDs(self,df_uniprot=None):
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			df_uniprot=pd.read_csv(self.databases+"neg_uniprot_og_cid_"+self.db_organism[1]+".tsv",sep="\t")
			df_uniprot=df_uniprot.set_index("OG",drop=True)
		except (OverwriteException,FileNotFoundError):
			if df_uniprot is None:
				print("creating DatFrame with IDs for different organisms of orthologous groups to query them for CIDs")
				df_uniprot=self.translate_string_df_to_uniprot()
			df_uniprot=df_uniprot.dropna()
			for og in df_uniprot.index:
				df_og=self.find_nonint_cids_in_pubchem(df_uniprot.loc[og,self.db_organism[1]].split(";"))
				cids_list=df_og.CIDs.values.tolist()
				if cids_list==[]: # if no inactive CIDs found, delete row
					df_uniprot=df_uniprot.drop(og)
				else:
					cids=set(cids_list[0])
					if len(cids)>1: # if more than one CID found, build intersection between all proteins (only CIDs considered, which are present over all UniProt IDs)
						for i in range(1,len(cids_list)):
							cids=set(cids & cids_list[i])
					if cids!=set(): # if intersection not empty, put it into dataframe
						#make set of cids to one string for dataframe
						cids_string=""
						for cid in cids:
							if cids_string!="":
								cids_string=cids_string+";"+cid
							else:
								cids_string=cid
						df_uniprot.loc[og,"CIDs"]=cids_string
					else: #else delete row for corresponding OG
						df_uniprot=df_uniprot.drop(og)
			df_uniprot=self.padd_df(df_uniprot)
			df_uniprot.to_csv(self.databases+"neg_uniprot_og_cid_"+self.db_organism[1]+".tsv",header=True,index=True,sep="\t")
			print("DataFrame with OGs, UniProtIDs per organism and CIDs in "+self.db_organism[1]+" created")
			print(str(len(df_uniprot))+" OGs with corresponding non-interactive (inactive) bioassaydata found")
		return(df_uniprot)
		
		
		
	#iterates through all CIDs and pads them with zeros in front to a length of 8, so that they match CIDs from database
	def padd_df(self,df_uniprot):
		for i in df_uniprot.index:
			cids_padded=list()
			for cid in df_uniprot.loc[i,"CIDs"].split(";"):
				if len(cid)<8:
					cids_padded.append("0"*(8-len(cid))+cid)
				else:
					cids_padded.append(cid)
			df_uniprot.loc[i,"CIDs"]=";".join(cids_padded)
		return(df_uniprot)
	
	
	def meta_intersect_bioassay_files(self,df_uniprot=None):
		if df_uniprot is None:
			df_uniprot=self.find_CIDs() #has pubchem CIDs padded with zeros
		metabolites=self.meta_intersect_files() #metabolite CIDs regarding STITCH (and Marcin's) database
		metabolite_names=self.translate_CIDs_offline(metabolites)
		#merge CIDs from STITCH and PubChem to one DataFrame for translation
		if "PubChem_CID" not in metabolite_names.columns:
			metabolite_cids_pubchem=self.translate_metabolites_pubchem(metabolite_names["Name"].tolist())
			metabolite_names["PubChem_CID"]=np.empty((len(metabolite_names),0)).tolist() 
			for i in metabolite_names.index:
				if metabolite_names.loc[i,"Name"] in metabolite_cids_pubchem.index: #if online translation worked
					metabolite_names.loc[i,"PubChem_CID"]=metabolite_cids_pubchem.loc[metabolite_names.loc[i,"Name"],"CIDs"]
				else: #if online translation on PubChem did not find a valid result, take the STITCH CID
					metabolite_names.loc[i,"PubChem_CID"]=[i]
			listwise=True
		else:
			#padd PubChem CIDs with zeros
			pubchem_cids_padded=self.padd_cids(metabolite_names["PubChem_CID"])
			metabolite_names["PubChem_CID"]=pubchem_cids_padded
			listwise=False
		pairs=0 #number of protein metabolite pairs
		organisms=list(df_uniprot)
		organisms.remove("CIDs")
		organisms.remove(self.db_organism[1])
		#exclude all metabolites in negative set from bioassay data which are not present in above metabolite intersection from experimental data
		for og in df_uniprot.index:
			if listwise:
				df_uniprot.loc[og,"CIDs"]=set(set(df_uniprot.loc[og,"CIDs"].split(";")) & set(itertools.chain.from_iterable(metabolite_names["PubChem_CID"].tolist()))) # set(itertools.chain.from_iterable(metabolite_cids_pubchem["CIDs"].tolist()))) #better: intersection between online PubChem CID translation and bioassay PubChem CIDs
			else:
				df_uniprot.loc[og,"CIDs"]=set(set(df_uniprot.loc[og,"CIDs"].split(";")) & set(metabolite_names["PubChem_CID"].tolist())) # set(itertools.chain.from_iterable(metabolite_cids_pubchem["CIDs"].tolist()))) #better: intersection between online PubChem CID translation and bioassay PubChem CIDs
			#retranslate CIDs from PubChem CIDs to STITCH CIDs
			stitch_cids=set()
			for cid in df_uniprot.loc[og,"CIDs"]:
				for cid_list in metabolite_names["PubChem_CID"]:
					if cid in cid_list:
						stitch_cids.add(metabolite_names.index[metabolite_names["PubChem_CID"].tolist().index(cid_list)])
			#if intersection empty, try to find if there exists discrepancy to CID used in STITCH and if so, take it, else delete row in dataframe for this OG
			if stitch_cids==set(): #for case of intersection between STITCH and PubChem CID, change here to df_uniprot.loc[og,"CIDs"]
				df_uniprot=df_uniprot.drop(og)
			else:
				df_uniprot.loc[og,"CIDs"]=";".join(stitch_cids) #use translated Stitch CIDs for set
				#df_uniprot.loc[og,"CIDs"]=stitch_cids
				for org in organisms: #count number of protein-metabolite interaction pairs
					#pairs=pairs+len(df_uniprot.loc[og,org].split(";"))*len(df_uniprot.loc[og,"CIDs"])
					pairs=pairs+len(df_uniprot.loc[og,org].split(";"))*len(stitch_cids)
		print(str(len(df_uniprot))+" OGs with metabolites in set")
		print(str(len(set(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),df_uniprot.CIDs.tolist()))))))+" metabolites present in set")
		#print(str(len(set(itertools.chain.from_iterable(df_uniprot.CIDs.values.tolist()))))+" metabolites present in set")
		print(str(len(list(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),df_uniprot.CIDs.tolist()))))))+" OG-metabolite interaction pairs present in set")
		#print(str(len(list(itertools.chain.from_iterable(df_uniprot.CIDs.values.tolist()))))+" OG-metabolite interaction pairs present in set")
		print(str(pairs)+" protein-metabolite interaction-pairs in set")
		return(df_uniprot)
		
		
		
	
	#trim negative set to proteinwise (UniProt-CID pairs) set, and exclude highly correlated profiles
	def proteinwise_x_set(self,x_set,method=np.mean,corr_thresh=0.5):
		protwise_x_set=pd.DataFrame(columns=["protein","CIDs"]).set_index("protein",drop=False)
		intersecting_prots=self.prot_intersect_files()
		exp_cols=list() #columns with experiments
		filehandlers=list() #list of filehandler instances
		for expfilename in self.expfiles:
			exp_cols.append(expfilename[:-5].split("/")[-1])
			filehandlers.append(FileHandler(expfilename,experimental=self.experimental,databases=self.databases))
		for og in x_set.index:
			proteins=set(set(";".join(x_set.loc[og,exp_cols].tolist()).split(";")) & intersecting_prots) #intersection of proteins to corresponding OG to intersecting proteins
			if len(proteins)>1: #if more than one protein, calculate profile correlation 
				low_corr_pairs=list() #list of protein pairs with low profile correlation
				for expdata in filehandlers:
					corr_df=expdata.profiles_correlation_proteins(proteins,profiles="raw",method=method)
					for i in range(len(corr_df.index)):
						for j in range(i,len(corr_df.columns)):
							if corr_df.iloc[i,j]<corr_thresh:
								low_corr_pairs.append((corr_df.index[i],corr_df.columns[j]))
				#filter to pairs having correlation below corr_thresh in all expdatas
				c=Counter(map(lambda x: str(x),low_corr_pairs))
				overall_low_corr_pairs=list(Counter(e for e in c.elements() if c[e]==len(self.expfiles)).elements())
				if len(overall_low_corr_pairs)>0:
					#create graph
					g=Graph.TupleList(list(map(lambda x: eval(x),overall_low_corr_pairs)),directed=False)
					largest_clique=g.largest_cliques()[0] #first largest clique
					proteins=list(map(lambda x: g.vs[x]["name"],largest_clique))
					for protein in proteins:
						protwise_x_set.loc[protein,"CIDs"]=x_set.loc[og,"CIDs"]
				else: #just take first uniprot ID
					protwise_x_set.loc[list(proteins)[0],"CIDs"]=x_set.loc[og,"CIDs"]
			elif len(proteins)==1:
				protwise_x_set.loc[list(proteins)[0],"CIDs"]=x_set.loc[og,"CIDs"]
		protwise_x_set["protein"]=protwise_x_set.index
		return(protwise_x_set)
	
	
	
	
	#saves dataframe to file
	def save_set(self,tabfile=None,df=None,method=None): #df is df_uniprot
		if tabfile is None:
			tabfile=self.databases+"negative_set_"+self.db_organism[1]+".tsv"
		if df is None:
			df=self.meta_intersect_bioassay_files()
		#converts sets in CIDs to string, not necessary if ";".join(stitch_cids) in self.meta_intersect_bioassay_files
		'''
		for i in df.index:
			cids_string=""
			for j in df.loc[i,"CIDs"]:
				if cids_string=="":
					cids_string=j
				else:
					cids_string=cids_string+";"+j
			df.loc[i,"CIDs"]=cids_string
		'''
		df.to_csv(tabfile,header=True,index=True,sep="\t")
		return(df)




###############################################################
###############################################################
###############################################################

# class for reading experimental data, interacting with database to find positive subset (+) of training set from STITCH database



class PositiveSet(XSet):
	
	
	#find orthologous groups for all experimental data files (different species)
	def get_orthologs_between_species(self): 
		prot_all_orgs=list()
		#for every excel sheet with experimental data (=organism)...
		for i in range(len(self.expfiles)):
			expfilename=self.expfiles[i]
			org=DBFileInteractions(expfilename,simulation=self.simulation,experimental=self.experimental,databases=self.databases)
			# ...find orthologs and save them as tsv file
			if self.simulation==False:
				pos=org.find_positive_candidates()
				org.extract_and_save_orthologs(pos,self.databases+org.database.organism+"_pos_string_orthologs.tsv",overwrite=self.overwrite,include_metas=True)
				df_org_x=pd.read_csv(self.databases+org.database.organism+"_pos_string_orthologs.tsv",sep="\t")
			else:
				pos=org.find_positive_candidates()
				org.extract_and_save_orthologs(pos,"../simulation_data/"+org.database.organism+"_pos_string_orthologs.tsv",overwrite=self.overwrite,include_metas=True)
				df_org_x=pd.read_csv("../simulation_data/"+org.database.organism+"_pos_string_orthologs.tsv",sep="\t")
			# ... put them into a set and determine intersect to previous ones
			'''
			#only if all COG per protein match
			if i==0:
				orthologs=set(df_org_x.COG)
			else:
				self.intersect(orthologs,set(df_org_x.COG))
			'''
			#for every COG on its own 
			list_org_x=df_org_x.COG.values.tolist()
			for listelement in list_org_x:
				if not pd.isna(listelement):
					cogs=str(listelement).split(";")
					for element in cogs:
						prot_all_orgs.append(element)
		orthologs=set(prot_all_orgs)
		return(orthologs)
	
	
	#get metabolites for orthologeous proteins, maybe more efficient to include get_orthologs_between_species, since the databases have to be openend once, but maybe faster if intersection between species gets smaller
	#creation of the gold standard: dataframe with orthologeous proteins and interacting metabolites
	def get_metabolites_for_all_orthologs(self,orthologs=None):
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			df_res=pd.read_csv(self.databases+"positive_set_raw.tsv",sep="\t")
			df_res=df_res.set_index(list(df_res)[0],drop=True) #set OG as index, coloumn "Unnamed: 0"
			print("dataframe with STRING IDs and not filtered metabolite IDs found and loaded")
		except (OverwriteException,FileNotFoundError):
			if orthologs==None:
				orthologs=self.get_orthologs_between_species()
			df_res=pd.DataFrame(index=orthologs,columns=["CID"])
			#for every excel sheet with experimental data import corresponding string_orthologs dataframe and merge them according to given orthologeous proteins (orthologs) 
			for i in range(len(self.expfiles)):
				expfilename=self.expfiles[i]
				organism=expfilename[:-5].split("/")[-1] #.split("_")[0]
				if self.simulation==False:
					df_org_x=pd.read_csv(self.databases+organism+"_pos_string_orthologs.tsv",sep="\t")
				else:
					df_org_x=pd.read_csv("../simulation_data/"+organism+"_pos_string_orthologs.tsv",sep="\t")
				df_org_x=df_org_x.dropna()
				df_org_x=df_org_x.set_index("COG",drop=True)
				for cogs in df_org_x.index: # since cogs are concatenated by ";"
					cogs_list=cogs.split(";")
					for cog in cogs_list:
						if cog in df_res.index: 
							try: #if one COG occurs only once
								df_res.loc[cog,"CID_"+organism]=df_org_x.loc[cogs,"CID"] 
								df_res.loc[cog,"STRING_ID_"+organism]=df_org_x.loc[cogs,"STRING_ID"] 
							except:
								#since there may be one COG represented several times, put them all to a list
								coglist=df_org_x.loc[cogs,"CID"].values.tolist() 
								stringidlist=df_org_x.loc[cogs,"STRING_ID"].values.tolist() 
								#filter from this list those CIDs which are present in every COG
								for j in range(len(coglist)):
									if j==0:
										cids=set(coglist[j].split(";"))
									else:
										cids=set(cids & set(coglist[j].split(";"))) #intersection
								#append STRING_IDs to dataframe, set to exclude double occurences
								stringids=set()
								for k in range(len(stringidlist)):
									stringids.add(stringidlist[k])
								df_res.loc[cog,"CID_"+organism]=";".join(list(cids))
								df_res.loc[cog,"STRING_ID_"+organism]=";".join(list(stringids))
							#get intersect in CID for all organisms as intersect from CID and CID_+organism
							if i==0: #pd.isna(df_res.ix[cog,"CID"]):
								df_res.loc[cog,"CID"]=df_res.loc[cog,"CID_"+organism]
							else:
								df_res.loc[cog,"CID"]=";".join(list(set(set(df_res.loc[cog,"CID_"+organism].split(";")) & set(str(df_res.loc[cog,"CID"]).split(";")))))
				#remove all OG where current organism can not contribute with data
				df_res=df_res[pd.isna(df_res["CID_"+organism])==False]
			#remove COG without metabolite CIDs
			df_res=df_res[df_res.CID!=""]
			#df_res=df_res[pd.isna(df_res.CID)==False] #to remove NaN
			df_res.to_csv(self.databases+"positive_set_raw.tsv",header=True,index=True,sep="\t")
			print("dataframe with STRING IDs and not filtered metabolite IDs created")
		return(df_res)
		
		
	
	
	#get proteinwise positve set raw (not filtered to stereospecific CIDs
	def protwise_pos_set(self):
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing positive_set_raw requested")
			df_res=pd.read_csv(self.databases+"positive_set_raw.tsv",sep="\t")
			df_res=df_res.set_index(list(df_res)[0],drop=True) #set protein as index, coloumn "Unnamed: 0"
			print("dataframe with STRING IDs and not filtered metabolite IDs found and loaded")
		except FileNotFoundError:
			#load intersecting proteins (UniProt IDs) and intersection with Stitch DB
			proteins=self.prot_intersect_files() #from files
			organism=self.expfiles[0].split("/")[-1].split("_")[0]+".xlsx" #S.cerevisiae or A.thaliana
			db=DBHandler(organism,score_cutoff=None,overwrite=self.overwrite,simulation=self.simulation,experimental=self.experimental,databases=self.databases,analyses=self.analyses)
			#db=DBHandler(organism,score_cutoff=None,overwrite=ps.overwrite,simulation=ps.simulation,experimental=ps.experimental,databases=ps.databases,analyses=ps.analyses)
			positive_candidates=set(db.get_proteins() & proteins)
			stitch_db=pd.read_csv(db.file,sep="\t").set_index("protein")
			stitch_proteins=db.proteins.set_index("ID")
			#initialize DataFrame with proteins and corresponding CIDs
			df_res=pd.DataFrame(index=positive_candidates,columns=["STRING_ID","CID"])
			#for every protein, find CIDs from Stitch Database
			for uniprot in df_res.index:
				df_res.loc[uniprot,"STRING_ID"]=stitch_proteins.loc[uniprot,"STRING_ID"]
				if type(stitch_db.loc[stitch_proteins.loc[uniprot,"STRING_ID"],"chemical"])==str: #one cid match
					df_res.loc[uniprot,"CID"]=stitch_db.loc[stitch_proteins.loc[uniprot,"STRING_ID"],"chemical"]
				else: #several cid matches
					df_res.loc[uniprot,"CID"]=";".join(stitch_db.loc[stitch_proteins.loc[uniprot,"STRING_ID"],"chemical"])
			df_res.to_csv(self.databases+"positive_set_raw.tsv",sep="\t",index=True,header=True)
		return(df_res)
		
		
	
	
	#filters for CIDs or CID0, which are stereo-specific compounds, CIDm/CID1 are with merged stereo-isomers and for metabolites in file intersection
	def filter_stereospecific_cids(self,metabolites=None,df=None):
		try:
			df=pd.read_csv(self.databases+"positive_set_string_all_cids.tsv",sep="\t")
		except FileNotFoundError:
			if df is None:
				if self.proteinwise==True:
					df=self.protwise_pos_set()
				else:
					df=self.get_metabolites_for_all_orthologs()
			if metabolites is None:
				metabolites=self.meta_intersect_files()
			metabolites_padded=self.padd_cids(metabolites)
			organisms=list(filter(lambda org: org[:9]=="STRING_ID", list(df)))
			pairs=0
			for i in df.index:
				filtered_cids=""
				cids=df.loc[i,"CID"].split(";")
				for cid in cids:
					if cid[3]=="s" or cid[3]=="0":
						if cid[4:] in metabolites_padded:  #only take CID number, not "CIDs-number"
							if filtered_cids=="":
								filtered_cids=cid[4:]
							else:
								filtered_cids=filtered_cids+";"+cid[4:]
				df.loc[i,"CID"]=filtered_cids
				if self.proteinwise==False:
					#count number of PMI-pairs
					prots=set()
					for org in organisms:
						prots=set(prots | set(df.loc[i,org].split(";")))###
					pairs=pairs+len(prots)*len(df.loc[i,"CID"].split(";"))
				else:
					pairs=pairs+len(df.loc[i,"CID"].split(";"))
			df=df[df.CID!=""]
			df.to_csv(self.databases+"positive_set_string_all_cids.tsv",header=True,index=True,sep="\t")
			print(str(len(df))+" OGs/proteins in positive set")
			print(str(len(set(itertools.chain.from_iterable(list(map(lambda cids: cids.split(";"),df.CID.values.tolist()))))))+" metabolites in positive set")
			if self.proteinwise==False:
				print(str(len(df.CID.values.tolist()))+" OG-metabolite interactions in positive set")
			print(str(pairs)+" protein-metabolite interaction-pairs in positive set")
		return(df)
	
	
	
	### GENERAL ###
	#trim BioLiP database to organisms in experimental data and translate metabolites
	def trim_biolip_to_orglist(self,organisms=["Arabidopsis thaliana","Saccharomyces cerevisiae","Escherichia coli"]):
		biolip_db=pd.read_csv("../databases/BindingDB_All_interactions_nonan.tsv",sep="\t")
		for org in organisms:
			biolip_db_org=biolip_db[biolip_db["Target Source Organism According to Curator or DataSource"]==org]
			#translate pubchem cids to metabolite names and then to stitch cids
			pubchem_to_names=self.translate_cids_pubchem(biolip_db_org["PubChem CID"].tolist())
			names_to_stitch=self.translate_metabolites_offline(pubchem_to_names["Name"].tolist())
			translation_df=pubchem_to-names.reset_index().set_index("Name").join(names_to_stitch)
			#stopped here, because most of ligands from biolip for arabidopsis were not biologically relevant
			
			
		
		
		
	# translate positive set to UniProtIDs and saves it
	def save_set(self,tabfile=None,df=None,method=None):
		df_res=pd.DataFrame()
		if tabfile is None:
			tabfile=self.databases+"positive_set.tsv"
		if df is None:
			df=self.filter_stereospecific_cids()
		if self.proteinwise==True:
			df_res=df.drop("STRING_ID",axis=1)
			df_res.reset_index(inplace=True)
			df_res.columns=["protein","CIDs"]
			df_res.to_csv(tabfile,header=True,index=False,sep="\t")
		else:
			organisms=list(filter(lambda org: org[:9]=="STRING_ID", list(df)))
			#organisms=list(map(lambda x: x.split("/")[-1][:-5],self.expfiles))
			for i in range(len(organisms)):
				org=organisms[i]#[10:]
				db=pd.read_csv("../databases/"+org.split("_")[0]+"_ID_collection.tsv",sep="\t") ### new: .split("_")[0], to make only one per organism
				db=db.set_index("STRING_ID",drop=True)
				for ind in df.index:
					string_ids=df.loc[ind,organisms[i]].split(";")
					uniprot_ids=""
					for string_id in string_ids:
						if uniprot_ids=="":
							uniprot_ids=db.loc[string_id,"ID"]
						else:
							uniprot_ids=uniprot_ids+";"+db.loc[string_id,"ID"]
					df_res.loc[ind,org]=uniprot_ids
					if i==0:
						df_res.loc[ind,"CIDs"]=df.loc[ind,"CID"]
			df_res.to_csv(tabfile,header=True,index=True,sep="\t")
		return(df_res)
		
		



###############################################################
###############################################################
###############################################################

# class to construct positive subset like negative set from PubChem


class PositiveSet_PubChem(NegativeSet):
	
	#query given list of proteins uniprot IDs at PubChem for bioassays and extract compounds (metabolites) that are inactive (=>interaction score =0)
	def find_int_cids_in_pubchem(self,uniprot_ids):
		df=pd.DataFrame(columns=["UniProtID","CIDs"])
		#query every protein at PubChem for metabolite CIDs where inactivity (no effect on activity of protein) were determined
		t1=0
		for uniprot_id in uniprot_ids:
			t2=time.time()
			if t2-t1>=60/400 or t1==0: #on PubChem, no more than 400 requests per minute
				t1=time.time()
				query_for_cids="https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/"+uniprot_id+"/cids/TXT?cids_type=active"
				resp_cids=requests.post(query_for_cids)
				out_cids=resp_cids.text
				#...if protein and corresponding "inactive" CIDs found on PubChem, eliminate doubled CIDs
				if out_cids[:11]!="Status: 404":
					cids=set(out_cids.split("\n"))
					cids.remove("")
					df=df.append(pd.DataFrame(data={"UniProtID":[uniprot_id],"CIDs":[cids]}))
			else:
				time.sleep(60/400-(t2-t1))
		return(df)
	
	
	
	#fill dataframe (df_uniprot) with CIDs for each OG
	def find_CIDs(self,df_uniprot=None):
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			df_uniprot=pd.read_csv(self.databases+"pos_uniprot_og_cid_"+self.db_organism[1]+".tsv",sep="\t")
			df_uniprot=df_uniprot.set_index("OG",drop=True)
		except (OverwriteException,FileNotFoundError):
			if df_uniprot is None:
				print("creating DatFrame with IDs for different organisms of orthologous groups to query them for CIDs")
				df_uniprot=self.translate_string_df_to_uniprot()
			df_uniprot=df_uniprot.dropna()
			for og in df_uniprot.index:
				df_og=self.find_int_cids_in_pubchem(df_uniprot.loc[og,self.db_organism[1]].split(";"))
				cids_list=df_og.CIDs.values.tolist()
				if cids_list==[]: # if no active CIDs found, delete row
					df_uniprot=df_uniprot.drop(og)
				else:
					cids=set(cids_list[0])
					if len(cids)>1: # if more than one CID found, build intersection between all proteins (only CIDs considered, which are present over all UniProt IDs)
						for i in range(1,len(cids_list)):
							cids=set(cids & cids_list[i])
					if cids!=set(): # if intersection not empty, put it into dataframe
						#make set of cids to one string for dataframe
						cids_string=""
						for cid in cids:
							if cids_string!="":
								cids_string=cids_string+";"+cid
							else:
								cids_string=cid
						df_uniprot.loc[og,"CIDs"]=cids_string
					else: #else delete row for corresponding OG
						df_uniprot=df_uniprot.drop(og)
			df_uniprot=self.padd_df(df_uniprot)
			df_uniprot.to_csv(self.databases+"pos_uniprot_og_cid_"+self.db_organism[1]+".tsv",header=True,index=True,sep="\t")
			print("DataFrame with OGs, UniProtIDs per organism and CIDs in "+self.db_organism[1]+" created")
			print(str(len(df_uniprot))+" OGs with corresponding interactive (active) bioassaydata found")
		return(df_uniprot)
		
	
	
	#saves dataframe to file
	def save_set(self,tabfile=None,df=None,method=None): #df is df_uniprot
		if tabfile is None:
			tabfile=self.databases+"positive_set_"+self.db_organism[1]+"_pubchem.tsv"
		if df is None:
			df=self.meta_intersect_bioassay_files()
		#converts sets in CIDs to string, not necessary if ";".join(stitch_cids) in self.meta_intersect_bioassay_files
		'''
		for i in df.index:
			cids_string=""
			for j in df.loc[i,"CIDs"]:
				if cids_string=="":
					cids_string=j
				else:
					cids_string=cids_string+";"+j
			df.loc[i,"CIDs"]=cids_string
		'''
		df.to_csv(tabfile,header=True,index=True,sep="\t")
		return(df)
		


###############################################################
###############################################################
###############################################################

# class to construct training set from subsets, and create gold standard
#optionally include sampled subsets if one subset is bigger
#always balances such that positive and negative subsets are of equal length by trimming
#unbalanced vs balanced set describes balancing such that metabolites have to occur in both subsets of the training set



class TrainingSet(XSet):
	def __init__(self,db_orglist,simulation=False,overwrite=False,experimental="../experimental_data/",databases="../databases/",analyses="../analyses/",feature="",methods=[np.mean],normalized="",normalized2=None,proteinwise=False):
		super().__init__(feature=feature,simulation=simulation,overwrite=overwrite,experimental=experimental,databases=databases,analyses=analyses,methods=methods,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
		self.db_orglist=db_orglist
		
	#############################################
	# Loading and merging sets
	#############################################
	
	#merges positive sets from databases to one
	def merge_positive_set_with_pubchems(self):
		db_orgs="_".join(list(map(lambda x: x[1],self.db_orglist)))
		for method in self.methods:
			positive_merged_set_normprofiles=pd.read_csv(self.databases+"positive_set_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",sep="\t")
			positive_merged_set_normprofiles=positive_merged_set_normprofiles.set_index([self.protcol,"metabolite"],drop=True).drop_duplicates(keep="first") #one protein sometimes belongs to several OGs
			for db_org in self.db_orglist:
				positive_set_pubchem_normprofiles=pd.read_csv(self.databases+"positive_set_"+db_org[1]+"_pubchem_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",sep="\t")
				positive_set_pubchem_normprofiles=positive_set_pubchem_normprofiles.set_index([self.protcol,"metabolite"],drop=True)
				previous_length=len(positive_merged_set_normprofiles)
				positive_merged_set_normprofiles=positive_merged_set_normprofiles.append(positive_set_pubchem_normprofiles).drop_duplicates(keep="first")
				print("organism "+str(db_org[1])+" added "+str(len(positive_merged_set_normprofiles)-previous_length)+" new interactions")
			positive_merged_set_normprofiles=positive_merged_set_normprofiles[~positive_merged_set_normprofiles.index.duplicated(keep=False)]
			positive_merged_set_normprofiles.dropna().to_csv(self.databases+"positive_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			positive_merged_set_normprofiles=self.refeaturize_x_set(positive_merged_set_normprofiles)
			positive_merged_set_normprofiles.dropna().to_csv(self.databases+"positive_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles"+feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		
		
	#merges negative sets from PubChem for different organisms to one
	def merge_negative_set_with_pubchems(self,feature=None):
		if feature is None:
			feature=self.feature
		db_orgs="_".join(list(map(lambda x: x[1],self.db_orglist)))
		for method in self.methods:
			negative_merged_set_normprofiles=pd.DataFrame()
			for db_org in self.db_orglist:
				negative_set_normprofiles=pd.read_csv(self.databases+"negative_set_"+db_org[1]+"_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",sep="\t")
				negative_set_normprofiles=negative_set_normprofiles.set_index([self.protcol,"metabolite"],drop=True)
				previous_length=len(negative_merged_set_normprofiles)
				negative_merged_set_normprofiles=negative_merged_set_normprofiles.append(negative_set_normprofiles).drop_duplicates() #one protein sometimes belongs to several OGs
				print("organism "+str(db_org[1])+" added "+str(len(negative_merged_set_normprofiles)-previous_length)+" new interactions")
			negative_merged_set_normprofiles=negative_merged_set_normprofiles[~negative_merged_set_normprofiles.index.duplicated(keep=False)]
			negative_merged_set_normprofiles.dropna().to_csv(self.databases+"negative_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			negative_merged_set_normprofiles=self.refeaturize_x_set(negative_merged_set_normprofiles)
			negative_merged_set_normprofiles.dropna().to_csv(self.databases+"negative_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles"+feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		
		
	
	#load xset and merge over different organisms, x={"positive","negative"}
	def load_merged_x_set(self,method,x,feature=None,):
		if feature is None:
			feature=self.feature
		db_orgs="_".join(list(map(lambda d: d[1],self.db_orglist)))
		try:
			x_merged_set=pd.read_csv(self.databases+x+"_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles"+feature+"_"+method.__name__+".tsv",sep="\t")
		except FileNotFoundError:
			eval("self.merge_"+x+"_set_with_pubchems()")
			x_merged_set=pd.read_csv(self.databases+x+"_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles"+feature+"_"+method.__name__+".tsv",sep="\t")
		x_merged_set=x_merged_set.set_index([self.protcol,"metabolite"],drop=True)
		return(x_merged_set)
	
	
	
	#load xset, in this method x_set is loaded and different expansions/trimmings are applied, which can be switched boolean
	# metabolite_extension: to extend with similar metabolites
	# manual_selection: to trim by manual selection of co-eluting pairs, please provide confidence level
	def load_x_set(self,method,x,metabolite_extension=False,manual_selection=False,confidence=""):
		x_set=self.load_merged_x_set(method=method,x=x)
		appendix=""
		if manual_selection:
			db_orgs="_".join(list(map(lambda d: d[1],self.db_orglist)))
			try:
				x_set=pd.read_csv(self.databases+x+"_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_trimmed_"+confidence+"_confidence_"+method.__name__+".tsv",sep="\t")
				x_set=x_set.set_index([self.protcol,"metabolite"])
			except FileNotFoundError:
				x_set=self.trim_x_set_by_manual_selection(x=x,x_set=x_set,confidence=confidence,method=method)
				x_set.to_csv(self.databases+x+"_merged_"+db_orgs+"_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_trimmed_"+confidence+"_confidence_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		if metabolite_extension:
			appendix+="_meta_extended"
			x_set=self.expand_xset_with_similar_metabolites(xset=x_set,method=method)
		return(x_set)
	
	
	
	#merge positive and negative set together without balancing
	#extend_with_sampling=True if you want to take whole set and have as many negatives as positives and filled them with instances from sampled set
	def get_unbalanced_training_set(self,method,extend_with_sampling=True,metabolite_extension=False,manual_selection=False,confidence=""):
		positive_set=self.load_x_set(x="positive",method=method,metabolite_extension=metabolite_extension,manual_selection=manual_selection,confidence=confidence)
		negative_set=self.load_x_set(x="negative",method=method,metabolite_extension=metabolite_extension,manual_selection=False,confidence=confidence)
		appendix=""
		if metabolite_extension:
			appendix+="_meta_extended"
		if manual_selection:
			appendix=appendix+"_"+confidence+"_confidence"
		if extend_with_sampling:
			#if one x_set is bigger than the other, fill with data from sampled sets
			if len(positive_set)>len(negative_set): 
				nss=XSet_sampled(experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,x="negative",y_set=positive_set,balanced=False,normalized=self.normalized,normalized2=self.normalized2,proteinwise=self.proteinwise)
				negative_sampled_set=nss.load_set(nss.databases+nss.x+"_sampled_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+appendix+"_"+method.__name__+".tsv",method=method)
				positive_set=positive_set.assign(interaction=True) #positive_set["interaction"]=True
				negative_set=negative_set.assign(interaction=False) #negative_set["interaction"]=False
				diff=len(positive_set)-len(negative_set)
				if len(negative_sampled_set)>0:
					if len(negative_sampled_set)>=diff: #if negative sampled set bigger/equal than difference between positive and negative set
						training_set=positive_set.append(negative_set)
					else:
						diff=len(negative_sampled_set)
						training_set=negative_set.append(positive_set.iloc[sample(range(len(positive_set)),len(negative_set)+len(negative_sampled_set))])
					negative_sampled_set=negative_sampled_set.reset_index().set_index([self.protcol,"metabolite"]) # reset index to multi-index
					negative_append_set=negative_sampled_set.iloc[sample(range(len(negative_sampled_set)),diff)]
					negative_append_set=negative_append_set.assign(interaction=False) #negative_append_set["interaction"]=False
					training_set=training_set.append(negative_append_set)
				else:
					training_set=negative_set.append(positive_set.iloc[sample(range(len(positive_set)),len(negative_set))])
			elif len(negative_set)>len(positive_set):
				pss=XSet_sampled(experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,x="positive",y_set=negative_set,balanced=False,normalized=self.normalized,normalized2=self.normalized2,proteinwise=self.proteinwise)
				positive_sampled_set=pss.load_set(pss.databases+pss.x+"_sampled_set_"+pss.normalized+pss.normalized2+"normprofiles"+pss.feature+appendix+"_"+method.__name__+".tsv",method=method)
				positive_set=positive_set.assign(interaction=True) #positive_set["interaction"]=True
				negative_set=negative_set.assign(interaction=False) #negative_set["interaction"]=False
				diff=len(negative_set)-len(positive_set) #how many sampled can be used to extend
				if len(positive_sampled_set)>0:
					if len(positive_sampled_set)>=diff: #if positive sampled set bigger/equal than difference between positive and negative set
						training_set=positive_set.append(negative_set)
					else: #elif positive sampled set smaller than difference between positive and negative set
						diff=len(positive_sampled_set)
						training_set=positive_set.append(negative_set.iloc[sample(range(len(negative_set)),len(positive_set)+len(positive_sampled_set))])
					positive_sampled_set=positive_sampled_set.reset_index().set_index([self.protcol,"metabolite"]) # reset index to multi-index
					positive_append_set=positive_sampled_set.iloc[sample(range(len(positive_sampled_set)),diff)]
					positive_append_set=positive_append_set.assign(interaction=True) #positive_append_set["interaction"]=True
					training_set=training_set.append(positive_append_set)
				else:
					training_set=positive_set.append(negative_set.iloc[sample(range(len(negative_set)),len(positive_set))])
			if self.normalized=="":
				training_set.dropna().to_csv(self.analyses+"training_set"+self.feature+"_unbalanced_sampled"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
				self.split_training_set_into_xsets(training_set,approach=self.feature+"_unbalanced_sampled"+appendix,method=method)
			else:
				training_set.dropna().to_csv(self.analyses+"training_set"+self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_unbalanced_sampled"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
				self.split_training_set_into_xsets(training_set,approach=self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_unbalanced_sampled"+appendix,method=method)
		else:
			neg_indices=sample(range(len(negative_set)),floor(min(len(negative_set),len(positive_set))))
			pos_indices=sample(range(len(positive_set)),floor(min(len(negative_set),len(positive_set))))
			positive_set=positive_set.assign(interaction=True) #positive_set["interaction"]=True
			negative_set=negative_set.assign(interaction=False) #negative_set["interaction"]=False
			training_set=positive_set.iloc[pos_indices,:].append(negative_set.iloc[neg_indices,:])
			training_set.dropna().to_csv(self.analyses+"training_set"+self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_unbalanced"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			self.split_training_set_into_xsets(training_set,approach=self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_unbalanced"+appendix,method=method)
		return(training_set)
		
		
		
	#############################################
	# Balancing sets
	#############################################
	
	#balances training set, balance_to={"positive","negative","both"}, means taking as many interaction pairs as in balance_to set, if both takes minimum of both sets for each metabolite
	def balance_training_set(self,positive_set,negative_set,balance_to="both"):
		training_set=pd.DataFrame()
		#reduce positive and negative set to intersecting metabolites
		intersect_metas=list(set(set(positive_set.index.get_level_values(1))&set(negative_set.index.get_level_values(1))))
		positive_set=positive_set.loc[list(map(lambda i: i in intersect_metas,positive_set.index.get_level_values(1)))]
		positive_set=positive_set.assign(interaction=True) #positive_set["interaction"]=True
		negative_set=negative_set.loc[list(map(lambda i: i in intersect_metas,negative_set.index.get_level_values(1)))]
		negative_set=negative_set.assign(interaction=False) #negative_set["interaction"]=False
		#for every metabolite in intersection, sample as many indices containing it from positive and negative set as are in xset with less occurences of this metabolite
		for m in intersect_metas:
			posset_m=positive_set.loc[list(map(lambda i: i==m,positive_set.index.get_level_values(1)))]
			negset_m=negative_set.loc[list(map(lambda i: i==m,negative_set.index.get_level_values(1)))]
			if balance_to=="both":
				samplerange=min(len(posset_m),len(negset_m))
			elif balance_to=="positive":
				if len(negset_m)>=len(posset_m):
					samplerange=len(posset_m)
				else:
					samplerange=len(negset_m)
			elif balance_to=="negative":
				if len(posset_m)>=len(negset_m):
					samplerange=len(negset_m)
				else:
					samplerange=len(posset_m)
			training_set=training_set.append(posset_m.iloc[sample(range(len(posset_m)),samplerange)])
			training_set=training_set.append(negset_m.iloc[sample(range(len(negset_m)),samplerange)])
		return(training_set)
		
	
	#balances training_set from databases without sampling
	def get_balanced_db_set(self,method,metabolite_extension=False,manual_selection=False,confidence=""):
		positive_set=self.load_x_set(x="positive",method=method,metabolite_extension=metabolite_extension,manual_selection=manual_selection,confidence=confidence)
		negative_set=self.load_x_set(x="negative",method=method,metabolite_extension=metabolite_extension,manual_selection=False,confidence=confidence)
		training_set_from_dbs=self.balance_training_set(positive_set,negative_set,balance_to="both")
		#delete pairs which are present doubled (e.g. found in positive and negative set)
		training_set_from_dbs=training_set_from_dbs[~training_set_from_dbs.index.duplicated(keep=False)].drop_duplicates()
		appendix=""
		if metabolite_extension:
			appendix+="_meta_extended"
		if manual_selection:
			appendix=appendix+"_"+confidence+"_confidence"
		if self.normalized=="":
			training_set_from_dbs.dropna().to_csv(self.analyses+"training_set"+self.feature+"_balanced_from_dbs"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			self.split_training_set_into_xsets(training_set_from_dbs,approach=self.feature+"_balanced_from_dbs"+appendix,method=method)
		else:
			training_set_from_dbs.dropna().to_csv(self.analyses+"training_set"+self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_balanced_from_dbs"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			self.split_training_set_into_xsets(training_set_from_dbs,approach=self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_balanced_from_dbs"+appendix,method=method)
		return(training_set_from_dbs)
	
	
	def get_balanced_sampled_set(self,method,metabolite_extension=False,manual_selection=False,confidence=""):
		#load positive and negative set
		positive_set=self.load_x_set(x="positive",method=method,metabolite_extension=metabolite_extension,manual_selection=manual_selection,confidence=confidence)
		negative_set=self.load_x_set(x="negative",method=method,metabolite_extension=metabolite_extension,manual_selection=False,confidence=confidence)
		appendix=""
		if metabolite_extension:
			appendix+="_meta_extended"
		if manual_selection:
			appendix=appendix+"_"+confidence+"_confidence"
		#get balanced negative sampled set
		nss=XSet_sampled(experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,x="negative",y_set=positive_set,normalized=self.normalized,normalized2=self.normalized2,proteinwise=self.proteinwise)
		negative_sampled_set=nss.load_set(nss.databases+nss.x+"_balanced_sampled_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+appendix+"_"+method.__name__+".tsv",method=method)
		negative_sampled_set=negative_sampled_set.reset_index().set_index([self.protcol,"metabolite"]) # reset index to multi-index
		training_set_neg_sampled=self.balance_training_set(positive_set,negative_sampled_set,balance_to="positive") 
		#get balanced positive sampled set
		pss=XSet_sampled(experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,x="positive",y_set=negative_set,normalized=self.normalized,normalized2=self.normalized2,proteinwise=self.proteinwise)
		positive_sampled_set=pss.load_set(pss.databases+pss.x+"_balanced_sampled_set_"+self.normalized+self.normalized2+"normprofiles"+pss.feature+appendix+"_"+method.__name__+".tsv",method=method)
		positive_sampled_set=positive_sampled_set.reset_index().set_index([self.protcol,"metabolite"]) # reset index to multi-index
		training_set_pos_sampled=self.balance_training_set(positive_sampled_set,negative_set,balance_to="negative") 
		#get balanced double sampled set (union of sampled sets), ##### length=2*len(positive_set)+2*len(negative_set)-3*len(training_set_from_dbs)
		training_set_sampled=training_set_pos_sampled.append(training_set_neg_sampled)
		#delete pairs which are present doubled (e.g. found in positive and negative set due to sampling)
		if len(training_set_sampled)>0:
			training_set_sampled=training_set_sampled[~training_set_sampled.index.duplicated(keep=False)].drop_duplicates()
		#combine balanced sampled set with balanced from dbs set
		training_set=self.get_balanced_db_set(method=method,metabolite_extension=metabolite_extension,manual_selection=manual_selection,confidence=confidence)
		training_set=training_set.append(training_set_sampled).drop_duplicates()
		if self.normalized=="":
			training_set.dropna().to_csv(self.analyses+"training_set"+self.feature+"_balanced_sampled"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			self.split_training_set_into_xsets(training_set,approach=self.feature+"_balanced_sampled"+appendix,method=method)
		else:
			training_set.dropna().to_csv(self.analyses+"training_set"+self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_balanced_sampled"+appendix+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
			self.split_training_set_into_xsets(training_set,approach=self.feature+"_"+self.normalized+self.normalized2+"profiles"+"_balanced_sampled"+appendix,method=method)
		return(training_set)
	
	
	#splits given training set into positive and negative set with normprofiles, saves them and constructs corresponding sets with OG-UniProt_IDs-CIDs 
	#approach={"_balanced_from_dbs","_balanced_sampled"}
	def split_training_set_into_xsets(self,training_set,approach,method):
		#split into positive and negative set with normprofiles and save them
		negative_set_normprofiles=training_set[training_set["interaction"]==False]
		negative_set_normprofiles.drop("interaction",axis=1).to_csv(self.databases+"negative_set"+approach+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		positive_set_normprofiles=training_set[training_set["interaction"]==True]
		positive_set_normprofiles.drop("interaction",axis=1).to_csv(self.databases+"positive_set"+approach+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		#construct positive and negative set with OGs, UniProt_IDs, CIDs
		if self.proteinwise==False:
			complete_set=pd.read_csv(self.databases+"complete_set.tsv",sep="\t").set_index("Unnamed: 0")
			positive_uniprot_og_cid=complete_set.loc[positive_set_normprofiles.index.get_level_values(0).unique()]
			negative_uniprot_og_cid=complete_set.loc[negative_set_normprofiles.index.get_level_values(0).unique()]
			for og in positive_uniprot_og_cid.index:
				positive_uniprot_og_cid.loc[og,"CIDs"]=reduce(lambda x,y: str(y)+";"+str(x),positive_set_normprofiles[positive_set_normprofiles.index.get_level_values(0)==og].index.get_level_values(1))
			for og in negative_uniprot_og_cid.index:
				negative_uniprot_og_cid.loc[og,"CIDs"]=reduce(lambda x,y: str(y)+";"+str(x),negative_set_normprofiles[negative_set_normprofiles.index.get_level_values(0)==og].index.get_level_values(1))
			positive_uniprot_og_cid.to_csv(self.databases+"positive_set"+approach+".tsv",header=True,index=True,sep="\t")
			negative_uniprot_og_cid.to_csv(self.databases+"negative_set"+approach+".tsv",header=True,index=True,sep="\t") 
		else:
			positive_uniprot_cid=pd.DataFrame(index=positive_set_normprofiles.index.get_level_values(0).unique(),columns=["protein","CIDs"])
			for prot in positive_set_normprofiles.index:
				positive_uniprot_cid.loc[prot[0],"protein"]=prot[0]
				positive_uniprot_cid.loc[prot[0],"CIDs"]=reduce(lambda x,y: str(y)+";"+str(x),positive_set_normprofiles[positive_set_normprofiles.index.get_level_values(0)==prot[0]].index.get_level_values(1))
			negative_uniprot_cid=pd.DataFrame(index=negative_set_normprofiles.index.get_level_values(0).unique(),columns=["protein","CIDs"])
			for prot in negative_set_normprofiles.index:
				negative_uniprot_cid.loc[prot[0],"protein"]=prot[0]
				negative_uniprot_cid.loc[prot[0],"CIDs"]=reduce(lambda x,y: str(y)+";"+str(x),negative_set_normprofiles[negative_set_normprofiles.index.get_level_values(0)==prot[0]].index.get_level_values(1))
			positive_uniprot_cid.to_csv(self.databases+"positive_set"+approach+".tsv",header=True,index=False,sep="\t")
			negative_uniprot_cid.to_csv(self.databases+"negative_set"+approach+".tsv",header=True,index=False,sep="\t") 
		return
	
	
	
	#for a given xset and training_set returns xset without pairs which are in training set
	def extract_antibalanced_set(self,xset,training_set):
		intersect_indices=list(set(set(xset.index.tolist())&set(training_set.index.tolist())))
		return(xset.drop(intersect_indices))
	
	
	
	#############################################
	#determine chemically similar metabolites that correlate in elution profiles and add them to corresponding OG
	#############################################
	
	
	#determines correlation between metabolite profiles to detect co-eluting metabolites and returns metabolite pairs with correlation above given threshold
	def metabolite_profile_correlation(self,method=np.mean,threshold=0.9,plot=True):
		complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",sep="\t")
		complete_set=complete_set.set_index([self.protcol,"metabolite"])
		metabolites=complete_set.index.get_level_values(1).unique().tolist()
		metabolite_df=complete_set.loc[complete_set.index.get_level_values(0)[0],]
		corr_df=metabolite_df.transpose().corr()
		corr_df.to_csv(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_all_metas_stitchcids.tsv",sep="\t",header=True,index=True)
		#plot heatmap of correlations
		if plot==True:
			f=plt.figure()
			sns.heatmap(corr_df.where(np.tril(np.ones(corr_df.shape),k=-1).astype(np.bool)),annot=True,vmin=-1,vmax=1,linewidths=0.5,cmap="hot_r")
			plt.show()
			f.savefig(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_map.png")
		#triangle=corr_df.where(np.tril(np.ones(corr_df.shape),k=-1).astype(np.bool))[corr_df>=threshold]
		meta_corr=pd.DataFrame(columns=["m1","m2","elution_profile_correlation","Tanimoto_coefficient"])
		metabolites_high_corr=set()
		for m1 in corr_df.index:
			for m2 in corr_df.columns:
				if m1==m2:
					continue
				if corr_df.loc[m1,m2]>=threshold:
					meta_corr=meta_corr.append(pd.DataFrame(columns=["m1","m2","elution_profile_correlation","Tanimoto_coefficient"],data=[[m1,m2,corr_df.loc[m1,m2],np.nan]]),ignore_index=True)
					metabolites_high_corr.add(m1)
					metabolites_high_corr.add(m2)
		#meta_corr["m1"].to_csv(self.analyses+"m1.csv",index=False,header=False)
		#meta_corr["m2"].to_csv(self.analyses+"m2.csv",index=False,header=False)
		#translate cids from stitch to pubchem
		metabolite_names=self.translate_CIDs_from_Stitch_to_PubChem(metabolites_high_corr) ###
		
		#construct dataframe with correlation of elution profiles and Tanimoto index for metabolites
		meta_corr=meta_corr.set_index(["m1","m2"])
		#get rid of B-A pairs for every A-B
		for i in meta_corr.index:
			if (i[0],i[1]) in meta_corr.index:
				meta_corr=meta_corr.drop((i[1],i[0]))
		metabolite_names["PubChem_CID"].to_csv(self.databases+"metabolites_"+self.normalized+self.normalized2+"profiles_high_corr.csv",index=False,header=False)
		metabolite_names=metabolite_names.set_index("Stitch_CID")
		#pd.DataFrame(metabolite_pubchem_high_corr).to_csv(self.databases+"metabolites_high_corr.csv",index=False,header=False)
		return([meta_corr,metabolite_names])
	
	
	#calculates profile correlation for all metabolites and prints list of all PubChem CIDs to file
	def metabolite_profile_correlation_all(self,colormap1,colormap2,method=np.mean):
		try:
			try:
				both_df=pd.read_csv(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_tanimoto_map.tsv",sep="\t").set_index("metabolite")
			except FileNotFoundError:
				tanimoto_df=pd.read_csv(self.analyses+"metabolite_tanimoto_coefficients_all_metas.csv")
				corr_df=pd.read_csv(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_all_metas.tsv",sep="\t")
				metabolite_names=self.translate_CIDs_from_Stitch_to_PubChem(corr_df["metabolite"]).set_index("PubChem_CID")
				#rename PubChem CIDs to Stitch CIDs in tanimoto_df
				for i in range(len(tanimoto_df.index)):
					pubchem_cid=tanimoto_df.iloc[i,0]
					if tanimoto_df.iloc[i,0]!=int(tanimoto_df.columns[i+1]):
						print("column metabolite does not match row metabolite!")
					else:
						tanimoto_df.iloc[i,0]=metabolite_names.loc[pubchem_cid,"Stitch_CID"]
						tanimoto_df.rename(columns={tanimoto_df.columns[i+1]:str(metabolite_names.loc[pubchem_cid,"Stitch_CID"])},inplace=True)
				#merge both, such that upper triangle is tanimoto coefficients, and lower profile correlation
				corr_df=corr_df.set_index("metabolite")
				corr_df=corr_df.where(np.tril(np.ones(corr_df.shape),k=-1).astype(np.bool)).replace(np.nan,0)
				tanimoto_df=tanimoto_df.set_index("CID")
				tanimoto_df=tanimoto_df.where(np.triu(np.ones(tanimoto_df.shape),k=1).astype(np.bool)).replace(np.nan,0)
				both_df=corr_df+tanimoto_df/100
				np.fill_diagonal(both_df.values,np.nan)
				#rename CIDs to names
				metabolite_names=metabolite_names.reset_index().set_index("Stitch_CID")
				for i in range(len(both_df)):
					both_df.rename(index={both_df.index[i]:metabolite_names.loc[both_df.index[i],"Name"]},inplace=True)
					both_df.rename(columns={both_df.columns[i]:metabolite_names.loc[int(both_df.columns[i]),"Name"]},inplace=True)
				both_df.to_csv(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_tanimoto_map.tsv",sep="\t",index=True,header=True)
			both_df.index.name=""
			#plot
			f,ax=plt.subplots(figsize=(10,10))
			sns.set(font_scale=0.6)
			#ax=sns.heatmap(both_df,annot=False,vmin=-1,vmax=1,linewidths=0.1,cmap="hot_r",xticklabels=True,yticklabels=True)
			lt=sns.heatmap(both_df.where(np.tril(np.ones(both_df.shape),k=-1).astype(np.bool)),annot=False,vmin=-1,vmax=1,linewidths=0.1,cmap=colormap1,xticklabels=True,yticklabels=True,ax=ax,cbar_ax=f.add_axes([0.91,0.15,0.03,0.4]))
			#f.colorbar(lt,shrink=0.6,location="left")
			ut=sns.heatmap(both_df.where(np.triu(np.ones(both_df.shape),k=1).astype(np.bool)),annot=False,vmin=-1,vmax=1,linewidths=0.1,cmap=colormap2,xticklabels=True,yticklabels=True,ax=ax,cbar_ax=f.add_axes([0.91,0.56,0.03,0.4])) #,cbar_kws={"anchor":(0.6,1.0),"panchor":(1.0,0.6),"shrink":0.4}
			cbar_lt=ax.collections[0].colorbar
			cbar_lt.ax.tick_params(labelsize=12)
			cbar_ut=ax.collections[1].colorbar
			cbar_ut.ax.tick_params(labelsize=12)
			f.tight_layout(rect=[0,0,0.9,1])
			plt.show()
			f.savefig(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_tanimoto_map.png")
		except FileNotFoundError:
			try:
				corr_df=pd.read_csv(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_all_metas.tsv",sep="\t").set_index("metabolite")
			except FileNotFoundError:
				complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles_"+method.__name__+".tsv",sep="\t")
				complete_set=complete_set.set_index([self.protcol,"metabolite"])
				metabolites=complete_set.index.get_level_values(1).unique().tolist()
				metabolite_df=complete_set.loc[complete_set.index.get_level_values(0)[0],]
				corr_df=metabolite_df.transpose().corr()
				corr_df.to_csv(self.analyses+"metabolite_"+self.normalized+self.normalized2+"profile_correlation_all_metas.tsv",sep="\t",header=True,index=True)
			#save pubchem cids to list for determination of tanimoto indices
			metabolite_names=self.translate_CIDs_from_Stitch_to_PubChem(corr_df.index)
			metabolite_names["PubChem_CID"].to_csv(self.analyses+"all_metabolite_"+self.normalized+self.normalized2+"profiles_pubchemCIDs.csv",sep="\t",header=False,index=False)
			print("check manually if CIDs correspond to same names as on PubChem") #e.g. Guanosine in Stitch library is 6802, on PubChem 135398635
			print("if discrepancy between Stitch CID and PubChem CID found, please adapt in IDTranslations.translate_CIDs_from_Stitch_to_PubChem()")
			print("upload saved list "+self.analyses+"all_metabolite_pubchemCIDs.csv as ID list to: https://pubchem.ncbi.nlm.nih.gov/score_matrix/score_matrix.cgi")
			print("select second ID list as None")
			print("save output to "+self.analyses+"metabolite_tanimoto_coefficients_all_metas.csv")
			print("then rerun this function TrainingSet.metabolite_profile_correlation_all")
	
	
	
	#determines the tanimoto index for a given list of metabolites, outside of script
	def metabolites_tanimoto_index(self):
		print("check manually if CIDs correspond to same names as on PubChem") #e.g. Guanosine in Stitch library is 6802, on PubChem 135398635
		print("if discrepancy between Stitch CID and PubChem CID found, please adapt in IDTranslations.translate_CIDs_from_Stitch_to_PubChem()")
		print("upload saved list "+self.databases+"metabolites_high_corr.csv as ID list to: https://pubchem.ncbi.nlm.nih.gov/score_matrix/score_matrix.cgi")
		print("select second ID list as None")
		print("save output to "+self.databases+"metas_high_corr_similarity_matrix.csv")
		print("then rerun TrainingSet.get_coeluting_similar_metabolites")
		print("if "+self.databases+"metabolites_high_corr.csv is missing, run TrainingSet.metabolite_profile_correlation")
		
		
	
	
	
	#extracts metabolites with high chemical similarity and high correlation in elution profiles 
	def get_coeluting_similar_metabolites(self,T_thresh=0.85,corr_thresh=0.9):
		meta_corr,metabolite_names=self.metabolite_profile_correlation(threshold=corr_thresh,plot=False)
		try:
			similarity_matrix=pd.read_csv(self.databases+"metas_high_corr_similarity_matrix.csv")
			similarity_matrix=similarity_matrix.set_index("CID",drop=True)
		except FileNotFoundError:
			self.metabolites_tanimoto_index()
			return
		meta_corr.reset_index(inplace=True)
		complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
		metabolites=set(complete_set["metabolite"].tolist())
		for i in meta_corr.index:
			if meta_corr.loc[i,"m1"] and meta_corr.loc[i,"m2"] in metabolites:
				if metabolite_names.loc[meta_corr.loc[i,"m1"],"PubChem_CID"]!=meta_corr.loc[i,"m1"]:
					m1=metabolite_names.loc[meta_corr.loc[i,"m1"],"PubChem_CID"]
				else:
					m1=meta_corr.loc[i,"m1"]
				if metabolite_names.loc[meta_corr.loc[i,"m2"],"PubChem_CID"]!=meta_corr.loc[i,"m2"]:
					m2=str(metabolite_names.loc[meta_corr.loc[i,"m2"],"PubChem_CID"])
				else:
					m2=str(meta_corr.loc[i,"m2"])
				try:
					meta_corr.loc[i,"Tanimoto_coefficient"]=similarity_matrix.loc[m1,m2]/100
				except KeyError:
					print("CID "+str(m1)+" or "+m2+" are not found on PubChem, please try to correct manually, else dropped (metabolite not found in metas_high_corr_similarity_matrix.csv)")
			else:
				meta_corr.loc[i,"Tanimoto_coefficient"]=0 #if metabolite not in complete set, set Tanimoto_coeff to zero to delete it from dataframe
		similar_metabolites=meta_corr[meta_corr["Tanimoto_coefficient"]>T_thresh]
		similar_metabolites.to_csv(self.databases+"similar_metabolites.tsv",index=False,header=True,sep="\t")
		
	
	
	#after already merging all the xsets, add interactions between OGs and their similar metabolites
	def expand_xset_with_similar_metabolites(self,xset,method=np.mean):
		try:
			similar_metabolites=pd.read_csv(self.databases+"similar_metabolites.tsv",sep="\t")
		except FileNotFoundError:
			self.get_coeluting_similar_metabolites(T_thresh=0.85,corr_thresh=0.9)
			try: #if similarity matrix was not downloaded yet
				similar_metabolites=pd.read_csv(self.databases+"similar_metabolites.csv",sep="\t")
			except FileNotFoundError:
				return #download similarity matrix with Tanimoto coefficients
		complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+self.normalized2+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
		complete_set=complete_set.set_index([self.protcol,"metabolite"])
		#complete_set=complete_set.dropna()
		m1s=similar_metabolites["m1"].tolist()
		m2s=similar_metabolites["m2"].tolist()
		for i in xset.index:
			if i[1] in m1s:
				intersect_metas=similar_metabolites[similar_metabolites["m1"]==i[1]]["m2"].tolist()
				for m in intersect_metas:
					xset=xset.append(complete_set.loc[(i[0],m)])
			if i[1] in m2s:
				intersect_metas=similar_metabolites[similar_metabolites["m2"]==i[1]]["m1"].tolist()
				for m in intersect_metas:
					xset=xset.append(complete_set.loc[(i[0],m)])
		return(xset.drop_duplicates()) #drop duplicates for the case, that interaction with similar metabolite was already in xset
	
	
	
	#plots protein and metabolite profile overlay for given index
	def plot_protein_metabolite_overlay_OLD(self,method,pos_index=0,neg_index=None,confidence="low",savefigs=False,plot_feature_engineered=False):
		#load sets
		positive_set=self.load_merged_x_set(x="positive",method=method,feature="")
		negative_set=self.load_merged_x_set(x="negative",method=method,feature="")
		
		#preliminary checks
		if len(positive_set)<=pos_index:
			pos_index=-1
			print("given index exceeds length of positive set, setting to last index")
		if len(negative_set)<=neg_index:
			neg_index=-1
			print("given index exceeds length of negative set, setting to last index")
		#get separations between experiments
		experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
		separators=list()
		num_fractions=0
		for exp in experiments[:-1]:
			num_fractions+=len(list(filter(lambda c: exp in c,positive_set.columns)))/2
			separators.append(num_fractions-0.5)
		#split sets into metabolite and protein profile
		metabolite_cols=list(filter(lambda col: "metabolites" in col,positive_set.columns))
		protein_cols=list(filter(lambda col: "proteins" in col,positive_set.columns))
		metabolite_positive_set=positive_set[metabolite_cols]
		metabolite_negative_set=negative_set[metabolite_cols]
		protein_positive_set=positive_set[protein_cols]
		protein_negative_set=negative_set[protein_cols]
		#plotting
		f,axn=plt.subplots(2,1,figsize=(10,10),num=pos_index)
		axn[0].set_title("Positive Set",weight="bold")
		pos_meta, =axn[0].plot(metabolite_positive_set.iloc[pos_index].tolist(),color="g")
		pos_prot, =axn[0].plot(protein_positive_set.iloc[pos_index].tolist(),color="b")
		ymin=min(metabolite_positive_set.iloc[pos_index].tolist()+protein_positive_set.iloc[pos_index].tolist())
		ymax=max(metabolite_positive_set.iloc[pos_index].tolist()+protein_positive_set.iloc[pos_index].tolist())
		axn[0].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		pos_meta.set_label("Metabolite profile: "+str(positive_set.iloc[pos_index].name[1]))
		pos_prot.set_label("Protein profile: "+positive_set.iloc[pos_index].name[0])
		axn[0].legend()
		if plot_feature_engineered:
			positive_set_featurized=self.load_merged_x_set(x="positive",method=method,feature=None)
			axn[1].set_title("Positive Set feature engineered",weight="bold")
			pos_prot_feat, =axn[1].plot(positive_set_featurized.iloc[pos_index].tolist(),color="r")
			ymin,ymax=axn[1].get_ylim()
			axn[1].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
			pos_prot_feat.set_label("profile"+self.feature+": "+positive_set_featurized.iloc[pos_index].name[0]+","+str(positive_set_featurized.iloc[pos_index].name[1]))
		else:
			axn[1].set_title("Negative Set",weight="bold")
			neg_meta, =axn[1].plot(metabolite_negative_set.iloc[neg_index].tolist(),color="g")
			neg_prot, =axn[1].plot(protein_negative_set.iloc[neg_index].tolist(),color="b")
			ymin=min(metabolite_negative_set.iloc[neg_index].tolist()+protein_negative_set.iloc[neg_index].tolist())
			ymax=max(metabolite_negative_set.iloc[neg_index].tolist()+protein_negative_set.iloc[neg_index].tolist())
			axn[1].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
			neg_meta.set_label("Metabolite profile: "+str(negative_set.iloc[neg_index].name[1]))
			neg_prot.set_label("Protein profile: "+negative_set.iloc[neg_index].name[0])
		axn[1].legend()
		plt.show()
		
	
	
	#plots protein and metabolite profile overlay for given index
	def plot_protein_metabolite_overlay(self,method,metabolite_positive_set=None,protein_positive_set=None,metabolite_negative_set=None,protein_negative_set=None,separators=None,savefigs=False,plot_feature_engineered=False):
		if metabolite_positive_set is None or protein_positive_set is None or metabolite_negative_set is None or protein_negative_set is None or separators is None:
			#load sets
			positive_set=self.load_merged_x_set(x="positive",method=method,feature="")
			negative_set=self.load_merged_x_set(x="negative",method=method,feature="")
			
			#get separations between experiments
			experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
			separators=list()
			num_fractions=0
			for exp in experiments[:-1]:
				num_fractions+=len(list(filter(lambda c: exp in c,positive_set.columns)))/2
				separators.append(num_fractions-0.5)
			#split sets into metabolite and protein profile
			metabolite_cols=list(filter(lambda col: "metabolites" in col,positive_set.columns))
			protein_cols=list(filter(lambda col: "proteins" in col,positive_set.columns))
			metabolite_positive_set=positive_set[metabolite_cols]
			metabolite_negative_set=negative_set[metabolite_cols]
			protein_positive_set=positive_set[protein_cols]
			protein_negative_set=negative_set[protein_cols]
		
		#plotting
		self.f,axn=plt.subplots(2,1,figsize=(10,10),num=self.i)
		axn[0].set_title("Positive Set",weight="bold")
		pos_meta, =axn[0].plot(metabolite_positive_set.loc[self.pos_index].tolist(),color="g")
		pos_prot, =axn[0].plot(protein_positive_set.loc[self.pos_index].tolist(),color="b")
		ymin=min(metabolite_positive_set.loc[self.pos_index].tolist()+protein_positive_set.loc[self.pos_index].tolist())
		ymax=max(metabolite_positive_set.loc[self.pos_index].tolist()+protein_positive_set.loc[self.pos_index].tolist())
		axn[0].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		pos_meta.set_label("Metabolite profile: "+str(self.pos_index[1]))
		pos_prot.set_label("Protein profile: "+self.pos_index[0])
		axn[0].legend()
		if plot_feature_engineered:
			positive_set_featurized=self.load_merged_x_set(x="positive",method=method,feature=None)
			axn[1].set_title("Positive Set feature engineered",weight="bold")
			pos_prot_feat, =axn[1].plot(positive_set_featurized.loc[self.pos_index].tolist(),color="r")
			ymin,ymax=axn[1].get_ylim()
			axn[1].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
			pos_prot_feat.set_label("profile"+self.feature+": "+self.pos_index[0]+","+str(self.pos_index[1]))
		else:
			axn[1].set_title("Negative Set",weight="bold")
			neg_meta, =axn[1].plot(metabolite_negative_set.loc[self.neg_index].tolist(),color="g")
			neg_prot, =axn[1].plot(protein_negative_set.loc[self.neg_index].tolist(),color="b")
			ymin=min(metabolite_negative_set.loc[self.neg_index].tolist()+protein_negative_set.loc[self.neg_index].tolist())
			ymax=max(metabolite_negative_set.loc[self.neg_index].tolist()+protein_negative_set.loc[self.neg_index].tolist())
			axn[1].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
			neg_meta.set_label("Metabolite profile: "+str(self.neg_index[1]))
			neg_prot.set_label("Protein profile: "+self.neg_index[0])
		axn[1].legend()
		plt.show()
		return
	
	
	
	#saves index from plot of protein_metabolite_overlay
	def save_pair(self,confidences,method,savefigs,plot_feature_engineered):
		for confidence in confidences:
			saved_pos_indices=open(self.analyses+"overlay_plots/positive_indices_"+method.__name__+"_"+self.normalized+self.normalized2+"profiles_"+confidence+"_confidence.txt","a")
			saved_pos_indices.write(self.pos_index[0]+"\t"+str(self.pos_index[1])+"\n")
			saved_pos_indices.close()
			if plot_feature_engineered==False:
				saved_neg_indices=open(self.analyses+"overlay_plots/negative_indices_"+method.__name__+self.normalized+self.normalized2+"profiles_"+"_"+confidence+"_confidence.txt","a")
				saved_neg_indices.write(self.neg_index[0]+"\t"+str(self.neg_index[1])+"\n")
				saved_neg_indices.close()
		if savefigs:
			self.f.savefig(self.analyses+"overlay_plots/overlay_protein_metabolite_"+self.normalized+self.normalized2+"profile_i_"+str(pos_index)+"_"+str(neg_index)+method.__name__+self.feature+".png")
		#self.interaction_selection.quit()
		plt.close()
		return
	
	
	# observe co-elution profiles of protein metabolite pairs and save selected ones
	# x={"positive","negative"}
	# confidence={"high","medium","low"} #or whatever
	# if xset is given, it will be returned trimmed, else the x_set from self.load_merged_x_set is taken
	def trim_x_set_manually(self,x,method,x_set=None,savefigs=False):
		if x_set is None:
			x_set=self.load_merged_x_set(x=x,method=method,feature="") #load x_set without feature engineering
		if not os.path.isdir(self.analyses+"overlay_plots/"):
			os.makedirs(self.analyses+"overlay_plots/")
		#load sets
		positive_set=self.load_merged_x_set(x="positive",method=method,feature="")
		negative_set=self.load_merged_x_set(x="negative",method=method,feature="")
		
		#get separations between experiments
		experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
		separators=list()
		num_fractions=0
		for exp in experiments[:-1]:
			num_fractions+=len(list(filter(lambda c: exp in c,positive_set.columns)))/2
			separators.append(num_fractions-0.5)
		#split sets into metabolite and protein profile
		metabolite_cols=list(filter(lambda col: "metabolites" in col,positive_set.columns))
		protein_cols=list(filter(lambda col: "proteins" in col,positive_set.columns))
		metabolite_positive_set=positive_set[metabolite_cols]
		metabolite_negative_set=negative_set[metabolite_cols]
		protein_positive_set=positive_set[protein_cols]
		protein_negative_set=negative_set[protein_cols]
		
		#Dialog to set confidence
		self.interaction_selection=Tk()
		#interaction_selection.withdraw()
		Button(self.interaction_selection,text="high",command=lambda: self.save_pair(["high","low"],method=method,savefigs=savefigs,plot_feature_engineered=True)).pack()
		Button(self.interaction_selection,text="low",command=lambda: self.save_pair(["low"],method=method,savefigs=savefigs,plot_feature_engineered=True)).pack()
		Button(self.interaction_selection,text="None",command=lambda: plt.close).pack()
		self.interaction_selection.after(0,lambda: self.plotloop(x_set,positive_set,negative_set,method,metabolite_positive_set,protein_positive_set,metabolite_negative_set,protein_negative_set,separators=separators,savefigs=False,plot_feature_engineered=True))
		self.interaction_selection.mainloop()
		#interaction_selection.quit()
		return
	
	
	def plotloop(self,x_set,positive_set,negative_set,method,metabolite_positive_set,protein_positive_set,metabolite_negative_set,protein_negative_set,separators,savefigs=False,plot_feature_engineered=True):
		for i in range(len(x_set)):
			self.i=i
			self.interaction_selection.update()
			#preliminary checks
			if len(positive_set)<=i:
				self.pos_index=positive_set.iloc[-1].name
				print("given index exceeds length of positive set, setting to last index")
			else:
				self.pos_index=positive_set.iloc[i].name
			if len(negative_set)<=i:
				self.neg_index=negative_set.iloc[-1].name
				print("given index exceeds length of negative set, setting to last index")
			else:
				self.neg_index=negative_set.iloc[i].name
			self.plot_protein_metabolite_overlay(method,metabolite_positive_set,protein_positive_set,metabolite_negative_set,protein_negative_set,separators=separators,savefigs=False,plot_feature_engineered=True)
			#self.interaction_selection.update() ###
		self.interaction_selection.destroy() #.quit()
		return
	
	
	#mean_as_default=True: if for given method no indices manual selected, take those from mean selection
	def trim_x_set_by_manual_selection(self,x,method,x_set=None,confidence="low",mean_as_default=True):
		if x_set is None:
			x_set=self.load_merged_x_set(x=x,method=method,feature=None)
		#load trimmed indices and feature engineered x_set
		try:
			trimmed_indices=pd.read_csv(self.analyses+"overlay_plots/"+x+"_indices_"+method.__name__+"_"+self.normalized+self.normalized2+"profiles_"+confidence+"_confidence.txt",sep="\t",header=None,names=[self.protcol,"metabolite"])
		except FileNotFoundError:
			#if indices not trimmed, take mean indices
			if mean_as_default==True:
				try:
					trimmed_indices=pd.read_csv(self.analyses+"overlay_plots/"+x+"_indices_"+np.mean.__name__+"_"+self.normalized+self.normalized2+"profiles_"+confidence+"_confidence.txt",sep="\t",header=None,names=[self.protcol,"metabolite"])
				except FileNotFoundError:
					self.trim_x_set_manually(x=x,x_set=None,method=method)
					trimmed_indices=pd.read_csv(self.analyses+"overlay_plots/"+x+"_indices_"+method.__name__+"_"+self.normalized+self.normalized2+"profiles_"+confidence+"_confidence.txt",sep="\t",header=None,names=[self.protcol,"metabolite"])
			else:
				self.trim_x_set_manually(x=x,x_set=None,method=method)
				trimmed_indices=pd.read_csv(self.analyses+"overlay_plots/"+x+"_indices_"+method.__name__+"_"+self.normalized+self.normalized2+"profiles_"+confidence+"_confidence.txt",sep="\t",header=None,names=[self.protcol,"metabolite"])
		#trim and return
		trimmed_xset=x_set.loc[trimmed_indices.set_index([self.protcol,"metabolite"]).index].dropna()
		return(trimmed_xset)
	
	
	
	#plots given number of profiles from positive and negative set
	#indices (optional) is list of indices to plot (integers), num_profiles has to be None
	def plot_training_set(self,method,indices=[0],num_profiles=None,save=False):
		positive_set=self.load_merged_x_set(x="positive",method=method)
		negative_set=self.load_merged_x_set(x="negative",method=method)
		#getting indices of profiles to plot via sampling
		if num_profiles is None:
			pos_profile_indices=indices
			neg_profile_indices=indices
		else:
			if num_profiles>=len(positive_set):
				pos_profile_indices=range(len(positive_set))
				print("plotting all profiles from positive_set")
			else:
				pos_profile_indices=sample(range(len(positive_set)),num_profiles)
			if num_profiles>=len(negative_set):
				neg_profile_indices=range(len(negative_set))
				print("plotting all profiles from negative_set")
			else:
				neg_profile_indices=sample(range(len(negative_set)),num_profiles)
		#get separations between experiments
		experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
		separators=list()
		num_fractions=0
		for exp in experiments[:-1]:
			num_fractions+=len(list(filter(lambda c: exp in c,positive_set.columns)))
			separators.append(num_fractions-0.5)
		#plotting
		f,axn=plt.subplots(2,1,figsize=(10,10))
		axn[0].set_title("Positive Set",weight="bold")
		axn[1].set_title("Negative Set",weight="bold")
		for i in pos_profile_indices:
			axn[0].plot(positive_set.iloc[i].tolist(),color="b")
		ymin,ymax=axn[0].get_ylim()
		axn[0].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		for j in neg_profile_indices:
			axn[1].plot(negative_set.iloc[j].tolist(),color="r")
		ymin,ymax=axn[1].get_ylim()
		axn[1].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		plt.show()
		if save:
			f.savefig(self.analyses+"plot_posneg_set_"+method.__name__+self.feature+".png")
	
	
	
	
	#show histogram of maximum of crosscorr-feature-vector for positive and negative set among experiments
	def distribution_of_maxcrosscorr(self,method):
		positive_set=self.load_merged_x_set(x="positive",method=method)
		negative_set=self.load_merged_x_set(x="negative",method=method)
		f,axn=plt.subplots(len(self.expfiles),2,figsize=(10,10))
		experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
		for i,exp in enumerate(experiments):
			exp_cols=list(filter(lambda c: exp in c,positive_set.columns))
			positive_set_exp=positive_set[exp_cols]
			negative_set_exp=negative_set[exp_cols]
			pos_columns=positive_set_exp.apply(lambda x: x.sort_values(ascending=False).index[0],axis=1).tolist()
			neg_columns=negative_set_exp.apply(lambda x: x.sort_values(ascending=False).index[0],axis=1).tolist()
			if self.mode=="full":
				middle=(len(exp_cols)+1)/2-1
			else:
				middle=len(exp_cols)-1
			pos_columns_numbers=list(map(lambda c: int(c.split("_")[-1])-middle,pos_columns))
			neg_columns_numbers=list(map(lambda c: int(c.split("_")[-1])-middle,neg_columns))
			#plot
			axn[i][0].hist(pos_columns_numbers,color="b")
			ymin,ymax=axn[i][0].get_ylim()
			axn[i][0].vlines(np.mean(pos_columns_numbers),ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
			axn[i][1].hist(neg_columns_numbers,color="r")
			ymin,ymax=axn[i][1].get_ylim()
			axn[i][1].vlines(np.mean(neg_columns_numbers),ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		axn[0][0].set_title("Positive Set",weight="bold")
		axn[0][1].set_title("Negative Set",weight="bold")
		plt.show()
		f.savefig(self.analyses+"histogram_optimal_shifts_"+method.__name__+self.feature+".png")
	
	
	#plot example of interacting and not interacting protein-metabolite profile +feature engineered profile
	def fig_pmrelation_example(self,pos_index,neg_index,method=np.min,feature="_same_crosscorr_not_normalized"):
		positive_set=self.load_merged_x_set(x="positive",method=method,feature="")
		negative_set=self.load_merged_x_set(x="negative",method=method,feature="")
		#get separations between experiments
		experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
		separators=list()
		num_fractions=0
		for exp in experiments[:-1]:
			num_fractions+=len(list(filter(lambda c: exp in c,positive_set.columns)))/2
			separators.append(num_fractions-0.5)
		#split sets into metabolite and protein profile
		metabolite_cols=list(filter(lambda col: "metabolites" in col,positive_set.columns))
		protein_cols=list(filter(lambda col: "proteins" in col,positive_set.columns))
		metabolite_positive_set=positive_set[metabolite_cols]
		metabolite_negative_set=negative_set[metabolite_cols]
		protein_positive_set=positive_set[protein_cols]
		protein_negative_set=negative_set[protein_cols]
		#plot
		f,axn=plt.subplots(4,1,sharex=True)
		#POSITIVE
		axn[0].set_title("elution profiles")#,weight="bold")
		pos_meta, =axn[0].plot(metabolite_positive_set.loc[pos_index].tolist(),color="tab:blue")
		pos_prot, =axn[0].plot(protein_positive_set.loc[pos_index].tolist(),color="tab:green")
		ymin=min(metabolite_positive_set.loc[pos_index].tolist()+protein_positive_set.loc[pos_index].tolist())
		ymax=max(metabolite_positive_set.loc[pos_index].tolist()+protein_positive_set.loc[pos_index].tolist())
		axn[0].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		axn[0].set_ylabel("norm. MS-intensity")
		#pos_meta.set_label("Metabolite profile: "+str(positive_set.loc[pos_index].name[1]))
		#pos_prot.set_label("Protein profile: "+positive_set.loc[pos_index].name[0])
		positive_set_featurized=self.load_merged_x_set(x="positive",method=method,feature=feature)
		axn[1].set_title("cross-correlation profile")#,weight="bold")
		pos_prot_feat, =axn[1].plot(positive_set_featurized.loc[pos_index].tolist(),color="tab:red")
		ymin,ymax=axn[1].get_ylim()
		axn[1].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		axn[1].set_ylabel("cross-corr.")
		#pos_prot_feat.set_label("profile"+self.feature+": "+positive_set_featurized.iloc[pos_index].name[0]+","+str(positive_set_featurized.iloc[pos_index].name[1]))
		#NEGATIVE
		axn[2].set_title("elution profiles")#,weight="bold")
		neg_meta, =axn[2].plot(metabolite_negative_set.loc[neg_index].tolist(),color="tab:blue")
		neg_prot, =axn[2].plot(protein_negative_set.loc[neg_index].tolist(),color="tab:green")
		ymin=min(metabolite_negative_set.loc[neg_index].tolist()+protein_negative_set.loc[neg_index].tolist())
		ymax=max(metabolite_negative_set.loc[neg_index].tolist()+protein_negative_set.loc[neg_index].tolist())
		axn[2].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		axn[2].set_ylabel("norm. MS-intensity")
		#pos_meta.set_label("Metabolite profile: "+str(positive_set.loc[pos_index].name[1]))
		#pos_prot.set_label("Protein profile: "+positive_set.loc[pos_index].name[0])
		negative_set_featurized=self.load_merged_x_set(x="negative",method=method,feature=feature)
		axn[3].set_title("cross-correlation profile")#,weight="bold")
		neg_prot_feat, =axn[3].plot(negative_set_featurized.loc[neg_index].tolist(),color="tab:red")
		ymin,ymax=axn[3].get_ylim()
		axn[3].vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		axn[3].set_ylabel("cross-corr.")
		axn[3].set_xlabel("feature (SEC-fraction per dataset)")
		f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/figures/methods/crosscorr_examples.png")
		f.savefig("../overall_analysis/figures/methods/crosscorr_examples.svg")



###############################################################
###############################################################
###############################################################

#class that constructs feature engineered profiles for all pairs in experimental data to make predictions on



class CompleteSet(XSet):
		
		
	#find orthologous groups for all experimental data files (different species)
	def get_orthologs_between_species(self): #nearly same as in PositiveSet, put more general form into XSet.py
		prot_all_orgs=list()
		#for every excel sheet with experimental data (=organism)...
		for i in range(len(self.expfiles)):
			expfilename=self.expfiles[i]
			org=DBFileInteractions(expfilename,simulation=self.simulation,experimental=self.experimental,databases=self.databases)
			expdata=FileHandler(expfilename,experimental=self.experimental,databases=self.databases)
			uniprot_set=set(expdata.get_proteins())
			# ...find orthologs and save them as tsv file
			if self.simulation==False:
				org.extract_and_save_orthologs(uniprot_set,self.databases+org.database.organism+"_complete_string_orthologs.tsv",overwrite=self.overwrite,include_metas=False)
				df_org_x=pd.read_csv(self.databases+org.database.organism+"_complete_string_orthologs.tsv",sep="\t")
			else:
				org.extract_and_save_orthologs(uniprot_set,"../simulation_data/"+org.database.organism+"_complete_string_orthologs.tsv",overwrite=self.overwrite,include_metas=False)
				df_org_x=pd.read_csv("../simulation_data/"+org.database.organism+"_complete_string_orthologs.tsv",sep="\t")
			# ... put them into a set and determine intersect to previous ones
			'''
			#only if all COG per protein match
			if i==0:
				orthologs=set(df_org_x.COG)
			else:
				self.intersect(orthologs,set(df_org_x.COG))
			'''
			#for every COG on its own 
			list_org_x=df_org_x.COG.values.tolist()
			for listelement in list_org_x:
				if not pd.isna(listelement):
					cogs=str(listelement).split(";")
					for element in cogs:
						prot_all_orgs.append(element)
		orthologs=set(prot_all_orgs)
		return(orthologs)
		
		
	#get metabolites for orthologous proteins, maybe more efficient to include get_orthologs_between_species, since the databases have to be openend once, but maybe faster if intersection between species gets smaller
	#creation of the gold standard: dataframe with orthologous proteins and interacting metabolites
	def get_proteins_for_all_orthologs(self,orthologs=None): #also similar to get_metabolites_for_all_orthologs in PositiveSet, check for putting it into XSet.py, just CID part deleted
		try:
			if self.overwrite==True:
				raise OverwriteException("overwriting of existing file requested")
			df_res=pd.read_csv(self.databases+"ortholog_mapping.tsv",sep="\t")
			df_res=df_res.set_index(list(df_res)[0],drop=True) #set OG as index, coloumn "Unnamed: 0"
			print("dataframe with STRING IDs and not filtered metabolite IDs found and loaded")
		except (OverwriteException,FileNotFoundError):
			if orthologs==None:
				orthologs=self.get_orthologs_between_species()
			df_res=pd.DataFrame(index=orthologs)
			#for every excel sheet with experimental data import corresponding string_orthologs dataframe and merge them according to given orthologeous proteins (orthologs) 
			for i in range(len(self.expfiles)):
				expfilename=self.expfiles[i]
				organism=expfilename[:-5].split("/")[-1] #.split("_")[0]
				if self.simulation==False:
					df_org_x=pd.read_csv(self.databases+organism+"_complete_string_orthologs.tsv",sep="\t")
				else:
					df_org_x=pd.read_csv("../simulation_data/"+organism+"_complete_string_orthologs.tsv",sep="\t")
				if "CID" in list(df_org_x):
					df_org_x=df_org_x.drop("CID",axis=1)
				df_org_x=df_org_x.dropna() #drop rows where there are NaN
				df_org_x=df_org_x.set_index("COG",drop=True)
				for cogs in df_org_x.index: # since cogs are concatenated by ";"
					cogs_list=cogs.split(";")
					for cog in cogs_list:
						if cog in df_res.index: 
							try: #if one COG occurs only once
								df_res.loc[cog,"STRING_ID_"+organism]=df_org_x.loc[cogs,"STRING_ID"] 
							except:
								#since there may be one COG represented several times, put them all to a list
								stringidlist=df_org_x.loc[cogs,"STRING_ID"].values.tolist() 
								#append STRING_IDs to dataframe, set to exclude double occurences
								stringids=set()
								for k in range(len(stringidlist)):
									stringids.add(stringidlist[k])
								df_res.loc[cog,"STRING_ID_"+organism]=";".join(list(stringids))
			df_res.to_csv(self.databases+"ortholog_mapping.tsv",header=True,index=True,sep="\t")
			print("dataframe with STRING IDs and not filtered metabolite IDs created")
		return(df_res)
		
		
		
	#gets metabolite intersection between files and maps all orthologous groups present in intersection of all experiments over all metabolites in intersection between files
	def create_complete_set_profiles(self,complete_set=None):
		metabolites=self.meta_intersect_files()
		metas_string=";".join(metabolites)
		if self.proteinwise==False:
			if complete_set is None:
				complete_set=self.load_set(self.databases+"complete_set.tsv")
				complete_set=complete_set.dropna() #drops rows, where not for all experiments proteins are found to corresponding OG
			#load metabolites and convert set to string
			metabolites=self.meta_intersect_files()
			metas_string=";".join(list(metabolites))
			complete_set["OG"]=complete_set.index
		else:
			if complete_set is None:
				complete_set=pd.DataFrame(columns=["protein","CIDs"])
			proteins=self.prot_intersect_files()
			complete_set["protein"]=list(proteins)
			complete_set.set_index("protein",drop=False,inplace=True)
		complete_set["CIDs"]=metas_string#*len(complete_set)
		#return(complete_set.set_index(["OG","CIDs"],drop=True)) #set OG,CID as index
		return(complete_set)
		
		
		
	# translate positive set to UniProtIDs and saves it
	def save_set(self,tabfile=None,df=None,method=None): #also copied from PositiveSet.py and CID-part deleted
		if tabfile is None:
			tabfile=self.databases+"complete_set.tsv"
		if df is None:
			df=self.get_proteins_for_all_orthologs()
		organisms=list(filter(lambda org: org[:9]=="STRING_ID", list(df)))
		df_res=pd.DataFrame()
		for i in range(len(organisms)):
			org=organisms[i][10:]
			string_proteins=list(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),df["STRING_ID_"+org].dropna().values.tolist()))))
			db=self.translate_proteins(string_proteins,in_id_type="STRING_ID",out_id_types=[""],save=False)
			for ind in df[pd.notnull(df["STRING_ID_"+org])].index:
				string_ids=df.loc[ind,organisms[i]].split(";")
				uniprot_ids=""
				for string_id in string_ids:
					if uniprot_ids=="":
						try:
							uniprot_ids=db.loc[string_id,""].iloc[0]
						except AttributeError:
							uniprot_ids=db.loc[string_id,""]
					else:
						try:
							uniprot_ids=uniprot_ids+";"+db.loc[string_id,""].iloc[0]
						except AttributeError:
							uniprot_ids=uniprot_ids+";"+db.loc[string_id,""]
				df_res.loc[ind,org]=uniprot_ids
		df_res.to_csv(tabfile,header=True,index=True,sep="\t")
		return(df_res)
	



###############################################################
###############################################################
###############################################################

#class to perform predictions using machine learning algorithms and evaluate their performance

class MLClassifiers(IDTranslations):
	#normalized = normalized+normalized2
	def __init__(self,experimental="../experimental_data/",databases="../databases/",analyses="../analyses/",approach="",feature="",normalized="",proteinwise=False):
		self.experimental=experimental
		self.databases=databases
		self.analyses=analyses
		self.approach=approach #defines training set: "_balanced_from_dbs" = balanced, without additional sampling, "_balanced_sampled" = balanced with sampling of positive and negative set
		self.feature=feature #"_product" if product of metabolite and protein profiles should be taken
		self.normalized=normalized
		self.proteinwise=proteinwise
		if proteinwise==True:
			self.protcol="protein"
		else:
			self.protcol="OG"
		self.expfiles=self.get_expfiles()
		
		
	
	
	
	#load positive and negative set, extract corresponding profile with provided method
	def load_known_interactions(self,overwrite=False,method=np.mean,profiles="raw",trim=1.0,random_indices=False,filtered_ogs=None):
		#load positive and negative set and extract profiles with method
		try:
			positive_set=pd.read_csv(self.databases+"positive_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
			positive_set=positive_set.set_index([self.protcol,"metabolite"],drop=True)
		except FileNotFoundError:
			print("construct positive_set"+self.feature+self.approach+"_"+method.__name__+".tsv with TrainingSet.py")
		try:
			negative_set=pd.read_csv(self.databases+"negative_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
			negative_set=negative_set.set_index([self.protcol,"metabolite"],drop=True)
		except FileNotFoundError:
			print("construct negative_set"+self.feature+self.approach+"_"+method.__name__+".tsv with Training_Set.py")
		#positive_set["interaction"]=[True]*len(positive_set)
		#negative_set["interaction"]=[False]*len(negative_set)
		if random_indices==True:#balance and trim negative and positive set randomly
			neg_indices=sample(range(len(negative_set)),floor(min(len(negative_set),len(positive_set))*trim))
			pos_indices=sample(range(len(positive_set)),floor(min(len(negative_set),len(positive_set))*trim))
			known_set=positive_set.iloc[pos_indices,:].append(negative_set.iloc[neg_indices,:])
		elif type(filtered_ogs)==list:
			positive_set["interaction"]=True
			negative_set["interaction"]=False
			pos_ogs=list(set(set(filtered_ogs) & set(positive_set.index.get_level_values(0))))
			neg_ogs=list(set(set(filtered_ogs) & set(negative_set.index.get_level_values(0))))
			print(min(len(negative_set.loc[neg_ogs]),len(positive_set.loc[pos_ogs]))) ### number of indices used per set
			pos_indices=sample(positive_set.loc[pos_ogs].index.tolist(),floor(min(len(negative_set.loc[neg_ogs]),len(positive_set.loc[pos_ogs]))*trim))
			neg_indices=sample(negative_set.loc[neg_ogs].index.tolist(),floor(min(len(negative_set.loc[neg_ogs]),len(positive_set.loc[pos_ogs]))*trim))
			known_set=positive_set.loc[pos_indices,:].append(negative_set.loc[neg_indices,:])
		else: 
			#balance and trim negative and positive set by taking first rows (needed for comparison of methods)
			# neg_indices=range(min(len(positive_set),len(negative_set)))
			# pos_indices=neg_indices
			#balance and trim sets by taking random indices
			neg_indices=sample(range(len(negative_set)),floor(min(len(negative_set),len(positive_set))*trim))
			pos_indices=sample(range(len(positive_set)),floor(min(len(negative_set),len(positive_set))*trim))
			known_set=positive_set.iloc[pos_indices,:].append(negative_set.iloc[neg_indices,:])
		#merge them to known_set dataframe
		return(known_set)
		
		
	#load training set
	#load known_set build in function load_known_interactions
	def load_known_set(self,method,randomized=False,complete_randomized=False,filtered_ogs=None):
		try: #try to load saved set because creation of known_set is really time consuming
			if complete_randomized==True:
				known_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
			elif type(filtered_ogs)==list: #if given list of ogs which have to be taken for training data
				raise OverwriteException("overwriting of existing file with filtered OGs requested")
				#known_set=pd.read_csv(self.databases+"complete_set_normprofiles_"+method.__name__+"_corr_filtered_OGs.tsv",sep="\t")
				#known_set=known_set.set_index(["OG","metabolite"])
			else:
				known_set=pd.read_csv(self.analyses+"training_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
				known_set=known_set.set_index([self.protcol,"metabolite"])
				#multiply given fraction pairs
				#calculate fractions for every experiment of protein and corresponding metabolite profile
			#known_set=self.load_known_interactions(method=method)
		except (OverwriteException,FileNotFoundError):
			if complete_randomized==True:
				print("Please construct complete set with CompleteSet.py")
			elif type(filtered_ogs)==list:
				known_set=self.load_known_interactions(overwrite=False,method=method,profiles="raw",trim=1.0,random_indices=False,filtered_ogs=filtered_ogs)
				known_set.to_csv(self.analyses+"known_set_"+method.__name__+"_corr_filtered_OGs.tsv",header=True,index=True,sep="\t")
			else:
				print("please construct training set with TrainingSet.py")
				print(self.analyses+"training_set"+self.feature+self.approach+"_"+method.__name__+".tsv")
				#known_set=self.load_known_interactions(overwrite=False,method=method,profiles="raw",trim=1.0,random_indices=False)
				#known_set.to_csv(self.analyses+"known_"+self.approach+"set_"+method.__name__+".tsv",header=True,index=True,sep="\t")        
				#print("set for "+method.__name__+" saved")
				return
		return(known_set)
	
	
	
	
	#Random Forest classifier
	def rf_clf(self):
		return(RandomForestClassifier(oob_score=True))
		#return(RandomForestClassifier())
	
	#SVM classifier, kernel={"linear","poly","rbf"}
	def svm_clf(self,C=1,kernel="linear",degree=3,param=None):
		return(svm.SVC(C=C,probability=True,kernel=kernel,degree=degree))
		
		
	#split data into training and test set
	def fit_with_training_and_test_data(self,clf,known_set,training_size):
		train_data,test_data,train_res,test_res=train_test_split(known_set[list(known_set)[:-1]],known_set["interaction"],train_size=training_size)
		clf.fit(train_data,train_res)
		test_pred=clf.predict(test_data)
		#confusion_mat=metrics.confusion_matrix(test_res,test_pred)
		return(metrics.precision_recall_fscore_support(test_res,test_pred,average=None)) #average= {"micro", "macro", "weighted", None}
		
		
	#fits classifier with training data and runs it on training data to estimate training error
	def fit_and_test_trainset(self,clf,known_set,training_size):
		#without cross validation
		train_data,test_data,train_res,test_res=train_test_split(known_set[list(known_set)[:-1]],known_set["interaction"],train_size=training_size)
		clf.fit(train_data,train_res)
		train_pred=clf.predict(train_data)
		#return(metrics.precision_recall_fscore_support(train_res,train_pred,average=None)) #average= {"micro", "macro", "weighted", None}
		return(metrics.accuracy_score(train_res,train_pred)) #change if other measure wished
		
		
	#split data into training, cross validation and test set
	def fit_with_cross_validation(self,clf,known_set,num_kfolds=5,train_size=0.6,test_size=None): #num_kfolds=cv
		cv=ShuffleSplit(n_splits=num_kfolds,train_size=train_size,test_size=test_size)
		accuracies=cross_val_score(clf,known_set[list(known_set)[:-1]],known_set["interaction"],cv=cv,scoring="accuracy")
		train_data,test_data,train_res,test_res=train_test_split(known_set[list(known_set)[:-1]],known_set["interaction"],train_size=train_size)
		f1s=cross_val_score(clf,known_set[list(known_set)[:-1]],known_set["interaction"],cv=cv,scoring="f1")
		precisions=cross_val_score(clf,known_set[list(known_set)[:-1]],known_set["interaction"],cv=cv,scoring="precision")
		recalls=cross_val_score(clf,known_set[list(known_set)[:-1]],known_set["interaction"],cv=cv,scoring="recall")
		clf.fit(train_data,train_res)
		test_pred=clf.predict(test_data)
		tn,fp,fn,tp=metrics.confusion_matrix(test_res,test_pred).ravel()
		sensitivity=tp/(tp+fn)
		specificity=tn/(tn+fp)
		accuracy=metrics.accuracy_score(test_res,test_pred)
		#print("accuracy: "+str(accuracies.mean())+"\t"+str(accuracy))
		return([precisions.mean(),recalls.mean(),f1s.mean(),accuracies.mean(),sensitivity,specificity])
	
	
	#creates summary table with some quality measures for used classifiers
	def measure_table(self,methods,classifiers,reps=100):
		df=pd.DataFrame(columns=["methods","classifiers","precision","recall","f_1 score","accuracy","sensitivity","specificity","precision_sd","recall_sd","f1_sd","accuracy_sd","sens_sd","spec_sd"])
		methods_idx=list(itertools.chain.from_iterable(list(map(lambda method: [method.__name__]*len(classifiers),methods))))
		classifiers_idx=classifiers*len(methods)
		df["classifiers"]=classifiers_idx
		df["methods"]=methods_idx
		df=df.set_index(["methods","classifiers"])
		for method in methods:
			known_set=self.load_known_set(method=method)
			for classifier in classifiers:
				summary=list()
				clf=eval("self."+classifier)
				for i in range(reps):
					summary.append(self.fit_with_cross_validation(clf,known_set))
				df.loc[(method.__name__,classifier)]=list(np.mean(summary,axis=0))+list(np.std(summary,axis=0))
		df.to_csv(self.analyses+"analysis_summary"+self.feature+self.approach+".tsv",index=True,header=True,sep="\t")
		return(df)
		
		
		
	#trims whole known_set
	def trim_known_set(self,known_set,trim):
		#split known_set into positive and negative set
		positive_set=known_set.loc[known_set["interaction"]==True]
		negative_set=known_set.loc[known_set["interaction"]==False]
		#sample indices from trimmed set
		neg_indices=sample(range(len(negative_set)),floor(min(len(negative_set),len(positive_set))*trim))
		pos_indices=sample(range(len(positive_set)),floor(min(len(negative_set),len(positive_set))*trim))
		#merge them to known_set dataframe
		return(positive_set.iloc[pos_indices,:].append(negative_set.iloc[neg_indices,:]))
		
	
	
	#methods=list of methods to combine profiles, classifiers=list of classifiers to create learning curves for (initialization has to be implemented)
	#training_sizes=list of training sizes for learning curve (values between (0,1] ), if randomized True, assign random True/False interaction to combination of positive and negative set,
	#if complete randomized True, draw randomly profiles to from intersection of OGs and metabolites and assign random interaction for training
	def create_learning_curves(self,methods,classifiers,training_sizes=np.linspace(0.1,1.0,5),plot=True,filtered_ogs=None,randomized=False,complete_randomized=False,reps=10,appendix="",legend=True,title=None,len_trainset=None,prob_int=None):
		#colour=["r","b","g","k","k","k","k"] 
		colour=["r","b","g","c","m","y","k","r","b","g","c","m","y","k","r","b","g","c","m","y","k","r","b","g","c","m","y","k"]
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		#classifierlabels=["Random Forest","SVM linear kernel","SVM Gaussian kernel"]
		classifierlabels=classifiers
		df=pd.DataFrame(columns=["methods","classifiers","train_score_mean","train_score_min","train_score_max","test_score_mean","test_score_min","test_score_max"])
		methods_idx=list(itertools.chain.from_iterable(list(map(lambda method: [method.__name__]*len(classifiers),methods))))
		classifiers_idx=classifiers*len(methods)
		df["classifiers"]=classifiers_idx
		df["methods"]=methods_idx
		df=df.set_index(["methods","classifiers"])
		cv=ShuffleSplit(n_splits=10,test_size=0.2)
		if complete_randomized==True:
			randomized=True
		f,axn=plt.subplots(len(methods),len(classifiers),sharey=True,sharex=True,figsize=(10,10))
		ax_ctr=0
		for method in methods:
			known_set=self.load_known_set(method=method,filtered_ogs=None,randomized=randomized,complete_randomized=complete_randomized)
			if type(filtered_ogs)==list:
				known_set=known_set.loc[filtered_ogs]
			if len_trainset is None:
				len_trainset=len(known_set)
				if complete_randomized==True:
					print("please provide length of training set as parameter len_trainset")
					return
			if complete_randomized==True:
				known_set_org=known_set.copy()
				test_set=self.load_known_set(method=method)
			if randomized==True:
				repetitions=reps #assign that many times random interactions and draw known_set 
				if prob_int is None:
					print("please give probability for an interaction as parameter prob_int")
					return
			else:
				repetitions=1
			#f=plt.figure()
			for i in range(len(classifiers)):
				clf=eval("self."+classifiers[i])
				train_scores_mean=list()
				test_scores_mean=list()
				train_scores_min=list()
				train_scores_max=list()
				test_scores_min=list()
				test_scores_max=list()
				for rep in range(repetitions): 
					if complete_randomized==True:
						indices=sample(range(len(known_set_org)),len_trainset)
						known_set=known_set_org.iloc[indices,:]
						known_set=known_set.set_index([self.protcol,"metabolite"])
					if randomized==True:
						known_set["interaction"]=choices([True,False],weights=[prob_int,1-prob_int],k=len(known_set)) #randomly assign if there is an interaction or not
					train_sizes,train_scores,test_scores=learning_curve(clf,known_set[list(known_set)[:-1]],known_set["interaction"],cv=cv,train_sizes=training_sizes,scoring="accuracy")
					train_scores_mean.append(np.mean(train_scores,axis=1))
					train_scores_sd=np.std(train_scores,axis=1)
					train_scores_min.append(np.mean(train_scores,axis=1)-train_scores_sd)
					train_scores_max.append(np.mean(train_scores,axis=1)+train_scores_sd)
					test_scores_mean.append(np.mean(test_scores,axis=1))
					test_scores_sd=np.std(test_scores,axis=1)
					test_scores_min.append(np.mean(test_scores,axis=1)-test_scores_sd)
					test_scores_max.append(np.mean(test_scores,axis=1)+test_scores_sd)
				df.loc[(method.__name__,classifiers[i])]=[";".join(str(n) for n in np.mean(train_scores_mean,axis=0)),";".join(str(n) for n in np.min(train_scores_min,axis=0)),";".join(str(n) for n in np.max(train_scores_max,axis=0)),";".join(str(n) for n in np.mean(test_scores_mean,axis=0)),";".join(str(n) for n in np.min(test_scores_min,axis=0)),";".join(str(n) for n in np.max(test_scores_max,axis=0))]
				if plot==True:
					if type(axn)==np.ndarray:
						ax=axn.flat[ax_ctr]
					else:
						ax=axn
					#plt.errorbar(train_sizes,test_scores_mean,yerr=test_scores_sd,ecolor=colour[i],color=colour[i],label=classifiers[i]+"_cv_score")
					#plt.errorbar(train_sizes,train_scores_mean,yerr=train_scores_sd,ecolor=colour[i+len(classifiers)],color=colour[i+len(classifiers)],label=classifiers[i]+"_train_score")
					ax.fill_between(train_sizes,np.min(train_scores_min,axis=0),np.max(train_scores_max,axis=0),alpha=0.1,color="k")
					ax.fill_between(train_sizes,np.min(test_scores_min,axis=0),np.max(test_scores_max,axis=0),alpha=0.1,color=colour[i])
					ax.plot(train_sizes,np.mean(train_scores_mean,axis=0),color="k",label=classifierlabels[i]+" training score")
					ax.plot(train_sizes,np.mean(test_scores_mean,axis=0),color=colour[i],label=classifierlabels[i]+" test score")
					ax.set_title(abc[ax_ctr],weight="bold")
					#if abc[ax_ctr] in "GHI":
					if type(axn)!=np.ndarray or ax_ctr>=len(axn)-len(classifiers):
						ax.set_xlabel("Training data size")
					if abc[ax_ctr] in "ADG":
						ax.set_ylabel("Accuracy")
					ax_ctr=ax_ctr+1
					print(ax_ctr)
			if plot==True:
				if legend==True:
					plt.legend(loc="lower right")
				#plt.xlabel("Training data size")
				#plt.ylabel("Accuracy")
		if plot==True:
			f.tight_layout(rect=[0,0,0.9,1])
			plt.show()
			if complete_randomized==True:
				f.savefig(self.analyses+"learning_curve_randomized"+self.feature+self.approach+appendix+".png")
			elif randomized==True:
				f.savefig(self.analyses+"learning_curve_randomized_known_set"+self.feature+self.approach+appendix+".png")
			else:
				f.savefig(self.analyses+"learning_curve"+self.feature+self.approach+appendix+".png")
			#plt.close()
		if complete_randomized==True:
			df.to_csv(self.analyses+"learning_curves_randomized"+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
		elif randomized==True:
			df.to_csv(self.analyses+"learning_curves_randomized_known_set"+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
		else:
			df.to_csv(self.analyses+"learning_curves"+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
		return
		
		
	
	#creates figure with learning and ROC curves based on same test set for random and trained classifiers
	def create_learning_ROC_curves_comparison(self,method,classifier,training_sizes=np.linspace(0.1,1.0,5),reps=1000,appendix="",test_size=0.2):
		#def create_learning_curves(self,),filtered_ogs=None,randomized=False,complete_randomized=False,reps=10,appendix="",legend=True,title=None,len_trainset=None,prob_int=None): 
		colour=["tab:red","tab:blue","tab:orange","tab:green","tab:purple","tab:brown","k"]
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		lc_df=pd.DataFrame(columns=["training","train_size","train_score_mean","train_score_min","train_score_max","test_score_mean","test_score_min","test_score_max"])
		lc_df["training"]=["True training","randomized training_set","randomized complete_set"] #,"randomized network"]
		lc_df.set_index("training",inplace=True)
		cv=ShuffleSplit(n_splits=10,test_size=test_size)
		clf=eval("self."+classifier)
		#load training sets: true training = test_set, randomized training_set =test_set, randomized_complete=random_set, network randomized=network_set
		test_set=self.load_known_set(method=method)
		random_set=self.load_known_set(method=method,complete_randomized=True).set_index([self.protcol,"metabolite"])
		#orig_network=pd.read_csv(,"\t") ###
		#orig_network=pd.read_csv(self.analyses+"predictions"+self.feature+self.approach+"_"+method+".tsv",sep="\t").set_index(["OG","metabolite"])[classifier]
		cv_iter=list()
		for train_indices,test_indices in cv.split(test_set[list(test_set)[:-1]]):
			cv_iter.append([train_indices,test_indices])
		#########################
		#learning curve
		#########################
		df_trainsizes=list()
		df_trainscore_norm=list()
		df_trainscore_sd_norm=list()
		df_testscore_norm=list()
		df_testscore_sd_norm=list()
		df_trainscore_rand_ts=list()
		df_trainscore_sd_rand_ts=list()
		df_testscore_rand_ts=list()
		df_testscore_sd_rand_ts=list()
		df_trainscore_rand_cs=list()
		df_trainscore_sd_rand_cs=list()
		df_testscore_rand_cs=list()
		df_testscore_sd_rand_cs=list()
		#df_trainscore_rand_net=list()
		#df_trainscore_sd_rand_net=list()
		#df_testscore_rand_net=list()
		#df_testscore_sd_rand_net=list()
		for train_size in training_sizes:
			len_trainset=int(np.floor(train_size*len(test_set)))
			#accuracies for given train_size
			accuracies_train_normal=list()
			accuracies_test_normal=list()
			accuracies_train_random_ts=list()
			accuracies_test_random_ts=list()
			accuracies_train_random_cs=list()
			accuracies_test_random_cs=list()
			#accuracies_train_random_net=list()
			#accuracies_test_random_net=list()
			for cv_i in cv_iter:
				test_data=test_set.iloc[cv_i[1],:-1]
				test_res=test_set.iloc[cv_i[1],-1]
				#normal
				train_data_normal=test_set.iloc[cv_i[0][:len_trainset],:-1]
				train_res_normal=test_set.iloc[cv_i[0][:len_trainset],-1]
				clf_normal=clone(clf)
				clf_normal.fit(train_data_normal,train_res_normal)
				test_pred_normal=clf_normal.predict(test_data)
				train_pred_normal=clf_normal.predict(train_data_normal)
				accuracies_test_normal.append(metrics.accuracy_score(test_res,test_pred_normal))
				accuracies_train_normal.append(metrics.accuracy_score(train_res_normal,train_pred_normal))
				for r in range(reps):
					#randomized training set
					#train_data_random_ts=train_data_normal
					train_res_random_ts=choices([True,False],k=len_trainset)
					clf_random_ts=clone(clf)
					try:
						clf_random_ts.fit(train_data_normal,train_res_random_ts)
						test_pred_random_ts=clf_random_ts.predict(test_data)
						train_pred_random_ts=clf_random_ts.predict(train_data_normal)
						accuracies_test_random_ts.append(metrics.accuracy_score(test_res,test_pred_random_ts))
						accuracies_train_random_ts.append(metrics.accuracy_score(train_res_random_ts,train_pred_random_ts))
					except ValueError:
						accuracies_test_random_ts.append(np.nan)
						accuracies_train_random_ts.append(np.nan)
					#randomized from complete set
					rand_indices=sample(range(len(random_set)),len_trainset)
					train_data_random_cs=random_set.iloc[rand_indices,:]
					train_res_random_cs=choices([True,False],k=len_trainset)
					clf_random_cs=clone(clf)
					try:
						clf_random_cs.fit(train_data_random_cs,train_res_random_cs)
						test_pred_random_cs=clf_random_cs.predict(test_data)
						train_pred_random_cs=clf_random_cs.predict(train_data_random_cs)
						accuracies_test_random_cs.append(metrics.accuracy_score(test_res,test_pred_random_cs))
						accuracies_train_random_cs.append(metrics.accuracy_score(train_res_random_cs,train_pred_random_cs))
					except:
						accuracies_test_random_cs.append(np.nan)
						accuracies_train_random_cs.append(np.nan)
					#randomized from network randomization
					#network_set=self.randomize_network(orig_network) #########
					#...
			df_trainsizes.append(len_trainset)
			df_trainscore_norm.append(np.nanmean(accuracies_train_normal))
			df_trainscore_sd_norm.append(np.nanstd(accuracies_train_normal))
			df_testscore_norm.append(np.nanmean(accuracies_test_normal))
			df_testscore_sd_norm.append(np.nanstd(accuracies_test_normal))
			df_trainscore_rand_ts.append(np.nanmean(accuracies_train_random_ts))
			df_trainscore_sd_rand_ts.append(np.nanstd(accuracies_train_random_ts))
			df_testscore_rand_ts.append(np.nanmean(accuracies_test_random_ts))
			df_testscore_sd_rand_ts.append(np.nanstd(accuracies_test_random_ts))
			df_trainscore_rand_cs.append(np.nanmean(accuracies_train_random_cs))
			df_trainscore_sd_rand_cs.append(np.nanstd(accuracies_train_random_cs))
			df_testscore_rand_cs.append(np.nanmean(accuracies_test_random_cs))
			df_testscore_sd_rand_cs.append(np.nanstd(accuracies_test_random_cs))
			#df_trainscore_rand_net.append(np.mean(accuracies_train_random_net))
			#df_trainscore_sd_rand_net.append(np.std(accuracies_train_random_net))
			#df_testscore_rand_net.append(np.mean(accuracies_test_random_net))
			#df_testscore_sd_rand_net.append(np.std(accuracies_test_random_net))
		lc_df["train_size"]=";".join(list(map(lambda x: str(x),df_trainsizes)))
		lc_df.loc["True training","train_score_mean"]=";".join(map(str,df_trainscore_norm))
		lc_df.loc["True training","train_score_min"]=";".join(list(map(lambda i: str(df_trainscore_norm[i]-df_trainscore_sd_norm[i]),range(len(df_trainscore_norm)))))
		lc_df.loc["True training","train_score_max"]=";".join(list(map(lambda i: str(df_trainscore_norm[i]+df_trainscore_sd_norm[i]),range(len(df_trainscore_norm)))))
		lc_df.loc["True training","test_score_mean"]=";".join(map(str,df_testscore_norm))
		lc_df.loc["True training","test_score_min"]=";".join(list(map(lambda i: str(df_testscore_norm[i]-df_testscore_sd_norm[i]),range(len(df_testscore_norm)))))
		lc_df.loc["True training","test_score_max"]=";".join(list(map(lambda i: str(df_testscore_norm[i]+df_testscore_sd_norm[i]),range(len(df_testscore_norm)))))
		lc_df.loc["randomized training_set","train_score_mean"]=";".join(map(str,df_trainscore_rand_ts))
		lc_df.loc["randomized training_set","train_score_min"]=";".join(list(map(lambda i: str(df_trainscore_rand_ts[i]-df_trainscore_sd_rand_ts[i]),range(len(df_trainscore_rand_ts)))))
		lc_df.loc["randomized training_set","train_score_max"]=";".join(list(map(lambda i: str(df_trainscore_rand_ts[i]+df_trainscore_sd_rand_ts[i]),range(len(df_trainscore_rand_ts)))))
		lc_df.loc["randomized training_set","test_score_mean"]=";".join(map(str,df_testscore_rand_ts))
		lc_df.loc["randomized training_set","test_score_min"]=";".join(list(map(lambda i: str(df_testscore_rand_ts[i]-df_testscore_sd_rand_ts[i]),range(len(df_testscore_rand_ts)))))
		lc_df.loc["randomized training_set","test_score_max"]=";".join(list(map(lambda i: str(df_testscore_rand_ts[i]+df_testscore_sd_rand_ts[i]),range(len(df_testscore_rand_ts)))))
		lc_df.loc["randomized complete_set","train_score_mean"]=";".join(map(str,df_trainscore_rand_cs))
		lc_df.loc["randomized complete_set","train_score_min"]=";".join(list(map(lambda i: str(df_trainscore_rand_cs[i]-df_trainscore_sd_rand_cs[i]),range(len(df_trainscore_rand_cs)))))
		lc_df.loc["randomized complete_set","train_score_max"]=";".join(list(map(lambda i: str(df_trainscore_rand_cs[i]+df_trainscore_sd_rand_cs[i]),range(len(df_trainscore_rand_cs)))))
		lc_df.loc["randomized complete_set","test_score_mean"]=";".join(map(str,df_testscore_rand_cs))
		lc_df.loc["randomized complete_set","test_score_min"]=";".join(list(map(lambda i: str(df_testscore_rand_cs[i]-df_testscore_sd_rand_cs[i]),range(len(df_testscore_rand_cs)))))
		lc_df.loc["randomized complete_set","test_score_max"]=";".join(list(map(lambda i: str(df_testscore_rand_cs[i]+df_testscore_sd_rand_cs[i]),range(len(df_testscore_rand_cs)))))
		lc_df.to_csv(self.analyses+"learning_curves_trained_and_random_csts"+self.feature+self.approach+appendix+".tsv",sep="\t",index=True,header=True)
		######################
		# ROC curve
		######################
		#all fprs and tprs to be calculated
		fprs_normal=list()
		tprs_normal=list()
		fprs_random_ts=list()
		tprs_random_ts=list()
		fprs_random_cs=list()
		tprs_random_cs=list()
		for cv_i in cv_iter:
			test_data=test_set.iloc[cv_i[1],:-1]
			test_res=test_set.iloc[cv_i[1],-1]
			#normal
			train_data_normal=test_set.iloc[cv_i[0],:-1]
			train_res_normal=test_set.iloc[cv_i[0],-1]
			clf_normal=clone(clf)
			clf_normal.fit(train_data_normal,train_res_normal)
			#test_pred_normal=clf_normal.predict(test_data)
			test_scores_normal=clf_normal.decision_function(test_data)
			fpr_normal,tpr_normal,thresholds_normal=self.my_roc_curve(test_res,test_scores_normal)
			fprs_normal+=fpr_normal
			tprs_normal+=tpr_normal
			for rep in range(reps):
				#randomized training set
				#train_data_random_ts=train_data_normal
				train_res_random_ts=choices([True,False],k=int(np.floor(len(test_set)*(1-test_size))))
				clf_random_ts=clone(clf)
				clf_random_ts.fit(train_data_normal,train_res_random_ts)
				#test_pred_random_ts=clf_tandom_ts.predict(test_data)
				test_scores_random_ts=clf_random_ts.decision_function(test_data)
				fpr_random_ts,tpr_random_ts,thresholds_random_ts=self.my_roc_curve(test_res,test_scores_random_ts)
				fprs_random_ts+=fpr_random_ts
				tprs_random_ts+=tpr_random_ts
				#randomized from complete set
				rand_indices=sample(range(len(random_set)),int(np.floor(len(test_set)*(1-test_size))))
				train_data_random_cs=random_set.iloc[rand_indices,:]
				train_res_random_cs=choices([True,False],k=int(np.floor(len(test_set)*(1-test_size))))
				clf_random_cs=clone(clf)
				clf_random_cs.fit(train_data_random_cs,train_res_random_cs)
				#test_pred_random_cs=clf_tandom_cs.predict(test_data)
				test_scores_random_cs=clf_random_cs.decision_function(test_data)
				fpr_random_cs,tpr_random_cs,thresholds_random_cs=self.my_roc_curve(test_res,test_scores_random_cs)
				fprs_random_cs+=fpr_random_cs
				tprs_random_cs+=tpr_random_cs
		#no doubled fpr-values, or else several fpr-tpr mappings in plot
		#if doubled fpr-value, take mean over tprs
		#better would be global fit of roc-curve
		normal_df=pd.DataFrame(columns=["fpr","tpr"])
		normal_df["fpr"]=list(set(fprs_normal))
		normal_df.set_index("fpr",inplace=True)
		random_ts_df=pd.DataFrame(columns=["fpr","tpr"])
		random_ts_df["fpr"]=list(set(fprs_random_ts))
		random_ts_df.set_index("fpr",inplace=True)
		random_cs_df=pd.DataFrame(columns=["fpr","tpr"])
		random_cs_df["fpr"]=list(set(fprs_random_cs))
		random_cs_df.set_index("fpr",inplace=True)
		#normal
		for fpr in normal_df.index:
			if fprs_normal.count(fpr)>1: #if fpr occuring more than once, average tpr
				indices=[k for k, x in enumerate(fprs_normal) if x==fpr]
				tprs=[tprs_normal[i] for i in indices]
				normal_df.loc[fpr,"tpr"]=np.mean(tprs)
			else:
				normal_df.loc[fpr,"tpr"]=tprs_normal[fprs_normal.index(fpr)]
		#random_ts
		for fpr in random_ts_df.index:
			if fprs_random_ts.count(fpr)>1: #if fpr occuring more than once, average tpr
				indices=[k for k, x in enumerate(fprs_random_ts) if x==fpr]
				tprs=[tprs_random_ts[i] for i in indices]
				random_ts_df.loc[fpr,"tpr"]=np.mean(tprs)
			else:
				random_ts_df.loc[fpr,"tpr"]=tprs_random_ts[fprs_random_ts.index(fpr)]
		#random_cs
		for fpr in random_cs_df.index:
			if fprs_random_cs.count(fpr)>1: #if fpr occuring more than once, average tpr
				indices=[k for k, x in enumerate(fprs_random_cs) if x==fpr]
				tprs=[tprs_random_cs[i] for i in indices]
				random_cs_df.loc[fpr,"tpr"]=np.mean(tprs)
			else:
				random_cs_df.loc[fpr,"tpr"]=tprs_random_cs[fprs_random_cs.index(fpr)]
		#calculate AUC
		normal_df=normal_df.reset_index().sort_values(["fpr"])
		auc_normal=metrics.auc(normal_df["fpr"].tolist(),normal_df["tpr"].tolist())
		random_ts_df=random_ts_df.reset_index().sort_values(["fpr"])
		auc_random_ts=metrics.auc(random_ts_df["fpr"].tolist(),random_ts_df["tpr"].tolist())
		random_cs_df=random_cs_df.reset_index().sort_values(["fpr"])
		auc_random_cs=metrics.auc(random_cs_df["fpr"].tolist(),random_cs_df["tpr"].tolist())
		#save everything to dataframe
		roc_df=pd.DataFrame(columns=["training","AUC","true_positive_rate","false_positive_rate"])
		roc_df["training"]=["True training","randomized training_set","randomized complete_set"] #,"randomized network"]
		roc_df.set_index("training",inplace=True)
		roc_df["AUC"]=[auc_normal,auc_random_ts,auc_random_cs]
		roc_df.loc["True training","true_positive_rate"]=";".join(str(n) for n in normal_df["tpr"].tolist())
		roc_df.loc["True training","false_positive_rate"]=";".join(str(n) for n in normal_df["fpr"].tolist())
		roc_df.loc["randomized training_set","true_positive_rate"]=";".join(str(n) for n in random_ts_df["tpr"].tolist())
		roc_df.loc["randomized training_set","false_positive_rate"]=";".join(str(n) for n in random_ts_df["fpr"].tolist())
		roc_df.loc["randomized complete_set","true_positive_rate"]=";".join(str(n) for n in random_cs_df["tpr"].tolist())
		roc_df.loc["randomized complete_set","false_positive_rate"]=";".join(str(n) for n in random_cs_df["fpr"].tolist())
		roc_df.to_csv(self.analyses+"ROC_curve_trained_and_random_csts"+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
		return
	
	
	
	
	
	#get insight on variance and bias
	def var_bias(self,methods,classifiers,len_trainset,trimvalues,measure="accuracy",training_size=0.8,reps=6):
		#initialize dataframes
		trainset_trimmed=list(map(lambda x: floor(x*len_trainset),trimvalues))
		classifiers_train=list(map(lambda clf: clf+"_train",classifiers))
		cv_df=pd.DataFrame(columns=["classifier","training_size","precision","recall","F_1","accuracy"])
		ind=classifiers+classifiers_train
		cv_df["classifier"]=list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(list(map(lambda clf: [clf]*len(trimvalues),ind))))))
		cv_df["training_size"]=trainset_trimmed*len(ind)
		cv_df=cv_df.set_index(["classifier","training_size"],drop=True)
		for method in methods:
			for i in range(len(list(cv_df))):
				cv_df[list(cv_df)[i]]=np.empty((len(cv_df),0)).tolist()
			known_set=self.load_known_set(method=method)
			for trim in trimvalues:
				for i in range(reps): #sample known_set several times
					known_set=self.trim_known_set(known_set_org,trim=trim)
					for clf_str in classifiers:
						clf=eval("self."+clf_str)
						accuracy=self.fit_and_test_trainset(clf,known_set,training_size=training_size)
						cv_df.loc[(clf_str+"_train",floor(trim*len_trainset)),measure].append(accuracy) 
						clf_cv=eval("self."+clf_str)
						quality_cv=self.fit_with_cross_validation(clf_cv,known_set,num_kfolds=5)
						for j in range(len(list(cv_df))):
							cv_df.loc[(clf_str,floor(trim*len_trainset)),list(cv_df)[j]].append(quality_cv[j])
					if trim==1.0: #if full set used, no need to repeat sampling
						break
			cv_df.to_csv(self.analyses+"training_sizes_var_bias_"+method.__name__+self.approach+".tsv",header=True,index=True,sep="\t")
			self.plot(cv_df,method,measure=measure,appendix="_bias_var_")
		
		
		
	#self implemented function to calculate roc curve
	def my_roc_curve(self,test_res,test_scores):
		tpr=list() #true positive rate
		fpr=list() #false positive rate
		#append scores to dataframe
		test_df=pd.DataFrame(data=test_res)
		test_df.columns=["res"]
		test_df["scores"]=test_scores
		#sort by scores
		test_df=test_df.sort_values(["scores"])
		thresholds=test_df["scores"].values.tolist()
		thresholds[0]=thresholds[0]-1 #set lowest threshold even lower to include all
		for thresh in thresholds:
			test_df["thresh"]=test_df["scores"]>thresh
			try:
				tn,fp,fn,tp=metrics.confusion_matrix(test_df["res"],test_df["thresh"]).ravel()
				tpr.append(tp/(tp+fn))
				fpr.append(fp/(fp+tn))
			except ValueError: 
				continue
		return([fpr,tpr,thresholds])
		
		
		
	#calculate ROC and AUC
	def roc_auc(self,methods,classifiers,training_size,len_trainset=None,reps=100,appendix="",plot=True,legend=False,randomized=False,complete_randomized=False,prob_int=None):
		colours=["r","b","g","c","m","y","r","b","g","c","m","y","k"]
		alpha=np.linspace(start=1,stop=0.2,num=len(methods))
		#alpha=[1,0.4,0.2]
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		classifierlabels_dict={"rf_clf()":"RF","svm_clf(C=0.01,kernel='linear')":"lin. SVM,C=0.01","svm_clf(C=0.1,kernel='linear')":"lin. SVM,C=0.1","svm_clf(C=1,kernel='linear')":"lin. SVM,C=1","svm_clf(C=10,kernel='linear')":"lin. SVM,C=10","svm_clf(C=100,kernel='linear')":"lin. SVM,C=100","svm_clf(kernel='rbf')":"RBF SVM", "svm_clf(kernel='poly',degree=8)":"poly. SVM", "svm_clf(kernel='sigmoid')":"sigm. SVM"}
		classifierlabels=list(map(lambda x: classifierlabels_dict[x],classifiers))
		if plot==True:
			f,axn=plt.subplots(len(methods),len(classifiers),sharey=False,sharex=True,figsize=(10,10))
		ax_ctr=0
		roc_df=pd.DataFrame(columns=["methods","classifiers","AUC","true_positive_rate","false_positive_rate"])
		methods_idx=list(itertools.chain.from_iterable(list(map(lambda method: [method.__name__]*len(classifiers),methods))))
		classifiers_idx=classifiers*len(methods)
		roc_df["classifiers"]=classifiers_idx
		roc_df["methods"]=methods_idx
		roc_df=roc_df.set_index(["methods","classifiers"])
		for m in range(len(methods)):
			#roc_df=pd.DataFrame(index=classifiers,columns=["AUC","true_positive_rate","false_positive_rate"])
			known_set=self.load_known_set(methods[m],complete_randomized=complete_randomized)
			if len_trainset is None:
				len_trainset=len(known_set)
				if complete_randomized==True:
					print("please provide length of training set as parameter len_trainset")
					return
			if complete_randomized==True:
				known_set_org=known_set.copy()
				if prob_int is None:
					print("please give probability for an interaction as parameter prob_int")
					return
			#f=plt.figure()
			#known_set["interaction"]=choices([True,False],k=len(known_set)) #ranomized => AUC=0.5
			rep_df=pd.DataFrame(index=classifiers,columns=["fprs","tprs","thresholds"])
			for column in rep_df.columns:
				rep_df[column]=np.empty((len(rep_df),0)).tolist() 
			for av_rep in range(reps): #100 repetitions to average roc_curve
				cv=ShuffleSplit(n_splits=100,test_size=0.2) #train_size=0.6)
				if complete_randomized==True:
					indices=sample(range(len(known_set_org)),len_trainset)
					known_set=known_set_org.iloc[indices,:]
					known_set=known_set.set_index([self.protcol,"metabolite"])
				if randomized==True:
					known_set["interaction"]=choices([True,False],weights=[prob_int,1-prob_int],k=len(known_set)) #randomly assign if there is an interaction or not
				##########################################
				train_data,test_data,train_res,test_res=train_test_split(known_set[list(known_set)[:-1]],known_set["interaction"],train_size=training_size)
				for i in range(len(classifiers)):
					clf=eval("self."+classifiers[i])
					clf.fit(train_data,train_res)
					test_pred=clf.predict(test_data)
					if classifiers[i]=="rf_clf()":
						#test_scores=clf.oob_decision_function_[:,1]
						test_scores=clf.predict_proba(test_data)[:,1]
					else:
						test_scores=clf.decision_function(test_data)
					false_positive_rate,true_positive_rate,thresholds=self.my_roc_curve(test_res,test_scores) ###,pos_label=True)
					'''
					print(classifiers[i]+", "+method.__name__)
					print("accuracy:\t"+str(metrics.accuracy_score(test_res,test_pred)))
					print(metrics.confusion_matrix(test_res,test_pred))
					print("thresh")
					print(thresholds)
					print("fpr")
					print(false_positive_rate)
					print("tpr")
					print(true_positive_rate)
					'''
					rep_df.loc[classifiers[i],"fprs"]=rep_df.loc[classifiers[i],"fprs"]+false_positive_rate
					rep_df.loc[classifiers[i],"tprs"]=rep_df.loc[classifiers[i],"tprs"]+true_positive_rate
					rep_df.loc[classifiers[i],"thresholds"]=rep_df.loc[classifiers[i],"thresholds"]+thresholds
					##false_positive_rate,true_positive_rate,thresholds=metrics.roc_curve(test_res,test_scores,pos_label=True)
					#r_auc=metrics.auc(false_positive_rate,true_positive_rate)
					#plt.plot(false_positive_rate,true_positive_rate,colours[i],label=classifiers[i]+"(area = %0.2f)" % r_auc)
			for i in range(len(classifiers)):
				fprs=list()
				tprs=list()
				#fpr_indices=list(set(rep_df.ix[classifiers[i],"fprs"]).remove(0))
				for j in range(len(rep_df.loc[classifiers[i],"fprs"])):#range(len(fpr_indices)):
					if rep_df.loc[classifiers[i],"fprs"][j] not in fprs:
						if rep_df.loc[classifiers[i],"fprs"].count(rep_df.loc[classifiers[i],"fprs"][j])>1:
							fpr=rep_df.loc[classifiers[i],"fprs"][j]
							indices=[k for k, x in enumerate(rep_df.loc[classifiers[i],"fprs"]) if x==fpr]
							tpr=np.nanmean([rep_df.loc[classifiers[i],"tprs"][index] for index in indices])
							fprs.append(fpr)
							tprs.append(tpr)
						else:
							fprs.append(rep_df.loc[classifiers[i],"fprs"][j])
							tprs.append(rep_df.loc[classifiers[i],"tprs"][j])
				auc_df=pd.DataFrame()
				auc_df["fprs"]=fprs
				auc_df["tprs"]=tprs
				auc_df=auc_df.sort_values(["fprs"])
				auc_df.dropna(inplace=True)
				r_auc=metrics.auc(auc_df["fprs"].tolist(),auc_df["tprs"].tolist())
				if plot==True:
					ax=axn.flat[ax_ctr]
					ax.plot([0,auc_df.iloc[0,0]],[0,auc_df.iloc[0,1]],color=colours[i],alpha=alpha[m])
					ax.plot(auc_df["fprs"].tolist(),auc_df["tprs"].tolist(),colours[i],alpha=alpha[m],label=classifierlabels[i]+" (AUC = %0.2f)" % r_auc)
				roc_df.loc[(methods[m].__name__,classifiers[i]),"AUC"]=r_auc
				roc_df.loc[(methods[m].__name__,classifiers[i]),"true_positive_rate"]=";".join(str(n) for n in auc_df["tprs"].tolist())
				roc_df.loc[(methods[m].__name__,classifiers[i]),"false_positive_rate"]=";".join(str(n) for n in auc_df["fprs"].tolist())
				#print(classifiers[i]+" "+method.__name__+" fprs")
				#print(fprs)
				#print("tprs")
				#print(tprs)
				if plot==True:
					ax.plot([0,1],[0,1],color="k",linestyle="--")
					ax.set_title(abc[ax_ctr],weight="bold")
					if ax_ctr>=len(axn.flat)-len(classifiers):
						ax.set_xlabel("False positive rate")
					if abc[ax_ctr] in "ADGJM":
						ax.set_ylabel("True positive rate")
					if legend==True:
						ax.legend(loc="lower right",fontsize="small")
					print(ax_ctr)
					ax_ctr=ax_ctr+1
				
			#plt.ylabel("True Positive Rate")
			#plt.xlabel("False Positive Rate")
		if plot==True:
			f.tight_layout(rect=[0,0,0.9,1])
			plt.show()
			f.savefig(self.analyses+"ROC_curve"+self.feature+self.approach+appendix+".png")
			plt.close()
		roc_df.to_csv(self.analyses+"ROC_curve"+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
		return
		
		
		
	#check whether repetitions correlate via rv-coefficients
	def corr_rep(self,overwrite=False):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		prots_mets=list()
		orgs=list()
		for i in range(len(self.expfiles)):
			expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
			organism=self.expfiles[i][:-5].split("/")[-1]
			try:
				if overwrite==True:
					raise OverwriteException("overwriting of existing file requested")
				rv_prot=pd.read_csv(self.analyses+organism+"_protein_heatmap_repetitions.tsv",sep="\t")
				rv_met=pd.read_csv(self.analyses+organism+"_metabolite_heatmap_repetitions.tsv",sep="\t")
			except (OverwriteException,FileNotFoundError):
				rv_prot,rv_met=expdata.profiles_correlation_reps(profiles="raw")
				#saving correlations to "tsv"
				rv_prot.to_csv(self.analyses+organism+"_protein_heatmap_repetitions.tsv",index=True,header=True,sep="\t")
				rv_met.to_csv(self.analyses+organism+"_metabolite_heatmap_repetitions.tsv",index=True,header=True,sep="\t")
			#print(organism)
			#print(rv_prot)
			#print(rv_met)
			#read tsv-files, because else TypeError when creating heatmap from dataframe
			rv_prot=pd.read_csv(self.analyses+organism+"_protein_heatmap_repetitions.tsv",sep="\t")
			rv_prot=rv_prot.set_index(list(rv_prot)[0])
			rv_prot.index.name=""
			rv_prot.columns=list(map(lambda x: x.replace("Experiment_","Repetition "),rv_prot.columns))
			rv_prot.index=list(map(lambda x: x.replace("Experiment_","Repetition "),rv_prot.index))
			rv_met=pd.read_csv(self.analyses+organism+"_metabolite_heatmap_repetitions.tsv",sep="\t")
			rv_met=rv_met.set_index(list(rv_met)[0])
			rv_met.index.name=""
			rv_met.columns=list(map(lambda x: x.replace("Experiment_","Repetition "),rv_met.columns))
			rv_met.index=list(map(lambda x: x.replace("Experiment_","Repetition "),rv_met.index))
			rv_prot=rv_prot.where(np.tril(np.ones(rv_prot.shape),k=-1).astype(np.bool)) #lower triangle belongs to protein
			rv_met=rv_met.where(np.triu(np.ones(rv_met.shape),k=1).astype(np.bool)) #upper triangle belongs to metabolite
			rv_org=rv_prot.multiply(rv_met,fill_value=1)#put all into lists
			prots_mets.append(rv_org)
			orgs.append(organism)
		#plotting
		f,axn=plt.subplots(ceil(len(self.expfiles)/2),2,sharey=False,sharex=False,figsize=(7,7))
		cbar_ax=f.add_axes([0.91,0.3,0.03,0.4])
		for i in range(len(self.expfiles)): #for i,ax in enumerate(axn.flat):
			ax=axn.flat[i]
			ax.set_title(abc[i],weight="bold")
			ax.xaxis.set_ticks_position("none")
			ax.yaxis.set_ticks_position("none")
			sns.heatmap(prots_mets[i],ax=ax,annot=True,cbar=i==0,vmin=0,vmax=1,cbar_ax=None if i else cbar_ax,linewidths=0.5,cmap="hot_r")
		#f.tick_params(which="both",bottom=False,top=False,labelbottom=True)
		f.tight_layout(rect=[0,0,0.9,1])
		plt.show()
		f.savefig(self.analyses+"heatmaps_repetitions.png")
		plt.close()
	
	
	#counts for every protein in every experiment how many OGs where assigned to it
	def count_og_assignments(self):
		try:
			num_ogs=pd.read_csv(self.analyses+"number_of_OGs_per_protein.tsv",sep="\t")
		except FileNotFoundError:
			for i in range(len(self.expfiles)):
				#load proteins from excel file
				expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
				prots=expdata.get_proteins()
				#load OG - UniProt library
				og_assignments=pd.read_csv(self.databases+expdata.organism+"_string_cog_from_expdata.tsv",sep="\t").set_index("UniProt_ID")
				#count for every uniprot ID, how often an OG was assigned
				df=pd.DataFrame(index=prots,columns=[expdata.organism])
				for uniprot_id in df.index:
					if uniprot_id in og_assignments.index:
						if type(og_assignments.loc[uniprot_id,"COG"])==str:
							df.loc[uniprot_id,expdata.organism]=1
						else:
							df.loc[uniprot_id,expdata.organism]=len(og_assignments.loc[uniprot_id,"COG"])
					else:
						df.loc[uniprot_id,expdata.organism]=0
				if i==0:
					num_ogs=df[expdata.organism].value_counts()
				else:
					num_ogs=pd.concat([num_ogs,df[expdata.organism].value_counts()],axis=1)
			num_ogs.fillna(0).to_csv(self.analyses+"number_of_OGs_per_protein.tsv",sep="\t",index=True,header=True)
		return(num_ogs)
		
		
		
		
	
	
	#creates figure with following plots:
	#1st column: histogram showing frequency of pearson correlation of protein profiles per OG
	#2nd column: pie chart with number of proteins per OG (how many OGs have 1 protein, how many two and forth)
	#3rd column: pie chart with number of OGs annotated to proteins
	#xsets is list of e.g. complete_set and training_set, or different complete_sets
	def corr_og_hist(self,xsets,method,profiles="raw",appendix=""):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		alphas=np.linspace(start=0.5,stop=1,num=len(xsets))
		f,axn=plt.subplots(len(xsets),3,sharex=False,sharey=False,figsize=(10,7))
		for x,xset in enumerate(xsets):
			try:
				correlationlists=pd.read_csv(self.analyses+"corrs_"+str(x)+"_"+method.__name__+".tsv",sep="\t")
				correlationlist=correlationlists[method.__name__].tolist()
				num_prots=pd.read_csv(self.analyses+"number_of_proteins_per_OG_"+str(x)+".tsv",sep="\t")
			except FileNotFoundError:
				correlationlists=pd.DataFrame()
				correlationlist=list()
				for i in range(len(self.expfiles)):
					expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
					xset_org=xset[[expdata.organism]]
					xset_org.columns=["UniProt_IDs"]
					xset_org=xset_org.dropna()
					correlations,num_prots_org=expdata.profiles_correlation_ogs(xset_org,profiles=profiles,method=method)
					correlationlist=correlationlist+list(itertools.chain.from_iterable(correlations["correlations"].tolist()))
					if i==0:
						num_prots=num_prots_org
					else:
						num_prots=num_prots.join(num_prots_org)
				correlationlists[method.__name__]=correlationlist
				correlationlists.to_csv(self.analyses+"corrs_"+str(x)+"_"+method.__name__+".tsv",sep="\t")
				num_prots=num_prots.fillna(0).sort_index()
				num_prots["sum"]=num_prots.sum(axis=1)
				num_prots.to_csv(self.analyses+"number_of_proteins_per_OG_"+str(x)+"_"+".tsv",index=True,header=True,sep="\t")
			#trim num_prots
			num_prots_trimmed=pd.DataFrame(columns=["num_prots"],index=list(map(lambda idx: str(idx),range(1,5))))
			num_prots_trimmed["num_prots"]=num_prots.loc[1:4,"sum"].tolist()
			num_prots_trimmed.loc["5-10","num_prots"]=num_prots.loc[5:10,"sum"].sum()
			num_prots_trimmed.loc[">10","num_prots"]=num_prots.loc[11:,"sum"].sum()
			num_ogs=self.count_og_assignments() 
			#first column: histogram
			axn.flat[x*3].set_title(abc[x*3],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			axn.flat[x*3].hist(correlationlist,bins=20,histtype="bar", alpha=alphas[x])
			axn.flat[x*3].set_ylabel("Frequency")
			axn.flat[x*3].set_xlabel("Pearson correlation coefficient")
			#y_align_titles=f.subplotpars.top
			#second column: pie chart #prots/OG
			axn.flat[x*3+1].set_title(abc[x*3+1],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			axn.flat[x*3+1].pie(num_prots_trimmed["num_prots"].tolist(),labels=num_prots_trimmed.index.tolist(),autopct=lambda p: "{:.0f}".format(int(round(p/100*num_prots_trimmed["num_prots"].sum()))))
			if x==0:
				#third column: pie chart #OGs/prot
				axn.flat[x*3+2].set_title(abc[x*3+2],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
				axn.flat[x*3+2].pie(num_ogs.sum(axis=1),labels=num_ogs.index,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_ogs.sum(axis=1).sum()))))
			#if x==len(xsets)-1:
			#	axn.flat[x*3].set_xlabel("Pearson correlation coefficient")
		f.tight_layout(rect=[0,0,0.9,1])
		plt.show()
		f.savefig(self.analyses+"corr_ogs_histograms_"+appendix+"og.png")
		plt.close()
	
	
	#merges uniprot IDs and OGs from positive and negative set
	def merge_pos_neg_orthos(self,method=np.mean):
		cs=CompleteSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,normalized=self.normalized,normalized2="",feature=self.feature,proteinwise=self.proteinwise)
		complete_set=cs.load_set(cs.databases+"complete_set.tsv") #.replace(np.nan,"")
		training_set=self.load_known_set(method=method)
		training_og_set=complete_set.loc[training_set.index.get_level_values(0).unique()]
		return(training_og_set)
		
	
	
	#calculates boxplots for given set
	def corr_boxplots(self,xset,profiles="raw",method=np.mean):
		try:
			df=pd.read_csv(self.analyses+"correlations_inside_og.tsv",sep="\t")
		except FileNotFoundError:
			df=pd.DataFrame(index=xset.index)
			#f=plt.figure(figsize=(20,10))
			for i in range(len(self.expfiles)):
				expdata=FileHandler(self.expfiles[i],experimental=self.experimental,databases=self.databases)
				organism=self.expfiles[i][:-5].split("/")[-1]
				xset_org=xset[[organism]]
				xset_org.columns=["UniProt_IDs"]
				correlations,num_prots=expdata.profiles_correlation_ogs(xset_org,profiles=profiles,method=method)
				if i==0:
					df=correlations.copy()
				else:
					for ii in df.index:
						if correlations.loc[ii,"correlations"]!=[]:
							df.loc[ii,"correlations"]=df.loc[ii,"correlations"]+correlations.loc[ii,"correlations"]
							#df.ix[ii,"correlations"]=list(df.ix[ii,"correlations"],correlations.ix[ii,"correlations"])
			df.to_csv(self.analyses+"correlations_inside_og.tsv",sep="\t",header=True,index=True)
		plotdata=list()
		ogs=list()
		for ii in df.index:
			if df.loc[ii,"correlations"]!=[]: #only take ogs where at least 2 proteins are present
				plotdata.append(df.loc[ii,"correlations"])
				ogs.append(ii)
				#df.ix[ii,"correlations"].boxplot(meanline=True,showmeans=True,showcaps=True,showbox=True,showfliers=False,ax=ax)
		f,axes=plt.subplots(nrows=10,ncols=1,figsize=(20,10))
		for i in range(1,11):
			ax=f.add_subplot(10,1,i)
			if i==10:
				bp=ax.boxplot(plotdata[((i-1)*10):])
				ax.set_xticklabels(ogs[((i-1)*10):])
			else:
				bp=ax.boxplot(plotdata[((i-1)*10):(i*10)])
				ax.set_xticklabels(ogs[((i-1)*10):(i*10)])
		f.tight_layout()
		f.savefig(self.analyses+"boxplots_ogs_reps"+self.approach+".pdf")
		plt.close()
		return(df)
		
		
		
	#filters orthologous groups for OGs which have non-consistent correlation between proteins (below given threshold) and returns OGs as list
	def filter_unimodal_ogs(self,corr_df=None,corr_thresh=0.5):
		if corr_df is None:
			training_og_set=self.merge_pos_neg_orthos()
			corr_df=self.corr_boxplots(training_og_set)
		for og in corr_df.index:
			if corr_df.loc[og,"correlations"]!=[] and min(corr_df.loc[og,"correlations"])>corr_thresh:
				corr_df=corr_df.drop(og)
		return(corr_df.index.tolist())
		
		
		
	#uses only OGs with consistent correlations between proteins (unimodal correlations distrubution) to create learning curve
	def train_with_unimodal_ogs(self,methods,classifiers,training_sizes=np.linspace(0.1,1.0,5),filtered_ogs=None,appendix="_filtered_ogs"):
		if filtered_ogs is None:
			filtered_ogs=self.filter_unimodal_ogs()
		self.create_learning_curves(methods,classifiers,training_sizes,filtered_ogs=filtered_ogs,appendix=appendix,legend=False)
	
	
	
	#loads set to make predictions on, either trimmed or not
	#which={"combined","metabolite","protein","both"} #determines on which profiles decide to exclude underrepresented instances
	def load_prediction_set(self,which=None,method=np.mean):
		if which is None:
			complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+"normprofiles"+self.feature+"_"+method.__name__+".tsv",sep="\t")
			complete_set=complete_set.set_index([self.protcol,"metabolite"])
		else:
			try:
				complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+"normprofiles"+self.feature+"_trimmed_to_"+which+"_in_ts_"+method.__name__+".tsv",sep="\t")
				complete_set=complete_set.set_index([self.protcol,"metabolite"])
			except FileNotFoundError:
				complete_set=self.exclude_underrepresented_profiles_from_prediction_data(which=which,method=method)
		return(complete_set.dropna())
		
		
		
	#makes prediction for orthologous groups occuring in all files (complete set)
	def make_predictions(self,methods,classifiers,known_set=None,reps=100,trim=""):
		for method in methods:
			#load training_data
			if known_set is None:
				known_set=self.load_known_set(method)
			#load complete set and drop rows which are already in training data
			complete_set=self.load_prediction_set(method=method)
			#positive_set=pd.read_csv(self.databases+"positive_set_normprofiles_"+method.__name__+".tsv",sep="\t")
			#negative_set=pd.read_csv(self.databases+"negative_"+self.neg_approach+"set_"+self.db_organism[1]+"_normprofiles_"+method.__name__+".tsv",sep="\t")
			#complete_set=complete_set.set_index(["OG","metabolite"])
			#positive_set=positive_set.set_index(["OG","metabolite"])
			#negative_set=negative_set.set_index(["OG","metabolite"])
			#drop indices that are in training set, considering that some OGs where deleted due to missing IDs
			#complete_set=complete_set.drop(list(filter(lambda x: x in complete_set.index,positive_set.index.tolist())))
			#complete_set=complete_set.drop(list(filter(lambda x: x in complete_set.index,negative_set.index.tolist())))
			#complete_set=complete_set.drop(list(filter(lambda x: x in complete_set.index,known_set.index.tolist())))
			predictions=pd.DataFrame(index=complete_set.index,columns=classifiers)
			for classifier in classifiers:
				if classifier=="rf_clf()": #average over several repetitions and do majority vote
					preds=list()
					for rep in range(reps):
						clf=eval("self."+classifier)
						clf.fit(known_set[list(known_set)[:-1]],known_set["interaction"])
						preds.append(clf.predict(complete_set))
					p=list(map(lambda x: bool(x),np.median(preds,axis=0))) ###
					predictions[classifier]=p
				else:
					clf=eval("self."+classifier)
					clf.fit(known_set[list(known_set)[:-1]],known_set["interaction"])
					predictions[classifier]=clf.predict(complete_set)
			predictions.to_csv(self.analyses+"predictions"+trim+"_"+method.__name__+self.feature+self.approach+".tsv",index=True,header=True,sep="\t")
		return
	
	
	#assortativity of a bipartite network relative to node number in other type
	#g=graph
	def my_assortativity(self,g):
		#define type (assign to group of bipartiteness)
		g.vs["type"]=list(map(lambda x: type(x)==int,g.vs["name"])) #1 for metabolites, 0 for proteins ###adapt
		max_possible_degree_prots=np.sum(g.vs["type"])
		max_possible_degree_metas=len(g.vs)-np.sum(g.vs["type"])
		avg_neighbour_degrees=list()
		degrees=list()
		for v in g.vs:
			if v["type"]==1: #metabolites
				avg_neighbour_degrees.append(np.mean(g.vs[g.neighbors(v)].degree())/max_possible_degree_prots)
				degrees.append(v.degree()/max_possible_degree_metas)
			else: #proteins
				avg_neighbour_degrees.append(np.mean(g.vs[g.neighbors(v)].degree())/max_possible_degree_metas)
				degrees.append(v.degree()/max_possible_degree_prots)
		avg_neighbour_degrees=[0 if np.isnan(x) else x for x in avg_neighbour_degrees]
		degrees=[0 if np.isnan(x) else x for x in degrees]
		assortativity=np.corrcoef(avg_neighbour_degrees,degrees)[0][1]
		return(assortativity)
	
	
	#make average random predictions for training from training set with random interaction labels and randomized from complete set
	#also calculate assortativity and for complete randomization accuracy on training data
	def make_random_predictions(self,method,classifier,reps=10000):
		prediction_set=self.load_prediction_set(method=method)
		df_ass=pd.DataFrame(columns=["assortativity","assortativity_sd","my_assortativity","my_assortativity_sd","accuracy","accuracy_sd"],index=["training","complete"])
		predictions_ts=pd.DataFrame(index=prediction_set.index,columns=["num_true_ts",classifier])
		predictions_cs=pd.DataFrame(index=prediction_set.index,columns=["num_true_cs",classifier])
		predictions_ts["num_true_ts"]=0 #counts number of "True" assignemnts
		predictions_cs["num_true_cs"]=0
		assortativities_ts=list()
		assortativities_cs=list()
		my_assortativities_ts=list()
		my_assortativities_cs=list()
		accuracies_cs=list()
		training_set=self.load_known_set(method)
		len_trainset=len(training_set)
		complete_set=prediction_set.copy()
		for rep in range(reps):
			#from training set
			clf_ts=eval("self."+classifier)
			clf_ts.fit(training_set[list(training_set)[:-1]],choices([True,False],k=len_trainset))
			predictions_ts[classifier]=clf_ts.predict(prediction_set)
			predictions_ts["num_true_ts"]+=predictions_ts[classifier].tolist()
			##create graph for assortativity
			edges_ts=[tuple(x) for x in predictions_ts[predictions_ts.loc[:,classifier]==True].index]
			g_ts=Graph.TupleList(edges_ts,directed=False)
			##include vertices (OGs and metabolites) which have no interacting partners
			vertices_g_ts=set(g_ts.vs["name"])
			vertices_df_ts=set(predictions_ts.index.get_level_values(0).tolist()+predictions_ts.index.get_level_values(1).tolist())
			remaining_vertices_ts=vertices_df_ts-vertices_g_ts
			g_ts.add_vertices(list(remaining_vertices_ts))
			assortativities_ts.append(g_ts.assortativity_degree(directed=False))
			my_assortativities_ts.append(mlc_yeast.my_assortativity(g_ts)) 
			#random profiles from complete set
			clf_cs=eval("self."+classifier)
			rand_indices=sample(range(len(complete_set)),len_trainset)
			train_data_random_cs=complete_set.iloc[rand_indices,:]
			clf_cs.fit(train_data_random_cs,choices([True,False],k=len_trainset))
			predictions_cs[classifier]=clf_cs.predict(prediction_set)
			predictions_cs["num_true_cs"]+=predictions_cs[classifier].tolist()
			accuracies_cs.append(metrics.accuracy_score(training_set["interaction"].tolist(),predictions_cs.loc[training_set.index,classifier].tolist()))
			##create graph for assortativity
			edges_cs=[tuple(x) for x in predictions_cs[predictions_cs.loc[:,classifier]==True].index]
			g_cs=Graph.TupleList(edges_cs,directed=False)
			##include vertices (OGs and metabolites) which have no interacting partners
			vertices_g_cs=set(g_cs.vs["name"])
			vertices_df_cs=set(predictions_cs.index.get_level_values(0).tolist()+predictions_cs.index.get_level_values(1).tolist())
			remaining_vertices_cs=vertices_df_cs-vertices_g_cs
			g_cs.add_vertices(list(remaining_vertices_cs))
			assortativities_cs.append(g_cs.assortativity_degree(directed=False))
			my_assortativities_cs.append(self.my_assortativity(g_cs))
		#if value in prediction >=reps/2, assign True
		predictions_ts[classifier]=predictions_ts["num_true_ts"]>=reps/2
		predictions_ts.to_csv(self.analyses+"predictions_random_ts_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",index=True,header=True,sep="\t")
		predictions_cs[classifier]=predictions_cs["num_true_cs"]>=reps/2
		predictions_cs.to_csv(self.analyses+"predictions_random_cs_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",index=True,header=True,sep="\t")
		#assortativities
		df_ass.loc["training","assortativity"]=np.nanmean(assortativities_ts)
		df_ass.loc["training","assortativity_sd"]=np.nanstd(assortativities_ts)
		df_ass.loc["complete","assortativity"]=np.nanmean(assortativities_cs)
		df_ass.loc["complete","assortativity_sd"]=np.nanstd(assortativities_cs)
		#use own assortativity relative to degree of other node type
		df_ass.loc["training","my_assortativity"]=np.nanmean(my_assortativities_ts)
		df_ass.loc["training","my_assortativity_sd"]=np.nanstd(my_assortativities_ts)
		df_ass.loc["complete","my_assortativity"]=np.nanmean(my_assortativities_cs)
		df_ass.loc["complete","my_assortativity_sd"]=np.nanstd(my_assortativities_cs)
		df_ass.loc["complete","accuracy"]=np.nanmean(accuracies_cs)
		df_ass.loc["complete","accuracy_sd"]=np.nanstd(accuracies_cs)
		df_ass.to_csv(self.analyses+"comparison_to_random_ts_cs_assortativity_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t",index=True,header=True)
		return
		
	
	
	
	#returns prediction probabilities for interactions
	def make_prediction_probabilities(self,methods,classifiers,reps=1000):
		for method in methods:
			training_set=self.load_known_set(method)
			prediction_set=self.load_prediction_set(method=method)
			predictions=pd.DataFrame(index=prediction_set.index,columns=classifiers)
			for classifier in classifiers:
				if classifier=="rf_clf()": #average over several repetitions and do majority vote
					preds=list()
					for rep in range(reps):
						clf=eval("self."+classifier)
						clf.fit(training_set[list(training_set)[:-1]],training_set["interaction"])
						preds.append(list(map(lambda x: x[1],clf.predict_proba(prediction_set)))) #probability for interaction
					predictions[classifier]=np.mean(preds,axis=0) ###
				else:
					clf=eval("self."+classifier)
					clf.fit(training_set[list(training_set)[:-1]],training_set["interaction"])
					predictions[classifier]=list(map(lambda x: x[1],clf.predict_proba(prediction_set))) #probability for interaction
					#predictions[classifier]=list(map(lambda x: x[1],clf.decision_function(prediction_set)))
			predictions.to_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",index=True,header=True,sep="\t")
	
	
	
	#plots histogram of prediction probabilities
	def plot_prediction_probabilities(self,methods,classifiers):
		f,axes=plt.subplots(nrows=len(methods),ncols=len(classifiers),sharex=True,sharey=False,figsize=(20,10))
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		i=0
		for method in methods:
			try: 
				probabilities=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
			except FileNotFoundError:
				probabilities=make_prediction_probabilities(methods,classifiers)
				probabilities=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
			for classifier in classifiers:
				axes[i].hist(probabilities[classifier].tolist(),bins=50)
				axes[i].set_title(abc[i])
				if method==methods[-1]:
					axes[i].set_xlabel("prediction probability")
				if classifier==classifiers[0]:
					axes[i].set_ylabel("frequency")
				i+=1
		plt.show()
		f.savefig(self.analyses+"prediction_probability_distribution.png")
	
	
	#rank orthologous groups by highest amount of interacting metabolites and save top "top"
	#classifier1= classifier to sort
	def rank_proteins(self,methods,classifiers,classifier1=None,top=10,predictions=None,appendix=""):
		if classifier1 is None:
			classifier1=classifiers[0]
		for method in methods:
			if predictions is None:
				predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
				predictions=predictions.set_index([self.protcol,"metabolite"])
			prot_rank=pd.DataFrame(index=predictions.index.get_level_values(0).unique(),columns=classifiers)
			for classifier in classifiers:
				for og in prot_rank.index:
					prot_rank.loc[og,classifier]=predictions.loc[og,classifier].sum()
				prot_rank=prot_rank.sort_values([classifier],ascending=False)
				prots=prot_rank.index[:top]
				top_highest_rank_prots=pd.DataFrame(index=prots,columns=["number of interacting metabolites","metabolites"])
				tip_highest_rank_prots=pd.DataFrame(index=prots,columns=["number of interacting metabolites","metabolites"])
				top_highest_rank_prots["metabolites"]=np.empty((len(top_highest_rank_prots),0)).tolist()
				for og in prots:
					for i in range(len(predictions.loc[og,classifier])):
						if predictions.loc[(og,predictions.loc[og,classifier].index[i]),classifier]==True:
							top_highest_rank_prots.loc[og,"metabolites"].append(predictions.loc[og,classifier].index[i])
					#return(tip_highest_rank_prots)
					tip_highest_rank_prots.loc[og,"metabolites"]=self.translate_CIDs_offline(self.padd_cids(top_highest_rank_prots.loc[og,"metabolites"]))["Name"].tolist()
				tip_highest_rank_prots["number of interacting metabolites"]=prot_rank.loc[prots,classifier].tolist()
				print(tip_highest_rank_prots)
				tip_highest_rank_prots.to_csv(self.analyses+str(top)+"_highest_ranked_"+self.protcol+"s_"+method.__name__+classifier+appendix+".tsv",sep="\t")
			prot_rank=prot_rank.sort_values([classifier1])### resort by first
			prot_rank.to_csv(self.analyses+self.protcol+"_ranks_"+method.__name__+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
			
			
	#rank metabolites by highest amount of interacting orthologous groups
	#classifier1 is classifier to sort final df on
	def rank_metabolites(self,methods,classifiers,classifier1=None,predictions=None,appendix=""):
		if classifier1 is None:
			classifier1=classifiers[0]
		for method in methods:
			if predictions is None:
				predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
				predictions=predictions.set_index(["metabolite",self.protcol])
			meta_rank=pd.DataFrame(index=predictions.index.get_level_values(0).unique(),columns=classifiers)
			for classifier in classifiers:
				for m in meta_rank.index:
					meta_rank.loc[m,classifier]=predictions.loc[m,classifier].sum()
				meta_rank=meta_rank.sort_values([classifier],ascending=False)
				metas=self.translate_CIDs_offline(self.padd_cids(meta_rank.index))
				top_highest_rank_metas=pd.DataFrame(index=metas["Name"],columns=["number of interacting "+self.protcol+"s",self.protcol+"s"])
				top_highest_rank_metas[self.protcol+"s"]=np.empty((len(top_highest_rank_metas),0)).tolist()
				top_highest_rank_metas["number of interacting "+self.protcol+"s"]=meta_rank.loc[meta_rank.index,classifier].tolist()
				for m in range(len(metas)):
					for i in range(len(predictions.loc[meta_rank.index[m],classifier])):
						if predictions.loc[(meta_rank.index[m],predictions.loc[meta_rank.index[m],classifier].index[i]),classifier]==True:
							top_highest_rank_metas.loc[metas.iloc[m,0],self.protcol+"s"].append(predictions.loc[meta_rank.index[m],classifier].index[i])
				print(top_highest_rank_metas)
				top_highest_rank_metas.to_csv(self.analyses+"highest_ranked_metas_"+method.__name__+"_"+classifier+appendix+".tsv",sep="\t")
			meta_rank=meta_rank.sort_values([classifier1])# resort by given classifier
			#translate metabolite cids to metabolite names 
			meta_translation=self.translate_CIDs_offline(self.padd_cids(meta_rank.index.tolist()))
			for m in meta_rank.index: #to make omit shifts due to missing translations
				meta_rank.loc[m,"Metabolite Name"]=meta_translation.loc[self.padd_cids([m])[0],"Name"]
			meta_rank.to_csv(self.analyses+"metabolite_ranks_"+method.__name__+self.feature+self.approach+appendix+".tsv",index=True,header=True,sep="\t")
	
	
	#checks for consistency of predictions (for random forest)
	def check_consistency(self,methods,classifiers):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		for method in methods:
			known_set=self.load_known_set(method)
			complete_set=self.load_prediction_set(method=method)
			f,axn=plt.subplots(len(classifiers),2,sharex=False,sharey=False,figsize=(10,10))
			for c in range(len(classifiers)):
				classifier=classifiers[c]
				reps_pred=pd.DataFrame()
				probs=pd.DataFrame()
				probslist=list()
				for i in range(1,11):
					while True:
						clf=eval("self."+classifier)
						clf.fit(known_set[list(known_set)[:-1]],known_set["interaction"])
						reps_pred[i]=clf.predict(complete_set)
						probs[i]=clf.predict_proba(complete_set)[:,1]
						probslist=probslist+probs[i].tolist()
						if reps_pred[i].sum()!=0 and reps_pred[i].sum()!=len(reps_pred[i]):
							break
				ax=axn.flat[c*2] #axis for histogram
				ax2=axn.flat[c*2+1] #axis for heatmap
				ax.hist(probslist,bins=10,histtype="bar")
				ax.set_title(abc[c*2],weight="bold")
				ax.set_xlabel("Predicted probability for interaction")
				ax.set_ylabel("Frequency")
				sns.heatmap(reps_pred.corr().where(np.tril(np.ones(reps_pred.corr().shape),k=-1).astype(np.bool)),ax=ax2,annot=True,vmin=-1,vmax=1,linewidths=0.5,cmap="hot_r")
				ax2.set_title(abc[2*c+1],weight="bold")
				ax2.set_xlabel("Repetition")
				ax2.set_ylabel("Repetition")
			f.tight_layout(rect=[0,0,0.9,1])
			plt.show()
			f.savefig(self.analyses+"consistency_check_"+method.__name__+self.approach+".png")
			
			
	
	#constructs barplot comparing occurences of metabolites in positive, negative and complete set
	def compare_metabolite_occur_training_pred(self,method=np.mean):
		#load complete set and drop rows which are already in training data
		complete_set=self.load_prediction_set(method=method)
		#positive_set=pd.read_csv(self.databases+"positive_set_normprofiles_"+method.__name__+".tsv",sep="\t")
		#negative_set=pd.read_csv(self.databases+"negative_"+self.neg_approach+"set_"+self.db_organism[1]+"_normprofiles_"+method.__name__+".tsv",sep="\t")
		#positive_set=positive_set.set_index(["OG","metabolite"])
		#negative_set=negative_set.set_index(["OG","metabolite"])
		#drop indices that are in training set, considering that some OGs where deleted due to missing IDs
		#complete_set=complete_set.drop(list(filter(lambda x: x in complete_set.index,positive_set.index.tolist())))
		#complete_set=complete_set.drop(list(filter(lambda x: x in complete_set.index,negative_set.index.tolist())))
		#load training_data
		known_set=self.load_known_set(method=method)
		complete_set=complete_set.drop(list(filter(lambda x: x in complete_set.index,known_set.index.tolist())))
		#create dataframe that stores occurences of every metabolite in training and to be predicted data
		df=pd.DataFrame(index=set(known_set.index.get_level_values(1).unique() | complete_set.index.get_level_values(1).unique()), columns=["prediction_data_occurence","training_data_occurence","positive_training_data_occurence","negative_training_data_occurence"])
		#loop through all metabolites in training set and count there occurences
		for metabolite in df.index:
			df.loc[metabolite,"training_data_occurence"]=len(list(filter(lambda m: m==metabolite,known_set.index.get_level_values(1))))
			df.loc[metabolite,"positive_training_data_occurence"]=len(list(filter(lambda m: m==metabolite,known_set[known_set["interaction"]==True].index.get_level_values(1))))
			df.loc[metabolite,"negative_training_data_occurence"]=len(list(filter(lambda m: m==metabolite,known_set[known_set["interaction"]==False].index.get_level_values(1))))
			df.loc[metabolite,"prediction_data_occurence"]=len(list(filter(lambda m: m==metabolite,complete_set.index.get_level_values(1))))
		#calculate percentage of occurence
		df["training_data_%_occurence"]=df["training_data_occurence"]/df["prediction_data_occurence"]*100
		df["positive_training_data_%_occurence"]=df["positive_training_data_occurence"]/df["prediction_data_occurence"]*100
		df["negative_training_data_%_occurence"]=df["negative_training_data_occurence"]/df["prediction_data_occurence"]*100
		df.to_csv(self.analyses+"metabolite_occurence_training_predictions"+self.feature+self.approach+".tsv",header=True,index=True,sep="\t")
		
		#plot occurences of metabolites in training and prediction data as barplot
		f=plt.figure(figsize=(20,10))
		ind=np.arange(len(df.index))
		p1=plt.bar(ind,df["positive_training_data_%_occurence"].tolist())
		p2=plt.bar(ind,df["negative_training_data_%_occurence"].tolist())
		plt.xticks(ind,df.index,rotation=60)
		plt.legend((p1[0],p2[0]),("positive set","negative set"))
		#plt.show()
		f.savefig(self.analyses+"metabolite_occurences_training_predictions"+self.feature+self.approach+".png")
		print(str(round(len(known_set.index.get_level_values(1).unique())/len(complete_set.index.get_level_values(1).unique())*100,2))+" % of metabolites in training data are represented in predictions")
		
		#same for OGs
		df_og=pd.DataFrame(index=set(known_set.index.get_level_values(0).unique() | complete_set.index.get_level_values(0).unique()), columns=["prediction_data_occurence","training_data_occurence","positive_training_data_occurence","negative_training_data_occurence"])
		#loop through all metabolites in training set and count there occurences
		for og in df_og.index:
			df_og.loc[og,"training_data_occurence"]=len(list(filter(lambda o: o==og,known_set.index.get_level_values(0))))
			df_og.loc[og,"positive_training_data_occurence"]=len(list(filter(lambda o: o==og,known_set[known_set["interaction"]==True].index.get_level_values(0))))
			df_og.loc[og,"negative_training_data_occurence"]=len(list(filter(lambda o: o==og,known_set[known_set["interaction"]==False].index.get_level_values(0))))
			df_og.loc[og,"prediction_data_occurence"]=len(list(filter(lambda o: o==og,complete_set.index.get_level_values(0))))
		#calculate percentage of occurence
		df_og["training_data_%_occurence"]=df_og["training_data_occurence"]/df_og["prediction_data_occurence"]*100
		df_og["positive_training_data_%_occurence"]=df_og["positive_training_data_occurence"]/df_og["prediction_data_occurence"]*100
		df_og["negative_training_data_%_occurence"]=df_og["negative_training_data_occurence"]/df_og["prediction_data_occurence"]*100
		print(str(round(len(known_set.index.get_level_values(0).unique())/len(complete_set.index.get_level_values(0).unique())*100,2))+" % of OGs occuring in training data are represented in predictions")
		df.to_csv(self.analyses+"OG_occurence_training_predictions"+self.feature+self.approach+".tsv",header=True,index=True,sep="\t")
	
	
	
	
	#checks correlations between predictions of different classifiers
	def corr_preds(self,methods,classifiers):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		f,axn=plt.subplots(ceil(len(classifiers)/2),2,sharey=True,sharex=True,figsize=(7,7))
		cbar_ax=f.add_axes([0.91,0.3,0.03,0.4])
		preds=list()
		for method in methods:
			predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
			predictions=predictions.set_index([self.protcol,"metabolite"])
			preds.append(predictions.corr(method="pearson"))
		for i in range(len(classifiers)):
			axn.flat[i].set_title(abc[i],weight="bold")
			axn.flat[i].xaxis.set_ticks_position("none")
			axn.flat[i].yaxis.set_ticks_position("none")
			sns.heatmap(preds[i],ax=axn.flat[i],annot=True,cbar=i == 0,vmin=0,vmax=1,cbar_ax=None if i else cbar_ax,linewidths=0.5,cmap="hot_r")
		if len(classifiers)/2<ceil(len(classifiers)/2):
			axn[-1,-1].axis("off")
		#f.tick_params(which="both",bottom=False,top=False,labelbottom=True)
		f.tight_layout(rect=[0,0,0.9,1])
		plt.show()
		f.savefig(self.analyses+"heatmaps_predictions"+self.feature+self.approach+".png")
		plt.close()
	
	
	
	#averages profiles of metabolites and OGs present in training and prediction data and compares them
	def compare_fractional_coverage(self,profiles="raw"):
		cs=CompleteSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,normalized=self.normalized)
		complete_set=cs.load_set(self.databases+"complete_set.tsv")
		complete_set=complete_set.dropna() #drops rows, where not for all experiments proteins are found to corresponding OG
		complete_set["CIDs"]=";".join(cs.meta_intersect_files())
		ps=PositiveSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,normalized=self.normalized)
		positive_set=ps.load_set(self.databases+"positive_set"+self.feature+self.approach+".tsv")
		ns=NegativeSet(overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,normalized=self.normalized)
		negative_set=ns.load_set(self.databases+"negative_set"+self.feature+self.approach+".tsv")
		
		#get average profile for complete (to be predicted) set (including training data!)
		ogs_pred,metabolites_pred=cs.average_meta_and_prot_profiles_per_file(xset=complete_set,profiles=profiles)
		#get average profile for positive set
		ogs_pos,metabolites_pos=ps.average_meta_and_prot_profiles_per_file(xset=positive_set,profiles=profiles)
		#get average profile for negative set
		ogs_neg,metabolites_neg=ns.average_meta_and_prot_profiles_per_file(xset=negative_set,profiles=profiles)
		ogs_pos=list(map(lambda x: x*0.5,ogs_pos))
		ogs_neg=np.sum([ogs_pos,list(map(lambda x: x*0.5,ogs_neg))],axis=0)
		metabolites_pos=list(map(lambda x: x*0.5,metabolites_pos))
		metabolites_neg=np.sum([metabolites_pos,list(map(lambda x: x*0.5,metabolites_neg))],axis=0)
		
		#plot
		f,axn=plt.subplots(nrows=2,ncols=1,figsize=(20,10))
		fractions=range(1,len(ogs_pred)+1)
		#OGs
		axn[0].plot(fractions,ogs_pred,color="k")
		axn[0].plot(fractions,ogs_neg,color="r")
		axn[0].plot(fractions,ogs_pos,color="b")
		axn[0].fill_between(fractions,ogs_neg,ogs_pred,alpha=0.1,color="k")
		axn[0].fill_between(fractions,ogs_pos,ogs_neg,alpha=0.1,color="r")
		axn[0].fill_between(fractions,[0]*len(fractions),ogs_pos,alpha=0.1,color="b")
		axn[0].set_title("A",weight="bold") #OGs
		#Metabolites
		axn[1].plot(fractions,metabolites_pred,color="k")
		axn[1].plot(fractions,metabolites_neg,color="r")
		axn[1].plot(fractions,metabolites_pos,color="b")
		axn[1].fill_between(fractions,metabolites_neg,metabolites_pred,alpha=0.1,color="k")
		axn[1].fill_between(fractions,metabolites_pos,metabolites_neg,alpha=0.1,color="r")
		axn[1].fill_between(fractions,[0]*len(fractions),metabolites_pos,alpha=0.1,color="b")
		axn[1].set_title("B",weight="bold")
		plt.show()
		f.savefig(self.analyses+"training_set_profile_coverage"+self.approach+".png")
		plt.close()
	
	
	
	def pearson_affinity(self,M):
		return(1-np.array([[pearsonr(a,b)[0] for a in M] for b in M]))
		
		
		
	#determine number of clusters using elbow criterion
	#clustertype={"hierarchical distance-based","hierarchical correlation-based","K-Means"}
	#UNDER CONSTRUCTION
	def determine_n_clusters_for_profile_clustering(self,method,n_clusters_max,clustertype):
		#load profiles for prediction data
		complete_set=self.load_prediction_set(method=method)
		#split into protein and metabolite profiles
		complete_set_metabolites=complete_set[list(filter(lambda col: "metabolite" in col,list(complete_set)))]
		complete_set_metabolites.index=complete_set_metabolites.index.get_level_values(1)
		complete_set_metabolites=complete_set_metabolites.reset_index().drop_duplicates(subset="metabolite",keep="first").set_index("metabolite")
		complete_set_metabolites.columns=range(len(list(complete_set_metabolites)))
		complete_set_proteins=complete_set[list(filter(lambda col: "proteins" in col,list(complete_set)))]
		complete_set_proteins.index=complete_set_proteins.index.get_level_values(0)
		complete_set_proteins=complete_set_proteins.reset_index().drop_duplicates(subset=self.protcol,keep="first").set_index("OG")
		complete_set_proteins.columns=range(len(list(complete_set_proteins)))
		#cluster them
		#based on distance to centroids
		distortions_comb=list()
		distortions_prot=list()
		distortions_meta=list()
		#silhouette indices
		si_comb=list()
		si_prot=list()
		si_meta=list()
		if clustertype=="hierarchical correlation-based":
			#corr=complete_set_proteins.corrwith(complete_set_metabolites,axis=0,method="pearson")
			corr=complete_set_proteins.append(complete_set_metabolites).transpose().corr()
		for k in range(2,n_clusters_max+1):
			#for combined profiles
			if clustertype=="K-Means":
				clusterer_comb=KMeans(n_clusters=k)
				clusterer_prot=KMeans(n_clusters=k)
				clusterer_meta=KMeans(n_clusters=k)
			else:
				clusterer_comb=AgglomerativeClustering(n_clusters=k,linkage="average",affinity="euclidean")
				clusterer_prot=AgglomerativeClustering(n_clusters=k,linkage="average",affinity="euclidean")
				clusterer_meta=AgglomerativeClustering(n_clusters=k,linkage="average",affinity="euclidean")
			if clustertype=="hierarchical correlation-based":
				cluster_labels_comb=clusterer_comb.fit_predict(corr)
				#distortions_comb.append(sum(np.min(cdist(corr,clusterer_comb.cluster_centers_,"euclidean"),axis=1))/corr.shape[0])
				si_comb.append(metrics.silhouette_score(corr,cluster_labels_comb))
			else:
				cluster_labels_comb=clusterer_comb.fit_predict(complete_set)
				#distortions_comb.append(sum(np.min(cdist(complete_set,clusterer_comb.cluster_centers_,"euclidean"),axis=1))/complete_set.shape[0])
				try:
					si_comb.append(metrics.silhouette_score(complete_set,cluster_labels_comb))
				except MemoryError:
					si_comb.append(np.nan)
				#for protein profiles
				cluster_labels_prot=clusterer_prot.fit_predict(complete_set_proteins)
				#distortions_prot.append(sum(np.min(cdist(complete_set_proteins,clusterer_prot.cluster_centers_,"euclidean"),axis=1))/complete_set_proteins.shape[0])
				si_prot.append(metrics.silhouette_score(complete_set_proteins,cluster_labels_prot))
				#for metabolite profiles
				cluster_labels_meta=clusterer_meta.fit_predict(complete_set_metabolites)
				#distortions_meta.append(sum(np.min(cdist(complete_set_metabolites,clusterer_meta.cluster_centers_,"euclidean"),axis=1))/complete_set_metabolites.shape[0])
				si_meta.append(metrics.silhouette_score(complete_set_metabolites,cluster_labels_meta))
				#mi_comb.append(metrics.cluster.normalized_mutual_info_score(#for every cluster ,#to every other cluster))
		#scree plot
		f,axn=plt.subplots(nrows=3,ncols=1,figsize=(20,10))
		#axn[0].plot(range(2,n_clusters_max+1),distortions_comb)
		axn[0].set_title("A - combined profiles",weight="bold")
		if clustertype!="hierarchical correlation-based":
			axn[1].plot(range(2,n_clusters_max+1),distortions_prot)
			axn[2].plot(range(2,n_clusters_max+1),distortions_meta)
			axn[1].set_title("B - protein profiles",weight="bold")
			axn[2].set_title("C - metabolite profiles",weight="bold")
			axn[2].set_xlabel("#clusters")
		list(map(lambda ax: ax.set_ylabel("distortion"),axn))
		plt.show()
		f.savefig(self.analyses+"scree_plot_clustering_predictions_"+clustertype+self.feature+self.approach+".png")
		plt.close()
		#silhouette index
		f,axn=plt.subplots(nrows=3,ncols=1,figsize=(20,10))
		axn[0].plot(range(2,n_clusters_max+1),si_comb)
		axn[0].set_title("A - combined profiles",weight="bold")
		if clustertype!="hierarchical correlation-based":
			axn[1].plot(range(2,n_clusters_max+1),si_prot)
			axn[2].plot(range(2,n_clusters_max+1),si_meta)
			axn[1].set_title("B - protein profiles",weight="bold")
			axn[2].set_title("C - metabolite profiles",weight="bold")
			axn[2].set_xlabel("#clusters")
		list(map(lambda ax: ax.set_ylabel("silhouette index"),axn))
		plt.show()
		f.savefig(self.analyses+"scree_plot_clustering_predictions_si_"+clustertype+self.feature+self.approach+".png")
		plt.close()
		
		
		
		
	#determine number of clusters for hierarchical clustering from dendrogram
	def determine_n_clusters_for_profile_hierarchical_clustering(self,method,trim_last=20):
		#load profiles for prediction data
		complete_set=self.load_prediction_set(method=method)
		#split into protein and metabolite profiles
		complete_set_metabolites=complete_set[list(filter(lambda col: "metabolite" in col,list(complete_set)))]
		complete_set_proteins=complete_set[list(filter(lambda col: "proteins" in col,list(complete_set)))]
		#cluster them
		#plot dendrogram
		f,axn=plt.subplots(nrows=3,ncols=1)
		dend_comb=hierarchy.dendrogram(hierarchy.linkage(complete_set,method="ward"),p=trim_last,truncate_mode="lastp",ax=axn[0])
		dend_prot=hierarchy.dendrogram(hierarchy.linkage(complete_set_proteins,method="ward"),p=trim_last,truncate_mode="lastp",ax=axn[1])
		dend_meta=hierarchy.dendrogram(hierarchy.linkage(complete_set_metabolites,method="ward"),p=trim_last,truncate_mode="lastp",ax=axn[2])
		axn[0].set_title("A - combined profiles",weight="bold")
		axn[1].set_title("B - protein profiles",weight="bold")
		axn[2].set_title("C - metabolite profiles",weight="bold")
		#plt.show()
		f.savefig(self.analyses+"dendrogram_hierachical_clustering_cutted"+self.approach+".png")
		
		
		
		
	#takes as input method and list of number of clusters for the kmeans clustering algorithm with list[0]=#clusters for combined profiles, list[1] for protein profiles and list[2] for metabolite profiles
	#returns dataframe with p-values from one-sided exact fisher test, pvalue under-representation<0.05 means that profile type is under-represented 
	#clustertype={"hierarchical distance-based","hierarchical correlation-based","K-Means"}
	def profile_clustering(self,method=np.mean,clustertype="K-Means",n_clusters=[10,10,29]):
		#initialize dataframe with pvalues for every cluster
		df=pd.DataFrame(columns=["cluster"])
		df["cluster"]=range(max( n_clusters))
		df=df.set_index("cluster")
		#load profiles for prediction data
		complete_set=self.load_prediction_set(method=method)
		training_set=self.load_known_set(method=method)
		#interaction=training_set["interaction"]
		training_set=training_set[list(training_set)[:-1]] #delete column with interaction
		#split profiles into OG and metabolite profiles, change indices
		complete_set_metabolites=complete_set[list(filter(lambda col: "metabolite" in col,list(complete_set)))]
		complete_set_metabolites.index=complete_set_metabolites.index.get_level_values(1)
		#complete_set_metabolites=complete_set_metabolites.reset_index().drop_duplicates(subset="metabolite",keep="first").set_index("metabolite")
		complete_set_proteins=complete_set[list(filter(lambda col: "proteins" in col,list(complete_set)))]
		complete_set_proteins.index=complete_set_proteins.index.get_level_values(0)
		#complete_set_proteins=complete_set_proteins.reset_index().drop_duplicates(subset="OG",keep="first").set_index("OG")
		training_set_metabolites=training_set[list(filter(lambda col: "metabolite" in col,list(training_set)))]
		training_set_metabolites.index=training_set_metabolites.index.get_level_values(1)
		#training_set_metabolites=training_set_metabolites.reset_index().drop_duplicates(subset="metabolite",keep="first").set_index("metabolite")
		training_set_proteins=training_set[list(filter(lambda col: "proteins" in col,list(training_set)))]
		training_set_proteins.index=training_set_proteins.index.get_level_values(0)
		#training_set_proteins=training_set_proteins.reset_index().drop_duplicates(subset="OG",keep="first").set_index("OG")
		#run fisher's exact test one-sided (equals hypergeometric test) or two-sided (equals sum of hypergeometric test on over- and under-representation of cluster)
		for i in range(3): #for combined, protein and metabolite profile
			if i==0:
				profilename="combined"
				dataset1=training_set
				dataset2=complete_set
			elif i==1:
				profilename="protein profiles"
				dataset1=training_set_proteins
				dataset2=complete_set_proteins
			elif i==2:
				profilename="metabolite profiles"
				dataset1=training_set_metabolites
				dataset2=complete_set_metabolites
			#cluster and predict labels
			if clustertype=="K-Means":
				kmeans=KMeans(n_clusters=n_clusters[i])
				dataset2["cluster"]=kmeans.fit_predict(dataset2)
				dataset1["cluster"]=kmeans.predict(dataset1)
			elif "hierarchical" in clustertype:
				if clustertype=="hierarchical correlation-based":
					pass
					#implement if neccessary
					#AgglomerativeClustering(distance_matrix or correlation_matrix)
			for clabel in range(n_clusters[i]):
				#contingency table
				c_in_train=sum(dataset1["cluster"]==clabel)
				c_in_complete=sum(dataset2["cluster"]==clabel)
				rest_in_train=sum(dataset1["cluster"]!=clabel)
				rest_in_complete=sum(dataset2["cluster"]!=clabel)
				oddsratio_u,pvalue_u=fisher_exact([[c_in_train,c_in_complete],[rest_in_train,rest_in_complete]],alternative="less") #alternative={"less","greater","two-sided"}
				oddsratio_o,pvalue_o=fisher_exact([[c_in_train,c_in_complete],[rest_in_train,rest_in_complete]],alternative="greater")
				df.loc[clabel,"p-value under-representation "+profilename]=pvalue_u
				df.loc[clabel,"p-value over-representation "+profilename]=pvalue_o
				df.loc[clabel,profilename.split(" ")[0]]=";".join(str(e) for e in set(dataset2.index[dataset2["cluster"]==clabel].tolist()))
		df.to_csv(self.analyses+"test_on_cluster_representatives"+self.approach+".tsv",header=True,index=True,sep="\t")
		#print(df)
	
	
	
	
	
	#train classifiers reps times on random profiles and look whether they make same predictions as classifiers trained on training data
	#if training and test sets are given, train and test on them, else train on standard training_set and compare predictions on to be predicted data
	def comparison_predictions_made_by_random_training(self,methods,classifiers,training_set=None,test_set=None,reps=1000):
		for method in methods:
			#initialize dataframe with empty lists
			df=pd.DataFrame(index=classifiers,columns=["random","training"])
			df["random"]=[[] for i in range(len(df))]
			df["training"]=[[] for i in range(len(df))]
			#load profiles for prediction data
			complete_set=self.load_prediction_set(method=method)
			if training_set is None: # train on standard training set
				training_set=self.load_known_set(method=method)
			if test_set=="under-represented":
				try:
					clusters=pd.read_csv(self.analyses+"test_on_cluster_representatives"+self.approach+".tsv",sep="\t")
				except FileNotFoundError:
					self.profile_clustering(method=method)
					clusters=pd.read_csv(self.analyses+"test_on_cluster_representatives.tsv",sep="\t")
				#select clusters with over-represented profiles and build intersecting indices with training data
				clusters=clusters[clusters["p-value under-representation "+"combined"]<=alpha]
				indices=list(map(lambda y: eval(y),list(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),clusters["combined"].dropna().tolist()))))))
				underrep_test=complete_set.loc[indices]
			frequencies=list()
			for classifier in classifiers:
				predictions=list()
				for i in range(reps):
					#draw random indices from complete set and train classifier with them and random labels (interactions)
					indices=sample(complete_set.index.tolist(),len(training_set))
					interactions=choices([True,False],k=len(indices))
					clf=eval("self."+classifier)
					clf.fit(complete_set.loc[indices],interactions)
					#make predictions
					if test_set is None: #compare predictions on whole complete set
						predictions.append(clf.predict(complete_set))
					elif test_set=="under-represented":
						predictions.append(clf.predict(underrep_test))
					else: #compare predictions on given test set
						predictions.append(clf.predict(test_set[list(test_set)[:-1]]))
				#average predictions
				df.at[classifier,"random"]=np.median(predictions,axis=0) #alternatively compute empirical probability
				df.at[classifier,"random: percentage of predicted interactions"]=np.sum(df.at[classifier,"random"])/len(df.at[classifier,"random"])*100
				#predictions based on training data
				clf=eval("self."+classifier)
				clf.fit(training_set[list(training_set)[:-1]],training_set["interaction"])
				if test_set is None:
					df.at[classifier,"training"]=clf.predict(complete_set)
				elif test_set=="under-represented":
					df.at[classifier,"training"]=clf.predict(underrep_test)
				else:
					df.at[classifier,"training"]=test_set["interaction"].tolist()
				df.at[classifier,"training: percentage of predicted interactions"]=np.sum(df.at[classifier,"training"])/len(df.at[classifier,"training"])*100 #determines percentage of interactions predicted as True
				#compare for given method how many times classifer and random predict same
				ctr=0
				for j in range(len(df.loc[classifier,"random"])):
					if df.at[classifier,"random"][j]==df.at[classifier,"training"][j]:
						ctr+=1
				frequencies.append(ctr/len(df.at[classifier,"random"]))
			df["random_pred=training_pred"]=frequencies
			if test_set is None:
				df.to_csv(self.analyses+"comparison_predictions_made_by_random_training_"+method.__name__+self.feature+self.approach+".tsv",index=True,header=True,sep="\t")
			elif type(test_set)==str:
				df.to_csv(self.analyses+"comparison_predictions_made_by_random_training_"+method.__name__+self.feature+self.approach+"_"+test_set+".tsv",index=True,header=True,sep="\t")
			else:
				df.to_csv(self.analyses+"comparison_predictions_made_by_random_training_"+method.__name__+self.feature+self.approach+"_selected_test_set.tsv",index=True,header=True,sep="\t")
	
	
	
	#count occurences of metabolites in training data and compare to metabolite ranks
	def count_occurences_of_metabolites_in_training(self,method,classifiers):
		positive_set=pd.read_csv(self.databases+"positive_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
		positive_set=positive_set.set_index([self.protcol,"metabolite"],drop=True)
		negative_set=pd.read_csv(self.databases+"negative_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
		negative_set=negative_set.set_index([self.protcol,"metabolite"],drop=True)
		training_set=self.load_known_set(method=method)
		positive_training_set=training_set[training_set["interaction"]==True]
		negative_training_set=training_set[training_set["interaction"]==False]
		#complete_set=pd.read_csv("../databases/complete_set_normprofiles_"+method.__name__+".tsv",sep="\t")
		try:
			rank_df=pd.read_csv(self.analyses+"metabolite_ranks_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
		except FileNotFoundError:
			self.make_predictions([method],classifiers,reps=1000)
			self.rank_metabolites([method],classifiers)
			rank_df=pd.read_csv(self.analyses+"metabolite_ranks_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
		rank_df=rank_df.set_index("metabolite")
		rank_df=pd.DataFrame(index=rank_df.index,data=rank_df["svm_clf(kernel='linear')"])
		rank_df=rank_df.sort_values(by="svm_clf(kernel='linear')")
		posset_counter=Counter(positive_set.index.get_level_values(1).tolist())
		negset_counter=Counter(negative_set.index.get_level_values(1).tolist())
		postrainset_counter=Counter(positive_training_set.index.get_level_values(1).tolist())
		negtrainset_counter=Counter(negative_training_set.index.get_level_values(1).tolist())
		for metabolite in posset_counter:
			rank_df.loc[metabolite,"occurence in positive set"]=posset_counter[metabolite]
		for metabolite in negset_counter:
			rank_df.loc[metabolite,"occurence in negative set"]=negset_counter[metabolite]
		for metabolite in postrainset_counter:
			rank_df.loc[metabolite,"occurence in positive training set"]=postrainset_counter[metabolite]
		for metabolite in negtrainset_counter:
			rank_df.loc[metabolite,"occurence in negative training set"]=negtrainset_counter[metabolite]
		rank_df.to_csv(self.analyses+"metabolite_ranks_"+method.__name__+"_train_preds"+self.feature+self.approach+".tsv",index=True,header=True,sep="\t")
	
	
	#count occurences of OGs in training data and compare to OG ranks
	def count_occurences_of_OGs_in_training(self,method,classifiers):
		positive_set=pd.read_csv(self.databases+"positive_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
		positive_set=positive_set.set_index([self.protcol,"metabolite"],drop=True)
		negative_set=pd.read_csv(self.databases+"negative_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
		negative_set=negative_set.set_index([self.protcol,"metabolite"],drop=True)
		training_set=self.load_known_set(method=method)
		positive_training_set=training_set[training_set["interaction"]==True]
		negative_training_set=training_set[training_set["interaction"]==False]
		#complete_set=pd.read_csv("../databases/complete_set_normprofiles_"+method.__name__+".tsv",sep="\t")
		try:
			rank_df=pd.read_csv(self.analyses+self.protcol+"_ranks_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
		except FileNotFoundError:
			self.make_predictions([method],classifiers,reps=1000)
			self.rank_metabolites([method],classifiers)
			rank_df=pd.read_csv(self.analyses+self.protcol+"_ranks_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
		rank_df=rank_df.set_index(self.protcol)
		rank_df=pd.DataFrame(index=rank_df.index,data=rank_df["svm_clf(kernel='linear')"])
		rank_df=rank_df.sort_values(by="svm_clf(kernel='linear')")
		posset_counter=Counter(positive_set.index.get_level_values(0).tolist())
		negset_counter=Counter(negative_set.index.get_level_values(0).tolist())
		postrainset_counter=Counter(positive_training_set.index.get_level_values(0).tolist())
		negtrainset_counter=Counter(negative_training_set.index.get_level_values(0).tolist())
		for og in posset_counter:
			rank_df.loc[og,"occurence in positive set"]=posset_counter[og]
		for og in negset_counter:
			rank_df.loc[og,"occurence in negative set"]=negset_counter[og]
		for og in postrainset_counter:
			rank_df.loc[og,"occurence in positive training set"]=postrainset_counter[og]
		for og in negtrainset_counter:
			rank_df.loc[og,"occurence in negative training set"]=negtrainset_counter[og]
		rank_df.to_csv(self.analyses+self.protcol+"_ranks_"+method.__name__+"_train_preds"+self.feature+self.approach+".tsv",index=True,header=True,sep="\t")
	
	
	#profilename={"combined","metabolite profiles", "protein profiles"}
	def train_without_overrepresented_profiles_and_eval(self,classifier,profilename="combined",method=np.mean,alpha=0.05):
		#load dataframe with overrepresented profile-clusters
		try:
			clusters=pd.read_csv(self.analyses+"test_on_cluster_representatives"+self.approach+".tsv",sep="\t")
		except FileNotFoundError:
			self.profile_clustering(method=method)
			clusters=pd.read_csv(self.analyses+"test_on_cluster_representatives"+self.approach+".tsv",sep="\t")
		#select clusters with over-represented profiles and build intersecting indices with training data
		clusters=clusters[clusters["p-value over-representation "+profilename]<=alpha]
		while True:
			#load training set
			training_set=self.load_known_set(method=method)
			intersect_df=pd.DataFrame() #training_set
			if profilename=="metabolite profiles":
				indices=clusters["metabolite"].dropna().tolist()
				intersect_indices=set(set(indices) & set(training_set.index.get_level_values(1).tolist()))
				list(map(lambda meta: print(Counter(training_set[training_set.index.get_level_values(1)==meta].index.get_level_values(1).tolist())),intersect_indices))
				selection=input("choose metabolites to exclude from training (if several,separate by comma without spaces): ")
				#generate test set from intersect_indices
				for meta in selection.split(","):
					intersect_df=intersect_df.append(training_set[training_set.index.get_level_values(1)==int(meta)])
			elif profilename=="protein profiles":
				indices=list(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),clusters["protein"].dropna().tolist()))))
				intersect_indices=set(set(indices) & set(training_set.index.get_level_values(0).tolist()))
				list(map(lambda og: print(Counter(training_set[training_set.index.get_level_values(1)==og].index.get_level_values(1).tolist())),intersect_indices))
				selection=input("choose OGs to exclude from training (if several,separate by comma without spaces): ")
				#generate test set from intersect_indices and exclude them from training set
				for og in selected_indices:
					intersect_df=intersect_df.append(training_set[training_set.index.get_level_values(0)==int(og)])
			elif profilename=="combined":
				indices=list(map(lambda y: eval(y),list(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),clusters["combined"].dropna().tolist()))))))
				intersect_indices=set(set(indices) & set(training_set.index.tolist()))
				for tup in intersect_indices:
					intersect_df=intersect_df.append(training_set[training_set.index==tup])
			#exclude test set from training data
			training_set=training_set.drop(intersect_df.index)
			#train and predict
			clf=eval("self."+classifier)
			clf.fit(training_set[list(training_set)[:-1]],training_set["interaction"])
			predictions=clf.predict(intersect_df[list(intersect_df)[:-1]])
			c=0
			for i in range(len(predictions)):
				if predictions[i]==intersect_df["interaction"][i]:
					c+=1
			c=c/len(predictions)
			print("percentage of correctly predicted interactions: "+str(c))
			loop=input("Continue? [m]etabolite profiles, [p]rotein profiles, [c]ombined, [b]reak: ")
			if loop[0]=="m":
				profilename="metabolite profiles"
			elif loop[0]=="p":
				profilename="protein profiles"
			elif loop[0]=="c":
				profilename="combined"
			else:
				break
		return([training_set,intersect_df])
	
	
	
	#compare predictions of Arabidopsis made by training with and without yeast
	def compare_preds_between_sets_of_experimental_data(self):
		preds_all_org=pd.read_csv("../analyses_all_org/predictions_mean"+self.feature+self.approach+".tsv",sep="\t")
		preds_all_org=preds_all_org.set_index([self.protcol,"metabolite"])
		preds_ara=pd.read_csv("../analyses_Ara_only/predictions_mean"+self.feature+self.approach+".tsv",sep="\t")
		preds_ara=preds_ara.set_index([self.protcol,"metabolite"])
		df=pd.DataFrame(index=set(preds_all_org.index & preds_ara.index))
		for i in df.index:
			df.loc[i,"comparison"]=preds_all_org.loc[i,list(preds_all_org)[1]]==preds_ara.loc[i,list(preds_ara)[1]]
		print(str(round(sum(df["comparison"])/len(df)*100,1))+" % of predictions coincide")
	
	
	
	#excludes profiles belonging to clusters underrepresented in training to make predictions
	#which={"combined","proteins","metabolites","both"}
	#if strict in which, only include profiles occuring in training set
	def exclude_underrepresented_profiles_from_prediction_data(self,which="both",method=np.mean):
		if "strict" not in which:
			try:
				cluster_represents=pd.read_csv(self.analyses+"test_on_cluster_representatives"+self.feature+self.approach+".tsv",sep="\t")
			except FileNotFoundError:
				print("please run self.profile_clustering()")
				return
		complete_set=self.load_prediction_set(method=method)
		training_set=self.load_known_set(method=np.mean)
		if "combined" in which:
			include_pairs=list(set(";".join(cluster_represents[cluster_represents["p-value under-representation combined"]>=0.05]["combined"].tolist()).split(";")))
			include_pairs=list(map(lambda x: eval(x),include_pairs))
			complete_set=complete_set.loc[include_pairs,:]
		if "protein" in which or "both" in which:
			#exclude OGs clustered to in training underrepresented profiles
			if "strict" not in which:
				include_ogs=set(";".join(cluster_represents[cluster_represents["p-value under-representation protein profiles"]>=0.05]["protein"].tolist()).split(";"))
				#exclude OGs not occuring in training set
				include_ogs=set(set(include_ogs) & set(training_set.index.get_level_values(0)))
			else:
				include_ogs=set(training_set.index.get_level_values(0))
			complete_set=complete_set.reset_index().set_index(self.protcol)
			for og in complete_set.index.unique():
				if og not in include_ogs:
					complete_set=complete_set.drop(og)
			complete_set=complete_set.reset_index()
		if "metabolites" in which or "both" in which:
			#exclude metabolites clustered to underrepresented profiles
			if "strict" not in which:
				include_metas=set(";".join(cluster_represents[cluster_represents["p-value under-representation metabolite profiles"]>=0.05]["metabolite"].tolist()).split(";"))
				# exclude metabolites not occuring in training set
				include_metas=set(set(map(lambda m: int(m),include_metas)) & set(training_set.index.get_level_values(1)))
			else:
				include_metas=set(training_set.index.get_level_values(1))
			complete_set=complete_set.reset_index().set_index("metabolite")
			for m in complete_set.index.unique():
				if m not in include_metas:
					complete_set=complete_set.drop(m)
		complete_set=complete_set.reset_index().set_index([self.protcol,"metabolite"])
		predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")
		predictions.set_index([self.protcol,"metabolite"]).loc[complete_set.index].to_csv(self.analyses+"predictions_underrep_"+which+"_excluded_"+method.__name__+self.feature+self.approach+".tsv",header=True,index=True,sep="\t")
		complete_set.to_csv(self.databases+"complete_set_"+self.normalized+"normprofiles"+self.feature+"_trimmed_to_"+which+"_in_ts_"+method.__name__+".tsv",sep="\t",header=True,index=True)
		return(complete_set)
	
	
	
	#plots values of feature attributes (for linear SVM only)
	def plot_feature_attributes(self,method,classifier):
		training_set=self.load_known_set(method=method)
		clf=eval("self."+classifier)
		clf.fit(training_set[list(training_set)[:-1]],training_set["interaction"])
		feature_weights=clf.coef_[0]
		#get separations between experiments
		experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
		separators=list()
		num_fractions=0
		for exp in experiments[:-1]:
			num_fractions+=len(list(filter(lambda c: exp in c,training_set.columns)))
			separators.append(num_fractions-0.5)
		#plot
		f,ax=plt.subplots()
		ax.plot(feature_weights,color="tab:red")
		ymin,ymax=ax.get_ylim()
		ax.vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
		ax.set_xlabel("feature coefficient")
		ax.set_ylabel("weight")
		f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row2_feature_attributes"+self.feature+"_"+method.__name__+"_"+classifier+".png")
		f.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row2_feature_attributes"+self.feature+"_"+method.__name__+"_"+classifier+".svg")
		return
	
	
	
	
	#boxplots with distribution of kappas for agreement of predictions from predicted network and random networks (from training and complete set)
	def boxplots_kappas_random_norm(self,method,classifier,trim="",appendix=""):
		predictions_norm=pd.read_csv(self.analyses+"predictions"+trim+"_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")[[self.protcol,"metabolite",classifier]].set_index([self.protcol,"metabolite"])
		predictions_norm.rename(columns={classifier:"prediction"},inplace=True)
		kappas_ts_cs=list()
		kappas_ts_norm=list()
		kappas_cs_norm=list()
		for r in range(reps):
			self.make_random_predictions(method=method,classifier=classifier,reps=1)
			predictions_training=pd.read_csv(self.analyses+"predictions_random_ts_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
			predictions_training.drop("num_true_ts",axis=1,inplace=True)
			predictions_training.columns=["prediction"]
			predictions_complete=pd.read_csv(self.analyses+"predictions_random_cs_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
			predictions_complete.drop("num_true_cs",axis=1,inplace=True)
			predictions_complete.columns=["prediction"]
			kappas_ts_cs.append(metrics.cohen_kappa_score(predictions_training["prediction"],predictions_complete["prediction"]))
			kappas_ts_norm.append(metrics.cohen_kappa_score(predictions_training["prediction"],predictions_norm["prediction"]))
			kappas_cs_norm.append(metrics.cohen_kappa_score(predictions_complete["prediction"],predictions_norm["prediction"]))
		#remove nans
		kappas_ts_norm=list(filter(lambda x: not np.isnan(x),kappas_ts_norm))
		kappas_cs_norm=list(filter(lambda x: not np.isnan(x),kappas_cs_norm))
		kappas_ts_cs=list(filter(lambda x: not np.isnan(x),kappas_ts_cs))
		#plot
		colours=["tab:red","tab:purple","tab:green"]
		f=plt.figure()
		bplot=plt.boxplot([kappas_ts_norm,kappas_cs_norm,kappas_ts_cs],labels=["pn-rts","pn-rp","rts-rp"],notch=True,patch_artist=True)
		#plt.setp(bplot["boxes"],color=colours[0])
		for patch, colour in zip(bplot['boxes'],colours):
				patch.set_facecolor(colour)
		plt.ylabel("Cohen's kappa")
		f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/figures/comparison_random/boxplots_agreement_random_predicted_"+appendix+".png")
		f.savefig("../overall_analysis/figures/comparison_random/boxplots_agreement_random_predicted_"+appendix+".svg")
		return
	
	
	
	
	# for every metabolite, give most likely (ranked) interacting proteins according to prediction probability greater equal p.
	def rank_interacting_proteins(self,method,classifier,p=0.8):
		predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index("metabolite")
		prediction_probabilities=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index("metabolite")
		#predictions=pd.read_csv(mlc.analyses+"predictions_"+method.__name__+mlc.feature+mlc.approach+".tsv",sep="\t").set_index("metabolite")
		#prediction_probabilities=pd.read_csv(mlc.analyses+"prediction_probabilities_"+method.__name__+mlc.feature+mlc.approach+".tsv",sep="\t").set_index("metabolite")
		#extract interactions above p and which are predicted as True
		prediction_probabilities_clf=pd.DataFrame(index=predictions.index)
		prediction_probabilities_clf[self.protcol]=predictions[self.protcol]
		#prediction_probabilities_clf[mlc.protcol]=predictions[mlc.protcol]
		prediction_probabilities_clf["probability"]=prediction_probabilities[classifier]
		prediction_probabilities_clf["probability >= "+str(p)]=prediction_probabilities[classifier]>=p
		prediction_probabilities_clf["prediction"]=predictions[classifier]
		prediction_probabilities_clf["prediction and probability >= "+str(p)]=prediction_probabilities_clf["probability >= "+str(p)]&prediction_probabilities_clf["prediction"]
		prediction_probabilities_clf=prediction_probabilities_clf[prediction_probabilities_clf["prediction and probability >= "+str(p)]==True]
		#ranking
		metabolites=prediction_probabilities_clf.index.unique()
		max_interactors=np.max(list(Counter(prediction_probabilities_clf.index).values()))
		rank_df=pd.DataFrame(index=range(max_interactors),columns=metabolites)
		for metabolite in metabolites:
			if type(prediction_probabilities_clf.loc[metabolite,self.protcol])==str:
				prot_list=[prediction_probabilities_clf.loc[metabolite,self.protcol]]
			else:
				#prot_list=prediction_probabilities_clf.loc[metabolite,self.protcol].tolist()
				prot_list=prediction_probabilities_clf.loc[metabolite].sort_values(by="probability",ascending=False)[self.protcol].tolist()
			rank_df[metabolite]=prot_list+[np.nan]*(max_interactors-len(prot_list))
		#translate metabolite CIDs to metabolite names
		#trans_df=mlc.translate_CIDs_offline(mlc.padd_cids(rank_df.columns.tolist()))
		trans_df=self.translate_CIDs_offline(self.padd_cids(rank_df.columns.tolist()))
		rank_df.columns=trans_df["Name"].tolist()
		rank_df.to_csv(self.analyses+"ranked_protein_targets_greaterequal_"+str(p)+"_prediction_probability_"+method.__name__+self.feature+self.approach+".tsv",sep="\t",index=False,header=True)
		return



###############################################################
###############################################################
###############################################################

#class to construct network from predictions and make analyses on networks

class NetworkAnalysis():
	#normalized=normalized+normalized2
	#trim= trimming of predictions, e.g. trim over-represented profiles from training
	#classifier is name of classifier used in MLanalysis
	def __init__(self,feature,approach,appendix="/",method=np.mean,trim="",classifier="svm_clf(kernel='linear')",normalized="",proteinwise=False):
		self.databases="../databases"+appendix
		self.analyses="../analyses"+appendix
		self.method=method
		self.feature=feature
		self.normalized=normalized
		self.approach=approach
		self.trim=trim
		self.classifier=classifier
		self.proteinwise=proteinwise
		if proteinwise==True:
			self.protcol="protein"
		else:
			self.protcol="OG"
		self.predictions=self.load_predictions(trim=trim) #predictions
		self.g=self.create_graph_from_preds() #graph
		if self.g.is_bipartite():
			self.g.vs["type"]=list(map(lambda x: type(x)==int,self.g.vs["name"])) #1 for metabolites, 0 for proteins
		else:
			print("WARNING: Graph is not bipartite")
		
	#load predictions, if prediction set trimmed etc, change here
	#e.g. trim="underrep_both_strict_excluded"
	def load_predictions(self,trim="",appendix=""):
		return(pd.read_csv(self.analyses+"predictions"+trim+"_"+self.method.__name__+self.feature+self.approach+appendix+".tsv",sep="\t")[[self.protcol,"metabolite",self.classifier]])
		
		
	#load prediction probabilities
	def load_prediction_probabilities(self,trim=""):
		return(pd.read_csv(self.analyses+"prediction_probabilities"+trim+"_"+self.method.__name__+self.feature+self.approach+".tsv",sep="\t")[[self.protcol,"metabolite",self.classifier]])
	
	
	#takes self predictions (og-wise) and transforms them to predictions for proteins (uniprot_ids)
	def ogwise_to_protwise_preds(self,trim=""):
		try:
			protwise_predictions=pd.read_csv(self.analyses+"predictions_protwise"+trim+"_"+self.method.__name__+self.feature+self.approach+"_"+self.classifier+".tsv",sep="\t")[["protein","metabolite",self.classifier]]
		except FileNotFoundError:
			complete_set=pd.read_csv(self.databases+"complete_set.tsv",sep="\t").dropna()
			complete_set.set_index(complete_set.columns[0],inplace=True)
			protwise_predictions=pd.DataFrame(columns=["protein","metabolite",self.classifier]).set_index(["protein","metabolite"])
			og_predictions=self.predictions.set_index([self.protcol,"metabolite"])
			for og in og_predictions.index:
				proteins=set(";".join(complete_set.loc[og[0]].tolist()).split(";"))
				for protein in proteins:
					protwise_predictions.loc[(protein,og[1]),self.classifier]=og_predictions.loc[og,self.classifier]
			#protwise_predictions=protwise_predictions.reset_index().rename({"Unnamed:0":"protein","index":"protein"},axis=1)
			protwise_predictions.to_csv(self.analyses+"predictions_protwise"+trim+"_"+self.method.__name__+self.feature+self.approach+"_"+self.classifier+".tsv",sep="\t",header=True,index=True)
		return(protwise_predictions.reset_index())
	
	
	#takes self predictions (prot-wise) and transforms them to og-wise predictions
	#other is og-wise experiment on same datasets
	#intersect={"all","majority"}, describes how an interaction of an OG is defined from its proteins, if all proteins or the majority have to interact
	def protwise_to_ogwise_preds(self,other,trim="",intersect="all"):
		try:
			ogwise_predictions=pd.read_csv(self.analyses+"predictions_ogwise_"+intersect+trim+"_"+self.method.__name__+self.feature+self.approach+"_"+self.classifier+".tsv",sep="\t")[["OG","metabolite",self.classifier]]
		except FileNotFoundError:
			complete_set=pd.read_csv(other.databases+"complete_set.tsv",sep="\t").dropna()
			complete_set.set_index(complete_set.columns[0],inplace=True)
			ogwise_predictions=other.predictions.copy().set_index(["OG","metabolite"])
			protwise_predictions=self.predictions.set_index([self.protcol,"metabolite"])
			metabolites=protwise_predictions.index.get_level_values(1).unique()
			for og in ogwise_predictions.index:
				proteins=list(set(";".join(complete_set.loc[og[0]].tolist()).split(";")) & set(protwise_predictions.index.get_level_values(0)))
				locs=list(map(lambda x: (x,og[1]), proteins))
				if intersect=="all" and protwise_predictions.loc[locs,self.classifier].sum()==len(proteins):
					ogwise_predictions.loc[og,self.classifier]=True
				elif intersect=="majority" and protwise_predictions.loc[locs,self.classifier].sum()>=len(proteins)/2:
					ogwise_predictions.loc[og,self.classifier]=True
				elif intersect=="any" and protwise_predictions.loc[locs,self.classifier].sum()>0:
					ogwise_predictions.loc[og,self.classifier]=True
				else:
					ogwise_predictions.loc[og,self.classifier]=False
			ogwise_predictions.to_csv(self.analyses+"predictions_ogwise_"+intersect+trim+"_"+self.method.__name__+self.feature+self.approach+"_"+self.classifier+".tsv",sep="\t",index=True,header=True)
		return(ogwise_predictions)
	
	
	
	
	#creates graph object from predictions
	def create_graph_from_preds(self):
		edges=[tuple(x) for x in self.predictions[self.predictions.iloc[:,-1]==True].iloc[:,:-1].values]
		g=Graph.TupleList(edges,directed=False)
		#include vertices (OGs and metabolites) which have no interacting partners
		vertices_g=set(g.vs["name"])
		vertices_df=set(self.predictions[self.protcol].tolist()+self.predictions["metabolite"].tolist())
		remaining_vertices=vertices_df-vertices_g
		g.add_vertices(list(remaining_vertices))
		return(g)
	
	
	# plot degree distribution
	def density_of_graph(self):
		num_metabolites=len(self.predictions["metabolite"].unique())
		og_degree=list(filter(lambda x: x<=num_metabolites, self.g.degree()))
		num_ogs=len(og_degree)
		num_fewint_ogs=len(list(filter(lambda x: x<num_metabolites and x>0,og_degree)))
		manyint_ogs=len(list(filter(lambda x: x==num_metabolites, og_degree)))
		solo_ogs=len(list(filter(lambda x: x==0, og_degree)))
		print(str(round(num_fewint_ogs/num_ogs*100,2))+"% of "+self.protcol+"s interact with several but not all metabolites")
		print(str(round(manyint_ogs/num_ogs*100,2))+"% of "+self.protcol+"s interact with all metabolites")
		print(str(round(solo_ogs/num_ogs*100,2))+"% of "+self.protcol+"s interact with no metabolite")
		print("mean: "+str(np.mean(og_degree)))
		print("median: "+str(np.median(og_degree)))
		f=plt.figure(figsize=(10,10))
		plt.hist(og_degree)
		plt.title(self.feature+self.approach)
		plt.xlabel("degree")
		plt.ylabel("frequency")
		plt.show()
		f.savefig(self.analyses+"degree_distribution"+self.feature+self.approach+self.trim+"_"+self.method.__name__+".png")
	
	
	#powerlaw function
	def powerlaw(self,x,c,m,c0):
		return(c*(x**m)+c0)
	
	
	def ln_ln_powerlaw(self,x,m,const):
		return(m*x+const)
	
	
	# plot cumulative degree distribution
	#random={"ts","cs"}
	def degree_distribution(self,random=""):
		if random!="":
			self.predictions=pd.read_csv(self.analyses+"predictions_random_"+random+"_"+self.method.__name__+self.feature+self.approach+"_"+self.classifier+".tsv",sep="\t")[[self.protcol,"metabolite",self.classifier]]
			self.g=self.create_graph_from_preds()
			self.g.vs["type"]=list(map(lambda x: type(x)==int,self.g.vs["name"])) #1 for metabolites, 0 for proteins
		f=plt.figure(figsize=(10,10))
		num_prots=len(self.g.vs)-np.sum(self.g.vs["type"])
		num_metas=np.sum(self.g.vs["type"])
		##for metabolites and proteins together
		degrees_comb=self.g.degree()
		cum_probabilities_comb=list(map(lambda d: len(list(filter(lambda f: f>=d,degrees_comb)))/len(degrees_comb),range(num_prots+num_metas+1)))
		probabilities_comb=list(map(lambda d: len(list(filter(lambda f: f==d,degrees_comb)))/len(degrees_comb),range(num_prots+num_metas+1))) #,range(np.max(degrees_comb)+1)
		df_comb=pd.DataFrame(columns=["degree","probability","cumulative_probability"])
		df_comb["degree"]=range(num_prots+num_metas+1)#range(np.max(degrees_comb)+1)
		df_comb["probability"]=probabilities_comb
		df_comb["cumulative_probability"]=cum_probabilities_comb
		##for proteins only
		degrees_prot=list(filter(lambda x: x!=None,map(lambda i: self.g.degree()[i] if self.g.vs["type"][i]==False else None,range(self.g.vcount()))))
		cum_probabilities_prot=list(map(lambda d: len(list(filter(lambda f: f>=d,degrees_prot)))/num_prots,range(num_metas+1)))
		probabilities_prot=list(map(lambda d: len(list(filter(lambda f: f==d,degrees_prot)))/num_prots,range(num_metas+1)))
		df_prot=pd.DataFrame(columns=["degree","probability","cumulative_probability"])
		df_prot["degree"]=range(num_metas+1)#range(np.max(degrees_prot)+1)
		df_prot["probability"]=probabilities_prot
		df_prot["cumulative_probability"]=cum_probabilities_prot
		##for metabolites only
		degrees_met=list(filter(lambda x: x!=None,map(lambda i: self.g.degree()[i] if self.g.vs["type"][i]==True else None,range(self.g.vcount()))))
		probabilities_met=list(map(lambda d: len(list(filter(lambda f: f>=d,degrees_met)))/len(degrees_met),range(num_prots+1)))
		cum_probabilities_met=list(map(lambda d: len(list(filter(lambda f: f==d,degrees_met)))/len(degrees_met),range(num_prots+1)))
		df_met=pd.DataFrame(columns=["degree","probability","cumulative_probability"])
		df_met["degree"]=range(num_prots+1)
		df_met["probability"]=probabilities_met
		df_met["cumulative_probability"]=probabilities_met
		##fit
		popt_comb,pcov_comb=curve_fit(self.powerlaw,df_comb["degree"],df_comb["probability"])
		lnln_popt_comb,lnln_pcov_comb=curve_fit(self.ln_ln_powerlaw,df_comb["degree"],df_comb["probability"])
		df_comb["fit"]=self.powerlaw(df_comb["degree"],*popt_comb)
		popt_prot,pcov_prot=curve_fit(self.powerlaw,df_prot["degree"],df_prot["probability"])
		lnln_popt_prot,lnln_pcov_prot=curve_fit(self.ln_ln_powerlaw,df_prot["degree"],df_prot["probability"])
		df_prot["fit"]=self.powerlaw(df_prot["degree"],*popt_prot)
		popt_met,pcov_met=curve_fit(self.powerlaw,df_met["degree"],df_met["probability"])
		lnln_popt_met,lnln_pcov_met=curve_fit(self.ln_ln_powerlaw,df_met["degree"],df_met["probability"])
		df_met["fit"]=self.powerlaw(df_met["degree"],*popt_met)
		##plot proteins
		plt.plot(df_prot["degree"],df_prot["probability"].tolist(),"bo")
		plt.plot(df_prot["degree"],df_prot["fit"],"b-",label="power law fit with P(k > K) = %0.2f $K^{%0.2f}$+%0.2f" %(popt_prot[0],popt_prot[1],popt_prot[2]))
		#plot metabolites
		plt.plot(df_met["degree"],df_met["cumulative_probability"].tolist(),"rd")
		plt.plot(df_met["degree"],df_met["fit"],"r-",label="power law fit with P(k > K) = %0.2f $K^{%0.2f}$+%0.2f" %(popt_met[0],popt_met[1],popt_met[2]))
		plt.yscale("log")
		plt.xscale("log")
		plt.xlabel("Degree of nodes K")
		plt.ylabel("P(k > K)")
		plt.legend()
		plt.show()
		#save
		f.savefig(self.analyses+"degree_distribution"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".png")
		powerlaw_comb="P(k > K) = %0.2f $K^{%0.2f}$+%0.2f" %(popt_comb[0],popt_comb[1],popt_comb[2])
		powerlaw_prot="P(k > K) = %0.2f $K^{%0.2f}$+%0.2f" %(popt_prot[0],popt_prot[1],popt_prot[2])
		powerlaw_met="P(k > K) = %0.2f $K^{%0.2f}$+%0.2f" %(popt_met[0],popt_met[1],popt_met[2])
		powerlaw_comb_file=open(self.analyses+"powerlaw_func_comb"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".txt","w")
		powerlaw_comb_file.write(powerlaw_comb)
		powerlaw_prot_file=open(self.analyses+"powerlaw_func_prot"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".txt","w")
		powerlaw_prot_file.write(powerlaw_prot)
		powerlaw_met_file=open(self.analyses+"powerlaw_func_met"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".txt","w")
		powerlaw_met_file.write(powerlaw_met)
		df_comb.to_csv(self.analyses+"degree_distribution_comb"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".tsv",sep="\t")
		df_prot.to_csv(self.analyses+"degree_distribution_prot"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".tsv",sep="\t")
		df_met.to_csv(self.analyses+"degree_distribution_met"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".tsv",sep="\t")
		file_degrees_prot=open(self.analyses+"protein_degrees"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".txt","w")
		file_degrees_prot.write("\n".join(str(x) for x in degrees_prot))
		file_degrees_prot.close()
		file_degrees_met=open(self.analyses+"metabolite_degrees"+random+self.feature+self.approach+self.trim+"_"+self.method.__name__+"_"+self.classifier+".txt","w")
		file_degrees_met.write("\n".join(str(x) for x in degrees_met))
		file_degrees_met.close()
		self.predictions=self.load_predictions(trim=self.trim)
		self.create_graph_from_preds()
		self.g.vs["type"]=list(map(lambda x: type(x)==int,self.g.vs["name"])) #1 for metabolites, 0 for proteins
		return
	
	
	
	#compare two prediction sets, the added value of other to self
	def compare_predictions(self,other):
		intersect=self.predictions.set_index([self.protcol,"metabolite"]).join(other.predictions.set_index([self.protcol,"metabolite"]),lsuffix="_first",rsuffix="_second")
		if len(intersect.dropna())!=len(intersect):
			print("CAREFUL!!!")
			print("set of predicted data is not same for both classifiers")
			intersect=intersect.dropna()
			print("only taking intersection of both predictions into account, length of intersection: "+str(len(intersect)))
		intersect[intersect.columns[0]]=intersect[intersect.columns[0]].astype(int)
		intersect[intersect.columns[1]]=intersect[intersect.columns[1]].astype(int)
		intersect["difference"]=list(map(lambda i: intersect.iloc[i,0]-intersect.iloc[i,1],range(len(intersect)))) #intersect.diff(axis=1)
		#how often do the two prediction sets give same result (coincide)
		agreement=len(intersect[intersect["difference"]==0])/len(intersect)
		print(str(round(agreement*100,2))+"% of predictions among both sets coincide")
		#how many of the interactions found by the first classifier were also found by the second? = intersection of interactions
		intersecting_interactions=(intersect[intersect.columns[0]].sum()-len(intersect[intersect["difference"]==1]))/intersect[intersect.columns[0]].sum()
		print(str(round(intersecting_interactions*100,2))+"% of the interactions found by first set are verified by second")
		#how often does the first classifier says true but the second false
		second_contradicts=len(intersect[intersect["difference"]==1])/intersect[intersect.columns[0]].sum()
		print("For "+str(round(second_contradicts*100,2))+"% of interactions the second predictions contradict the first (the second says 'False' where the first says 'True')")
		#how much does the number of interactions increase by adding the second classifier
		#added_interactions=(intersect[intersect.columns[1]].sum()-intersect[intersect.columns[0]].sum())/intersect[intersect.columns[0]].sum()
		added_interactions=len(intersect[intersect["difference"]==-1])/intersect[intersect.columns[0]].sum()  #/intersecting_interactions
		print("The second prediction set added "+str(round(added_interactions*100,2))+"% new interactions")
		
	
	
	#plots venn diagramm for 2-6 sets
	#nas is list of type [self,other1,other2,...]
	#labels is list with names of sets
	def venn(self,nas,labels=None,appendix="",colors=None):
		#if prediction sets from different organisms, trim to same prediction sets (intersection between experiments)
		intersect_predictions=set(nas[0].predictions.set_index([self.protcol,"metabolite"]).index.tolist())
		#intersect_predictions=list(list(map(lambda na: set(set(na.predictions.set_index([self.protcol,"metabolite"]).index.tolist())&intersect_predictions),nas))[-1])
		for na in nas[1:]:
			intersect_predictions=intersect_predictions & set(na.predictions.set_index([na.protcol,"metabolite"]).index.tolist())
		intersect_predictions=list(intersect_predictions)
		ints=list()
		for na in nas:
			preds=na.predictions.set_index([self.protcol,"metabolite"])[na.classifier].loc[intersect_predictions]
			ints.append(set(preds[preds==True].index.tolist()))### preds[preds==True] ### preds[preds[preds.columns[0]]==True].index.tolist())
		#matplotlib_venn.venn3(ints,labels)
		data=venn.get_labels(ints)
		if colors is None:
			f,axn=eval("venn.venn"+str(len(nas))+"(data,names=labels)")
		else:
			f,axn=eval("venn.venn"+str(len(nas))+"(data,names=labels,colors=colors)")
		plt.show()
		f.savefig(self.analyses+"venn_diagram_"+"_".join(labels)+appendix+".png")
		f.savefig(self.analyses+"venn_diagram_"+"_".join(labels)+appendix+".svg")
		
	


###############################################################
###############################################################
###############################################################

#class to make analyses on experimental data and predictions
#includes network analyses on classes and construction of figures from thesis
class DataAnalysis(IDTranslations):
	
	def __init__(self,approach="",posset=None,negset=None,appendix="/",feature="",normalized="",proteinwise=False):
		self.feature=feature
		self.normalized=normalized
		self.approach=approach
		self.analyses="../analyses"+appendix
		self.experimental="../experimental_data"+appendix
		self.databases="../databases"+appendix
		self.expfiles=self.get_expfiles()
		self.proteinwise=proteinwise
		if proteinwise==False:
			self.protcol="OG"
		else:
			self.protcol="protein"
		
	
	#get files with experimental data
	def get_expfiles(self):
		expfiles=glob.glob(self.experimental+"*.xlsx") # if you want to predict data different from training data, adjust this path to prediction data 
		if expfiles==[]:
			root=Tk()
			root.withdraw()
			expfiles=glob.glob(askdirectory(initialdir="../")+"/*.xlsx")
		return(expfiles)
		
	
	################################
	#Analysis of experimental data
	################################
	
	#counts metabolites in file of experimental data (including ambiguous metabolites)
	def count_metabolites_file(self,expdata):
		return(len(expdata.get_metabolites()))
	
	#count proteins in file of experimental data (including ambiguous and not translated proteins)
	def count_proteins_file(self,expdata):
		return(len(expdata.get_proteins()))
	
	#counts metabolites in given xset
	def count_metabolites_set(self,xset):
		'''
		metabolites=set()
		for i in xset.index:
			for j in xset.loc[i,"CIDs"].split(";"):
				metabolites.add(j)
		'''
		try:
			metabolites=set(";".join(list(map(lambda x: str(x),xset["CIDs"].tolist()))).split(";"))
		except: 
			metabolites=set(";".join(list(map(lambda x: str(x),xset["metabolite"].tolist()))).split(";"))
		return(len(metabolites))
		
	#counts proteins in given xset
	def count_proteins_set(self,xset):
		if self.proteinwise==False:
			proteins=set()
			for expfilename in self.expfiles:
				organism=expfilename[:-5].split("/")[-1]
				for i in xset.index:
					if not pd.isna(xset.loc[i,organism]):
						for j in xset.loc[i,organism].split(";"):
							proteins.add(j)
		else:
			proteins=xset.reset_index()["protein"].unique()
		return(len(proteins))
	
	#counts proteins-meta pairs for protwise, og-meta-pairs for ogwise analysis
	def count_og_meta_pairs(self,xset):
		#return(len(list(itertools.chain.from_iterable(list(map(lambda x: x.split(";"),xset["CIDs"].tolist()))))))
		try:
			pairs=len(";".join(list(map(lambda x: str(x),xset["CIDs"].tolist()))).split(";"))
		except:
			pairs=len(";".join(list(map(lambda x: str(x),xset["metabolite"].tolist()))).split(";"))
		return(pairs)
	
	#counts protein-metabolite pairs for og-wise analysis
	def count_prot_meta_pairs(self,xset):
		if self.proteinwise==False:
			pairs=0
			for expfilename in self.expfiles:
				organism=expfilename[:-5].split("/")[-1]
				for i in xset.index:
					if not pd.isna(xset.loc[i,organism]):
						for j in xset.loc[i,organism].split(";"):
							pairs+=len(xset.loc[i,"CIDs"].split(";"))
		else:
			print("use self.count_og_meta_pairs(xset) to count protein-metabolite pairs for proteinwise analysis")
		return(pairs)
		
	
	#construct summary dataframe about all experimental datasets used
	def construct_df_files(self,method=np.min):
		organisms=list()
		for expfilename in self.expfiles:
			organisms.append(expfilename[:-5].split("/")[-1])
		cs=CompleteSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,normalized=self.normalized)
		complete_set=cs.load_set(cs.databases+"complete_set.tsv")
		complete_set_profiles=cs.load_set(cs.databases+"complete_set_"+cs.normalized+"normprofiles"+cs.feature+"_"+method.__name__+".tsv")
		df=pd.DataFrame(index=organisms+["intersection"],columns=["proteins","metabolites","orthologous groups"])
		proteins=set()
		metabolites=set()
		for expfilename in self.expfiles:
			organism=expfilename[:-5].split("/")[-1]
			expdata=FileHandler(expfilename)
			df.loc[organism,"proteins"]=self.count_proteins_file(expdata)
			df.loc[organism,"metabolites"]=self.count_metabolites_file(expdata)
			proteins=proteins|expdata.get_proteins()
			metabolites=metabolites|expdata.get_metabolites()
			if self.proteinwise==False:
				df.loc[organism,"orthologous groups"]=len(complete_set[organism].dropna())
		#df.ix["total"]=df.sum()
		df.loc["intersection","metabolites"]=len(complete_set_profiles["metabolite"].unique())
		df.loc["union","proteins"]=len(proteins)
		df.loc["union","metabolites"]=len(metabolites)
		#df.loc["intersection","metabolites"]=len(cs.create_complete_set_profiles().iloc[0]["CIDs"].split(";"))
		if self.proteinwise==False:
			df.loc["intersection","proteins"]=self.count_proteins_set(complete_set.loc[complete_set_profiles.index.unique()])
			df.loc["intersection","orthologous groups"]=len(complete_set_profiles.index.unique())
			df.loc["union","orthologous groups"]=len(complete_set)
		else:
			df.loc["intersection","proteins"]=self.count_proteins_set(complete_set_profiles)
			df.drop("orthologous groups",axis=1,inplace=True)
		return(df)
	
	
	#constructs summary dataframe for training sets used
	def construct_df_sets(self,method):
		df=pd.DataFrame(index=["positive set","negative set","training set"],columns=["proteins","metabolites","orthologous groups",self.protcol+"-metabolite-interaction-pairs"])
		ps=PositiveSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,normalized=self.normalized,normalized2="",proteinwise=self.proteinwise)
		ns=NegativeSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,normalized=self.normalized,normalized2="",proteinwise=self.proteinwise)
		if self.proteinwise==True:
			positive_set=ps.load_set(ps.databases+"positive_set"+ps.feature+self.approach+"_"+method.__name__+".tsv")
			negative_set=ns.load_set(self.databases+"negative_set"+self.feature+self.approach+"_"+method.__name__+".tsv")
		else:
			positive_set=ps.load_set(ps.databases+"positive_set"+ps.feature+self.approach+".tsv")
			negative_set=ns.load_set(ns.databases+"negative_set"+ns.feature+self.approach+".tsv")
		pos_neg_set=positive_set.append(negative_set)
		pos_neg_set.drop_duplicates(inplace=True)
		training_set=pd.read_csv(self.analyses+"training_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t")
		training_set=training_set.set_index([self.protcol,"metabolite"])
		sets=[positive_set,negative_set]
		for i in range(len(sets)):
			df.loc[df.index[i],"proteins"]=self.count_proteins_set(sets[i])
			df.loc[df.index[i],"metabolites"]=self.count_metabolites_set(sets[i])
			df.loc[df.index[i],"orthologous groups"]=len(sets[i])
			df.loc[df.index[i],self.protcol+"-metabolite-interaction-pairs"]=self.count_og_meta_pairs(sets[i])
			if self.proteinwise==False:
				df.loc[df.index[i],"protein-metabolite-interaction-pairs"]=self.count_prot_meta_pairs(sets[i])
		df.loc["training set","proteins"]=self.count_proteins_set(pos_neg_set.loc[training_set.index.get_level_values(0).unique()])
		df.loc["training set","metabolites"]=len(training_set.index.get_level_values(1).unique())
		if self.proteinwise==False:
			df.loc["training set","orthologous groups"]=len(training_set.index.get_level_values(0).unique())
			prot_met_pairs_ps_in_ts=self.count_prot_meta_pairs(positive_set.loc[training_set[training_set["interaction"]==True].index.get_level_values(0).unique()])
			prot_met_pairs_ns_in_ts=self.count_prot_meta_pairs(negative_set.loc[training_set[training_set["interaction"]==False].index.get_level_values(0).unique()])
			df.loc["training set","protein-metabolite-interaction-pairs"]=prot_met_pairs_ps_in_ts+prot_met_pairs_ns_in_ts
		else:
			df.drop("orthologous groups",axis=1,inplace=True)
		df.loc["training set",self.protcol+"-metabolite-interaction-pairs"]=len(training_set)
		return(df)
	
	
	
	
	#creates figure with following plots for different experiments (Ara_only,Eucaryotes,...):
	#1st column: histogram showing frequency of pearson correlation of protein profiles per OG
	#2nd column: pie chart with number of proteins per OG (how many OGs have 1 protein, how many two and forth)
	#3rd column: pie chart with number of OGs annotated to protein
	def corr_ogs_hist_pies(self,mlcs,method,profiles="raw"):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		colors=["tab:green","tab:blue","tab:red","tab:blue","tab:gray"]
		titles=["A.thaliana","S.cerevisiae"]
		#alphas=np.linspace(start=0.5,stop=1,num=len(xsets))
		alphas=np.linspace(1,0.3,6)
		f=plt.figure(figsize=(10,10))
		r=len(mlcs)*2
		for m,mlc in enumerate(mlcs):
			#create grid for histrogram and three pie charts for classifier
			ax1=plt.subplot2grid((r,3),(2*m,0),colspan=1,rowspan=2)
			ax2=plt.subplot2grid((r,3),(2*m,1))
			ax3=plt.subplot2grid((r,3),(2*m+1,1))
			ax4=plt.subplot2grid((r,3),(2*m,2),colspan=1,rowspan=2)
			#load dataframes with correlations between proteins in one OG
			try:
				correlations_df=pd.read_csv(mlc.analyses+"corrs_"+method.__name__+".tsv",sep="\t")
				num_prots=pd.read_csv(mlc.analyses+"number_of_proteins_per_OG.tsv",sep="\t")
			except FileNotFoundError:
				#read complete set and training set
				cs=CompleteSet(simulation=False,overwrite=False,databases=mlc.databases,analyses=mlc.analyses,experimental=mlc.experimental,normalized=self.normalized)
				complete_set=cs.load_set(mlc.databases+"complete_set.tsv")
				training_set=mlc.merge_pos_neg_orthos(method=method)
				correlations_df=pd.DataFrame()
				cs_correlationlist=list()
				ts_correlationlist=list()
				for i in range(len(mlc.expfiles)):
					expdata=FileHandler(mlc.expfiles[i],experimental=mlc.experimental,databases=mlc.databases)
					cs_org=complete_set[[expdata.organism]]
					cs_org.columns=["UniProt_IDs"]
					cs_org=cs_org.dropna()
					cs_correlations,cs_num_prots_org=expdata.profiles_correlation_ogs(cs_org,profiles=profiles,method=method)
					cs_correlationlist=cs_correlationlist+list(itertools.chain.from_iterable(cs_correlations["correlations"].tolist()))
					ts_org=training_set[[expdata.organism]]
					ts_org.columns=["UniProt_IDs"]
					ts_org=ts_org.dropna()
					ts_correlations,ts_num_prots_org=expdata.profiles_correlation_ogs(ts_org,profiles=profiles,method=method)
					ts_correlationlist=ts_correlationlist+list(itertools.chain.from_iterable(ts_correlations["correlations"].tolist()))
					if i==0:
						cs_num_prots=cs_num_prots_org
						ts_num_prots=ts_num_prots_org
					else:
						cs_num_prots=cs_num_prots.join(cs_num_prots_org)
						ts_num_prots=ts_num_prots.join(ts_num_prots_org) 
				correlations_df["complete_set"]=cs_correlationlist
				correlations_df.loc[:len(ts_correlationlist)-1,"training_set"]=ts_correlationlist
				correlations_df.to_csv(mlc.analyses+"corrs_"+method.__name__+".tsv",sep="\t")
				cs_num_prots=cs_num_prots.fillna(0).sort_index()
				ts_num_prots=ts_num_prots.fillna(0).sort_index()
				cs_num_prots["complete_set"]=cs_num_prots.sum(axis=1)
				ts_num_prots["training_set"]=ts_num_prots.sum(axis=1)
				num_prots=pd.concat([cs_num_prots["complete_set"],ts_num_prots["training_set"]],axis=1)
				num_prots.to_csv(mlc.analyses+"number_of_proteins_per_OG.tsv",index=True,header=True,sep="\t")
			#trim num_prots
			num_prots=num_prots.replace(np.nan,0)
			num_prots_trimmed=pd.DataFrame(columns=["complete_set","training_set"],index=list(map(lambda idx: str(idx),range(1,5))))
			num_prots_trimmed["complete_set"]=num_prots.loc[1:4,"complete_set"].tolist()
			num_prots_trimmed["training_set"]=num_prots.loc[1:4,"training_set"].tolist()
			num_prots_trimmed.loc["5-10","complete_set"]=num_prots.loc[5:10,"complete_set"].sum()
			num_prots_trimmed.loc[">10","complete_set"]=num_prots.loc[11:,"complete_set"].sum()
			num_prots_trimmed.loc["5-10","training_set"]=num_prots.loc[5:10,"training_set"].sum()
			num_prots_trimmed.loc[">10","training_set"]=num_prots.loc[11:,"training_set"].sum()
			num_prots_trimmed["alpha"]=alphas
			num_prots_trimmed=num_prots_trimmed.replace(0.0,np.nan)
			#first column: histogram
			#ax1.set_title(abc[m*4],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			ax1.hist(correlations_df["complete_set"],bins=20,histtype="bar", alpha=1,color=colors[2*m])
			ax1.hist(correlations_df["training_set"],bins=20,histtype="bar", alpha=1,color=colors[2*m+1])
			ax1.set_ylabel("Frequency")
			ax1.set_xlabel("Pearson correlation coefficient")
			#second column: pie charts #prots/OG
			complete_colors=[colors[2*m]]*6
			training_colors=[colors[2*m+1]]*6
			#ax2.set_title(abc[m*4+1],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			ax2.set_title(titles[m],fontweight="bold",fontstyle="italic")
			cn=ax2.pie(num_prots_trimmed["complete_set"].dropna().tolist(),labels=num_prots_trimmed["complete_set"].dropna().index.tolist(),colors=complete_colors[:len(num_prots_trimmed["complete_set"].dropna())]) #,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_prots_trimmed["complete_set"].sum()))))
			for i in range(len(cn[0])):
				cn[0][i].set_alpha(num_prots_trimmed.drop("training_set",axis=1).dropna().iloc[i,1])
			#ax3.set_title(abc[m*4+2],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			tn=ax3.pie(num_prots_trimmed["training_set"].dropna().tolist(),labels=num_prots_trimmed["training_set"].dropna().index.tolist(),colors=training_colors[:len(num_prots_trimmed["training_set"].dropna())]) #,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_prots_trimmed["training_set"].sum()))))
			for i in range(len(tn[0])):
				tn[0][i].set_alpha(num_prots_trimmed.drop("complete_set",axis=1).dropna().iloc[i,1])
			#third column: pie chart proteins with assigned OG
			num_ogs=mlc.count_og_assignments() 
			#ax4.set_title(abc[m*4+3],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			ax4.pie(num_ogs.sum(axis=1),labels=num_ogs.index,colors=[colors[4],colors[2*m]]) #,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_ogs.sum(axis=1).sum()))))
		f.tight_layout(rect=[0,0,0.9,1])
		#f.tight_layout()
		#f.subplots_adjust()
		plt.show()
		f.savefig(self.analyses+"corr_ogs_histograms_pies.png")
		
	
	#plot heatmaps for correlation over repetitions
	def heatmaps_repetitions(self):
		pass
	
	#creates figure with following plots for different experiments (Ara_only,Eucaryotes,...):
	#1st column: heatmaps with rv-coefficients for proteins and metabolites among repetitions
	#2nd column: histogram showing frequency of pearson correlation of protein profiles per OG
	#3rd column: pie chart with number of proteins per OG (how many OGs have 1 protein, how many two and forth)
	def corr_ogs_heat_hist_pies(self,mlcs,method,profiles="raw",first_col=True,sec_col=True,third_col=True,colors=None,titles=None):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		if colors is None:
			colors=["tab:green","tab:blue","tab:orange","tab:blue","tab:gray"]
		tableau_red=(214/255,39/255,40/255) #eucaryotes
		tableau_blue=(31/255,119/255,180/255) #metabolites
		tableau_green=(44/255,160/255,44/255) #arabidopsis
		tableau_orange=(255/255,127/255,14/255) #yeast
		tableau_yellow=(255/255,221/255,113/255)
		colormap1=LinearSegmentedColormap.from_list("tabblue",[tableau_yellow,tableau_blue],N=100)
		colormap2=LinearSegmentedColormap.from_list("tabred",[tableau_yellow,tableau_red],N=100)
		if titles is None:
			titles=["A.thaliana","S.cerevisiae"]
		org_dict={"A.thaliana_cell_cultures":"cell culture","A.thaliana_rosettes":"rosettes","A.thaliana_seedlings_4_rep":"seedlings","S.cerevisiae_cell_cultures":"cell culture","S.cerevisiae_growthphases_condition1":"log phase","S.cerevisiae_growthphases_condition2":"diauxic growth","S.cerevisiae_growthphases_condition3":"stat. phase"}
		#alphas=np.linspace(start=0.5,stop=1,num=len(xsets))
		alphas=np.linspace(1,0.3,6)
		#count number of rows for figure
		numfiles=0
		for mlc in mlcs:
			if len(mlc.expfiles)>2:
				numfiles+=len(mlc.expfiles)
			else:
				numfiles+=2
		r=2*numfiles+2
		axn1=list()
		#f,axn=plt.subplots(figsize=(20,10),sharex=True)
		f=plt.figure(figsize=(15,10))
		gridspec.GridSpec(nrows=r,ncols=6)#,width_ratios=[4,4,1]) #width ratio not working
		#axn.set_global()
		for m,mlc in enumerate(mlcs):
			if m==0:
				ax0a=plt.subplot2grid((r,6),(0,0),rowspan=1,colspan=2)
				ax0b=plt.subplot2grid((r,6),(1,0),rowspan=1,colspan=2)
				if first_col==False:
					ax0a.axis("off")
					ax0b.axis("off")
				axn1.append(ax0a)
				#axn1.append(ax0b)
				ax2=plt.subplot2grid((r,6),(0,2),colspan=3,rowspan=len(mlc.expfiles)*2+2) #histograms
				ax3=plt.subplot2grid((r,6),(0,5),rowspan=len(mlc.expfiles)+1) #pie expdata
				ax4=plt.subplot2grid((r,6),(len(mlc.expfiles)+1,5),rowspan=len(mlc.expfiles)+1)
			else:
				ax2=plt.subplot2grid((r,6),(len(axn1)*2,2),colspan=3,rowspan=len(mlc.expfiles)*2) #histograms
				ax3=plt.subplot2grid((r,6),(len(axn1)*2,5),colspan=1,rowspan=len(mlc.expfiles)) #pie expdata
				ax4=plt.subplot2grid((r,6),(len(axn1)*2+len(mlc.expfiles),5),rowspan=len(mlc.expfiles)) #pie training data
				#create grid for heatmaps, histrogram and pie charts for classifier
			last=len(axn1)
			for i in range(len(mlc.expfiles)):
				axn1.append(plt.subplot2grid((r,6),(len(axn1)*2,0),colspan=2,rowspan=2)) #heatmaps
			#load dataframes with heatmaps
			prots_mets=list()
			orgs=list()
			mlc.expfiles.sort()
			for i in range(len(mlc.expfiles)):
				expdata=FileHandler(mlc.expfiles[i],experimental=mlc.experimental,databases=mlc.databases)
				organism=mlc.expfiles[i][:-5].split("/")[-1]
				try:
					rv_prot=pd.read_csv(mlc.analyses+organism+"_protein_heatmap_repetitions.tsv",sep="\t").set_index("Unnamed: 0")
					rv_met=pd.read_csv(mlc.analyses+organism+"_metabolite_heatmap_repetitions.tsv",sep="\t").set_index("Unnamed: 0")
				except FileNotFoundError:
					rv_prot,rv_met=expdata.profiles_correlation_reps(profiles="raw")
					#saving correlations to "tsv"
					rv_prot.to_csv(mlc.analyses+organism+"_protein_heatmap_repetitions.tsv",index=True,header=True,sep="\t")
					rv_met.to_csv(mlc.analyses+organism+"_metabolite_heatmap_repetitions.tsv",index=True,header=True,sep="\t")
				#join them
				rv_prot=rv_prot.where(np.tril(np.ones(rv_prot.shape),k=-1).astype(np.bool)) #lower triangle belongs to protein
				rv_met=rv_met.where(np.triu(np.ones(rv_met.shape),k=1).astype(np.bool)) #upper triangle belongs to metabolite
				rv_org=rv_prot.multiply(rv_met,fill_value=1)
				prots_mets.append(rv_org)
				orgs.append(org_dict[organism])
			#load dataframes with correlations between proteins in one OG
			try:
				correlations_df=pd.read_csv(mlc.analyses+"corrs_"+method.__name__+".tsv",sep="\t")
				num_prots=pd.read_csv(mlc.analyses+"number_of_proteins_per_OG.tsv",sep="\t")
			except FileNotFoundError:
				#read complete set and training set
				cs=CompleteSet(simulation=False,overwrite=False,databases=mlc.databases,analyses=mlc.analyses,experimental=mlc.experimental,normalized=mlc.normalized)
				complete_set=cs.load_set(mlc.databases+"complete_set.tsv")
				training_set=mlc.merge_pos_neg_orthos(method=method)
				correlations_df=pd.DataFrame()
				cs_correlationlist=list()
				ts_correlationlist=list()
				for i in range(len(mlc.expfiles)):
					expdata=FileHandler(mlc.expfiles[i],experimental=mlc.experimental,databases=mlc.databases)
					cs_org=complete_set[[expdata.organism]]
					cs_org.columns=["UniProt_IDs"]
					cs_org=cs_org.dropna()
					cs_correlations,cs_num_prots_org=expdata.profiles_correlation_ogs(cs_org,profiles=profiles,method=method)
					cs_correlationlist=cs_correlationlist+list(itertools.chain.from_iterable(cs_correlations["correlations"].tolist()))
					ts_org=training_set[[expdata.organism]]
					ts_org.columns=["UniProt_IDs"]
					ts_org=ts_org.dropna()
					ts_correlations,ts_num_prots_org=expdata.profiles_correlation_ogs(ts_org,profiles=profiles,method=method)
					ts_correlationlist=ts_correlationlist+list(itertools.chain.from_iterable(ts_correlations["correlations"].tolist()))
					if i==0:
						cs_num_prots=cs_num_prots_org
						ts_num_prots=ts_num_prots_org
					else:
						cs_num_prots=cs_num_prots.join(cs_num_prots_org)
						ts_num_prots=ts_num_prots.join(ts_num_prots_org) 
				correlations_df["complete_set"]=cs_correlationlist
				correlations_df.loc[:len(ts_correlationlist)-1,"training_set"]=ts_correlationlist
				correlations_df.to_csv(mlc.analyses+"corrs_"+method.__name__+".tsv",sep="\t")
				cs_num_prots=cs_num_prots.fillna(0).sort_index()
				ts_num_prots=ts_num_prots.fillna(0).sort_index()
				cs_num_prots["complete_set"]=cs_num_prots.sum(axis=1)
				ts_num_prots["training_set"]=ts_num_prots.sum(axis=1)
				num_prots=pd.concat([cs_num_prots["complete_set"],ts_num_prots["training_set"]],axis=1)
				num_prots.to_csv(mlc.analyses+"number_of_proteins_per_OG.tsv",index=True,header=True,sep="\t")
			#trim num_prots
			num_prots=num_prots.replace(np.nan,0)
			num_prots_trimmed=pd.DataFrame(columns=["complete_set","training_set"],index=list(map(lambda idx: str(idx),range(1,5))))
			num_prots_trimmed["complete_set"]=num_prots.loc[1:4,"complete_set"].tolist()
			num_prots_trimmed["training_set"]=num_prots.loc[1:4,"training_set"].tolist()
			num_prots_trimmed.loc["5-10","complete_set"]=num_prots.loc[5:10,"complete_set"].sum()
			num_prots_trimmed.loc[">10","complete_set"]=num_prots.loc[11:,"complete_set"].sum()
			num_prots_trimmed.loc["5-10","training_set"]=num_prots.loc[5:10,"training_set"].sum()
			num_prots_trimmed.loc[">10","training_set"]=num_prots.loc[11:,"training_set"].sum()
			num_prots_trimmed["alpha"]=alphas
			num_prots_trimmed=num_prots_trimmed.replace(0.0,np.nan)
			#first column: heatmaps
			if first_col==True:
				for i,ax1 in enumerate(axn1[last:]):
					if i==0 and m==0:
						#cbar_ax=ax0
						cbar_ax_red=ax0a
						cbar_ax_blue=ax0b
						#from mpl_toolkits.axes_grid1 import make_axes_locatable
						#divider=make_axes_locatable(ax0)
						#cbar_ax_red=divider.append_axes("top",size="50%",pad=0.05)
						#cbar_ax_blue=divider.append_axes("bottom",size="50%",pad=0.05)
						#f.colorbar(ax0,cax=cbar_ax,orientation="horizontal")
					ax1.xaxis.set_ticks_position("none")
					ax1.yaxis.set_ticks_position("none")
					prots_mets[i].index.name=orgs[i]
					sns.heatmap(prots_mets[i].where(np.tril(np.ones(prots_mets[i].shape),k=-1).astype(np.bool)),ax=ax1,annot=False,xticklabels=False,yticklabels=False,cbar=i==0 and m==0,vmin=0,vmax=1,cbar_ax=None if i!=0 or m!=0 else cbar_ax_red,linewidths=0.1,cmap=colormap2,cbar_kws={"orientation":"horizontal"})#,"shrink":0.5}) #"anchor":(0.5,1.0),"panchor":(0.5,0.0)
					sns.heatmap(prots_mets[i].where(np.triu(np.ones(prots_mets[i].shape),k=1).astype(np.bool)),ax=ax1,annot=False,xticklabels=False,yticklabels=False,cbar=i==0 and m==0,vmin=0,vmax=1,cbar_ax=None if i!=0 or m!=0 else cbar_ax_blue,linewidths=0.1,cmap=colormap1,cbar_kws={"orientation":"horizontal"})#,"shrink":0.5}) #"anchor":(0.5,1.0),"panchor":(0.5,0.0)
					#f.colorbar(hm,cax=cbar_ax,orientation="horizontal")
					#second column: histogram
				#ax1.set_title(abc[m*4],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
			else:
				for ax1 in axn1[last:]:
					ax1.axis("off")
			if sec_col==True:
				ax2.hist(correlations_df["complete_set"],bins=20,histtype="bar", alpha=1,color=colors[2*m])
				ax2.hist(correlations_df["training_set"].dropna(),bins=20,histtype="bar", alpha=1,color=colors[2*m+1])
				ax2.set_ylabel("abs. frequency")
				if m==len(mlcs)-1:
					ax2.set_xlabel("Pearson correlation coefficient")
				#ax2.set_title(abc[m*4+1],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
				if titles[m]!="":
					ax2.set_title(titles[m],fontweight="bold",fontstyle="italic")
			else:
				ax2.axis("off")
			#third column: pie charts #prots/OG
			complete_colors=[colors[2*m]]*6
			training_colors=[colors[2*m+1]]*6
			if third_col==True:
				cn=ax3.pie(num_prots_trimmed["complete_set"].dropna().tolist(),labels=num_prots_trimmed["complete_set"].dropna().index.tolist(),colors=complete_colors[:len(num_prots_trimmed["complete_set"].dropna())]) #,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_prots_trimmed["complete_set"].sum()))))
				for i in range(len(cn[0])):
					cn[0][i].set_alpha(num_prots_trimmed.drop("training_set",axis=1).dropna().iloc[i,1])
				#ax3.set_title(abc[m*4+2],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
				tn=ax4.pie(num_prots_trimmed["training_set"].dropna().tolist(),labels=num_prots_trimmed["training_set"].dropna().index.tolist(),colors=training_colors[:len(num_prots_trimmed["training_set"].dropna())]) #,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_prots_trimmed["training_set"].sum()))))
				for i in range(len(tn[0])):
					tn[0][i].set_alpha(num_prots_trimmed.drop("complete_set",axis=1).dropna().iloc[i,1])
				#third column: pie chart proteins with assigned OG
				#num_ogs=mlc.count_og_assignments() 
				#ax4.set_title(abc[m*4+3],weight="bold",horizontalalignment="left",verticalalignment="center_baseline")
				#ax4.pie(num_ogs.sum(axis=1),labels=num_ogs.index,colors=[colors[4],colors[2*m]]) #,autopct=lambda p: "{:.0f}".format(int(round(p/100*num_ogs.sum(axis=1).sum()))))
			else:
				ax3.axis("off")
				ax4.axis("off")
		#f.tight_layout(rect=[0,0,0.9,1])
		#f.subplots_adjust()
		#plt.show()
		f.tight_layout()
		f.subplots_adjust(bottom=0.2,right=0.8)
		plt.show()
		app="_cols"
		if first_col==True:
			app=app+"1"
		if sec_col==True:
			app=app+"2"
		if third_col==True:
			app=app+"3"
		f.savefig("../overall_analysis/corr_ogs_heatmaps_histograms_pies"+app+".png")
		f.savefig("../overall_analysis/corr_ogs_heatmaps_histograms_pies"+app+".svg")
	
	
	########################################
	#Annotate pathways and classes
	########################################
	
	#dictionary for KEGG pathways
	def get_pathway_dict(self):
		#original
		'''
		pathway_dict={
		"00332":"Carbapenem biosynthesis",
		"00250":"Alanine, aspartate and glutamate metabolism",
		"04113":"Meiosis - yeast",
		"00190":"Oxidative phosphorylation",
		"00330":"Arginine and proline metabolism",
		"00380":"Tryptophan metabolism",
		"00410":"beta-Alanine metabolism",
		"00760":"Nicotinate and nicotinamide metabolism",
		"00750":"Vitamin B6 metabolism",
		"00270":"Cysteine and methionine metabolism",
		"01230":"Biosynthesis of amino acids",
		"00230":"Purine metabolism",
		"00740":"Riboflavin metabolism",
		"00770":"Panthotenate and CoA biosynthesis",
		"04213":"Longevity regulating pathway - multiple species",
		"00480":"Glutathione metabolism",
		"00400":"Phenylalanine, tyrosine and tryptophan biosynthesis",
		"00970":"Aminoacyl-tRNA biosynthesis",
		"00730":"Thiamine metabolism",
		"00350":"Tyrosine metabolism",
		"00680":"Methane metabolism",
		"00130":"Ubiquinone and other terpenoid-quinone biosynthesis",
		"01100":"Metabolic pathways",
		"01110":"Biosynthesis of secondary metabolites",
		"01210":"2-Oxocarboxylic acid metabolism",
		"00460":"Cyanoamino acid metabolism",
		"04111":"Cell cycle -yeast",
		"00220":"Arginine biosynthesis",
		"00240":"Pyrimidine metabolism",
		"00260":"Glycine, serine and threonine metabolism",
		"00261":"Monobactam biosynthesis",
		"00910":"Nitrogen metabolism"}
		'''
		#abbreviations
		pathway_dict={
		"00195":"Photosynthesis",
		"00332":"Carbapenem b.",
		"00250":"Ala, Asp, Glu m.",
		"04113":"Meiosis - yeast",
		"00190":"Ox. phosphorylation",
		"00330":"Arg and Pro m.",
		"00380":"Trp m.",
		"00410":"ß-Ala m.",
		"00760":"Nicotinate and NAM m.",
		"00750":"Vitamin B6 m.",
		"00270":"Cys and Met m.",
		"01230":"B. of amino acids",
		"00230":"Purine m.",
		"00740":"Riboflavin m.",
		"00770":"Panthotenate and CoA b.",
		"04213":"Longevity regulating",
		"00480":"Glutathione m.",
		"00400":"Phe, Tyr and Trp b.",
		"00970":"Aminoacyl-tRNA b.",
		"00730":"Thiamine m.",
		"00350":"Tyr m.",
		"00680":"Methane m.",
		"00130":"Ubiquinone, terp.-qu. b.",
		"01100":"Metabolic pathways",
		"01110":"B. of sec. metabolites",
		"01210":"2-Oxocarboxylic acid m.",
		"00460":"Cyanoamino acid m.",
		"04111":"Cell cycle",
		"00220":"Arg b.",
		"00240":"Pyrimidine m.",
		"00260":"Gly, Ser and Thr m.",
		"00261":"Monobactam b.",
		"00910":"Nitrogen m.",
		"00360":"Phenylalanine m.",
		"00908":"Zeatin b.",
		"00940":"Phenylpropanoid b.",
		"00966":"Glucosinolate b.",
		"00780":"Biotin m.",
		"00960":"Tropane, piperidine, pyridine alkaloid b.",
		"00950":"Isoquinoline alkaloid b."}
		return(pathway_dict)
	
	
	
	#creates a dataframe which counts interactions and non-interactions between pathways
	#exclude_pathways is set of pathways to be excluded (too general terms)
	def make_pathway_dataframe(self,method,classifier):
		#load pathways
		protein_kegg_ec_df=pd.read_csv(self.analyses+"UniProt_BRITEclasses_KEGGpathways.tsv",sep="\t").set_index("UniProt ID")
		metabolite_kegg_df=pd.read_csv("../databases/metabolite_kegg_classes_pathways.tsv",sep="\t").set_index("CID")
		pathway_dict=self.get_pathway_dict()
		#load predictions
		predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")[[self.protcol,"metabolite",classifier]].set_index([self.protcol,"metabolite"])
		predictions.columns=["prediction"]
		##predictions=predictions[predictions["prediction"]==True]
		#load complete_set
		if self.proteinwise==False:
			cs=CompleteSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature,normalized=self.normalized,proteinwise=self.proteinwise)
			complete_set=cs.load_set(cs.databases+"complete_set.tsv").dropna()
		#intersecting pathways
		metabolite_pathways=set(map(lambda x: x[8:],";".join(metabolite_kegg_df["KEGG pathways"].dropna()).split(";")))
		protein_pathways=set(map(lambda x: x[8:],";".join(protein_kegg_ec_df["KEGG pathways"].dropna()).split(";")))
		intersect_pathways=set(metabolite_pathways & protein_pathways)
		intersect_pathway_names=list(map(lambda x: pathway_dict[x],intersect_pathways))
		#initialize pathway dataframe 
		pathway_df=pd.DataFrame(columns=["protein pathway","metabolite pathway","number of interactions","number of non-interactions"])
		prot_pathways=len(intersect_pathway_names)*intersect_pathway_names
		prot_pathways.sort()
		pathway_df["protein pathway"]=prot_pathways
		pathway_df["metabolite pathway"]=len(intersect_pathway_names)*intersect_pathway_names
		pathway_df["number of interactions"]=0
		pathway_df["number of non-interactions"]=0
		pathway_df.set_index(["protein pathway","metabolite pathway"],inplace=True)
		pathway_df_trimmed=pathway_df.copy() #trimmed case excludes proteins and metabolites without interactions (important for significance)
		#append pathways to predictions
		predictions["protein pathways"]=protein_kegg_ec_df.loc[predictions.index.get_level_values(0),"KEGG pathways"].tolist()
		predictions["metabolite pathways"]=metabolite_kegg_df.loc[predictions.index.get_level_values(1),"KEGG pathways"].tolist()
		predictions.dropna(inplace=True)
		for i in predictions.index:
			for prot_pathway in list(map(lambda x: pathway_dict[x[8:]] if x[8:] in pathway_dict else x,predictions.loc[i,"protein pathways"].split(";"))):
				if prot_pathway in intersect_pathway_names:
					for meta_pathway in list(map(lambda x: pathway_dict[x[8:]] if x[8:] in pathway_dict else x,predictions.loc[i,"metabolite pathways"].split(";"))):
						if meta_pathway in intersect_pathway_names:
							if self.proteinwise==True:
								pathway_df.loc[(prot_pathway,meta_pathway),"number of interactions"]+=predictions.loc[i,"prediction"]
								#trimmed case
								if predictions.xs(key=i[0],level=0)["prediction"].sum()!=0 or predictions.xs(key=i[1],level=1)["prediction"].sum()!=0:
									pathway_df_trimmed.loc[(prot_pathway,meta_pathway),"number of non-interactions"]+=abs(predictions.loc[i,"prediction"]-1)
								pathway_df.loc[(prot_pathway,meta_pathway),"number of non-interactions"]+=abs(predictions.loc[i,"prediction"]-1)
							else:
								if predictions.loc[i,"prediction"]==True:
									pathway_df.loc[(prot_pathway,meta_pathway),"number of interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
								else:
									pathway_df.loc[(prot_pathway,meta_pathway),"number of non-interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
									#trimmed case:
									if predictions.xs(key=i[0],level=0)["prediction"].sum()!=0 or predictions.xs(key=i[1],level=1)["prediction"].sum()!=0:
										pathway_df_trimmed.loc[(prot_pathway,meta_pathway),"number of non-interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
		pathway_df_trimmed["number of interactions"]=pathway_df["number of interactions"]
		pathway_df.reset_index(inplace=True)
		pathway_df_trimmed.reset_index(inplace=True)
		pathway_df.to_csv(self.analyses+"interactions_kegg_pathways_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t",index=False,header=True)
		pathway_df_trimmed.to_csv(self.analyses+"interactions_kegg_pathways_trimmed_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t",index=False,header=True)
		return(pathway_df)
	
	
	
	#dictionary for BRITE classes
	def get_brite_dict(self):
		brite_dict={"00535":"Proteoglycans",
		"00536":"Glycosaminoglycan binding proteins",
		"00537":"GPI-anchored proteins",
		"01000":"Enzymes",
		"01001":"Protein kinases",
		"01002":"Peptidases and inhibitors",
		"01003":"Glycosyltransferases",
		"01004":"Lipid biosynthesis proteins",
		"01005":"Lipopolysaccharide biosyntesis proteins",
		"01006":"Prenyltransferases",
		"01007":"Amino acid related enzymes",
		"01009":"Protein phosphatases and associated proteins",
		"00194":"Photosynthesis proteins",
		"00199":"Cytochrome P450",
		"02000":"Transporters",
		"02022":"Two-component system",
		"02044":"Secretion system",
		"02048":"Prokaryotic defense system",
		"03000":"Transcription factors",
		"03009":"Ribosome biogenesis",
		"03011":"Ribosome",
		"03012":"Translation factors",
		"03016":"Transfer RNA biogenesis",
		"03019":"Messenger RNA biogenesis",
		"03021":"Transcription machinery",
		"03029":"Mitochondrial biogenesis",
		"03032":"DNA replication proteins",
		"03036":"Chromosome and associated proteins",
		"03041":"Spliceosome",
		"03051":"Proteasome",
		"03110":"Chaperones and folding catalysts",
		"03400":"DNA repair and recombination proteins",
		"04031":"GTP-binding proteins",
		"04040":"Ion channels",
		"04090":"CD molecules",
		"04091":"Lectins",
		"04121":"Ubiquitin system",
		"04131":"Membrane trafficking",
		"04147":"Exosome",
		"04812":"Cytoskeleton proteins"}
		return(brite_dict)
		
	
	
	#create dataframes with #interactions among classes (EC and BRITE)
	#better not give predictions and method but approach and feature, or deepcopy predictions beforehand
	#random={"complete","training"}
	def make_graph_data_among_classes(self,classifier,predictions_raw=None,complete_set=None,method=np.mean,approach=None,feature=None,random=False):
		if feature is None:
			feature=self.feature
		if approach is None:
			approach=self.approach
		if predictions_raw is None:
			if random=="complete":
				predictions_raw=pd.read_csv(self.analyses+"predictions_random_cs_"+method.__name__+feature+approach+"_"+classifier+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
				predictions_raw.rename(columns={classifier:"prediction"},inplace=True)
			elif random=="training":
				predictions_raw=pd.read_csv(self.analyses+"predictions_random_ts_"+method.__name__+feature+approach+"_"+classifier+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
				predictions_raw.rename(columns={classifier:"prediction"},inplace=True)
			else: #if not random
				na=NetworkAnalysis(appendix=self.databases[12:],method=method,approach=approach,feature=feature,trim="",classifier=classifier)
				predictions_raw=na.predictions.set_index([self.protcol,"metabolite"])
				predictions_raw.columns=["prediction"]
		if complete_set is None:
			if self.proteinwise==True:
				try:
					complete_set=pd.read_csv(self.analyses+"UniProt_BRITEclasses_KEGGpathways.tsv",sep="\t",dtype={"EC class":str}).set_index("UniProt ID")
				except FileNotFoundError:
					print("assign classes to proteins with XSet.classify_proteins() and XSet.assign_EC_numbers_to_proteins()")
					return
			else:
				try:
					complete_set=pd.read_csv(self.analyses+"OG_UniProt_BRITEclasses_KEGGpathways.tsv",sep="\t",dtype={"intersect EC class":str}).set_index("Unnamed: 0").replace(np.nan,"")
				except FileNotFoundError:
					print("assign classes to OGs")
					#cs=CompleteSet(simulation=False,overwrite=False,experimental=self.experimental,databases=self.databases,analyses=self.analyses,feature=self.feature)
					#complete_set=cs.load_set(cs.databases+"complete_set.tsv").dropna()
					return
		metabolite_kegg_df=pd.read_csv("../databases/metabolite_kegg_classes_pathways.tsv",sep="\t").set_index("CID")
		#antipredictions=predictions[predictions["prediction"]==False]
		#predictions=predictions[predictions["prediction"]==True]
		#ANNOTATION OF CLASSES TO PREDICTIONS
		try:
			if random=="complete":
				predictions=pd.read_csv(self.analyses+"predictions_random_cs_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",sep="\t",dtype={"EC":str}).set_index([self.protcol,"metabolite"])
			elif random=="training":
				predictions=pd.read_csv(self.analyses+"predictions_random_ts_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",sep="\t",dtype={"EC":str}).set_index([self.protcol,"metabolite"])
			else: #not random
				predictions=pd.read_csv(self.analyses+"predictions_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",sep="\t",dtype={"EC":str}).set_index([self.protcol,"metabolite"])
			print("annotated predictions found and loaded")
		except FileNotFoundError:
			predictions=predictions_raw.copy()
			if self.proteinwise==False:
				for i in predictions.index:
					predictions.loc[i,"BRITE classes"]=complete_set.loc[i[0],"intersect BRITE classes"]
					predictions.loc[i,"EC"]=complete_set.loc[i[0],"intersect EC class"]
					predictions.loc[i,"metabolite class"]=metabolite_kegg_df.loc[i[1],"My class"]
			else:
				for i in predictions.index:
					predictions.loc[i,"BRITE classes"]=complete_set.loc[i[0],"BRITE classes"]
					predictions.loc[i,"EC"]=complete_set.loc[i[0],"EC class"]
					predictions.loc[i,"metabolite class"]=metabolite_kegg_df.loc[i[1],"My class"]
			if random=="complete":
				predictions.to_csv(self.analyses+"predictions_random_cs_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",index=True,header=True,sep="\t")
			elif random=="training":
				predictions.to_csv(self.analyses+"predictions_random_ts_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",index=True,header=True,sep="\t")
			else: #not random
				predictions.to_csv(self.analyses+"predictions_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",index=True,header=True,sep="\t")
		predictions.replace(np.nan,"",inplace=True)
		#all classes
		brite_classes=set(";".join(predictions["BRITE classes"]).split(";"))
		ec_numbers=set(";".join(predictions["EC"]).split(";"))
		metabolite_classes=set(";".join(predictions["metabolite class"]).split(";"))
		#dataframe for brite classes
		brite_df=pd.DataFrame(columns=["BRITE class","metabolite class","number of interactions","number of non-interactions"])
		m_brite_classes=len(metabolite_classes)*list(brite_classes)
		m_brite_classes.sort()
		brite_df["BRITE class"]=m_brite_classes
		brite_df["metabolite class"]=len(brite_classes)*list(metabolite_classes)
		brite_df=brite_df.set_index(["BRITE class","metabolite class"])
		brite_df["number of interactions"]=0
		brite_df["number of non-interactions"]=0
		brite_df_trimmed=brite_df.copy() #excludes metabolites and proteins without any found interactions
		#dataframe for ec numbers
		ec_df=pd.DataFrame(columns=["EC number","metabolite class","number of interactions","number of non-interactions"])
		m_ec_numbers=len(metabolite_classes)*list(ec_numbers)
		m_ec_numbers.sort()
		ec_df["EC number"]=m_ec_numbers
		ec_df["metabolite class"]=len(ec_numbers)*list(metabolite_classes)
		ec_df=ec_df.set_index(["EC number","metabolite class"])
		ec_df["number of interactions"]=0
		ec_df["number of non-interactions"]=0
		ec_df_trimmed=ec_df.copy()
		#count pairs
		for i in predictions.index:
			if predictions.loc[i,"prediction"]==True:
				if not pd.isna(predictions.loc[i,"BRITE classes"]):
					for brite_class in predictions.loc[i,"BRITE classes"].split(";"):
						if self.proteinwise==True:
							brite_df.loc[(brite_class,predictions.loc[i,"metabolite class"]),"number of interactions"]+=1
						else:
							brite_df.loc[(brite_class,predictions.loc[i,"metabolite class"]),"number of interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
				if not pd.isna(predictions.loc[i,"EC"]):
					if self.proteinwise==True:
						ec_df.loc[(predictions.loc[i,"EC"],predictions.loc[i,"metabolite class"]),"number of interactions"]+=1
					else:
						ec_df.loc[(predictions.loc[i,"EC"],predictions.loc[i,"metabolite class"]),"number of interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
			elif predictions.loc[i,"prediction"]==False:
				if not pd.isna(predictions.loc[i,"BRITE classes"]):
					for brite_class in predictions.loc[i,"BRITE classes"].split(";"):
						if self.proteinwise==True:
							brite_df.loc[(brite_class,predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=1
						else:
							brite_df.loc[(brite_class,predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
					#trimmed case:
					if predictions.xs(key=i[0],level=0)["prediction"].sum()!=0 or predictions.xs(key=i[1],level=1)["prediction"].sum()!=0:
						for brite_class in predictions.loc[i,"BRITE classes"].split(";"):
							if self.proteinwise==True:
								brite_df_trimmed.loc[(brite_class,predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=1
							else:
								brite_df_trimmed.loc[(brite_class,predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
				if not pd.isna(predictions.loc[i,"EC"]):
					if self.proteinwise==True:
						ec_df.loc[(predictions.loc[i,"EC"],predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=1
					else:
						ec_df.loc[(predictions.loc[i,"EC"],predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
					#trimmed case
					if predictions.xs(key=i[0],level=0)["prediction"].sum()!=0 or predictions.xs(key=i[1],level=1)["prediction"].sum()!=0:
						if self.proteinwise==True:
							ec_df_trimmed.loc[(predictions.loc[i,"EC"],predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=1
						else:
							ec_df_trimmed.loc[(predictions.loc[i,"EC"],predictions.loc[i,"metabolite class"]),"number of non-interactions"]+=len(complete_set.loc[i[0],"union UniProt IDs"].split(";"))
		brite_df_trimmed["number of interactions"]=brite_df["number of interactions"]
		ec_df_trimmed["number of interactions"]=ec_df["number of interactions"]
		brite_df=brite_df.reset_index()
		ec_df=ec_df.reset_index()
		brite_df_trimmed.reset_index(inplace=True)
		ec_df_trimmed.reset_index(inplace=True)
		if random=="complete":
			brite_df.to_csv(self.analyses+"random_cs_interactions_briteclass_metabolite_class_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			ec_df.to_csv(self.analyses+"random_cs_interactions_EC_metabolite_class_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			brite_df_trimmed.to_csv(self.analyses+"random_cs_interactions_briteclass_metabolite_class_trimmed_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			ec_df_trimmed.to_csv(self.analyses+"random_cs_interactions_EC_metabolite_class_trimmed_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
		elif random=="training":
			brite_df.to_csv(self.analyses+"random_ts_interactions_briteclass_metabolite_class_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			ec_df.to_csv(self.analyses+"random_ts_interactions_EC_metabolite_class_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			brite_df_trimmed.to_csv(self.analyses+"random_ts_interactions_briteclass_metabolite_class_trimmed_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			ec_df_trimmed.to_csv(self.analyses+"random_ts_interactions_EC_metabolite_class_trimmed_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
		else:
			brite_df.to_csv(self.analyses+"interactions_briteclass_metabolite_class_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			ec_df.to_csv(self.analyses+"interactions_EC_metabolite_class_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			brite_df_trimmed.to_csv(self.analyses+"interactions_briteclass_metabolite_class_trimmed_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
			ec_df_trimmed.to_csv(self.analyses+"interactions_EC_metabolite_class_trimmed_"+method.__name__+feature+approach+".tsv",sep="\t",index=False,header=True)
		return([brite_df,ec_df])
		
		
		
	#count number of interactions for every metabolite class
	def count_interactions_per_metabolite_class(self,predictions):
		metabolite_kegg_df=pd.read_csv("../databases/metabolite_kegg_classes_pathways.tsv",sep="\t").set_index("CID")
		predictions["metabolite_class"]=list(map(lambda x: metabolite_kegg_df.loc[x[1],"My class"], predictions.index))
		df=pd.DataFrame(index=predictions["metabolite_class"].unique(),columns=["num_interactions","num_metas","num_prots"])
		for i in df.index:
			df.loc[i,"num_interactions"]=predictions[predictions["metabolite_class"]==i]["prediction"].sum()
			df.loc[i,"num_metas"]=len(predictions[predictions["metabolite_class"]==i].index.get_level_values(1).unique())
			subpredictions=predictions[predictions["metabolite_class"]==i]
			df.loc[i,"num_prots"]=len(subpredictions[subpredictions["prediction"]==True].index.get_level_values(0).unique())
		df.to_csv(self.analyses+"interactions_proteins_metaboliteclass"+self.feature+self.approach+".tsv",sep="\t",header=True,index=True)
		return
	
	
	
	#count number of proteins and metabolites per class for complete and training set, for given DataAnalysis instances
	def count_class_assignments(self,das,colors,approaches,names,method=np.mean,feature=None,classifier="svm_clf(kernel='linear')"):
		if feature is None:
			feature=self.feature
		f,axn=plt.subplots(figsize=(10,10),nrows=len(das),ncols=2,sharex=False)
		#f.canvas.draw()
		for d,da in enumerate(das):
			#for prediction set
			ec_names={"":"Other Proteins", "1":"Oxidorecutases","2":"Transferases","3":"Hydrolases","4":"Lyases","5":"Isomerases","6":"Ligases","7":"Translocases"}
			predictions=pd.read_csv(da.analyses+"predictions_annotated_"+method.__name__+feature+approaches[d]+"_"+classifier+".tsv",sep="\t",dtype={"EC":str}).set_index([da.protcol,"metabolite"]).replace(np.nan,"")
			predictions=predictions.replace({"EC":ec_names})
			if d==0:
				#construct dataframe for plot
				meta_df=pd.DataFrame(index=predictions["metabolite class"].unique(),columns=names+list(map(lambda x: x+"_training",names)))
				og_df=pd.DataFrame(index=predictions["EC"].unique(),columns=names+list(map(lambda x: x+"_training",names)))
			num_ecs=Counter(predictions["EC"])
			num_mcs=Counter(predictions["metabolite class"])
			totnum_ogs=len(predictions.index.get_level_values(0).unique())
			totnum_metas=len(predictions.index.get_level_values(1).unique())
			for i in num_ecs:
				#num_ecs[i]=num_ecs[i]/totnum_metas
				og_df.loc[i,names[d]]=num_ecs[i]/totnum_metas
			for i in num_mcs:
				#num_mcs[i]=num_mcs[i]/totnum_ogs
				meta_df.loc[i,names[d]]=num_mcs[i]/totnum_ogs
			#for training set
			mlc=MLClassifiers(experimental=da.experimental,databases=da.databases,analyses=da.analyses,approach=approaches[d],feature=feature,normalized=da.normalized,proteinwise=da.proteinwise)
			training_set=mlc.load_known_set(method=method)
			ts_num_ecs=Counter(predictions.reset_index().set_index(da.protcol).loc[training_set.index.get_level_values(0).unique(),"EC"])
			ts_num_mcs=Counter(predictions.reset_index().set_index("metabolite").loc[training_set.index.get_level_values(1).unique(),"metabolite class"])
			if np.nan in ts_num_ecs:
				del ts_num_ecs[np.nan]
			if np.nan in ts_num_mcs:
				del ts_num_mcs[np.nan]
			for i in ts_num_ecs:
				#ts_num_ecs[i]=ts_num_ecs[i]/totnum_metas
				og_df.loc[i,names[d]+"_training"]=ts_num_ecs[i]/totnum_metas
			for i in ts_num_mcs:
				#ts_num_mcs[i]=ts_num_mcs[i]/totnum_ogs
				meta_df.loc[i,names[d]+"_training"]=ts_num_mcs[i]/totnum_ogs
			'''
			ec_lists=sorted(num_ecs.items())
			mc_lists=sorted(num_mcs.items())
			ts_ec_lists=sorted(ts_num_ecs.items())
			ts_mc_lists=sorted(ts_num_mcs.items())
			x_ecs,y_ecs=zip(*ec_lists)
			x_mcs,y_mcs=zip(*mc_lists)
			x_tsecs,y_tsecs=zip(*ts_ec_lists)
			x_tsmcs,y_tsmcs=zip(*ts_mc_lists)
			axn[d][0].bar(x_ecs,y_ecs,color=colors[d])
			axn[d][0].bar(x_tsecs,y_tsecs,color="tab:blue")
			axn[d][0].set_ylabel("absolute frequency")
			axn[d][1].bar(x_mcs,y_mcs,color=colors[d],alpha=0.7)
			axn[d][1].bar(x_tsmcs,y_tsmcs,color="tab:blue",alpha=0.7)
			'''
			axn[d][0].bar(og_df.index,og_df[names[d]],color=colors[d])
			axn[d][0].bar(og_df.index,og_df[names[d]+"_training"],color="tab:blue")
			if len(das)%2==0:
				axn[d][0].set_ylabel("abs. frequency")
			else:
				if d==np.median(range(len(das))):
					axn[d][0].set_ylabel("abs. frequency")
			axn[d][1].bar(meta_df.index,meta_df[names[d]],color=colors[d],alpha=0.7)
			axn[d][1].bar(meta_df.index,meta_df[names[d]+"_training"],color="tab:blue",alpha=0.7)
			if d==len(das)-1: #xticklabels for last plot
				axn[d][0].set_xticks(np.arange(len(og_df.index)),og_df.index.tolist())
				axn[d][0].set_xticklabels(og_df.index.tolist(),rotation=45,ha="right")
				axn[d][0].set_xlabel("EC class")
				axn[d][1].set_xticks(np.arange(len(meta_df.index)),meta_df.index.tolist())
				axn[d][1].set_xticklabels(meta_df.index.tolist(),rotation=45,ha="right")
				axn[d][1].set_xlabel("Metabolite class")
			else: #no xticklabels
				axn[d][0].set_xticklabels(len(og_df)*[""])
				axn[d][1].set_xticklabels(len(meta_df)*[""])
		#axn[d][1].set_xticks(np.arange(len(meta_df.index)),meta_df.index.tolist())
		#axn[d][1].set_xticklabels(meta_df.index.tolist(),rotation=45,ha="right")
		#f.text(0.02,0.4,"abs. frequency",ha="center",va="center",rotation="vertical")
		f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/"+self.protcol+"_class_distribution.png")
		f.savefig("../overall_analysis/"+self.protcol+"_class_distribution.svg")
		meta_df.to_csv("../overall_analysis/metabolite_"+self.protcol+"_class_distribution.tsv",sep="\t",index=True,header=True)
		og_df.to_csv("../overall_analysis/"+self.protcol+"_class_distribution.tsv",sep="\t",index=True,header=True)
		return
	
	
	
	#all class assignments merged over experiments
	def class_assignments_merged(self,color1=(44/255,160/255,44/255),color2=(255/255,127/255,14/255),color3=(214/255,39/255,40/255)):
		try:
			meta_df=pd.read_csv("../overall_analysis/metabolite_"+self.protcol+"_class_distribution.tsv",sep="\t").set_index("Unnamed: 0")
			og_df=pd.read_csv("../overall_analysis/"+self.protcol+"_class_distribution.tsv",sep="\t").set_index("Unnamed: 0")
		except FileNotFoundError:
			print("generate class distribution dataframes with self.count_class_assignments()")
			return
		if "Euc" in og_df.columns: #og_df["Euc"].sum()!=0: #if with both
			train_label=None
		else:
			train_label="Gold standard"
		f,axn=plt.subplots(figsize=(10,10),nrows=1,ncols=2)
		xlocs=np.arange(len(og_df))
		bar_width=0.25
		spacing=1.1
		axn[0].bar(xlocs-spacing*bar_width,og_df["Ath"],bar_width,label="A.thaliana",color=color1)
		axn[0].bar(xlocs-spacing*bar_width,og_df["Ath_training"],bar_width,label=None,color="tab:blue",alpha=1)
		axn[0].bar(xlocs,og_df["Sce"],bar_width,label="S.cerevisisae",color=color2)
		axn[0].bar(xlocs,og_df["Sce_training"],bar_width,label=None,color="tab:blue",alpha=1)
		axn[0].set_xticks(xlocs)
		axn[0].set_xticklabels(og_df.index.tolist(),rotation=45,ha="right")
		axn[0].set_xlabel("EC class")
		axn[0].set_ylabel("Abs. frequency")
		
		axn[1].bar(xlocs-spacing*bar_width,meta_df["Ath"],bar_width,label="A.thaliana",color=color1)
		axn[1].bar(xlocs-spacing*bar_width,meta_df["Ath_training"],bar_width,label=None,color="tab:blue",alpha=1)
		axn[1].bar(xlocs,meta_df["Sce"],bar_width,label="S.cerevisiae",color=color2)
		axn[1].bar(xlocs,meta_df["Sce_training"],bar_width,label=train_label,color="tab:blue",alpha=1)
		axn[1].set_xticks(xlocs)
		axn[1].set_xticklabels(meta_df.index.tolist(),rotation=45,ha="right")
		axn[1].set_xlabel("Metabolite class")
		
		if "Euc" in og_df.columns: #og_df["Euc"].sum()!=0:
			axn[0].bar(xlocs+spacing*bar_width,og_df["Euc"],bar_width,label="Both",color=color3)
			axn[0].bar(xlocs+spacing*bar_width,og_df["Euc_training"],bar_width,label=None,color="tab:blue",alpha=1)
			axn[1].bar(xlocs+spacing*bar_width,meta_df["Euc"],bar_width,label="Both",color=color3)
			axn[1].bar(xlocs+spacing*bar_width,meta_df["Euc_training"],bar_width,label="Gold standard",color="tab:blue",alpha=1)
		axn[1].legend(loc="upper right")
		f.tight_layout()
		plt.subplots_adjust(bottom=0.4)
		plt.show()
		f.savefig("../overall_analysis/"+self.protcol+"_class_distribution_merged.png")
		f.savefig("../overall_analysis/"+self.protcol+"_class_distribution_merged.svg")
	
	
	
	##########################################
	#Analysis of pathway and classes networks
	##########################################
	'''
		brite_edges=[tuple(brite_df.loc[x]) for x in brite_df.index]
		brite_g=Graph.TupleList(brite_edges,directed=False,weights=True)
		sum_ints=brite_df["number of interactions"].sum() #sum of all interactions
		num_ints={}
		#for brite classes
		brite_df=brite_df.set_index("BRITE class")
		for bc in brite_df.index:
			num_ints[bc]=brite_df.loc[bc,"number of interactions"].sum()
		#for Metabolite classes
		brite_df=brite_df.reset_index().set_index("metabolite class")
		for mc in brite_df.index:
			num_ints[mc]=brite_df.loc[mc,"number of interactions"].sum()
		#setting vertex size on number of interactions relative to sum of all interactions
		brite_df=brite_df.reset_index()
		for vertex in brite_g.vs:
			vertex["vertex_size"]=num_ints[vertex["name"]]
		brite_g.vs["vertex_size"]=list(map(lambda x: x/sum_ints*200,brite_g.vs["vertex_size"]))
		#name brite classes #########stopped here
		brite_names={"":"No Enzyme", "1":"Oxidorecutases","2":"Transferases","3":"Hydrolases","4":"Lyases","5":"Isomerases","6":"Ligases","7":"Translocases"}
		#to be continued....
	'''
	
	
	#create weighted EC graph
	#random={"complete","training"}
	def make_EC_graph(self,method=np.mean,approach="_unbalanced_low_confidence",classifier="svm_clf(kernel='linear')",feature=None,layout=None,color1=None,color2=None,random=False,font_size=14,plot_all_vertices=False,trim=""):
		if feature is None:
			feature=self.feature
		size_x=1000
		size_y=1000
		margin=60
		labeldistance=5
		initial_edge_width=0.1
		vertex_size_scale=0.2
		vertex_type1="circle"
		vertex_type2="circle" #"rectangle"
		if color1 is None:
			color1=(31/255,119/255,180/255) #metabolites
		if color2 is None:
			color2=(214/255,39/255,40/255) #proteins
		grey=(65/255,68/255,81/255)
		font="DejaVu Sans"
		sig_level=0.05
		ec_names={"":"Other Proteins", "1":"Oxidorecutases","2":"Transferases","3":"Hydrolases","4":"Lyases","5":"Isomerases","6":"Ligases","7":"Translocases"}
		inv_ec_names={v:k for k,v in ec_names.items()}
		#load number of interactions between classes
		try:
			if random=="complete":
				ec_df=pd.read_csv(self.analyses+"random_cs_interactions_EC_metabolite_class"+trim+"_"+method.__name__+feature+approach+".tsv",sep="\t",dtype={"EC number":str})
			elif random=="training":
				ec_df=pd.read_csv(self.analyses+"random_ts_interactions_EC_metabolite_class"+trim+"_"+method.__name__+feature+approach+".tsv",sep="\t",dtype={"EC number":str})
			else:
				ec_df=pd.read_csv(self.analyses+"interactions_EC_metabolite_class"+trim+"_"+method.__name__+feature+approach+".tsv",sep="\t",dtype={"EC number":str})
			ec_df=ec_df.replace(np.nan,"")
		except FileNotFoundError:
			print("creating dataframe with interactions between protein and metabolite classes")
			#trim not taken into consideration
			if random==False:
				brite_df,ec_df=self.make_graph_data_among_classes(method=method,approach=approach,classifier=classifier,feature=feature)
			else:
				brite_df,ec_df=self.make_graph_data_among_classes(method=method,approach=approach,classifier=classifier,feature=feature,random=random)
		#construct weighted graph
		g_ec_df=ec_df.drop("number of non-interactions",axis=1)
		ec_edges=[tuple(g_ec_df.loc[x]) for x in g_ec_df.replace(0,np.nan).dropna().index]
		ec_g=Graph.TupleList(ec_edges,directed=False,weights=True)
		if type(layout) is pd.core.frame.DataFrame:
			#rename Proteins to EC numbers
			layout_classes=list(map(lambda x: inv_ec_names[x] if x in inv_ec_names else x,layout.index.tolist()))
			missing_vertices=list(set(layout_classes)-set(ec_g.vs["name"]))
			if plot_all_vertices==True:
				ec_g.add_vertices(missing_vertices)
			else: layout.drop(missing_vertices,axis=0,inplace=True)
		else:
			missing_vertices=list()
		ec_g.vs["type"]=list(map(lambda x: len(x)>1,ec_g.vs["name"])) #1 for metabolites, 0 for proteins
		#rename ec classes
		ec_df=ec_df.replace({"EC number":ec_names})
		for vertex in ec_g.vs:
			if vertex["name"] in ec_names:
				vertex["name"]=ec_names[vertex["name"]]
		
		##create dictionary with number of interactions per vertex
		sum_ints=ec_df["number of interactions"].sum() #sum of all interactions
		sum_nonints=ec_df["number of non-interactions"].sum()
		num_ints={}
		num_nonints={}
		#for EC
		ec_df=ec_df.set_index("EC number")
		for ec in ec_df.index:
			num_ints[ec]=ec_df.loc[ec,"number of interactions"].sum()
			num_nonints[ec]=ec_df.loc[ec,"number of non-interactions"].sum()
		#for Metabolite classes
		ec_df=ec_df.reset_index().set_index("metabolite class")
		for mc in ec_df.index:
			num_ints[mc]=ec_df.loc[mc,"number of interactions"].sum()
			num_nonints[mc]=ec_df.loc[mc,"number of non-interactions"].sum()
		for missing_class in missing_vertices:
			num_ints[missing_class]=0
		#setting vertex size on number of interactions relative to sum of all interactions
		ec_df=ec_df.reset_index()
		for vertex in ec_g.vs:
			vertex["vertex_size"]=num_ints[vertex["name"]]
		ec_g.vs["vertex_size"]=list(map(lambda x: x/sum_ints*vertex_size_scale*size_x,ec_g.vs["vertex_size"]))
		#setting vertex shape
		for vertex in ec_g.vs:
			vertex["vertex_shape"]=vertex_type1 if vertex["type"] else vertex_type2
		#setting edge width 
		#for edge in ec_g.es:
		#	edge["edge_width"]=edge["weight"]/sum_ints*100
		
		#setting color, the darker the more instances in class
		# count class occurences:
		predictions=pd.read_csv(self.analyses+"predictions_annotated_"+method.__name__+feature+approach+"_"+classifier+".tsv",sep="\t",dtype={"EC":str}).set_index([self.protcol,"metabolite"]).replace(np.nan,"")
		predictions=predictions.replace({"EC":ec_names})
		num_classes=Counter(predictions["EC"])+Counter(predictions["metabolite class"])
		totnum_ogs=len(predictions.index.get_level_values(0).unique())
		print("number OGs:\t"+str(totnum_ogs))
		totnum_metas=len(predictions.index.get_level_values(1).unique())
		print("number metabolites:\t"+str(totnum_metas))
		color_dict={True: color1, False: color2}
		#color_dict={True: "blue", False:"red"}
		#ec_g.vs["vertex_color"]=[color_dict[moltype] for moltype in ec_g.vs["type"]]
		colormap1=LinearSegmentedColormap.from_list("tabblue",[color1,grey],N=100)
		colormap2=LinearSegmentedColormap.from_list("tabred",[color2,grey],N=100)
		colormap_dict={True: colormap1, False:colormap2}
		#set color due to colormap on number of instances
		for vertex in ec_g.vs:
			colorratio=num_classes[vertex["name"]]/len(predictions)
			vertex["vertex_color"]=colormap_dict[vertex["type"]](colorratio)
		
		#setting plotting parameters
		#setting plot arguments
		visual_style={}
		visual_style["bbox"]=(size_x,size_y)
		visual_style["margin"]=margin
		if type(layout) is pd.core.frame.DataFrame:
			#check if vertex order is same, if not bring in same order
			if ec_g.vs["name"]!=layout.index.tolist():
				sorter=dict(zip(ec_g.vs["name"],range(len(ec_g.vs["name"]))))
				layout["rank"]=layout.index.map(sorter)
				layout.sort_values("rank",axis=0,inplace=True)
				layout.drop("rank",axis=1,inplace=True)
			if ec_g.vs["name"]!=layout.index.tolist():
				print("vertex order in given layout still does not correspond vertex order of interaction dataframe!")
			visual_style["layout"]=[tuple(layout.loc[x]) for x in layout.index]
		else:
			if layout is None:
				visual_style["layout"]=ec_g.layout() #"large"
			else:
				visual_style["layout"]=ec_g.layout(layout)
			visual_style["layout"].fit_into((visual_style["margin"],visual_style["margin"],size_x-visual_style["margin"],size_y-visual_style["margin"]),keep_aspect_ratio=False)
		#visual_style["vertex_label"]=ec_g.vs["name"]
		visual_style["vertex_color"]=ec_g.vs["vertex_color"] 
		visual_style["vertex_size"]=ec_g.vs["vertex_size"]
		visual_style["vertex_shape"]=ec_g.vs["vertex_shape"]
		#visual_style["vertex_label_dist"]=1.2
		#visual_style["edge_width"]=ec_g.es["edge_width"]
		visual_style["edge_width"]=initial_edge_width
		visual_style["vertex_frame_color"]=visual_style["vertex_color"]
		
		#get coordinates of vertices
		if type(layout) is not pd.core.frame.DataFrame:
			x,y=np.array(visual_style["layout"].coords).T
		else: 
			x=layout["x"].tolist()
			y=layout["y"].tolist()
		coordinate_df=pd.DataFrame(index=ec_g.vs["name"],columns=["x","y"])
		coordinate_df["x"]=x
		coordinate_df["y"]=y
		
		#make dataframe with all classes against all classes
		#ec_df=ec_df.reset_index().set_index(["EC number","metabolite class"])
		ec_df_1=ec_df.copy()
		ec_df_1.columns=["class_1","class_2","number of interactions","number of non-interactions"]
		ec_df_2=ec_df_1.copy()
		ec_df_2.columns=["class_2","class_1","number of interactions","number of non-interactions"]
		class_df=ec_df_1.set_index(["class_1","class_2"]).append(ec_df_2.set_index(["class_1","class_2"]))
		
		p=plot(ec_g,**visual_style,opacity=1)
		ctx=cairo.Context(p.surface)
		#context=drawing.coord.CoordinateSystem
		for edge in ec_g.es:
			#height=np.sqrt((x[edge.source]-x[edge.target])**2+(y[edge.source]-y[edge.target])**2)/2
			#coordinates of central point B on edge (trapez variant with B1 and B2 for less peaky triangles)
			center_x=x[edge.source]+0.5*(x[edge.target]-x[edge.source])
			center_y=y[edge.source]+0.5*(y[edge.target]-y[edge.source])
			B=(center_x,center_y)
			width_source=class_df.loc[(ec_g.vs[edge.source]["name"],ec_g.vs[edge.target]["name"]),"number of interactions"]/num_ints[ec_g.vs[edge.source]["name"]]*ec_g.vs[edge.source]["vertex_size"]
			width_target=class_df.loc[(ec_g.vs[edge.target]["name"],ec_g.vs[edge.source]["name"]),"number of interactions"]/num_ints[ec_g.vs[edge.target]["name"]]*ec_g.vs[edge.target]["vertex_size"]
			###draw triangle as conical shaped edge from source to target
			#calculate coordinates of points on triangle baseline
			perpend_vec=(1,-(x[edge.target]-x[edge.source])/(y[edge.target]-y[edge.source])*1)
			norm_perpend_vec=np.sqrt(perpend_vec[0]**2+perpend_vec[1]**2)
			perpend_vec_normlen_source=(perpend_vec[0]/norm_perpend_vec*0.5*width_source,perpend_vec[1]/norm_perpend_vec*0.5*width_source)
			#check for perpendicularity (==0): perpend_vec_normlen[0]*(x[edge.target]-x[edge.source])+perpend_vec_normlen[1]*(y[edge.target]-y[edge.source])
			A=(x[edge.source]+perpend_vec_normlen_source[0],y[edge.source]+perpend_vec_normlen_source[1])
			C=(x[edge.source]-perpend_vec_normlen_source[0],y[edge.source]-perpend_vec_normlen_source[1])
			#check for perpendicularity (==0): (C[0]-A[0])*(x[edge.target]-x[edge.source])+(C[1]-A[1])*(y[edge.target]-y[edge.source])
			#trapez variant for less peaky triangle at B
			if width_source>2 and width_target>2:
				b_width=2
			else:
				b_width=width_source
			B1=(B[0]+b_width/2*perpend_vec[0]/norm_perpend_vec,B[1]+b_width/2*perpend_vec[1]/norm_perpend_vec)
			B2=(B[0]-b_width/2*perpend_vec[0]/norm_perpend_vec,B[1]-b_width/2*perpend_vec[1]/norm_perpend_vec)
			#draw triangle
			ctx.move_to(A[0],A[1])
			ctx.line_to(B1[0],B1[1])
			ctx.line_to(B2[0],B2[1])
			ctx.line_to(C[0],C[1])
			ctx.line_to(A[0],A[1])
			ctx.close_path()
			color_source=ec_g.vs[edge.source]["vertex_color"] #color1 if ec_g.vs[edge.source]["type"]==1 else color2
			ctx.set_source_rgb(color_source[0],color_source[1],color_source[2])
			ctx.fill()
			###draw triangle as conical shaped edge from target to source
			#calculate coordinates of points on triangle baseline
			perpend_vec_normlen_target=(perpend_vec[0]/norm_perpend_vec*0.5*width_target,perpend_vec[1]/norm_perpend_vec*0.5*width_target)
			D=(x[edge.target]+perpend_vec_normlen_target[0],y[edge.target]+perpend_vec_normlen_target[1])
			E=(x[edge.target]-perpend_vec_normlen_target[0],y[edge.target]-perpend_vec_normlen_target[1])
			#draw triangle
			ctx.move_to(D[0],D[1])
			#ctx.line_to(B[0],B[1])
			ctx.line_to(B1[0],B1[1])
			ctx.line_to(B2[0],B2[1])
			ctx.line_to(E[0],E[1])
			ctx.line_to(D[0],D[1])
			ctx.close_path()
			color_target=ec_g.vs[edge.target]["vertex_color"] #color1 if ec_g.vs[edge.target]["type"]==1 else color2
			ctx.set_source_rgb(color_target[0],color_target[1],color_target[2])
			ctx.fill()
			#if significantly many interactions from source to target, draw * on triangle
			ints_source_target=class_df.loc[(ec_g.vs[edge.source]["name"],ec_g.vs[edge.target]["name"]),"number of interactions"]
			nonints_source_target=class_df.loc[(ec_g.vs[edge.source]["name"],ec_g.vs[edge.target]["name"]),"number of non-interactions"]
			contingency_table_source=[[ints_source_target,nonints_source_target],[num_ints[ec_g.vs[edge.source]["name"]]-ints_source_target,num_nonints[ec_g.vs[edge.source]["name"]]-nonints_source_target]]
			oddsratio_source,pvalue_source=fisher_exact(contingency_table_source,alternative="greater")
			if pvalue_source<sig_level*1/len(ec_g.es): #bonferroni correction for multiple testing
				ctx.set_source_rgb(0,0,0)
				##draw *
				#quarter_x=x[edge.source]+0.25*(x[edge.target]-x[edge.source])
				#quarter_y=y[edge.source]+0.25*(y[edge.target]-y[edge.source])
				#ctx.select_font_face(font)
				#ctx.set_font_size(font_size)
				#ctx.move_to(quarter_x,quarter_y)
				#ctx.show_text(u"\u2732")
				##draw black edge
				norm_edge=np.sqrt((x[edge.target]-x[edge.source])**2+(y[edge.target]-y[edge.source])**2)
				normlen_edge=((x[edge.target]-x[edge.source])/norm_edge,(y[edge.target]-y[edge.source])/norm_edge)
				v_source_edge=(x[edge.source]+ec_g.vs[edge.source]["vertex_size"]*0.5*normlen_edge[0],y[edge.source]+ec_g.vs[edge.source]["vertex_size"]*0.5*normlen_edge[1])
				#ctx.move_to(x[edge.source],y[edge.source]) #from center of vertex
				ctx.move_to(v_source_edge[0],v_source_edge[1]) #from edge of vertex
				#line:
				ctx.line_to(B[0],B[1])
				#ctx.set_line_width(0.08)
				ctx.stroke()
			#if significantly many interactions from target to source, draw * on triangle
			ints_target_source=class_df.loc[(ec_g.vs[edge.target]["name"],ec_g.vs[edge.source]["name"]),"number of interactions"]
			nonints_target_source=class_df.loc[(ec_g.vs[edge.target]["name"],ec_g.vs[edge.source]["name"]),"number of non-interactions"]
			contingency_table_target=[[ints_target_source,nonints_target_source],[num_ints[ec_g.vs[edge.target]["name"]]-ints_target_source,num_nonints[ec_g.vs[edge.target]["name"]]-nonints_target_source]]
			oddsratio_target,pvalue_target=fisher_exact(contingency_table_target,alternative="greater")
			if pvalue_target<sig_level*1/len(ec_g.es):
				ctx.set_source_rgb(0,0,0)
				##draw *
				#quarter_x=x[edge.target]+0.25*(x[edge.source]-x[edge.target])
				#quarter_y=y[edge.target]+0.25*(y[edge.source]-y[edge.target])
				#ctx.select_font_face(font)
				#ctx.set_font_size(font_size)
				#ctx.move_to(quarter_x,quarter_y)
				#ctx.show_text(u"\u2732")
				##draw black edge
				norm_edge=np.sqrt((x[edge.source]-x[edge.target])**2+(y[edge.source]-y[edge.target])**2)
				normlen_edge=((x[edge.source]-x[edge.target])/norm_edge,(y[edge.source]-y[edge.target])/norm_edge)
				v_target_edge=(x[edge.target]+ec_g.vs[edge.target]["vertex_size"]*0.5*normlen_edge[0],y[edge.target]+ec_g.vs[edge.target]["vertex_size"]*0.5*normlen_edge[1])
				#ctx.move_to(x[edge.target],y[edge.target]) #center of vertex
				ctx.move_to(v_target_edge[0],v_target_edge[1]) #edge of vertex
				ctx.line_to(B[0],B[1])
				#ctx.set_line_width(0.08)
				ctx.stroke()
		#add vertex labels
		for v,vertex in enumerate(ec_g.vs):
			ctx.set_source_rgb(0,0,0)
			ctx.select_font_face(font)
			ctx.set_font_size(font_size)
			(x_l,y_l,width_l,height_l,dx_l,dy_l)=ctx.text_extents(vertex["name"])
			x_label=x[v]-width_l/2
			if y[v]<=size_y/2: #labels above node for upper half, below for lower half
				y_label=y[v]-vertex["vertex_size"]/2-labeldistance #labeldistance pixels above vertex
			else:
				y_label=y[v]+vertex["vertex_size"]/2+labeldistance+height_l #labeldistance pixels below vertex
			ctx.move_to(x_label,y_label)
			ctx.show_text(vertex["name"])
		#display and save
		p.show()
		if random==False:
			p.save(self.analyses+"EC_MC_network"+trim+feature+approach+".png")
			p.save(self.analyses+"EC_MC_network"+trim+feature+approach+".svg")
			coordinate_df.to_csv(self.analyses+"EC_MC_network"+trim+feature+approach+"_coordinates.tsv",sep="\t",header=True,index=True)
		else:
			p.save(self.analyses+"EC_MC_network_"+random+"_randomized"+trim+feature+approach+".png")
			p.save(self.analyses+"EC_MC_network_"+random+"_randomized"+trim+feature+approach+".svg")
			coordinate_df.to_csv(self.analyses+"EC_MC_network_"+random+"_randomized"+trim+feature+approach+"_coordinates.tsv",sep="\t",header=True,index=True)
		return
	
	
	
	#saves given colormap to a colorbar plot
	def save_cbar(self,colormap,savename,vmax=1):
		f,ax=plt.subplots()
		cbar=plt.cm.ScalarMappable(cmap=colormap,norm=plt.Normalize(vmin=0,vmax=vmax))
		cbar._A=[]
		plt.colorbar(cbar)
		#cbar2=plt.cm.ScalarMappable(cmap=colormap2,norm=plt.Normalize(vmin=0,vmax=totnum_ogs))
		#cbar2._A=[]
		#plt.colorbar(cbar2)
		ax.remove()
		f.savefig("../overall_analysis/figures/comparison_random/"+savename+".png")
		f.savefig("../overall_analysis/figures/comparison_random/"+savename+".svg")
		plt.close()
		return
	
	
	
	######################################
	#Comparison of various classifiers 
	######################################
	
	#creates figure to compare classifiers on their learning and ROC curves, degree distributions and Venn Diagrams of predictions.
	# predictions need to be made beforehand as well as degree distributions with NetworkAnalysis.py
	# experiment={"_Ara_only2/","_eucaryotes/",...}
	# train_sizes is list, same as used for MLClassifiers.create_learning_curves()
	#no trimming included yet!
	
	def classifier_comparison(self,experiments,features,approaches,classifiers,methods,appendix=None,colour=None,train_sizes=[0.3,0.4,0.5,0.6,0.7,0.8],legend_roc=True,legend_dd=True,title=False):
		#initialize figure
		num_rows=len(experiments)*len(features)*len(approaches)*len(classifiers)*len(methods)
		f,axn=plt.subplots(nrows=num_rows,ncols=3)
		r=-1 #:= row of figure
		if colour is None:
			colour=["tab:red","tab:blue","tab:green","tab:orange","tab:purple","tab:brown","k"]
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		#get maximal degrees for degree distributions:
		upper_ylim_mets=0
		upper_ylim_prots=0
		for experiment in experiments:
			for feature in features:
				for approach in approaches:
					for classifier in classifiers:
						for method in methods:
							prot_degrees_norm=open("../analyses"+experiment+"protein_degrees"+feature+approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
							met_degrees_norm=open("../analyses"+experiment+"metabolite_degrees"+feature+approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
							prot_degrees=[int(x) for x in prot_degrees_norm]
							met_degrees=[int(x) for x in met_degrees_norm]
							upper_ylim_mets=np.max([upper_ylim_mets,np.max(met_degrees)])
							upper_ylim_prots=np.max([upper_ylim_prots,np.max(prot_degrees)])
		for experiment in experiments:
			for feature in features:
				for approach in approaches:
					#1. read learning curves data
					lc_df=pd.read_csv("../analyses"+experiment+"learning_curves"+feature+approach+".tsv",sep="\t")
					lc_df=lc_df.set_index(["classifiers","methods"])
					#2. read ROC curve data
					roc_df=pd.read_csv("../analyses"+experiment+"ROC_curve"+feature+approach+".tsv",sep="\t")
					roc_df=roc_df.set_index(["classifiers","methods"])
					for classifier in classifiers:
						for method in methods:
							r+=1
							#plot learning curve: first column
							axn.flat[r*3].fill_between(train_sizes,list(map(lambda x: eval(x),lc_df.loc[(classifier,method.__name__),"train_score_min"].split(";"))),list(map(lambda x: eval(x),lc_df.loc[(classifier,method.__name__),"train_score_max"].split(";"))),alpha=0.1,color="k")
							axn.flat[r*3].fill_between(train_sizes,list(map(lambda x: eval(x),lc_df.loc[(classifier,method.__name__),"test_score_min"].split(";"))),list(map(lambda x: eval(x),lc_df.loc[(classifier,method.__name__),"test_score_max"].split(";"))),alpha=0.1,color=colour[r])
							axn.flat[r*3].plot(train_sizes,list(map(lambda x: eval(x),lc_df.loc[(classifier,method.__name__),"train_score_mean"].split(";"))),color="k",label="training score")
							axn.flat[r*3].plot(train_sizes,list(map(lambda x: eval(x),lc_df.loc[(classifier,method.__name__),"test_score_mean"].split(";"))),color=colour[r],label="test score")
							if title:
								axn.flat[r*3].set_title(abc[r*3],weight="bold")
							axn.flat[r*3].set_ylabel("accuracy")
							#plot ROC curve: second column
							axn.flat[r*3+1].plot([0,1],[0,1],color="k",linestyle="--")
							#axn.flat[r*3+1].plot([0,roc_df.loc[(classifier,method.__name__),"false_positive_rate"][0]],[0,roc_df.loc[(classifier,method),"true_positive_rate"][0]],color=colours[r])
							axn.flat[r*3+1].plot(list(map(lambda x: eval(x),roc_df.loc[(classifier,method.__name__),"false_positive_rate"].split(";"))),list(map(lambda x: eval(x),roc_df.loc[(classifier,method.__name__),"true_positive_rate"].split(";"))),colour[r],label="(AUC = %0.2f)" % roc_df.loc[(classifier,method.__name__),"AUC"])
							if title:
								axn.flat[r*3+1].set_title(abc[r*3+1],weight="bold")
							axn.flat[r*3+1].set_ylabel("tpr")
							if legend_roc:
								axn.flat[r*3+1].legend()
							#3. read and plot degree distributions: third column
							prot_degrees_norm=open("../analyses"+experiment+"protein_degrees"+feature+approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
							met_degrees_norm=open("../analyses"+experiment+"metabolite_degrees"+feature+approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
							prot_degrees=[int(x) for x in prot_degrees_norm]
							met_degrees=[int(x) for x in met_degrees_norm]
							bins=np.linspace(0,upper_ylim_mets,100)
							'''
							ddp_df=pd.read_csv("../analyses"+experiment+"degree_distribution_prot"+feature+approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
							ddm_df=pd.read_csv("../analyses"+experiment+"degree_distribution_met"+feature+approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
							eq_prot=open("../analyses"+experiment+"powerlaw_func_prot"+feature+approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
							eq_met=open("../analyses"+experiment+"powerlaw_func_met"+feature+approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
							axn.flat[r*3+2].plot(ddp_df["degree"],ddp_df["cumulative_probability"].tolist(),"o",color=colour[r])
							axn.flat[r*3+2].plot(ddm_df["degree"],ddm_df["cumulative_probability"].tolist(),".",color=colour[r],alpha=0.1)
							axn.flat[r*3+2].plot(ddp_df["degree"],ddp_df["fit"],linestyle="-",color=colour[r],label=eq_prot)
							axn.flat[r*3+2].plot(ddm_df["degree"],ddm_df["fit"],linestyle="-",color=colour[r],alpha=0.1,label=eq_met)
							if title:
								axn.flat[r*3+2].set_title(abc[r*3+2],weight="bold")
							axn.flat[r*3+2].set_ylabel("P(k > K)")
							axn.flat[r*3+2].set_yscale("log")
							axn.flat[r*3+2].set_xscale("log")
							'''
							axn.flat[r*3+2].hist(prot_degrees,color=colour[r],alpha=0.5)
							axn.flat[r*3+2].set_ylabel("abs. freq.",color=colour[r])
							#axn.flat[r*3+2].set_ylim(0,upper_ylim_mets)
							axn.flat[r*3+2].tick_params(axis="y",color=colour[r],labelcolor=colour[r])
							axn.flat[r*3+2].tick_params(axis="x",color=colour[r])
							axn.flat[r*3+2].set_xlim(0,upper_ylim_prots)
							ax1=axn.flat[r*3+2].twinx()
							ax2=ax1.twiny()
							ax2.hist(met_degrees,color="tab:blue",alpha=0.5)
							#ax1.set_ylim(0,upper_ylim_prots)
							ax1.tick_params(axis="y",color="tab:blue",labelcolor="tab:blue")
							ax2.set_xlim(0,upper_ylim_mets)
							ax2.tick_params(axis="x",color="tab:blue",labelcolor="tab:blue")
							ax1.set_ylabel("abs. freq.",color="tab:blue")
							if r<num_rows-1:
								axn.flat[r*3+2].set_xticklabels([])
							if r!=0:
								ax2.set_xticklabels([],color="tab:blue")
							if legend_dd:
								axn.flat[r*3+2].legend()
		axn.flat[r*3].set_xlabel("size of training set")
		axn.flat[r*3+1].set_xlabel("fpr")
		axn.flat[r*3+2].set_xlabel("Degree of nodes K")
		plt.subplots_adjust(hspace=0.4,wspace=0.3) #adjust vertical space between subplots,default=0.2
		#f.tight_layout()
		plt.show()
		if appendix is None:
			appendix="_"+"_".join(experiments+features+approaches+classifiers+list(map(lambda m: m.__name__,methods))).replace("/","")
		f.savefig("../overall_analysis/figures/classifier_selection/classifier_comparison"+appendix+".png")
		f.savefig("../overall_analysis/figures/classifier_selection/classifier_comparison"+appendix+".svg")
	
	
	
	#makes figure for classifier selection => barplots
	#methods_notnorm is list of method.__name__
	def fig_selection_of_classifier(self,classifier_oi,method_norm="mean normalized",appendix=""):
		performance_table=pd.read_csv(self.analyses+"analysis_summary"+self.feature+self.approach+".tsv",sep="\t")
		performance_table_norm=pd.read_csv(self.analyses+"analysis_summary"+self.feature[:-15]+self.approach+".tsv",sep="\t").set_index("methods")
		num_rows=int((len(performance_table.columns)-2)/2)
		num_cols=len(performance_table["methods"].unique())+1
		performance_table=performance_table.set_index("methods")
		methods=performance_table.index.unique().tolist()+[method_norm]
		method_dict={"mean":"mean","median":"median","amax":"max","amin":"min","mean normalized":"mean normalized","amin normalized":"min normalized"}
		performance_measure_dict={"precision":"precision","recall":"recall","f_1 score":"F$_1$ score","accuracy":"accuracy","sensitivity":"sensitivity","specificity":"specificity"}
		classifier_dict={"svm_clf(C=0.01,kernel='linear')":"lin. SVM, C=0.01","svm_clf(C=0.1,kernel='linear')":"lin. SVM, C=0.1","svm_clf(C=1,kernel='linear')":"lin. SVM, C=1","svm_clf(C=10,kernel='linear')":"lin. SVM, C=10","svm_clf(C=100,kernel='linear')":"lin. SVM, C=100","svm_clf(kernel='rbf')":"Gauss. SVM","svm_clf(kernel='poly',degree=8)":"poly. SVM, d=8","svm_clf(kernel='sigmoid')":"sigmoid SVM","rf_clf()":"Random Forest"}
		classifier_names=list(map(lambda x: classifier_dict[x],performance_table.loc[methods[0],"classifiers"].tolist()))
		#set red colour for selected classifier of interest, blue for others
		#colors=["tab:blue"]*2+["tab:red"]+["tab:blue"]*(len(classifier_names)-3)
		clf_oi_index=classifier_names.index(classifier_dict[classifier_oi])
		colors=["tab:blue"]*len(classifier_names)
		colors[clf_oi_index]="tab:red"
		alphas=np.linspace(1,0.4,num_cols)
		f,axn=plt.subplots(figsize=(10,10),nrows=num_rows,ncols=num_cols,sharey=True)
		#f,axn=plt.subplots(nrows=2,ncols=2)
		for c in range(num_cols):
			method=methods[c]
			for r in range(num_rows):
				if method!=method_norm:
					xdata=performance_table.loc[method,"classifiers"]
					ydata=performance_table.loc[method,performance_table.columns[r+1]]
					upper_errors=performance_table.loc[method,performance_table.columns[r+1+num_rows]] #/performance_table.loc[method,performance_table.columns[r+1]]
				else:
					xdata=performance_table_norm.loc[method.split(" ")[0],"classifiers"]
					ydata=performance_table_norm.loc[method.split(" ")[0],performance_table_norm.columns[r+1]]
					upper_errors=performance_table_norm.loc[method.split(" ")[0],performance_table_norm.columns[r+1+num_rows]] #/performance_table_norm.loc[method.split(" ")[0],performance_table_norm.columns[r+1]].tolist()
				lower_errors=list(map(lambda x: ydata.iloc[x] if upper_errors.iloc[x]>ydata.iloc[x] else upper_errors.iloc[x],range(len(upper_errors.tolist()))))
				errors=[lower_errors,upper_errors.tolist()]
				axn[r][c].bar(xdata,ydata,yerr=errors,color=colors,alpha=alphas[c],tick_label=classifier_names)
				if r==0:
					axn[r][c].set_title(method_dict[method])
				if c==0:
					axn[r][c].set_ylabel(performance_measure_dict[performance_table.columns[r+1]])
				if r==num_rows-1:
					#axn[r][c].set_xticks(np.arange(len(classifier_names)),classifier_names)
					axn[r][c].set_xticklabels(classifier_names,rotation=65,ha="right")
				else:
					axn[r][c].set_xticklabels(len(classifier_names)*[""])
			#axn[r+1][c].hist()
		#f.tight_layout()
		plt.subplots_adjust(bottom=0.2)
		plt.show()
		f.savefig("../overall_analysis/figures/classifier_selection/selection_of_clf_and_pooling"+appendix+".png")
		f.savefig("../overall_analysis/figures/classifier_selection/selection_of_clf_and_pooling"+appendix+".svg")
	
	
	#methods_notnorm is list of method.__name__
	def fig_selection_of_classifier_1PM(self,classifier_oi,performance_measure,methods_notnorm=["amin","mean","amax"],method_norm="amin normalized",appendix=""):
		performance_table=pd.read_csv(self.analyses+"analysis_summary"+self.feature+self.approach+".tsv",sep="\t").set_index("methods")
		method_dict={"mean":"mean","median":"median","amax":"max","amin":"min","mean normalized":"mean normalized","amin normalized":"min normalized"}
		performance_measure_dict={"precision":"precision","recall":"recall","f_1 score":"F$_1$ score","accuracy":"accuracy","sensitivity":"sensitivity","specificity":"specificity"}
		classifier_dict={"svm_clf(C=0.01,kernel='linear')":"lin. SVM, C=0.01","svm_clf(C=0.1,kernel='linear')":"lin. SVM, C=0.1","svm_clf(C=1,kernel='linear')":"lin. SVM, C=1","svm_clf(C=10,kernel='linear')":"lin. SVM, C=10","svm_clf(C=100,kernel='linear')":"lin. SVM, C=100","svm_clf(kernel='rbf')":"Gauss. SVM","svm_clf(kernel='poly',degree=8)":"poly. SVM, d=8","svm_clf(kernel='sigmoid')":"sigmoid SVM","rf_clf()":"Random Forest"}
		classifier_names=list(map(lambda x: classifier_dict[x],performance_table.loc[methods_notnorm[0],"classifiers"].tolist()))
		methods=methods_notnorm+[method_norm]
		clf_oi_index=classifier_names.index(classifier_dict[classifier_oi])
		colors=["tab:blue"]*len(classifier_names)
		colors[clf_oi_index]="tab:red"
		alphas=np.linspace(1,0.6,4)
		num_rows=2
		num_cols=2
		f,axn=plt.subplots(figsize=(10,8),nrows=num_rows,ncols=num_cols,sharey=True)
		m=0
		pm_index=performance_table.columns.tolist().index(performance_measure)
		for c in range(num_cols):
			for r in range(num_rows):
				method=methods[m]
				if method!=method_norm:
					xdata=performance_table.loc[method,"classifiers"]
					ydata=performance_table.loc[method,performance_measure]
					upper_errors=performance_table.loc[method,performance_table.columns[pm_index+int((len(performance_table.columns)-1)/2)]] #/performance_table.loc[method,performance_table.columns[r+1]]
				else:
					xdata=performance_table_norm.loc[method.split(" ")[0],"classifiers"]
					ydata=performance_table_norm.loc[method.split(" ")[0],performance_measure]
					upper_errors=performance_table_norm.loc[method.split(" ")[0],performance_table.columns[pm_index+int((len(performance_table.columns)-1)/2)]] #/performance_table_norm.loc[method.split(" ")[0],performance_table_norm.columns[r+1]].tolist()
				lower_errors=list(map(lambda x: ydata.iloc[x] if upper_errors.iloc[x]>ydata.iloc[x] else upper_errors.iloc[x],range(len(upper_errors.tolist()))))
				errors=[lower_errors,upper_errors.tolist()]
				axn[r][c].bar(xdata,ydata,yerr=errors,color=colors,alpha=alphas[m],tick_label=classifier_names)
				axn[r][c].set_title(method_dict[method])
				if c==0:
					axn[r][c].set_ylabel(performance_measure_dict[performance_measure])
				if r==1:
					#axn[r][c].set_xticks(np.arange(len(classifier_names)),classifier_names)
					axn[r][c].set_xticklabels(classifier_names,rotation=65,ha="right")
				else:
					axn[r][c].set_xticklabels(len(classifier_names)*[""])
				m+=1
			#axn[r+1][c].hist()
		#f.tight_layout()
		plt.subplots_adjust(bottom=0.2)
		plt.show()
		f.savefig("../overall_analysis/figures/classifier_selection/selection_of_clf_and_pooling_"+performance_measure+appendix+".png")
		f.savefig("../overall_analysis/figures/classifier_selection/selection_of_clf_and_pooling_"+performance_measure+appendix+".svg")
		return
	
	
	
	#https://fcpython.com/visualisation/radar-charts-matplotlib
	#methods_notnorm is list of method.__name__
	def radar_chart_selection_classifier(self,methods_notnorm,method_norm="amin normalized",colour="tab:blue",classifier=None,appendix="",title=False):
		performance_table=pd.read_csv(self.analyses+"analysis_summary"+self.feature+self.approach+".tsv",sep="\t").set_index("methods")
		method_dict={"mean":"mean","median":"median","amax":"max","amin":"min","mean normalized":"mean normalized","amin normalized":"min normalized"}
		performance_measure_dict={"precision":"precision","recall":"recall","f_1 score":"F$_1$ score","accuracy":"accuracy","sensitivity":"sensitivity","specificity":"specificity"}
		performance_measures=list(performance_measure_dict)
		performance_measures_names=list(map(lambda x: performance_measure_dict[x],performance_measures))
		classifier_dict={"svm_clf(C=0.01,kernel='linear')":"linear SVM, C=0.01","svm_clf(C=0.1,kernel='linear')":"linear SVM, C=0.1","svm_clf(C=1,kernel='linear')":"linear SVM, C=1","svm_clf(C=10,kernel='linear')":"linear SVM, C=10","svm_clf(C=100,kernel='linear')":"linear SVM, C=100","svm_clf(kernel='rbf')":"Gaussian SVM","svm_clf(kernel='poly',degree=8)":"polynomial SVM, d=8","svm_clf(kernel='sigmoid')":"sigmoid SVM","rf_clf()":"Random Forest"}
		classifier_names=list(map(lambda x: classifier_dict[x],performance_table.loc[methods_notnorm[0],"classifiers"].tolist()))
		#one method and one classifier only, classifier needs to be specified
		if len(methods_notnorm)==1:
			method=methods_notnorm[0]
			performance_selection=performance_table.loc[method,performance_measures+[performance_measures[0],"classifiers"]].set_index("classifiers")
			f,axn=plt.subplots(figsize=(8,8),nrows=1,ncols=1)
			angles=[n/len(performance_measures)*2*np.pi for n in range(len(performance_measures))]
			angles+=angles[:1]
			ax=plt.subplot(1,1,1,polar=True)
			plt.xticks(angles[:-1],performance_measures_names)
			ax.set_xticklabels(performance_measures_names)
			ax.plot(angles,performance_selection.loc[classifier],label=classifier_dict[classifier],color=colour)
			ax.fill(angles,performance_selection.loc[classifier],color=colour,alpha=0.1)
			if title:
				ax.set_title(classifier_dict[classifier])
			plt.show()
			f.savefig("../overall_analysis/figures/classifier_selection/radar_charts"+appendix+".png")
			f.savefig("../overall_analysis/figures/classifier_selection/radar_charts"+appendix+".svg")
		else:
			# several methods
			performance_table_norm=pd.read_csv(self.analyses+"analysis_summary"+self.feature[:-15]+self.approach+".tsv",sep="\t").set_index("methods")
			num_rows=int(np.ceil((len(methods_notnorm)+1)/2))
			num_cols=2
			methods=methods_notnorm+[method_norm]
			f,axn=plt.subplots(figsize=(10,10),nrows=num_rows,ncols=num_cols)
			for m,method in enumerate(methods):
				#trim and combine performance_table for analysis
				if method==method_norm:
					performance_selection=performance_table_norm.loc[method.split(" ")[0],performance_measures+[performance_measures[0],"classifiers"]].set_index("classifiers")
				else:
					performance_selection=performance_table.loc[method,performance_measures+[performance_measures[0],"classifiers"]].set_index("classifiers")
				angles=[n/len(performance_measures)*2*np.pi for n in range(len(performance_measures))]
				angles+=angles[:1]
				ax=plt.subplot(num_rows,num_cols,m+1,polar=True)
				plt.xticks(angles[:-1],performance_measures_names)
				ax.set_xticklabels(performance_measures_names)
				for classifier in performance_selection.index:
					ax.plot(angles,performance_selection.loc[classifier],label=classifier_dict[classifier])
					ax.fill(angles,performance_selection.loc[classifier],alpha=0.1)
					ax.set_title(method_dict[method])
			f.tight_layout()
			#plt.subplots_adjust(wspace=0.3)
			plt.show()
			f.savefig("../overall_analysis/figures/classifier_selection/radar_charts"+appendix+".png")
			f.savefig("../overall_analysis/figures/classifier_selection/radar_charts"+appendix+".svg")
		return
	
	
	
	#histogram for frequency of prediction probabilities for different 
	def hist_predprobs(self,method,classifier):
		predprobs=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index(["OG","metabolite"])
		f=plt.figure()
		plt.hist(predprobs[classifier])
		plt.xlabel("P(Interaction)")
		plt.ylabel("abs. frequency")
		f.savefig(self.analyses+"histogram_prob(interaction).png")
		f.savefig(self.analyses+"histogram_prob(interaction).svg")
		
	
	
	#################################
	#COMPARISON TO EXPERIMENTAL DATA
	#################################
	
	#Comparison from Yeast Ser-Leu targets to AP and TPP experimental data
	#predictions either from eucaryotes or yeast experiment
	def comparison_to_ap_tpp_pcc(self,exp_predictions,appendix=""):
		#function arguments
		tableau_red=(214/255,39/255,40/255,0.5) 
		tableau_blue=(31/255,119/255,180/255,0.5)
		tableau_green=(44/255,160/255,44/255,0.6)
		tableau_orange=(255/255,127/255,14/255,0.4)
		metabolite=13919048 #Ser-Leu ###
		inputfile="../overall_analysis/Ser-Leu targets - putative.xlsx" ###
		sheet="Sheet1" ###
		sce_idtrans_df=pd.read_csv("../databases/S.cerevisiae_ID_collection.tsv",sep="\t").set_index("GENENAME")
		all_union_probs=pd.read_csv("../overall_analysis/union_prediction_probabilities_no_me2.tsv",sep="\t")
		#extract experimental targets and protein sets
		target_excel=pd.ExcelFile(inputfile)
		target_df=target_excel.parse(sheet)
		tpp_all=set(target_df["All proteins TPP"].dropna().tolist())
		tpp_all_uniprot=sce_idtrans_df.loc[tpp_all,"ID"].dropna().tolist()
		tpp_targets=set(target_df["Targets TPP"].dropna().tolist())
		tpp_targets_uniprot=sce_idtrans_df.loc[tpp_targets,"ID"].dropna().tolist()
		ap_targets=set(target_df["Targets AP"].dropna().tolist())
		ap_targets_uniprot=sce_idtrans_df.loc[ap_targets,"ID"].dropna().tolist()
		ap_all=set(target_df["All proteins AP"].dropna().tolist())
		ap_all_uniprot=sce_idtrans_df.loc[ap_all,"ID"].dropna().tolist()
		pcc_targets=set(target_df["Targets PCC"].dropna().tolist())
		pcc_targets_uniprot=sce_idtrans_df.loc[pcc_targets,"ID"].dropna().tolist()
		#extract ML results
		exp_probs=all_union_probs[~pd.isna(all_union_probs["prediction_probability_Saccharomyces"])]
		exp_metabolite_predictions=exp_predictions[exp_predictions.index.get_level_values(1)==metabolite]
		if self.proteinwise==False:
			exp_ogs=exp_metabolite_predictions.index.get_level_values(0).unique()
			exp_prots=set(";".join(all_union_probs.set_index("OG").loc[exp_ogs,"Saccharomyces_UniProt_IDs"].dropna().tolist()).split(";"))
			exp_int_ogs=exp_metabolite_predictions[exp_metabolite_predictions["prediction"]==True].index.get_level_values(0).unique()
			exp_int_prots=set(";".join(all_union_probs.set_index("OG").loc[exp_int_ogs,"Saccharomyces_UniProt_IDs"].dropna().tolist()).split(";"))
		else:
			exp_prots=set(exp_metabolite_predictions.index.get_level_values(0).unique().tolist())
			exp_int_prots=set(exp_metabolite_predictions[exp_metabolite_predictions["prediction"]==True].index.get_level_values(0).unique().tolist())
		ml_targets_uniprot=list(exp_int_prots)
		#trim to intersection
		tpp_targets_uniprot_in_ml=set(set(tpp_targets_uniprot) & exp_prots)
		ap_targets_uniprot_in_ml=set(set(ap_targets_uniprot) & exp_prots)
		pcc_targets_uniprot_in_ml=set(set(pcc_targets_uniprot) & exp_prots)
		tpp_targets_uniprot_in_ap=set(set(tpp_targets_uniprot) & set(ap_all_uniprot))
		ap_targets_uniprot_in_tpp=set(set(ap_targets_uniprot) & set(tpp_all_uniprot))
		tpp_targets_uniprot_in_ml_and_ap=set(set(tpp_targets_uniprot) & exp_prots & set(ap_all_uniprot))
		ap_targets_uniprot_in_ml_and_tpp=set(set(ap_targets_uniprot) & exp_prots & set(tpp_all_uniprot))
		ap_tpp_targets_uniprot_in_ml=set(set(tpp_targets_uniprot) & set(ap_targets_uniprot) & exp_prots)
		ml_targets_uniprot_in_tpp=set(exp_int_prots & set(tpp_all_uniprot))
		ml_targets_uniprot_in_ap=set(exp_int_prots & set(ap_all_uniprot))
		ml_targets_uniprot_in_ap_or_tpp=set(exp_int_prots & set(set(ap_all_uniprot) | set(tpp_all_uniprot)))
		ml_targets_uniprot_in_ap_and_tpp=set(exp_int_prots & set(set(ap_all_uniprot) & set(tpp_all_uniprot)))
		#prediction dataframes for cohen's kappa calculation
		exp_metabolite_predictions=exp_metabolite_predictions.reset_index().set_index("protein")
		tpp_df=pd.DataFrame(index=tpp_all_uniprot,columns=["prediction"])
		tpp_df["prediction"]=False
		tpp_df.loc[tpp_targets_uniprot,"prediction"]=True
		ap_df=pd.DataFrame(index=ap_all_uniprot,columns=["prediction"])
		ap_df["prediction"]=False
		try:
			ap_df.loc[ap_targets_uniprot,"prediction"]=True
		except KeyError:
			for p in ap_targets_uniprot:
				ap_df.loc[p,"prediction"]=True
		tpp_ml_prots_intersect=set(tpp_df.index) & set(exp_metabolite_predictions.index)
		ap_ml_prots_intersect=set(ap_df.index) & set(exp_metabolite_predictions.index)
		ap_tpp_prots_intersect=set(ap_df.index) & set(tpp_df.index)
		#VENN DIAGRAMS
		#ap, tpp, ml
		#interactions of proteins present in all three ways
		data=venn.get_labels([tpp_targets_uniprot_in_ml_and_ap,ap_targets_uniprot_in_ml_and_tpp,ml_targets_uniprot_in_ap_and_tpp])
		f,axn=venn.venn3(data,names=["TPP","AP","ML"],colors=[tableau_blue,tableau_red,tableau_orange])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ap_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ap_ml"+appendix+".svg")
		#ap, tpp, pcc, ml, for proteins in ML
		data=venn.get_labels([tpp_targets_uniprot_in_ml,ap_targets_uniprot_in_ml,pcc_targets_uniprot_in_ml,ml_targets_uniprot_in_ap_and_tpp])
		f,axn=venn.venn4(data,names=["TPP","AP","PCC","ML"],colors=[tableau_blue,tableau_red,tableau_green,tableau_orange])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ap_pcc_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ap_pcc_ml"+appendix+".svg")
		
		## venn diagrams showing how many of the with ML found targets are verified by experimental approach
		#intersect tpp - ml
		data=venn.get_labels([tpp_targets_uniprot_in_ml,ml_targets_uniprot_in_tpp])
		tpp_ml_df=pd.DataFrame(index=tpp_ml_prots_intersect,columns=["TPP","ML"])
		tpp_ml_df.loc[tpp_ml_prots_intersect,"TPP"]=tpp_df.loc[tpp_ml_prots_intersect,"prediction"]
		tpp_ml_df.loc[tpp_ml_prots_intersect,"ML"]=exp_metabolite_predictions.loc[tpp_ml_prots_intersect,"prediction"]
		print("kappa TPP-ML: "+str(metrics.cohen_kappa_score(tpp_ml_df["TPP"],tpp_ml_df["ML"])))
		f,axn=venn.venn2(data,names=["TPP","ML"],colors=[tableau_blue,tableau_orange])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ml"+appendix+".svg")
		#intersect ap - ml
		data=venn.get_labels([ap_targets_uniprot_in_ml,ml_targets_uniprot_in_ap])
		ap_ml_df=pd.DataFrame(index=ap_ml_prots_intersect,columns=["AP","ML"])
		ap_ml_df.loc[ap_ml_prots_intersect,"AP"]=ap_df.loc[ap_ml_prots_intersect,"prediction"]
		ap_ml_df.loc[ap_ml_prots_intersect,"ML"]=exp_metabolite_predictions.loc[ap_ml_prots_intersect,"prediction"]
		print("kappa AP-ML: "+str(metrics.cohen_kappa_score(ap_ml_df["AP"],ap_ml_df["ML"])))
		f,axn=venn.venn2(data,names=["AP","ML"],colors=[tableau_red,tableau_orange])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_targets_ap_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_targets_ap_ml"+appendix+".svg")
		#intersect pcc-ml
		data=venn.get_labels([pcc_targets_uniprot_in_ml,ml_targets_uniprot])
		f,axn=venn.venn2(data,names=["PCC","ML"],colors=[tableau_green,tableau_orange])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_targets_pcc_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_targets_pcc_ml"+appendix+".svg")
		#intersect tpp-ap
		data=venn.get_labels([tpp_targets_uniprot_in_ap,ap_targets_uniprot_in_tpp])
		ap_tpp_df=pd.DataFrame(index=ap_tpp_prots_intersect,columns=["AP","TPP"])
		ap_tpp_df.loc[ap_tpp_prots_intersect,"AP"]=ap_df.loc[ap_tpp_prots_intersect,"prediction"]
		ap_tpp_df.loc[ap_tpp_prots_intersect,"TPP"]=tpp_df.loc[ap_tpp_prots_intersect,"prediction"]
		print("kappa AP-TPP: "+str(metrics.cohen_kappa_score(ap_tpp_df["AP"],ap_tpp_df["TPP"])))
		f,axn=venn.venn2(data,names=["TPP","AP"],colors=[tableau_blue,tableau_red])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ap.png")
		f.savefig("../overall_analysis/venn_diagram_targets_tpp_ap.svg")
		return
		
		
	
	
	#comparison of yeast Ser-Leu targets to PCC results
	def comparison_to_pcc(self,exp_predictions,appendix=""):
		metabolite=13919048 #Ser-Leu ###
		inputfile="../overall_analysis/Ser-Leu targets - putative.xlsx" ###
		sheet="Sheet1" ###
		#extract ML results
		sce_idtrans_df=pd.read_csv("../databases/S.cerevisiae_ID_collection.tsv",sep="\t").set_index("GENENAME")
		all_union_probs=pd.read_csv("../overall_analysis/union_prediction_probabilities_no_me2.tsv",sep="\t")
		exp_metabolite_predictions=exp_predictions[exp_predictions.index.get_level_values(1)==metabolite]
		if self.proteinwise==False:
			exp_probs=all_union_probs[~pd.isna(all_union_probs["prediction_probability_Saccharomyces"])]
			exp_ogs=exp_metabolite_predictions.index.get_level_values(0).unique()
			exp_prots=set(";".join(all_union_probs.set_index("OG").loc[exp_ogs,"Saccharomyces_UniProt_IDs"].dropna().tolist()).split(";"))
			exp_int_ogs=exp_metabolite_predictions[exp_metabolite_predictions["prediction"]==True].index.get_level_values(0).unique()
			exp_int_prots=set(";".join(all_union_probs.set_index("OG").loc[exp_int_ogs,"Saccharomyces_UniProt_IDs"].dropna().tolist()).split(";"))
		else:
			exp_prots=set(exp_metabolite_predictions.index.get_level_values(0).unique().tolist())
			exp_int_prots=set(exp_metabolite_predictions[exp_metabolite_predictions["prediction"]==True].index.get_level_values(0).unique().tolist())
		ml_targets_uniprot=list(exp_int_prots)
		#extract AP and TPP results
		target_excel=pd.ExcelFile(inputfile)
		target_df=target_excel.parse(sheet)
		tpp_all=set(target_df["All proteins TPP"].dropna().tolist())
		tpp_all_uniprot=sce_idtrans_df.loc[tpp_all,"ID"].dropna().tolist()
		tpp_targets=set(target_df["Targets TPP"].dropna().tolist())
		tpp_targets_uniprot=sce_idtrans_df.loc[tpp_targets,"ID"].dropna().tolist()
		ap_targets=set(target_df["Targets AP"].dropna().tolist())
		ap_targets_uniprot=sce_idtrans_df.loc[ap_targets,"ID"].dropna().tolist()
		ap_all=set(target_df["All proteins AP"].dropna().tolist())
		ap_all_uniprot=sce_idtrans_df.loc[ap_all,"ID"].dropna().tolist()
		### comparison to PCC, no information about all proteins, so interactions not trimmed to intersecting proteins
		pcc_targets_inputfile="../overall_analysis/Ser-Leu_yeast_PCC.xlsx"
		cc_inputfile=self.experimental+"S.cerevisiae_cell_cultures.xlsx"
		gp1_inputfile=self.experimental+"S.cerevisiae_growthphases_condition1.xlsx"
		gp2_inputfile=self.experimental+"S.cerevisiae_growthphases_condition2.xlsx"
		gp3_inputfile=self.experimental+"S.cerevisiae_growthphases_condition3.xlsx"
		pcc_targets_excel=pd.ExcelFile(pcc_targets_inputfile)
		cc_excel=pd.ExcelFile(cc_inputfile)
		gp1_excel=pd.ExcelFile(gp1_inputfile)
		gp2_excel=pd.ExcelFile(gp1_inputfile)
		gp3_excel=pd.ExcelFile(gp1_inputfile)
		cc_protsheet=cc_excel.parse("Protein_profiles_rep_raw")
		gp1_protsheet=gp1_excel.parse("Protein_profiles_rep_raw")
		gp2_protsheet=gp2_excel.parse("Protein_profiles_rep_raw")
		gp3_protsheet=gp3_excel.parse("Protein_profiles_rep_raw")
		cc_all_proteins=set(cc_protsheet["UniProt_IDs"][1:].dropna().tolist())
		gp1_all_proteins=set(gp1_protsheet["UniProt_IDs"][1:].dropna().tolist())
		gp2_all_proteins=set(gp2_protsheet["UniProt_IDs"][1:].dropna().tolist())
		gp3_all_proteins=set(gp3_protsheet["UniProt_IDs"][1:].dropna().tolist())
		intersecting_prots_pcc=set(cc_all_proteins & gp1_all_proteins & gp2_all_proteins & gp3_all_proteins) #intersect_pcc
		intersecting_prots_ml=set(cc_all_proteins & gp1_all_proteins & gp2_all_proteins & gp3_all_proteins & exp_prots) #intersect_pcc_ml
		all_prots_df=pcc_targets_excel.parse("Ser-Leu")
		cc_target_df=pcc_targets_excel.parse("Log Phase - 1st dataset")
		gp1_target_df=pcc_targets_excel.parse("Log Phase - 2nd dataset")
		gp2_target_df=pcc_targets_excel.parse("Diauxic - 3rd dataset")
		gp3_target_df=pcc_targets_excel.parse("Stationary - 4th dataset")
		all_pcc_prots=set(all_prots_df["Protein_ID"].tolist())
		gp1_targets=set(gp1_target_df.iloc[:,0].tolist())
		gp2_targets=set(gp2_target_df.iloc[:,0].tolist())
		gp3_targets=set(gp3_target_df.iloc[:,0].tolist())
		cc_targets=set(cc_target_df.iloc[:,0].tolist())
		gp1_targets_in_intersect_pcc=set(gp1_targets & intersecting_prots_pcc)
		gp2_targets_in_intersect_pcc=set(gp2_targets & intersecting_prots_pcc)
		gp3_targets_in_intersect_pcc=set(gp3_targets & intersecting_prots_pcc)
		cc_targets_in_intersect_pcc=set(cc_targets & intersecting_prots_pcc)
		gp1_targets_in_intersect=set(gp1_targets & intersecting_prots_ml)
		gp2_targets_in_intersect=set(gp2_targets & intersecting_prots_ml)
		gp3_targets_in_intersect=set(gp3_targets & intersecting_prots_ml)
		cc_targets_in_intersect=set(cc_targets & intersecting_prots_ml)
		ml_targets_in_intersect=set(set(ml_targets_uniprot) & intersecting_prots_ml)
		#Venn diagrams
		#ml with pcc intersect
		data=venn.get_labels([ml_targets_in_intersect,set(gp1_targets_in_intersect & gp2_targets_in_intersect & gp3_targets_in_intersect & cc_targets_in_intersect)])
		f,axn=venn.venn2(data,names=["ml","pcc"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_intersect_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_intersect_ml"+appendix+".svg")
		#ml with pcc solo
		data=venn.get_labels([ml_targets_in_intersect,gp1_targets_in_intersect,gp2_targets_in_intersect,gp3_targets_in_intersect,cc_targets_in_intersect])
		f,axn=venn.venn5(data,names=["ml","gp1","gp2","gp3","cc"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_ml"+appendix+".svg")
		#pcc solo
		data=venn.get_labels([gp1_targets_in_intersect_pcc,gp2_targets_in_intersect_pcc,gp3_targets_in_intersect_pcc,cc_targets_in_intersect_pcc])
		f,axn=venn.venn4(data,names=["gp1","gp2","gp3","cc"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc.png")
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc.svg")
		#comparison pcc to ap and tpp
		pcc_intersect=set(cc_all_proteins & gp1_all_proteins & gp2_all_proteins & gp3_all_proteins & set(ap_all_uniprot) & set(tpp_all_uniprot))
		pcc_intersect_ap=set(cc_all_proteins & gp1_all_proteins & gp2_all_proteins & gp3_all_proteins & set(ap_all_uniprot))
		pcc_intersect_tpp=set(cc_all_proteins & gp1_all_proteins & gp2_all_proteins & gp3_all_proteins & set(tpp_all_uniprot))
		cc_targets_pcc_intersect=set(cc_targets & pcc_intersect)
		gp1_targets_pcc_intersect=set(gp1_targets & pcc_intersect)
		gp2_targets_pcc_intersect=set(gp2_targets & pcc_intersect)
		gp3_targets_pcc_intersect=set(gp3_targets & pcc_intersect)
		ap_targets_pcc_intersect=set(set(ap_targets_uniprot) & pcc_intersect)
		tpp_targets_pcc_intersect=set(set(tpp_targets_uniprot) & pcc_intersect)
		cc_targets_pcc_intersect_ap=set(cc_targets & pcc_intersect_ap)
		gp1_targets_pcc_intersect_ap=set(gp1_targets & pcc_intersect_ap)
		gp2_targets_pcc_intersect_ap=set(gp2_targets & pcc_intersect_ap)
		gp3_targets_pcc_intersect_ap=set(gp3_targets & pcc_intersect_ap)
		ap_targets_pcc_intersect_ap=set(set(ap_targets_uniprot) & pcc_intersect_ap)
		cc_targets_pcc_intersect_tpp=set(cc_targets & pcc_intersect_tpp)
		gp1_targets_pcc_intersect_tpp=set(gp1_targets & pcc_intersect_tpp)
		gp2_targets_pcc_intersect_tpp=set(gp2_targets & pcc_intersect_tpp)
		gp3_targets_pcc_intersect_tpp=set(gp3_targets & pcc_intersect_tpp)
		tpp_targets_pcc_intersect_tpp=set(set(tpp_targets_uniprot) & pcc_intersect_tpp)
		#pcc to ap
		data=venn.get_labels([cc_targets_pcc_intersect_ap,gp1_targets_pcc_intersect_ap,gp2_targets_pcc_intersect_ap,gp3_targets_pcc_intersect_ap,ap_targets_pcc_intersect_ap])
		f,axn=venn.venn5(data,names=["cc","gp1","gp2","gp3","ap"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_ap.png")
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_ap.svg")
		#pcc to tpp
		data=venn.get_labels([cc_targets_pcc_intersect_tpp,gp1_targets_pcc_intersect_tpp,gp2_targets_pcc_intersect_tpp,gp3_targets_pcc_intersect_tpp,tpp_targets_pcc_intersect_tpp])
		f,axn=venn.venn5(data,names=["cc","gp1","gp2","gp3","tpp"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_tpp.png")
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_tpp.svg")
		#pcc to ap and tpp
		data=venn.get_labels([cc_targets_pcc_intersect,gp1_targets_pcc_intersect,gp2_targets_pcc_intersect,gp3_targets_pcc_intersect,ap_targets_pcc_intersect,tpp_targets_pcc_intersect])
		f,axn=venn.venn6(data,names=["cc","gp1","gp2","gp3","ap","tpp"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_ap_tpp.png")
		f.savefig("../overall_analysis/venn_diagram_yeast_serleu_targets_pcc_ap_tpp.svg")
		return
	
	
	
	
	#comparison tyr-asp to arabidopsis predictions
	def comparison_to_ap_tpp_pcc_tyrasp_targets(self,predictions,appendix=""):
		tableau_red=(214/255,39/255,40/255,0.5) 
		tableau_blue=(31/255,119/255,180/255,0.5)
		tableau_green=(44/255,160/255,44/255,0.6)
		tableau_orange=(255/255,127/255,14/255,0.4)
		metabolite=19816752 #Tyr-Asp
		# IDs given in Araport (or TAIR) identifiers
		tyr_asp_targets_filename="../overall_analysis/TyrAsp_ara_targets.xlsx"
		tyr_asp_tpp_prots_filename="../overall_analysis/TPP_all_prots.xlsx"
		trans_db=pd.read_csv("../databases/A.thaliana_ID_collection.tsv",sep="\t").set_index("ID")
		ap_sheet="AP"
		tpp_sheet="TPP"
		pcc_sheet="cc"
		#load tyr-asp targets from comparing approaches
		target_excel=pd.ExcelFile(tyr_asp_targets_filename)
		ap_target_df=target_excel.parse(ap_sheet)
		tpp_target_df=target_excel.parse(tpp_sheet)
		pcc_target_df=target_excel.parse(pcc_sheet)
		ap_targets=set(ap_target_df["UniProt_ID"].tolist())
		tpp_targets=set(tpp_target_df["UniProt_ID"].tolist())
		seedlings_targets=set(pcc_target_df["Seedlings_UniProt"].tolist())
		cellcultures_targets=set(pcc_target_df["Cell Cultures_UniProt"].tolist())
		rosettes_targets=set(pcc_target_df["Rosettes_UniProt"].tolist())
		#pcc_targets as intersection
		pcc_targets=set(seedlings_targets & cellcultures_targets & rosettes_targets)
		#translate ML interactions to TAIR
		ml_preds_metabolite=predictions.reset_index().set_index("metabolite").loc[metabolite]
		ml_ints_metabolite=ml_preds_metabolite[ml_preds_metabolite["prediction"]==True]
		if self.proteinwise==True:
			ml_targets=set(ml_ints_metabolite[self.protcol].tolist())
		else:
			print("transform predictions from og-wise to proteinwise")
		#ml_targets=set(trans_db.loc[ml_targets_uniprot,"ENSEMBLGENOME_PRO_ID"].tolist())
		#load all proteins for TPP, AP is for whole proteome
		tpp_all_prots_excel=pd.ExcelFile(tyr_asp_tpp_prots_filename)
		tpp_all_df=tpp_all_prots_excel.parse("Sheet1")
		tpp_all=set(tpp_all_df["UniProt_ID"].tolist())
		#get intersection of all proteins
		ml_prots=set(ml_preds_metabolite["protein"].tolist())
		#ml_prots=set(trans_db.loc[ml_prots_uniprot,"ENSEMBLGENOME_PRO_ID"].tolist())
		all_prots=set(ml_prots & tpp_all)
		print(str(len(all_prots))+" proteins in intersection of experiments")
		#trim proteins to intersecting prots
		ap_targets_intersect=set(ap_targets & all_prots)
		ap_targets_in_tpp=set(ap_targets & tpp_all)
		tpp_targets_intersect=set(tpp_targets & all_prots)
		ml_targets_intersect=set(ml_targets & all_prots)
		pcc_targets_intersect=set(pcc_targets & all_prots)
		#prediction dataframes for cohen's kappa calculation
		ml_preds_metabolite=ml_preds_metabolite.reset_index().set_index("protein")
		tpp_df=pd.DataFrame(index=tpp_all,columns=["prediction"])
		tpp_df["prediction"]=False
		try:
			tpp_df.loc[tpp_targets,"prediction"]=True
		except KeyError:
			for p in tpp_targets:
				tpp_df.loc[p,"prediction"]=True
		ap_df=pd.DataFrame(index=set(ml_prots | tpp_all),columns=["prediction"])
		ap_df["prediction"]=False
		try:
			ap_df.loc[ap_targets,"prediction"]=True
		except KeyError:
			for p in ap_targets:
				ap_df.loc[p,"prediction"]=True
		tpp_ml_prots_intersect=set(tpp_df.index) & set(ml_preds_metabolite.index)
		ap_ml_prots_intersect=set(ap_df.index) & set(ml_preds_metabolite.index)
		ap_tpp_prots_intersect=set(ap_df.index) & set(tpp_df.index)
		#venn diagram ap, tpp, and ml
		data=venn.get_labels([tpp_targets_intersect,ap_targets_intersect,ml_targets_intersect])
		f,axn=venn.venn3(data,names=["TPP","AP","ML"],colors=[tableau_blue,tableau_red,tableau_green])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_tpp_ap_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_tpp_ap_ml"+appendix+".svg")
		#venn diagram ap, tpp, pcc and ml
		data=venn.get_labels([tpp_targets_intersect,ap_targets_intersect,ml_targets_intersect,pcc_targets_intersect])
		f,axn=venn.venn4(data,names=["TPP","AP","ML","PCC"],colors=[tableau_blue,tableau_red,tableau_green,tableau_orange])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_tpp_ap_ml_pcc"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_tpp_ap_ml_pcc"+appendix+".svg")
		#venn diagram ap and ml
		data=venn.get_labels([ap_targets,ml_targets])
		ap_ml_df=pd.DataFrame(index=ap_ml_prots_intersect,columns=["AP","ML"])
		ap_ml_df.loc[ap_ml_prots_intersect,"AP"]=ap_df.loc[ap_ml_prots_intersect,"prediction"]
		ap_ml_df.loc[ap_ml_prots_intersect,"ML"]=ml_preds_metabolite.loc[ap_ml_prots_intersect,"prediction"]
		print("kappa AP-ML: "+str(metrics.cohen_kappa_score(ap_ml_df["AP"],ap_ml_df["ML"])))
		f,axn=venn.venn2(data,names=["AP","ML"],colors=[tableau_red,tableau_green])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_ap_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_ap_ml"+appendix+".svg")
		#venn diagram tpp and ml
		data=venn.get_labels([tpp_targets_intersect,ml_targets_intersect])
		tpp_ml_df=pd.DataFrame(index=tpp_ml_prots_intersect,columns=["TPP","ML"])
		tpp_ml_df.loc[tpp_ml_prots_intersect,"TPP"]=tpp_df.loc[tpp_ml_prots_intersect,"prediction"]
		tpp_ml_df.loc[tpp_ml_prots_intersect,"ML"]=ml_preds_metabolite.loc[tpp_ml_prots_intersect,"prediction"]
		print("kappa TPP-ML: "+str(metrics.cohen_kappa_score(tpp_ml_df["TPP"],tpp_ml_df["ML"])))
		f,axn=venn.venn2(data,names=["TPP","ML"],colors=[tableau_blue,tableau_green])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_tpp_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_tpp_ml"+appendix+".svg")
		#venn diagram pcc and ml
		data=venn.get_labels([pcc_targets,ml_targets])
		f,axn=venn.venn2(data,names=["PCC","ML"],colors=[tableau_orange,tableau_green])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_pcc_ml"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_pcc_ml"+appendix+".svg")
		#venn diagram tpp and ap
		data=venn.get_labels([ap_targets_in_tpp,tpp_targets])
		ap_tpp_df=pd.DataFrame(index=ap_tpp_prots_intersect,columns=["AP","TPP"])
		ap_tpp_df.loc[ap_tpp_prots_intersect,"AP"]=ap_df.loc[ap_tpp_prots_intersect,"prediction"]
		ap_tpp_df.loc[ap_tpp_prots_intersect,"TPP"]=tpp_df.loc[ap_tpp_prots_intersect,"prediction"]
		print("kappa AP-TPP: "+str(metrics.cohen_kappa_score(ap_tpp_df["AP"],ap_tpp_df["TPP"])))
		f,axn=venn.venn2(data,names=["AP","TPP"],colors=[tableau_red,tableau_blue])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_ap_tpp"+appendix+".png")
		f.savefig("../overall_analysis/venn_diagram_tyrasp_targets_ap_tpp"+appendix+".svg")
		return
	
	
	#compares SLIMP to AP and TPP on Ser-Leu targets in yeast (venn diagrams)
	def comparison_serleu_AP_TPP_SLIMP(self,predictions_yeast):
		metabolite=13919048 #Ser-Leu ###
		#load sheet/proteins
		inputfile=pd.ExcelFile("../ForBoris.xlsx")
		serleu=inputfile.parse("SerLeu")
		slimp_proteins_uniprot=predictions_yeast.index.get_level_values(0).unique()
		##slimp targets threshold 80% probability
		slimp80_targets_uniprot=set(serleu["SLIMP80 (Uniprot)"].tolist())
		##slimp targets no threshold
		serleu_predictions=predictions_yeast.reset_index().set_index("metabolite").loc[metabolite]
		serleu_targets=serleu_predictions[serleu_predictions["prediction"]==True]
		if self.proteinwise==True:
			slimp_targets_uniprot=set(serleu_targets[self.protcol].tolist())
		else:
			print("transform predictions from og-wise to proteinwise")
		ap_targets=set(serleu["AP"].tolist())
		tpp_targets=set(serleu["TPP (differential)"].tolist())
		tpp_proteins=set(serleu["TPP (all measured proteins)"].tolist())
		tpp_proteins=tpp_targets | tpp_proteins # targets not in all measured proteins
		#translate IDs
		trans_df=da_ara.translate_proteins(slimp_proteins_uniprot.tolist(),in_id_type="ID",out_id_types=["GENENAME"],save=False)
		slimp_proteins=set(trans_df["GENENAME"].tolist())
		slimp_targets=set(trans_df.loc[slimp_targets_uniprot,"GENENAME"].tolist())
		slimp80_targets=set(trans_df.loc[slimp80_targets_uniprot,"GENENAME"].tolist())
		#remove np.nans
		if np.nan in ap_targets:
			ap_targets.remove(np.nan)
		if np.nan in tpp_targets:
			tpp_targets.remove(np.nan)
		if np.nan in tpp_proteins:
			tpp_proteins.remove(np.nan)
		if np.nan in slimp_targets:
			slimp_targets.remove(np.nan)
		if np.nan in slimp80_targets:
			slimp80_targets.remove(np.nan)
		if np.nan in slimp_proteins:
			slimp_proteins.remove(np.nan)
		#get proteins in intersection among experiments
		tpp_targets_intersect=tpp_targets & slimp_proteins
		ap_targets_intersect=ap_targets & slimp_proteins & tpp_proteins #
		slimp_targets_intersect=slimp_targets & tpp_proteins
		slimp80_targets_intersect=slimp80_targets & tpp_proteins
		data=venn.get_labels([ap_targets_intersect,tpp_targets_intersect,slimp_targets_intersect])
		data80=venn.get_labels([ap_targets_intersect,tpp_targets_intersect,slimp80_targets_intersect])
		#data=venn.get_labels([ap_targets_intersect,tpp_targets_intersect,slimp_targets])
		#create dataframes to calculate aggreements (Cohen's kappa)
		tpp_df=pd.DataFrame(index=tpp_proteins,columns=["interaction"])
		tpp_df["interaction"]=False
		tpp_df.loc[tpp_targets,"interaction"]=True
		ap_df=pd.DataFrame(index=tpp_proteins|slimp_proteins,columns=["interaction"])
		ap_df["interaction"]=False
		try:
			ap_df.loc[ap_targets,"interaction"]=True
		except KeyError:
			for p in ap_targets:
				ap_df.loc[p,"interaction"]=True
		slimp_df=pd.DataFrame(index=slimp_proteins,columns=["interaction"])
		slimp_df["interaction"]=False
		slimp_df.loc[slimp_targets,"interaction"]=True
		slimp80_df=pd.DataFrame(index=slimp_proteins,columns=["interaction"])
		slimp80_df["interaction"]=False
		slimp80_df.loc[slimp80_targets,"interaction"]=True
		#kappa df
		kappa_df=pd.DataFrame(index=["AP","TPP","SLIMP","SLIMP80"],columns=["AP","TPP","SLIMP","SLIMP80"])
		#pairwise dataframes
		##AP-SLIMP
		ap_slimp_df=pd.DataFrame(index=set(ap_df.index)&set(slimp_df.index),columns=["AP","SLIMP","SLIMP80"])
		ap_slimp_df.loc[ap_slimp_df.index,"AP"]=ap_df.loc[ap_slimp_df.index,"interaction"]
		ap_slimp_df.loc[ap_slimp_df.index,"SLIMP"]=slimp_df.loc[ap_slimp_df.index,"interaction"]
		ap_slimp_df.loc[ap_slimp_df.index,"SLIMP80"]=slimp80_df.loc[ap_slimp_df.index,"interaction"]
		kappa_df.loc["AP","SLIMP"]=metrics.cohen_kappa_score(ap_slimp_df["AP"],ap_slimp_df["SLIMP"])
		kappa_df.loc["AP","SLIMP80"]=metrics.cohen_kappa_score(ap_slimp_df["AP"],ap_slimp_df["SLIMP80"])
		##TPP-SLIMP
		tpp_slimp_df=pd.DataFrame(index=set(tpp_df.index)&set(slimp_df.index),columns=["TPP","SLIMP","SLIMP80"])
		tpp_slimp_df.loc[tpp_slimp_df.index,"TPP"]=tpp_df.loc[tpp_slimp_df.index,"interaction"]
		tpp_slimp_df.loc[tpp_slimp_df.index,"SLIMP"]=slimp_df.loc[tpp_slimp_df.index,"interaction"]
		tpp_slimp_df.loc[tpp_slimp_df.index,"SLIMP80"]=slimp80_df.loc[tpp_slimp_df.index,"interaction"]
		kappa_df.loc["TPP","SLIMP"]=metrics.cohen_kappa_score(tpp_slimp_df["TPP"],tpp_slimp_df["SLIMP"])
		kappa_df.loc["TPP","SLIMP80"]=metrics.cohen_kappa_score(tpp_slimp_df["TPP"],tpp_slimp_df["SLIMP80"])
		##AP-TPP
		ap_tpp_df=pd.DataFrame(index=set(ap_df.index)&set(tpp_df.index),columns=["AP","TPP"])
		ap_tpp_df.loc[ap_tpp_df.index,"AP"]=ap_df.loc[ap_tpp_df.index,"interaction"]
		ap_tpp_df.loc[ap_tpp_df.index,"TPP"]=tpp_df.loc[ap_tpp_df.index,"interaction"]
		kappa_df.loc["AP","TPP"]=metrics.cohen_kappa_score(ap_tpp_df["AP"],ap_tpp_df["TPP"])
		##SLIMP-SLIMP80
		slimp_df.loc[slimp_df.index,"SLIMP80"]=slimp80_df.loc[slimp_df.index,"interaction"]
		kappa_df.loc["SLIMP","SLIMP80"]=metrics.cohen_kappa_score(slimp_df["interaction"],slimp_df["SLIMP80"])
		kappa_df.to_csv("../overall_analysis/cohens_kappa_SerLeu_targets.tsv",sep="\t",index=True)
		#fisher's exact test whether prediction of two methods is significantly different from random
		##				SLIMP int	SLIMP not int
		# OTHER int		
		# OTHER not int	
		#contingency_table=[[intersect,other-intersect],[slimp-intersect,both not int]]
		#oddsratio = unconditional Maximum Likelihood Estimate
		oddsratio_df=pd.DataFrame(index=["AP","TPP","SLIMP","SLIMP80"],columns=["AP","TPP","SLIMP","SLIMP80"])
		pvalue_df=pd.DataFrame(index=["AP","TPP","SLIMP","SLIMP80"],columns=["AP","TPP","SLIMP","SLIMP80"])
		##SLIMP-AP
		ap_slimp80_df=ap_slimp_df[["AP","SLIMP80"]]
		ap_slimp_df=ap_slimp_df[["AP","SLIMP"]]
		ap_slimp_df["sum"]=ap_slimp_df.sum(axis=1)
		ap_slimp80_df["sum"]=ap_slimp80_df.sum(axis=1)
		contingency_table_slimp_ap=[[len(ap_slimp_df[ap_slimp_df["sum"]==2]),(ap_slimp_df[ap_slimp_df["AP"]==1]["SLIMP"]==0).sum()],[(ap_slimp_df[ap_slimp_df["SLIMP"]==1]["AP"]==0).sum(),len(ap_slimp_df[ap_slimp_df["sum"]==0])]]
		oddsratio_slimp_ap,pvalue_slimp_ap=fisher_exact(contingency_table_slimp_ap,alternative="greater")
		oddsratio_df.loc["SLIMP","AP"]=oddsratio_slimp_ap
		pvalue_df.loc["SLIMP","AP"]=pvalue_slimp_ap
		contingency_table_slimp80_ap=[[len(ap_slimp80_df[ap_slimp80_df["sum"]==2]),(ap_slimp80_df[ap_slimp80_df["AP"]==1]["SLIMP80"]==0).sum()],[(ap_slimp80_df[ap_slimp80_df["SLIMP80"]==1]["AP"]==0).sum(),len(ap_slimp80_df[ap_slimp80_df["sum"]==0])]]
		oddsratio_slimp80_ap,pvalue_slimp80_ap=fisher_exact(contingency_table_slimp80_ap,alternative="greater")
		oddsratio_df.loc["SLIMP80","AP"]=oddsratio_slimp80_ap
		pvalue_df.loc["SLIMP80","AP"]=pvalue_slimp80_ap
		##SLIMP-TPP
		tpp_slimp80_df=tpp_slimp_df[["TPP","SLIMP80"]]
		tpp_slimp_df=tpp_slimp_df[["TPP","SLIMP"]]
		tpp_slimp_df["sum"]=tpp_slimp_df.sum(axis=1)
		tpp_slimp80_df["sum"]=tpp_slimp80_df.sum(axis=1)
		contingency_table_slimp_tpp=[[len(tpp_slimp_df[tpp_slimp_df["sum"]==2]),(tpp_slimp_df[tpp_slimp_df["TPP"]==1]["SLIMP"]==0).sum()],[(tpp_slimp_df[tpp_slimp_df["SLIMP"]==1]["TPP"]==0).sum(),len(tpp_slimp_df[tpp_slimp_df["sum"]==0])]]
		oddsratio_slimp_tpp,pvalue_slimp_tpp=fisher_exact(contingency_table_slimp_tpp,alternative="two-sided")
		oddsratio_df.loc["SLIMP","TPP"]=oddsratio_slimp_tpp
		pvalue_df.loc["SLIMP","TPP"]=pvalue_slimp_tpp
		contingency_table_slimp80_tpp=[[len(tpp_slimp80_df[tpp_slimp80_df["sum"]==2]),(tpp_slimp80_df[tpp_slimp80_df["TPP"]==1]["SLIMP80"]==0).sum()],[(tpp_slimp80_df[tpp_slimp80_df["SLIMP80"]==1]["TPP"]==0).sum(),len(tpp_slimp80_df[tpp_slimp80_df["sum"]==0])]]
		oddsratio_slimp80_tpp,pvalue_slimp80_tpp=fisher_exact(contingency_table_slimp80_tpp,alternative="two-sided")
		oddsratio_df.loc["SLIMP80","TPP"]=oddsratio_slimp80_tpp
		pvalue_df.loc["SLIMP80","TPP"]=pvalue_slimp80_tpp
		#AP-TPP
		ap_tpp_df["sum"]=ap_tpp_df.sum(axis=1)
		contingency_table_ap_tpp=[[len(ap_tpp_df[ap_tpp_df["sum"]==2]),(ap_tpp_df[ap_tpp_df["AP"]==1]["TPP"]==0).sum()],[(ap_tpp_df[ap_tpp_df["TPP"]==1]["AP"]==0).sum(),len(ap_tpp_df[ap_tpp_df["sum"]==0])]]
		oddsratio_ap_tpp,pvalue_ap_tpp=fisher_exact(contingency_table_ap_tpp,alternative="greater")
		oddsratio_df.loc["AP","TPP"]=oddsratio_ap_tpp
		pvalue_df.loc["AP","TPP"]=pvalue_ap_tpp
		pvalue_df.to_csv("../overall_analysis/pvalues_SerLeu_targets.tsv",sep="\t",index=True)
		oddsratio_df.to_csv("../overall_analysis/oddsratios_SerLeu_targets.tsv",sep="\t",index=True)
		#plot
		f,axn=venn.venn3(data,names=["AP","TPP","SLIMP"])
		plt.show()
		f.savefig("../overall_analysis/SerLeu_targets_SLIMP.svg")
		f.savefig("../overall_analysis/SerLeu_targets_SLIMP.png")
		g,axn=venn.venn3(data80,names=["AP","TPP","SLIMP"])
		plt.show()
		g.savefig("../overall_analysis/SerLeu_targets_SLIMP80.svg")
		g.savefig("../overall_analysis/SerLeu_targets_SLIMP80.png")
		#prepare dataframe
		maxprots=np.max([len(tpp_proteins),len(slimp_proteins),len(ap_targets)])
		serleu_res=pd.DataFrame(index=range(maxprots),columns=["AP targets","AP targets in intersection","TPP targets","TPP all measured proteins","TPP targets in intersection","SLIMP targets","SLIMP all measured proteins","SLIMP targets in intersection"])
		serleu_res.loc[range(len(ap_targets)),"AP targets"]=list(ap_targets)
		serleu_res.loc[range(len(ap_targets_intersect)),"AP targets in intersection"]=list(ap_targets_intersect)
		serleu_res.loc[range(len(tpp_targets)),"TPP targets"]=list(tpp_targets)
		serleu_res.loc[range(len(tpp_targets_intersect)),"TPP targets in intersection"]=list(tpp_targets_intersect)
		serleu_res.loc[range(len(tpp_proteins)),"TPP all measured proteins"]=list(tpp_proteins)
		serleu_res.loc[range(len(slimp_targets)),"SLIMP targets"]=list(slimp_targets)
		serleu_res.loc[range(len(slimp_targets_intersect)),"SLIMP targets in intersection"]=list(slimp_targets_intersect)
		serleu_res.loc[range(len(slimp_proteins)),"SLIMP all measured proteins"]=list(slimp_proteins)
		#save to excel file
		outputfile=pd.ExcelWriter("../ForOla_SerLeu.xlsx")
		serleu_res.to_excel(outputfile,"SerLeu",index=False)
		outputfile.save()
		return
		
		
	
	#compares TyrAsp targets between AP, TPP and SLIMP for A.thaliana and for Tyr, Phe containing dipeptides
	def comparison_TyrPhe_TyrAsp_AP_TPP_SLIMP(self,predictions_ara):
		metabolite=19816752 #Tyr-Asp
		#load sheet/proteins
		inputfile=pd.ExcelFile("../ForBoris.xlsx")
		tyrasp=inputfile.parse("TyrAsp")
		slimp_proteins_uniprot=predictions_ara.index.get_level_values(0).unique()
		slimp_tyrphe_targets=set(tyrasp["SLIMP (Phe, Tyr)"].tolist())
		tyrasp_predictions=predictions_ara.reset_index().set_index("metabolite").loc[metabolite]
		tyrasp_targets=tyrasp_predictions[tyrasp_predictions["prediction"]==True]
		if self.proteinwise==True:
			slimp_targets_uniprot=set(tyrasp_targets[self.protcol].tolist())
		else:
			print("transform predictions from og-wise to proteinwise")
		ap_targets=set(tyrasp["AP"].tolist())
		tpp_targets=set(tyrasp["TPP (differential)"].tolist())
		tpp_proteins=set(tyrasp["TPP (all measured proteins)"].tolist())
		#translate IDs
		trans_df=da_ara.translate_proteins(slimp_proteins_uniprot.tolist(),in_id_type="ID",out_id_types=["ARAPORT_ID"],save=False)
		slimp_proteins=set(trans_df["ARAPORT_ID"].tolist())
		slimp_targets=set(trans_df.loc[slimp_targets_uniprot,"ARAPORT_ID"].tolist())
		#remove np.nans
		if np.nan in ap_targets:
			ap_targets.remove(np.nan)
		if np.nan in tpp_targets:
			tpp_targets.remove(np.nan)
		if np.nan in tpp_proteins:
			tpp_proteins.remove(np.nan)
		if np.nan in slimp_targets:
			slimp_targets.remove(np.nan)
		if np.nan in slimp_proteins:
			slimp_proteins.remove(np.nan)
		if np.nan in slimp_tyrphe_targets:
			slimp_tyrphe_targets.remove(np.nan)
		#get proteins in intersection among experiments
		tpp_targets_intersect=tpp_targets & slimp_proteins
		ap_targets_intersect=ap_targets & slimp_proteins & tpp_proteins #
		slimp_targets_intersect=slimp_targets & tpp_proteins
		slimp_tyrphe_targets_intersect=slimp_tyrphe_targets & tpp_proteins
		data=venn.get_labels([ap_targets_intersect,tpp_targets_intersect,slimp_targets_intersect])
		data_tyrphe=venn.get_labels([ap_targets_intersect,tpp_targets_intersect,slimp_tyrphe_targets_intersect])
		#create dataframes to calculate aggreements (Cohen's kappa)
		tpp_df=pd.DataFrame(index=tpp_proteins,columns=["interaction"])
		tpp_df["interaction"]=False
		tpp_df.loc[tpp_targets,"interaction"]=True
		ap_df=pd.DataFrame(index=tpp_proteins|slimp_proteins,columns=["interaction"])
		ap_df["interaction"]=False
		try:
			ap_df.loc[ap_targets,"interaction"]=True
		except KeyError:
			for p in ap_targets:
				ap_df.loc[p,"interaction"]=True
		slimp_df=pd.DataFrame(index=slimp_proteins,columns=["interaction"])
		slimp_df["interaction"]=False
		slimp_df.loc[slimp_targets,"interaction"]=True
		slimp_tyrphe_df=pd.DataFrame(index=slimp_proteins,columns=["interaction"])
		slimp_tyrphe_df["interaction"]=False
		slimp_tyrphe_df.loc[slimp_tyrphe_targets & slimp_proteins,"interaction"]=True
		#kappa df
		kappa_df=pd.DataFrame(index=["AP","TPP","SLIMP","SLIMP(Tyr,Phe)"],columns=["AP","TPP","SLIMP","SLIMP(Tyr,Phe)"])
		#pairwise dataframes
		##AP-SLIMP
		ap_slimp_df=pd.DataFrame(index=set(ap_df.index)&set(slimp_df.index),columns=["AP","SLIMP","SLIMP(Tyr,Phe)"])
		ap_slimp_df.loc[ap_slimp_df.index,"AP"]=ap_df.loc[ap_slimp_df.index,"interaction"]
		ap_slimp_df.loc[ap_slimp_df.index,"SLIMP"]=slimp_df.loc[ap_slimp_df.index,"interaction"]
		ap_slimp_df.loc[ap_slimp_df.index,"SLIMP(Tyr,Phe)"]=slimp_tyrphe_df.loc[ap_slimp_df.index,"interaction"]
		kappa_df.loc["AP","SLIMP"]=metrics.cohen_kappa_score(ap_slimp_df["AP"],ap_slimp_df["SLIMP"])
		kappa_df.loc["AP","SLIMP(Tyr,Phe)"]=metrics.cohen_kappa_score(ap_slimp_df["AP"],ap_slimp_df["SLIMP(Tyr,Phe)"])
		##TPP-SLIMP
		tpp_slimp_df=pd.DataFrame(index=set(tpp_df.index)&set(slimp_df.index),columns=["TPP","SLIMP","SLIMP(Tyr,Phe)"])
		tpp_slimp_df.loc[tpp_slimp_df.index,"TPP"]=tpp_df.loc[tpp_slimp_df.index,"interaction"]
		tpp_slimp_df.loc[tpp_slimp_df.index,"SLIMP"]=slimp_df.loc[tpp_slimp_df.index,"interaction"]
		tpp_slimp_df.loc[tpp_slimp_df.index,"SLIMP(Tyr,Phe)"]=slimp_tyrphe_df.loc[tpp_slimp_df.index,"interaction"]
		kappa_df.loc["TPP","SLIMP"]=metrics.cohen_kappa_score(tpp_slimp_df["TPP"],tpp_slimp_df["SLIMP"])
		kappa_df.loc["TPP","SLIMP(Tyr,Phe)"]=metrics.cohen_kappa_score(tpp_slimp_df["TPP"],tpp_slimp_df["SLIMP(Tyr,Phe)"])
		##AP-TPP
		ap_tpp_df=pd.DataFrame(index=set(ap_df.index)&set(tpp_df.index),columns=["AP","TPP"])
		ap_tpp_df.loc[ap_tpp_df.index,"AP"]=ap_df.loc[ap_tpp_df.index,"interaction"]
		ap_tpp_df.loc[ap_tpp_df.index,"TPP"]=tpp_df.loc[ap_tpp_df.index,"interaction"]
		kappa_df.loc["AP","TPP"]=metrics.cohen_kappa_score(ap_tpp_df["AP"],ap_tpp_df["TPP"])
		##SLIMP-SLIMP(Tyr,Phe)
		slimp_df.loc[slimp_df.index,"SLIMP(Tyr,Phe)"]=slimp_tyrphe_df.loc[slimp_df.index,"interaction"]
		kappa_df.loc["SLIMP","SLIMP(Tyr,Phe)"]=metrics.cohen_kappa_score(slimp_df["interaction"],slimp_df["SLIMP(Tyr,Phe)"])
		kappa_df.to_csv("../overall_analysis/cohens_kappa_TyrAsp_targets.tsv",sep="\t",index=True)
		#fisher's exact test whether prediction of two methods is significantly different from random
		##				SLIMP int	SLIMP not int
		# OTHER int		
		# OTHER not int	
		#contingency_table=[[intersect,other-intersect],[slimp-intersect,both not int]]
		#oddsratio = unconditional Maximum Likelihood Estimate
		oddsratio_df=pd.DataFrame(index=["AP","TPP","SLIMP","SLIMP(Tyr,Phe)"],columns=["AP","TPP","SLIMP","SLIMP(Tyr,Phe)"])
		pvalue_df=pd.DataFrame(index=["AP","TPP","SLIMP","SLIMP(Tyr,Phe)"],columns=["AP","TPP","SLIMP","SLIMP(Tyr,Phe)"])
		##SLIMP-AP
		ap_slimp_tyrphe_df=ap_slimp_df[["AP","SLIMP(Tyr,Phe)"]]
		ap_slimp_df=ap_slimp_df[["AP","SLIMP"]]
		ap_slimp_df["sum"]=ap_slimp_df.sum(axis=1)
		ap_slimp_tyrphe_df["sum"]=ap_slimp_tyrphe_df.sum(axis=1)
		contingency_table_slimp_ap=[[len(ap_slimp_df[ap_slimp_df["sum"]==2]),(ap_slimp_df[ap_slimp_df["AP"]==1]["SLIMP"]==0).sum()],[(ap_slimp_df[ap_slimp_df["SLIMP"]==1]["AP"]==0).sum(),len(ap_slimp_df[ap_slimp_df["sum"]==0])]]
		oddsratio_slimp_ap,pvalue_slimp_ap=fisher_exact(contingency_table_slimp_ap,alternative="greater")
		oddsratio_df.loc["SLIMP","AP"]=oddsratio_slimp_ap
		pvalue_df.loc["SLIMP","AP"]=pvalue_slimp_ap
		contingency_table_slimp_tyrphe_ap=[[len(ap_slimp_tyrphe_df[ap_slimp_tyrphe_df["sum"]==2]),(ap_slimp_tyrphe_df[ap_slimp_tyrphe_df["AP"]==1]["SLIMP(Tyr,Phe)"]==0).sum()],[(ap_slimp_tyrphe_df[ap_slimp_tyrphe_df["SLIMP(Tyr,Phe)"]==1]["AP"]==0).sum(),len(ap_slimp_tyrphe_df[ap_slimp_tyrphe_df["sum"]==0])]]
		oddsratio_slimp_tyrphe_ap,pvalue_slimp_tyrphe_ap=fisher_exact(contingency_table_slimp_tyrphe_ap,alternative="greater")
		oddsratio_df.loc["SLIMP(Tyr,Phe)","AP"]=oddsratio_slimp_tyrphe_ap
		pvalue_df.loc["SLIMP(Tyr,Phe)","AP"]=pvalue_slimp_tyrphe_ap
		##SLIMP-TPP
		tpp_slimp_tyrphe_df=tpp_slimp_df[["TPP","SLIMP(Tyr,Phe)"]]
		tpp_slimp_df=tpp_slimp_df[["TPP","SLIMP"]]
		tpp_slimp_df["sum"]=tpp_slimp_df.sum(axis=1)
		tpp_slimp_tyrphe_df["sum"]=tpp_slimp_tyrphe_df.sum(axis=1)
		contingency_table_slimp_tpp=[[len(tpp_slimp_df[tpp_slimp_df["sum"]==2]),(tpp_slimp_df[tpp_slimp_df["TPP"]==1]["SLIMP"]==0).sum()],[(tpp_slimp_df[tpp_slimp_df["SLIMP"]==1]["TPP"]==0).sum(),len(tpp_slimp_df[tpp_slimp_df["sum"]==0])]]
		oddsratio_slimp_tpp,pvalue_slimp_tpp=fisher_exact(contingency_table_slimp_tpp,alternative="greater")
		oddsratio_df.loc["SLIMP","TPP"]=oddsratio_slimp_tpp
		pvalue_df.loc["SLIMP","TPP"]=pvalue_slimp_tpp
		contingency_table_slimp_tyrphe_tpp=[[len(tpp_slimp_tyrphe_df[tpp_slimp_tyrphe_df["sum"]==2]),(tpp_slimp_tyrphe_df[tpp_slimp_tyrphe_df["TPP"]==1]["SLIMP(Tyr,Phe)"]==0).sum()],[(tpp_slimp_tyrphe_df[tpp_slimp_tyrphe_df["SLIMP(Tyr,Phe)"]==1]["TPP"]==0).sum(),len(tpp_slimp_tyrphe_df[tpp_slimp_tyrphe_df["sum"]==0])]]
		oddsratio_slimp_tyrphe_tpp,pvalue_slimp_tyrphe_tpp=fisher_exact(contingency_table_slimp_tyrphe_tpp,alternative="greater")
		oddsratio_df.loc["SLIMP(Tyr,Phe)","TPP"]=oddsratio_slimp_tyrphe_tpp
		pvalue_df.loc["SLIMP(Tyr,Phe)","TPP"]=pvalue_slimp_tyrphe_tpp
		#AP-TPP
		ap_tpp_df["sum"]=ap_tpp_df.sum(axis=1)
		contingency_table_ap_tpp=[[len(ap_tpp_df[ap_tpp_df["sum"]==2]),(ap_tpp_df[ap_tpp_df["AP"]==1]["TPP"]==0).sum()],[(ap_tpp_df[ap_tpp_df["TPP"]==1]["AP"]==0).sum(),len(ap_tpp_df[ap_tpp_df["sum"]==0])]]
		oddsratio_ap_tpp,pvalue_ap_tpp=fisher_exact(contingency_table_ap_tpp,alternative="greater")
		oddsratio_df.loc["TPP","AP"]=oddsratio_ap_tpp
		pvalue_df.loc["TPP","AP"]=pvalue_ap_tpp
		pvalue_df.to_csv("../overall_analysis/pvalues_TyrAsp_targets.tsv",sep="\t",index=True)
		oddsratio_df.to_csv("../overall_analysis/oddsratios_TyrAsp_targets.tsv",sep="\t",index=True)
		#plot
		f,axn=venn.venn3(data,names=["AP","TPP","SLIMP"])
		plt.show()
		f.savefig("../overall_analysis/TyrAsp_targets_SLIMP.svg")
		f.savefig("../overall_analysis/TyrAsp_targets_SLIMP.png")
		g,axn=venn.venn3(data_tyrphe,names=["AP","TPP","SLIMP"])
		plt.show()
		g.savefig("../overall_analysis/TyrPhe_targets_SLIMP.svg")
		g.savefig("../overall_analysis/TyrPhe_targets_SLIMP.png")
		#prepare dataframe
		maxprots=np.max([len(tpp_proteins),len(slimp_proteins),len(ap_targets)])
		tyrasp_res=pd.DataFrame(index=range(maxprots),columns=["AP targets","AP targets in intersection","TPP targets","TPP all measured proteins","TPP targets in intersection","SLIMP targets","SLIMP all measured proteins","SLIMP targets in intersection","SLIMP targets (Tyr,Phe)","SLIMP targets (Tyr,Phe) in intersection"])
		tyrasp_res.loc[range(len(ap_targets)),"AP targets"]=list(ap_targets)
		tyrasp_res.loc[range(len(ap_targets_intersect)),"AP targets in intersection"]=list(ap_targets_intersect)
		tyrasp_res.loc[range(len(tpp_targets)),"TPP targets"]=list(tpp_targets)
		tyrasp_res.loc[range(len(tpp_targets_intersect)),"TPP targets in intersection"]=list(tpp_targets_intersect)
		tyrasp_res.loc[range(len(tpp_proteins)),"TPP all measured proteins"]=list(tpp_proteins)
		tyrasp_res.loc[range(len(slimp_targets)),"SLIMP targets"]=list(slimp_targets)
		tyrasp_res.loc[range(len(slimp_targets_intersect)),"SLIMP targets in intersection"]=list(slimp_targets_intersect)
		tyrasp_res.loc[range(len(slimp_proteins)),"SLIMP all measured proteins"]=list(slimp_proteins)
		tyrasp_res.loc[range(len(slimp_tyrphe_targets)),"SLIMP targets (Tyr,Phe)"]=list(slimp_tyrphe_targets)
		tyrasp_res.loc[range(len(slimp_tyrphe_targets_intersect)),"SLIMP targets (Tyr,Phe) in intersection"]=list(slimp_tyrphe_targets_intersect)
		#save to excel file
		outputfile=pd.ExcelWriter("../ForOla_TyrAsp.xlsx")
		tyrasp_res.to_excel(outputfile,"TyrAsp",index=False)
		outputfile.save()
		return
	
	
	
	#calculate how many of Ewkas experimental verified interactions were found in Arabidopsis predictions
	def comparison_to_Ewkas_set(self,normalized=None,feature=None,method=np.mean,approach=None,classifier="svm_clf(C=1,kernel='linear')",confidence=None):
		db_orgs="human_yeast_ecoli"
		if feature is None:
			feature=self.feature
		if approach is None:
			approach=self.approach
		if normalized is None:
			normalized=self.normalized
		if confidence is None:
			if "high" in approach:
				confidence="high"
			else:
				confidence="low"
		#load Ewkas set
		known_positives=pd.read_csv("../databases/known_protein-metabolites-interactions_for_atha_prots.txt",sep="\t")
		known_positives["CID_int"]=list(map(lambda x: int(x[4:]),known_positives["CID"]))
		#rename column protein to STRING_ID
		known_positives.rename({"protein":"STRING_ID"},inplace=True)
		if self.proteinwise==False:
			#load string-OG database and known interactions
			string_OG_db=pd.read_csv(self.databases+"ortholog_mapping.tsv",sep="\t").set_index("Unnamed: 0")
			string_OG_db=string_OG_db.replace(np.nan,"")
			string_OG_db["string_ids"]=string_OG_db[string_OG_db.columns[1:]].sum(axis=1)
			known_positives["OG"]=[[] for i in range(len(known_positives))]
			for i in known_positives.index:
				for j in string_OG_db.index:
					if known_positives.loc[i,"STRING_ID"] in string_OG_db.loc[j,"string_ids"].split("."): #known_positives["protein"] is String ID
						known_positives.loc[i,"OG"].append(j)
		else:
			id_db=pd.read_csv("../databases/A.thaliana_ID_collection.tsv",sep="\t").set_index("STRING_ID")
			for i in known_positives.index:
				known_positives.loc[i,"protein"]=id_db.loc[known_positives.loc[i,"STRING_ID"].split(".")[-1],"ID"]
		#extract known OG-metabolite interactions
		known_ints=pd.DataFrame(columns=[self.protcol,"CID"])
		for i in known_positives.index:
			if len(known_positives.loc[i,self.protcol])>0:
				for og in known_positives.loc[i,self.protcol]:
					known_ints=known_ints.append(pd.DataFrame(data=[[og,known_positives.loc[i,"CID_int"]]],columns=[self.protcol,"CID"]))
		known_ints=known_ints.set_index([self.protcol,"CID"]) #51 interactions
		#load positive set, complete set and predictions => Ara_only2
		positive_set=pd.read_csv(self.databases+"positive_merged_"+db_orgs+"_set_"+normalized+"normprofiles"+feature+"_"+method.__name__+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		positive_set_trimmed=pd.read_csv(self.databases+"positive_merged_"+db_orgs+"_set_"+normalized+"normprofiles"+feature+"_trimmed_"+confidence+"_confidence_"+method.__name__+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		complete_set=pd.read_csv(self.databases+"complete_set_"+normalized+"normprofiles"+feature+"_"+method.__name__+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+feature+approach+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		prediction_probabilities=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+feature+approach+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		###get intersection to data sets
		#intersection between ewkas interactions and experiment
		intersecting_ints=set(set(known_ints.index)&set(predictions.index)) 
		print(str(predictions.loc[list(intersecting_ints),classifier].sum())+" out of "+str(len(intersecting_ints))+" interactions found experimentally predicted correctly, equalling "+str(np.round(predictions.loc[list(intersecting_ints),classifier].sum()/len(intersecting_ints)*100,2))+"%") #including training data
		#intersection between ewkas interactions and positive set not in training
		pos_intersecting_ints=set(set(known_ints.index)&set(set(positive_set.index.tolist()))-set(positive_set_trimmed.index.tolist())) 
		print(str(predictions.loc[list(pos_intersecting_ints),classifier].sum())+" out of "+str(len(pos_intersecting_ints))+" interactions occuring in positive set, but not in training set were predicted correctly, equalling "+str(np.round(predictions.loc[list(pos_intersecting_ints),classifier].sum()/len(pos_intersecting_ints)*100,2))+"%") #interactions occuring in positive set, but not in training
		postrimmed_intersecting_ints=set(set(known_ints.index.tolist()) & set(positive_set_trimmed.index.tolist())) #interactions seen in training data
		#intersection between ewkas interactions and interactions not already seen in training
		not_seen_ints=intersecting_ints-postrimmed_intersecting_ints # interactions not seen by classifier
		#predictions.loc[list(pos_intersecting_ints),classifier].sum()/len(pos_intersecting_ints) #positive set not in training
		#predictions.loc[list(postrimmed_intersecting_ints),classifier].sum()/len(postrimmed_intersecting_ints) #interactions seen in training
		print(str(predictions.loc[list(not_seen_ints),classifier].sum())+" out of "+str(len(not_seen_ints))+" interactions not seen in training were predicted correctly, equalling "+str(np.round(predictions.loc[list(not_seen_ints),classifier].sum()/len(not_seen_ints)*100,2))+"%") #interactions not in training
		print(prediction_probabilities.loc[list(intersecting_ints),classifier]) #how sure about decision
		#figure
		ml_not_seen_intersecting_predictions=predictions.loc[list(not_seen_ints),classifier]
		ml_not_seen_intersecting_ints=set(ml_not_seen_intersecting_predictions[ml_not_seen_intersecting_predictions==True].index)
		data=venn.get_labels([ml_not_seen_intersecting_ints,not_seen_ints])
		f,axn=venn.venn2(data,names=["ML","Ewkas Set"])
		plt.show()
		f.savefig("../overall_analysis/venn_diagram_ara_comparison_to_Ewka_"+confidence+"_confidence.png")
		f.savefig("../overall_analysis/venn_diagram_ara_comparison_to_Ewka_"+confidence+"_confidence.svg")
		#pie
		pf=plt.figure()
		plt.pie([len(ml_not_seen_intersecting_ints),len(not_seen_ints)-len(ml_not_seen_intersecting_ints)],labels=["interaction predicted","no interaction predicted"],colors=["tab:blue","tab:red"],autopct="%1.1f%%")
		plt.show()
		pf.savefig("../overall_analysis/pie_chart_ara_comparison_to_Ewka_"+confidence+"_confidence.png")
		pf.savefig("../overall_analysis/pie_chart_ara_comparison_to_Ewka_"+confidence+"_confidence.svg")
		return
	
	
	
	#compares all metabolite targets between Promis, slimp and promis on slimp dataset
	def comparison_to_PROMIS(self,org):
		if org=="ath":
			promis_raw=pd.read_csv("../PROMIS/Raw_overlapping_interactions.txt",sep="\t")#names
			promis_concat=pd.read_csv("../PROMIS/Concat_overlapping_interactions.txt",sep="\t") #stitch
			slimp_concat=pd.read_csv("../PROMIS/Ath_ranked_protein_targets_greaterequal_0.5_prediction_probability_amin_same_crosscorr_binsize4_not_normalized_maxmaxprofiles_unbalanced_low_confidence.tsv",sep="\t")
			#trans_df=pd.read_csv("../PROMIS/Metabolite_names_overlapping.txt",sep="\t",dtype={"PubchemID":str}).set_index("PubchemID")
		elif org=="sce":
			promis_raw=pd.read_csv("../PROMIS/Raw_Overlapping_Interactions_Yeast.txt",sep="\t")#names
			promis_concat=pd.read_csv("../PROMIS/Concat_Overlapping_Interactions_Yeast.txt",sep="\t") #stitch
			slimp_concat=pd.read_csv("../PROMIS/Sce_ranked_protein_targets_greaterequal_0.5_prediction_probability_amin_same_crosscorr_binsize4_not_normalized_maxmaxprofiles_unbalanced_low_confidence.tsv",sep="\t")
		#get dataframe with names for stitch cids
		slimp_concat_cols_pubchem=self.translate_CIDs_from_Stitch_to_PubChem(slimp_concat.columns.tolist())
		slimp_trans_df=slimp_concat_cols_pubchem.set_index("Stitch_CID")
		#slimp_concat.columns=list(map(lambda x: str(x),slimp_concat_cols_pubchem["PubChem_CID"].tolist()))
		#translate from names to stitch cids
		promis_raw_cols_stitch=self.translate_metabolites_offline(promis_raw.columns.tolist())
		promis_raw_columns=pd.DataFrame(columns=["Stitch_CID"],index=promis_raw.columns.tolist())
		for m in promis_raw.columns:
			try:
				promis_raw_columns.loc[m,"Stitch_CID"]=str(int(promis_raw_cols_stitch.loc[m,"CID"]))
			except KeyError:
				promis_raw_columns.loc[m,"Stitch_CID"]=m
		promis_raw.columns=promis_raw_columns["Stitch_CID"].tolist()
		intersect=set(promis_raw.columns.tolist()) & set(promis_concat.columns.tolist()) & set(slimp_concat.columns.tolist())
		labels=["promis_raw","promis_concat","slimp"]
		for metabolite in intersect:
			ints=[set(promis_raw[metabolite].dropna().tolist()),set(promis_concat[metabolite].dropna().tolist()),set(slimp_concat[metabolite].dropna().tolist())]
			data=venn.get_labels(ints)
			f,axn=eval("venn.venn"+str(len(ints))+"(data,names=labels)")
			plt.title(slimp_trans_df.loc[int(metabolite),"Name"])
			#plt.show()
			f.savefig("../PROMIS/venns_sce/sce_venn_"+metabolite+".svg")
		return
	
	
	
	#calculate accuracy on promis
	def calculate_accuracy_on_PROMIS(self,org,method=np.min):
		training_set=pd.read_csv(self.analyses+"training_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t",dtype={"metabolite":str}).set_index([self.protcol,"metabolite"])
		#training_set=pd.read_csv(mlc.analyses+"training_set"+mlc.feature+mlc.approach+"_"+method.__name__+".tsv",sep="\t",dtype={"metabolite":str}).set_index([mlc.protcol,"metabolite"])
		if org=="ath":
			promis_raw=pd.read_csv("../PROMIS/Raw_overlapping_interactions.txt",sep="\t")#names
			promis_concat=pd.read_csv("../PROMIS/Concat_overlapping_interactions.txt",sep="\t") #stitch
			slimp_concat=pd.read_csv("../PROMIS/Ath_ranked_protein_targets_greaterequal_0.5_prediction_probability_amin_same_crosscorr_binsize4_not_normalized_maxmaxprofiles_unbalanced_low_confidence.tsv",sep="\t")
			#trans_df=pd.read_csv("../PROMIS/Metabolite_names_overlapping.txt",sep="\t",dtype={"PubchemID":str}).set_index("PubchemID")
		elif org=="sce":
			promis_raw=pd.read_csv("../PROMIS/Raw_Overlapping_Interactions_Yeast.txt",sep="\t")#names
			promis_concat=pd.read_csv("../PROMIS/Concat_Overlapping_Interactions_Yeast.txt",sep="\t") #stitch
			slimp_concat=pd.read_csv("../PROMIS/Sce_ranked_protein_targets_greaterequal_0.5_prediction_probability_amin_same_crosscorr_binsize4_not_normalized_maxmaxprofiles_unbalanced_low_confidence.tsv",sep="\t")
		#translate from names to stitch cids
		promis_raw_cols_stitch=self.translate_metabolites_offline(promis_raw.columns.tolist())
		promis_raw_columns=pd.DataFrame(columns=["Stitch_CID"],index=promis_raw.columns.tolist())
		for m in promis_raw.columns:
			try:
				promis_raw_columns.loc[m,"Stitch_CID"]=str(int(promis_raw_cols_stitch.loc[m,"CID"]))
			except KeyError:
				promis_raw_columns.loc[m,"Stitch_CID"]=m
		promis_raw.columns=promis_raw_columns["Stitch_CID"].tolist()
		tp_raw=0
		tn_raw=0
		fp_raw=0
		tn_notpred_raw=0
		fn_raw=0
		fn_notpred_raw=0
		tp_concat=0
		tn_concat=0
		fp_concat=0
		tn_notpred_concat=0
		fn_concat=0
		fn_notpred_concat=0
		missing_metas_raw=set()
		missing_metas_concat=set()
		#for every instance in training set, check if it is present in PROMIS targets
		for pair in training_set.index:
			#PROMIS_raw
			if pair[1] in promis_raw.columns and pair[0] in promis_raw[pair[1]].tolist(): #Positive test result
				if training_set.loc[pair,"interaction"]==True:
					tp_raw+=1
				elif training_set.loc[pair,"interaction"]==False:
					fp_raw+=1
			elif pair[1] in promis_raw.columns and pair[0] not in promis_raw[pair[1]].tolist(): #Negative test result
				if training_set.loc[pair,"interaction"]==True:
					fn_raw+=1
				elif training_set.loc[pair,"interaction"]==False:
					tn_raw+=1
			elif pair[1] not in promis_raw.columns: #either not predicted or not tested
				print(pair[1]+" not in promis_raw")
				missing_metas_raw.add(pair[1])
				#case "not predicted"
				if training_set.loc[pair,"interaction"]==True:
					fn_notpred_raw+=1
				elif training_set.loc[pair,"interaction"]==False:
					tn_notpred_raw+=1
			else:
				print("dum dum")
			#PROMIS_concat
			if pair[1] in promis_concat.columns and pair[0] in promis_concat[pair[1]].tolist(): #Positive test result
				if training_set.loc[pair,"interaction"]==True:
					tp_concat+=1
				elif training_set.loc[pair,"interaction"]==False:
					fp_concat+=1
			elif pair[1] in promis_concat.columns and pair[0] not in promis_concat[pair[1]].tolist(): #Negative test result
				if training_set.loc[pair,"interaction"]==True:
					fn_concat+=1
				elif training_set.loc[pair,"interaction"]==False:
					tn_concat+=1
			elif pair[1] not in promis_concat.columns: #either not predicted or not tested
				print(pair[1]+" not in promis_concat")
				missing_metas_concat.add(pair[1])
				#case "not predicted"
				if training_set.loc[pair,"interaction"]==True:
					fn_notpred_concat+=1
				elif training_set.loc[pair,"interaction"]==False:
					tn_notpred_concat+=1
			else:
				print("dum dum concat")
		#fn_notpred_raw for the case, that missing metabolites are treated as prediction of absence of interaction, fn_raw as not tested (=> excluded)
		fn_np_raw=fn_raw+fn_notpred_raw
		fn_np_concat=fn_concat+fn_notpred_concat
		tn_np_raw=tn_raw+tn_notpred_raw
		tn_np_concat=tn_concat+tn_notpred_concat
		#calculate performance measures
		acc_raw=(tp_raw+tn_raw)/(tp_raw+tn_raw+fp_raw+fn_raw)
		acc_np_raw=(tp_raw+tn_np_raw)/(tp_raw+tn_np_raw+fp_raw+fn_np_raw)
		acc_concat=(tp_concat+tn_concat)/(tp_concat+tn_concat+fp_concat+fn_concat)
		acc_np_concat=(tp_concat+tn_np_concat)/(tp_concat+tn_np_concat+fp_concat+fn_np_concat)
		sens_raw=tp_raw/(tp_raw+fn_raw)
		sens_concat=tp_concat/(tp_concat+fn_concat)
		sens_np_raw=tp_raw/(tp_raw+fn_np_raw)
		sens_np_concat=tp_concat/(tp_concat+fn_np_concat)
		spec_raw=tn_raw/(tn_raw+fp_raw)
		spec_concat=tn_concat/(tn_concat+fp_concat)
		spec_np_raw=tn_np_raw/(tn_np_raw+fp_raw)
		spec_np_concat=tn_np_concat/(tn_np_concat+fp_concat)
		f1_raw=(2*tp_raw)/(2*tp_raw+fp_raw+fn_raw)
		f1_concat=(2*tp_concat)/(2*tp_concat+fp_concat+fn_concat)
		f1_np_raw=(2*tp_raw)/(2*tp_raw+fp_raw+fn_np_raw)
		f1_np_concat=(2*tp_concat)/(2*tp_concat+fp_concat+fn_np_concat)
		prec_raw=tp_raw/(tp_raw+fp_raw)
		prec_concat=tp_concat/(tp_concat+fp_concat)
		#results into dataframe
		df=pd.DataFrame(index=["accuracy","f1_score","precision","sensitivity","specificity"],columns=["PROMIS","PROMIS_excluded","PROMIS_concat","PROMIS_concat_excluded"])
		df.loc["accuracy"]=[acc_np_raw,acc_raw,acc_np_concat,acc_concat]
		df.loc["f1_score"]=[f1_np_raw,f1_raw,f1_np_concat,f1_concat]
		df.loc["precision"]=[prec_raw,prec_raw,prec_concat,prec_concat]
		df.loc["sensitivity"]=[sens_np_raw,sens_raw,sens_np_concat,sens_concat]
		df.loc["specificity"]=[spec_np_raw,spec_raw,spec_np_concat, spec_concat]
		df.to_csv("../PROMIS/performance_measures_"+org+".tsv",sep="\t",index=True,header=True)
		return
	
	
	
	#retrieve prediction results for given (df) interactions
	def verify_selected_interactions_og(self,predictions,df=None):
		if df is None:
			protein_abbrev=["GAPDH","PEPCK","ALD","PNP1","FBA6","PGK1","KPHMT1","RBP47B","MTI"]
			metabolite=["Tyr-Asp","Ala-Ile","Gly-Phe","Xanthine","Gly-Phe","Ser-Leu","Panthotenic acid","2'-3'-cAMP","MTA"]
			cid=[19816752,417358,92953,1188,92953,13919048,6613,101812,439176]
			uniprot_ath=["Q9SAJ6;Q5E924;P25858;P25857;P25856;Q9FX54;Q9LPW0","Q9T074","Q8W033;Q70DU8;Q9SU63;Q56YU0;Q8S528;Q0WSF1;Q9SYG7;Q70E96;F4JC27","","Q9SJQ9","Q9LD57;P50318;Q9SAJ4","O82357;Q9M315","Q0WW84;P42731;Q9FXA2;Q05196;Q9SAB3;O22173;Q9LX90;O64380;Q93VI4;F4I3B3;Q9SX79","Q9ZUG4"]
			uniprot_sce=["P00359;P00360;P00358","P10963","P22281;P32872;E9P9H2","Q05788","P14540","P00560","P38122","P04147","Q06489;A6ZWZ9;C8ZJE0;B3LK82;C76XT0;B5VTQ8"]
			df=pd.DataFrame(columns=["protein abbreviation","metabolite","OG","CID","UniProt A.thaliana","UniProt S.cerevisiae"])
			df["protein abbreviation"]=protein_abbrev
			df["metabolite"]=metabolite
			df["CID"]=cid
			df["UniProt A.thaliana"]=uniprot_ath
			df["UniProt S.cerevisiae"]=uniprot_sce
		#annotate OGs for Arabidopsis and yeast proteins
		ath_og_trans=pd.read_csv("../databases/A.thaliana_string_uniprot_cog.tsv",sep="\t").set_index("UniProt_ID")
		sce_og_trans=pd.read_csv("../databases/S.cerevisiae_string_uniprot_cog.tsv",sep="\t").set_index("UniProt_ID")
		ath_ogs=list()
		sce_ogs=list()
		for i in df.index:
			#print(df.loc[i,"protein abbreviation"])
			#print("Arabidopsis")
			ath_og=list()
			ath_uniprots=list()
			for ath_uniprot in df.loc[i,"UniProt A.thaliana"].split(";"):
				try:
					ath_og.append(ath_og_trans.loc[ath_uniprot,"COG"])
					ath_uniprots.append(ath_uniprot)
				except:
					pass
					#print(ath_uniprot)
			ath_ogs.append(";".join(set(ath_og)))
			df.loc[i,"UniProt A.thaliana"]=";".join(ath_uniprots)
			#print("Saccharomyces")
			sce_og=list()
			sce_uniprots=list()
			for sce_uniprot in df.loc[i,"UniProt S.cerevisiae"].split(";"):
				try:
					sce_og.append(sce_og_trans.loc[sce_uniprot,"COG"])
					sce_uniprots.append(sce_uniprot)
				except:
					pass
					#print(sce_uniprot)
			sce_ogs.append(";".join(set(sce_og)))
			df.loc[i,"UniProt S.cerevisiae"]=";".join(sce_uniprots)
		df["OG A.thaliana"]=ath_ogs
		df["OG S.cerevisiae"]=sce_ogs
		#get union of OGs over Ath and Sce
		for i in df.index:
			ogs=set(set(df.loc[i,"OG A.thaliana"].split(";")) | set(df.loc[i,"OG S.cerevisiae"].split(";")))
			if "" in ogs:
				ogs.remove("")
			df.loc[i,"OG"]=";".join(ogs)
		#prediction result
		for i in df.index:
			predictions_og=list()
			cid=df.loc[i,"CID"]
			for og in df.loc[i,"OG"].split(";"):
				if og in predictions.index.get_level_values(0) and cid in predictions.index.get_level_values(1):
					predictions_og.append(predictions.loc[(og,cid),"prediction"])
			df.loc[i,"predictions"]=";".join(list(map(lambda x: str(x),predictions_og)))
		df.to_csv(self.analyses+"prediction_verification"+self.feature+self.approach+".tsv",sep="\t",header=True,index=True)
		return
	
	
	
	#extracts predictions of interest
	#predictions is df with "prediction column",organism ={"A.thaliana","S.cerevisiae"}
	def verify_selected_interactions_protwise(self,predictions,df=None,organism=None,appendix=""):
		if df is None:
			protein_abbrev=["GAPDH","PEPCK","ALD","PNP1","FBA6","PGK1","KPHMT1","RBP47B","MTI"]
			metabolite=["Tyr-Asp","Ala-Ile","Gly-Phe","Xanthine","Gly-Phe","Ser-Leu","Panthotenic acid","2'-3'-cAMP","MTA"]
			cid=[19816752,417358,92953,1188,92953,13919048,6613,101812,439176]
			if organism=="A.thaliana":
				uniprot=["Q9SAJ6;Q5E924;P25858;P25857;P25856;Q9FX54;Q9LPW0","Q9T074","Q8W033;Q70DU8;Q9SU63;Q56YU0;Q8S528;Q0WSF1;Q9SYG7;Q70E96;F4JC27","","Q9SJQ9","Q9LD57;P50318;Q9SAJ4","O82357;Q9M315","Q0WW84;P42731;Q9FXA2;Q05196;Q9SAB3;O22173;Q9LX90;O64380;Q93VI4;F4I3B3;Q9SX79","Q9ZUG4"]
			elif organism=="S.cerevisiae":
				uniprot=["P00359;P00360;P00358","P10963","P22281;P32872;E9P9H2","Q05788","P14540","P00560","P38122","P04147","Q06489;A6ZWZ9;C8ZJE0;B3LK82;C76XT0;B5VTQ8"]
			protein_abbrev_protwise=list(itertools.chain.from_iterable(map(lambda x: [protein_abbrev[x]]*len(uniprot[x].split(";")),range(len(protein_abbrev)))))
			metabolite_protwise=list(itertools.chain.from_iterable(map(lambda x: [metabolite[x]]*len(uniprot[x].split(";")),range(len(metabolite)))))
			cid_protwise=list(itertools.chain.from_iterable(map(lambda x: [cid[x]]*len(uniprot[x].split(";")),range(len(cid)))))
			df=pd.DataFrame(columns=["protein abbreviation","metabolite","protein","CID","prediction"])
			df["protein abbreviation"]=protein_abbrev_protwise
			df["metabolite"]=metabolite_protwise
			df["protein"]=";".join(uniprot).split(";")
			df["CID"]=cid_protwise
		df.set_index(["protein","CID"],inplace=True)
		for prot_cid in df.index:
			try:
				df.loc[prot_cid,"prediction"]=predictions.loc[prot_cid,"prediction"]
			except KeyError:
				continue
		df.to_csv(self.analyses+"prediction_verification_protwise"+self.feature+self.approach+appendix+".tsv",sep="\t",header=True,index=True)
		return
		
	
	
	#decides which verification funciton to use (OG or proteinwise)
	def verify_selected_interactions(self,predictions,df=None,organism=None,appendix=""):
		if self.proteinwise==False:
			self.verify_selected_interactions_og(predictions=predictions,df=df)
		else:
			self.verify_selected_interactions_protwise(predictions=predictions,df=df,organism=organism,appendix=appendix)
		return
	
	
	##########################
	#Trimming of predicitons
	#########################
	
	#select set of metabolites from profiles
	def select_metabolite_profiles(self,method=np.mean):
		try:
			metabolites=pd.read_csv(self.analyses+"metabolite_selection.txt",sep="\t")
		except FileNotFoundError:
			#load and extract metabolite profiles
			complete_set=pd.read_csv(self.databases+"complete_set_"+self.normalized+"normprofiles_"+method.__name__+".tsv",sep="\t").set_index("metabolite")
			complete_set=complete_set[~complete_set.index.duplicated()]
			non_metabolite_columns=list(filter(lambda col: "_metabolites" not in col,complete_set.columns))
			metabolite_profiles=complete_set.drop(non_metabolite_columns,axis=1)
			#plot
			#get separations between experiments
			experiments=list(map(lambda expfilename: expfilename[:-5].split("/")[-1],self.expfiles)) #experiments
			separators=list()
			num_fractions=0
			for exp in experiments[:-1]:
				num_fractions+=len(list(filter(lambda c: exp in c,metabolite_profiles.columns)))
				separators.append(num_fractions-0.5)
			#figure
			for i,cid in enumerate(metabolite_profiles.index):
				f,axn=plt.subplots(1,1,figsize=(10,10),num=i)
				axn.set_title(cid,weight="bold")
				pos_meta, =axn.plot(metabolite_profiles.loc[cid].tolist(),color="g")
				ymin=0
				ymax=max(metabolite_profiles.loc[cid].tolist())
				axn.vlines(separators,ymin=ymin,ymax=ymax,color="k",linestyles="dashed")
				plt.show()
				#save dialog
				save=messagebox.askyesno("Save?","Do you want to save the plot?")
				if save:
					metabolite_selection=open(self.analyses+"metabolite_selection.txt","a")
					metabolite_selection.write(str(cid)+"\n")
					metabolite_selection.close()
		return
	
	
	
	#exclude metabolites from metabolite_selection.txt interacting most
	def include_given_metabolites_only(self,predictions,method=np.mean):
		metabolites_df=pd.read_csv(self.analyses+"metabolite_selection.txt",header=None)
		metabolites=metabolites_df[metabolites_df.columns[0]].tolist()
		for prot_met in predictions.index:
			if prot_met[1] not in metabolites:
				predictions.drop(prot_met,axis=0,inplace=True)
		return(predictions)
		
	
	
	#exclude proteins interacting with majority of metabolites
	#threshold is threshold for trimming, e.g. trim proteins which interact with 75% or more of metabolites
	def trim_highest_ranked_proteins(self,predictions,classifier,method=np.mean,threshold=0.75):
		prot_ranks=pd.read_csv(self.analyses+self.protcol+"_ranks_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index(self.protcol)
		num_metas=len(predictions.index.get_level_values(1).unique())
		trim_prots=prot_ranks[prot_ranks[classifier]>=(threshold*num_metas)].index.tolist()
		for i in predictions.index:
			if i[0] in trim_prots:
				predictions.drop(i,axis=0,inplace=True)
		return(predictions)
	
	
	#exclude metabolites interacting with majority of metabolites
	#threshold is threshold for trimming, e.g. trim metabolites which interact with 75% or more of proteins
	def trim_highest_ranked_metabolites(self,predictions,classifier,method=np.mean,threshold=0.75):
		meta_ranks=pd.read_csv(self.analyses+"metabolite_ranks_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index("metabolite")
		num_prots=len(predictions.index.get_level_values(0).unique())
		trim_metas=meta_ranks[meta_ranks[classifier]>=(threshold*num_prots)].index.tolist()
		for i in predictions.index:
			if i[0] in trim_metas:
				predictions.drop(i,axis=0,inplace=True)
		return(predictions)
	
	
	#executes trimming of predictions and saves result, for metabolites and proteins only and together
	#trim_by_sel determines whether metabolites are trimmed by rank (=False) or by selection (=True)
	def trim_predictions_prots_metas_both(self,classifier_sel,predictions_sel,threshold,trim_by_sel=True,method=np.mean):
		print("interaction frequencies")
		for i in range(len(classifier_sel)):
			#trim by selected metabolites
			if trim_by_sel==True:
				try:
					predictions_trimmed_meta=pd.read_csv(self.analyses+"predictions_trimmed_less_interacting_metabolites_"+classifier_sel[i]+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
				except FileNotFoundError:
					predictions_trimmed_meta=self.include_given_metabolites_only(predictions=predictions_sel[i],method=method)
					predictions_trimmed_meta.to_csv(self.analyses+"predictions_trimmed_less_interacting_metabolites_"+classifier_sel[i]+".tsv",sep="\t",index=True,header=True)
			else:
				#trim highest ranked metabolites
				try:
					predictions_trimmed_meta=pd.read_csv(self.analyses+"predictions_trimmed_less_interacting_metas_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
				except FileNotFoundError:
					predictions_trimmed_meta=self.trim_highest_ranked_metabolites(predictions=predictions_sel[i],classifier=classifier_sel[i],threshold=threshold,method=method)
					predictions_trimmed_meta.to_csv(self.analyses+"predictions_trimmed_less_interacting_metas_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t",index=True,header=True)
			print(classifier_sel[i])
			print("metabolites trimmed: "+str(len(predictions_trimmed_meta[predictions_trimmed_meta["prediction"]==True])/len(predictions_trimmed_meta)))
			#trim highest ranked proteins
			try:
				predictions_trimmed_prot=pd.read_csv(self.analyses+"predictions_trimmed_less_interacting_prots_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
			except FileNotFoundError:
				predictions_trimmed_prot=self.trim_highest_ranked_proteins(predictions=predictions_sel[i],classifier=classifier_sel[i],threshold=threshold,method=method)
				predictions_trimmed_prot.to_csv(self.analyses+"predictions_trimmed_less_interacting_prots_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t",index=True,header=True)
			print(len(predictions_trimmed_prot[predictions_trimmed_prot["prediction"]==True])/len(predictions_trimmed_prot))
			#trim highest ranked proteins and metabolites
			try:
				if trim_by_sel==True:
					predictions_trimmed_prot_meta=pd.read_csv(self.analyses+"predictions_trimmed_less_interacting_prots_and_metas_selection_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
				else:
					predictions_trimmed_prot_meta=pd.read_csv(self.analyses+"predictions_trimmed_less_interacting_prots_and_metas_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
			except FileNotFoundError:
				predictions_trimmed_prot_meta=self.trim_highest_ranked_proteins(predictions=predictions_trimmed_meta,classifier=classifier_sel[i],threshold=threshold,method=method)
				if trim_by_sel==True:
					predictions_trimmed_prot_meta.to_csv(self.analyses+"predictions_trimmed_less_interacting_prots_and_metas_selection_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t",index=True,header=True)
				else:
					predictions_trimmed_prot_meta.to_csv(self.analyses+"predictions_trimmed_less_interacting_prots_and_metas_"+classifier_sel[i]+"_thresh_"+str(threshold)+".tsv",sep="\t",index=True,header=True)
			print(len(predictions_trimmed_prot_meta[predictions_trimmed_prot_meta["prediction"]==True])/len(predictions_trimmed_prot_meta))
		return
	
	
	#trim predictions on given prediction probability threshold
	def trim_predictions_on_predprobs(self,predictions,classifier,threshold,method=np.mean):
		prediction_probabilities=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		predprobs_trimmed=prediction_probabilities[prediction_probabilities[classifier]>=threshold]
		intersect_predprobs_preds=set(set(predprobs_trimmed.index.tolist()) & set(predictions.index.tolist()))
		predictions_trimmed=predictions.copy()
		predictions_trimmed["prediction"]=False
		predictions_trimmed.loc[intersect_predprobs_preds,"prediction"]=True
		return(predictions_trimmed)
	
	
	
	#trim predictions on proteins interacting with more than met_thresh metabolites and probability higher prob_thresh
	def trim_predictions_on_predprobs_and_protrank(self,predictions,met_thresholds,prob_thresholds,classifier,est_int_prob,method=np.mean):
		try:
			df=pd.read_csv(self.analyses+"estimation_interaction_probability_"+classifier+".tsv",sep="\t")
		except FileNotFoundError:
			df=pd.DataFrame(columns=["metabolite threshold","prediction probability threshold","interaction probability"])
			df.set_index(["metabolite threshold","prediction probability threshold"],inplace=True)
			for met_thresh in met_thresholds:
				predictions_trimmed_prot=pd.read_csv(self.analyses+"predictions_trimmed_less_interacting_prots_"+classifier+"_thresh_"+str(met_thresh)+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
				for prob_thresh in prob_thresholds:
					predictions_trimmed_prot_probs=self.trim_predictions_on_predprobs(predictions=predictions_trimmed_prot,classifier=classifier,threshold=prob_thresh,method=method)
					df.loc[(met_thresh,prob_thresh),"interaction probability"]=len(predictions_trimmed_prot_probs[predictions_trimmed_prot_probs["prediction"]==True])/len(predictions_trimmed_prot_probs)
			df.to_csv(self.analyses+"estimation_interaction_probability_"+classifier+".tsv",sep="\t")
			df.reset_index(inplace=True)
		#prepare data for plotting
		#surface showing estimated interaction probability
		x,y=np.meshgrid(np.linspace(df["metabolite threshold"].min(),df["metabolite threshold"].max(),5),np.linspace(df["prediction probability threshold"].min(),df["prediction probability threshold"].max(),5))
		z=x*0+y*0+est_int_prob
		#3D plot
		f=plt.figure()
		ax=f.gca(projection="3d")
		ax.plot_surface(x,y,z,color="tab:red",alpha=0.3)
		ax.plot_trisurf(df["metabolite threshold"],df["prediction probability threshold"],df["interaction probability"],color="tab:blue")
		ax.set_xlim(df["metabolite threshold"].min(),df["metabolite threshold"].max())
		ax.set_ylim(df["prediction probability threshold"].min(),df["prediction probability threshold"].max())
		ax.set_zlim(0,0.5)
		ax.set_xlabel("metabolite interactors threshold")
		ax.set_ylabel("prediction probability threshold")
		ax.set_zlabel("interaction probability")
		plt.show()
		f.savefig("../overall_analysis/figures/interactome_estimation/estimation_interaction_probability_"+classifier+".png",sep="\t")
		f.savefig("../overall_analysis/figures/interactome_estimation/estimation_interaction_probability_"+classifier+".svg",sep="\t")
		return(df)
		
	
	
	#trim predictions by prediction probability
	#thresholds is list of prediction probability thresholds
	def interaction_probability_from_prediction_prob(self,classifier,thresholds,method=np.mean):
		prediction_probability=pd.read_csv(self.analyses+"prediction_probabilities_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		predprobs=prediction_probability[classifier].tolist()
		interaction_probabilities=list()
		for thresh in thresholds:
			predprobs_trimmed=list(filter(lambda x: x>=thresh,predprobs))
			interaction_probabilities.append(len(predprobs_trimmed)/len(prediction_probability))
		#plot
		f=plt.figure()
		plt.plot(thresholds,interaction_probabilities)
		plt.show()
		return
	
	
	
	##############################
	#figure: comparison to random
	##############################
	
	#compares learning curve and ROC curve and degree distributions to random training
	def comparison_to_random_training(self,classifier,appendix="",method=np.min,colours=["tab:orange","tab:green","tab:red","tab:blue","tab:purple","tab:brown","k"],legend=True):
		#load data
		try:
			lc_df=pd.read_csv(self.analyses+"learning_curves_trained_and_random_csts"+self.feature+self.approach+appendix+".tsv",sep="\t")
			roc_df=pd.read_csv(self.analyses+"ROC_curve_trained_and_random_csts"+self.feature+self.approach+appendix+".tsv",sep="\t")
			ylim_lc=(np.min([np.round(float(x),3) for x in ";".join(lc_df["test_score_min"].tolist()).split(";")]),1)
		except FileNotFoundError:
			print("please run MLanalysis.create_learning_ROC_curves_comparison(appendix="+appendix+")")
			return
		try:
			#degree distributions
			'''
			##true network
			ddp_df_true=pd.read_csv(self.analyses+"degree_distribution_prot"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
			ddm_df_true=pd.read_csv(self.analyses+"degree_distribution_met"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
			eq_prot_true=open(self.analyses+"powerlaw_func_prot"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
			eq_met_true=open(self.analyses+"powerlaw_func_met"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
			##random ts and cs
			ddp_df_ts=pd.read_csv(self.analyses+"degree_distribution_protts"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
			ddm_df_ts=pd.read_csv(self.analyses+"degree_distribution_metts"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
			eq_prot_ts=open(self.analyses+"powerlaw_func_protts"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
			eq_met_ts=open(self.analyses+"powerlaw_func_metts"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
			ddp_df_cs=pd.read_csv(self.analyses+"degree_distribution_protcs"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
			ddm_df_cs=pd.read_csv(self.analyses+"degree_distribution_metcs"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".tsv",sep="\t")
			eq_prot_cs=open(self.analyses+"powerlaw_func_protcs"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
			eq_met_cs=open(self.analyses+"powerlaw_func_metcs"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read()
			#concatenate them to list
			ddp_dfs=[ddp_df_true,ddp_df_ts,ddp_df_cs]
			ddm_dfs=[ddm_df_true,ddm_df_ts,ddm_df_cs]
			eq_prots=[eq_prot_true,eq_prot_ts,eq_prot_cs]
			eq_mets=[eq_met_true,eq_met_ts,eq_met_cs]
			'''
			prot_degrees_norm=open(self.analyses+"protein_degrees"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
			met_degrees_norm=open(self.analyses+"metabolite_degrees"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
			prot_degrees_ts=open(self.analyses+"protein_degreests"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
			met_degrees_ts=open(self.analyses+"metabolite_degreests"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
			prot_degrees_cs=open(self.analyses+"protein_degreescs"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
			met_degrees_cs=open(self.analyses+"metabolite_degreescs"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".txt","r").read().split("\n")
			prot_degrees=[[int(x) for x in prot_degrees_norm],[int(y) for y in prot_degrees_ts],[int(z) for z in prot_degrees_cs]]
			met_degrees=[[int(x) for x in met_degrees_norm],[int(y) for y in met_degrees_ts],[int(z) for z in met_degrees_cs]]
			upper_ylim_mets=np.max([int(x) for x in met_degrees_norm+met_degrees_ts+met_degrees_cs])
			upper_ylim_prots=np.max([int(x) for x in prot_degrees_norm+prot_degrees_ts+prot_degrees_cs])
			bins=np.linspace(0,upper_ylim_mets,100)
		except FileNotFoundError:
			print("create data for degree distributions with NetworkAnalysis.degree_distribution()")
			return
		f,axn=plt.subplots(figsize=(10,10),nrows=len(lc_df),ncols=3)
		#colours=["tab:orange","tab:green","tab:red","tab:blue","tab:purple","tab:brown","k"]
		for r in range(len(lc_df)):
			#learning curve
			axn[r][0].fill_between(lc_df.loc[lc_df.index[r],"train_size"].split(";"),list(map(lambda x: eval(x),lc_df.loc[lc_df.index[r],"train_score_min"].split(";"))),list(map(lambda x: eval(x),lc_df.loc[lc_df.index[r],"train_score_max"].split(";"))),alpha=0.1,color="k")
			axn[r][0].fill_between(lc_df.loc[lc_df.index[r],"train_size"].split(";"),list(map(lambda x: eval(x),lc_df.loc[lc_df.index[r],"test_score_min"].split(";"))),list(map(lambda x: eval(x),lc_df.loc[lc_df.index[r],"test_score_max"].split(";"))),alpha=0.1,color=colours[r])
			axn[r][0].plot(lc_df.loc[lc_df.index[r],"train_size"].split(";"),list(map(lambda x: eval(x),lc_df.loc[lc_df.index[r],"train_score_mean"].split(";"))),color="k",label="training score")
			axn[r][0].plot(lc_df.loc[lc_df.index[r],"train_size"].split(";"),list(map(lambda x: eval(x),lc_df.loc[lc_df.index[r],"test_score_mean"].split(";"))),color=colours[r],label="test score")
			axn[r][0].set_ylabel("accuracy")
			axn[r][0].set_ylim(ylim_lc)
			if r<len(lc_df)-1:
				axn[r][0].set_xticklabels([])
			#ROC curve
			axn[r][1].plot([0,1],[0,1],color="k",linestyle="--")
			axn[r][1].plot([0,eval(roc_df.loc[roc_df.index[r],"false_positive_rate"].split(";")[0])],[0,eval(roc_df.loc[roc_df.index[r],"true_positive_rate"].split(";")[0])],color=colours[r])
			axn[r][1].plot(list(map(lambda x: eval(x),roc_df.loc[roc_df.index[r],"false_positive_rate"].split(";"))),list(map(lambda x: eval(x),roc_df.loc[roc_df.index[r],"true_positive_rate"].split(";"))),color=colours[r],label="(AUC = %0.2f)" % roc_df.loc[roc_df.index[r],"AUC"])
			axn[r][1].set_ylabel("true positive rate")
			axn[r][1].legend(loc="lower right")
			if r<len(lc_df)-1:
				axn[r][1].set_xticklabels([])
			#Degree distribution
			'''
			axn[r][2].plot(ddp_dfs[r]["degree"],ddp_dfs[r]["cumulative_probability"].tolist(),"o",color=colours[r])
			axn[r][2].plot(ddm_dfs[r]["degree"],ddm_dfs[r]["cumulative_probability"].tolist(),".",color=colours[r],alpha=0.1)
			axn[r][2].plot(ddp_dfs[r]["degree"],ddp_dfs[r]["fit"],linestyle="-",color=colours[r],label=eq_prots[r])
			axn[r][2].plot(ddm_dfs[r]["degree"],ddm_dfs[r]["fit"],linestyle="-",color=colours[r],alpha=0.1,label=eq_mets[r])
			axn[r][2].set_ylabel("P(k > K)")
			axn[r][2].set_yscale("log")
			axn[r][2].set_xscale("log")
			if legend:
				axn[r][2].legend(prop={"size":8})
			'''
			axn[r][2].hist(prot_degrees[r],color=colours[r],alpha=0.5)
			axn[r][2].set_ylabel("abs. frequency",color=colours[r])
			axn[r][2].set_ylim(0,upper_ylim_mets)
			axn[r][2].tick_params(axis="y",color=colours[r],labelcolor=colours[r])
			axn[r][2].tick_params(axis="x",color=colours[r])
			ax1=axn[r][2].twinx()
			ax2=ax1.twiny()
			ax2.hist(met_degrees[r],color="tab:blue",alpha=0.5)
			ax1.set_ylim(0,upper_ylim_prots)
			ax1.tick_params(axis="y",color="tab:blue",labelcolor="tab:blue")
			ax2.set_xlim(0,upper_ylim_mets)
			ax2.tick_params(axis="x",color="tab:blue",labelcolor="tab:blue")
			ax1.set_ylabel("abs. frequency",color="tab:blue")
			if r<len(lc_df)-1:
				axn[r][2].set_xticklabels([])
			if r!=0:
				ax2.set_xticklabels([],color="tab:blue")
		axn[r][0].set_xlabel("size of training set")
		axn[r][1].set_xlabel("false positive rate")
		axn[r][2].set_xlabel("Degree of nodes K")
		f.tight_layout(rect=(0,0.1,0.9,0.9))
		#plt.subplots_adjust(bottom=0.3,wspace=0.9)
		plt.show()
		f.savefig("../overall_analysis/figures/comparison_random/comparison_to_random_training_LC_ROC_DD"+self.feature+self.approach+appendix+".png")
		f.savefig("../overall_analysis/figures/comparison_random/comparison_to_random_training_LC_ROC_DD"+self.feature+self.approach+appendix+".svg")
		return
	
	
	
	#assortativity of a bipartite network relative to node number in other type
	#g=graph
	def my_assortativity(self,g):
		#define type (assign to group of bipartiteness)
		g.vs["type"]=list(map(lambda x: type(x)==int,g.vs["name"])) #1 for metabolites, 0 for proteins ###adapt
		max_possible_degree_prots=np.sum(g.vs["type"])
		max_possible_degree_metas=len(g.vs)-np.sum(g.vs["type"])
		avg_neighbour_degrees=list()
		degrees=list()
		for v in g.vs:
			if v["type"]==1: #metabolites
				avg_neighbour_degrees.append(np.mean(g.vs[g.neighbors(v)].degree())/max_possible_degree_prots)
				degrees.append(v.degree()/max_possible_degree_metas)
			else: #proteins
				avg_neighbour_degrees.append(np.mean(g.vs[g.neighbors(v)].degree())/max_possible_degree_metas)
				degrees.append(v.degree()/max_possible_degree_prots)
		avg_neighbour_degrees=[0 if np.isnan(x) else x for x in avg_neighbour_degrees]
		degrees=[0 if np.isnan(x) else x for x in degrees]
		assortativity=np.corrcoef(avg_neighbour_degrees,degrees)[0][1]
		return(assortativity)
	
	
	
	#compares accuracy, assortativity between random networks as bar plots
	#and edge and node weights correlations as heatmaps
	def comparison_to_random_network(self,classifier,int_probs,method=np.mean,reps=1000):
		#load network and test set
		#ec_names={"":"Other Proteins", "1":"Oxidorecutases","2":"Transferases","3":"Hydrolases","4":"Lyases","5":"Isomerases","6":"Ligases","7":"Translocases"}
		predictions=pd.read_csv(self.analyses+"predictions_annotated_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t",dtype={"EC":str,"metabolite":int}).set_index([self.protcol,"metabolite"]).replace(np.nan,"")
		#predictions=predictions.replace({"EC":ec_names})
		test_set=pd.read_csv(self.analyses+"training_set"+self.feature+self.approach+"_"+method.__name__+".tsv",sep="\t").set_index([self.protcol,"metabolite"])
		performance_table=pd.read_csv(self.analyses+"analysis_summary"+self.feature+self.approach+".tsv",sep="\t").set_index(["methods","classifiers"])
		df=pd.DataFrame(columns=["accuracy","accuracy_sd","assortativity","assortativity_sd","my_assortativity","my_assortativity_sd"],index=int_probs)
		#create graph
		edges=[tuple(x) for x in predictions[predictions.loc[:,"prediction"]==True].index]
		g=Graph.TupleList(edges,directed=False)
		#trim non-interacting molecules
		df.loc["trimmed network","accuracy"]=performance_table.loc[(method.__name__,classifier),"accuracy"]
		df.loc["trimmed network","accuracy_sd"]=performance_table.loc[(method.__name__,classifier),"accuracy_sd"]
		df.loc["trimmed network","assortativity"]=g.assortativity_degree(directed=False)
		df.loc["trimmed network","my_assortativity"]=self.my_assortativity(g)
		#include vertices (OGs and metabolites) which have no interacting partners
		vertices_g=set(g.vs["name"])
		vertices_df=set(predictions.index.get_level_values(0).tolist()+predictions.index.get_level_values(1).tolist())
		remaining_vertices=vertices_df-vertices_g
		g.add_vertices(list(remaining_vertices))
		#normal
		df.loc["normal network","accuracy"]=performance_table.loc[(method.__name__,classifier),"accuracy"]
		df.loc["normal network","accuracy_sd"]=performance_table.loc[(method.__name__,classifier),"accuracy_sd"]
		df.loc["normal network","assortativity"]=g.assortativity_degree(directed=False)
		df.loc["normal network","my_assortativity"]=self.my_assortativity(g)
		#drop nan
		test_set.drop(set(test_set.index.tolist())-set(predictions.index.tolist()),inplace=True)
		accuracy_training_normal=metrics.accuracy_score(test_set["interaction"].tolist(),predictions.loc[test_set.index,"prediction"].tolist()) #accuracy on training data!
		#random
		for p in int_probs:
			predictions_rand=predictions.copy()
			accuracies_p=list()
			assortativities_p=list()
			my_assortativities_p=list()
			for r in range(reps):
				predictions_rand["prediction"]=choices([True,False],weights=[p,1-p],k=len(predictions_rand))
				accuracies_p.append(metrics.accuracy_score(test_set["interaction"].tolist(),predictions_rand.loc[test_set.index,"prediction"].tolist()))
				#create graph
				edges=[tuple(x) for x in predictions_rand[predictions_rand.loc[:,"prediction"]==True].index]
				g_rand=Graph.TupleList(edges,directed=False)
				#include vertices (OGs and metabolites) which have no interacting partners
				vertices_g_rand=set(g_rand.vs["name"])
				vertices_preds_rand=set(predictions_rand.index.get_level_values(0).tolist()+predictions_rand.index.get_level_values(1).tolist())
				remaining_vertices_rand=vertices_preds_rand-vertices_g_rand
				g_rand.add_vertices(list(remaining_vertices_rand))
				assortativities_p.append(g_rand.assortativity_degree(directed=False))
				my_assortativities_p.append(self.my_assortativity(g_rand))
			df.loc[p,"accuracy"]=np.mean(accuracies_p)
			df.loc[p,"accuracy_sd"]=np.std(accuracies_p)
			df.loc[p,"assortativity"]=np.mean(assortativities_p)
			df.loc[p,"assortativity_sd"]=np.std(assortativities_p)
			#use own assortativity relative to degree of other node type
			df.loc[p,"my_assortativity"]=np.mean(my_assortativities_p)
			df.loc[p,"my_assortativity_sd"]=np.std(my_assortativities_p)
		#3. compare edge correlations
		#...
		#4. compare vertex correlations
		df.to_csv(self.analyses+"comparison_to_random_intprobs_accuracy_assortativity_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t",index=True,header=True)
		return
		
	
	
	#plot histograms of random accuracies and assortativities
	#first column: accuracies, second column: assortativities
	def hist_acc_ass(self,method,classifier):
		df1=pd.read_csv(self.analyses+"comparison_to_random_intprobs_accuracy_assortativity_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t")
		df2=pd.read_csv(self.analyses+"comparison_to_random_ts_cs_assortativity_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t")
		df=df1.append(df2).set_index("Unnamed: 0")
		yaxis=["accuracy","my_assortativity"]
		yaxis_dict={"accuracy":"accuracy","my_assortativity":"rel. assortativity"}
		label_dict={"complete":"rand. profiles","training":"rand. train. set","normal network":"pred. network"}
		f,axn=plt.subplots(nrows=1,ncols=2)
		#do not plot accuracy for training
		for i in range(2):
			labels=list(map(lambda x: label_dict[x] if x in label_dict else x,df.index))
			colors=["tab:blue"]*len(df)
			colors[labels.index("pred. network")]="tab:red"
			xdata=labels
			ydata=df.loc[:,yaxis[i]].tolist()
			errors=df.loc[:,yaxis[i]+"_sd"].tolist()
			if yaxis[i]=="accuracy":
				ind=xdata.index(label_dict["training"])
				del(xdata[ind])
				del(ydata[ind])
				del(errors[ind])
				del(colors[ind])
			axn[i].bar(xdata,ydata,yerr=errors,color=colors,tick_label=labels)
			axn[i].set_ylabel(yaxis_dict[yaxis[i]])
			axn[i].set_xlabel("P(interaction)")
			axn[i].set_xticklabels(labels,rotation=65,ha="right")
		f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/figures/comparison_random/accuracy_assortativity_hist"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".png")
		f.savefig("../overall_analysis/figures/comparison_random/accuracy_assortativity_hist"+self.feature+self.approach+"_"+method.__name__+"_"+classifier+".svg")
		return
	
	
	#venn diagram over random predictions and made predictions
	def venn_diagram_randoms(self,predictions_exp,predictions_training,predictions_complete,colors=None,labels=["pred. network","rand. train. set","rand. profiles"],appendix="",fontsize=14):
		predictions=[predictions_exp,predictions_training,predictions_complete]
		ints=list()
		for preds in predictions:
			ints.append(set(preds[preds["prediction"]==True].index.tolist()))
		data=venn.get_labels(ints)
		if colors is None:
			f,axn=venn.venn3(data,names=labels,fontsize=fontsize)
		else:
			f,axn=venn.venn3(data,names=labels,fontsize=fontsize,colors=colors)
		plt.show()
		f.savefig("../overall_analysis/figures/comparison_random/venn_diagram_random"+appendix+".png")
		f.savefig("../overall_analysis/figures/comparison_random/venn_diagram_random"+appendix+".svg")
		return
	
	
	
	
	
	#########################
	#figure feature engineering 
	#########################
	
	#creates figure feature engineering panel example (1st row)
	#mode={"same","full"}
	def fig_feature_engineering_example(self,mode="full"):
		#init profiles
		prot_profile=[0.02,0.2,0.5,0.9,1,0.9,0.5,0.2,0.02,0.02,0.02,0.02,0.02,0.02,0.02]
		int_met_profile=[0.02,0.2,0.5,0.9,1,0.9,0.5,0.2,0.02,0.02,0,0.02,0.02,0.02,0.02]
		nonint_met_profile=[0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.2,0.5,0.9,1,0.9,0.5,0.2,0.02]
		crosscorr_int=np.correlate(prot_profile,int_met_profile,mode=mode)
		crosscorr_nonint=np.correlate(prot_profile,nonint_met_profile,mode=mode)
		if mode=="same":
			fractional_shift=[-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7]
		else:
			fractional_shift=[-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
		xlim_max=15.5
		ylim_max=1.05
		#figure protein profile
		fig_prot=plt.figure()
		ax=SubplotZero(fig_prot,111)
		fig_prot.add_subplot(ax)
		for direction in ["xzero", "yzero"]:
			ax.axis[direction].set_axisline_style("-|>")
			ax.axis[direction].set_visible(True)
		for direction in ["left", "right", "bottom", "top"]:
			ax.axis[direction].set_visible(False)
		ax.plot(prot_profile,color="tab:green",linewidth=5)
		ax.set_xlim(0,xlim_max)
		ax.set_ylim(0,ylim_max)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_yticks([])
		#ax.set_xlabel("frac.")
		#ax.set_ylabel("norm. int.")
		#plt.show()
		fig_prot.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row1_protein_profile.svg")
		plt.close()
		
		#figure interacting metabolite profile
		fig_int_met=plt.figure()
		ax=SubplotZero(fig_int_met,111)
		fig_int_met.add_subplot(ax)
		for direction in ["xzero", "yzero"]:
			ax.axis[direction].set_axisline_style("-|>")
			ax.axis[direction].set_visible(True)
		for direction in ["left", "right", "bottom", "top"]:
			ax.axis[direction].set_visible(False)
		ax.plot(int_met_profile,color="tab:blue",linewidth=5)
		ax.set_xlim(0,xlim_max)
		ax.set_ylim(0,ylim_max)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_yticks([])
		#ax.set_xlabel("frac.")
		#ax.set_ylabel("norm. int.")
		#plt.show()
		fig_int_met.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row1_int_metabolite_profile.svg")
		plt.close()
		
		#figure not interacting metabolite profile
		fig_nonint_met=plt.figure()
		ax=SubplotZero(fig_nonint_met,111)
		fig_nonint_met.add_subplot(ax)
		for direction in ["xzero", "yzero"]:
			ax.axis[direction].set_axisline_style("-|>")
			ax.axis[direction].set_visible(True)
		for direction in ["left", "right", "bottom", "top"]:
			ax.axis[direction].set_visible(False)
		ax.plot(nonint_met_profile,color="tab:blue",linewidth=5)
		ax.set_xlim(0,xlim_max)
		ax.set_ylim(0,ylim_max)
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_yticks([])
		#ax.set_xlabel("frac.")
		#ax.set_ylabel("norm. int.")
		#plt.show()
		fig_nonint_met.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row1_nonint_metabolite_profile.svg")
		plt.close()
		
		#figure crosscorr_interaction
		fig_corr_int=plt.figure()
		ax=SubplotZero(fig_corr_int,111)
		fig_corr_int.add_subplot(ax)
		for direction in ["xzero", "yzero"]:
			ax.axis[direction].set_axisline_style("-|>")
			ax.axis[direction].set_visible(True)
		for direction in ["left", "right", "bottom", "top"]:
			ax.axis[direction].set_visible(False)
		ax.plot(fractional_shift,crosscorr_int,color="tab:red",linewidth=5)
		ax.set_xlim(np.min(fractional_shift),np.max(fractional_shift))
		ax.set_ylim(0,1.05*np.max(crosscorr_int))
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_yticks([])
		#ax.set_xlabel("frac. shift")
		#ax.set_ylabel("corr")
		#plt.show()
		fig_corr_int.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row1_crosscorr_int_profile_"+mode+".svg")
		plt.close()
		
		#figure crosscorr no interaction
		fig_corr_nonint=plt.figure()
		ax=SubplotZero(fig_corr_nonint,111)
		fig_corr_nonint.add_subplot(ax)
		for direction in ["xzero", "yzero"]:
			ax.axis[direction].set_axisline_style("-|>")
			ax.axis[direction].set_visible(True)
		for direction in ["left", "right", "bottom", "top"]:
			ax.axis[direction].set_visible(False)
		ax.plot(fractional_shift,crosscorr_nonint,color="tab:red",linewidth=5)
		ax.set_xlim(np.min(fractional_shift),np.max(fractional_shift))
		ax.set_ylim(0,1.05*np.max(crosscorr_int))
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.set_xticks([])
		ax.set_yticks([])
		#ax.set_xlabel("frac. shift")
		#ax.set_ylabel("corr")
		#plt.show()
		fig_corr_nonint.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row1_crosscorr_nonint_profile_"+mode+".svg")
		plt.close()
		return
	
	
	
	#creates figure feature engineering panel weights (2nd row)
	def fig_feature_engineering_weights(self):
		print("create venn diagram between different features with NetworkAnalysis.py")
		print("create feature weighting for selected with MLanalysis.plot_feature_attributes()")
		return
	
	
	#creates figure feature engineering panel performance measures (3rd row)
	#feature_oi = feature of interest, the one which shall be highlighted
	def fig_feature_engineering_performance(self,features,method,classifier,feature_oi="_same_crosscorr_not_normalized",appendix=""):
		#collect tables with performance measures
		performance_measure_dict={"precision":"precision","recall":"recall","f_1 score":"F$_1$ score","accuracy":"accuracy","sensitivity":"sensitivity","specificity":"specificity"}
		feature_dict={"_same_crosscorr":"not normalized","_same_crosscorr_not_normalized":"profile length","_same_crosscorr_binsize10_not_normalized":"binsize 10","_same_crosscorr_binsize4_not_normalized":"binsize 4","_same_crosscorr_binsize2_not_normalized":"binsize 2"}
		for feature in features:
			performance_table=pd.read_csv(self.analyses+"analysis_summary"+feature+self.approach+".tsv",sep="\t").set_index(["methods","classifiers"])
			if feature==features[0]:
				pm_df=pd.DataFrame(columns=performance_table.columns)
			pm_df.loc[feature]=performance_table.loc[(method.__name__,classifier)]
		#parmeters for figure
		num_cols=int(len(pm_df.columns)/2)
		##set red colour for selected feature of interest, blue for others
		feature_oi_index=features.index(feature_oi)
		colors=["tab:blue"]*len(features)
		colors[feature_oi_index]="tab:red"
		#plot
		f,axn=plt.subplots(nrows=1,ncols=num_cols,sharey=False)
		for c in range(num_cols):
			xdata=pm_df.index
			ydata=pm_df.loc[:,pm_df.columns[c]]
			upper_errors=pm_df.loc[:,pm_df.columns[num_cols+c]]
			lower_errors=list(map(lambda x: ydata.iloc[x] if upper_errors.iloc[x]>ydata.iloc[x] else upper_errors.iloc[x],range(len(upper_errors.tolist()))))
			errors=[lower_errors,upper_errors.tolist()]
			axn[c].bar(xdata,ydata,yerr=errors,color=colors,tick_label=list(map(lambda x: feature_dict[x],features)))
			axn[c].set_ylabel(performance_measure_dict[pm_df.columns[c]])
			#axn[c].set_xticks(np.arange(len(features)),list(map(lambda x: feature_dict[x],features)))
			axn[c].set_xticklabels(list(map(lambda x: feature_dict[x],features)),rotation=65,ha="right")
		plt.subplots_adjust(bottom=0.3,wspace=0.9)
		#f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row3_"+classifier+"_"+method.__name__+appendix+".png")
		f.savefig("../overall_analysis/figures/feature_engineering/feature_engineering_row3_"+classifier+"_"+method.__name__+appendix+".svg")
		return
		
	
	
	#plots figure heading letters
	def plot_letters(self):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		f,axn=plt.subplots(nrows=4,ncols=6)
		for a,ax in enumerate(axn.flat):
			ax.axis("off")
			ax.set_title(abc[a],weight="bold")
		plt.show()
		f.savefig("../overall_analysis/figures/letters.svg")
		return
	
	
	
	#plots figure heading headers
	def plot_headers(self,headers,plot_with_letters=False,weight="bold"):
		abc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		f,axn=plt.subplots(nrows=len(headers),ncols=1)
		for a,header in enumerate(headers):
			axn[a].axis("off")
			if plot_with_letters:
				axn[a].set_title(abc[a]+" "+header,weight=weight)
			else:
				axn[a].set_title(header,weight=weight)
		plt.show()
		f.savefig("../overall_analysis/figures/headers.svg")
		return
	
	
	###########################################
	# Figure: Enrichment
	###########################################
	
	
	#heatmap brite class vs metabolite class
	def heatmap_bc_mc(self,method,significance_level=0.05,trim="",square=False):
		brite_df=pd.read_csv(self.analyses+"interactions_briteclass_metabolite_class"+trim+"_"+method.__name__+self.feature+self.approach+".tsv",sep="\t").set_index("BRITE class")
		exclude_brite_classes=[np.nan,"br:sce00001","br:ath00001"]
		brite_dict=self.get_brite_dict()
		brite_df.drop(exclude_brite_classes,axis=0,inplace=True)
		brite_df["BRITE class name"]=list(map(lambda x: brite_dict[x[6:]],brite_df.index.tolist()))
		#rearrange df for heatmap
		heatmap=pd.DataFrame(index=brite_df["BRITE class name"].unique(),columns=brite_df["metabolite class"].unique())
		brite_df=brite_df.reset_index().set_index(["BRITE class name","metabolite class"])
		for bc in heatmap.index:
			for mc in heatmap.columns:
				heatmap.loc[bc,mc]=brite_df.loc[(bc,mc),"number of interactions"]
		#filter out classes with zero interactions
		mc_zeros=heatmap.columns[heatmap.sum(axis=0)==0]
		bc_zeros=heatmap.index[heatmap.sum(axis=1)==0]
		heatmap.drop(mc_zeros,axis=1,inplace=True)
		heatmap.drop(bc_zeros,axis=0,inplace=True)
		#significance
		significance_df=heatmap.copy()
		#number of interactions and non-interactions per BRITE class for significance
		sum_ints_df=pd.DataFrame(index=brite_df.index.get_level_values(0).unique(),columns=["num_ints","num_nonints"])
		sum_ints_df["num_ints"]=list(map(lambda x: brite_df.loc[x,"number of interactions"].sum(),sum_ints_df.index))
		sum_ints_df["num_nonints"]=list(map(lambda x: brite_df.loc[x,"number of non-interactions"].sum(),sum_ints_df.index))
		for bc in significance_df.index:
			for mc in significance_df.columns:
				ints=brite_df.loc[(bc,mc),"number of interactions"]
				nonints=brite_df.loc[(bc,mc),"number of non-interactions"]
				contingency_table=[[ints,nonints],[sum_ints_df.loc[bc,"num_ints"]-ints,sum_ints_df.loc[bc,"num_nonints"]-nonints]]
				oddsratio,pvalue=fisher_exact(contingency_table,alternative="greater")
				if pvalue<significance_level*1/len(sum_ints_df):
					significance_df.loc[bc,mc]="*"
				else:
					significance_df.loc[bc,mc]=""
		#plot heatmap
		f=plt.figure()#figsize=(25,21))
		sns.heatmap(heatmap,annot=significance_df,fmt="",square=square,robust=True,xticklabels=True,yticklabels=True,cmap="hot_r") #vmin=0,vmax=brite_df["number of interactions"].max(),linewidths=0.5
		f.tight_layout()
		plt.show()
		f.savefig("../overall_analysis/figures/enrichment/heatmap_BC_MC"+trim+self.approach+self.feature+"_"+method.__name__+"_"+classifier+".png")
		f.savefig("../overall_analysis/figures/enrichment/heatmap_BC_MC"+trim+self.approach+self.feature+"_"+method.__name__+"_"+classifier+".svg")
		return
	
	
	#plots heatmap with number of interactions between different pathways in intersection
	def heatmap_pathways(self,method,classifier,significance_level=0.05,trim="",square=False):
		try:
			pathway_df=pd.read_csv(self.analyses+"interactions_kegg_pathways"+trim+"_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t")
		except FileNotFoundError:
			pathway_df=self.make_pathway_dataframe(method=method,classifier=classifier)
			pathway_df=pd.read_csv(self.analyses+"interactions_kegg_pathways"+trim+"_"+method.__name__+self.feature+self.approach+"_"+classifier+".tsv",sep="\t")
		pathway_df.set_index(["protein pathway","metabolite pathway"],inplace=True)
		pathway_dict=self.get_pathway_dict()
		exclude_pathways={"01100"}
		exclude_pathway_names=set(map(lambda x: pathway_dict[x],exclude_pathways))
		#initialize dataframe for heatmap
		heatmap=pd.DataFrame(index=set(pathway_df.index.get_level_values(0))-exclude_pathway_names,columns=set(pathway_df.index.get_level_values(1))-exclude_pathway_names)
		for prot_pathway in heatmap.index:
			for meta_pathway in heatmap.columns:
				heatmap.loc[prot_pathway,meta_pathway]=pathway_df.loc[(prot_pathway,meta_pathway),"number of interactions"]
		#filter out pathways with zero interactions
		metapaths_zeros=heatmap.columns[heatmap.sum(axis=0)==0]
		protpaths_zeros=heatmap.index[heatmap.sum(axis=1)==0]
		heatmap.drop(metapaths_zeros,axis=1,inplace=True)
		heatmap.drop(protpaths_zeros,axis=0,inplace=True)
		#significance
		significance_df=heatmap.copy()
		#number of interactions and non-interactions per BRITE class for significance
		sum_ints_df=pd.DataFrame(index=pathway_df.index.get_level_values(0).unique(),columns=["num_ints","num_nonints"])
		sum_ints_df["num_ints"]=list(map(lambda x: pathway_df.loc[x,"number of interactions"].sum(),sum_ints_df.index))
		sum_ints_df["num_nonints"]=list(map(lambda x: pathway_df.loc[x,"number of non-interactions"].sum(),sum_ints_df.index))
		for protpath in significance_df.index:
			for metapath in significance_df.columns:
				ints=pathway_df.loc[(protpath,metapath),"number of interactions"]
				nonints=pathway_df.loc[(protpath,metapath),"number of non-interactions"]
				contingency_table=[[ints,nonints],[sum_ints_df.loc[protpath,"num_ints"]-ints,sum_ints_df.loc[protpath,"num_nonints"]-nonints]] ###last entry was metpath
				oddsratio,pvalue=fisher_exact(contingency_table,alternative="greater")
				if pvalue<significance_level*1/len(sum_ints_df):
					significance_df.loc[protpath,metapath]="*"
				else:
					significance_df.loc[protpath,metapath]=""
		#plot heatmap
		f,ax=plt.subplots()#figsize=25,21))
		sns.heatmap(heatmap,ax=ax,annot=significance_df,fmt="",square=square,robust=True,xticklabels=True,yticklabels=True,cmap="hot_r") #vmin=0,vmax=brite_df["number of interactions"].max(),linewidths=0.5
		f.tight_layout()
		#plt.subplots_adjust(bottom=0.5)
		plt.show()
		#save +adjust
		f.savefig("../overall_analysis/figures/enrichment/heatmap_pathways"+trim+self.approach+self.feature+"_"+method.__name__+"_"+classifier+".png")
		f.savefig("../overall_analysis/figures/enrichment/heatmap_pathways"+trim+self.approach+self.feature+"_"+method.__name__+"_"+classifier+".svg")
		return
	
	
	
	#wordcloud showing most abundand GO terms (2 words each) for interacting proteins
	def wordcloud_GOterms(self,method,classifier,predictions=None,width=450,height=150,min_font_size=6,appendix=""):
		#load predictions and annotations
		prot_annotations=pd.read_csv(self.analyses+"all_protein_annotations.tsv",sep="\t").set_index("Entry")
		if predictions is None:
			predictions=pd.read_csv(self.analyses+"predictions_"+method.__name__+self.feature+self.approach+".tsv",sep="\t")[[self.protcol,"metabolite",classifier]].set_index(self.protcol)
			predictions.rename(columns={classifier:"prediction"},inplace=True)
			#trim predictions to interactions only
			predictions=predictions[predictions["prediction"]==True]
		GO_columns=["Gene ontology (biological process)","Gene ontology (cellular component)","Gene ontology (GO)","Gene ontology (molecular function)"]
		#change here for OG-wise approach
		#change here for different GO_columns
		go_terms="; ".join(prot_annotations.loc[predictions.index,"Gene ontology (GO)"].dropna().tolist())
		stopwords={" of "," to "," by "," in "," from "," and "}
		go_terms_trimmed=go_terms.replace("GO","")
		for s in stopwords:
			go_terms_trimmed=go_terms_trimmed.replace(s," ")
		wordcloud=WordCloud(width=width,height=height,background_color='white',min_font_size=min_font_size).generate(go_terms_trimmed) 
		f=plt.figure(figsize=(8,8),facecolor=None) 
		plt.imshow(wordcloud) 
		plt.axis("off") 
		f.tight_layout() 
		plt.show() 
		f.savefig("../overall_analysis/figures/enrichment/wordcloud_GOterms"+self.feature+self.approach+method.__name__+classifier+appendix+".png")
		f.savefig("../overall_analysis/figures/enrichment/wordcloud_GOterms"+self.feature+self.approach+method.__name__+classifier+appendix+".svg")
		#f.savefig("../overall_analysis/figures/enrichment/wordcloud_GOterms"+da_yeast.feature+da_yeast.approach+method.__name__+classifier+"_broad.png")
		return
	
	
	
	#transforms predictions to uniprot-pubchemCID-metabolite_name-interaction dataframe
	def collect_predictions_with_ids(self,predictions,predictions_name,prediction_probabilities=None):
		predictions.reset_index(inplace=True)
		#translate metabolites
		string_pubchem_transdf=self.translate_CIDs_from_Stitch_to_PubChem(cids=predictions["metabolite"].unique()).set_index("Stitch_CID")
		predictions["metabolite_name"]=list(map(lambda x: string_pubchem_transdf.loc[x,"Name"],predictions["metabolite"]))
		predictions["PubChem_CID"]=list(map(lambda x: string_pubchem_transdf.loc[x,"PubChem_CID"],predictions["metabolite"]))
		if self.proteinwise==True:
			predictions.rename(columns={"metabolite":"Stitch_CID","prediction":"interaction","protein":"UniProtID"},inplace=True)
			#translate proteins
			organism=self.expfiles[0].split("/")[-1].split("_")[0]
			prot_trans_df=pd.read_csv("../databases/"+organism+"_ID_collection.tsv",sep="\t").set_index("ID")
			predictions["STRING_ID"]=list(map(lambda x: prot_trans_df.loc[x,"STRING_ID"] if x in prot_trans_df.index else "",predictions["UniProtID"]))
			predictions["GENENAME"]=list(map(lambda x: prot_trans_df.loc[x,"GENENAME"] if x in prot_trans_df.index else "",predictions["UniProtID"]))
			#OGs
			og_mappingdf=pd.read_csv("../databases/"+organism+"_string_uniprot_cog.tsv",sep="\t").set_index("UniProt_ID")
			predictions["Orthologous_Group"]=list(map(lambda x: og_mappingdf.loc[x,"COG"] if x in og_mappingdf.index else "",predictions["UniProtID"]))
			predictions=predictions.reindex(columns=["UniProtID","STRING_ID","GENENAME","Orthologous_Group","metabolite_name","PubChem_CID","Stitch_CID","interaction"])
			if prediction_probabilities is not None:
				predictions["probability_interaction"]=list(map(lambda x: prediction_probabilities.loc[x,prediction_probabilities.columns[0]],[tuple(y) for y in predictions[["UniProtID","Stitch_CID"]].values]))
		else:
			predictions.rename(columns={"metabolite":"Stitch_CID","prediction":"interaction","OG":"Orthologous_Group"},inplace=True)
			complete_set=pd.read_csv(self.databases+"complete_set.tsv",sep="\t").dropna()
			complete_set.set_index(complete_set.columns[0],inplace=True)
			predictions["UniProt_IDs"]=list(map(lambda x: ";".join(set((";".join(complete_set.loc[x])).split(";"))),predictions["Orthologous_Group"]))
			predictions=predictions.reindex(columns=["Orthologous_Group","UniProt_IDs","metabolite_name","PubChem_CID","Stitch_CID","interaction"])
		#save
		predictions.to_csv("../overall_analysis/"+predictions_name+".tsv",sep="\t",index=False,header=True)
		return
	
	
	
	#adds Stitch CIDs to metabolite raw data sheets extracted manually to supplemental_datasets_metabolites_raw.xlsx
	def add_CIDs_to_raw_sheets(self):
		inputfile="../supplemental_datasets_metabolites_raw.xlsx"
		excelfile=pd.ExcelFile(inputfile)
		sheets=excelfile.sheet_names
		writer=pd.ExcelWriter(inputfile)
		for sheet in sheets:
			original=excelfile.parse(sheet)
			ids=self.translate_metabolites_offline(original["Metabolite Names"].tolist(),online_translations=False)
			new=original.set_index("Metabolite Names")
			for metabolite in new.index:
				if metabolite in ids.index and ids.loc[metabolite,"CID"]!="":
					new.loc[metabolite,"Stitch_CID"]=int(ids.loc[metabolite,"CID"])
				else:
					new.loc[metabolite,"Stitch_CID"]=""
			#change order of columns
			new=new.reset_index()[["Stitch_CID"]+original.columns.tolist()]
			#save dataframes to sheets
			new.to_excel(writer,sheet+"_new",index=False)
		writer.save()
		print("Change back column headers and formatting manually in Excel!")
		return
	






###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################


if __name__=='__main__':
	###################################################
	#Programmatic workflow for one analysis experiment (e.g. S.cerevisiae)
	###################################################
	
	################################################
	#Initialize parameters
	proteinwise=True
	method=np.min
	methods=[np.min]
	appendix="_yeast/"
	organism="S.cerevisiae"
	feature="_same_crosscorr_binsize4_not_normalized"
	normalized="max" #normalize profiles to "normalized" , e.g. sum, max, min; needs to be pandas method as string, default="" (equals "sum")
	normalized2="max"
	organisms=[("9606","human"),("4932","yeast"),("511145","ecoli")] #for orthogonal query on PubChem
	confidence="low" #confidence from own manual selection of gold standard for training set
	approach="_maxmaxprofiles_unbalanced_low_confidence"
	classifier="svm_clf(C=0.1,kernel='linear')"
	classifiers=[classifier]
	score_cutoff=800 #score threshold for the STITCH database (for creation of positive set)
	################################################
	
	################################################
	#Prepare Data (since in the excel file with supplemental datasets are many datasets, loading it every time may overcharge RAM and be time consuming, all datasets are put into separate excel files.)
	supplemental_datasets="../supplemental_datasets_ids.xlsx" #adapt here were supplemental datasets are saved (default: parent folder)
	PD=PrepareData()
	PD.create_folders(appendix=appendix)
	PD.prepare_expdata(organism=organism,supplemental_datasets=supplemental_datasets)
	PD.prepare_databases(supplemental_datasets)
	################################################
	
	################################################
	#construct complete set to make predictions on
	
	cs=CompleteSet(simulation=False,overwrite=False,experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix,feature=feature,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
	complete_set=cs.create_complete_set_profiles()
	##stacking profiles
	complete_set_profiles=cs.stack_profiles(complete_set,method=method,profiles="raw")
	complete_set_profiles.to_csv(cs.databases+"complete_set_"+cs.normalized+"profiles_"+method.__name__+".tsv",sep="\t")
	##normalizing profiles
	complete_set_normprofiles=cs.normalize_concat_profiles_per_org(complete_set_profiles).dropna() #dropping NaN, e.g. in case of multiple profiles per molecule
	complete_set_normprofiles.to_csv(cs.databases+"complete_set_"+cs.normalized+cs.normalized2+"normprofiles_"+method.__name__+".tsv",header=True,index=True,sep="\t")
	##feature engineering on profiles
	cs.refeaturize_x_set(complete_set_normprofiles).to_csv(cs.databases+"complete_set_"+cs.normalized+cs.normalized2+"normprofiles"+cs.feature+"_"+method.__name__+".tsv",header=True,index=True,sep="\t")
	print("set of metabolite protein pairs to make predictions on created")
	
	
	################################################
	#provide ID translations as well as mapping to COGs for creation of gold standard
	db=DBHandler(organism+".xlsx",score_cutoff=score_cutoff,overwrite=False,simulation=False,experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix)
	db.get_orthologs(save=True)
	################################################
	#Creation of gold standard (positive and negative training subset can be created simultaneously)
	################################################
	#build positive subset of training set
	ps=PositiveSet(simulation=False,overwrite=False,experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix,feature=feature,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
	##construct positive set from STITCH
	if ps.proteinwise==True:
		positive_set=ps.load_set(ps.databases+"positive_set.tsv").reset_index().set_index(ps.protcol,drop=False)
	else:
		positive_set=ps.load_set(ps.databases+"positive_set.tsv")##stacking profiles
	positive_set_profiles=ps.stack_profiles(positive_set,method=method,profiles="raw")
	positive_set_profiles.to_csv(ps.databases+"positive_set_"+ps.normalized+"profiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
	##normalizing profiles
	positive_set_normprofiles=ps.normalize_concat_profiles_per_org(positive_set_profiles)
	positive_set_normprofiles.to_csv(ps.databases+"positive_set_"+ps.normalized+ps.normalized2+"normprofiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
	##feature engineering on profiles
	positive_set_refprofiles=ps.refeaturize_x_set(positive_set_normprofiles)
	positive_set_refprofiles.to_csv(ps.databases+"positive_set_"+ps.normalized+ps.normalized2+"normprofiles"+ps.feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
	
	##construct positive set from PubChem
	for organism in organisms:
		ps_pubchem=PositiveSet_PubChem(overwrite=False,experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix,db_organism=organism,feature=feature,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
		positive_set_pubchem=ps_pubchem.load_set(ps_pubchem.databases+"positive_og_set_"+ps_pubchem.db_organism[1]+"_pubchem.tsv")
		if ps_pubchem.proteinwise==True:
			positive_set_pubchem=ps_pubchem.proteinwise_x_set(positive_set_pubchem,method=method,corr_thresh=0.5)
			positive_set_pubchem.to_csv(ps_pubchem.databases+"positive_set_"+ps_pubchem.db_organism[1]+"_pubchem.tsv",sep="\t",index=False,header=True)
		##stack profiles
		positive_set_pubchem_profiles=ps_pubchem.stack_profiles(positive_set_pubchem,method=method,profiles="raw")
		positive_set_pubchem_profiles.to_csv(ps_pubchem.databases+"positive_set_"+ps_pubchem.db_organism[1]+"_pubchem_"+ps_pubchem.normalized+"profiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		##normalize profiles
		positive_set_pubchem_normprofiles=ps_pubchem.normalize_concat_profiles_per_org(positive_set_pubchem_profiles)
		positive_set_pubchem_normprofiles.to_csv(ps_pubchem.databases+"positive_set_"+ps_pubchem.db_organism[1]+"_pubchem_"+ps_pubchem.normalized+ps_pubchem.normalized2+"normprofiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
	print("positive training subset created")
	################################################
	
	################################################
	#Build negative subset of training set
	for organism in organisms:
		ns=NegativeSet(overwrite=False,experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix,db_organism=organism,feature=feature,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
		negative_set=ns.load_set(ns.databases+"negative_og_set_"+ns.db_organism[1]+".tsv").dropna() ### dropna() new due to missing uniprot
		if ns.proteinwise==True:
			negative_set=ns.proteinwise_x_set(negative_set,method=method,corr_thresh=0.5)
			negative_set.to_csv(ns.databases+"negative_set_"+ns.db_organism[1]+".tsv",sep="\t",index=False,header=True)
		##stack profiles
		negative_set_profiles=ns.stack_profiles(negative_set,method=method,profiles="raw")
		negative_set_profiles.to_csv(ns.databases+"negative_set_"+ns.db_organism[1]+"_"+ns.normalized+"profiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		#normalize profiles
		negative_set_normprofiles=ns.normalize_concat_profiles_per_org(negative_set_profiles)
		negative_set_normprofiles.to_csv(ns.databases+"negative_set_"+ns.db_organism[1]+"_"+ns.normalized+ns.normalized2+"normprofiles_"+method.__name__+".tsv",index=True,header=True,sep="\t")
		#feature engineering
		negative_set_refnormprofiles=ns.refeaturize_x_set(negative_set_normprofiles)
		negative_set_refnormprofiles.to_csv(ns.databases+"negative_set_"+ns.db_organism[1]+"_"+ns.normalized+ns.normalized2+"normprofiles"+ns.feature+"_"+method.__name__+".tsv",index=True,header=True,sep="\t")
	print("negative training subset created")
	################################################
	
	################################################
	#Merge subsets to training set
	ts=TrainingSet(methods=methods,db_orglist=organisms,experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix,feature=feature,normalized=normalized,normalized2=normalized2,proteinwise=proteinwise)
	#merge subsets from stitch and pubchem + feature engineering
	positive_set=ts.load_merged_x_set(x="positive",method=method)
	negative_set=ts.load_merged_x_set(x="negative",method=method)
	training_set_unbalanced=ts.get_unbalanced_training_set(method=method,extend_with_sampling=False,metabolite_extension=False,manual_selection=True,confidence=confidence)
	print("gold standard created")
	################################################
	
	################################################
	#machine learning 
	mlc=MLClassifiers(experimental="../experimental_data"+appendix,databases="../databases"+appendix,analyses="../analyses"+appendix,approach=approach,feature=feature,normalized=normalized+normalized2,proteinwise=proteinwise)
	#make predictions
	mlc.make_predictions([method],classifiers,reps=1)
	mlc.make_prediction_probabilities([method],classifiers,reps=1)
	print("predictions made")
	'''
	#performance evaluation on measures
	summary_table=mlc.measure_table(methods,classifiers)
	#evaluation on learning and ROC curve
	mlc.create_learning_ROC_curves_comparison(method,classifier,[0.3,0.4,0.5,0.6,0.7,0.8],reps=1000,appendix="_svmlinreg_c0p1")
	#rank proteins and metabolites
	mlc.rank_proteins(methods=[method],classifiers=classifiers)
	mlc.rank_metabolites(methods=[method],classifiers=classifiers)
	#comparison to random training
	mlc_yeast.make_random_predictions(method=method,classifier=classifier,reps=1000)
	mlc.boxplots_kappas_random_norm(method=method,classifier=classifier,appendix="_"+method.__name__+"_svmlinreg0p1")
	'''
	################################################
	
	################################################
	#network analysis
	na=NetworkAnalysis(appendix=appendix,method=method,approach=approach,feature=feature,trim="",classifier=classifier,normalized=normalized,proteinwise=proteinwise)
	#degree distributions
	#na.degree_distribution()
	#comparison of predictions from different classifiers, feature engineerings, pooling methods, normalizations, protein-wise/OG-wise with na.venn()
	
	
	################################################
	'''
	FIGURES and TABLES in manuscript
	################################################
	#data analysis
	approach_ara="_maxmaxprofiles_unbalanced_low_confidence"
	approach_yeast="_maxmaxprofiles_unbalanced_low_confidence"
	approach_euc="_maxmaxprofiles_unbalanced_low_confidence"
	appendix_ara="_ara_protwise/"
	appendix_ara_og="_Ara_only2/"
	appendix_yeast_og="_yeast_only2/"
	appendix_yeast="_yeast_protwise/"
	appendix_euc="_eucaryotes2/"
	
	##################################
	# Initialize Objects and load Sets and predictions
	##################################
	#initialize DataAnalysis objects
	da_ara=DataAnalysis(appendix=appendix_ara,approach=approach_ara,feature=feature,normalized=normalized,proteinwise=proteinwise)
	da_ara_og=DataAnalysis(appendix=appendix_ara_og,approach=approach_ara,feature=feature,normalized=normalized,proteinwise=False)
	da_yeast=DataAnalysis(appendix=appendix_yeast,approach=approach_yeast,feature=feature,normalized=normalized,proteinwise=proteinwise)
	da_yeast_og=DataAnalysis(appendix=appendix_yeast_og,approach=approach_yeast,feature=feature,normalized=normalized,proteinwise=False)
	da_euc=DataAnalysis(appendix=appendix_euc,approach=approach_euc,feature=feature,normalized=normalized,proteinwise=False)
	#initialize Complete Set objects
	ara_cs=CompleteSet(simulation=False,overwrite=False,experimental="../experimental_data"+appendix_ara,databases="../databases"+appendix_ara,analyses="../analyses"+appendix_ara,feature=feature,normalized=normalized,proteinwise=proteinwise)
	yeast_cs=CompleteSet(simulation=False,overwrite=False,experimental="../experimental_data"+appendix_yeast,databases="../databases"+appendix_yeast,analyses="../analyses"+appendix_yeast,feature=feature,normalized=normalized,proteinwise=proteinwise)
	#initialize NetworkAnalysis instances
	na_ara=NetworkAnalysis(appendix=appendix_ara,method=method,approach=approach_ara,feature=feature,trim="",classifier=classifier,normalized=normalized,proteinwise=proteinwise)
	na_ara_og=NetworkAnalysis(appendix=appendix_ara_og,method=method,approach=approach_ara,feature=feature,trim="",classifier=classifier,normalized=normalized,proteinwise=False)
	na_euc=NetworkAnalysis(appendix=appendix_euc,method=method,approach=approach_euc,feature=feature,trim="",classifier=classifier,normalized=normalized,proteinwise=False)
	na_yeast=NetworkAnalysis(appendix=appendix_yeast,method=method,approach=approach_yeast,feature=feature,trim="",classifier=classifier,normalized=normalized,proteinwise=True)
	na_yeast_og=NetworkAnalysis(appendix=appendix_yeast_og,method=method,approach=approach_yeast,feature=feature,trim="",classifier=classifier,normalized=normalized,proteinwise=False)
	#load predictions
	ara_predictions=na_ara.predictions.set_index([na_ara.protcol,"metabolite"])
	ara_predictions.columns=["prediction"]
	ara_og_predictions=na_ara_og.predictions.set_index([na_ara_og.protcol,"metabolite"])
	ara_og_predictions.columns=["prediction"]
	euc_predictions=na_euc2.predictions.set_index([na_euc2.protcol,"metabolite"])
	euc_predictions.columns=["prediction"]
	yeast_predictions=na_yeast.predictions.set_index([na_yeast.protcol,"metabolite"])
	yeast_predictions.columns=["prediction"]
	yeast_og_predictions=na_yeast_og.predictions.set_index([na_yeast_og.protcol,"metabolite"])
	yeast_og_predictions.columns=["prediction"]
	
	##################################
	# summary table experimental data files
	##################################
	df_euc_files=da_euc.construct_df_files().sort_index()
	df_euc_files.to_csv(da_euc.analyses+"data_summary_files.tsv",header=True,index=True,sep="\t")
	df_euc_files=df_euc_files.rename(columns={"proteins": "Proteins","metabolites": "Metabolites", "orthologous groups": "Orthologous Groups"},index={"intersection": "Intersection", "union":"Union"})
	df_ara_files=da_ara.construct_df_files().sort_index()
	df_ara_files.to_csv(da_ara.analyses+"data_summary_files.tsv",header=True,index=True,sep="\t")
	df_ara_files=df_ara_files.rename(columns={"proteins": "Proteins","metabolites": "Metabolites", "orthologous groups": "Orthologous Groups"},
	index={"A.thaliana_cell_cultures": "A. thaliana cell culture","A.thaliana_rosettes": "A. thaliana rosettes", "A.thaliana_seedlings_4_rep": "A. thaliana seedlings","intersection": "A. thaliana intersection","union":"A. thaliana union"})
	df_ara_og_files=da_ara_og.construct_df_files().sort_index()
	df_ara_og_files.to_csv(da_ara_og.analyses+"data_summary_files.tsv",header=True,index=True,sep="\t")
	df_ara_og_files=df_ara_og_files.rename(columns={"proteins": "Proteins","metabolites": "Metabolites", "orthologous groups": "Orthologous Groups"},
	index={"A.thaliana_cell_cultures": "A. thaliana cell culture","A.thaliana_rosettes": "A. thaliana rosettes", "A.thaliana_seedlings_4_rep": "A. thaliana seedlings","intersection": "A. thaliana intersection","union":"A. thaliana union"})
	df_yeast_files=da_yeast.construct_df_files().sort_index()
	df_yeast_files.to_csv(da_yeast.analyses+"data_summary_files.tsv",header=True,index=True,sep="\t")
	df_yeast_files=df_yeast_files.rename(columns={"proteins": "Proteins","metabolites": "Metabolites", "orthologous groups": "Orthologous Groups"},
	index={"S.cerevisiae_cell_cultures": "S. cerevisiae cell culture","S.cerevisiae_growthphases_condition1": "S. cerevisiae growing condition 1", "S.cerevisiae_growthphases_condition2": "S. cerevisiae growing condition 2", "S.cerevisiae_growthphases_condition3": "S. cerevisiae growing condition 3", "intersection": "S. cerevisiae intersection","union":"S. cerevisisae union"})
	df_yeast_og_files=da_yeast_og.construct_df_files().sort_index()
	df_yeast_og_files.to_csv(da_yeast_og.analyses+"data_summary_files.tsv",header=True,index=True,sep="\t")
	df_yeast_og_files=df_yeast_og_files.rename(columns={"proteins": "Proteins","metabolites": "Metabolites", "orthologous groups": "Orthologous Groups"},
	index={"S.cerevisiae_cell_cultures": "S. cerevisiae cell culture","S.cerevisiae_growthphases_condition1": "S. cerevisiae growing condition 1", "S.cerevisiae_growthphases_condition2": "S. cerevisiae growing condition 2", "S.cerevisiae_growthphases_condition3": "S. cerevisiae growing condition 3", "intersection": "S. cerevisiae intersection","union":"S. cerevisisae union"})
	file_summary=df_ara_files.append(df_yeast_files).append(df_euc_files.loc["Intersection"]).append(df_euc_files.loc["Union"])
	file_summary=file_summary[["Orthologous Groups", "Proteins","Metabolites"]]
	file_summary.loc[df_ara_og_files.index,"Orthologous Groups"]=df_ara_og_files.loc[:,"Orthologous Groups"]
	file_summary.loc[df_yeast_og_files.index,"Orthologous Groups"]=df_yeast_og_files.loc[:,"Orthologous Groups"]
	file_summary.to_csv("../overall_analysis/file_summary.tsv",sep="\t",header=True,index=True)
	
	##################################
	#summary table for training sets
	##################################
	df_ara_sets=da_ara.construct_df_sets(method=method)
	df_ara_sets.to_csv(da_ara.analyses+"data_summary_training_sets"+da_ara.approach+".tsv",header=True,index=True,sep="\t")
	df_ara_og_sets=da_ara_og.construct_df_sets(method=method)
	df_ara_og_sets.to_csv(da_ara_og.analyses+"data_summary_training_sets"+da_ara_og.approach+".tsv",header=True,index=True,sep="\t")
	df_yeast_sets=da_yeast.construct_df_sets(method=method)
	df_yeast_sets.to_csv(da_yeast.analyses+"data_summary_training_sets"+da_yeast.approach+".tsv",header=True,index=True,sep="\t")
	df_yeast_og_sets=da_yeast_og.construct_df_sets(method=method)
	df_yeast_og_sets.to_csv(da_yeast_og.analyses+"data_summary_training_sets"+da_yeast_og.approach+".tsv",header=True,index=True,sep="\t")
	df_euc_sets=da_euc.construct_df_sets(method=method)
	df_euc_sets.to_csv(da_euc.analyses+"data_summary_training_sets"+da_euc.approach+".tsv",header=True,index=True,sep="\t")
	#concatenate them
	df_sets_index=3*["A.thaliana"]+3*["S.cerevisiae"]+3*["Both"]
	df_sets=df_ara_sets.append(df_yeast_sets).append(df_euc_sets)
	df_sets["experiment"]=df_sets_index
	df_sets=df_sets.reset_index().set_index("experiment")
	df_sets.to_csv("../overall_analysis/sets_summary.tsv",sep="\t",header=True,index=True)
	
	##################################
	# description of OG annotations
	##################################
	#Heatmaps are independent of method, approach, feature
	#initialize Arabidopsis MLClassifier
	approach_og="_unbalanced_sampled_low_confidence"
	feature_og="_same_crosscorr_not_normalized"
	normalized_og=""
	appendix_ara_og="_Ara_only2/"
	mlc_ara=MLClassifiers(experimental="../experimental_data"+appendix_ara_og,databases="../databases"+appendix_ara_og,analyses="../analyses"+appendix_ara_og,approach=approach_og,feature=feature_og,normalized=normalized_og)
	#initialize Yeast MLClassifier
	appendix_yeast_og="_yeast_only2/"
	mlc_yeast=MLClassifiers(experimental="../experimental_data"+appendix_yeast_og,databases="../databases"+appendix_yeast_og,analyses="../analyses"+appendix_yeast_og,approach=approach_og,feature=feature_og,normalized=normalized_og)
	#both classifier
	mlc_euc=MLClassifiers(experimental="../experimental_data"+appendix_euc,databases="../databases"+appendix_euc,analyses="../analyses"+appendix_euc,approach=approach_euc,feature=feature,normalized=normalized)
	#create figure
	mlcs=[mlc_ara,mlc_yeast]
	mlcs_cols23=[mlc_euc]
	colors_euc=["tab:red","tab:blue"]
	#da_ara.corr_ogs_heat_hist_pies(mlcs,method=np.min)
	da_ara.corr_ogs_heat_hist_pies(mlcs,method=np.min,sec_col=False,third_col=False)
	da_euc.corr_ogs_heat_hist_pies(mlcs_cols23,method=np.min,first_col=False,sec_col=True,third_col=True,colors=colors_euc,titles=[""])
	plt.rcParams.update({'font.size': 8})
	da_euc.corr_ogs_heat_hist_pies(mlcs_cols23,method=np.min,first_col=False,sec_col=False,third_col=True,colors=colors_euc,titles=[""])
	#da_ara.corr_ogs_heat_hist_pies(mlcs_cols23,method=np.min,first_col=False,sec_col=False)
	
	##############################################
	# FREQUENCY OF INTERACTIONS
	###########################################
	len(ara_predictions[ara_predictions["prediction"]==True])/len(ara_predictions)
	len(ara_og_predictions[ara_og_predictions["prediction"]==True])/len(ara_og_predictions)
	len(yeast_predictions[yeast_predictions["prediction"]==True])/len(yeast_predictions)
	len(yeast_og_predictions[yeast_predictions["prediction"]==True])/len(yeast_og_predictions)
	len(euc_predictions[euc_predictions["prediction"]==True])/len(euc_predictions)
	
	###########################################
	# COMPARISON TO EXPERIMENTAL DATA
	###########################################
	#Comparison of Ara and Euc predictions to Ewkas set
	da_ara.comparison_to_Ewkas_set(normalized=normalized,feature=feature,method=method,approach=approach_ara,classifier=classifier_svmlin_reg0p1,confidence="low")
	da_euc.comparison_to_Ewkas_set(normalized=normalized,feature=feature,method=method,approach=approach_euc,classifier=classifier,confidence="low")
	#Comparison of Yeast and Euc predictions to AP and TPP
	da_yeast.comparison_to_ap_tpp_pcc(yeast_predictions,appendix="_yeast_maxmax_protwise_min_svmlin_no_sampling")
	da_yeast.comparison_to_ap_tpp_pcc(yeast_linreg0p1_predictions,appendix="_yeast_maxmax_protwise_min_svmlinreg0p1_no_sampling")
	da_yeast.comparison_to_ap_tpp_pcc(yeast_rbf_predictions,appendix="_yeast_maxmax_protwise_min_svmrbf_no_sampling")
	#Comparison Tyr-Asp targets
	da_ara.comparison_to_ap_tpp_pcc_tyrasp_targets(ara_predictions)
	ara_og_protwise_predictions=na_ara_og.ogwise_to_protwise_preds().set_index(["protein","metabolite"])
	ara_og_protwise_predictions.columns=["prediction"]
	da_ara.comparison_to_ap_tpp_pcc_tyrasp_targets(ara_og_protwise_predictions,appendix="_OGwise")
	
	#################################################################
	# CLASSIFICATION OF METABOLITES AND PATHWAY ANNOTATION WITH KEGG => DONE
	#################################################################
	#PUBCHEM
	#get metabolites per experiment
	ara_metabolites=ara_predictions.index.get_level_values(1).unique().tolist()
	yeast_metabolites=yeast_predictions.index.get_level_values(1).unique().tolist()
	euc_metabolites=euc_predictions.index.get_level_values(1).unique().tolist()
	#translate them to PubChem CIDs
	ara_metabolite_names=ara_cs.translate_CIDs_from_Stitch_to_PubChem(ara_metabolites)
	yeast_metabolite_names=yeast_cs.translate_CIDs_from_Stitch_to_PubChem(yeast_metabolites)
	yeast_only_metabolite_names=yeast_cs.translate_CIDs_from_Stitch_to_PubChem(set(yeast_metabolites)-set(euc_metabolites))
	euc_metabolite_names=euc_cs.translate_CIDs_from_Stitch_to_PubChem(euc_metabolites)
	ara_cids=ara_metabolite_names["PubChem_CID"].tolist()
	yeast_cids=yeast_metabolite_names["PubChem_CID"].tolist()
	euc_cids=euc_metabolite_names["PubChem_CID"].tolist()
	#classify them
	ara_metabolite_df=ara_cs.retrieve_CID_classes(ara_cids)
	yeast_metabolite_df=yeast_cs.retrieve_CID_classes(yeast_cids)
	euc_metabolite_df=euc_cs.retrieve_CID_classes(euc_cids)
	yeast_only_metabolite_names.to_csv("../analyses_yeast_only2/metabolites_yeast_only.tsv",sep="\t",header=True,index=True)
	
	#KEGG
	metabolite_kegg_df=pd.read_csv("../databases/metabolite_kegg_classes.tsv",sep="\t")
	metabolite_kegg_ids=metabolite_kegg_df["KEGG_ID"].dropna().tolist()
	pathway_df=ara_cs.find_pathways_metabolites(metabolite_kegg_ids)
	for i in pathway_df.index:
		metabolite_kegg_df.loc[metabolite_kegg_df["KEGG_ID"]==i,"KEGG pathways"]=pathway_df.loc[i,"KEGG pathways"]
	metabolite_kegg_df.to_csv("../databases/metabolite_kegg_classes_pathways.tsv",sep="\t",index=False,header=True)
	metabolite_kegg_df.set_index("CID",inplace=True)
	ara_metabolite_kegg_df=metabolite_kegg_df.loc[ara_metabolites]
	ara_metabolite_pathways=set(";".join(ara_metabolite_kegg_df["KEGG pathways"].dropna()).split(";"))
	len(ara_metabolite_pathways)
	#ara: metabolites from 87 pathways
	yeast_metabolite_kegg_df=metabolite_kegg_df.loc[yeast_metabolites]
	yeast_metabolite_pathways=set(";".join(yeast_metabolite_kegg_df["KEGG pathways"].dropna()).split(";"))
	len(yeast_metabolite_pathways)
	#yeast: metabolites from 124 pathways
	
	#################################################################
	#CLASSIFICATION OF PROTEINS AND PATHWAY ANNOTATION WITH KEGG
	#################################################################
	##get proteins per experiment, OG-wise
	euc_ogs=euc_predictions.index.get_level_values(0).unique().tolist()
	#get proteins per experiment protwise
	ara_proteins=ara_predictions.index.get_level_values(0).unique().tolist()
	yeast_proteins=yeast_predictions.index.get_level_values(0).unique().tolist()
	ara_kegg_df=ara_cs.classify_proteins(ara_proteins,"ath")
	yeast_kegg_df=yeast_cs.classify_proteins(yeast_proteins,"sce")
	#for euc_kegg_df, delete organismprefixes
	for i in euc_kegg_df.index:
		if not pd.isna(euc_kegg_df.loc[i,"BRITE classes"]):
			euc_kegg_df.loc[i,"BRITE classes"]=";".join(list(map(lambda x: "br:"+x[6:],euc_kegg_df.loc[i,"BRITE classes"].split(";"))))
		if not pd.isna(euc_kegg_df.loc[i,"KEGG pathways"]):
			euc_kegg_df.loc[i,"KEGG pathways"]=";".join(list(map(lambda x: "path:"+x[8:],euc_kegg_df.loc[i,"KEGG pathways"].split(";"))))
	
	#decide which classes to take
	all_classes=";".join(ara_kegg_df["BRITE classes"].replace(np.nan,"").tolist()).split(";")
	all_classes_trimmed=list(filter(lambda x: x!="br:ath00001" and x!="br:ath01000",all_classes))
	plt.hist(all_classes_trimmed)
	plt.show()
	c=Counter(all_classes)
	
	#assign EC classes proteinwise
	ara_cs.collect_uniprot_IDs(xset=None,predictions=ara_predictions)
	ara_kegg_ec_df=ara_cs.assign_EC_numbers_to_proteins(df=ara_kegg_df)
	ara_kegg_ec_df.to_csv(ara_cs.analyses+"UniProt_BRITEclasses_KEGGpathways.tsv",sep="\t",index=True,header=True)
	yeast_cs.collect_uniprot_IDs(xset=None,predictions=yeast_predictions)
	yeast_kegg_ec_df=yeast_cs.assign_EC_numbers_to_proteins(df=yeast_kegg_df)
	yeast_kegg_ec_df.to_csv(yeast_cs.analyses+"UniProt_BRITEclasses_KEGGpathways.tsv",sep="\t",index=True,header=True)
	
	#ara: 119 pathways
	ara_protein_pathways=set(";".join(ara_kegg_ec_df["KEGG pathways"].dropna()).split(";"))
	len(ara_protein_pathways)
	#yeast: 99 pathways annotated to proteins
	yeast_protein_pathways=set(";".join(yeast_kegg_ec_df["KEGG pathways"].dropna()).split(";"))
	len(yeast_protein_pathways)
	
	#intersecting pathways proteins and metabolites
	ara_metabolite_pathway_numbers=set(map(lambda x: x[8:],ara_metabolite_pathways))
	ara_protein_pathway_numbers=set(map(lambda x: x[8:],ara_protein_pathways))
	len(set(ara_metabolite_pathway_numbers & ara_protein_pathway_numbers))
	#ara: 29 pathways in intersection
	yeast_metabolite_pathway_numbers=set(map(lambda x: x[8:],yeast_metabolite_pathways))
	yeast_protein_pathway_numbers=set(map(lambda x: x[8:],yeast_protein_pathways))
	len(set(yeast_metabolite_pathway_numbers & yeast_protein_pathway_numbers))
	#yeast: 20 pathways in intersection
	
	#count BRITE class assignments
	metabolite_kegg_df=pd.read_csv("../databases/metabolite_kegg_classes_pathways.tsv",sep="\t").set_index("CID")
	#yeast
	yeast_protein_kegg_df=yeast_cs.classify_proteins(uniprot_ids=None,org=None)
	len(set(";".join(yeast_protein_kegg_df["BRITE classes"].dropna()).split(";"))) #37 protein brite classes
	yeast_metabolites=yeast_predictions.index.get_level_values(1).unique().tolist()
	yeast_metabolite_classes=metabolite_kegg_df.loc[yeast_metabolites,"My class"].unique() #8 metabolite classes
	#ara
	ara_protein_kegg_df=ara_cs.classify_proteins(uniprot_ids=None,org=None)
	len(set(";".join(ara_protein_kegg_df["BRITE classes"].dropna()).split(";"))) #39 protein brite classes
	ara_metabolites=ara_predictions.index.get_level_values(1).unique()
	ara_metabolite_classes=metabolite_kegg_df.loc[ara_metabolites,"My class"].unique() #8 metabolite classes
	
	#count GO Terms
	#yeast
	prot_annotations=pd.read_csv(yeast_cs.analyses+"all_protein_annotations.tsv",sep="\t").set_index("Entry")
	go_terms="; ".join(prot_annotations.loc[yeast_predictions.index.get_level_values(0),"Gene ontology (GO)"].dropna().tolist())
	len(set(go_terms.split(";"))) #3264 GO Terms
	#ara
	prot_annotations=pd.read_csv(ara_cs.analyses+"all_protein_annotations.tsv",sep="\t").set_index("Entry")
	go_terms="; ".join(prot_annotations.loc[ara_predictions.index.get_level_values(0),"Gene ontology (GO)"].dropna().tolist())
	len(set(go_terms.split(";"))) #3997 GO Terms
	
	#################################################################
	#classify OGs and annotate pathways with KEGG
	#################################################################
	#assign classes and pathways to OGs (KEGG)
	ara_kegg_df_OG=pd.read_csv("../analyses_Ara_only2/uniprot_kegg_briteclasses_keggpathways.tsv",sep="\t").set_index("UniProt ID")
	yeast_kegg_df_OG=pd.read_csv("../analyses_yeast_only2/uniprot_kegg_briteclasses_keggpathways.tsv",sep="\t").set_index("UniProt ID")
	euc_kegg_df=ara_kegg_df_OG.append(yeast_kegg_df_OG) #4435 proteins
	euc_complete_set_classes_pathways=euc_cs.assign_classes_and_pathways_to_OGs(complete_set=euc_complete_set.copy(),kegg_df=euc_kegg_df)
	#assign EC numbers to OGs
	euc_complete_set_classes_pathways_ec=euc_cs.assign_EC_numbers_to_OGs(euc_complete_set_classes_pathways)
	#save dfs
	euc_complete_set_classes_pathways_ec.to_csv(euc_cs.analyses+"OG_UniProt_BRITEclasses_KEGGpathways.tsv",sep="\t",index=True,header=True)
	
	#################################################################
	#Make graph data of interactions among OG and metabolite classes
	#################################################################
	ara_predictions_annot=ara_predictions.copy()
	da_ara.make_graph_data_among_classes(method=method,predictions_raw=ara_predictions_annot,classifier=classifier)
	yeast_predictions_annot=yeast_linreg0p1_predictions.copy()
	da_yeast.make_graph_data_among_classes(method=method,predictions_raw=yeast_predictions_annot,classifier=classifier_svmlin_reg0p1)
	euc_predictions_annot=euc_predictions.copy()
	da_euc.make_graph_data_among_classes(method=method,predictions_raw=euc_predictions_annot,classifier=classifier)
	
	#################################
	#plot class distribution barplot
	#################################
	#das=[da_ara,da_yeast,da_euc]
	das=[da_ara,da_yeast]
	colors=[tableau_green,tableau_orange,tableau_red]
	approaches=[approach_ara,approach_yeast,approach_euc]
	names=["Ath","Sce","Euc"]
	da_ara.count_class_assignments(das=das,colors=colors,approaches=approaches,names=names,method=method,feature=None,classifier=classifier)
	da_ara.class_assignments_merged(color1=tableau_green,color2=tableau_orange,color3=tableau_red)
	#with both
	da_euc.count_class_assignments(das=[da_ara,da_yeast,da_euc],colors=colors,approaches=approaches,names=names,method=method,feature=None,classifier=classifier)
	plt.rcParams.update({'font.size': 12})
	da_euc.class_assignments_merged(color1=tableau_green,color2=tableau_orange,color3=tableau_red)
	
	#################################################################
	#Selection of classifier and pooling method
	#################################################################
	plt.rcParams.update({'font.size': 12})
	da_yeast.fig_selection_of_classifier(classifier_oi=classifier_svmlin,appendix="_svmlin")
	da_yeast.fig_selection_of_classifier(classifier_oi=classifier_svmlin_reg,appendix="_svmlinreg")
	da_yeast.fig_selection_of_classifier(classifier_oi=classifier_svmlin_reg0p1,method_norm="amin normalized",appendix="_svmlinreg0p1_bigsize")
	da_yeast.fig_selection_of_classifier(classifier_oi=classifier_svmrbf,appendix="_svmrbf")
	da_yeast.fig_selection_of_classifier_1PM(classifier_oi=classifier_svmlin_reg0p1,methods_notnorm=["amin","mean","amax"],method_norm="amin normalized",performance_measure="f_1 score")
	#comparison of different methods:
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=["svm_clf(C=1,kernel='linear')"],methods=[np.mean,np.median,np.min,np.max],appendix="_diff_methods_svmlin"+feature)
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=[classifier_svmlin_reg],methods=[np.mean,np.median,np.min,np.max],appendix="_diff_methods_svmlin_reg"+feature)
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=[classifier_svmrbf],methods=[np.mean,np.median,np.min,np.max],appendix="_diff_methods_svmrbf"+feature)
	#comparison of different classifiers
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=["svm_clf(C=1,kernel='linear')","svm_clf(kernel='rbf')","svm_clf(kernel='poly',degree=8)","svm_clf(kernel='sigmoid')","rf_clf()"],methods=[np.min],appendix="_diff_classifiers_min"+feature)
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=["svm_clf(C=1,kernel='linear')","svm_clf(kernel='rbf')","svm_clf(kernel='poly',degree=8)","svm_clf(kernel='sigmoid')","rf_clf()"],methods=[np.mean],appendix="_diff_classifiers_mean"+feature)
	#comparison linear kernel regularization
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=["svm_clf(C=0.01,kernel='linear')","svm_clf(C=0.1,kernel='linear')","svm_clf(C=1,kernel='linear')","svm_clf(C=10,kernel='linear')","svm_clf(C=100,kernel='linear')"],methods=[np.min],appendix="_svmlin_reg_min"+feature)
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=["svm_clf(C=0.01,kernel='linear')","svm_clf(C=0.1,kernel='linear')","svm_clf(C=1,kernel='linear')","svm_clf(C=10,kernel='linear')","svm_clf(C=100,kernel='linear')"],methods=[np.mean],appendix="_svmlin_reg_mean"+feature)
	#suppfig:
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=[feature],approaches=[approach_yeast],classifiers=["svm_clf(C=0.01,kernel='linear')","svm_clf(C=0.1,kernel='linear')","svm_clf(C=1,kernel='linear')","svm_clf(C=10,kernel='linear')","svm_clf(C=100,kernel='linear')","rf_clf()"],methods=[np.min],appendix="_min_selected",legend_roc=True,legend_dd=False)
	headers=["Learning Curves","ROC Curves","Degree Distributions","lin. SVM", "C=0.01","C=0.1","C=1","C=10","C=100","RF"]
	da_yeast.plot_headers(headers, plot_with_letters=False)
	
	#methods_notnorm=["amin","mean","amax"]
	methods_notnorm=["amin"]
	da_yeast.radar_chart_selection_classifier(methods_notnorm=methods_notnorm,method_norm="amin normalized",classifier=classifier_svmlin_reg0p1,colour="tab:red",appendix="_svmlinreg0p1",title="linear SVM, C=0.1")
	
	#all experiments
	plt.rcParams.update({'font.size': 16})
	methods_notnorm_og=["mean"]
	da_yeast.radar_chart_selection_classifier(methods_notnorm=methods_notnorm,method_norm="amin normalized",classifier=classifier_svmlin_reg0p1,colour="tab:orange",appendix="_yeast_svmlinreg0p1")
	da_yeast_og.radar_chart_selection_classifier(methods_notnorm=methods_notnorm_og,method_norm="amin normalized",classifier=classifier_svmlin_reg0p1,colour="tab:orange",appendix="_yeastOG_svmlinreg0p1_mean")
	da_ara.radar_chart_selection_classifier(methods_notnorm=methods_notnorm,method_norm="amin normalized",classifier=classifier_svmlin_reg0p1,colour="tab:green",appendix="_ara_svmlinreg0p1")
	da_ara_og.radar_chart_selection_classifier(methods_notnorm=methods_notnorm_og,method_norm="amin normalized",classifier=classifier_svmlin_reg0p1,colour="tab:green",appendix="_araOG_svmlinreg0p1_mean")
	da_euc.radar_chart_selection_classifier(methods_notnorm=methods_notnorm_og,method_norm="amin normalized",classifier=classifier_svmlin_reg0p1,colour="tab:red",appendix="_euc_svmlinreg0p1_mean")
	
	
	headers=["Protein-wise","OG-wise","S. cerevisiae","A. thaliana","Both"]
	da_yeast.plot_headers(headers, plot_with_letters=False)
	
	#################################################################
	#figure feature engineering
	#################################################################
	
	plt.rcParams.update({'font.size': 18})
	da_yeast.fig_feature_engineering_example(mode="full")
	da_yeast.fig_feature_engineering_weights()
	features=["_same_crosscorr_not_normalized","_same_crosscorr_binsize10_not_normalized","_same_crosscorr_binsize4_not_normalized","_same_crosscorr_binsize2_not_normalized"]
	da_yeast.fig_feature_engineering_performance(features=features,method=method,classifier=classifier_svmlin_reg0p1,feature_oi=feature,appendix="_bin4")
	#svmlin
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=features,approaches=[approach_yeast],classifiers=["svm_clf(C=1,kernel='linear')"],methods=[np.min],appendix="_diff_features_svmlin_min")
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=features,approaches=[approach_yeast],classifiers=["svm_clf(C=1,kernel='linear')"],methods=[np.mean],appendix="_diff_features_svmlin_mean")
	#svmlinreg
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=features,approaches=[approach_yeast],classifiers=["svm_clf(C=0.1,kernel='linear')"],methods=[np.min],appendix="_diff_features_svmlinreg0p1_min")
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=features,approaches=[approach_yeast],classifiers=["svm_clf(C=0.01,kernel='linear')"],methods=[np.mean],appendix="_diff_features_svmlinreg_mean")
	#rbf
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=features,approaches=[approach_yeast],classifiers=[classifier_svmrbf],methods=[np.min],appendix="_diff_features_svmrbf_min")
	da_yeast.classifier_comparison(experiments=[appendix_yeast],features=features,approaches=[approach_yeast],classifiers=[classifier_svmrbf],methods=[np.mean],appendix="_diff_features_svmrbf_mean")
	
	
	#################################################################
	#comparison to random
	#################################################################
	#LC and ROC and DD 
	da_yeast.comparison_to_random_training(appendix="_svmlinreg_c0p1")
	da_yeast.comparison_to_random_training(classifier=classifier_svmlin_reg0p1,appendix="_svmlinreg_c0p1",method=method,colours=["tab:red","tab:blue","tab:green","k"])
	colours=["tab:red","tab:purple","tab:green","k"]
	da_yeast.comparison_to_random_training(classifier=classifier_svmlin_reg0p1,appendix="_svmlinreg_c0p1",method=method,colours=colours)

	headers=["A Learning Curves","B ROC Curves","C Degree Distributions","E Accuracy","Assortativity","D Agreement"]
	da_yeast.plot_headers(headers, plot_with_letters=False)
	#acc and ass
	interaction_probabilities=[0.01,0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.5]
	da_yeast.comparison_to_random_network(classifier=classifier_svmlin_reg0p1,int_probs=interaction_probabilities,method=method,reps=1000)
	da_yeast.hist_acc_ass(method=method,classifier=classifier_svmlin_reg0p1)
	#venn
	tableau_red=(214/255,39/255,40/255,0.5) 
	tableau_blue=(31/255,119/255,180/255,0.5)
	tableau_green=(44/255,160/255,44/255,0.6)
	tableau_purple=(148/255,103/255,189/255,0.5)
	venn_colours=[tableau_red,tableau_purple,tableau_green]
	
	predictions_norm=yeast_linreg0p1_predictions.copy()
	predictions_training=pd.read_csv(da_yeast.analyses+"predictions_random_ts_"+method.__name__+da_yeast.feature+da_yeast.approach+"_"+classifier_svmlin_reg0p1+".tsv",sep="\t").set_index([da_yeast.protcol,"metabolite"])
	predictions_training.drop("num_true_ts",axis=1,inplace=True)
	predictions_training.columns=["prediction"]
	predictions_complete=pd.read_csv(da_yeast.analyses+"predictions_random_cs_"+method.__name__+da_yeast.feature+da_yeast.approach+"_"+classifier_svmlin_reg0p1+".tsv",sep="\t").set_index([da_yeast.protcol,"metabolite"])
	predictions_complete.drop("num_true_cs",axis=1,inplace=True)
	predictions_complete.columns=["prediction"]
	train_exp_intersect=set(set(yeast_linreg0p1_predictions[yeast_linreg0p1_predictions["prediction"]==True].index.tolist()) & set(predictions_training[predictions_training["prediction"]==True].index.tolist()))
	da_yeast.venn_diagram_randoms(predictions_exp=predictions_norm,predictions_training=predictions_training,predictions_complete=predictions_complete,colors=venn_colours,appendix="_linreg0p1_min",fontsize=22)
	
	predictions_yeast_trimmed=yeast_linreg0p1_predictions.drop(train_exp_intersect,axis=0)
	len(predictions_yeast_trimmed[predictions_yeast_trimmed["prediction"]==True])/len(predictions_yeast_trimmed)
	
	#networks
	da_yeast.make_graph_data_among_classes(method=method,predictions_raw=None,classifier=classifier_svmlin_reg0p1,random="training")
	da_yeast.make_graph_data_among_classes(method=method,predictions_raw=None,classifier=classifier_svmlin_reg0p1,random="complete")
	
	da_yeast.make_EC_graph(layout=net_layout,method=method,approach=approach_yeast,classifier=classifier_svmlin_reg0p1,feature=None,color1=tableau_blue,color2=tableau_green,random="training")
	da_yeast.make_EC_graph(layout=net_layout,method=method,approach=approach_yeast,classifier=classifier_svmlin_reg0p1,feature=None,color1=tableau_blue,color2=tableau_red,random="complete")
	#network for interaction probability=prediction frequency
	#to do or not
	
	
	
	######################################
	# Enrichment analysis
	######################################
	
	#networks
	net_layout=pd.read_csv("../analyses_Ara_only2/EC_MC_network_same_crosscorr_not_normalized_unbalanced_sampled_high_confidence_coordinates_1.tsv",sep="\t").set_index("Unnamed: 0")
	net_layout_large=pd.read_csv("../overall_analysis/figures/enrichment/EC_MC_network_same_crosscorr_binsize4_not_normalized_maxmaxprofiles_unbalanced_low_confidence_coordinates_large1.tsv",sep="\t").set_index("Unnamed: 0")
	yeast_predictions_annot=yeast_linreg0p1_predictions.copy()
	da_yeast.make_graph_data_among_classes(method=method,predictions_raw=yeast_predictions_annot,classifier=classifier_svmlin_reg0p1)
	
	#color1=metabolites
	#color2=proteins
	#normal
	da_yeast.make_EC_graph(layout=net_layout_large,method=method,approach=approach_yeast,classifier=classifier_svmlin_reg0p1,feature=feature,color1=tableau_blue,color2=tableau_orange,font_size=18,plot_all_vertices=True,trim="_trimmed")
	
	#save colormap bars
	colormap1=LinearSegmentedColormap.from_list("tabblue",[tableau_blue,tableau_grey],N=100)
	colormap2=LinearSegmentedColormap.from_list("tabgreen",[tableau_green,tableau_grey],N=100)
	colormap3=LinearSegmentedColormap.from_list("taborange",[tableau_orange,tableau_grey],N=100)
	colormap4=LinearSegmentedColormap.from_list("tabred",[tableau_red,tableau_grey],N=100)
	num_metas=len(yeast_linreg0p1_predictions.index.get_level_values(1).unique())
	num_prots=len(yeast_linreg0p1_predictions.index.get_level_values(0).unique())
	da_yeast.save_cbar(colormap1,"blue_to_grey",num_metas)
	da_yeast.save_cbar(colormap3,"orange_to_grey",num_prots)
	da_yeast.save_cbar(colormap2,"green_to_grey",num_prots)
	da_yeast.save_cbar(colormap4,"red_to_grey",num_prots)
	
	
	#heatmap brite classes - metabolite classes
	da_yeast.heatmap_bc_mc(method=method,trim="_trimmed",square=True)
	
	#heatmap kegg pathways
	da_yeast.heatmap_pathways(method=method,classifier=classifier_svmlin_reg0p1,trim="_trimmed",square=True)
	
	
	###################
	headers=["A Presence of interaction","B Absence of interaction"]
	da_yeast.plot_headers(headers, plot_with_letters=False)
	
	####################
	#suppfig: learning, ROC curves and degree distributions
	####################
	colours=["tab:green","tab:green","tab:orange","tab:orange","tab:red"]
	headers=["A. thaliana","S. cerevisiae","Both","OG-wise","protein-wise","A learning curves", "B ROC curves", "C Degree distributions"]
	da_yeast.plot_headers(headers, plot_with_letters=False)
	
	da_yeast.classifier_comparison(experiments=[appendix_ara,appendix_ara_og,appendix_yeast,appendix_yeast_og,appendix_euc],features=[feature],approaches=[approach_yeast],classifiers=[classifier_svmlin_reg0p1],methods=[method],appendix="_diff_experiments_svmlinreg0p1_min",legend_roc=True,legend_dd=False,colour=colours)
	'''
