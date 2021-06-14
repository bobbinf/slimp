#!/usr/bin/python3
import computer_scripts

if __name__=='__main__':
	###################################################
	#Programmatic workflow for one analysis experiment (e.g. S.cerevisiae, approximate duration: 25 hours)
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
