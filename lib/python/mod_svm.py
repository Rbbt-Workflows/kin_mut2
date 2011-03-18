#!/usr/bin/python

import sys, re, os, string, urllib, time, math, random, tempfile

		
from mod_utils import div_sec,isNumber,frange

from math import exp,log

import shutil

from multiprocessing import Pool
from multiprocessing import cpu_count

from subprocess import Popen , PIPE

from numpy import array,zeros,dot,flipud
from scipy import reshape, sqrt, identity
from numpy.matlib import repmat, repeat
from numpy import median as numpy_median 
from numpy import sum as sum_numpy
from numpy import exp as exp_numpy

from operator import itemgetter
from itertools import product
from itertools import groupby
from itertools import chain

from mod_categorical_data_analysis import get_log_likelihood
from mod_categorical_data_analysis import p_value
from mod_utils import log_sec,div_sec

min_increment = 5e-03

###############################################################

def train_test_model_selection_grid(svm):
	
	svm.train()
	
	svm.test(skip_performance_eval=True)
		
	return svm

#############################################################

def parse_pattern_file(filename):
		
		l_pattern = []
		l_labels  = []
		
		if filename == None:
			raise Exception('ERROR','mod_svm.parse_pattern_file: The file to parse has not been provided')
		
		if not os.path.exists(filename):
			raise Exception('ERROR','mod_svm.parse_pattern_file: The file to parse does not exist')
		
		fpat = file(filename,'r')
		l_lines = map(lambda l: filter(lambda x: x[0]<>'#', l.split('\t')), map(string.strip,fpat.readlines()))
		fpat.close()
		
		l_labels  = map(int,map(itemgetter(0),l_lines))
		l_pattern = map(lambda line: line[1:], l_lines)
		
		return l_pattern,l_labels
	
#############################################################

def parse_pattern_file_v2(filename):
		
		l_pattern = []
		l_labels  = []
		
		if filename == None:
			raise Exception('ERROR','mod_svm.parse_pattern_file: The file to parse has not been provided')
		
		if not os.path.exists(filename):
			raise Exception('ERROR','mod_svm.parse_pattern_file: The file to parse does not exist')
		
		fpat = file(filename,'r')
		l_lines = map(lambda l: l.split('\t'), map(string.strip,fpat.readlines()))
		fpat.close()
		
		l_labels  = map(int,map(itemgetter(0),l_lines))
		l_pattern = map(lambda line: filter(lambda x: x[0]<>'#', line[1:]), l_lines)
		l_comment = list(chain(*map(lambda line: filter(lambda x: x[0]=='#', line[1:]), l_lines)))
		
		return l_pattern,l_labels,l_comment

#############################################################

class c_feature:
	
	def __init__(self,feature=None,index=0):
		
		self.feature = feature
		self.index   = index
		self.no_positives = 0
		self.no_negatives = 0
		self.per_pos      = 0
		self.per_neg      = 0
		self.positive_distb = []
		self.negative_distb = []
		self.l_positive     = []
		self.l_negative     = []
		self.significance   = 0
		self.rnk_assoc_target = -1
		self.rnk_svm = -1
		
	def set_feature(self,feat):
		
		self.feature = feat
		
	def get_feature(self):
		
		return self.feature
	
	def set_index(self,i):
		
		self.index = i

	def get_index(self):
		
		return self.index
	
	def add_count_positive(self,score):
		
		if score == 0:
			return 
		
		self.no_positives += 1
		
	def add_count_negative(self,score):
		
		if score == 0:
			return
		
		self.no_negatives += 1
	
	def add_positives(self,score):
		
		if score == 0:
			return 
		
		self.l_positive.append(score)
		
	def add_negatives(self,score):
		
		if score == 0:
			return
		
		self.l_negative.append(score)
		
	def set_positives_negatives_percent(self,per_pos,per_neg):
		
		self.per_pos = per_pos
		self.per_neg = per_neg
		
	def get_positive_percent(self):
		
		return self.per_pos
	
	def get_negative_percent(self):
		
		return self.per_neg
		
	def get_number_positives(self):
		
		return self.no_positives
	
	def get_number_negatives(self):
		
		return self.no_negatives
	
	def __value2bin(self,value,bins):
		
		l_lower_bound = filter(lambda (i,b): b >= value, enumerate(bins))

		if l_lower_bound == []:
			raise Exception('ERROR','mod_svm.c_feature: __value2bin: The value %1.2f is not in the range of the bins' % (value))

		return l_lower_bound[0][0]
	
	def calculate_distributions(self,Npos,Nneg,bins=None):
		
		if bins == None:
			self.positive_distb = [(k,float(len(list(g)))/Npos) for k,g in groupby(sorted(self.l_positive))]
			self.negative_distb = [(k,float(len(list(g)))/Nneg) for k,g in groupby(sorted(self.l_negative))]
		else:
			self.positive_distb = [(k,float(len(list(g)))/Npos) for k,g in groupby(sorted(self.l_positive,key=lambda x: self.__value2bin(x, bins)))]
			self.negative_distb = [(k,float(len(list(g)))/Nneg) for k,g in groupby(sorted(self.l_negative,key=lambda x: self.__value2bin(x, bins)))]

	def get_positive_distribution(self):
		
		return self.positive_distb
	
	def get_negative_distribution(self):
		
		return self.negative_distb
	
	def calculate_significance(self,N_pos,N_neg):
		
		"""
            The contingency table is:
            
                            Positive_set       Negative_set
                            ----------------------------
            feature           n11=no_pos       n12=no_neg
            not_feature       n21=Npos-no_pos  n22=Nneg-no_neg
            
            no_pos = Number of items of the positive set where the feature appears.
            no_neg = Number of items of the negative set where the feature appears.
            Npos   = Total number of items in the positive set
            Nneg   = Total number of items in the negative set
        """
		
		n11 = self.no_positives
		n12 = self.no_negatives
		n21 = N_pos-n11
		n22 = N_neg-n12 

		contingency_table = (n11,n12,n21,n22)
		
		self.significance = p_value(get_log_likelihood(contingency_table))  
	
	def calculate_log_likelihood(self,N_pos,N_neg):
		
		"""
			The calculus is taken from the paper:
				"Comparing Corpora using Fequency Profiling". Paul Rayson and Roger Garside.
				WCC '00 Proceedings of the workshop on Comparing corpora - Volume 9. 2000
		
            The contingency table is:
            
                            Positive_set       Negative_set
                            ----------------------------
            feature           n11=no_pos       n12=no_neg
            not_feature       n21=Npos-no_pos  n22=Nneg-no_neg
            
            no_pos = Number of items of the positive set where the feature appears.
            no_neg = Number of items of the negative set where the feature appears.
            Npos   = Total number of items in the positive set
            Nneg   = Total number of items in the negative set
            
            The log-likelihood (LL) measures the relative frequency difference between the positive and negative
            sets. The higher the value the more significative the difference is.  
            
            On-line calculator: http://ucrel.lancs.ac.uk/llwizard.html
             
        """
		
		n11 = float(self.no_positives)
		n12 = float(self.no_negatives)
		n21 = N_pos-n11
		n22 = N_neg-n12
		
		coeff = div_sec((n11+n12),(N_pos+N_neg))
		
		E1 = N_pos*coeff
		E2 = N_neg*coeff 
		
		try:
			LL = 2*(n11*log_sec(div_sec(n11,E1))+n12*log_sec(div_sec(n12,E2)))
		except:
			print "aqui"
		
		self.significance = LL
	
	def set_significance(self,sig):
		
		self.significance = sig
	
	def get_significance(self):
		
		return self.significance
	
	def set_rank_association_target(self,rnk_index):
		
		self.rnk_assoc_target = rnk_index
		
	def get_rank_association_target(self):
		
		return self.rnk_assoc_target
	
	def set_rank_svm(self,rnk_index):
		
		self.rnk_svm = rnk_index

	def get_rank_svm(self):
		
		return self.rnk_svm

#############################################################

class c_result:
	
	def __init__(self,*arg):
		
		__class__ = "c_result"
		
		self.confusion_matrix = {}
		self.accuracy  = 0
		self.precision = 0
		self.recall    = 0
		self.fscore    = 0
		self.roc_area  = 0
		self.roc_values = []
		
		try:
		
			if isinstance(arg,tuple):
			
				if isinstance(arg[0],c_result):
				
					self = arg[0]
		
		except IndexError:
			
			pass
		
	def get_accuracy(self):
		
		return self.accuracy
	
	def set_accuracy(self,ac):
		
		self.accuracy = ac
		
	def get_precision(self):
		
		return self.precision
	
	def set_precision(self,prec):
		
		self.precision = prec
		
	def get_recall(self):
		
		return self.recall
		
	def set_recall(self,rec):
		
		self.recall = rec
		
	def get_fscore(self):
		
		return self.fscore
	
	def set_fscore(self,fscore):
		
		self.fscore = fscore
		
	def get_confusion_matrix(self):
		
		return self.confusion_matrix
	
	def set_confusion_matrix(self,conf_m):
		
		self.confusion_matrix = conf_m
		
	def get_roc_area(self):
		
		return self.roc_area
	
	def set_roc_area(self,roc_a):
		
		self.roc_area = roc_a
		
	def set_roc_values(self,roc_v):
		
		self.roc_values = roc_v
		
	def get_roc_values(self):
		
		return self.roc_values
		
	def __repr__(self):
		
		return "TP:%d FP:%d FN:%d TN:%d\tacc:%1.2f%% prec:%1.2f%% rec:%1.2f%% fscore:%1.2f%% roc_a:%1.4f\n" % (self.confusion_matrix.get('TP',0),self.confusion_matrix.get('FP',0),self.confusion_matrix.get('FN',0),self.confusion_matrix.get('TN',0),self.accuracy,self.precision,self.recall,self.fscore,self.roc_area)
	
	def sum(self,other):
		
		if not isinstance(other,c_result):
			
			raise self.__class__, "The object to be added is not a c_result instance"
		
		self.accuracy  = self.accuracy  + other.accuracy
		self.precision = self.precision + other.precision
		self.recall    = self.recall    + other.recall
		self.fscore    = self.fscore    + other.fscore
		self.roc_area  = self.roc_area  + other.roc_area
		self.confusion_matrix['TP'] = self.confusion_matrix.get('TP',0) + other.confusion_matrix.get('TP',0)
		self.confusion_matrix['TN'] = self.confusion_matrix.get('TN',0) + other.confusion_matrix.get('TN',0)
		self.confusion_matrix['FP'] = self.confusion_matrix.get('FP',0) + other.confusion_matrix.get('FP',0)
		self.confusion_matrix['FN'] = self.confusion_matrix.get('FN',0) + other.confusion_matrix.get('FN',0)
		
	def div(self,N):
		
		if not isNumber(N):
			raise Exception("ERROR","mod_svm.%s.div: The object to divide is not a number" % (self.__class_)) 
			
		if N == 0:
			raise Exception("ERROR", "mod_svm.%s.div: The number to divide is zero" % (self.__class_))
		
		self.accuracy  = self.accuracy/N
		self.precision = self.precision/N
		self.recall    = self.recall/N
		self.fscore    = self.fscore/N
		self.roc_area  = self.roc_area/N
		self.confusion_matrix['TP'] = self.confusion_matrix['TP']/N
		self.confusion_matrix['TN'] = self.confusion_matrix['TN']/N
		self.confusion_matrix['FP'] = self.confusion_matrix['FP']/N
		self.confusion_matrix['FN'] = self.confusion_matrix['FN']/N
		
	def cmp_roc_fscore(self,other):
		
		selfKey  = (float(self.get_roc_area()),float(self.get_fscore()))
		otherKey = (float(other.get_roc_area()),float(other.get_fscore()))
		
		return cmp(selfKey,otherKey)
	
	def cmp_roc_a(self,other):
		
		selfKey  = self.get_roc_area()
		otherKey = other.get_roc_area()
		
		return cmp(selfKey,otherKey)
	
	def cmp_fscore(self,other):
		
		selfKey  = self.get_fscore()
		otherKey = other.get_fscore()
		
		return cmp(selfKey,otherKey)
	
	def cmp_acc(self,other):
		
		selfKey  = self.get_accuracy()
		otherKey = other.get_accuracy()
		
		return cmp(selfKey,otherKey)
	
		
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

def get_confusion_matrix(args):
	
	if len(args) == 2: 
		l_scores = args[0]
		l_labels = args[1]
	else:
		l_scores = args[0]
		l_labels = args[1]	
		margin_threshold = args[2]
	
	l_tuples_label_score = map(lambda x,y: (x,y), l_labels,l_scores)
		
	l_P = filter(lambda x:  x[0] > 0, l_tuples_label_score)
	l_N = filter(lambda x:  x[0] <= 0, l_tuples_label_score)
		
	P = len(l_P)
	N = len(l_N)
		
	if P+N <> len(l_labels):
		raise Exception('ERROR','mod_svm.get_confusion_matrix: Problems in the binary labels')
		
	TP = len(filter(lambda x: x[1] > margin_threshold, l_P))
	TN = len(filter(lambda x: x[1] < margin_threshold, l_N))
	FN = P-TP
	FP = N-TN
		
	conf_matrix = {'TP':TP,'TN':TN,'FP':FP,'FN':FN,'P':P,'N':N}
		
	return conf_matrix

###########################################################################
				
class c_svm_light():
	
	def __init__(self,path_svm_light,kernel=None):
				
		self.beta = 1
		self.beta_pow_2 = math.pow(self.beta, 2)
		self.kernel = None
		if kernel <> None:
			self.configure_kernel(kernel)
		
		self.l_params = []
		self.dict_params = {} 
		self.modelfilename  = None
		self.train_filename = None
		self.test_filename  = None
		self.alpha_filename = None
		self.log_filename   = None
		self.prediction_filename = None
								
		self.path_svm_light = path_svm_light
		self.path_output_files = ""
					
		self.__support_vectors = None
		self.__alpha_Y         = None
		
	#######################################################
		
	def train(self,filename_training=None,model_filename=None,alpha_filename=None):
		
		if filename_training == None and model_filename == None and alpha_filename == None:
			filename_training = self.train_filename
			model_filename    = self.modelfilename
			alpha_filename    = self.alpha_filename
			
			if filename_training == None and model_filename == None and alpha_filename == None:
				raise Exception('ERROR','mod_svm.c_svm_light.train: The files needed for training had not been provided')
		
		model_file = os.path.join(self.path_output_files,model_filename)
		alpha_file = os.path.join(self.path_output_files,alpha_filename)
		
		svm_learn = os.path.join(self.path_svm_light,"svm_learn")
		
		l_call = self.__get_svm_light_call(self.kernel,self.l_params)

		train_sal = Popen([svm_learn]+l_call+["-m",'100',"-a",alpha_file,filename_training,model_file], stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata) = train_sal.communicate()
		train_sal.wait()
		
		if logdata <> '':
			raise Exception('ERROR','mod_svm.c_svm_light.train: %s' % (logdata))
		
		f_r = file(os.path.join(self.path_output_files,"svm_train.log"),'w')
		f_r.write(output)
		f_r.flush()
		f_r.close()
		
		return output

	#######################################################
	
	def train_K_fold(self,hash_filename_training=None,model_filename=None,alpha_filename=None):
		
		"""
			In K-fold training the input is not a single pattern file but a set of K files.
			The pattern file is splitted into K pieces and the training is composed by K sub-trainings. 
			In each round the Kth piece is retained and the training is performed with the remaining pieces. The Kth piece is used to validate the model
			
			The hash structure that is needed as input is composed by K keys each of one corresponds to a differerent partition.
			The entry with k key is composed by a tuple. The first element is a file with the training patterns and the other with the validation pattern
		"""
		
		if hash_filename_training==None and model_filename==None and alpha_filename==None:
			hash_filename_training = self.train_filename
			model_filename    = self.modelfilename
			alpha_filename    = self.alpha_filename
			
			if hash_filename_training == None and model_filename == None and alpha_filename == None:
				raise Exception('ERROR','mod_svm.c_svm_light.train_K_fold: The files needed for training had not been provided')
		
		l_folds = sorted(hash_filename_training.keys()) 
		
		if filter(lambda k: not isNumber(k), l_folds) <> []:
			raise Exception('ERROR','c_svm_light.train_K_fold: The input hash for training is not properly built')
		
		(root,ext) = os.path.splitext(model_filename)
		
		log_filename        = root + '.log'
		
		l_results        = []
		l_results_fscore = []
		
		for k in l_folds:
			
			(filename_training,filename_validation) = hash_filename_training[k]
			
			model_filename_k      = root + '.%d.svm' % (k)
			filename_prediction_k = root + '.%d.out' % (k)
			
			self.train(filename_training, model_filename_k, alpha_filename)
			
			val_result = self.test(model_filename_k, filename_validation, filename_prediction_k, log_filename)
			
			l_results_fscore.append((k,val_result.get_fscore()))
			l_results.append((k,val_result))
			
		median_fscore = numpy_median(map(itemgetter(1),l_results_fscore))
		
		k_median = filter(lambda (k,fs): fs==median_fscore, l_results_fscore)[0][0]
		
		result = filter(lambda (k,res): k==k_median, l_results)[0][1]
		
		(filename_training,filename_validation) = hash_filename_training[k_median]
		
		self.train(filename_training, model_filename, alpha_filename)
		
		filename_prediction          =  root + '.out'
		filename_prediction_k_median =  root + '.%d.out' % (k_median)
		
		shutil.copyfile(filename_prediction_k_median, filename_prediction)
						
		return result 		
	
	#######################################################
	
	def model_selection_grid_linear(self,filename_training,model_filename,alpha_filename,l_c,filename_valid=None):
		
		if filename_valid == None:
			filename_valid = filename_training
			
		(root,ext) = os.path.splitext(model_filename)
		
		log_filename = root + '.log'
		
		l_svm = []
		
		path_svm_light = self.get_path_svm_light()
		kernel         = self.get_kernel()
		
		for c in l_c:
			
			svm = c_svm_light(path_svm_light,kernel)
			
			model_filename      = root + '%1.2f.svm' % (c)
			filename_prediction = root + '%1.2f.out' % (c)
			
			svm.configure_train_filename(filename_training)
			svm.configure_test_filename(filename_valid)
			svm.configure_model(model_filename)
			svm.configure_alpha_filename(alpha_filename)
			svm.configure_log_filename(log_filename)
			svm.configure_prediction_filename(filename_prediction)
			
			svm.configure_params([c])
			
			l_svm.append(svm)
		
		pool = Pool(cpu_count())
		
		results = pool.map_async(train_test_model_selection_grid,l_svm)
		
		results.wait()
		
		l_svm = results.get()
		
		l_results = []
		
		l_pattern,l_labels = parse_pattern_file(filename_valid)
		
		os.sys.stdout.write('\nModel Selection for linear kernel\n\nsoft.margin\tgamma\tresults\n')
		
		for svm in l_svm:
			
			[c] = svm.get_params()
			
			filename_prediction = root + '%1.2f.out' % (c)
			
			result = svm.evaluate_performance(filename_prediction,l_labels)
			
			l_results.append(((c,),result))
			
			os.sys.stdout.write('%1.1f\t%s' % (c,result.__repr__()))
			os.sys.stdout.flush()
			
		return l_results	
	
	#######################################################
	
	def model_selection_grid_linear_Kfold(self,hash_filename_training,model_filename,alpha_filename,l_c):
			
		(root,ext) = os.path.splitext(model_filename)
		
		log_filename = root + '.log'
		
		l_svm = []
		
		path_svm_light = self.get_path_svm_light()
		kernel         = self.get_kernel()
		
		K = len(hash_filename_training.keys())
		
		for c in l_c:
			
			svm = c_svm_light(path_svm_light,kernel)
			
			model_filename      = root + '%1.2f.svm' % (c)
			filename_prediction = root + '%1.2f.out' % (c)
			
			svm.configure_train_filename(hash_filename_training)
			svm.configure_model(model_filename)
			svm.configure_alpha_filename(alpha_filename)
			svm.configure_log_filename(log_filename)
			svm.configure_prediction_filename(filename_prediction)
			
			svm.configure_params([c])
			
			l_svm.append(svm)
		
		l_results = []
		
		os.sys.stdout.write('\nModel Selection for linear kernel and K-fold cross-validation (%d folds)\n\nsoft.margin\tgamma\tresults\n' % (K))
		
		for svm in l_svm:
			
			[c] = svm.get_params()
			
			filename_prediction = root + '%1.2f.out' % (c)
			
			result = svm.train_K_fold()
			
			l_results.append(((c,),result))
			
			os.sys.stdout.write('%1.1f\t%s' % (c,result.__repr__()))
			os.sys.stdout.flush()
			
		return l_results
	
	#######################################################
	
	def __get_best_results(self,l_results,opt_criteria):

		if l_results == []:
			raise Exception('ERROR','c_svm_light.__get_best_results: The list of results is empty')
		
				
		if opt_criteria == 'max_acc':
			(best_params,best_result) = sorted(l_results,reverse=True,key=lambda (param,res): res.get_accuracy())[0]
		elif opt_criteria == 'max_fscore':
			(best_params,best_result) = sorted(l_results,reverse=True,key=lambda (param,res): res.get_fscore())[0]
		elif opt_criteria == 'max_roc_a':
			(best_params,best_result) = sorted(l_results,reverse=True,key=lambda (param,res): res.get_roc_area())[0]
					
		return best_params,best_result
		
	#######################################################
		
	def model_selection_grid_radial(self,filename_training,model_filename,alpha_filename,l_c,l_g,filename_valid=None):
		
		if filename_valid == None:
			filename_valid = filename_training
			
		(root,ext) = os.path.splitext(model_filename)
		
		log_filename = root + '.log'
				
		l_pairs_c_g = map(lambda pair_c_g: list(pair_c_g), product(l_c,l_g))
		
		l_svm = []
		
		path_svm_light = self.get_path_svm_light()
		kernel         = self.get_kernel()
		
		for pair_c_g in l_pairs_c_g:
			
			svm = c_svm_light(path_svm_light,kernel)
			
			model_filename      = root + '%1.2f_%1.4f.svm' % (pair_c_g[0],pair_c_g[1])
			filename_prediction = root + '%1.2f_%1.4f.out' % (pair_c_g[0],pair_c_g[1])
			
			svm.configure_train_filename(filename_training)
			svm.configure_test_filename(filename_valid)
			svm.configure_model(model_filename)
			svm.configure_alpha_filename(alpha_filename)
			svm.configure_log_filename(log_filename)
			svm.configure_prediction_filename(filename_prediction)
			
			svm.configure_params(pair_c_g)
			
			l_svm.append(svm)
		
		pool = Pool(cpu_count())
		
		results = pool.map_async(train_test_model_selection_grid,l_svm)
		
		results.wait()
		
		l_svm = results.get()
		
		l_results = []
		
		l_pattern,l_labels = parse_pattern_file(filename_valid)
		
		os.sys.stdout.write('\nModel Selection for radial kernel\n\nsoft.margin\tgamma\tresults\n')
		
		for svm in l_svm:
			
			(c,g) = svm.get_params()
			
			filename_prediction = root + '%1.2f_%1.4f.out' % (c,g)
			
			result = svm.evaluate_performance(filename_prediction,l_labels)
			
			l_results.append(((c,g),result))
			
			os.sys.stdout.write('%1.2f\t%1.4f\t%s' % (c,g,result.__repr__()))
			os.sys.stdout.flush()
			
		return l_results
	
	#######################################################
	
	def model_selection_grid_radial_Kfold(self,hash_filename_training,model_filename,alpha_filename,l_c,l_g):
			
		(root,ext) = os.path.splitext(model_filename)
		
		log_filename = root + '.log'
				
		l_pairs_c_g = map(lambda pair_c_g: list(pair_c_g), product(l_c,l_g))
		
		l_svm = []
		
		path_svm_light = self.get_path_svm_light()
		kernel         = self.get_kernel()
		K              = len(hash_filename_training.keys())
		
		for pair_c_g in l_pairs_c_g:
			
			svm = c_svm_light(path_svm_light,kernel)
			
			model_filename      = root + '.%1.2f_%1.4f.svm' % (pair_c_g[0],pair_c_g[1])
			filename_prediction = root + '.%1.2f_%1.4f.out' % (pair_c_g[0],pair_c_g[1])
			
			svm.configure_train_filename(hash_filename_training)
			svm.configure_model(model_filename)
			svm.configure_alpha_filename(alpha_filename)
			svm.configure_log_filename(log_filename)
			svm.configure_prediction_filename(filename_prediction)
			
			svm.configure_params(pair_c_g)
			
			l_svm.append(svm)
		
		l_results = []
		
		os.sys.stdout.write('\nModel Selection for radial kernel and K-fold cross-validation (%d folds)\n\nsoft.margin\tgamma\tresults\n' % (K))
		
		for svm in l_svm:
			
			(c,g) = svm.get_params()
						
			result = svm.train_K_fold()
			
			l_results.append(((c,g),result))
			
			os.sys.stdout.write('%1.2f\t%1.4f\t%s' % (c,g,result.__repr__()))
			os.sys.stdout.flush()
			
		return l_results
	
	#######################################################
	
	def model_selection_grid(self,filename_training,model_filename,alpha_filename,hash_grid,**kwargs):
		
		#### Grids
		l_grid_soft_margin = hash_grid.get('soft_margin',[])
		l_grid_gamma       = hash_grid.get('gamma',[])
		
		if l_grid_soft_margin == [] and l_grid_gamma == []:
			raise Exception('ERROR','mod_svm.model_selection_grid: The list of parameters to be inspected had not been provided')
		elif l_grid_soft_margin == []:
			raise Exception('ERROR','mod_svm.model_selection_grid: The list of soft-margin parameters to be inspected had not been provided')
		elif l_grid_gamma == [] and self.kernel == 'radial':
			raise Exception('ERROR','mod_svm.model_selection_grid: The list of gamma parameters to be inspected had not been provided')
		
		#### Training parameters
		K = 0
		if kwargs.has_key('k_fold'):
			K = kwargs['k_fold']
		
		filename_valid = None
		if kwargs.has_key('filename_valid'):
			filename_valid = kwargs['filename_valid']
			
		opt_criteria = 'max_roc_a'
		if kwargs.has_key('opt_criteria'):
			opt_criteria = kwargs['opt_criteria']
		
		if K == None or K == 0 or K == 1: # No K-fold cross-validation
			
			if self.kernel == 'linear':
				
				l_results_grid = self.model_selection_grid_linear(filename_training, model_filename, alpha_filename, l_grid_soft_margin, filename_valid)
				
			elif self.kernel == 'radial':
								
				l_results_grid = self.model_selection_grid_radial(filename_training, model_filename, alpha_filename, l_grid_soft_margin, l_grid_gamma, filename_valid)
				
			(best_params,best_result) = self.__get_best_results(l_results_grid, opt_criteria)
		else:
			
			if not isinstance(filename_training,dict):
				raise Exception('ERROR','c_svm_light.model_selection_grid: The input for training with K-folds must be a hash')
			
			if self.kernel == 'linear':
				
				l_results_grid = self.model_selection_grid_linear_Kfold(filename_training, model_filename, alpha_filename, l_grid_soft_margin)
				
			elif self.kernel == 'radial':
								
				l_results_grid = self.model_selection_grid_radial_Kfold(filename_training, model_filename, alpha_filename, l_grid_soft_margin, l_grid_gamma)
				
			(best_params,best_result) = self.__get_best_results(l_results_grid, opt_criteria)
					
		if self.kernel  == 'linear':
			os.sys.stdout.write('\n%s\nOptimization criteria: %s\nBest parameter: c=%1.2f\nBest result: %s%s\n\n' % ('*'*50,opt_criteria,best_params[0],best_result.__repr__(),'*'*50))
		elif self.kernel == 'radial':
			os.sys.stdout.write('\n%s\nOptimization criteria: %s\nBest parameter: c=%1.2f\tg=%1.4f\nBest result: %s%s\n\n' % ('*'*50,opt_criteria,best_params[0],best_params[1],best_result.__repr__(),'*'*50))
		
		os.sys.stdout.flush()
		
		return best_params,best_result
	
	#######################################################
	
	def __enrich_prediction_file(self,filename_prediction,l_labels,l_patterns):
		
		f = open(filename_prediction,'r')
		l_out = map(string.strip,f.readlines())
		f.close()
		
		f = open(filename_prediction,'w')
		f.write('\n'.join(map(lambda l,out,pat: "%s\t%s\t%s" % (l,out,pat), l_labels, l_out, l_patterns)))
		f.close()
	
	#######################################################
	
	def test(self,filename_model=None,filename_pattern_test=None,filename_prediction=None,log_filename=None,**kwargs):
		
		if filename_model == None: # This is when the model is configured once and it can be tested with different patterns 
			
			filename_model = self.modelfilename
			
			if filename_model == None:
				raise Exception('ERROR','mod_svm.c_svm_light.test: The model filename has not been provided')
			
		if filename_pattern_test == None:
			
			filename_pattern_test = self.test_filename
			
			if filename_pattern_test == None:
				raise Exception('ERROR','mod_svm.c_svm_light.test: The test pattern filename has not been provided')
			
		if filename_prediction == None:
			
			filename_prediction = self.prediction_filename
			
			if filename_prediction == None:
				prediction_tmp_file = tempfile.NamedTemporaryFile("w")
				filename_prediction = prediction_tmp_file.name
		
		svm_classify = os.path.join(self.path_svm_light,"svm_classify")
		
		if log_filename == None:
			log_filename = self.log_filename
		
		if log_filename == None:
			test_sal = Popen([svm_classify,filename_pattern_test,filename_model, filename_prediction],stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		else:
			test_sal = Popen([svm_classify,filename_pattern_test,filename_model, filename_prediction, log_filename],stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
	
		(output,logdata) = test_sal.communicate()
		test_sal.wait()
		
		if logdata.lower() <> '':
			raise Exception('ERROR','mod_svm.c_svm_light.test: %s' % (logdata))
		
		if kwargs.has_key('skip_performance_eval'):
			if kwargs['skip_performance_eval'] == True:
				return None 
		
		l_labels = []
		
		if kwargs.has_key('labels'):
			l_labels = kwargs['labels']
		else:
			l_pattern,l_labels,l_comments = parse_pattern_file_v2(filename_pattern_test)
			
		if len(set(l_labels)) > 2:
			raise Exception('WARNING','mod_svm.c_svm_light.test: Multiclass classification problem is not implemented')
			
		result = self.evaluate_performance(filename_prediction,l_labels)
		
		self.__enrich_prediction_file(filename_prediction,l_labels,l_comments)
				
		return result
						
	#######################################################					
	
	def configure_beta(self,beta):
		
		self.beta = beta
		self.beta_pow_2 =  math.pow(self.beta, 2)
		
	def configure_params(self,l_params):
		
		self.l_params = l_params
		
	def get_params(self):
		
		return self.l_params
		
	def configure_model(self,svm_filename):
		
		self.modelfilename = svm_filename
		
	def configure_train_filename(self,train_filename):
		
		self.train_filename = train_filename
		
	def configure_test_filename(self,test_filename):
		
		self.test_filename = test_filename
		
	def configure_alpha_filename(self,alpha_filename):
		
		self.alpha_filename = alpha_filename
		
	def configure_log_filename(self,log_filename):
		
		self.log_filename = log_filename
		
	def configure_prediction_filename(self,pred_filename):
		
		self.prediction_filename = pred_filename 
		
	def configure_kernel(self,kernel):
		
		self.kernel = kernel
		
	def get_kernel(self):
		
		return self.kernel
						
	def get_path_svm_light(self):
				
		return self.path_svm_light
	
	def __get_svm_light_call(self,kernel,l_param):

		"""
			Kernel options:
			-t int      -> type of kernel function:
				#0: linear (default)
				#1: polynomial (s a*b+c)^d
				#2: radial basis function exp(-gamma ||a-b||^2)
				#3: sigmoid tanh(s a*b + c)
				#4: user defined kernel from kernel.h
			-d int      -> parameter d in polynomial kernel
			-g float    -> parameter gamma in rbf kernel
			-s float    -> parameter s in sigmoid/poly kernel
			-r float    -> parameter c in sigmoid/poly kernel
			-u string   -> parameter of user defined kernel
		"""
		
		l_call = []
		
		if kernel == 'radial':
					
			if len(l_param) <> 2:
		
				raise Exception('ERROR','c_svm_light.__get_svm_light_string_kernel: The parameter list for radial basis svm is empty or incomplete')
		
			l_call = ['-t','2','-c','%f' % (l_param[0]),'-g','%f' % (l_param[1])]
		
		elif kernel == 'linear':
		
			string_kernel = ' -t 0'
			
			if len(l_param) <> 1:
				
				raise Exception('ERROR','c_svm_light.__get_svm_light_string_kernel: The parameter list for linear basis svm is empty or incomplete')
			
			l_call = ['-t','0','-c','%f' % (l_param[0])]	
					
		elif kernel == 'polynomial':
		
			if len(l_param) < 4:
		
				raise Exception('ERROR','c_svm_light.__get_svm_light_string_kernel: The parameter list for polynomial svm is empty or incomplete')
			
			l_call = ['-t','1','-c','%f' % (l_param[0]),'-d','%d' % (l_param[1]),'-s','%f' % (l_param[2]),'-r','%f' % (l_param[3])]
		
		elif kernel == 'sigmoid':
		
			if len(l_param) < 3:
		
				raise Exception('ERROR','c_svm_light.__get_svm_light_string_kernel: The parameter list for sigmoid svm is empty or incomplete')

			l_call = ['-t','3','-c','%f' % (l_param[0]),'-s','%d' % (l_param[1]),'-r','%f' % (l_param[2])]
		
		return l_call
	
	def evaluate_performance(self,filename_out,l_labels):
		
		file_out = file(filename_out,'r')
		l_scores = map(float,map(string.strip,file_out.readlines()))
		file_out.close()
		
		conf_matrix = self.get_confusion_matrix(l_scores,l_labels)
			
		sensitivity = div_sec(conf_matrix['TP'],conf_matrix['P'])
		specificity = div_sec(conf_matrix['TN'],(conf_matrix['TN']+conf_matrix['TP'])) 
		acc         = 100*div_sec((conf_matrix['TP']+conf_matrix['TN']),(conf_matrix['P']+conf_matrix['N']))
		rec         = 100*sensitivity
		prec        = 100*div_sec(conf_matrix['TP'],(conf_matrix['TP']+conf_matrix['FP']))  
			
		l_fpr_tpr = self.calculate_roc_values(l_scores, l_labels)
		roc_area  = self.calculate_roc_area(l_fpr_tpr)
			
		fscore = self.get_f_score(acc, prec, rec)
	
		result = c_result()
		
		result.set_confusion_matrix(conf_matrix)
		result.set_accuracy(acc)
		result.set_precision(prec)
		result.set_recall(rec)
		result.set_fscore(fscore)
		result.set_roc_values(l_fpr_tpr)
		result.set_roc_area(roc_area)
			
		return result


	def get_p_n(self,l_labels):
		
		if len(set(l_labels)) > 2:
			
			raise Exception('ERROR','c_svm_light.get_p_n: Multiclass classification problem is not implemented')
		
		P = len(filter(lambda x:  x > 0, l_labels))
		N = len(filter(lambda x:  x < 0, l_labels))
		
		return {'P':P,'N':N}
	
	def get_confusion_matrix(self,l_scores,l_labels,margin_threshold=0):
				
		l_tuples_label_score = map(lambda x,y: (x,y), l_labels,l_scores)
		
		l_P = filter(lambda x:  x[0] > 0, l_tuples_label_score)
		l_N = filter(lambda x:  x[0] <= 0, l_tuples_label_score)
		
		P = len(l_P)
		N = len(l_N)
		
#		if P+N <> len(l_labels):
#			raise Exception('ERROR','c_svm_light.get_confusion_matrix: Problems in the binary labels')
		
		TP = len(filter(lambda x: x[1] > margin_threshold, l_P))
		TN = len(filter(lambda x: x[1] < margin_threshold, l_N))
		FN = P-TP
		FP = N-TN
		
		conf_matrix = {'TP':TP,'TN':TN,'FP':FP,'FN':FN,'P':P,'N':N}
		
		return conf_matrix
	
	def __get_TPR(self,TP,P):
		
		try:
			return float(TP)/P
		except ZeroDivisionError:
			return float(0)
		
	def __get_FPR(self,FP,N):
		
		try:
			return float(FP)/N
		except ZeroDivisionError:
			return float(0)
	
	def calculate_roc_values(self,l_scores,l_labels):
		
		"""
			Bibliography:
			ROC Graphs: Notes and Practical Considerations for Researchers
			Tom Fawcett (tom.fawcett@hp.com)
			HP Laboratories, MS 1143, 1501 Page Mill Road, Palo Alto, CA 94304
		"""
		
		if len(l_scores) <= 1:
			return []
				
		l_scores_sort = sorted(l_scores)
		mini = l_scores_sort[0]
		maxi = l_scores_sort[-1]
		increment  = max(l_scores_sort[1]-l_scores_sort[0],min_increment)
				
		l_margin_thres = frange(mini,maxi,increment) 
		l_margin_thres.append(1)
		
		num_cpu = 2
		
		if cpu_count() > 1:
			num_cpu = cpu_count() 
		
		pool = Pool(num_cpu-1)
		
		result = pool.map_async(get_confusion_matrix, map(lambda x: (l_scores,l_labels,x) , l_margin_thres))
		
		result.wait()
		
		l_results = result.get()
		
		l_roc_fpr = map(lambda x: self.__get_FPR(x['FP'],x['N']),l_results)
		l_roc_tpr = map(lambda x: self.__get_TPR(x['TP'],x['P']),l_results)
				
		return sorted(zip(l_roc_fpr,l_roc_tpr))
				
	def calculate_roc_area(self,l_fpr_tpr):
						
		area = 0
		
		if l_fpr_tpr == []:
			return area
		
		fpr_ini = l_fpr_tpr[0][0]
		tpr_ini = l_fpr_tpr[0][1]
		
		for i in range(1,len(l_fpr_tpr)):
			
			fpr_fin = l_fpr_tpr[i][0]
			tpr_fin = l_fpr_tpr[i][1]
			
			base      = fpr_fin-fpr_ini
			av_height = float(tpr_fin+tpr_ini)/2 
			
			area += base*av_height
			
			fpr_ini = fpr_fin
			tpr_ini = tpr_fin
							
		return area 			
			
		
	def get_performance(self,filename_test_log):

		file_test_log = file(filename_test_log,'r')

		l_lines = file_test_log.readlines()
		
		l_lines = map(lambda x: x.replace('%',''),l_lines)

		accuracy  = 0
		precision = 0
		recall    = 0

		for line in l_lines:

			if line.find("Accuracy on test set:") >= 0:

				accuracy = float(string.split(line[22:])[0])

			if line.find("Precision/recall on test set:") >= 0:

				precision = float(string.split(string.split(line[29:])[0],'/')[0])
				recall    = float(string.split(string.split(line[29:])[0],'/')[1])
		
		return (accuracy,precision,recall,self.get_f_score(accuracy,precision,recall))
		
	
		
#	def test(self,filename_model,filename_pattern_test,filename_prediction,l_labels,log_file=None):
#		
#		os.system(self.get_svm_light_string_test(filename_pattern_test, filename_model, filename_prediction, log_file))
#		
#		if l_labels <> []:
#		
#			result = self.evaluate_performance(filename_prediction,l_labels)
#		
#			return result
#		
#		return c_result()
			
	def get_f_score(self,acc,prec,rec):
	
		try:
			f_beta = float((1+self.beta_pow_2)*prec*rec)/((self.beta_pow_2*prec)+rec)
		except ZeroDivisionError:
			f_beta=float(0)

		return f_beta
	
	def __fillRow(self,i,d,Dist):
		
		#d = sum((Points[i] - Points[i+1:])**2,axis=1)
		Dist[i,i+1:] = d
		Dist[i+1:,i] = d
	
	def __calcDistanceMatrix(self,Points):
		
		"""
			Points is a list of arrays
		"""
		
		num_points = len(Points)
		
		dim_axis = Points.ndim-1
								
		l_d = map(lambda i: sum_numpy((Points[i] - Points[i+1:])**2,axis=dim_axis),range(num_points))
					
		Dist = zeros((num_points,num_points),float)
		
		map(lambda i,d: self.__fillRow(i,d,Dist),range(num_points),l_d)
			
		return Dist
	
	def calculate_gram_matrix_linear(self,X,gamma=None,degree=None):
		
		return dot(X,X.T)
	
	def calculate_gram_matrix_radial(self,X,gamma,degree=None):
		
		K = self.__calcDistanceMatrix(X)
		
		return exp_numpy(-gamma*K)
	
	def calculate_gram_matrix_poly(self,X,gamma,degree,coefficient):
		
		return (gamma * (dot(X,X.T)) + coefficient) ** degree
				
	def __parse_pattern_line(self,line):
	
		alpha = float(line[0])
	
		comment = line[-1]
				
		l_ = map(lambda x: x.split(':'),line[1:-1])
		
		l_weight = map(lambda x: (int(x[0]),float(x[1])),l_)
		
		return (alpha,l_weight,comment)

	def __assign_value(self,X,i,w):

		X[i-1] = w

	def __build_pattern(self,X,Nfeat):
		
		"""
			X is a list of tuples the tuples (feat,weight). The missing tuples are included with 0 weight 
		"""
		array = zeros(Nfeat,float)
		
		map(lambda (feat,w): self.__assign_value(array, feat, w), X)
		
		return array
		
	def load_model(self,filename,Nfeat_tot):
		
		#2 # kernel type
		#3 # kernel parameter -d 
		#0.2 # kernel parameter -g 
		#1 # kernel parameter -s 
		#1 # kernel parameter -r 
		#empty# kernel parameter -u 
		#9780 # highest feature index 
		#5588 # number of training documents 
		#1206  
		#0.25008607 # threshold b
		
		file_model = file(filename,"r")
	
		l_lines = file_model.readlines()
	
		num_sv = 0
		Nfeat  = 0
		kernel = ""
		l_kernel = ['linear','polynomial','radial'] 
				
		dict_params = {}
	
		for pos_line in range(1,len(l_lines)):
		
			line = l_lines[pos_line]
		
			if string.find(line,"# kernel type") > 0:
				kernel = l_kernel[int(line.split()[0])]
				pass
			elif string.find(line,"# kernel parameter -d") > 0 and kernel == 'polynomial':
				dict_params['gamma']  = 0
				dict_params['degree'] = float(line.split()[0]) 
			elif string.find(line,"# kernel parameter -g") > 0 and kernel <> 'linear':
				dict_params['gamma']  = float(line.split()[0])
				dict_params['degree'] = 0
			elif string.find(line,"# kernel parameter -s") > 0:
				pass
			elif string.find(line,"# kernel parameter -r") > 0:
				pass
			elif string.find(line,"# kernel parameter -u") > 0:
				pass
			elif string.find(line,"# highest feature index") > 0:
				Nfeat = int(line.split()[0])
				dict_params['nfeat'] = Nfeat
			elif string.find(line,"# number of training documents") > 0:
				pass
			elif string.find(line,"# number of support vectors plus 1") > 0:
				num_sv = int(string.split(line)[0])-1
			elif string.find(line,"# threshold b") > 0:
				dict_params['bias'] = float(line.split()[0])
				break
		
		if Nfeat > Nfeat_tot:
			raise Exception('ERROR','mod_svm.c_svm_light.load_model: There is an inconsistency in the model file: the number of features read from the model file overcomes the total number of features')
		
		l_alpha_patterns = map(string.split,l_lines[pos_line+1:])
		
		if len(l_alpha_patterns) <> num_sv:
			raise Exception("ERROR","mod_svm.c_svm_light.load_model: The number of support vectors is not correct")
	
		l_pattern = []
		l_alpha   = []
		
		l_patterns = map(lambda line: self.__parse_pattern_line(line),l_alpha_patterns)
	
		l_alpha_Y    = array(map(itemgetter(0),l_patterns))
		support_vect = array(map(lambda x: self.__build_pattern(x,Nfeat_tot), map(itemgetter(1),l_patterns)))
			
		return (l_alpha_Y,support_vect,kernel,dict_params)
	
	def rank_features(self,filename,l_features):
		
		"""
			The algorithm for ranking is taken from 'Variable Selection using SVM-based criteria'
			Alain Rakotomamonjy
			Journal of Machine Learning Research Vol. 3. pp. 1357-1370. 2004
		"""
		
		Nfeat = len(l_features)
		
		#l_pos_features = filter(lambda feat: feat.get_number_positives()>0, l_features)
		#Nfeat = len(l_pos_features)		
				
		(alpha_Y,supp_vect,kernel,dict_params) = self.load_model(filename,Nfeat)
		
		if kernel not in ['linear','radial','polynomial']:
			raise Exception("ERROR","mod_svm.c_svm_light.rank_features: The kernel %s is not correct" % (kernel))
		
		gamma  = dict_params.get('gamma',None)
		degree = dict_params.get('degree',None)
		
		if gamma == None or degree == None:
			raise Exception('ERROR','mod_svm.c_svm_light.rank_features: The parameters of the kernel had not been parsed correctly')
		
		if kernel == "linear":
			calculate_gram_matrix = self.calculate_gram_matrix_linear
		if kernel == "radial":
			calculate_gram_matrix = self.calculate_gram_matrix_radial
		if kernel == "polinomial":
			calculate_gram_matrix = self.calculate_gram_matrix_poly 
		
								
		#############################################################
		# Contribution of each feature i to the decision function:
		# DJ(i)=(1/2)(alpha_y.T*K*alpha_y - alpha_y.T*K(-i)*alpha_y)
		# K(-i) is K with the i feature removed
		#############################################################
		
		# Calculate Gram-Matrix
		K = calculate_gram_matrix(supp_vect,gamma,degree)
	
		# Calculate the first term in the sum (fixed term)
		DJ_1 = dot(dot(alpha_Y, K),alpha_Y)
				
		# Calculate the second term in the sum (the term varies according to the eliminated feature)
		DJ_2 = map(lambda i: dot(alpha_Y,dot(K-calculate_gram_matrix(supp_vect[:,i],gamma,degree),alpha_Y)),range(Nfeat))
						
		# Calculate the weight of each feature
		DJ = (0.5*(DJ_1 - DJ_2))

		weight_feat = DJ**2
	
		# Ranking considering the weight of each feature
		l_feat_weight = sorted(zip(l_features,weight_feat),key=itemgetter(1),reverse=True)
		
		map(lambda (feat,w),i: feat.set_rank_svm(i), l_feat_weight,range(len(l_feat_weight)))
						
		l_features = map(itemgetter(0),l_feat_weight)
		
		return l_features
	
	