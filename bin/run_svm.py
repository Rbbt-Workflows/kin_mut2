#!/usr/bin/python
'''
Created on Sep 16, 2010

@author: adelpozo
'''


from __future__ import division

import sys, re, os, string, urllib, time, math, random

from os import path as osp

root_dir = osp.dirname(osp.dirname(osp.realpath(__file__)))

modulesRelDirs = ["lib/python/"]

for moduleRelDir in modulesRelDirs:
        sys.path.insert(0,osp.join(root_dir,moduleRelDir))

model_dir = osp.join(root_dir,'share/model/')

import optparse

from ConfigParser import *

from math import log,pow

from stat import *

from operator import itemgetter
from itertools import groupby
from itertools import chain
from itertools import product as product_iter

from mod_utils import frange

import mod_svm

##########################################################################

class OptionParser(optparse.OptionParser):

    def check_required (self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)
        
        if atrib is None:
            self.error("%s option not supplied" % option)

##########################################################################

def read_config_file(cfg_filename):
    
    if not os.path.exists(cfg_filename):
        raise Exception("ERROR","run_svm.read_config_file: The file cfg were not found")
    
    config_parser = ConfigParser()
        
    config_parser.read(cfg_filename)
   
    hash_kernel = {}
        
    #########################
    if not config_parser.has_section('SVM_LIGHT'):
        raise Exception("ERROR","run_svm.read_config_file: The path of svm_light is not provided in the cfg file")
    if not config_parser.has_section('KERNEL'):
        raise Exception("ERROR","run_svm.read_config_file: The kernel information is not provided in the cfg file")
    if not config_parser.has_section('MODEL_SELECTION'):
        raise Exception("ERROR","run_svm.read_config_file: The parameters for model selection is not provided in the cfg file")
    
    path_svm_light = config_parser.get('SVM_LIGHT','path')
    
    if not os.path.exists(path_svm_light):
        raise Exception("ERROR","run_svm.read_config_file: The path of svm_light does not exist")

    hash_kernel['path_svm_light'] = path_svm_light
    
    kernel_type = config_parser.get('KERNEL','kernel')
    
    l_kernel = ['linear','polynomial','radial','sigmoid']
    
    if kernel_type not in l_kernel:
        raise Exception("ERROR","run_svm.read_config_file: The kernel is not properly configured. Options: %s" % (','.join(l_kernel)))
    
    hash_kernel['kernel_type'] = kernel_type
    
    try:
        soft_margin_param = config_parser.getfloat('KERNEL','soft_margin_param') # C parameter
    except: 
        soft_margin_param = 0
        
    try:
        gamma = config_parser.getfloat('KERNEL','gamma_rbf') # is the wide of the gaussian
    except:
        gamma = 0
    
    try:
        polynomial_exp = config_parser.getfloat('KERNEL','polynomial_exp') # d parameter. Exponent of the polynomial
    except:
        polynomial_exp = 0
    
    try:
        slope = config_parser.getfloat('KERNEL','slope') # s parameter. Slope of polynomial or sigmoid
    except:
        slope = 0
    
    try:
        cross_y = config_parser.getfloat('KERNEL','cross_y') # c parameter. Point where the function crosses the y axis in polynomial or sigmoid.
    except:
        cross_y = 0
        
    l_param = []
    
    if kernel_type == 'linear':
        l_param = [soft_margin_param]
    elif kernel_type == 'polynomial':
        l_param = [soft_margin_param,polynomial_exp,slope,cross_y]
    elif kernel_type == 'radial':
        l_param = [soft_margin_param,gamma]
    elif kernel_type == 'sigmoid':
        l_param = [soft_margin_param,slope,cross_y]
    
    hash_kernel['parameters'] = l_param
    
    K = config_parser.get('CROSS_VALIDATION','k_fold')
    
    if K == '':
        K = 0
    
    K = int(K)
    
    hash_kernel['k_fold'] = K
    
    #################################################################
    
    hash_kernel['model_selection'] = {}
    
    opt_criteria = config_parser.get('MODEL_SELECTION','opt_citeria')
    
    if opt_criteria <> '':
        hash_kernel['model_selection']['opt_criteria'] = opt_criteria 
     
        
    soft_margin1 = config_parser.get('MODEL_SELECTION','soft_margin1')
    
    if soft_margin1 == '':
        soft_margin1 = 0
    
    soft_margin1 = float(soft_margin1)
    
    soft_margin2 = config_parser.get('MODEL_SELECTION','soft_margin2')
    
    if soft_margin2 == '':
        soft_margin2 = 0
    
    soft_margin2 = float(soft_margin2)
    
    step_soft_margin = config_parser.get('MODEL_SELECTION','step_soft_margin')
    
    if step_soft_margin == '':
        step_soft_margin = 0
    
    step_soft_margin = float(step_soft_margin)
    
    gamma1 = config_parser.get('MODEL_SELECTION','gamma1')
    
    if gamma1 == '':
        gamma1 = 0
    
    gamma1 = float(gamma1)
    
    gamma2 = config_parser.get('MODEL_SELECTION','gamma2')
    
    if gamma2 == '':
        gamma2 = 0
    
    gamma2 = float(gamma2)
    
    step_gamma = config_parser.get('MODEL_SELECTION','step_gamma')
    
    if step_gamma == '':
        step_gamma = 0
    
    step_gamma = float(step_gamma)

    hash_kernel['model_selection']['soft_margin1'] = soft_margin1 
    hash_kernel['model_selection']['soft_margin2'] = soft_margin2
    hash_kernel['model_selection']['step_soft_margin'] = step_soft_margin
    hash_kernel['model_selection']['gamma1'] = gamma1 
    hash_kernel['model_selection']['gamma2'] = gamma2
    hash_kernel['model_selection']['step_gamma'] = step_gamma
    
    return hash_kernel 

##########################################################################

def generate_grid_search(hash_model_selection,kernel):

    c1     = hash_model_selection['soft_margin1']
    c2     = hash_model_selection['soft_margin2']
    step_c = hash_model_selection['step_soft_margin']
    
    if c1 >= c2:
        raise Exception('ERROR','run_svm.generate_grid_search: The limit points of soft_margin parameter are not properly configured')
    
    if step_c <= 0:
        raise Exception('ERROR','run_svm.generate_grid_search: The incremental amount for soft_margin scanning is not properly configured')
    
    grid_soft_margin = frange(c1,c2,step_c)
            
    if len(grid_soft_margin) <= 1:
        raise Exception('ERROR','run_svm.generate_grid_search: The list of soft_margin parameter has one element or is empty') 

    if kernel == 'linear':    
        
        return grid_soft_margin,[]
    
    g1     = hash_model_selection['gamma1']
    g2     = hash_model_selection['gamma2']
    step_g = hash_model_selection['step_gamma']
    
    if g1 >= g2:
        raise Exception('ERROR','run_svm.generate_grid_search: The limit points of gamma parameter are not properly configured')
    
    if step_g <= 0:
        raise Exception('ERROR','run_svm.generate_grid_search: The incremental amount for gamma scanning is not properly configured')
    
    grid_gamma = frange(g1,g2,step_g)
            
    if len(grid_gamma) <= 1:
        raise Exception('ERROR','run_svm.generate_grid_search: The list of gamma parameter has one element or is empty') 
             
    return grid_soft_margin,grid_gamma

##########################################################################

def __are_blocked(filename,output_path,K):
    
    (root,ext) = os.path.splitext(os.path.basename(filename))
    
    root = os.path.join(output_path,root)
    
    are_blocked = True
    
    for k in range(K):
        
        filename_train_k = root+ '.%d.train' % (k) + ext
        filename_valid_k = root+ '.%d.valid' % (k) + ext
        
        if not os.path.exists(filename_train_k) and not os.path.exists(filename_valid_k):
            are_blocked = False
            break
        else:
            if os.access(filename_train_k, os.W_OK) or os.access(filename_valid_k, os.W_OK):
                are_blocked = False
                break
            
    return are_blocked
    

def split_file_into_K_fold(filename,output_path,K):
    
    if __are_blocked(filename,output_path,K):
        
        hash_filenames_K = {}
        
        (root,ext) = os.path.splitext(os.path.basename(filename))
    
        root = os.path.join(output_path,root)
    
        for k in range(K):
        
            filename_train_k = root+ '.%d.train' % (k) + ext
            filename_valid_k = root+ '.%d.valid' % (k) + ext
        
            hash_filenames_K[k] = (filename_train_k,filename_valid_k)
            
        return hash_filenames_K
    
    l_pattern,l_labels = mod_svm.parse_pattern_file(filename)
    
    l_index_labels = zip(range(len(l_labels)),l_labels)
    
    neg_patterns = map(itemgetter(0),filter(lambda (i,lb): int(lb)<0, l_index_labels))
    pos_patterns = map(itemgetter(0),filter(lambda (i,lb): int(lb)>0, l_index_labels))
            
    map(lambda i: random.shuffle(pos_patterns),range(5))
    map(lambda i: random.shuffle(neg_patterns),range(5))
    
    size_K_pos = int(0.5+float(len(pos_patterns))/K)
    size_K_neg = int(0.5+float(len(neg_patterns))/K)
    
    f = open(filename,'r')
    l_lines = map(string.strip,f.readlines())
    f.close()
    
    hash_K_fold = {}
       
    for k in range(K-1):
        
        l_index_k = pos_patterns[size_K_pos*k:size_K_pos*(k+1)]   
        l_index_k.extend(neg_patterns[size_K_neg*k:size_K_neg*(k+1)])
        
        hash_K_fold[k] = map(lambda i: l_lines[i], l_index_k)
        
    k = K-1
     
    l_index_k = pos_patterns[size_K_pos*k:]   
    l_index_k.extend(neg_patterns[size_K_neg*k:])
       
    hash_K_fold[k] = map(lambda i: l_lines[i], l_index_k)
    
    (root,ext) = os.path.splitext(os.path.basename(filename))
    
    root = os.path.join(output_path,root)
    
    hash_filenames_K = {}
    
    for k in sorted(hash_K_fold.keys()):
        
        filename_train_k = root+ '.%d.train' % (k) + ext
        filename_valid_k = root+ '.%d.valid' % (k) + ext
        
        l_folds_train = list(set(range(K)).difference(set([k])))
        l_pat_train = []
        for kt in l_folds_train:
            l_pat_train.extend(hash_K_fold[kt])
        fk = open(filename_train_k,'w')
        fk.write('\n'.join(l_pat_train))
        fk.close()
        
        fk = open(filename_valid_k,'w')
        fk.write('\n'.join(hash_K_fold[k]))
        fk.close()
                
        hash_filenames_K[k] = (filename_train_k,filename_valid_k)
          
    return hash_filenames_K 
    
##########################################################################

def run(argv=None,**kwargs):
 
    if argv is None: argv = sys.argv    
       
    parser = OptionParser(add_help_option=True,description="SVM training and/or testing with svm_light code\nThe options with an (*) are obligatory\n")
    
    parser.add_option("--m",default=None,help="Performance Mode: t(training)/e(evaluation) (*)",dest="mode")
    parser.add_option("--tr",default=None,help="Pattern file with training examples (*)",dest="pattern_train_filename")
    parser.add_option("--ts",default=None,help="Pattern file with testing examples",dest="pattern_test_filename")
    parser.add_option("--o",default=None,help="Path of output files (*)",dest="output_path")
    parser.add_option("--svm",default=None,help="Model svm file. It is the support vector trained (*)",dest="model_filename")
    parser.add_option("--cfg",default=None,help="Configuration File (*)",dest="kernel_cfg")
    parser.add_option("--s",action="store_true",help="Model Selection Option",dest="model_selection")
    parser.add_option("--k",action="store_true",help="K-fold cross validation",dest="K_fold")
                        
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:

        sys.exit(0)
    
    parser.check_required("--m")
    parser.check_required("--svm")
    parser.check_required("--o")
    parser.check_required("--cfg")
    
    if options.mode == 't':
        parser.check_required("--tr")
    elif options.mode == 'e':
        parser.check_required("--ts")
    
    try:
        
        if not os.path.exists(options.output_path):
            os.mkdir(options.output_path)
            sys.stdout.write('\nThe path %s is being created\n\n')
            sys.stdout.flush()
        
        hash_kernel = read_config_file(options.kernel_cfg)
        
        SVM = mod_svm.c_svm_light(hash_kernel.get('path_svm_light'),hash_kernel.get('kernel_type'))
        
        SVM.configure_params(hash_kernel['parameters'])
        
        if options.mode == 't' and not options.model_selection: # Train without model selection
            
            input_file = options.pattern_train_filename
            model_file = os.path.join(model_dir,options.model_filename)
            (root,ext) = os.path.splitext(options.model_filename)
            alpha_file = os.path.join(options.output_path,root+'.alpha')
            log_file   = os.path.join(options.output_path,root+'.log')
            prediction_file = os.path.join(options.output_path,root+'.results')
            
            t1 = time.time()
            
            if options.K_fold:
                
                K = hash_kernel['k_fold']
                
                if K <= 1:
                    raise Exception('ERROR','run_svm.split_file_into_K_fold: The number of folds must be greater than or equal to 2')
                
                sys.stdout.write("\n\nCross-validation with K-fold\n\n")
                sys.stdout.write('Spitting the file %s into %d folds\n\nTraining... ' % (input_file,K))
                sys.stdout.flush()
                
                hash_filenames_K = split_file_into_K_fold(input_file,options.output_path,K)
            
                train_result = SVM.train_K_fold(hash_filenames_K, model_file, alpha_file)
                
                t2 = time.time()
                sys.stdout.write("done [%0.3fs].\n\n" % (t2-t1))
                sys.stdout.write("K-fold result:\n%s\n%s\n" % ('*'*50,train_result.__repr__()))
                sys.stdout.flush()
            
            else:
                
                sys.stdout.write('\nLearning the features from file %s... ' % (input_file))
                sys.stdout.flush()
                
                train_result = SVM.train(input_file,model_file,alpha_file)
            
                t2 = time.time()
                sys.stdout.write("done [%0.3fs].\n\n" % (t2-t1))
                sys.stdout.write("Train result: %s\n" % (train_result.__repr__()))
                sys.stdout.flush()
            
                if options.pattern_test_filename <> None:
                
                    input_file = os.path.join(options.output_path,options.pattern_test_filename)
            
                    t1 = time.time()
                    sys.stdout.write('\nTesting the model %s with the file %s... ' % (model_file,input_file))
                    test_result = SVM.test(model_file, input_file, prediction_file, log_file)
                    t2 = time.time()
                    sys.stdout.write("done [%0.3fs].\n\n" % (t2-t1))
                    sys.stdout.write('Test result:\n%s\n%s\n' % ('*'*50,test_result.__repr__()))
                    sys.stdout.flush()
        
        elif options.mode == 't' and options.model_selection: # Train with model selection
            
            input_file = options.pattern_train_filename
            model_file = os.path.join(model_dir,options.model_filename)
            (root,ext) = os.path.splitext(options.model_filename)
            alpha_file = os.path.join(options.output_path,root+'.alpha')
            log_file   = os.path.join(options.output_path,root+'.log')
            prediction_file = os.path.join(options.output_path,root+'.results')
            
            if hash_kernel.get('kernel_type') == 'sigmoid' or hash_kernel.get('kernel_type') == 'polynomial':
                raise Exception('WARNING','run_svm: The model selection option has been tuned for linear or radial-basis kernels')
            
            grid_soft_margin,grid_gamma = generate_grid_search(hash_kernel['model_selection'],hash_kernel.get('kernel_type'))
            
            hash_grid = {}
            hash_grid['soft_margin']  = grid_soft_margin
            hash_grid['gamma'] = grid_gamma
            
            if options.K_fold:
                
                K = hash_kernel['k_fold']
                
                if K <= 1:
                    raise Exception('ERROR','run_svm.split_file_into_K_fold: The number of folds must be greater than or equal to 2')
                
                hash_filenames_K = split_file_into_K_fold(input_file,options.output_path,K)
                
                best_param,best_result = SVM.model_selection_grid(hash_filenames_K, model_file, alpha_file, hash_grid, k_fold = K, opt_criteria = hash_kernel['model_selection']['opt_criteria'])
                
                SVM.configure_params(best_param)
                
                SVM.train_K_fold(hash_filenames_K,model_file,alpha_file)
                        
            else:
                
                best_param,best_result = SVM.model_selection_grid(input_file, model_file, alpha_file, hash_grid, opt_criteria = hash_kernel['model_selection']['opt_criteria'])
                
                SVM.configure_params(best_param)
            
                SVM.train(input_file,model_file,alpha_file)    
                
        elif options.mode == 'e': # Test
            
            input_file = options.pattern_test_filename
            model_file = os.path.join(model_dir,options.model_filename)
            #(root,ext) = os.path.splitext(options.model_filename)
            (root,ext) = os.path.splitext(osp.basename(input_file))
            prediction_file = os.path.join(options.output_path,osp.basename(input_file))
            log_file = os.path.join(options.output_path,root+'.log')
            
            t1 = time.time()
            sys.stdout.write('\nTesting the model %s with the file %s... ' % (model_file,input_file))
            test_result = SVM.test(model_file, input_file, prediction_file, log_file)
            t2 = time.time()
            sys.stdout.write("done [%0.3fs].\n" % (t2-t1))
            sys.stdout.write('Results:\n')
            sys.stdout.write(test_result.__repr__()+'\n')
            sys.stdout.flush()
            
        else:
            raise Exception('ERROR','run_svm: The mode option is not properly configured')
        
    except:
        print >> sys.stderr , '\n%s\t%s' % (sys.exc_info()[0],sys.exc_info()[1])
        sys.exit(2)

##########################################################################

if __name__=='__main__':
    
    run()
