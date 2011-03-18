#!/usr/bin/python

import sys, re, os, string, urllib, time, math, random, tempfile

from os import path as osp

from scipy.stats import chisqprob

from math import log

from mod_utils import log_sec

############################################################################

#### The contingency table is:
####	   	   	      	   	   
####	   	   	         X
####	   	         ----------
####	   	         w    not_w  total
####	  ---------------------------
####	Y | w        n11  n12     n1+
####	  | not_w    n21  n22     n2+
####      ---------------------------
####        total    n+1  n+2     n


def get_p_value_pearson_chi_squared(contingency_table):
	
	(n11,n12,n21,n22) = contingency_table

	n = n11+n12+n21+n22
	
	if n == 0:
		
		raise "The contingency table is empty"
	
	n_1_plus = float(n11 + n12)
	n_2_plus = float(n21 + n22)
	n_plus_1 = float(n11 + n21)
	n_plus_2 = float(n12 + n22)
	
	if n == n_1_plus:

		return float(1)
	
	elif n == n_2_plus:
		
		return float(1)
		
	# eij = (n_i_plus)(n_plus_j)/n
	e11 = (n_1_plus)*(n_plus_1)/n
	e12 = (n_1_plus)*(n_plus_2)/n
	e21 = (n_2_plus)*(n_plus_1)/n
	e22 = (n_2_plus)*(n_plus_2)/n
					
	chi2 = (math.pow(n11-e11,2)/e11) + (math.pow(n12-e12,2)/e12) + (math.pow(n21-e21,2)/e21) + (math.pow(n22-e22,2)/e22)

	p_value = chisqprob(chi2,1)

	return p_value

def get_log_likelihood(contingency_table):

	(n11,n12,n21,n22) = contingency_table

	g2 = n11*log_sec(n11)+n21*log_sec(n21)+n12*log_sec(n12)+n22*log_sec(n22)-(n11+n21)*log_sec(n11+n21)-(n11+n12)*log_sec(n11+n12)-(n21+n22)*log_sec(n21+n22)-(n12+n22)*log_sec(n12+n22)+(n11+n21+n12+n22)*log_sec(n11+n21+n12+n22)

	return 2*abs(g2)

def p_value(value,n=2): 

	p_value = chisqprob(value, n-1)

	return p_value


def get_p_value_likelihood_ratio_chi_squared(contingency_table):
	
	(n11,n12,n21,n22) = contingency_table
	
	n = n11+n12+n21+n22
	
	if n == 0:
		
		raise "The contingency table is empty"
		
	n_1_plus = float(n11 + n12)
	n_2_plus = float(n21 + n22)
	n_plus_1 = float(n11 + n21)
	n_plus_2 = float(n12 + n22)
	
	if n == n_1_plus:

		return float(1)
	
	elif n == n_2_plus:
		
		return float(1)
		
	# eij = (n_i_plus)(n_plus_j)/n
	e11 = (n_1_plus)*(n_plus_1)/n
	e12 = (n_1_plus)*(n_plus_2)/n
	e21 = (n_2_plus)*(n_plus_1)/n
	e22 = (n_2_plus)*(n_plus_2)/n
					
	try:
		
		g2 =  2*(n11*math.log(n11/e11) + n22*math.log(n12/e12) + n21*math.log(n21/e21) + n22*math.log(n22/e22))
	
	except ZeroDivisionError:
		return float(1)
	
	p_value = chisqprob(g2,1)

	return p_value

def get_odds_ratio(contingency_table):
	
	(n11,n12,n21,n22) = contingency_table
	
	try:
		odd = float(n11*n22)/n12*n21
	except ZeroDivisionError:
		odd = 0
		
	return odd
