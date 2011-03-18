#!/usr/bin/python

import re, string,math

import sys,os


def isInt(s):

	try:
		i = int(s)
	except ValueError:
		return False
	return True

def isFloat(s):
	try:
		i = float(s)
	except ValueError:
		return False
	return True

def isNumber(s):
	return isInt(s) or isFloat(s)

def frange(start, end, inc=None):

		"""
			start and end are includes in the returned list
		"""
	
		if inc == None:
			inc = 1.0
	
		L = [start]
	
		while 1:
	
			next = start + len(L) * inc
			
			if inc > 0 and next == end:
				L.append(next)
				break
			elif inc > 0 and next > end:
				break
			elif inc < 0 and next <= end:
				break
			
			L.append(next)
			
		return L

def sigmoid(inpt,alpha=1):
	
	return 1/(1+math.exp(-alpha*inpt))

def div_sec(num,deno):
	
	try:
		return float(num)/deno
	except ZeroDivisionError:
		return float(0)

def log_sec(value):

	try:
		log_value = math.log(value)
	except:
		log_value = float(0) 

	return log_value

############################################################################    