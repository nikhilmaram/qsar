# -*- coding: utf-8 -*-
'''
Last modified on March 8, 2017

@author: Yiting Ju
'''
import json
from pprint import pprint

		
def readJSON(jsonFilePath):
	with open(jsonFilePath) as jsonFile:
		jsonData = json.load(jsonFile)
		# pprint(jsonData)
	return jsonData




if __name__ == '__main__':
	print "TEST--Parsing"
	# print readJSON("/home/yiting/Dropbox/Spring2016/CLiCC/ModuleIntegration/clicc-flask-master/modules/qsar/parsing/TEST_summary_sample.json")



