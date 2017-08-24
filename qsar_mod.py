# -*- coding: utf-8 -*-
'''
Last modified on March 9, 2017

@author: Yiting Ju
'''

import json
import logging
import os
import shutil
import sys
from subprocess import Popen
from parsing import parse_episuite
#from chem_spider_api import ChemSpiderAPI
#from batch_files import Call_EPI
#from batch_files import Call_VEGA
#from batch_files import Call_TEST

#from parsing import parse_vega
#from parsing import parsing
import inspect
import csv
import time
# import datetime

# get the directory where this "qsar.py" file is
class_directory = os.path.dirname(os.path.abspath(
	inspect.getfile(inspect.currentframe()))
)

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

"""
	class_directory:
	<qsar>
	|__[smiles.txt]		# the file which smiles list will be write to
	|__configuration.json	# configuration file specifying some static parameters
	|__[process]
	|  |__[*currentTimeAsID*]	# the time in seconds since the epoch as a floating
	|							#   point number
	|__[results]
	|  |__[*currentTimeAsID*]
	|
	|
	|
	|__test_seeds
	|__batch_files
	|__parsing
	|__sikuli_scripts
	|__chem_spider_api	
"""

"""
	self._directory
	self.id
	self.logging_file
	self.smiles_path
	self.batch_folder
	self.results_folder
	self.sikuli_scripts
	self.seed_folder
	self.process_folder
	self.config
	self.script_list
	self.epi_script
"""
class QSARmod:
	def __init__(self):
		# define paths
		self._directory = class_directory
		# locate batch folder, sikuli folder and seed folder
		self.batch_folder = os.path.join(self._directory, 'batch_files')
		self.sikuli_scripts = os.path.join(self._directory, 'sikuli_scripts')
		self.seed_folder = os.path.join(self._directory, 'test_seeds')

		self.id = ""
		self.process_folder = ""
		self.results_folder = ""
		self.smiles_path = ""
		self.logging_file = ""

		# open config file and create dict from contents
		config_file = open(os.path.join(self._directory, 'configuration.json'), 'r')
		self.config = json.loads(config_file.read())
		config_file.close()

		# construct batch files
		#self.script_list = ["epi", "test", "vega"]  
		
		self.script_list = ["epi"]

		# only epi needs sikulix process
		# self.config['sikuli_cmd'] defines the path to "/sikulix/runsikulix" which trigger sikulix
		
		self.epi_script = (
				"{0} -r {1}".format(self.config['sikuli_cmd'],
								os.path.join(self.sikuli_scripts, 'epi' + '_script.sikuli')
							))
		print("Epi Script: \n",self.epi_script)
		logging.info(self.epi_script)


	def run(self,input_hash = {}):
		logging.info("Running the QSAR module...")
		start_time = time.time()
		self.id = str(int(round(start_time)))	# get current time in second, used for id
		print('input_hash:',input_hash)
		pass

		# Define paths for the folders and files	
		self.logging_file = os.path.join(self._directory, 'logging', self.id, 'logging.txt')
		self.process_folder = os.path.join(self._directory, 'process', self.id)
		self.results_folder = os.path.join(self._directory, 'results', self.id)
		self.smiles_path = os.path.join(self.process_folder, 'smiles.txt')

		# Create folders if they do not exist
		if not os.path.exists(self.process_folder):
			os.makedirs(self.process_folder)
		if not os.path.exists(self.results_folder):
			os.makedirs(self.results_folder)

		# Clear out old results and process files. Disabling may cause
		#   parsing of old results in the case of script failure
		if self.config['clear_results'] == True:
			for the_file in os.listdir(self.results_folder):
				file_path = os.path.join(self.results_folder, the_file)
				try:
					if os.path.isfile(file_path):
						os.unlink(file_path)		# delete the file
					elif os.path.isdir(file_path):
						shutil.rmtree(file_path)	# delete the directory
				except Exception, e:
					print e
		if self.config['clear_process'] == True:
			for the_file in os.listdir(self.process_folder):
				file_path = os.path.join(self.process_folder, the_file)
				try:
					if os.path.isfile(file_path):
						os.unlink(file_path)		# delete the file
					elif os.path.isdir(file_path):
						shutil.rmtree(file_path)	# delete the directory
				except Exception, e:
					print e

		smiles_file = ""
		epi_smiles_file = ""
		# allow varying input types for testing purposes by using an options hash
		#   take {'smiles_in': smiles_value} or {'file_in': results_text_file}
		if not input_hash:	# no arguments => run with seed files
			smiles_file = os.path.join(self.seed_folder, "smiles.txt")
			epi_smiles_file = os.path.join(self.process_folder, "smiles.txt")

		elif "smiles_in" in input_hash:
			print("smiles_in")
			smiles_file = os.path.join(self.process_folder, "smiles.txt")
			# EPI Suite is running in batch mode using a txt file as input,
			#   so create a text file of smiles
			write_smiles_to_file(input_hash["smiles_in"], smiles_file)

			epi_smiles_file = os.path.join(self.process_folder, "epi_smiles.txt")
			# as a placeholder given epi one valid chemical to process, 
			#    so program abortion could be avoided
			write_smiles_to_file("CC," + input_hash["smiles_in"], epi_smiles_file)

		elif 'file_in' in input_hash:
			print('file_in')
			smiles_file = input_hash['file_in']
			"""
				To do: add CC as a place holder...
					smiles_file got here is actually a file which should be further processed
			"""

		# Batch call for processing
		# self.config --> /qsar/configuration.txt
		# print self.script_list
		for script in self.script_list:		# self.script_list = ["epi", "test", "vega"]
			if self.config['run_{0}'.format(script)]:
				if script == "epi" and self.config['enable_sikuli']:
					epi_start_time = time.time()
					logging.info("Call EPI suites...")
					epiOutputResultsFilePath = os.path.join(self.process_folder, 
															"epi_result",
															"epi_results.txt")
					Call_EPI.EPI_batch_allEndpoints(epi_smiles_file,
													self.process_folder,
													epiOutputResultsFilePath,
													self.epi_script)
					logging.info("--- EPI call finishes in %s seconds ---" % (time.time() - epi_start_time))
		# parse results
		if self.config['parse_results']:
			parse_episuite.read_epi_result_toJson(os.path.join(self.process_folder,
																"epi_result",
																"epi_results.txt"),
																os.path.join(self.process_folder,
																"epi_result",
																"epi_results.json"))
'''	
				elif script == "test":
					test_start_time = time.time()
					logging.info("Call TEST...")
					testResultsFolderPath = os.path.join(self.process_folder,
														"test_result")
					testCMDFilePath = self.config['test_cmd']
					# print smiles_file
					# print testResultsFolderPath
					Call_TEST.TEST_batch_allEndpoints(testCMDFilePath,
														smiles_file,
														testResultsFolderPath)
					# Call_TEST.TEST_batch_allEndpoints(testCMDFilePath, smiles_file, testResultsFolderPath, endpointList = [3,4])
					logging.info("--- TEST call finishes in %s seconds ---" % (time.time() - test_start_time))

				elif script == "vega":
					vega_start_time = time.time()
					logging.info("Call VEGA...")
					vegaJarFilePath = self.config['vega_cmd']
					vegaProcessFolder = os.path.join(self.process_folder,
													"vega_result")
					vegaScriptFilePath = os.path.join(vegaProcessFolder,
													"vega_script") 	# it's a file, not a folder
					vegaSmilesSourceFilePath = smiles_file
					outputResultFilePath = os.path.join(vegaProcessFolder,
														"vega_results.txt")
					Call_VEGA.VEGA_batch_allEndpoints(vegaJarFilePath,
														vegaScriptFilePath,
														vegaSmilesSourceFilePath,
														outputResultFilePath)
					logging.info("--- VEGA call finishes in %s seconds ---" % (time.time() - vega_start_time))

				else:
					print "[Error] something is wrong"
	
		# parse results
		if self.config['parse_results']:
			parse_episuite.read_epi_result_toJson(os.path.join(self.process_folder,
																"epi_result",
																"epi_results.txt"),
																os.path.join(self.process_folder,
																"epi_result",
																"epi_results.json"))
		
			parse_vega.read_vega_result_toJSON(os.path.join(vegaProcessFolder, 
															"vega_results.txt"),
												os.path.join(vegaProcessFolder,
															"vega_results.json"))

			qsar_parsed_summary = parsing.parse(parsing.readJSON(os.path.join(self.process_folder,
																		"epi_result",
																		"epi_results.json")),
											parsing.readJSON(os.path.join(vegaProcessFolder,
																		"vega_results.json")),
											parsing.readJSON(os.path.join(self.process_folder,
																		"test_result",
																		"test_results.json")),
											os.path.join(self.results_folder,
														"qsar_summary.json"))
			logging.info("--- finish in %s seconds ---" % (time.time() - start_time))
			return self.id, qsar_parsed_summary

		return None, None


	def yes(self):
		print "yes"
'''

def write_smiles_to_file(smiles_to_write, smiles_path):
	smiles_list_to_write = split_smiles_string(smiles_to_write)
	print(smiles_path,"\n")
	smiles_file = open(smiles_path, 'w+')
	for smiles in smiles_list_to_write:
		smiles_file.write(smiles + "\n")
	smiles_file.close()


def split_smiles_string(smiles_string):
	smiles_list = smiles_string.split(',')
	list_to_return = []
	for smiles in smiles_list:
		list_to_return.append(smiles.strip())
	return list_to_return


if __name__ == '__main__':

	q = QSARmod()
	#QSARmod().run()
	#q.run({'file_in': q.smiles_path})
	q.run(input_hash={"smiles_in":"c1ccccc1"})


