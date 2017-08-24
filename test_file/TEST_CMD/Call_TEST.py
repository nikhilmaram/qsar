from subprocess import Popen
import shutil, os, csv, json
from multiprocessing import Pool
from TESTHeaderRow import TEST_Headers 


default_inputFilePath = "./test_chemicals/Sample_Smiles_List.smi"
default_outputFolderPath = "./test_chemicals/cmd_results/"
default_fileType = 2	# smiles
default_method = 10
# DEFAULT_ENDPOINT_LIST = [1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28]
DEFAULT_ENDPOINT_LIST = [1, 2]

# folers always end with a '/'


'''
	Call the TEST to get QSAR of all endpoints for a batch of test_chemicals
	  The results are categorized for each chemicals in a single folder
	  Everything will be summarized in a single summary.json
	inputFilePath: Path to a .smi file of a list of chemicals' smiles (one smiles for each line)
	outputFolderPath: Path to a folder where the result will be stored
	fileType: 1 for SDF; 2 for SMILES (.smi); [We use 2]
	method: 10 for Consensus
	endpointList: The list of endpoints to be calculated

	Folder structure:
		cmd_results
		-smile1
		|-result_endpoint_1.csv
		|-result_endpoint_2.csv
		|-...
		|-summary.csv
		-smile2
		|-...
		...
		...
		-summary.json

	summary.json:
	[
		{
			"No.": "1",
 			"SMILES": "OCCOCCOCCOCCOCCOCC(C)OCC(C)OCC(C)OCC(C)OCC(C)O",
 			"Bioaccumulation factor  Exp_Value:dyn/cm": "N/A",
        	"Bioaccumulation factor  Pred_Value:dyn/cm": "7204.49",
 			...
 			...
		},
		{
			"No.": "2",
 			"SMILES": "KO",
 			"Bioaccumulation factor  Exp_Value:dyn/cm": "N/A",
        	"Bioaccumulation factor  Pred_Value:dyn/cm": "N/A",
 			...
 			...
		},
		...
		...
		{
			...
			...
		}
	]
'''
def TEST_batch_allEndpoints(inputFilePath, outputFolderPath, 
					fileType = 2, method = 10, endpointList = DEFAULT_ENDPOINT_LIST):
	if not os.path.exists(outputFolderPath):
			os.makedirs(outputFolderPath)
	folderPath = inputFilePath[:inputFilePath.rfind('/')+1]
	# print folderPath   	# ./test_chemicals/
	numSmiles, smilesList = listToFiles(inputFilePath, folderPath + "tempSMIFolder/")
	testResultList = []
	for i in range(numSmiles):
		tempInputPath = folderPath + "tempSMIFolder/" + "smile_" + str(i+1) + ".smi"
		outputFolderPathIns = outputFolderPath + "smile" + str(i+1) + "/"
		if not os.path.exists(outputFolderPathIns):
			os.makedirs(outputFolderPathIns)
		TEST_allEndpoints(tempInputPath, outputFolderPathIns, fileType, method, endpointList)
		resultList = organizeResultToSingleCSV(outputFolderPathIns, outputFolderPathIns+"summary.csv", endpointList)
		testResultList.append(resultList)
	writeJSONSummary(testResultList, smilesList, outputFolderPath)
	shutil.rmtree(folderPath + "tempSMIFolder/")
	print "temp folder " + folderPath + "tempSMIFolder/" + " is removed"


"""
	Call the TEST cmd to get QSAR results for (all) endpoints in endpointList
		multi-processing is employed in this function
"""
def TEST_allEndpoints(inputFilePath, outputFolderPath, 
						fileType = 2, method = 10, endpointList = DEFAULT_ENDPOINT_LIST):
	commands = []
	for endpoint in endpointList:
		# command = "java -jar ../VEGA-CLI.jar -script ./script_allModule_test"
		command = 'java -Xmx512m -cp "TEST.jar" ToxPredictor.Application.runTEST_From_Command_Line "' + \
			inputFilePath + '" ' + str(fileType) + \
			' "' + (outputFolderPath+"result_endpoint_{0}.csv").format(endpoint) + \
			'" ' + str(endpoint) + ' ' + str(method)
		# commandList = command.split(' ')
		commands.append(command)
	pool = Pool(processes = 5)
	pool.map(callPopen, commands)


def callPopen(command):
	print "Calling", command
	try:
	  	e = Popen(
	      	command,
	      	# cwd="/home/yiting/Downloads/Vega-1.1.1-binaries/vega-cli beta 1_cmd",
	      	cwd = "/home/yiting/Dropbox/Spring2016/CLiCC/ModuleIntegration/clicc-flask-master/modules/qsar/batch_files/TEST_CMD",
			shell=True
		)
		stdout, stderr = e.communicate()

	except IOError as (errno,strerror):
		print "I/O error({0}): {1}".format(errno, strerror)


def listToFiles(inputListFilePath, outputFolder):
	smilesList = []
	with open(inputListFilePath, 'r') as f:
		for row in f:
			smilesList.append(row)
	if not os.path.exists(outputFolder):
		os.makedirs(outputFolder)
		print "temp folder " + outputFolder + " is created"
	for i in range(len(smilesList)):
		with open(outputFolder + "smile_" + str(i+1) + ".smi", 'w') as f:
			f.write(smilesList[i])
	return len(smilesList), smilesList



def organizeResultToSingleCSV(smilesResultFolderPath, summaryFilePath, endpointList = DEFAULT_ENDPOINT_LIST):
	headerRowList = []
	dataRowList = []
	for i in endpointList:
		headerRow, dataRow = readTESTResult(smilesResultFolderPath + "result_endpoint_"+str(i)+".csv")
		headerRowList.append(headerRow)
		dataRowList.append(dataRow)
	with open(summaryFilePath, 'wb') as csvfile:
		spamwriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		for dataRowIns in dataRowList:
			spamwriter.writerow(dataRowIns)
	return dataRowList


def readTESTResult(resultFilePath):
	headerRow = []
	dataRow = []
	with open(resultFilePath, 'rb') as csvfile:
		csvReader = csv.reader(csvfile, delimiter='\t')
		row_count = 0
		for row in csvReader:
			if row_count == 0:
				headerRow = row[2:-1]	# the first and last row are discarded
				row_count += 1
			elif row_count == 1:
				dataRow = row[2:-1]
				row_count += 1
		if row_count != 2:
			print "error: the result file should only have two lines"
	return headerRow, dataRow


def writeJSONSummary(testResultList, smilesList, outputFolderPath):
	resultJsonObject = []
	resultCounter = 0
	for testResult in testResultList:		
		itemJsonObject = {}
		itemJsonObject["No."] = resultCounter+1
		itemJsonObject["Smiles"] = (smilesList[resultCounter]).strip()
		itemCount = 0
		for itemList in testResult:
			for item in itemList:
				itemJsonObject[TEST_Headers[itemCount]] = item
				itemCount += 1
		resultJsonObject.append(itemJsonObject)
		resultCounter += 1

	# write JSON to file
	jsonOutputPath = outputFolderPath + "summary.json"
	with open(jsonOutputPath, "w") as outputFile:
		json.dump(resultJsonObject, outputFile, sort_keys=True, indent= 4, separators=(',', ': '))



if __name__ == "__main__":
	print "TEST"
	TEST_batch_allEndpoints("/home/yiting/Dropbox/Spring2016/CLiCC/ModuleIntegration/clicc-flask-master/modules/qsar/test_seeds/smiles.txt",
							"/home/yiting/Dropbox/Spring2016/CLiCC/ModuleIntegration/clicc-flask-master/modules/qsar/process/test_result/")
	# TEST_batch_allEndpoints("./test_chemicals/source_test_onlyOne.txt", default_outputFolderPath)
	# for i in default_endpointList:
	# 	print readTESTResult(default_outputFolderPath + "simle_1/result_endpoint_"+str(i)+".txt")

	# organizeResultToSingleCSV(default_outputFolderPath + "simle_2/", default_outputFolderPath + "simle_1/summary.csv")