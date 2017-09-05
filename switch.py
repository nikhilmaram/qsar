#!/usr/bin/python
from __future__ import print_function
import os,sys
sys.path.insert(1, os.path.join(sys.path[0], 'qsar'))
#print(sys.path)
import json, math
import numpy as np
import pprint
import time
#from qsar_mod import QSARmod
#import qsar_mod
from multiprocessing import cpu_count
from parse import parse
from vega_file import Call_VEGA
import subprocess
import re
import uuid 
import hashlib
from episuite_file.parse_episuite import read_epi_result_toJson
from vega_file.parse_vega import read_vega_result_toJSON
from test_file.Call_TEST import TEST_batch_allEndpoints

DIR_PATH = os.path.dirname(os.path.realpath(__file__))

def readJSON(jsonFilePath):
    with open(jsonFilePath) as jsonFile:
        jsonData = json.load(jsonFile)
        # pprint(jsonData)
    return jsonData

def make_dir_if_necessary(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

epinone=    {
        "BP  C  est": "N/A",
        "BP  C  exp": "N/A",
        "CAS": "N/A",
        "HLC  Pa-m3/mole  Bond":"N/A" ,
        "HLC  Pa-m3/mole  Group":"N/A" ,
        "HLC  Pa-m3/mole  exp":"N/A",
        "HLDegAero  hour": "N/A",
        "HLDegAir  hour":"N/A" ,
        "HLDegSed  hour": "N/A",
        "HLDegSoil  hour": "N/A",
        "HLDegSsed  hour": "N/A",
        "HLDegWater  hour": "N/A",
        "Kb_HL_pH7  days":"N/A" ,
        "Kb_HL_pH8  days": "N/A",
        "Kb_rateC  L/mol-sec": "N/A",
        "Km_10  /day": "N/A",
        "LogBCF  L/kg wet-wt  Arnot-Gobas":"N/A" ,
        "LogBCF  L/kg wet-wt  Regression": "N/A",
        "MP  C  est":"N/A" ,
        "MP  C  exp":"N/A" ,
        "MW  g/mol": "N/A",
        "Mol Formula": "N/A",
        "OH_HL  days": "N/A",
        "SMILES": "N/A",
        "VP  mmHg  est": "N/A",
        "VP  mmHg  exp":"N/A" ,
        "WS  mg/L  WATERNT  est":"N/A" ,
        "WS  mg/L  WATERNT  exp":"N/A" ,
        "WS  mg/L  WSKOW  est":"N/A" ,
        "WS  mg/L  WSKOW  exp":"N/A" ,
        "WWTair  %  10000hr":"N/A",
        "WWTair  %  Biowin/EPA":"N/A" ,
        "WWTbio  %  10000hr": "N/A",
        "WWTbio  %  Biowin/EPA": "N/A",
        "WWTremoval  %  10000hr": "N/A",
        "WWTremoval  %  Biowin/EPA": "N/A",
        "WWTslu  %  10000hr":"N/A" ,
        "WWTslu  %  Biowin/EPA": "N/A",
        "algaeChV_ecosar  mg/L": "N/A",
        "algaeEC50_96hr_ecosar  mg/L": "N/A",
        "aquaTox_acute  unitless": "N/A",
        "bioHC_HL  days": "N/A",
        "bio_HL  days": "N/A",
        "biodeg_MITIlinear  unitless":"N/A" ,
        "biodeg_MITInonlinear  unitless": "N/A",
        "biodeg_anaerobic  unitless": "N/A",
        "biodeg_linear  unitless":"N/A" ,
        "biodeg_nonlinear  unitless": "N/A",
        "biodeg_primary  unitless": "N/A",
        "biodeg_ready  unitless":"N/A" ,
        "biodeg_ultimate  unitless":"N/A" ,
        "biotrans_HL  days":"N/A" ,
        "dmChV_ecosar  mg/L":"N/A" ,
        "dmLC50_48hr_ecosar  mg/L":"N/A" ,
        "earthworm_14day_ecosar  mg/L":"N/A" ,
        "fishChVSW_ecosar  mg/L":"N/A" ,
        "fishChV_ecosar  mg/L":"N/A" ,
        "fishLC50SW_96hr_ecosar  mg/L":"N/A",
        "fishLC50_96hr_ecosar  mg/L":"N/A" ,
        "kAerAir  m3/ug  Koa":"N/A" ,
        "kAerAir  m3/ug  Mackay":"N/A" ,
        "kAirWater  unitless":"N/A",
        "kNO3  cm3/molecule-sec": "N/A",
        "kO3  cm3/molecule-sec": "N/A",
        "kOH  cm3/molecule-sec": "N/A",
        "kOctAir  unitless  est": "N/A",
        "kOctAir  unitless  exp": "N/A",
        "kOctWater  unitless  est":"N/A" ,
        "kOctWater  unitless  exp": "N/A",
        "kOrgWater  L/kg  Kow": "N/A",
        "kOrgWater  L/kg  MCI": "N/A",
        "kOrgWater  L/kg  exp": "N/A",
        "shrimpLC50_96hr_ecosar  mg/L": "N/A",
        "shrimpSWChV_ecosar  mg/L": "N/A"
    }




# empty TEST json component
TEST_Headers = [
    'Fathead minnow LC50 (96 hr)  Exp_Value:-Log10(mol/L)',
    'Fathead minnow LC50 (96 hr)  Pred_Value:-Log10(mol/L)',
    'Fathead minnow LC50 (96 hr)  Exp_Value:mg/L',
    'Fathead minnow LC50 (96 hr)  Pred_Value:mg/L',
    'Daphnia magna LC50 (48 hr)  Exp_Value:-Log10(mol/L)',
    'Daphnia magna LC50 (48 hr)  Pred_Value:-Log10(mol/L)',
    'Daphnia magna LC50 (48 hr)  Exp_Value:mg/L',
    'Daphnia magna LC50 (48 hr)  Pred_Value:mg/L',
    'T. pyriformis IGC50 (48 hr)  Exp_Value:-Log10(mol/L)',
    'T. pyriformis IGC50 (48 hr)  Pred_Value:-Log10(mol/L)',
    'T. pyriformis IGC50 (48 hr)  Exp_Value:mg/L',
    'T. pyriformis IGC50 (48 hr)  Pred_Value:mg/L',
    'Oral rat LD50  Exp_Value:-Log10(mol/kg)  HC',
    'Oral rat LD50  Pred_Value:-Log10(mol/kg)  HC',
    'Oral rat LD50  Exp_Value:mg/kg  HC',
    'Oral rat LD50  Pred_Value:mg/kg  HC',
    'Oral rat LD50  Exp_Value:-Log10(mol/kg)  FDA',
    'Oral rat LD50  Pred_Value:-Log10(mol/kg)  FDA',
    'Oral rat LD50  Exp_Value:mg/kg  FDA',
    'Oral rat LD50  Pred_Value:mg/kg  FDA',
    'Oral rat LD50  Exp_Value:-Log10(mol/kg)  NN',
    'Oral rat LD50  Pred_Value:-Log10(mol/kg)  NN',
    'Oral rat LD50  Exp_Value:mg/kg  NN',
    'Oral rat LD50  Pred_Value:mg/kg  NN',
    'Oral rat LD50  Exp_Value:-Log10(mol/kg)  C',
    'Oral rat LD50  Pred_Value:-Log10(mol/kg)  C',
    'Oral rat LD50  Exp_Value:mg/kg  C',
    'Oral rat LD50  Pred_Value:mg/kg  C',
    'Bioaccumulation factor  Exp_Value:Log10',
    'Bioaccumulation factor  Pred_Value:Log10',
    'Bioaccumulation factor  Exp_Value:',
    'Bioaccumulation factor  Pred_Value:',
    'Developmental Toxicity  Exp_Value',
    'Developmental Toxicity  Pred_Value',
    'Developmental Toxicity  Exp_Result',
    'Developmental Toxicity  Pred_Result',
    'Mutagenicity  Exp_Value',
    'Mutagenicity  Pred_Value',
    'Mutagenicity  Exp_Result',
    'Mutagenicity  Pred_Result',
    'Normal boiling point  Exp_Value:\xc2\xb0C',
    'Normal boiling point  Pred_Value:\xc2\xb0C',
    'Vapor pressure at 25\xc2\xb0C  Exp_Value:Log10(mmHg)',
    'Vapor pressure at 25\xc2\xb0C  Pred_Value:Log10(mmHg)',
    'Vapor pressure at 25\xc2\xb0C  Exp_Value:mmHg',
    'Vapor pressure at 25\xc2\xb0C  Pred_Value:mmHg',
    'Melting point  Exp_Value:\xc2\xb0C',
    'Melting point  Pred_Value:\xc2\xb0C',
    'Flash point  Exp_Value:\xc2\xb0C',
    'Flash point  Pred_Value:\xc2\xb0C',
    'Density  Exp_Value:g/cm\xc2\xb3  HC',
    'Density  Pred_Value:g/cm\xc2\xb3  HC',
    'Density  Exp_Value:g/cm\xc2\xb3  FDA',
    'Density  Pred_Value:g/cm\xc2\xb3  FDA',
    'Density  Exp_Value:g/cm\xc2\xb3  NN',
    'Density  Pred_Value:g/cm\xc2\xb3  NN',
    'Density  Exp_Value:g/cm\xc2\xb3  GC',
    'Density  Pred_Value:g/cm\xc2\xb3  GC',
    'Density  Exp_Value:g/cm\xc2\xb3  C',
    'Density  Pred_Value:g/cm\xc2\xb3  C',
    'Surface tension at 25\xc2\xb0C  Exp_Value:dyn/cm',
    'Surface tension at 25\xc2\xb0C  Pred_Value:dyn/cm',
    'Thermal conductivity at 25\xc2\xb0C  Exp_Value:mW/mK',
    'Thermal conductivity at 25\xc2\xb0C  Pred_Value:mW/mK',
    'Viscosity at 25\xc2\xb0C  Exp_Value:Log10(cP)',
    'Viscosity at 25\xc2\xb0C  Pred_Value:Log10(cP)',
    'Viscosity at 25\xc2\xb0C  Exp_Value:cP',
    'Viscosity at 25\xc2\xb0C  Pred_Value:cP',
    'Water solubility at 25\xc2\xb0C  Exp_Value:-Log10(mol/L)',
    'Water solubility at 25\xc2\xb0C  Pred_Value:-Log10(mol/L)',
    'Water solubility at 25\xc2\xb0C  Exp_Value:mg/L',
    'Water solubility at 25\xc2\xb0C  Pred_Value:mg/L'
]


#empty vega component
veganone={
        "BCF model (CAESAR) - ADI": "N/A",
        "BCF model (CAESAR) - assessment": "N/A",
        "BCF model (CAESAR) - experimental value": "N/A",
        "BCF model (CAESAR) - prediction [log(L/kg)]": "N/A",
        "BCF model (KNN/Read-Across) - ADI": "N/A",
        "BCF model (KNN/Read-Across) - assessment": "N/A",
        "BCF model (KNN/Read-Across) - experimental value": "N/A",
        "BCF model (KNN/Read-Across) - prediction [log(L/kg)]": "N/A",
        "BCF model (Meylan) - ADI": "N/A",
        "BCF model (Meylan) - assessment": "N/A",
        "BCF model (Meylan) - experimental value": "N/A",
        "BCF model (Meylan) - prediction [log(L/kg)]": "N/A",
        "Carcinogenicity model (CAESAR) - ADI": "N/A",
        "Carcinogenicity model (CAESAR) - assessment": "N/A",
        "Carcinogenicity model (CAESAR) - experimental value": "N/A",
        "Carcinogenicity model (CAESAR) - prediction": "N/A",
        "Carcinogenicity model (IRFMN/Antares) - ADI": "N/A",
        "Carcinogenicity model (IRFMN/Antares) - assessment": "N/A",
        "Carcinogenicity model (IRFMN/Antares) - experimental value": "N/A",
        "Carcinogenicity model (IRFMN/Antares) - prediction": "N/A",
        "Carcinogenicity model (IRFMN/ISSCAN-CGX) - ADI": "N/A",
        "Carcinogenicity model (IRFMN/ISSCAN-CGX) - assessment": "N/A",
        "Carcinogenicity model (IRFMN/ISSCAN-CGX) - experimental value": "N/A",
        "Carcinogenicity model (IRFMN/ISSCAN-CGX) - prediction": "N/A",
        "Carcinogenicity model (ISS) - ADI": "N/A",
        "Carcinogenicity model (ISS) - assessment": "N/A",
        "Carcinogenicity model (ISS) - experimental value": "N/A",
        "Carcinogenicity model (ISS) - prediction": "N/A",
        "Daphnia Magna LC50 48h (DEMETRA) - ADI": "N/A",
        "Daphnia Magna LC50 48h (DEMETRA) - assessment": "N/A",
        "Daphnia Magna LC50 48h (DEMETRA) - experimental value": "N/A",
        "Daphnia Magna LC50 48h (DEMETRA) - prediction [-log(mol/l)]": "N/A",
        "Daphnia Magna LC50 48h (EPA) - ADI": "N/A",
        "Daphnia Magna LC50 48h (EPA) - assessment": "N/A",
        "Daphnia Magna LC50 48h (EPA) - experimental value": "N/A",
        "Daphnia Magna LC50 48h (EPA) - prediction [-log(mol/l)]": "N/A",
        "Developmental Toxicity model (CAESAR) - ADI": "N/A",
        "Developmental Toxicity model (CAESAR) - assessment": "N/A",
        "Developmental Toxicity model (CAESAR) - experimental value": "N/A",
        "Developmental Toxicity model (CAESAR) - prediction": "N/A",
        "Developmental/Reproductive Toxicity library (PG) - ADI": "N/A",
        "Developmental/Reproductive Toxicity library (PG) - assessment": "N/A",
        "Developmental/Reproductive Toxicity library (PG) - experimental value": "N/A",
        "Developmental/Reproductive Toxicity library (PG) - prediction": "N/A",
        "Estrogen Receptor Relative Binding Affinity model (IRFMN) - ADI": "N/A",
        "Estrogen Receptor Relative Binding Affinity model (IRFMN) - assessment": "N/A",
        "Estrogen Receptor Relative Binding Affinity model (IRFMN) - experimental value": "N/A",
        "Estrogen Receptor Relative Binding Affinity model (IRFMN) - prediction": "N/A",
        "Estrogen Receptor-mediated effect (IRFMN/CERAPP) - ADI": "N/A",
        "Estrogen Receptor-mediated effect (IRFMN/CERAPP) - assessment": "N/A",
        "Estrogen Receptor-mediated effect (IRFMN/CERAPP) - experimental value": "N/A",
        "Estrogen Receptor-mediated effect (IRFMN/CERAPP) - prediction": "N/A",
        "Fathead Minnow LC50 96h (EPA) - ADI": "N/A",
        "Fathead Minnow LC50 96h (EPA) - assessment": "N/A",
        "Fathead Minnow LC50 96h (EPA) - experimental value": "N/A",
        "Fathead Minnow LC50 96h (EPA) - prediction [-log(mol/l)]": "N/A",
        "Fish Acute (LC50) Toxicity classification (SarPy/IRFMN) - ADI": "N/A",
        "Fish Acute (LC50) Toxicity classification (SarPy/IRFMN) - assessment": "N/A",
        "Fish Acute (LC50) Toxicity classification (SarPy/IRFMN) - experimental value": "N/A",
        "Fish Acute (LC50) Toxicity classification (SarPy/IRFMN) - prediction": "N/A",
        "Fish Acute (LC50) Toxicity model (KNN/Read-Across) - ADI": "N/A",
        "Fish Acute (LC50) Toxicity model (KNN/Read-Across) - assessment": "N/A",
        "Fish Acute (LC50) Toxicity model (KNN/Read-Across) - experimental value": "N/A",
        "Fish Acute (LC50) Toxicity model (KNN/Read-Across) - prediction [-log(mg/L)]": "N/A",
        "Fish Acute (LC50) Toxicity model (NIC) - ADI": "N/A",
        "Fish Acute (LC50) Toxicity model (NIC) - assessment": "N/A",
        "Fish Acute (LC50) Toxicity model (NIC) - experimental value": "N/A",
        "Fish Acute (LC50) Toxicity model (NIC) - prediction [log(1/(mmol/L))]": "N/A",
        "Id": "Molecule 1",
        "LogP model (ALogP) - ADI": "N/A",
        "LogP model (ALogP) - assessment": "N/A",
        "LogP model (ALogP) - experimental value": "N/A",
        "LogP model (ALogP) - prediction": "N/A",
        "LogP model (MLogP) - ADI": "N/A",
        "LogP model (MLogP) - assessment": "N/A",
        "LogP model (MLogP) - experimental value": "N/A",
        "LogP model (MLogP) - prediction": "N/A",
        "LogP model (Meylan/Kowwin) - ADI": "N/A",
        "LogP model (Meylan/Kowwin) - assessment": "N/A",
        "LogP model (Meylan/Kowwin) - experimental value": "N/A",
        "LogP model (Meylan/Kowwin) - prediction": "N/A",
        "Mutagenicity (Ames test) model (CAESAR) - ADI": "N/A",
        "Mutagenicity (Ames test) model (CAESAR) - assessment": "N/A",
        "Mutagenicity (Ames test) model (CAESAR) - experimental value": "N/A",
        "Mutagenicity (Ames test) model (CAESAR) - prediction": "N/A",
        "Mutagenicity (Ames test) model (ISS) - ADI": "N/A",
        "Mutagenicity (Ames test) model (ISS) - assessment": "N/A",
        "Mutagenicity (Ames test) model (ISS) - experimental value": "N/A",
        "Mutagenicity (Ames test) model (ISS) - prediction": "N/A",
        "Mutagenicity (Ames test) model (KNN/Read-Across) - ADI": "N/A",
        "Mutagenicity (Ames test) model (KNN/Read-Across) - assessment": "N/A",
        "Mutagenicity (Ames test) model (KNN/Reasmilesd-Across) - experimental value": "N/A",
        "Mutagenicity (Ames test) model (KNN/Read-Across) - prediction": "N/A",
        "Mutagenicity (Ames test) model (SarPy/IRFMN) - ADI": "N/A",
        "Mutagenicity (Ames test) model (SarPy/IRFMN) - assessment": "N/A",
        "Mutagenicity (Ames test) model (SarPy/IRFMN) - experimental value": "N/A",
        "Mutagenicity (Ames test) model (SarPy/IRFMN) - prediction": "N/A",
        "No.": "1",
        "Persistence (sediment) model (IRFMN) - ADI": "N/A",
        "Persistence (sediment) model (IRFMN) - assessment": "N/A",
        "Persistence (sediment) model (IRFMN) - experimental value": "N/A",
        "Persistence (sediment) model (IRFMN) - prediction": "N/A",
        "Persistence (soil) model (IRFMN) - ADI": "N/A",
        "Persistence (soil) model (IRFMN) - assessment": "N/A",
        "Persistence (soil) model (IRFMN) - experimental value": "N/A",
        "Persistence (soil) model (IRFMN) - prediction": "N/A",
        "Persistence (water) model (IRFMN) - ADI": "N/A",
        "Persistence (water) model (IRFMN) - assessment": "N/A",
        "Persistence (water) model (IRFMN) - experimental value": "N/A",
        "Persistence (water) model (IRFMN) - prediction": "N/A",
        "Ready Biodegradability model (IRFMN) - ADI": "N/A",
        "Ready Biodegradability model (IRFMN) - assessment": "N/A",
        "Ready Biodegradability model (IRFMN) - experimental value": "N/A",
        "Ready Biodegradability model (IRFMN) - prediction": "N/A",
        "Remarks": "N/A",
        "SMILES": "N/A",
        "Skin Sensitization model (CAESAR) - ADI": "N/A",
        "Skin Sensitization model (CAESAR) - assessment": "N/A",
        "Skin Sensitization model (CAESAR) - experimental value": "N/A",
        "Skin Sensitization model (CAESAR) - prediction": "N/A"}

def change_vega_script(vega_script_template,vega_file_folder,vega_script_path):
    script = None
    vega_singel_result_path = os.path.abspath(os.path.normpath('{0}/vega_result.txt').format(vega_file_folder))
    vega_multiple_result_path = os.path.abspath(os.path.normpath('{0}/results').format(vega_file_folder))    
    with open(vega_script_template,'r') as f:
        script = f.read()

        script = re.sub(r'<source>.*</source>','<source>'+
                        os.path.abspath(os.path.normpath('{0}/vega_source.txt').format(vega_file_folder))
                        +'</source>',script)

        script = re.sub(r'<singleTXT>.*</singleTXT>','<singleTXT>' +
                        vega_singel_result_path + '</singleTXT>',script)
        
        script = re.sub(r'<multipleTXT>.*</multipleTXT>','<multipleTXT>' +
                        vega_multiple_result_path + '</multipleTXT>',script)
    # print(script)
    make_dir_if_necessary(os.path.abspath(os.path.normpath('{0}/results').format(vega_file_folder)))
    with open(vega_script_path,'w') as f:
        f.write("".join(script))
    
    return vega_singel_result_path
    # print("script changed")
    # <singleTXT>/home/awsgui/Desktop/qsar/vega_file/result_test.txt</singleTXT>
    # <multipleTXT>/home/awsgui/Desktop/qsar/vega_file/results</multipleTXT>

def serialize_smiles_and_generate_scripts(smiles,temp_dir_path,epi,vega,test):
    VEGA_SCRIPT_PATH = ""
    VEGA_RESULT_PATH = ""
    TEST_SMILES_PATH = ""
    TEST_RESULT_PATH = ""
    '''
    Our episuite script always runs in batch mode. In batch mode, we need
    feed at least 1 'correct' smiles to episuite so that it can generate output.
    Thus, to handle potential 'wrong' smiles, we put a place holder in epi_smiles.txt
    '''
    if test:
        # serialize smiles for EPI-suite
        epi_file_path = os.path.join(temp_dir_path, "episuite_file")
        make_dir_if_necessary(epi_file_path)
        epi_smiles = open(os.path.join(epi_file_path,"epi_smiles.txt"), "w")
        epi_smiles.write("CC\n"+smiles+"\n")
        epi_smiles.close()
        # modify sikulix script

    '''
        We could let epi, vega and test read the same input file. However,
        we think seperating input files and making them more indenpent
        is a better choice.
        Also, epi input file has an additional placeholder, 
        but vega/test does NOT need a placeholder)
    '''
    if vega:
        # serialize smiles for VEGA
        vega_file_folder = os.path.join(temp_dir_path, "vega_file")
        make_dir_if_necessary(vega_file_folder)
        vega_smiles = open(os.path.join(vega_file_folder,'vega_source.txt'), "w")
        vega_smiles.write(smiles+"\n")
        vega_smiles.close()
        # generate VEGA script from template and save the script in temp folder
        VEGA_SCRIPT_PATH = os.path.join(temp_dir_path,'vega_file/vega_script')
        VEGA_RESULT_PATH = change_vega_script(os.path.normpath(DIR_PATH + "/vega_file/script_shaoyi"),
                           vega_file_folder, VEGA_SCRIPT_PATH)

    if test:
        # serialize smiles for TEST
        test_file_path = os.path.join(temp_dir_path, "test_file")
        make_dir_if_necessary(test_file_path)
        TEST_RESULT_PATH = test_file_path
        TEST_SMILES_PATH = os.path.join(test_file_path,"test_smiles.txt")
        test_smiles = open(TEST_SMILES_PATH, "w")
        test_smiles.write(smiles+"\n")
        test_smiles.close()

    return {"VEGA_SCRIPT_PATH":VEGA_SCRIPT_PATH,
            "VEGA_RESULT_PATH":VEGA_RESULT_PATH,
            "TEST_SMILES_PATH":TEST_SMILES_PATH,
            "TEST_RESULT_PATH":TEST_RESULT_PATH}

# smile: string, epi,vega,test: boolean, 
def switch(smiles,epi,vega,test,UUID,testopt="1"):

    # Given a MD5 hash for a smiles, create a folder that stores epi_script, vega_script
    # and temperory result for all 3 models. A MD5 is needed because smiles contain special
    # We can't use smiles as folder identifier directly because:
    #   1. We are not certain that smiles have a limited length (<255)
    #   2. EPI might modifies smiles, making it harder to combine results
    # Thus, we apply MD5 function upon smiles string to obtain a 'unique' identifier

    md5 = hashlib.md5()
    md5.update(sys.argv[1])
    SMILES_MD5 = md5.hexdigest()
    
    print(epi,vega,test)
    print(testopt,type(testopt))
    # create temporary directory for current smiles
    # TEMP_DIR_PATH: runtime temp files
    # RESULT_JSON_FOLDER: 1. json for each model 2. combined json
    TEMP_DIR_PATH = os.path.join(DIR_PATH,'history',SMILES_MD5)
    make_dir_if_necessary(TEMP_DIR_PATH)
    RESULT_JSON_FOLDER = os.path.join(TEMP_DIR_PATH,'json')
    make_dir_if_necessary(RESULT_JSON_FOLDER)
    PATH_DICT = serialize_smiles_and_generate_scripts(smiles,TEMP_DIR_PATH,epi,vega,test)

    if epi:
        # run sikulix script to operate epi
        epi_time = time.time()
        # This command assumes a symlink to runsikulix have been created
        os.system("{0} -r {1}".format(os.path.join(DIR_PATH,"runsikulix"),
                                      os.path.join(DIR_PATH,"sikuli_scripts/epi_script.sikuli")))
        read_epi_result_toJson(EPI_SUITE_SAMPLE_RESULTS_JSON_FILEPATH)
        # os.system("python " +DIR_PATH+ "/episuite_file/parse_episuite.py")
        #os.system("rm "+dir_path+"/episuite_file/epibat.out")
        print("EPI used {} seconds to complete.".format(time.time()-epi_time))
        #pass
    else:
        currentepi=epinone
        currentepi["SMILES"]=smiles

        resultJsonObject=[currentepi]

        jsonOutputPath = os.path.normpath(os.path.join(DIR_PATH,"episuite_file/epibat.json"))
        with open(jsonOutputPath, "w") as outputFile:
            json.dump(resultJsonObject, outputFile,
                    sort_keys=True, indent= 4, separators=(',', ': '))
    
    #vega switch to turn on or off
    if vega:
        vega_time = time.time()
        #os.system("java -jar ./vega_file/VEGA_CMD/VEGA-CLI.jar -script ./vega_file/script_allModule_test")
        # vega_script_path = os.path.normpath(DIR_PATH + "/vega_file/script_shaoyi")
        # change_vega_script(vega_script_path)    
        java_command = os.path.normpath("{0}/vega_file/VEGA_CMD/VEGA-CLI.jar -script {1}"
                                        .format(DIR_PATH,PATH_DICT['VEGA_SCRIPT_PATH']))
        # change the hard coded path in vega script
        os.system("java -jar " + java_command)
        read_vega_result_toJSON(PATH_DICT["VEGA_RESULT_PATH"],os.path.join(RESULT_JSON_FOLDER,'vega_result.json'))
        # os.system("python " + os.path.normpath(DIR_PATH+ "/vega_file/parse_vega.py"))
        print("VEGA used {} seconds to complete.".format(time.time()-vega_time))
        #print("{} process used".format(cpu_count()))
    else:
        #create the empty vage component if switch is off
        currentvega=veganone
        currentvega["SMILES"]=smiles

        resultJsonObject = [currentvega]
        
        #output the result
        jsonOutputPath = os.path.normpath(os.path.join(DIR_PATH,"vega_file/result_test.json"))
        with open(jsonOutputPath, "w") as outputFile:
            json.dump(resultJsonObject, outputFile,
                    sort_keys=True, indent= 4, separators=(',', ': '))


    #test switch to turn on or off
    test_time = time.time()

    if test:
        # print(test_time)
        # print(os.system("rm -rf {0}".format(os.path.normpath(DIR_PATH+"/test_file/for_testing/temp_result2"))))
        # try:
        # command = "python " + os.path.normpath(DIR_PATH + "/test_file/Call_TEST.py ") \
        #           + str(testopt)
        TEST_batch_allEndpoints(PATH_DICT["TEST_SMILES_PATH"],PATH_DICT["TEST_RESULT_PATH"],testopt)
        # print(command)
        # # subprocess.call(command,shell=True)
        # try:
        #     print(subprocess.check_output(command, stderr=subprocess.STDOUT,shell=True))
        # except subprocess.CalledProcessError as exc:
        #     print(exc.output)
        print("{0} process used".format(cpu_count()))
        print("TEST used {0} seconds to complete.".format(time.time()-test_time))
    else:
    #create the empty test component if switch is off
        currenttest={}

        for a in TEST_Headers:
            currenttest[a]="N/A"

        currenttest["Smiles"]=smiles

        resultJsonObject = [currenttest]

        #output the result
        jsonOutputPath = os.path.normpath(os.path.join(DIR_PATH,"test_file/for_testing/temp_result2/test_results.json"))
        with open(jsonOutputPath, "w") as outputFile:
            # json.dump(resultJsonObject, outputFile,
            #         sort_keys=True, indent= 4, separators=(',', ': '))
            # To save space on clusters, we will condense json on one line
            json.dump(resultJsonObject, outputFile)

    #os.system("python parsing.py")

    #epiJSON = readJSON(EPI_SUITE_SAMPLE_RESULTS_JSON_FILEPATH)
    #vegaJSON = readJSON(VEGA_SAMPLE_RESULTS_JSON_FILEPATH)
    print(smiles)
    # testJSON = readJSON(TEST_SAMPLE_RESULTS_JSON_FILEPATH)
    # print(testJSON)
    # return testJSON

    #outputFilePath = DEFAULT_JSON_OUTPUT_FILEPATH
    #qsar_dict = parse(epiJSON,vegaJSON,testJSON,outputFilePath)
    #print(qsar_dict)
    #return qsar_dict

def save_test_result(outfile_path):
    # append output json of current smiles to the large json of 500 smiles
    with open(outfile_path,'a') as fp_out:
        with open(TEST_SAMPLE_RESULTS_JSON_FILEPATH,'r') as fp_in:
            fp_out.write(fp_in.read())

def save_epi_result(outfile_path):
    with open(outfile_path,'a') as fp_out:
        with open(EPI_SUITE_SAMPLE_RESULTS_JSON_FILEPATH,'r') as fp_in:
            fp_out.write(fp_in.read())        

if __name__ == '__main__':
    #switch("C(Cl)Cl",True,True,False)
    #testopt, 1:all,0:density and orat
    # print(len(sys.argv),sys.argv)    
    # UUID = sys.argv[5]


    EPI_SUITE_SAMPLE_RESULTS_JSON_FILEPATH = os.path.normpath(DIR_PATH + "/episuite_file/epibat.json")
    # VEGA_SAMPLE_RESULTS_JSON_FILEPATH = os.path.normpath(DIR_PATH + "/vega_file/result_test.json")
    # TEST_SAMPLE_RESULTS_JSON_FILEPATH =  os.path.normpath(os.path.join(DIR_PATH + "/test_file/for_testing/temp_result_" + UUID + "/test_results.json"))
    # DEFAULT_JSON_OUTPUT_FILEPATH =  os.path.normpath(DIR_PATH + "/QSAR_summay_sample.json")

    if len(sys.argv) == 7:
        switch(sys.argv[1],eval(sys.argv[2]),eval(sys.argv[3]),eval(sys.argv[4]),sys.argv[5],sys.argv[6])
        print("switch finished")

    if len(sys.argv) == 8:
        switch(sys.argv[1],eval(sys.argv[2]),eval(sys.argv[3]),eval(sys.argv[4]),sys.argv[5],sys.argv[6])
        print("switch finished")
        save_test_result(sys.argv[7])
    
    # production
    if len(sys.argv) == 5:
        switch(sys.argv[1],eval(sys.argv[2]),eval(sys.argv[3]),eval(sys.argv[4]),'uuid',"1")
        # python switch.py CC False True False
    #switch("CC",False,True,False,test_opt)
    #switch("C(Cl)Cl",True,False,False,test_opt)

