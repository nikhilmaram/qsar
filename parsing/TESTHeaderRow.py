# -*- coding: utf-8 -*-
import os, csv


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



TEST_Endpoints = [
    "Fathead minnow LC50 (96 hr)",
    "Daphnia magna LC50 (48 hr)",
    "T. pyriformis IGC50 (48 hr)",
    "Oral rat LD50",
    "Bioaccumulation factor",
    "Developmental Toxicity",
    "Mutagenicity",
    "Normal boiling point",
    "Vapor pressure at 25°C",
    "Melting point",
    "Flash point",
    "Density",
    "Surface tension at 25°C",
    "Thermal conductivity at 25°C",
    "Viscosity at 25°C",
    "Water solubility at 25°C"
]


def getHeaders():
    headerList = []
    smilesFilesList = os.listdir("/home/yiting/Downloads/CLiCC/test run tox models from command line/test_chemicals/cmd_results/smile3")
    print sorted(smilesFilesList)
    for smilesFiles in sorted(smilesFilesList):
        if "result" in smilesFiles:
            with open("/home/yiting/Downloads/CLiCC/test run tox models from command line/test_chemicals/cmd_results/smile3/"+smilesFiles, 'rb') as csvfile:
                csvReader = csv.reader(csvfile, delimiter='\t')
                counter = 0
                for row in csvReader:
                    if counter == 0:
                        headerEntry = row[2:-1]
                        headerList.append(headerEntry)
                        print headerEntry
                    counter += 1
    print headerList
    print len(headerList)
    print len(TEST_Endpoints)
    return headerList

def combineHeaders(unitsHeader, endpoint):
    index = 0
    headerFinalList = []
    for endpointEntry in endpoint:
        for unit in unitsHeader[index]:
            tempHeader = endpointEntry + "  " + unit
            headerFinalList.append(tempHeader)
        index += 1
    print headerFinalList
    print len(headerFinalList)


if __name__ == '__main__':
    combineHeaders(getHeaders(), TEST_Endpoints)