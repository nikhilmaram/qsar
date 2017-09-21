import csv
import json
import sys,os
reload(sys)
sys.setdefaultencoding('utf8')

# The following code is not written by me(Shaoyi)
def convert(infile_path):
    with open(infile_path) as json_data:
        json_str = json.load(json_data)
        out_csv_path = os.path.join("csv",infile_path.split(os.sep)[-1].split('.')[-2] + ".csv")
        
        if not os.path.exists("csv"):
            os.makedirs("csv")
        f = csv.writer(open( out_csv_path, "wb+" ))

        keyarray=list()
        key2array=list()
        for key,value in json_str[0].items():
            keyarray.append(key)
            fl=0

            if key=="CAS" or key=="No." or key=="Smiles":
                key2array.append(" ")
            elif len(value)>1:
                for k,v in value.items():
                    key2array.append(k)
                    if fl<0:
                        keyarray.append(" ")
                    fl=fl-1

        # Write CSV Header, If you dont need that, remove this line
        f.writerow(keyarray)
        f.writerow(key2array)
        # for i in x:
        #     thisone=list()
        #     for j in keyarray:
        #         thisone.append(i[j])
        #     f.writerow(thisone)

        for i in json_str:
            thisone=list()
            for key,value in i.items():
                if key=="CAS" or key=="No." or key=="Smiles":
                    thisone.append(value)
                else:
                    if value!=None:
                        for k,v in value.items():
                            thisone.append(v)
                    else:
                        for k,v in json_str[0][key].items():
                            thisone.append(" ")
            f.writerow(thisone)





