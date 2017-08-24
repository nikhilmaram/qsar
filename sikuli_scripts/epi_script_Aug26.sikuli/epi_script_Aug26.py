import sys
import datetime
import os
import inspect

script_dir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe()))
)
# import locater script
path_array = script_dir.split(os.sep)
del path_array[-1]
base_script_dir = (os.sep).join(path_array)
sys.path.insert(0, base_script_dir)
from locater import get_locations
sys.path.pop(0)

locs = get_locations()

smiles_location = locs['smiles']
destination_folder = locs['results']
log_file = locs['log']

if log_file:
    try: 
       log = open(log_file, 'a')
       log.write("Started Epi Script. Log Loaded\n")
    except:
       log = False

App.focus("EPI Suite")
click("1454457927406.png")
click("Selection_032.png")
wait(Pattern("1454109316523.png").similar(0.80),10)
click(Pattern("1454109316523.png").similar(0.80).targetOffset(0,-6))
wait(Pattern("Selection_033.png").targetOffset(70,0), 10)
type(Pattern("Selection_033-1.png").targetOffset(70,0), "C:\users\yiting\Desktop\smiles.txt")
click("Selection_034.png")
click("Selection_036.png")
click("Selection_037.png")
wait("Selection_039.png", 3600)
click("Selection_038.png")
sleep(1)
App.focus("Note")
#wait(Pattern("Selection_046.png").similar(0.60), 10)
#click("Selection_044.png")
type(Key.ENTER)
#wait("Selection_045.png", 10)
sleep(3)
App.focus("Select an output file name for this batch data")
#type(destination_folder + r"\EPI_results.txt")
type("C:\users\yiting\Desktop\epi_results.txt")
sleep(1)
#click(Pattern("Selection_042.png").similar(0.6))
type(Key.ENTER)
if log:
    log.write("script reached end.\n\n")
    log.close()