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


import sys
def StopSikuli(event):
    sys.exit()

Env.addHotkey(Key.F2, KEY_CTRL, StopSikuli)


App.focus("EPI Suite")
click("Selection_001.png")
click("Selection_002.png")
wait("Selection_003.png",10)
click("Selection_003-1.png")
wait("Selection_004.png", 10)
type("Selection_004-1.png", smiles_location)
click("Selection_005.png")
click("Selection_036.png")
click("Selection_037.png")
while exists("Selection_008.png",1):
    sleep(1)
    if wait("Selection_006.png", 3600):
        click("Selection_007.png")
        sleep(5)
    if exists("Selection_009.png"):
        break
    
wait("Selection_039.png", 2)
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
type(destination_folder)
sleep(1)
#click(Pattern("Selection_042.png").similar(0.60))
type(Key.ENTER)
if log:
    log.write("script reached end.\n\n")
    log.close()