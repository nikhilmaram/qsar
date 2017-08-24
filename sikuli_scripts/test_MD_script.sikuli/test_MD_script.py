import os
import sys
import inspect
print sys.argv

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

switchApp("T.E.S.T (Toxicity Estimation Software Tool)")
close_batch_button = exists("1444861040674.png")
if close_batch_button:
    click(close_batch_button)

click("1444850471970.png")
click("1444850488942.png")
type(smiles_location)
click("1444850546160.png")

click("1444850586082.png")
click("1444850597386.png")
text_field = find("1444851328866.png").left(100)
click(text_field)
type('a',KEY_CTRL)
type(destination_folder)
click("1444850700783.png")

endpoint_menu = find("1444851023626.png").right(100)
endpoint_list = endpoint_menu.above(50)



click(endpoint_menu)
wheel(endpoint_list, WHEEL_UP, 5)
wheel(endpoint_list, WHEEL_DOWN, 2)
MD_button = find("1455318981716.png")
if MD_button:
    click(MD_button)
click("1444850761650.png")
wait("1445639543422.png", 3600)
type('w',KEY_CTRL)