import os
home = os.path.expanduser("~/")

os.chdir(home + "prog/Python/salome")

from salome import ImportComponentGUI
gg = ImportComponentGUI("GEOM")
from GeomTests import unit_tests
