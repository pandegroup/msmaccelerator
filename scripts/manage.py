import sys
import IPython as ip
import numpy as np

from msmbuilder import Serializer
import msmbuilder.Trajectory


from msmaccelerator import Project
from msmaccelerator.database import Session
from msmaccelerator import models

project = Project(sys.argv[1])

ip.embed()

