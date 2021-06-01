import sys
import os

sys.path.insert(0, ".." + os.sep + ".." + os.sep + ".." + os.sep)
from benchmarking import RotatedSurfaceBatchAnalysis


batch_analysis = RotatedSurfaceBatchAnalysis(dirname="")
batch_analysis.plot()
