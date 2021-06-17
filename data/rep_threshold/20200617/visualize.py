import sys
import os

sys.path.insert(0, ".." + os.sep + ".." + os.sep + ".." + os.sep)
from benchmarking import TopologicalBatchAnalysis


batch_analysis = TopologicalBatchAnalysis(dirname="")
batch_analysis.plot()
