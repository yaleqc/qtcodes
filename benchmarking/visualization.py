# -*- coding: utf-8 -*-
"""
visualization script for benchmark data
"""

#%%
import sys
import os

sys.path.insert(0, ".." + os.sep)
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 14, "pdf.fonttype": 42, "ps.fonttype": 42})

from benchmarking.benchmarking_tools import SurfaceCodeBenchmarkingTool

import glob

data_files = glob.glob("*.npz")
benchmarking_tools = []
for file in data_files:
    benchmarking_tools.append(SurfaceCodeBenchmarkingTool(filename=file))
    benchmarking_tools[-1].load_data()

# Plotting
fig = plt.figure(figsize=(3.5, 2.5), dpi=200)
ax = fig.subplots()

for benchmarking_tool in benchmarking_tools:
    benchmarking_tool.plot_benchmark_data(
        fig=fig,
        ax=ax,
        label="d={},T={}".format(benchmarking_tool.d, benchmarking_tool.T),
    )

plt.legend(loc="lower right", prop={"size": 6})
ax.set_xlabel("Physical Error Rate", size=10)
ax.set_ylabel("Logical Error Rate", size=10)
ax.set_title("Comparison of Surface Codes", size=10)
fig.tight_layout()
plt.savefig("comparison.png")
plt.show()


# %%
