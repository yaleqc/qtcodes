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
import numpy as np

dir = "data_identity_noise/"
data_files = glob.glob(dir + "*.npz")
benchmarking_tools = []
for file in data_files:
    benchmarking_tools.append(SurfaceCodeBenchmarkingTool(filename=file))
    benchmarking_tools[-1].load_data()

sorted_indxs = np.argsort(
    np.array([benchmarking_tool.d for benchmarking_tool in benchmarking_tools])
)

for log_plot in [True, False]:
    # Plotting
    fig = plt.figure(figsize=(3.5, 2.5), dpi=200)
    ax = fig.subplots()
    for i in sorted_indxs:
        benchmarking_tool = benchmarking_tools[i]
        benchmarking_tool.plot_benchmark_data(
            fig=fig,
            ax=ax,
            label="d={},T={}".format(benchmarking_tool.d, benchmarking_tool.T),
            log=log_plot,
        )

    plt.plot(
        benchmarking_tools[sorted_indxs[0]].benchmark_data["noise"],
        0.5 * benchmarking_tools[sorted_indxs[0]].benchmark_data["noise"],
        "--",
        label="breakeven",
    )
    plt.legend(loc="lower right", prop={"size": 6})
    ax.set_xlabel("Physical Error Rate", size=10)
    ax.set_ylabel("Logical Error Rate", size=10)
    ax.set_title("Comparison of Surface Codes", size=10)
    fig.tight_layout()
    plt.savefig(dir + "comparison" + ("_log" if log_plot else "") + ".png")
    plt.show()


# %%
