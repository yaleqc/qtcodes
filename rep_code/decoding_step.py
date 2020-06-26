import sys, os

import networkx as nx
import matplotlib.pyplot as plt

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "qiskit-ignis",
            "qiskit",
            "ignis",
            "verification",
            "topological_codes",
        )
    ),
)
from qiskit.ignis.verification.topological_codes import RepetitionCode
from fitters import GraphDecoder

d = 4
T = 2
code = RepetitionCode(d, T)
decoder = GraphDecoder(code)
# graph = decoder.make_error_graph("1 1  01")
nx.draw_networkx(decoder.S)
plt.plot()
plt.show()
