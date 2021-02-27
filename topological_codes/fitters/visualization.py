import networkx as nx
import matplotlib.pyplot as plt


class VisualizationMixin:
    def graph_2D(self, G, edge_label):
        pos = nx.get_node_attributes(G, "pos")
        nx.draw_networkx(G, pos)
        labels = nx.get_edge_attributes(G, edge_label)
        labels = {x: round(y, 3) for (x, y) in labels.items()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
        plt.show()
