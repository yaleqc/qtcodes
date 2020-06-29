# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 23:29:19 2020

@author: WillSun
"""


import networkx as nx
import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def generate_random_3Dgraph(n_nodes, radius, seed=None):
    if seed is not None:
        random.seed(seed)
    
    # Generate a dict of positions
    pos = {i: (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for i in range(n_nodes)}
    
    # Create random 3D network
    G = nx.random_geometric_graph(n_nodes, radius, pos=pos)
    return G

def generate_surface_3Dgraph(D, K, seed=None):
    if seed is not None:
        random.seed(seed)
    
    # Generate a dict of positions
    pos = {i: ((i%(D*D))//D, i%D, i//(D*D)) for i in range(D*D*K)}
    
    # Create random 3D network
    G = nx.random_geometric_graph(D*D*K, 2, pos=pos)
    return G

def network_plot_3D(G, angle, save=False):
    # Get node positions
    pos = nx.get_node_attributes(G, 'pos')
    
    # Get number of nodes
    n = G.number_of_nodes()
    # Get the maximum number of edges adjacent to a single node
    edge_max = max([G.degree(i) for i in range(n)])
    # Define color range proportional to number of edges adjacent to a single node
    colors = [plt.cm.plasma(G.degree(i)/edge_max) for i in range(n)] 
    # 3D network plot
    with plt.style.context(('ggplot')):
        
        fig = plt.figure(figsize=(10,7))
        ax = Axes3D(fig)
        
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        for key, value in pos.items():
            xi = value[0]
            yi = value[1]
            zi = value[2]
            
            # Scatter plot
            ax.scatter(xi, yi, zi, c=colors[key], s=20+20*G.degree(key), edgecolors='k', alpha=0.7)
            
            ax.text(xi, yi, zi, f'({zi}, {xi}, {yi})')
        
        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i,edge in enumerate(G.edges()):
            print(i, edge)
            x_1, x_2 = pos[edge[0]][0], pos[edge[1]][0]
            y_1, y_2 = pos[edge[0]][1], pos[edge[1]][1]
            z_1, z_2 = pos[edge[0]][2], pos[edge[1]][2]
            
            x = np.array((x_1, x_2))
            y = np.array((y_1, y_2))
            z = np.array((z_1, z_2))
            
            x_mid = (x_1 + x_2)/2
            y_mid = (y_1 + y_2)/2
            z_mid = (z_1 + z_2)/2
            ax.text(x_mid, y_mid, z_mid, i)
        
        # Plot the connecting lines
            ax.plot(x, y, z, c='black', alpha=0.5)
    
    # Set the initial view
    ax.view_init(30, angle)
    # Hide the axes
    ax.set_axis_off()
    
    if save is not False:
        plt.savefig("C:/Users/WillSun/Dropbox/Yale/IBM Hackathon/"+str(angle).zfill(3)+".png")
        plt.close('all')
    else:
        plt.show()
         
    return

n=200
#G = generate_random_3Dgraph(n_nodes=n, radius=0.25, seed=1)
G = generate_surface_3Dgraph(3, 2, seed=1)
for node in G.nodes():
    print(node)
'''
G = nx.Graph()
G.add_node("cool", pos=(0,1,2))
G.add_node("nice", pos=(2,1,-1))
G.add_edge("cool", "nice", weight=2)
for src, tgt in G.edges():
    print(src, tgt, G[src][tgt]['weight'])
'''
network_plot_3D(G, 30, save=False)