import numpy as np
import math
import scipy.special
import ctypes
from random import choices
from itertools import combinations, product
import networkx as nx
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

from pathlib import Path

path = Path(__file__).parent / "./libpypm.so"

_libpypm = ctypes.cdll.LoadLibrary(path)

_libpypm.infty.argtypes = None
_libpypm.infty.restype = ctypes.c_int
_libpypm.mwpm.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
_libpypm.mwpm.restype = None

# Integer that represents infinity (Blossom V algorithm). (Can be useful when converting float to int weights for MWPM)
INFTY = _libpypm.infty()


class XZZXDecoder:
    """
    Class to construct XZZX code on a rotated lattice. The current noise model is 
    one of Z-biased noise on data qubits only. 
    """
    def __init__(self, d, pz, eta):
        self.d = d
        self.pz = pz
        self.eta = eta # Define px=py=pz/eta, use eta <= 0 for infinite bias
        self.px = self._get_px()
        self.virtual = self._specify_virtual()
        self.S = self._specify_syndrome()


    def _get_px(self):
        if self.eta > 0:
            px = self.pz/self.eta
        else:
            px = 0

        return px


    def _specify_virtual(self):
        """
        Returns the virtual syndrome nodes for subgraph A and B separately in dictionary
        Returns: 
            virtual = {"A": [(-1, r, c), ... ], "B": [(-1, r, c), ...]}
        """
        virtual = {"A": [], "B": []}
        for i in range(0, self.d, 2):
            virtual["A"].append((-1, -0.5, i - 0.5))
            virtual["A"].append((-1, self.d - 0.5, i + 0.5))

            virtual["B"].append((-1, i + 0.5, -0.5))
            virtual["B"].append((-1, i - 0.5, self.d - 0.5))

        return virtual


    def _specify_syndrome(self):
        """
        Returns the real syndrome nodes for subgraph A and B separately in dictionary.
        S["A"] corresponds to S["Z"] for GraphDecoder
        Returns: 
            (syndrome) nodes = {"A": [(0, r, c), ... ], "B": [(0, r, c), ...]}
        """
        nodes = {"A": [], "B": []}
        for i in range(1, self.d, 2):
            nodes["A"].append((0, i - 0.5, self.d - 0.5))
            nodes["A"].append((0, i + 0.5, -0.5))

            nodes["B"].append((0, -0.5, i - 0.5))
            nodes["B"].append((0, self.d - 0.5, i + 0.5))

        for r in range(self.d - 1):
            for c in range(1, self.d, 2):
                if r % 2 == 0:
                    nodes["A"].append((0, r + 0.5, c - 0.5))
                    nodes["B"].append((0, r + 0.5, c + 0.5))
                else:
                    nodes["A"].append((0, r + 0.5, c + 0.5))
                    nodes["B"].append((0, r + 0.5, c - 0.5))

        return nodes  


    def get_error_syndrome(self):
        """
        Populate error on each data qubits, propagate errors to neighboring syndomes.
        Returns the actual logical error and the syndrome flips for error correction. 
        Returns:
            err_syndrome: Dictionary containing the nodes for subgraph A and B separately
            xL (int): The actual total logical X flips 
            zL (int): The actual total logical Z flips
        """
        pz = self.pz
        px = py = self.px
        dat_qubits = list(product(range(self.d), range(self.d)))
        anc_error = {"A": {}, "B": {}}
        err_syndrome = {"A": {}, "B": {}}
        xL = zL = 0

        # flipped = {"X": [], "Y": [], "Z": []}
        for pos in dat_qubits:
            # Determine error or not according to error probabilityF
            err_z = choices([0,1], [1-pz, pz])[0]
            err_x = choices([0,1], [1-px, px])[0]
            err_y = choices([0,1], [1-py, py])[0]

#             if err_x == 1:
#                 flipped["X"].append(pos)
#             if err_y == 1:
#                 flipped["Y"].append(pos)
#             if err_z == 1:
#                 flipped["Z"].append(pos)

            # Count actual logical error from data qubits
            if pos[0] == self.d - 1:
                if pos[1] % 2 == 0 and err_z == 1:
                    xL += 1
                elif pos[1] % 2 != 0 and err_x == 1:
                    xL += 1
                if err_y == 1:
                    xL += 1
            if pos[1] == self.d - 1:
                if pos[0] % 2 == 0 and err_x == 1:
                    zL += 1
                elif pos[0] % 2 != 0 and err_z == 1:
                    zL += 1
                if err_y == 1:
                    zL += 1                 

            left_up = (0, pos[0]-0.5, pos[1]-0.5)
            right_down = (0, pos[0]+0.5, pos[1]+0.5)
            left_down = (0, pos[0]+0.5, pos[1]-0.5)
            right_up = (0, pos[0]-0.5, pos[1]+0.5)

            # Propagate error from qubit to ancilla
            if left_up in self.S["A"]:
                try: anc_error["A"][left_up] += err_z + err_y
                except: anc_error["A"][left_up] = err_z + err_y
            elif left_up in self.S["B"]:
                try: anc_error["B"][left_up] += err_z + err_y
                except: anc_error["B"][left_up] = err_z + err_y    

            if right_down in self.S["A"]:
                try: anc_error["A"][right_down] += err_z + err_y
                except: anc_error["A"][right_down] = err_z + err_y
            elif right_down in self.S["B"]:
                try: anc_error["B"][right_down] += err_z + err_y
                except: anc_error["B"][right_down] = err_z + err_y

            if left_down in self.S["A"]:
                try: anc_error["A"][left_down] += err_x + err_y
                except: anc_error["A"][left_down] = err_x + err_y
            elif left_down in self.S["B"]:
                try: anc_error["B"][left_down] += err_x + err_y
                except: anc_error["B"][left_down] = err_x + err_y

            if right_up in self.S["A"]:
                try: anc_error["A"][right_up] += err_x + err_y
                except: anc_error["A"][right_up] = err_x + err_y
            elif right_up in self.S["B"]:
                try: anc_error["B"][right_up] += err_x + err_y
                except: anc_error["B"][right_up] = err_x + err_y

        # Determine error syndromes in subgraph A and B
        for subgraph in ["A", "B"]:
            err_syndrome[subgraph] = [
                pos
                for pos in anc_error[subgraph]
                if anc_error[subgraph][pos] % 2 != 0
            ]
        
        return err_syndrome, xL, zL   #, flipped


    def make_error_graph(self, err_syndrome, multi=False):
        """
        Construct error graph of subgraph A and B consisting of (node_a, node_n, distance)
        for MWPM. Multi-path summation optional. Distance calculated in self._path_summation(...),
        degeneracy calculated in self._degeneracy_cal.
        Args:
            err_syndrome: Dictionary containing the nodes for subgraph A and B separately
            multi: True for multi-path summation, False ignores such degeneracy
        Returns:
            error_graph: Dictionary containing the error graph ready for MWPM for subgraph A and B separately
        """
        error_graph = {"A": [], "B": []}
        for subgraph in ["A", "B"]:
            nodes = err_syndrome[subgraph] + self.virtual[subgraph] # All the nodes that can be matched
            if len(nodes) % 2 != 0:
                nodes.append((-1, self.d + 0.5, self.d + 0.5)) # Total nodes that can be matched must be even

            for source, target in combinations(nodes, 2):
                diff_r = target[1] - source[1]
                diff_c = target[2] - source[2]
                x, y = diff_c, -diff_r
                
                # Distance between virtual nodes is 0
                if source[0] < 0 and target[0] < 0:
                    distance = 0
                # The extra node only connects with virtual nodes
                elif source == (-1, self.d + 0.5, self.d + 0.5) or target == (-1, self.d + 0.5, self.d + 0.5):
                    continue
                elif y >= x:
                    if y >= -x:
                        distance = self._path_summation(source, target, (-1, 1), (-1, -1), subgraph, multi)
                    else:
                        distance = self._path_summation(source, target, (1, -1), (-1, -1), subgraph, multi)  
                else:
                    if y >= -x:
                        distance = self._path_summation(source, target, (-1, 1), (1, 1), subgraph, multi)
                    else:
                        distance = self._path_summation(source, target, (1, -1), (1, 1), subgraph, multi)

                if distance != None:
                    error_graph[subgraph].append((source, target, distance))

            # Scale the distance for MWPM function requires int(distance)
            max_abs_wt = max(abs(x[2]) for x in error_graph[subgraph]) 
            if max_abs_wt!=0:
                scaling=INFTY/10/max_abs_wt
            else:
                scaling=INFTY

            error_graph[subgraph] = [(x[0], x[1], int(x[2] * scaling)) for x in error_graph[subgraph]]
        
        return error_graph


    def _path_summation(self, source, target, vec_a, vec_b, subgraph, multi=False):
        """
        Calculates the distance between nodes. If multi-path summation is not implemented, 
        distance equals to edge number times edge weight (X/Z separately). If multi-path
        summation is implemented, degeneracy will be added as an extra term. 
        Args:
            source (tuple): The starting syndrome/virtual node
            target (tuple): The ending syndrome/virtual node
            vec_a (tuple): One of the two mutually orthogonal basis for describing the path
            vec_b (tuple): one of the two mutually orthogonal basis for describing the path
            subgraph ("A" or "B"): Specifying which subgraph the nodes are supposed to be on
            multi: True for multi-path summation, False ignores such degeneracy
        Returns:
            distance (float): distance between the source and target node for MWPM
        """
        # Solve for the edge number as two linear equations
        vec_lst = [vec_a, vec_b]
        RHS = np.array([target[1] - source[1], target[2] - source[2]])
        M = np.array(vec_lst).T
        result = np.linalg.inv(M).dot(RHS)
        result = result.astype(int)

        # Get edge weights
        weight = []
        for i in range(2):
            if vec_lst[i][0] != vec_lst[i][1]:
                # For pure Z noise, no possibility of X error
                if self.px == 0 and result[i] != 0:
                    return None
                elif self.px > 0:
                    weight.append(math.log1p(-2*self.px-self.pz) - math.log1p(2*self.px - 1))
                # No X_edge for pure Z noise, trivial case
                else:
                    weight.append(0)
            else:
                weight.append(math.log1p(-2*self.px-self.pz) - math.log1p(self.px + self.pz - 1))
        
        distance = result[0]*weight[0] + result[1]*weight[1]

        # Consider multi-path summation
        if multi:
            degeneracy = self._degeneracy_cal(result, vec_lst, source, subgraph)

            # TODO OPTIMIZE: CAN AT LEAST CUT ITERATION BY PASSING AN EXTRA ARGUMENT OF MATCHED VIRTUAL AND SYNDROME PAIR
            # If one of the nodes is virtual, check degeneracy of the other with all virtual nodes
            match = None
            if source[0] < 0:
                match = target
                matched = source
            elif target[0] < 0:
                match = source
                matched = target
            
            if match:
                vec_c = (-vec_a[0], -vec_a[1])
                vec_d = (-vec_b[0], -vec_b[1])
                vec_pairs = [(vec_a, vec_d), (vec_c, vec_b), (vec_c, vec_d)]
                
                self.virtual[subgraph].remove(matched)
                for virtual in self.virtual[subgraph]:
                    RHS = np.array([virtual[1] - match[1], virtual[2] - match[2]])
                    for pair in vec_pairs:
                        M = np.array(pair).T
                        result_n = np.linalg.inv(M).dot(RHS)
                        result_n = result_n.astype(int)

                        # In biased noise, Z_edges and X_edges are weighted differently
                        if self.eta != 1 and (result_n == result).all():
                            degeneracy += self._degeneracy_cal(result, pair, match, subgraph)
                            break
                        # In depolarizing noise, Z_edges and X_edges are weighted the same
                        elif self.eta == 1 and (result_n[0] + result_n[1]) == (result[0] + result[1]):
                            degeneracy += self._degeneracy_cal(result, pair, match, subgraph)
                            break
                        # -break- is used to avoid redundancy among new vec_pairs
                self.virtual[subgraph].append(matched)

            distance -= math.log(degeneracy)

        return distance


    def _degeneracy_cal(self, result, vec_lst, source, subgraph):
        """
        Calculates degeneracy given source and target node. Degeneracy is taken
        to be amongst shortest paths only. 
        Args:
            result (list): the separate number of basis vectors needed for the path
            vec_lst (list): the corresponding basis vectors
            source (tuple): the source node
            subgraph ("A" or "B"): Specifying which subgraph the nodes are supposed to be on
        """
        total = result[0] + result[1]
        edge_a = (vec_lst[0][0]*result[0] + source[1], vec_lst[0][1]*result[0] + source[2])
        edge_b = (vec_lst[1][0]*result[1] + source[1], vec_lst[1][1]*result[1] + source[2])
        edge_lst = [edge_a, edge_b] # Contains the two edge nodes besides source, target

        # Check for boundary conditions
        fault_node = []
        # Check both edge nodes
        for i in range(2):
            while not ((0,) + edge_lst[i] in self.S[subgraph]) and not ((-1,) + edge_lst[i] in self.virtual[subgraph]):
                if not fault_node:
                    fault_node.append(edge_lst[i])
                else:
                    tmp = []
                    for node in fault_node:
                        tmp.append((node[0] - vec_lst[i][0], node[1] - vec_lst[i][1]))
                        tmp.append((node[0] + vec_lst[(i+1)%2][0], node[1] + vec_lst[(i+1)%2][1]))
                    fault_node.extend(tmp)
                    fault_node = list(set(fault_node))
                edge_lst[i] = (edge_lst[i][0] - vec_lst[i][0], edge_lst[i][1] - vec_lst[i][1])

        if fault_node:
            vec_comb = list(product(*[vec_lst for i in range(total)]))
            vec_comb = [
                vecs
                for vecs in vec_comb
                if vecs.count(vec_lst[0]) == result[0]
            ] # Path with the correct number for each type of vectors
            degeneracy = 0
            for vec_path in vec_comb:
                node = (source[1], source[2])
                for vec in vec_path:
                    node = (node[0]+vec[0], node[1]+vec[1])
                    # Get rid of invalid path due to boundary conditions
                    if node in fault_node:
                        degeneracy -= 1
                        break
                degeneracy += 1
        else:
            degeneracy = int(scipy.special.binom(total, result[0]))

        return degeneracy


    def error_correct(self, matches):
        """
        Error correct according to syndromes, returned values are compared with
        actual logical error to determine the logical error rates. 
        Args:
            matches ([(node_a, node_b, edge), ...]): A list of all the matches from MWPM
        Retuns:
            xL (int): The calculated total logical X flips 
            zL (int): The calculated total logical Z flips
        """
        xL = zL = 0
        for match in matches:
            if match[0][0] < 0 and match[1][0] < 0:
                continue

            for node in match:
                if node in self.virtual["A"] and node[1] == self.d - 0.5:
                    xL += 1
                elif node in self.virtual["B"] and node[2] == self.d - 0.5:
                    zL += 1
            
        return xL, zL


    def _valid_syndrome(self, node, subgraph):
        """
        Supplementary function. Check whether a syndrome node is valid on our lattice.
        Args: 
        node ((r, c)): The node to be checked
        subgraph ("A" or "B"): Specifying which subgraph the node is supposed to be on
        Returns:
        Boolean
        """
        r, c = node[0], node[1]
        if subgraph == "A":
            if r > 0 and r < self.d - 1 and c > -1 and c < self.d:
                return True
            else: return False
        elif subgraph == "B":
            if c > 0 and c < self.d - 1 and r > -1 and r < self.d:
                return True
            else: return False


    def make_syndrome_graph(self):
        """
        Supplementary function. Make a complete NetworkX graph of the lattice with node position
        and edge weight. For visualization purposes. 
        Returns:
        node_graph: The complete NetworkX graph ready for plotting
        """
        node_graph = {"A": nx.Graph(), "B": nx.Graph()}

        start_nodes = {"A": (0.5, 0.5), "B": (0.5, 1.5)}
        for subgraph in ["A", "B"]:   
            start_node = start_nodes[subgraph]
            node_graph[subgraph].add_node(
                (0,) + start_node,
                virtual=0,
                pos=(start_node[1], -start_node[0]),
                time=0,
                pos_3D=(
                    start_node[1],
                    -start_node[0],
                    0,
                ),
            )
            self._populate_syndrome_graph(subgraph, (0,) + start_nodes[subgraph], [], node_graph)

        return node_graph


    def _populate_syndrome_graph(self, subgraph, curr_node, visited_nodes, node_graph):
        """
        Recursive function to populate syndrome subgraph at time 0 with syndrome_key A/B. The current_node
        is connected to neighboring nodes without revisiting a node.

        Args:
            subgraph ("A" or "B"): Which A/B syndrome subgraph these nodes are from.
            current_node ((0, r, c)): Current syndrome node to be connected with neighboring nodes.
            visited_nodes ([(0, r, c),]): List of syndrome nodes which have already been traver.
            node_graph (dictionary of two nx graphs): For appending nodes and edges to the complete graph

        Returns:
            None: function is to traverse the syndrome nodes and connect neighbors
        """
        neighbors = []
        r, c = curr_node[1], curr_node[2]
        neighbors.append((r - 1, c - 1, "Z"))
        neighbors.append((r + 1, c + 1, "Z"))
        neighbors.append((r - 1, c + 1, "X"))
        neighbors.append((r + 1, c - 1, "X"))

        normal_neighbors = [
            n 
            for n in neighbors
            if self._valid_syndrome(n, subgraph)
            and (0, n[0], n[1]) not in visited_nodes
        ]

        virtual_neighbors = [
            n
            for n in neighbors
            if (-1, n[0], n[1]) in self.virtual[subgraph]
            and (-1, n[0], n[1]) not in visited_nodes
        ]

        if not normal_neighbors and not virtual_neighbors:
            return
        
        for target in normal_neighbors:
            target_node = (0,) + target[:-1]
            if not node_graph[subgraph].has_node(target_node):
                node_graph[subgraph].add_node(
                    target_node,
                    virtual=0,
                    pos=(target[1], -target[0]),
                    time=0,
                    pos_3D=(target[1], -target[0], 0),
                )
            if target[2] == "Z":
                weight = round(math.log((1- 2*self.px - self.pz)/(self.pz + self.px)), 10)
            elif target[2] == "X":
                if self.px > 0:
                    weight = round(math.log((1- 2*self.px - self.pz)/(2*self.px)), 10)
                else:
                    weight = 0
            node_graph[subgraph].add_edge(
                curr_node, target_node, distance=weight
            )

        for target in virtual_neighbors:
            target_node = (-1,) + target[:-1]
            if not node_graph[subgraph].has_node(target_node):
                node_graph[subgraph].add_node(
                    target_node,
                    virtual=1,
                    pos=(target[1], -target[0]),
                    time=-1,
                    pos_3D=(target[1], -target[0], -1/2),
                )
            if target[2] == "Z":
                weight = round(math.log((1- 2*self.px - self.pz)/(self.pz + self.px)), 10)
            elif target[2] == "X":
                if self.px > 0:
                    weight = round(math.log((1- 2*self.px - self.pz)/(2*self.px)), 10)
                else:
                    weight = 0
            node_graph[subgraph].add_edge(
                curr_node, target_node, distance=weight
            )

        visited_nodes.append(curr_node)

        for target in normal_neighbors:
            self._populate_syndrome_graph(
                subgraph, (0,) + target[:-1], visited_nodes, node_graph
            )
        
        for target in virtual_neighbors:
            self._populate_syndrome_graph(
                subgraph, (-1,) + target[:-1], visited_nodes, node_graph
            )    


def mwpm_ids(edges):
    """
    Minimum Weight Perfect Matching using node ids (Blossom V algorithm).
    * Node ids are assumed to form a contiguous set of non-negative integers starting at zero, e.g.  {0, 1, ...}.
    * All nodes are assumed to participate in at least one edge.
    :param edges: Edges as [(node_id, node_id, weight), ...].
    :type edges: list of (int, int, int)
    :return: Set of matches as {(node_id, node_id), ...}. (Each tuple is sorted.)
    :rtype: set of (int, int)
    """
    node_ids = sorted(set(id for (id_a, id_b, _) in edges for id in (id_a, id_b)))
    n_nodes = len(node_ids)
    # Check node ids form contiguous set of non-negative integers starting at zero
    assert n_nodes == 0 or (node_ids[0] == 0 and node_ids[-1] == n_nodes - 1), (
        'Node ids are not a contiguous set of non-negative integers starting at zero.')
    n_edges = len(edges)
    nodes_a, nodes_b, weights = zip(*edges) if n_edges else ([], [], [])
    # Prepare arguments
    mates_array_type = ctypes.c_int * n_nodes
    edges_array_type = ctypes.c_int * n_edges
    mates_array = mates_array_type()
    # Call C interface
    _libpypm.mwpm(ctypes.c_int(n_nodes), mates_array, ctypes.c_int(n_edges),
                  edges_array_type(*nodes_a), edges_array_type(*nodes_b), edges_array_type(*weights))
    # Structure of mates: mates_array[i] = j means ith node matches jth node
    # Convert to more useful format: e.g. convert [1, 0, 3, 2] to {(0, 1), (2, 3)}
    mates = {tuple(sorted((a, b))) for a, b in enumerate(mates_array)}
    return mates


def mwpm(edges):
    """
    Convert to ID-sorted result for MWPM and back. 
    :param edges: Edges as [(node, node, weight), ...].
    :type edges: list of (object, object, int)
    :return: Set of matches as {(node, node), ...}.
    :rtype: set of (object, object)
    """
    nodes = list(set(node for (node_a, node_b, _) in edges for node in (node_a, node_b)))
    node_to_id = dict((n, i) for i, n in enumerate(nodes))
    edge_ids = [(node_to_id[node_a], node_to_id[node_b], weight) for node_a, node_b, weight in edges]
    mate_ids = mwpm_ids(edge_ids)
    mates = {(nodes[node_id_a], nodes[node_id_b]) for node_id_a, node_id_b in mate_ids}
    return mates


def graph_2D(G, edge_label):
    pos = nx.get_node_attributes(G, "pos")
    nx.draw_networkx(G, pos)
    labels = nx.get_edge_attributes(G, edge_label)
    labels = {x: round(y, 3) for (x, y) in labels.items()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    plt.show()

# %%
