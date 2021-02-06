import numpy as np
import math
import ctypes
from random import choices
from itertools import combinations, product
# import networkx as nx
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

_libpypm = ctypes.cdll.LoadLibrary("./libpypm.so")

_libpypm.infty.argtypes = None
_libpypm.infty.restype = ctypes.c_int
_libpypm.mwpm.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
_libpypm.mwpm.restype = None

# Integer that represents infinity (Blossom V algorithm). (Can be useful when converting float to int weights for MWPM)
INFTY = _libpypm.infty()


class XZZXDecoder:
    def __init__(self, d, pz, eta):
        self.d = d
        self.pz = pz
        self.eta = eta
        self.virtual = self._specify_virtual()
        self.S = self._specify_syndrome()


    def _specify_virtual(self):
        virtual = {"A": [], "B": []}
        for i in range(0, self.d, 2):
            virtual["A"].append((-1, -0.5, i - 0.5))
            virtual["A"].append((-1, self.d - 0.5, i + 0.5))

            virtual["B"].append((-1, i + 0.5, -0.5))
            virtual["B"].append((-1, i - 0.5, self.d - 0.5))

        return virtual


    def _specify_syndrome(self):
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


    def _valid_syndrome(self, node, subgraph):
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
        graph = {}
        node_graph = []
        start_nodes = {"A": (0.5, 0.5), "B": (0.5, 1.5)}
        for subgraph in ["A", "B"]:    
            self._populate_syndrome_graph(subgraph, (0,) + start_nodes[subgraph], [], node_graph)
            graph[subgraph] = node_graph
            node_graph = []
        return graph


    def _populate_syndrome_graph(self, subgraph, curr_node, visited_nodes, node_graph):
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
            if target[2] == "Z":
                weight = round(math.log((1-(self.pz/self.eta)-self.pz)/self.pz), 5)
            elif target[2] == "X":
                weight = round(math.log((1-(self.pz/self.eta)-self.pz)/(self.pz/self.eta)), 5)

            node_pair = [curr_node]
            node_pair.append((0,) + target[:-1])
            node_pair.append(weight)
            if node_pair not in node_graph and [node_pair[1], node_pair[0], node_pair[2]] not in node_graph:
                node_graph.append(node_pair)

        for target in virtual_neighbors:
            if target[2] == "Z":
                weight = round(math.log((1-(self.pz/self.eta)-self.pz)/self.pz), 5)
            elif target[2] == "X":
                weight = round(math.log((1-(self.pz/self.eta)-self.pz)/(self.pz/self.eta)), 5)

            node_pair = [curr_node]
            node_pair.append((-1,) + target[:-1])
            node_pair.append(weight)
            if node_pair not in node_graph and [node_pair[1], node_pair[0], node_pair[2]] not in node_graph:
                node_graph.append(node_pair)

        visited_nodes.append(curr_node)

        for target in normal_neighbors:
            self._populate_syndrome_graph(
                subgraph, (0,) + target[:-1], visited_nodes, node_graph
            )
        
        for target in virtual_neighbors:
            self._populate_syndrome_graph(
                subgraph, (-1,) + target[:-1], visited_nodes, node_graph
            )      


    def get_error_syndrome(self):
        pz = self.pz
        px = py = pz/self.eta
        dat_qubits = list(product(range(self.d), range(self.d)))
        anc_error = {"A": {}, "B": {}}
        err_syndrome = {"A": {}, "B": {}}

        for pos in dat_qubits:
            err_z = choices([0,1], [1-pz, pz])[0]
            err_x = choices([0,1], [1-px, px])[0]
            err_y = choices([0,1], [1-py, py])[0]

            left_up = (0, pos[0]-0.5, pos[1]-0.5)
            right_down = (0, pos[0]+0.5, pos[1]+0.5)
            left_down = (0, pos[0]+0.5, pos[1]-0.5)
            right_up = (0, pos[0]-0.5, pos[1]+0.5)

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

        for subgraph in ["A", "B"]:
            err_syndrome[subgraph] = [
                pos
                for pos in anc_error[subgraph]
                if anc_error[subgraph][pos] % 2 != 0
            ]
        
        return err_syndrome


    def make_error_graph(self, err_syndrome, multi=False):
        error_graph = {"A": [], "B": []}
        for subgraph in ["A", "B"]:
            nodes = err_syndrome[subgraph] + self.virtual[subgraph]

            for source, target in combinations(nodes, 2):
                diff_r = target[1] - source[1]
                diff_c = target[2] - source[2]
                x, y = diff_c, -diff_r
                
                if source[0] < 0 and target[0] < 0:
                    distance = 0
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

                error_graph[subgraph].append((source, target, distance))

            # NEED TO LOOK OUT FOR THE SIGN DEPENDING ON THE PATH SUMMATION RETURN
            max_abs_wt = max(abs(x[2]) for x in error_graph[subgraph]) 
            if max_abs_wt!=0:
                scaling=INFTY/10/max_abs_wt
            else:
                scaling=INFTY

            error_graph[subgraph] = [(x[0], x[1], int(x[2] * scaling)) for x in error_graph[subgraph]]
        
        return error_graph


    def _path_summation(self, source, target, vec_a, vec_b, subgraph, multi=False):
        vec_lst = [vec_a, vec_b]
        RHS = np.array([target[1] - source[1], target[2] - source[2]])
        M = np.array(vec_lst).T
        result = np.linalg.inv(M).dot(RHS)
        result = result.astype(int)

        print(source, target)
        print("vec_lst:", vec_lst)
        print("RHS:", RHS)
        print("result:", result)

        weight = []
        for vec in vec_lst:
            if vec[0] != vec[1]:
                weight.append(round(math.log((1-(self.pz/self.eta)-self.pz)/(self.pz/self.eta)), 5))
            else:
                weight.append(round(math.log((1-(self.pz/self.eta)-self.pz)/self.pz), 5))
        
        distance = result[0]*weight[0] + result[1]*weight[1]
        if multi:
            degeneracy = self._degeneracy_cal(result, vec_lst, source, subgraph)
            print("Before virtual:", degeneracy)

            # CAN OPTIMIZE BY PASSING AN EXTRA ARGUMENT OF MATCHED VIRTUAL AND SYNDROME PAIR
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
                        if (result_n == result).all():
                            degeneracy += self._degeneracy_cal(result, pair, match, subgraph)
                            print("Added virtual:". degeneracy)
                self.virtual[subgraph].append(matched)
            print("Final:", degeneracy)
            distance = round(distance - math.log(degeneracy), 5)

        return distance


    def _degeneracy_cal(self, result, vec_lst, source, subgraph):
        total = result[0] + result[1]
        degeneracy = math.comb(total, result[0])
        print("Before loss:", degeneracy)
        edge_a = (vec_lst[0][0]*result[0] + source[1], vec_lst[0][1]*result[0] + source[2])
        edge_b = (vec_lst[1][0]*result[1] + source[1], vec_lst[1][1]*result[1] + source[2])
        edge_lst = [edge_a, edge_b]
        print("edge_lst:", edge_lst)

        for i in range(2):
            cut = 0
            while not ((0,) + edge_lst[i] in self.S[subgraph]) and not ((-1,) + edge_lst[i] in self.virtual[subgraph]):
                cut += 1
                edge_lst[i] = (edge_lst[i][0] - vec_lst[i][0], edge_lst[i][1] - vec_lst[i][1])
                print("edge_lst[i]:", i, edge_lst[i], cut)
            if cut != 0:
                # THERE IS A BUG HERE
                degeneracy -= math.comb(total - result[i] + cut - 1, total - result[i])
                print("After loss:", degeneracy)
        
        return degeneracy


    def error_correct(self, matches):
        xL = zL = 0
        for match in matches:
            if match[0][0] < 0 and match[1][0] < 0:
                continue

            if match[0][1] == self.d - 0.5:
                zL += 1
            if match[1][1] == self.d - 0.5:
                zL += 1

            if match[0][2] == self.d - 0.5:
                xL += 1
            if match[1][2] == self.d - 0.5:
                xL += 1
            
        return xL, zL


def mwpm_ids(edges):
    """Minimum Weight Perfect Matching using node ids (Blossom V algorithm).

    Notes:

    * Node ids are assumed to form a contiguous set of non-negative integers starting at zero, e.g.  {0, 1, ...}.
    * All nodes are assumed to participate in at least one edge.

    :param edges: Edges as [(node_id, node_id, weight), ...].
    :type edges: list of (int, int, int)
    :return: Set of matches as {(node_id, node_id), ...}. (Each tuple is sorted.)
    :rtype: set of (int, int)
    """
    # extract and sort node ids
    node_ids = sorted(set(id for (id_a, id_b, _) in edges for id in (id_a, id_b)))
    # count n_nodes
    n_nodes = len(node_ids)
    # check node ids form contiguous set of non-negative integers starting at zero
    assert n_nodes == 0 or (node_ids[0] == 0 and node_ids[-1] == n_nodes - 1), (
        'Node ids are not a contiguous set of non-negative integers starting at zero.')
    # count n_edges
    n_edges = len(edges)
    # unzip edges
    nodes_a, nodes_b, weights = zip(*edges) if n_edges else ([], [], [])
    # prepare array types
    mates_array_type = ctypes.c_int * n_nodes
    edges_array_type = ctypes.c_int * n_edges
    # prepare empty mates
    mates_array = mates_array_type()
    # call C interface
    _libpypm.mwpm(ctypes.c_int(n_nodes), mates_array, ctypes.c_int(n_edges),
                  edges_array_type(*nodes_a), edges_array_type(*nodes_b), edges_array_type(*weights))
    # structure of mates: mates_array[i] = j means ith node matches jth node
    # convert to more useful format: e.g. convert [1, 0, 3, 2] to {(0, 1), (2, 3)}
    mates = {tuple(sorted((a, b))) for a, b in enumerate(mates_array)}
    return mates


def mwpm(edges):
    """Minimum Weight Perfect Matching using node objects (Blossom V algorithm).

    :param edges: Edges as [(node, node, weight), ...].
    :type edges: list of (object, object, int)
    :return: Set of matches as {(node, node), ...}.
    :rtype: set of (object, object)
    """
    # list of nodes without duplicates
    nodes = list(set(node for (node_a, node_b, _) in edges for node in (node_a, node_b)))
    # dict of node to id
    node_to_id = dict((n, i) for i, n in enumerate(nodes))
    # edges using ids
    edge_ids = [(node_to_id[node_a], node_to_id[node_b], weight) for node_a, node_b, weight in edges]
    # mwpm using ids
    mate_ids = mwpm_ids(edge_ids)
    # matches using objects
    mates = {(nodes[node_id_a], nodes[node_id_b]) for node_id_a, node_id_b in mate_ids}
    return mates


# def graph_2D(pairs):
#     G = nx.Graph()
#     for pair in pairs:
       
