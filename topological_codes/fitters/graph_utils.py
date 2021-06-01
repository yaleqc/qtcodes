"""
Common graph utilities.
"""

from collections import deque
import retworkx as rx
from typing import List, Optional


class GraphUtils:
    """
    Common Graph Utils
    """

    @staticmethod
    def num_shortest_paths(graph: rx.PyGraph, s: int) -> List[int]:
        """
        Returns a list of the number of distinct shortest paths
        from the source node s to every other node on graph.

        Arguments:
            graph (rx.PyGraph): the graph
            s (int): the index of the source node

        Returns:
            paths (List[int]):
                value: number of distinct shortest paths from source node s to target nodes
                index: target node graph index
        """
        n = len(graph.nodes())
        distance: List[Optional[int]] = [None] * n
        paths = [0] * n
        q = deque()
        q.append(s)
        distance[s] = 0
        paths[s] = 1

        visited = [False] * n
        while len(q) > 0:
            curr = q.popleft()
            for neighbor in graph.neighbors(curr):
                if not visited[neighbor]:
                    q.append(neighbor)
                    visited[neighbor] = True
                if (
                    distance[neighbor] is None
                    or distance[neighbor] > distance[curr] + 1
                ):
                    distance[neighbor] = distance[curr] + 1
                    paths[neighbor] = paths[curr]
                elif distance[neighbor] == distance[curr] + 1:
                    paths[neighbor] += paths[curr]

        return paths
