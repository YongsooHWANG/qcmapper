# -*-coding:utf-8-*-

# programmed by Y. Hwang (ETRI, Quantum Creative Research Lab.)
# 01 April 2016

# version 0.6

import numpy as np
import library.shortestpath as shortestpath

from pprint import pprint

qubit_connectivity_graph = None


def initialize_graph(qchip_data, **kwargs):
    """
        initialize graph for shortest path over nodes in graph
    """
    qubit_connectivity_graph = shortestpath.Graph()

    weight_type = kwargs.get("weight_type")
    for i in qchip_data: 
        qubit_connectivity_graph.add_vertex(i, weight_type=weight_type)
    
    if weight_type in [None, "depth"]:
        for k, v in qchip_data.items():
            for j in v:
                if k in qchip_data[j]:
                    qubit_connectivity_graph.add_edge(k, j, True)  # bi-directional
                else:
                    qubit_connectivity_graph.add_edge(k, j, False)  # one-directional

    elif weight_type in ["time", "fidelity"]:
        # set 큐빗 a 와 b 간의 swap 소요 시간
        for ctrl, trgt_list in qchip_data.items():
            for trgt, weight in trgt_list.items():
                qubit_connectivity_graph.add_edge(ctrl, trgt, weight, weight_type=weight_type)
                qubit_connectivity_graph.add_edge(trgt, ctrl, weight, weight_type=weight_type)
        
    return qubit_connectivity_graph


def findQubitConnectionPath(qchip_data, qubitA, qubitB, **kwargs):

    # initialize a graph adjacency matrix from the input relation
    qubit_connectivity_graph = initialize_graph(qchip_data, weight_type=kwargs.get("weight_type")) 

    source = qubit_connectivity_graph.get_vertex(qubitA)
    destination = qubit_connectivity_graph.get_vertex(qubitB)

    # find a path from the node i to the node j on the graph
    shortestpath.dijkstra(qubit_connectivity_graph, source, destination, weight_type=kwargs.get("weight_type"))

    target = qubit_connectivity_graph.get_vertex(qubitB)

    path = [qubitB]
    shortestpath.shortest(target, path)
    
    path.reverse()
    temp_path_list = []

    for k in range(len(path)-1):
        qubit_a = qubit_connectivity_graph.get_vertex(path[k])
        qubit_b = qubit_connectivity_graph.get_vertex(path[k+1])

        temp_path_list.append((path[k], path[k+1], qubit_a.get_weight(qubit_b)))

    return temp_path_list

    
def findALLConnectionPath(qubit_connectivity):

    qubit_connection_path = {}
    for i in qubit_connectivity.keys():
        for j in qubit_connectivity.keys():

            # initialize a graph adjacency matrix from the input relation
            qubit_connectivity_graph = initialize_graph(qubit_connectivity)
            # find a path from the node i to the node j on the graph
            shortestpath.dijkstra(qubit_connectivity_graph, qubit_connectivity_graph.get_vertex(i), qubit_connectivity_graph.get_vertex(j))
            target = qubit_connectivity_graph.get_vertex(j)

            path = [j]
            shortestpath.shortest(target, path)

            path.reverse()
            temp_path_list = []

            for k in range(len(path)-1):
                qubit_a = qubit_connectivity_graph.get_vertex(path[k])
                qubit_b = qubit_connectivity_graph.get_vertex(path[k+1])

                temp_path_list.append((path[k], path[k+1], qubit_a.get_weight(qubit_b)))

            qubit_connection_path[(i, j)] = temp_path_list

    return qubit_connection_path


