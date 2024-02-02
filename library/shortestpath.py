# -*-coding:utf-8-*-

import heapq
import math

class Vertex:

    def __lt__(self, other):
        return self.distance < other.distance

    def __init__(self, node, **kwargs):
        self.id = node
        self.adjacent = {}
        # Set distance to infinity for all nodes

        self.distance = math.inf

        # Mark all nodes unvisited
        self.visited = False
        # Predecessor
        self.previous = None

    def add_neighbor(self, neighbor, weight):
        self.adjacent[neighbor] = weight
        # self.adjacent[neighbor] = int(weight)

    def get_connections(self):
        return self.adjacent.keys()

    def get_adjacent(self):
        return self.adjacent

    def get_id(self):
        return self.id

    def get_weight(self, neighbor):
        try:
            return self.adjacent[neighbor]["weight"]
        except:
            return self.adjacent[neighbor]

    def set_distance(self, dist):
        self.distance = dist
        # self.distance = int(dist)

    def get_distance(self):
        return self.distance
        # return int(self.distance)

    def set_previous(self, prev):
        self.previous = prev

    def set_visited(self):
        self.visited = True

    def __str__(self):
        return str(self.id) + ' adjacent: ' + str([x.id for x in self.adjacent])


class Graph:
    def __init__(self):
        self.vert_dict = {}
        self.num_vertices = 0

    def __iter__(self):
        return iter(self.vert_dict.values())

    def add_vertex(self, node, **kwargs):
        self.num_vertices += 1 
        self.vert_dict[node] = Vertex(node, weight_type=kwargs.get("weight_type"))
        
    def get_vertex(self, n):
        return self.vert_dict.get(n)

    def add_edge(self, frm, to, both, weight_type=None):
        # edge 삽입하기에 앞서, 해당 노드들이 graph 에 있는지 확인
        # if frm not in self.vert_dict:
        #     self.add_vertex(frm)

        # if to not in self.vert_dict:
        #     self.add_vertex(to)

        if weight_type in [None, "depth"]:
            if both:
                # CNOT direction
                self.vert_dict[frm].add_neighbor(self.vert_dict[to], 4)
                # CNOT reverse direction
                self.vert_dict[to].add_neighbor(self.vert_dict[frm], 4)

            else:
                # CNOT direction
                self.vert_dict[frm].add_neighbor(self.vert_dict[to], 1)
                # CNOT reverse direction
                self.vert_dict[to].add_neighbor(self.vert_dict[frm], 3)

        elif weight_type in ["time", "fidelity"]:
            self.vert_dict[frm].add_neighbor(self.vert_dict[to], both)
            self.vert_dict[to].add_neighbor(self.vert_dict[frm], both)

    def get_vertices(self):
        return self.vert_dict.keys()

    def set_previous(self, current):
        self.previous = current

    def get_previous(self, current):
        return self.previous

    def get_size(self):
        return len(self.vert_dict.keys())

    def display(self):
        print(" Qubit connectivity")
        print(" --------------------------------------------------------- ")

        for v in self:
            for w in v.get_connections():
                vid = v.get_id()
                wid = w.get_id()
                print(" {0} -> {1} ({2})".format(vid, wid, v.get_weight(w)))
        print(" --------------------------------------------------------- ")


def shortest(v, path):
    ''' make shortest path from v.previous'''
    if v.previous:
        path.append(v.previous.get_id())
        shortest(v.previous, path)
    
    return


def dijkstra(aGraph, start, target, **kwargs):
    # print '''Dijkstra's shortest path'''
    # Set the distance for the start node to zero

    weight_type = kwargs.get("weight_type")
    
    start.set_distance(0)

    # Put tuple pair into the priority queue
    unvisited_queue = [(v.get_distance(), v) for v in aGraph]
    
    heapq.heapify(unvisited_queue)

    # weight_type is None -> # number of swap gates
    if weight_type in [None, "depth", "time"]:
        while unvisited_queue:
            # while len(unvisited_queue):
            # Pops a vertex with the smallest distance
            uv = heapq.heappop(unvisited_queue)
            current = uv[1]
            current.set_visited()

            # for next in v.adjacent:
            for next in current.adjacent:
                # if visited, skip
                if next.visited: continue

                new_dist = current.get_distance() + current.get_weight(next)

                if new_dist < next.get_distance():
                    next.set_distance(new_dist)
                    next.set_previous(current)

            # Rebuild heap
            # 1. Pop every item
            while unvisited_queue: heapq.heappop(unvisited_queue)

            # 2. Put all vertices not visited into the queue
            unvisited_queue = [(v.get_distance(), v) for v in aGraph if not v.visited]

            heapq.heapify(unvisited_queue)


    elif weight_type == "fidelity":
        while unvisited_queue:
            # while len(unvisited_queue):
            # Pops a vertex with the smallest distance
            uv = heapq.heappop(unvisited_queue)
            current = uv[1]
            current.set_visited()

            # for next in v.adjacent:
            for next in current.adjacent:
                # if visited, skip
                if next.visited: continue

                new_dist = current.get_distance() * current.get_weight(next)
                
                if new_dist < next.get_distance():
                    next.set_distance(new_dist)
                    next.set_previous(current)

            # Rebuild heap
            # 1. Pop every item
            while unvisited_queue: heapq.heappop(unvisited_queue)

            # 2. Put all vertices not visited into the queue
            unvisited_queue = [(v.get_distance(), v) for v in aGraph if not v.visited]

            heapq.heapify(unvisited_queue)


def get_shortest_path(qubit_connectivity_graph, srcs=0, dest=0):
    """
        find the shortest weight path over qubits
    """

    dijkstra(qubit_connectivity_graph, qubit_connectivity_graph.get_vertex(srcs), qubit_connectivity_graph.get_vertex(dest))
    target = qubit_connectivity_graph.get_vertex(dest)

    path = [dest]
    shortest(target, path)

    print("{0} -> {1} : {2}".format(srcs, dest, path[::-1]))

    list_path = []
    reverse_path = path[::-1]

    for i in range(len(reverse_path)-1):
        list_path.append(reverse_path[i:i+1])

    return list_path
