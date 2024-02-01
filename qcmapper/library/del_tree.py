# -*-coding:utf-8-*-

import collections
from pprint import pprint

class Node:
    def __init__(self, data):
        self.data = data
        self.children = []

    def addChildren(self, child):
    	self.children.append(child)

    def setParent(self, parent):
        self.parent = parent

    def getData(self):
        return self.data

    def getChildren(self):
        return self.children

    def getParent(self):
        return self.parent

    def setData(self, item, value):
        self.data.update({item: value})


class Tree:
    def __init__(self, all_modules):
        self.length_module = 0
        self.leaf_nodes = collections.defaultdict(list)
        self.all_modules = all_modules
        self.visited_modules = collections.defaultdict(int)
        # self.current_module = root
        

    def get_leaf_nodes(self):
        return self.leaf_nodes


    def display_tree(self, node, **kwargs):
        _node_data = node.getData()
        _node_name = _node_data["name"]
        _node_level = _node_data["level"]

        for i in range(_node_level):
            print("\t\t", end='')

        print("({}) {} {} {} {}".format(_node_level, _node_name, _node_data["Leaf"], _node_data["LeafParent"],_node_data["argument_size"]))

        for i in node.getChildren():
            self.display_tree(i, length=len(_node_name))
            

    def make_tree(self, **kwargs):
        current_module = "main"
        level = 0
        node = Node({"name": current_module, 
                     "Leaf": self.all_modules[current_module]["Leaf"], 
                     "argument_size": self.all_modules[current_module]["argument_size"],
                     "LeafParent": self.all_modules[current_module]["LeafParent"],
                     "level": level})

        node.setParent(None)
        
        if "duplicate" in kwargs:
            self.connect_children(node, level, duplicate=kwargs["duplicate"])
        else:
            self.connect_children(node, level)

        return node


    def connect_children(self, current_node, level, **kwargs):
        current_module = current_node.getData()["name"]
        

        if len(current_module) > self.length_module:
            self.length_module = len(current_module)
        
        # 현재 모듈의 child 모듈
        _called_modules = []
        for command in self.all_modules[current_module]["QASM"]:
            function = command[0]
            if function in list(self.all_modules.keys()):
                _called_modules.append(function)

        if "duplicate" in kwargs:
            if kwargs["duplicate"]:
                _called_modules = set(_called_modules)

        # 피호출 모듈 (child) 이 없으면, 해당 모듈은 leaf 모듈에 포함됨
        if not len(_called_modules):
            self.leaf_nodes[level].append(current_node)
            # self.leaf_nodes.append(current_node)

        else:
            level+=1
            # 피 호출 모듈에 대해서...
            for i in _called_modules:
                # 피 호출 모듈에 대해 node 생성
                child_node = Node({"name": i, 
                                   "Leaf": self.all_modules[i]["Leaf"], 
                                   "argument_size": self.all_modules[i]["argument_size"],
                                   "LeafParent": self.all_modules[current_module]["LeafParent"],
                                   "level": level})
                
                # 현재 모듈 노드의 child 로 피호출 모듈 추가
                current_node.addChildren(child_node)
                # 피 호출 모듈의 부모로 현재 노드 추가
                child_node.setParent(current_node)
                # 피호출 모듈 노드를 기준으로 재귀적으로 다시 child 찾고, connect
                if "duplicate" in kwargs:
                    self.connect_children(child_node, level, duplicate=kwargs["duplicate"])


    def display_statistics(self):
        pprint(self.visited_modules)


    def check_node(self, current_node, data, **kwargs):
        # self.visited_modules[current_node.getData()["name"]] += 1

        for child in current_node.getChildren():
            _data = child.getData()
            _flag = True

            for key, value in data.items():
                if not (key in _data.keys() and data[key] == _data[key]):
                    _flag = False
                    break

            if not _flag:
                self.check_node(child, data)
            else:
                if "mode" in kwargs:
                    if kwargs["mode"] == "count":
                        self.solution+=1
                    elif kwargs["mode"] == "solution":
                        self.solution = child
                        return solution

        if "mode" in kwargs:
            if kwargs["mode"] == "count":
                return self.solution


    # data 의 구조는 dictionary 로 {key: value}
    def search_tree(self, root, data, **kwargs):
        self.solution = None

        if "mode" in kwargs:
            if kwargs["mode"] == "count":
                self.solution = 0

        self.check_node(root, data, mode=kwargs["mode"])

        return self.solution