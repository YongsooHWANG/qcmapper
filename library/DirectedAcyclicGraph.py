# -*-coding:utf-8-*-

import networkx as nx
import collections
import itertools
from pprint import pprint

from library.gatelist import *

def createDAG(list_qasm, **kwargs):

	DAG = nx.DiGraph()

	# 노드 증가할 때 마다 1씩 increment
	node_index = 0
	list_nodes_associated_with_qubit = collections.defaultdict(list)
	list_nodes_connection = []
	list_root_nodes = []
	list_qubits = []

	if "goal" in kwargs: goal = kwargs["goal"]
	else: goal = "dag"

	algorithm_qubits = kwargs.get("algorithm_qubits")

	for tokens in list_qasm:
		if tokens[0] in list_2q_gates + list_move:
			if tokens[0] in list_2q_rotations:
				angle, ctrl, trgt = tokens[1:4]
				DAG.add_node(node_index, gate=tokens[0], angle=angle, ctrl=ctrl, trgt=trgt, id=node_index)

			else:
				# DAG 노드 생성
				ctrl, trgt = tokens[1:3]
				DAG.add_node(node_index, gate=tokens[0], ctrl=ctrl, trgt=trgt, id=node_index)
				
			# current 노드가 DAG의 root 인지 확인하기 위한 flag
			flag_children = False

			# current 노드로 연결되는 부모 노드가 있는지 확인하는 단계
			# 공유하는 큐빗의 유무로 부모-자식 관계 결정됨
			if len(list_nodes_associated_with_qubit[ctrl]):
				flag_children = True
				parent_node_id = list_nodes_associated_with_qubit[ctrl][-1]
				list_nodes_connection.append((parent_node_id, node_index))

			# ctrl/trgt 큐빗을 공유하는 그래프 상의 마지막 노드에 현재 노드 추가
			list_nodes_associated_with_qubit[ctrl].append(node_index)

			# move 의 경우, trgt 가 "measurement_qubit" 이면, trgt 에 대해서는 무시함
			if tokens[0] not in [move] and trgt != "measurement_qubit":
				if len(list_nodes_associated_with_qubit[trgt]):
					flag_children = True
					parent_node_id = list_nodes_associated_with_qubit[trgt][-1]
					list_nodes_connection.append((parent_node_id, node_index))

				list_nodes_associated_with_qubit[trgt].append(node_index)

			# 부모 노드 없으면, 현재 노드가 root 노드에 해당함
			if not flag_children:
				list_root_nodes.append(DAG.nodes[node_index])
	
		# measurement
		elif tokens[0] in list_1q_gates:
			# measure first
			if tokens[0] in list_measure:
				trgt, *cbit = tokens[1::2]
				DAG.add_node(node_index, gate=tokens[0], trgt=trgt, cbit=cbit[0], id=node_index)
				
			# rotational gate
			elif tokens[0] in [u]:
				# u 게이트 (IBM QX 경우) 이면, 세 각도가 모두 입력됨
				angle_x, angle_y, angle_z, trgt = tokens[1:]
				DAG.add_node(node_index, gate=tokens[0], trgt=trgt,
										angle={"x":angle_x, "y":angle_y, "z": angle_z}, id=node_index)
			
			elif tokens[0] in list_1q_rotations:
				angle, trgt = tokens[1:]
				# QASM 구조가 Gate qubit angle 순 (이전 버전) 이면, Gate angle qubit 순으로 바꿔 해석함
				# 즉, qubit 과 angle 을 서로 교환함
				try:
					if type(eval(trgt)) == float: angle, trgt = trgt, angle
				except: pass
				finally:
					DAG.add_node(node_index, gate=tokens[0], trgt=trgt, angle=angle, id=node_index)

				
			# else Hadamard, Pauli..
			else:
				trgt = tokens[1]
				DAG.add_node(node_index, gate=tokens[0], trgt=trgt, id=node_index)	

			flag_children = False

			if len(list_nodes_associated_with_qubit[trgt]):
				flag_children = True
				parent_node_id = list_nodes_associated_with_qubit[trgt][-1]
				list_nodes_connection.append((parent_node_id, node_index))

			if not flag_children:
				list_root_nodes.append(DAG.nodes[node_index])

			list_nodes_associated_with_qubit[trgt].append(node_index)


		elif tokens[0] in [barrier]:
			list_qubits = tokens[1:]
			DAG.add_node(node_index, gate=g.str_barrier, trgt=list_qubits, id=node_index)

			# 큐빗 목록상의 개별 큐빗들에 대해서...
			# 해당 큐빗에 대해서 선행 inst 가 존재한다면.. current node 를 선행 노드에 연결
			# 없으면.. 없는 부분에 대해서는 root 로 넣어줘야 함

			temp_counter = 0
			for qubit in list_qubits:
				if len(list_nodes_associated_with_qubit[qubit]):
					parent_node_id = list_nodes_associated_with_qubit[qubit][-1]
					list_nodes_connection.append((parent_node_id, node_index))
					temp_counter += 1

				list_nodes_associated_with_qubit[qubit].append(node_index)

			if not temp_counter:
				list_root_nodes.append(DAG.nodes[node_index])


		elif tokens[0] in [barrier_all]:
			DAG.add_node(node_index, gate=tokens[0], id=node_index)

			temp_counter = 0
			for qubit in algorithm_qubits:
				if len(list_nodes_associated_with_qubit[qubit]):
					parent_node_id = list_nodes_associated_with_qubit[qubit][-1]
					# 이전 노드와 지금 노드를 연결시켜주는 부분
					list_nodes_connection.append((parent_node_id, node_index))
					temp_counter+=1

				list_nodes_associated_with_qubit[qubit].append(node_index)

			if not temp_counter:
				list_root_nodes.append(DAG.nodes[node_index])


		elif tokens[0] in ["Qubit"]: 
			trgt = tokens[1]
			DAG.add_node(node_index, gate=tokens[0], trgt=trgt, id=node_index)
			list_root_nodes.append(DAG.nodes[node_index])

		elif tokens[0] in ["Cbit"]: pass			
		elif tokens[0] in ["Release"]:
			# qubit array 이름은 tokens[1]
			# 해당 qubit array 의 모든 인덱스에 대해서, Release 하는 명령을 추가함
			
			# array 이름
			target_qubit_array_name = tokens[1]
			# 현재까지 사용된 큐빗들 목록 : list_nodes_associated_with_qubit.keys()
			
			if len(list_nodes_associated_with_qubit.keys()):
				for qubit in list_nodes_associated_with_qubit.keys():
					if target_qubit_array_name in qubit:
						DAG.add_node(node_index, gate="Release", trgt=qubit, id=node_index)
					
						flag_children = False
						if len(list_nodes_associated_with_qubit[qubit]):
							flag_children = True
							parent_node_id = list_nodes_associated_with_qubit[qubit][-1]
							list_nodes_connection.append((parent_node_id, node_index))
					
						if not flag_children:
							list_root_nodes.append(DAG.nodes[node_index])
						
						list_nodes_associated_with_qubit[qubit].append(node_index)
	
						node_index+=1
			else:
				DAG.add_node(node_index, gate="Release", trgt=target_qubit_array_name, id=node_index)
				list_root_nodes.append(DAG.nodes[node_index])
				node_index+=1

		else:
			# raise error.Error("Error Happened : Not recognized instruction -> {}".format(tokens))
			raise Exception("Error Happened : Not recognized instruction -> {}".format(tokens))

		node_index+=1

	# directed acyclic graph 그래프 연결...
	DAG.add_edges_from(list_nodes_connection)

	# nx.draw_shell(DAG, with_labels=True, font_weight='bold')
	# nx.draw(DAG, cmap=plt.get_cmap('jet'))
	# plt.show()

	return {"DAG": DAG, "roots": list_root_nodes}



def get_parent_from_node(DAG, node, depth):
	'''
		function to return ancestors in depth steps from the current node
	'''
	parents = collections.defaultdict(list)
	parents[0] = [node]

	for i in range(depth):
		if i in parents.keys() and len(parents[i]):
			for j in parents[i]:
				parents[i+1].extend(list(DAG.predecessors(j)))

	list_parents = list(itertools.chain.from_iterable(parents.values()))
	list_parents.remove(node)

	return list_parents



def get_children_from_node(DAG, node, depth):
	'''
		DAG의 특정 노드로 부터 depth 만큼 아래의 노드 return
		Extended Set E 생성하기 위해 필요한 함수
	'''

	children = collections.defaultdict(list)
	children[0] = [node["id"]]

	for i in range(depth):
		if i in children.keys() and len(children[i]):
			for j in children[i]:
				children[i+1].extend(list(DAG.successors(j)))

	list_children = list(itertools.chain.from_iterable(children.values()))
	list_children.remove(node["id"])

	return list_children
