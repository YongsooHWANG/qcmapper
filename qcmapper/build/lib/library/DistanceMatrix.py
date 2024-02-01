# -*-coding:utf-8-*-

import collections
import math
import pandas
import itertools

from pprint import pprint

get_smaller = lambda a, b: a if a < b else b
get_bigger = lambda a, b: a if a > b else b


def calculate_swap_matrix(qchip_data, target_criterion):
	'''
		function to calculate a swap cost (time, fidelity) between two adjacent qubits
		args : cnot matrix of cnot links in quantum chip
	'''

	qchip_size = len(qchip_data.get("qubit_connectivity"))
	connectivity_matrix = qchip_data.get("qubit_connectivity")

	# cnot cost (time/error rate) 기반으로 swap cost 계산함
	# path 1: cnot a,b - cnot b,a - cnot a,b
	# path 2: cnot b,a - cnot a,b - cnot b,a
	# 상기 두 경로 비용을 비교해서, 더 좋은 값을 swap cost 로 결정

	if target_criterion == "time":
		matrix = [[0 if i==j else math.inf for i in range(qchip_size)] for j in range(qchip_size)]
		direc = [[None if i==j else None for i in range(qchip_size)] for j in range(qchip_size)]
		
		cnot_cost = qchip_data.get("cnot_gate_time")

		for i in connectivity_matrix:
			for j in connectivity_matrix[i]:

				cost_for = cnot_cost[i][j]*2 + cnot_cost[j][i]
				cost_back = cnot_cost[i][j] + cnot_cost[j][i]*2

				if cost_for < cost_back:
					matrix[i][j] = cost_for
					direc[i][j] = "{}>{}".format(i, j)

				else:
					matrix[i][j] = cost_back
					direc[i][j] = "{}>{}".format(j, i)
		
	
	elif target_criterion == "fidelity":
		matrix = [[1.0 if i==j else 0 for i in range(qchip_size)] for j in range(qchip_size)]
		direc = [[None if i==j else None for i in range(qchip_size)] for j in range(qchip_size)]
		
		cnot_cost = qchip_data.get("cnot_error_rate")

		for i in connectivity_matrix:
			for j in connectivity_matrix[i]:

				cost_for = (1.0-cnot_cost[i][j])**2 * (1-cnot_cost[j][i])
				cost_back = (1.0-cnot_cost[i][j]) * (1-cnot_cost[j][i])**2

				if cost_for < cost_back:
					matrix[i][j] = cost_back
					direc[i][j] = "{}>{}".format(j, i)

				else:
					matrix[i][j] = cost_for
					direc[i][j] = "{}>{}".format(i, j)
				
	return matrix, direc


def floyd_warshall(matrix, criterion):
	V = len(matrix)
	if not V:
		raise Exception("The input information of the graph is wrong")

	# path 찾기 위해 요구되는 
	pi = [[i if matrix[i][j] not in [0, math.inf] else -1 for j in range(V)] for i in range(V)]	

	if criterion in ["distance", "time"]:
		for k in range(V):
			next_matrix = [list(row) for row in matrix]
			for idx in itertools.product(range(V), range(V)):
				if matrix[idx[0]][idx[1]] > matrix[idx[0]][k] + matrix[k][idx[1]]:
					next_matrix[idx[0]][idx[1]] = matrix[idx[0]][k] + matrix[k][idx[1]]
					pi[idx[0]][idx[1]] = pi[k][idx[1]]

			matrix = next_matrix
					
	elif criterion == "fidelity":
		for k in range(V):
			next_matrix = [list(row) for row in matrix]
			for idx in itertools.product(range(V), range(V)):
				if matrix[idx[0]][idx[1]] < matrix[idx[0]][k] * matrix[k][idx[1]]:
					next_matrix[idx[0]][idx[1]] = matrix[idx[0]][k] * matrix[k][idx[1]]
					pi[idx[0]][idx[1]] = pi[k][idx[1]]

			matrix = next_matrix

	return matrix, pi


def find_paths(pi, qchip_size):
	"""
		function to find all the shortest path over multiple nodes
	"""
	paths = collections.defaultdict(list)
	
	for idx in itertools.product(range(qchip_size), range(qchip_size)):
		if idx[0] == idx[1]: continue
		paths[(idx[0], idx[1])].extend([idx[0], idx[1]])
		
		while True:
			current_node = paths[(idx[0], idx[1])][1]
			if pi[idx[0]][current_node] == idx[0]: break
			paths[(idx[0], idx[1])].insert(1, pi[idx[0]][current_node])

	return paths


def post_processing(matrix, qchip_data, path, target_criterion):
	"""
	"""
	qchip_size = len(qchip_data.get("qubit_connectivity"))

	if target_criterion == "time":
		cnot_cost = qchip_data.get("cnot_gate_time")
	
		for idx in itertools.product(range(qchip_size), range(qchip_size)):
			if idx[0] == idx[1]: continue
			final_edge = path[(idx[0], idx[1])][-2:]
			matrix[idx[0]][idx[1]] -= (cnot_cost[final_edge[0]][final_edge[1]] + cnot_cost[final_edge[1]][final_edge[0]])

	elif target_criterion == "fidelity":
		cnot_cost = qchip_data.get("cnot_error_rate")

		for idx in itertools.product(range(qchip_size), range(qchip_size)):
			if idx[0] == idx[1]: continue
			final_edge = path[(idx[0], idx[1])][-2:]
			matrix[idx[0]][idx[1]] /= (1-cnot_cost[final_edge[0]][final_edge[1]]) * (1-cnot_cost[final_edge[1]][final_edge[0]])

	return matrix


def generateDM(qchip_data, target_criterion="distance", **kwargs):
	connectivity_matrix = qchip_data.get("qubit_connectivity")
	qchip_size = len(qchip_data.get("qubit_connectivity"))
	
	if target_criterion == "distance":
		matrix = [[0 if i==j else math.inf for i in range(qchip_size)] for j in range(qchip_size)]
		for i in connectivity_matrix:
			for j in connectivity_matrix[i]:
				matrix[i][j] = 1
		
		answer_matrix, pi = floyd_warshall(matrix, target_criterion)
		# FT Circuit 회로 합성의 경우 path 가 필요없음
		flag_path = kwargs.get("flag_path")
		if flag_path:
			paths = find_paths(pi, qchip_size)
		else:
			paths = None
		
	else:
		# 먼저 cnot cost (time, fidelity) 기준으로 두 큐빗간 swap cost 계산함
		swap_cost_matrix, swap_direc = calculate_swap_matrix(qchip_data, target_criterion)

		# swap 비용 기준으로 최적 비용 행렬 계산
		# return
		# 	- matrix : 최적 비용 행렬
		# 	- pi : 최적 경로 계산에 필요한 
		matrix, pi = floyd_warshall(swap_cost_matrix, target_criterion)
		# 최적 비용 행렬 계산 하는 과정에서 함께 생성한 중간 노드 지정 행렬을 기반으로 경로 계산
		paths = find_paths(pi, qchip_size)
		# 경로의 마지막 구간을 기준으로 swap 비용 행렬에 대한 후 보정
		# 임의 두 큐빗간 cnot 비용 생성
		answer_matrix = post_processing(matrix, qchip_data, paths, target_criterion)

	# 임의 두 큐빗간 cnot 비용 및 경로 리턴
	return answer_matrix, paths
		
				
