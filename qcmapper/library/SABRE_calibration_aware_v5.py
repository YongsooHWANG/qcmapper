# -*-coding:utf-8-*-
# time/fidelity tracking 하지 않음 --> 따라서 사전에 계산한 최적 cost와 비교하지 않음
import os
import re
import sys

import simplejson as json
import collections
import numpy as np
import math
import copy
import multiprocessing

import library.DirectedAcyclicGraph as DirectedAcyclicGraph
import library.DistanceMatrix as DM

import library.parse_qasm as parse_qasm
import qubitmapping.qubitmapping as qubitmapping

import library.formatconversion as formatconversion
from library.gatelist import *

import library.checkup as checkup
from progress.bar import Bar
import time

decay = 0.1
decay_increment = 1 + decay

extended_set_weight = 0.5

list_measured_qubit = []
lap_depth = 0


def calculate_NNC_cost(FL, SWAP_candidate_list, DM, MT, **kwargs):
	'''
		cost function based on nearest neighbor cost
	'''
	flag_calibration_aware_mapping = False

	calibration_type = kwargs.get("calibration_type")
	
	if calibration_type == "time": TM = kwargs.get("TM")
	elif calibration_type == "fidelity": EM = kwargs.get("EM")

	qchip_data = kwargs.get("qchip")

	cost = {}
	if calibration_type == "time":
		qchip_performance = qchip_data.get("cnot_gate_time")
		for swap in SWAP_candidate_list:
			temp_sum = 0

			# SWAP 효과 기반 비용 계산하기 위해서, 먼저 SWAP
			MT[swap[1]], MT[swap[2]] = MT[swap[2]], MT[swap[1]]
			for node in FL:
				# swap first
				if node["gate"] in [swap]: continue
				if node["gate"] in list_2q_gates:
					temp_sum += TM["matrix"][MT[node["ctrl"]]][MT[node["trgt"]]]
				
				elif node["gate"] in [move]:
					# move 경우, 큐빗 a 에서 큐빗 b 까지 이동하는 비용
					# move a, b 계산 방식
					# 기존의 TM 행렬에서는 CNOT a, b 비용이 포함되어 있음. 
					# 이 값을 기반으로 move a, b 비용을 계산하려면, path over a to b 의 마지막 edge 의 비용을 cnot -> swap 으로 하면 된다.
					temp_swap_cost = TM["matrix"][MT[node["ctrl"]]][node["trgt"]]
					final_edge = TM["path"][(MT[node["ctrl"]], node["trgt"])][-2:]

					if len(final_edge):
						temp_swap_cost += qchip_performance[final_edge[0]][final_edge[1]] +\
											qchip_performance[final_edge[1]][final_edge[0]]
					temp_sum += temp_swap_cost
			
			# 다음 SWAP candidate 에 대한 효과 평가하기 위해 SWAP 복원
			MT[swap[1]], MT[swap[2]] = MT[swap[2]], MT[swap[1]]
			cost[swap] = temp_sum
				
	elif calibration_type == "fidelity":
		qchip_performance = qchip_data.get("cnot_error_rate")
		# 개별 swap candidate 에 대해서, FL 내 게이트 실행을 위해 필요한 비용 합산 함
		for swap in SWAP_candidate_list:
			temp_sum = 1
			MT[swap[1]], MT[swap[2]] = MT[swap[2]], MT[swap[1]]

			for node in FL:
				# swap first
				if node["gate"] in [swap]: continue
				
				elif node["gate"] in list_2q_gates:
					temp_sum *= EM["matrix"][MT[node["ctrl"]]][MT[node["trgt"]]]
				
				elif node["gate"] in [move]:
					# move 경우, 큐빗 a 에서 큐빗 b 까지 이동하는 비용
					
					# move a, b 비용 (fidelity) 계산 방식
					# 기존 EM 행렬에서 CNOT a, b 실행 비용이 포함되어 있음
					# 이 값을 기준으로 move a, b 비용을 계산하려면, path (a->b) 상의 마지막 edge 비용을 CNOT -> SWAP 으로 하면 된다.
					temp_swap_cost = EM["matrix"][MT[node["ctrl"]]][node["trgt"]]
					final_edge = EM["path"][(MT[node["ctrl"]], node["trgt"])][-2:]

					if len(final_edge):
						temp_swap_cost *= (1-qchip_performance[final_edge[0]][final_edge[1]]) *\
											(1-qchip_performance[final_edge[1]][final_edge[0]])
					temp_sum *= temp_swap_cost

			MT[swap[1]], MT[swap[2]] = MT[swap[2]], MT[swap[1]]
			cost[swap] = temp_sum

	else:
		for swap in SWAP_candidate_list:
			temp_sum = 0
			MT[swap[1]], MT[swap[2]] = MT[swap[2]], MT[swap[1]]
			
			for node in FL:
				# swap first
				if node["gate"] in [swap]: continue
				elif node["gate"] in list_2q_gates:
					temp_sum += DM["matrix"][MT[node["ctrl"]]][MT[node["trgt"]]]
			
				elif node["gate"] in [move]:
					temp_sum += DM["matrix"][MT[node["ctrl"]]][node["trgt"]]

			MT[swap[1]], MT[swap[2]] = MT[swap[2]], MT[swap[1]]
			cost[swap] = temp_sum
	
	return cost



def calculate_LAP_cost(FL, SWAP, DAG, DM, MT, listDecay, **kwargs):
	'''
		cost function based on Look-Ahead Ability and Parallelism
	'''

	qchip_data = kwargs.get("qchip")

	calibration_type = kwargs.get("calibration_type")
	if calibration_type == "time": TM = kwargs.get("TM")
	elif calibration_type == "fidelity": EM = kwargs.get("EM")

	decay = max(listDecay[SWAP[0]], listDecay[SWAP[1]])
	cost = {}
	extended_set = []

	if calibration_type == "time":
		temp_cost_F = 0
		temp_cost_E = 0

		for node in FL:
			if nodes["gate"] in [swap]: 
				# 알고리즘 상에 SWAP 이 있는 경우, classically 큐빗 relabelling 하는 것으로 처리함...
				# 실제 양자 연산을 수행하지는 않음 -> 따라서 비용은 0
				continue
			
			elif node["gate"] in list_2q_gates:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = MT[node["trgt"]]
				temp_cost_F += TM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]

			elif node["gate"] in list_barrier: continue

			elif node["gate"] in [move]:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = node["trgt"]

				temp_swap_cost = TM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]
				final_edge = TM["path"][(associated_physical_qubit_ctrl, associated_physical_qubit_trgt)][-2:]

				if len(final_edge):
					temp_swap_cost += qchip_data.get("cnot_gate_time")[final_edge[0]][final_edge[1]] +\
										qchip_data.get("cnot_gate_time")[final_edge[1]][final_edge[0]]
				temp_cost_F += temp_swap_cost


			partial_extended_set = DirectedAcyclicGraph.get_children_from_node(DAG, node, lap_depth)
			extended_set.extend(partial_extended_set)


		for node in extended_set:
			if DAG.nodes[node]["gate"] in [swap]:
				# 알고리즘 상에 SWAP 이 있는 경우, classically 큐빗 relabelling 하는 것으로 처리함...
				# 실제 양자 연산을 수행하지는 않음 -> 따라서 비용은 0
				pass

			elif DAG.nodes[node]["gate"] in list_2q_gates:
				associated_physical_qubit_ctrl = MT[DAG.nodes[node]["ctrl"]]
				associated_physical_qubit_trgt = MT[DAG.nodes[node]["trgt"]]
				temp_cost_E += TM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]

			elif DAG.nodes[node]["gate"] in list_barrier: continue

			elif DAG.nodes[node]["gate"] in [move]:
				associated_physical_qubit_ctrl = MT[DAG.nodes[node]["ctrl"]]
				associated_physical_qubit_trgt = DAG.nodes[node]["trgt"]

				# extended set 에서 trgt 큐빗이 "measurement_qubit" 으로 되어 있으면, 
				# 앞으로 어떻게 될지 정확히 얘측하기 어렵기 때문에, ctrl 큐빗의 현재 위치에서 가장 먼 거리에 위치한 측정 가능 큐빗까지의 거리를 반영함 
				if associated_physical_qubit_trgt == "measurement_qubit": 
					# 현재 물리 큐빗을 기준으로 swap 비용 행렬 확인, 없으면 계산
					list_measurable_qubits = qchip_data.get("measurable_qubits")
					
					temp_cost = {qubit: TM["matrix"][associated_physical_qubit_ctrl][qubit] for qubit in list_measurable_qubits}
					for m_qubit in temp_cost:
						final_edge = TM["path"][(associated_physical_qubit_ctrl, m_qubit)][-2:]
						if len(final_edge):
							temp_cost[m_qubit] += qchip_data.get("cnot_gate_time")[final_edge[0]][final_edge[1]] +\
												qchip_data.get("cnot_gate_time")[final_edge[1]][final_edge[0]]
					
					# 현재 큐빗을 기준으로 측정가능 큐빗까지의 이동 시간 + 측정 소요시간 합산하고, 가장 큰 비용의 큐빗 결정
					distance = {qubit: temp_cost[qubit] + qchip_data.get("measure_time")[qubit] for qubit in list_measurable_qubits}
					associated_physical_qubit_trgt = max(distance, key=distance.get)
					
					# 해당 큐빗까지의 비용 합산
					temp_cost_E += distance[associated_physical_qubit_trgt]
				
				else:
					temp_swap_cost = TM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]
					final_edge = TM["path"][(associated_physical_qubit_ctrl, associated_physical_qubit_trgt)][-2:]

					if len(final_edge):
						temp_swap_cost += qchip_data.get("cnot_gate_time")[final_edge[0]][final_edge[1]] +\
											qchip_data.get("cnot_gate_time")[final_edge[1]][final_edge[0]]
					temp_cost_E += temp_swap_cost


	elif calibration_type == "fidelity":
		temp_cost_E = 1
		temp_cost_F = 1

		for node in FL:
			if node["gate"] in [swap]: continue
			
			elif node["gate"] in list_2q_gates:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = MT[node["trgt"]]
				temp_cost_F *= EM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]

			elif node["gate"] in list_barrier: continue

			elif node["gate"] in [move]:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = node["trgt"]

				temp_swap_cost = EM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]
				final_edge = EM["path"][(associated_physical_qubit_ctrl, associated_physical_qubit_trgt)][-2:]

				if len(final_edge):
					temp_swap_cost *= (1-qchip_data.get("cnot_error_rate")[final_edge[0]][final_edge[1]]) *\
										(1-qchip_data.get("cnot_error_rate")[final_edge[1]][final_edge[0]])
				temp_cost_F *= temp_swap_cost


			partial_extended_set = DirectedAcyclicGraph.get_children_from_node(DAG, node, lap_depth)
			extended_set.extend(partial_extended_set)

		for node in extended_set:
			if DAG.nodes[node]["gate"] in [swap]: pass

			if DAG.nodes[node]["gate"] in list_2q_gates:
				associated_physical_qubit_ctrl = MT[DAG.nodes[node]["ctrl"]]
				associated_physical_qubit_trgt = MT[DAG.nodes[node]["trgt"]]
				temp_cost_E *= EM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]

			elif DAG.nodes[node]["gate"] in list_barrier: continue

			elif DAG.nodes[node]["gate"] in [move]:
				associated_physical_qubit_ctrl = MT[DAG.nodes[node]["ctrl"]]
				associated_physical_qubit_trgt = DAG.nodes[node]["trgt"]

				if associated_physical_qubit_trgt == "measurement_qubit":
					list_measurable_qubits = qchip_data.get("measurable_qubits")

					temp_cost = {qubit: EM["matrix"][associated_physical_qubit_ctrl][qubit] for qubit in list_measurable_qubits}

					for m_qubit in temp_cost:
						final_edge = EM["path"][(associated_physical_qubit_ctrl, m_qubit)][-2:]
						if len(final_edge):
							temp_cost[m_qubit] *= (1-qchip_data.get("cnot_error_rate")[final_edge[0]][final_edge[1]]) *\
												(1-qchip_data.get("cnot_error_rate")[final_edge[1]][final_edge[0]])
					
					distance = {qubit : temp_cost[qubit] * (1-qchip_data.get("measure_error")[qubit]) for qubit in list_measurable_qubits}
					associated_physical_qubit_trgt = min(distance, key=distance.get)

					temp_cost_E *= distance[associated_physical_qubit_trgt]
				
				else:
					temp_swap_cost = EM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]
					final_edge = EM["path"][(associated_physical_qubit_ctrl, associated_physical_qubit_trgt)][-2:]

					if len(final_edge):
						temp_swap_cost *= (1-qchip_data.get("cnot_error_rate")[final_edge[0]][final_edge[1]]) *\
											(1-qchip_data.get("cnot_error_rate")[final_edge[1]][final_edge[0]])
					temp_cost_E *= temp_swap_cost

	else:
		temp_cost_E = 0
		temp_cost_F = 0

		for node in FL:
			if node["gate"] in [swap]: continue

			if node["gate"] in list_2q_gates:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = MT[node["trgt"]]
				temp_cost_F += DM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]

			elif node["gate"] in [move]:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = node["ctrl"]
				temp_cost_F += DM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]	

			elif node["gate"] in list_barrier: continue

			partial_extended_set = DirectedAcyclicGraph.get_children_from_node(DAG, node, lap_depth)
			extended_set.extend(partial_extended_set)

		for node in extended_set:
			if DAG.nodes[node]["gate"] in [swap]: continue

			if DAG.nodes[node]["gate"] in list_2q_gates:
				associated_physical_qubit_ctrl = MT[DAG.nodes[node]["ctrl"]]
				associated_physical_qubit_trgt = MT[DAG.nodes[node]["trgt"]]
				temp_cost_E += DM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]

			elif DAG.nodes[node]["gate"] in list_barrier: continue

			elif DAG.nodes[node]["gate"] in [move]:
				associated_physical_qubit_ctrl = MT[node["ctrl"]]
				associated_physical_qubit_trgt = node["trgt"]
				temp_cost_E += DM["matrix"][associated_physical_qubit_ctrl][associated_physical_qubit_trgt]


	cost = decay * float(temp_cost_F/len(FL))
	if len(extended_set):
		cost += extended_set_weight * float(temp_cost_E/len(extended_set))

	return cost


def SABRE(DAG, FL, MT, DM, qchip_data, **kwargs):
	'''
		SWAP-based bidirectional search 알고리즘 core 파트
		args:
			DM: distance matrix from qubit connectivity
			DAG: directed acyclic graph from algorithm
			FL: front layer from DAG
			MT: random qubit mapping table
	'''
	list_syscode_commands = []

	physical_connectivity = qchip_data.get("qubit_connectivity")

	# cost function type 설정: default -> nnc
	cost_function = kwargs.get("cost")
	if cost_function is None:
		cost_function = "nnc"

	# calibration aware mapping 수행 관련 flag
	flag_calibration_aware_mapping = True
	
	direction = kwargs.get("direction")
	# Time matrix or Fidelity matrix 가 주어지면, calibration aware 합성 수행 
	# --> flag_calibration_aware_mapping = True
	
	calibration_type = kwargs.get("calibration_type")
	if calibration_type == "time":
		TM = kwargs.get("TM")
		# DTM = qchip_data.get("cnot_gate_time")
	
	elif calibration_type == "fidelity":
		EM = kwargs.get("EM")
		# DEM = qchip_data.get("cnot_error_rate")

	else:
		flag_calibration_aware_mapping = False
	
	# flag for writing system code 
	# default : False
	flag_write_syscode = kwargs.get("write_syscode")
	if flag_write_syscode is None: flag_write_syscode = False

	# flag for using swap gate, default : False
	flag_swap = kwargs.get("swap")
	if flag_swap is None: flag_swap = False
	

	inverse_MT = {v: k for k, v in MT.items()}
	listDecay = collections.defaultdict(lambda: 0)

	list_executed_nodes = set([])
	previous_best_SWAP = None
	list_waiting_operations = collections.defaultdict(list)

	list_measure_qubits = qchip_data.get("measurable_qubits")
	
	if list_measure_qubits is None or not len(list_measure_qubits):
		flag_measurement_qubits_specified = False
	else:
		flag_measurement_qubits_specified = True

	# measurement 큐빗은 고정이기 때문에, 현재 위치만을 기준으로 측정큐빗 까지 시간 저장 
	# TM_measurement = collections.defaultdict(float)

	while len(FL):
		Execute_gate_list = []
		# 1. front layer 에 있는 게이트들 중 직접 수행 가능한 연산들을 Execute_gate_list 에 추가
		for node in FL:
			# one-qubit gate: 항상 실행 가능함
			if node["gate"] in list_measure:
				
				# SABRE 에서 first and second traversal 은 초기 큐빗 배치를 결정하기 위한 과정임
				# 따라서, forward direction 에서는 측정할 큐빗 (논리 큐빗) 의 현재 위치가 양자칩이 제공하는 측정 가능 큐빗 (measurement_qubits) 이면 측정함
				# backward direction 에서는 이러한 내용 무시 : 무조건 실행 가능
				if flag_measurement_qubits_specified:
					if direction == "forward":
						if MT[node["trgt"]] in list_measure_qubits:
							Execute_gate_list.append(node)
					else: 
						Execute_gate_list.append(node)
				else:
					Execute_gate_list.append(node)

			# one-qubit gates + qubit/cbit declaration
			elif node["gate"] in list_1q_gates + ["Qubit", "Cbit"]:
				Execute_gate_list.append(node)
			
			# move : place a logical state a physical qubit holds to another physical qubit
			elif node["gate"] in [move]:
				
				# move 명령 for measurement operation
				if node["trgt"] == "measurement_qubit":
					# 측정할 큐빗 (논리 큐빗)의 현재 물리적 위치
					current_qubit = MT[node["ctrl"]]

					# current_qubit 으로 부터 거리가 가까운 측정 가능 큐빗을 확인하기 위한 
					if calibration_type == "time":
						# 대상 큐빗의 측정 게이트 fidelity 까지 반영할 필요...
						# 현재 TM matrix 는 CNOT 시간만 반영하고 있으니까..
						# 해야 할 일 : TM 에서 마지막 경로의 CNOT cost 를 제외하고, SWAP cost 로 변경한 다음, measurement performance 반영할 것
						
						# 1. current_qubit -> qubit 까지 이동 비용이 필요하므로 (원래 처음 데이터 만들어 질 때의 값)
						temp_cost = {qubit: TM["matrix"][current_qubit][qubit] for qubit in list_measure_qubits}
						# final edge 가 cnot 만 반영되어 있으니까, SWAP 반영하기
						# 측정 큐빗에서 measurement 성능 반영하기
						for m_qubit in temp_cost:
							final_edge = TM["path"][(current_qubit, m_qubit)][-2:]
							if len(final_edge):
								temp_cost[m_qubit] += qchip_data.get("cnot_gate_time")[final_edge[0]][final_edge[1]] +\
														qchip_data.get("cnot_gate_time")[final_edge[1]][final_edge[0]]

						distance = {qubit: temp_cost[qubit] + qchip_data.get("measure_time")[qubit] for qubit in list_measure_qubits}
						nearest_qubit = min(distance, key=distance.get)
						
					elif calibration_type == "fidelity":
						
						temp_cost = {qubit: EM["matrix"][current_qubit][qubit] for qubit in list_measure_qubits}

						for m_qubit in temp_cost:
							final_edge = EM["path"][(current_qubit, m_qubit)][-2:]
							if len(final_edge):
								temp_cost[m_qubit] *= (1-qchip_data.get("cnot_error_rate")[final_edge[0]][final_edge[1]]) *\
														(1-qchip_data.get("cnot_error_rate")[final_edge[1]][final_edge[0]])

						distance = {qubit: temp_cost[qubit] * (1-qchip_data.get("measure_error")[qubit]) for qubit in list_measure_qubits}
						nearest_qubit = max(distance, key=distance.get)
					
					else:
						# distance 의 경우에는 현재 큐빗에서 측정 가능 큐빗까지의 거리만 확인하므로,
						# 특별히 추가 보정해주어야 할 값이 없다.
						distance = {qubit: DM[current_qubit][qubit] for qubit in list_measure_qubits}
						nearest_qubit = max(distacne, key=distance.get)

					# 비용 분석을 통해서, 측정 명령을 실행할 가장 가까운 물리 큐빗
					node["trgt"] = nearest_qubit
					if current_qubit == node["trgt"]: Execute_gate_list.append(node)

				else:
					if MT[node["ctrl"]] == node["trgt"]: Execute_gate_list.append(node)

			# barrier-All
			elif node["gate"] == barrier_all:
				if all(node["gate"] == barrier_all for node in FL):
					Execute_gate_list.append(node)

			# selective barrier
			elif node["gate"] == barrier:
				if len(FL) == 1 and FL[0] == node: 
					Execute_gate_list.append(node)

				else:
					# barrier 명령으로 blocked 된 큐빗 목록
					barrier_blocked_qubits = node["trgt"]

					working_qubits = []
					for checknode in [item for item in FL]:
						if checknode["gate"] == barrier_all:
							working_qubits.extend(list(MT.keys()))

						elif checknode["gate"] == barrier:
							working_qubits.extend(checknode["trgt"])

						elif checknode["gate"] in list_2q_gates:
							working_qubits.extend([checknode["ctrl"], checknode["trgt"]])

						else:
							working_qubits.append(checknode["trgt"])

					if not len(set(working_qubits).intersection(set(barrier_blocked_qubits))):
						Execute_gate_list.append(node)


			elif node["gate"] == swap:
				# swap 게이트의 경우.. 알고리즘에 포함된 SWAP (e.g., QFT)
				# 회로 합성에 두가지 방식이 존재함
				# 1. 실제 물리적으로 큐빗 위치 교환 : Move 연산 활용하자...
				#    Swap a,b : a -> pi(b), b -> pi(a)
				
				# 2. 큐빗 association 교환 : 양자 연산 없이...
				Execute_gate_list.append(node)

			# two-qubit gate: two-qubit 게이트의 대상이 되는 두 큐빗이 인접해 있거나, 
			# 또는 calibration data 기준으로 최적으로 인접한 상태일 때 실행 가능함
			else:
				logical_ctrl = node["ctrl"]
				logical_trgt = node["trgt"]

				mapped_physical_ctrl = MT[logical_ctrl]
				mapped_physical_trgt = MT[logical_trgt]

				# conditionA : mapping 결과 trgt 큐빗이 ctrl 큐빗에 인접한가?
				if mapped_physical_trgt in qchip_data.get("qubit_connectivity")[mapped_physical_ctrl]:
					Execute_gate_list.append(node)

		# 2. Execute_gate_list 가 empty 가 아니면, 
		#	1) front layer 에서 해당 게이트 삭제
		# 	2) DAG 에서 해당 게이트의 후속 게이트 확인
		#	3) 해당 후속 게이트의 큐빗 종속성 확인해서, 바로 실행 가능하면, front layer 에 추가
		#		(1) 큐빗 종속성을 확인하는 방법: 해당 게이트의 선행자가 여전히 front layer에 존재하는지 확인
		#		(2) 만약 선행자가 존재하면, 해당 게이트의 실행은 아직 불가능하고
		#		(3) 선행자가 없으면, 해당 게이트는 바로 실행 가능함
		
		if len(Execute_gate_list):
			# list_FL_nodes = set([node["id"] for node in FL])

			for gate in Execute_gate_list:
				if flag_write_syscode:
					# measurement
					if gate["gate"] in list_measure:
						# index of cbit is equal to qubit : cbit = trgt
						list_syscode_commands.append([gate["gate"], MT[gate["trgt"]], gate["cbit"]])
						
						# 측정된 큐빗은 측정 큐빗 목록에 삽입: 추후 해당 큐빗 사용 금지 위해 사용함
						list_measured_qubit.append(MT[gate["trgt"]])

					# one qubit gate 타입에 따라 qasm instruction 의 포맷이 약간씩 다름
					elif gate["gate"] in list_1q_gates:
						# rotational gate
						if gate["gate"] in list_1q_rotations:
							list_syscode_commands.append([gate["gate"], gate["angle"], MT[gate["trgt"]]])
						
						# other H, Pauli, T, Tdag gates
						else:
							list_syscode_commands.append([gate["gate"], MT[gate["trgt"]]])
					
					elif gate["gate"] == "Qubit":
						list_syscode_commands.append([gate["gate"], MT[gate["trgt"]]])
					
					elif gate["gate"] == "Cbit":
							list_syscode_commands.append([gate["gate"], gate["trgt"]])
					
					elif gate["gate"] == barrier_all:
						list_syscode_commands.append([node["gate"]])
						list_instructions = list_waiting_operations.get("all")
						if list_instructions is not None:
							FL.extend(list_instructions)
							list_waiting_operations["all"] = []

					elif gate["gate"] == barrier:
						list_qubits = [MT[qubit] for qubit in node["trgt"]]
						list_syscode_commands.append([node["gate"], list_qubits])
						key = gate["id"]

						list_instructions = list_waiting_operations.get(key)
						if list_instructions is not None:
							FL.extend(list_instructions)
							del list_waiting_operations[key]

					elif gate["gate"] in [move]: pass

					elif gate["gate"] in list_2q_rotations:	
						list_syscode_commands.append([gate["gate"], gate["angle"], MT[gate["ctrl"]], MT[gate["trgt"]]])	

					elif gate["gate"] in list_2q_gates:
						list_syscode_commands.append([gate["gate"], MT[gate["ctrl"]], MT[gate["trgt"]]])

						if gate["gate"] == swap:
							MT[gate["ctrl"]], MT[gate["trgt"]] = MT[gate["trgt"]], MT[gate["ctrl"]]


				# 실행할 수 있는 system code 명령은 FL 에서 삭제
				FL.remove(gate)
				list_executed_nodes.add(gate["id"])

				# FL 에 존재하는 barrier 노드 : not barrier-all
				list_barrier_nodes_in_FL = [temp_node for temp_node in FL if temp_node["gate"] == barrier]
				# barrier 노드에 의해서 blocked 된 노드들
				list_barrier_blocked_qubits = []
				# barrier 노드에 의해 blocked 된 노드들을 기준으로 barrier 노드의 id 를 지정
				inverse_barrier_nodes = collections.defaultdict(list)

				for barrier_node in list_barrier_nodes_in_FL:
					list_barrier_blocked_qubits.extend(barrier_node["trgt"])
					for qubit in barrier_node["trgt"]:
						inverse_barrier_nodes[qubit].append(barrier_node["id"])
				list_barrier_blocked_qubits = set(list_barrier_blocked_qubits)
				
				# 삭제할 노드의 후속 노드들에 대해서....
				# 노드의 선행자들이 모두 실행완료되었으면.. 해당 노드가 실행 가능하다
				
				for j in DAG.successors(gate["id"]):
					flag_j_append = False

					# 조상들과 FL 내 원소들간의 비교
					# 단순히 현재 노드의 부모 노드만을 고려하면 안된다.
					ancestors = set(DAG.predecessors(j))

					# 후행자 j 의 선배들이 모두 실행된 상태라면... -> 기본적으로 후행자 j 를 FL 로 당겨올 수 있음
					if ancestors.issubset(list_executed_nodes):
						# 만약 FL 에 barrier_all 이 있으면.. FL 로 당겨오지 못하고
						# 임시 저장소에 저장...
						if barrier_all in [temp_node["gate"] for temp_node in FL]:
							list_waiting_operations["all"].append(DAG.nodes[j])


						# FL 에 barrier 노드들이 존재하면... 후행자 j 의 게이트/대상 큐빗에 따라 결정
						elif len(list_barrier_nodes_in_FL):
							# 후행자 j 노드의 양자 게이트가 2-큐비트 게이트이면...
							
							if DAG.nodes[j]["gate"] in list_2q_gates:
								# 혹시나 2-큐비트 게이트의 ctrl/trgt 큐비트가 barrier blocked 이면...
								if DAG.nodes[j]["ctrl"] in list_barrier_blocked_qubits:
									for barrier_node_id in inverse_barrier_nodes[DAG.nodes[j]["ctrl"]]:
										list_waiting_operations[barrier_node_id].append(DAG.nodes[j])

								elif DAG.nodes[j]["trgt"] in list_barrier_blocked_qubits:
									for barrier_node_id in inverse_barrier_nodes[DAG.nodes[j]["trgt"]]:
										list_waiting_operations[barrier_node_id].append(DAG.nodes[j])

								# 2-큐비트 게이트의 ctrl/trgt 큐비트가 barrier blocked 가 아니면..
								# 그냥 FL 로 당겨오면 됨
								else:
									FL.append(DAG.nodes[j])

							# 후행자 j 노드의 양자 게이트가 2-큐비트 게이트가 아니면..
							else:
								FL.append(DAG.nodes[j])
						
						# FL 에 barrier 노드가 없으면.. FL 로 당겨오면 됨
						else:
							FL.append(DAG.nodes[j])


		# 3. Execute_gate_list 가 empty 이면
		#	1) score 자료구조 초기화 (cost value 저장)
		#	2) 대상 SWAP 게이트들을 -> SWAP_candidate_list 에 추가
		#	3) SWAP_candidate_list 의 각 SWAP 에 대해서,	
		#		(1) 해당 SWAP을 반영한 일시적 매핑 확보
		#		(2) 해당 매핑 기준 cost value 계산
		#	4) minimal cost value 갖는 SWAP 선택
		#	5) 해당 SWAP 을 실제로 선택해서, 매핑 업데이트
		else:
			cost = {}
			SWAP_candidate_list = []

			# make SWAP candidate list
			for node in FL:
				if node["gate"] in list_barrier: continue

				if node["gate"] in list_2q_gates:
					logical_ctrl = node["ctrl"]
					logical_trgt = node["trgt"]

					mapped_physical_ctrl = MT[logical_ctrl]
					mapped_physical_trgt = MT[logical_trgt]

				elif node["gate"] in [move]:
					logical_ctrl = node["ctrl"]
					logical_trgt = node["trgt"]

					mapped_physical_ctrl = MT[logical_ctrl]
					mapped_physical_trgt = logical_trgt

				for neighbor in physical_connectivity[mapped_physical_ctrl]:
					SWAP_candidate_list.append(((logical_ctrl, logical_trgt), logical_ctrl, inverse_MT[neighbor]))

				# for move gate, the trgt qubit is the destination for moving
				# therefore, it must not be moved to other place
				if node["gate"] != move:
					for neighbor in physical_connectivity[mapped_physical_trgt]:
						SWAP_candidate_list.append(((logical_ctrl, logical_trgt), logical_trgt, inverse_MT[neighbor]))

			if not len(SWAP_candidate_list): continue

			# evaluate all the SWAP gates in the candidate list
			if cost_function == "nnc":
				if calibration_type == "time":
					cost = calculate_NNC_cost(FL, SWAP_candidate_list, DM, MT, 
											calibration_type=calibration_type,
											TM=TM, 
											qchip=qchip_data)

				elif calibration_type == "fidelity":
					cost = calculate_NNC_cost(FL, SWAP_candidate_list, DM, MT, 
											calibration_type=calibration_type,
											EM=EM, 
											qchip=qchip_data)

			elif cost_function == "lap":
				if calibration_type == "time":
					for SWAP in SWAP_candidate_list:
						MT[SWAP[1]], MT[SWAP[2]] = MT[SWAP[2]], MT[SWAP[1]]
						listDecay[SWAP[1]] += decay_increment
						listDecay[SWAP[2]] += decay_increment

						cost[SWAP] = calculate_LAP_cost(FL, SWAP, DAG, DM, MT, listDecay, 
												calibration_type=calibration_type,
												TM=TM, 
												qchip=qchip_data)

						listDecay[SWAP[1]] -= decay_increment
						listDecay[SWAP[2]] -= decay_increment
						MT[SWAP[1]], MT[SWAP[2]] = MT[SWAP[2]], MT[SWAP[1]]
				
				elif calibration_type == "fidelity":
					for SWAP in SWAP_candidate_list:
						MT[SWAP[1]], MT[SWAP[2]] = MT[SWAP[2]], MT[SWAP[1]]
						listDecay[SWAP[1]] += decay_increment
						listDecay[SWAP[2]] += decay_increment
						
						cost[SWAP] = calculate_LAP_cost(FL, SWAP, DAG, DM, MT, listDecay, 
											calibration_type=calibration_type,
											EM=EM, 
											# EM_swap=EM_swap, 
											qchip=qchip_data)

						listDecay[SWAP[1]] -= decay_increment
						listDecay[SWAP[2]] -= decay_increment
						MT[SWAP[1]], MT[SWAP[2]] = MT[SWAP[2]], MT[SWAP[1]]

			if calibration_type in ["time", "distance"]:
				best_SWAP = min(cost, key=cost.get)

			elif calibration_type in ["fidelity"]:
				best_SWAP = max(cost, key=cost.get)
			
			if len(cost)>1:
				while True:
					if best_SWAP != previous_best_SWAP: break
					else:
						del cost[best_SWAP]
						best_SWAP = min(cost, key=cost.get)

			listDecay[best_SWAP[1]] += (1+decay)
			listDecay[best_SWAP[2]] += (1+decay)

			associated_physical_qubit_ctrl = MT[best_SWAP[1]]
			associated_physical_qubit_trgt = MT[best_SWAP[2]]

			# swapping qubit mapping table
			MT[best_SWAP[1]], MT[best_SWAP[2]] = MT[best_SWAP[2]], MT[best_SWAP[1]]
			inverse_MT[associated_physical_qubit_ctrl], inverse_MT[associated_physical_qubit_trgt] =\
				inverse_MT[associated_physical_qubit_trgt], inverse_MT[associated_physical_qubit_ctrl]

			previous_best_SWAP = best_SWAP
			
			if flag_write_syscode:
				if flag_swap:
					# update qubit mapping by performing the best SWAP
					list_syscode_commands.append([swap, MT[best_SWAP[1]], MT[best_SWAP[2]]])
				
				else:
					if calibration_type == "time":
						# SWAP -> CNOT 어떻게 분해하는 게 좋은지...
						cost_forward = TM["matrix"][MT[best_SWAP[1]]][MT[best_SWAP[2]]]*2 + TM["matrix"][MT[best_SWAP[2]]][MT[best_SWAP[1]]]
						cost_backward = TM["matrix"][MT[best_SWAP[2]]][MT[best_SWAP[1]]]*2 + TM["matrix"][MT[best_SWAP[1]]][MT[best_SWAP[2]]]

						if cost_forward <= cost_backward: swap_decompose_direction = "forward"
						else: swap_decompose_direction = "backward"

					elif calibration_type == "fidelity":
						cost_forward = EM["matrix"][MT[best_SWAP[1]]][MT[best_SWAP[2]]]**2 + EM["matrix"][MT[best_SWAP[2]]][MT[best_SWAP[1]]]
						cost_backward = EM["matrix"][MT[best_SWAP[2]]][MT[best_SWAP[1]]]**2 + EM["matrix"][MT[best_SWAP[1]]][MT[best_SWAP[2]]]

						if cost_forward >= cost_backward: swap_decompose_direction = "forward"
						else: swap_decompose_direction = "backward"

					else:
						swap_decompose_direction = "forward"
					
					if swap_decompose_direction == "forward":
						list_syscode_commands.extend([[cnot, MT[best_SWAP[1]], MT[best_SWAP[2]]],
													  [cnot, MT[best_SWAP[2]], MT[best_SWAP[1]]],
													  [cnot, MT[best_SWAP[1]], MT[best_SWAP[2]]]])
					else:
						list_syscode_commands.extend([[cnot, MT[best_SWAP[2]], MT[best_SWAP[1]]],
													  [cnot, MT[best_SWAP[1]], MT[best_SWAP[2]]],
													  [cnot, MT[best_SWAP[2]], MT[best_SWAP[1]]]])


	if flag_write_syscode:
		return list_syscode_commands

	else:
		return True


def evaluate_syscode(system_code, **kwargs):
	''' 
		cnot gate 시간 기준으로 본 시스템 코드의 실행 시간 평가
	'''
	p = re.compile(r"[\{a-zA-Z0-9_.*/\->\+}]+")

	if "criterion" in kwargs:
		performance_criterion = kwargs["criterion"]

	if "qubit_mapping" in kwargs:
		inverse_qubit_mapping = {int(v): k for k, v in kwargs["qubit_mapping"].items()}

	physical_measure_time = None
	physical_measure_error = None
	gate_time = None
	error_rate = None

	flag_single_qubit = kwargs.get("single_qubit")
	
	qchip_data = kwargs.get("qchip_data")
	if qchip_data is None:
		raise Exception("qchip data is not provided.")
	

	if performance_criterion == "depth":
		import depth_analysis
		depth = depth_analysis.evaluate_circuit_depth(system_code)
		return depth


	elif performance_criterion == "time":
		circuit_time = collections.defaultdict(int)

		# qchip data 로 부터 필요한 정보만 추출함
		physical_cnot_time = qchip_data.get("cnot_gate_time")
		physical_measure_time = qchip_data.get("measure_time")
		gate_time = qchip_data.get("gate_time")

		for inst in system_code:
			# swap first
			if inst[0] in [swap]:
				ctrl, trgt = map(int, inst[1:])
				conditionA = ctrl in physical_cnot_time.keys() and trgt in physical_cnot_time.keys()
				conditionB = trgt in physical_cnot_time[ctrl] and ctrl in physical_cnot_time[trgt]

				if conditionA and conditionB:
					reference_time = max(circuit_time[ctrl], circuit_time[trgt])
					circuit_time[trgt] = circuit_time[ctrl] = reference_time +\
						(2*physical_cnot_time[ctrl][trgt]+physical_cnot_time[trgt][ctrl])
				else:
					raise Exception("physical swap is not possible for {} and {}".format(ctrl, trgt))

			elif inst[0] in list_2q_gates:
				ctrl, trgt = map(int, inst[1:])

				if ctrl in physical_cnot_time.keys() and trgt in physical_cnot_time[ctrl]:
					reference_time = max(circuit_time[ctrl], circuit_time[trgt])
					circuit_time[trgt] = circuit_time[ctrl] = reference_time + physical_cnot_time[ctrl][trgt]
				else:
					raise Exception("physical cnot connection is not possible for {} and {}".format(ctrl, trgt))


			elif inst[0] in list_measure:
				if physical_measure_time is not None:
					trgt = int(inst[1])
					circuit_time[trgt] += physical_measure_time[trgt]

			elif inst[0] in list_barrier: 
				continue

			else:
				if flag_single_qubit and gate_time is not None:

					if inst[0] in [u]:
						_, _, _, trgt = inst[1:]

					elif inst[0] in list_1q_rotations:
						_, trgt = inst[1:]
						
					else:
						trgt = inst[1]

					trgt = int(trgt)
					circuit_time[trgt] += gate_time[trgt]

		return max(list(circuit_time.values()))


	# circuit fidelity 측면에서 시스템 코드 평가
	elif performance_criterion == "fidelity":
		physical_cnot_error = qchip_data.get("cnot_error_rate")
		physical_measure_error = qchip_data.get("measure_error")
		error_rate = qchip_data.get("error_rate")

		circuit_fidelity = {k: 1 for k in inverse_qubit_mapping.keys()}
		
		fidelity_calculation = kwargs.get("fidelity_calculation")
		if fidelity_calculation is None: fidelity_calculation = "product"

		if fidelity_calculation == "propagation":
			time_ordered_syscode = formatconversion.transform_ordered_syscode(system_code, 
				qchip_data=qchip_data)
			list_qubits = formatconversion.extract_list_qubits(time_ordered_syscode)

			for time_idx, list_instructions in time_ordered_syscode.items():
				flag_qubits = {k: False for k in list_qubits}
				for instruction in list_instructions:
					token = p.findall(instruction)
					if not len(token): continue

					# swap first
					if token[0] in [swap]:
						ctrl, trgt = map(int, token[1:])
						conditionA = ctrl in physical_cnot_error.keys() and trgt in physical_cnot_error.keys()
						conditionB = trgt in physical_cnot_error[ctrl] and ctrl in physical_cnot_error[trgt]

						if conditionA and conditionB:
							circuit_fidelity[ctrl] *= (1-physical_cnot_error[ctrl][trgt])**2 * (1-physical_cnot_error[trgt][ctrl]) * circuit_fidelity[trgt]
						else:
							raise Exception("physical swap is not possible for {} and {}".format(ctrl, trgt))
						
						flag_qubits[ctrl] = flag_qubits[trgt] = True
						circuit_fidelity[trgt] = circuit_fidelity[ctrl]

					elif token[0] in list_2q_gates:
						*angle, ctrl, trgt = token[1:]
						ctrl, trgt = map(int, [ctrl, trgt])

						if ctrl in physical_cnot_error.keys() and trgt in physical_cnot_error[ctrl]:
							circuit_fidelity[ctrl] *= (1-physical_cnot_error[ctrl][trgt]) * circuit_fidelity[trgt]
						else:
							raise Exception("physical cnot connection is not possible for {} and {}".format(ctrl, trgt))

						flag_qubits[ctrl] = flag_qubits[trgt] = True
						circuit_fidelity[trgt] = circuit_fidelity[ctrl]

					elif token[0] in list_measure:
						if physical_measure_error is not None:
							trgt = int(token[1])
							circuit_fidelity[trgt] *= (1-physical_measure_error[trgt])
							flag_qubits[trgt] = True

					elif inst[0] in list_barrier: continue

					else: 
						if flag_single_qubit and error_rate is not None:
							if token[0] in [u]:
								_, _, _, trgt = token[1:]

							if token[0] in list_1q_rotations:
								_, trgt = token[1:]
								
							else:
								trgt = token[1]

							trgt = int(trgt)
							circuit_fidelity[trgt] *= (1-error_rate[trgt])
							flag_qubits[trgt] = True


		elif fidelity_calculation == "product":
			for inst in system_code:
				if inst[0] in [swap]:
					ctrl, trgt = map(int, inst[1:])
					conditionA = ctrl in physical_cnot_error.keys() and trgt in physical_cnot_error.keys()
					conditionB = trgt in physical_cnot_error[ctrl] and ctrl in physical_cnot_error[trgt]

					if conditionA and conditionB:
						circuit_fidelity[ctrl] *= (1-physical_cnot_error[ctrl][trgt])**2 * (1-physical_cnot_error[trgt][ctrl])
						circuit_fidelity[trgt] *= (1-physical_cnot_error[ctrl][trgt])**2 * (1-physical_cnot_error[trgt][ctrl])

					else:
						raise Exception("physical swap is not possible for {} and {}".format(ctrl, trgt))

				elif inst[0] in list_2q_gates:
					*angle, ctrl, trgt = inst[1:]
					ctrl, trgt = map(int, [ctrl, trgt])

					if ctrl in physical_cnot_error.keys() and trgt in physical_cnot_error[ctrl]:
						circuit_fidelity[ctrl] *= (1-physical_cnot_error[ctrl][trgt])
						circuit_fidelity[trgt] *= (1-physical_cnot_error[ctrl][trgt])
					else:
						raise Exception("physical CNOT connection is not possible for {} and {}".format(ctrl, trgt))


				elif inst[0] in list_measure:
					if physical_measure_error is not None:
						trgt = int(inst[1])
						circuit_fidelity[trgt] *= physical_measure_error[trgt]

				
				elif inst[0] in list_barrier: continue

				else:
					if flag_single_qubit and error_rate is not None:
						if inst[0] in [u]:
							_, _, _, trgt = inst[1:]

						elif inst[0] in list_1q_rotations:
							_, trgt = inst[1:]

						else:
							trgt = inst[1]

						trgt = int(trgt)
						circuit_fidelity[trgt] *= error_rate[trgt]

		return min(list(circuit_fidelity.values()))



def manage_forward_traversal(args, conn):
	'''
		first forward traversal (with random initial mapping) 관리 함수	
	'''
	qchip_size = len(args["QChip"]["qubit_connectivity"])
	
	qubit_mapping = qubitmapping.initialize_qubit_mapping(args["algorithm_qubits"], qchip_size, 
																option="random", seed=args.get("seed"))

	calibration_type = args.get("calibration_type")

	# calibration_type : time
	if calibration_type == "time":
		while True:
			ret = SABRE(args["DAG"], args["FL"], qubit_mapping, args["DM"], args["QChip"], 
					cost=args["cost"], 
					calibration_type=calibration_type, 
					write_syscode=True, 
					TM=args["TM"], 
					allow_swap=args["swap"],
					direction="forward")
			if ret: break

	# calibration_type : fidelity
	elif calibration_type == "fidelity":
		while True:
			# qubit_mapping = qubitmapping.initialize_qubit_mapping(args["algorithm_qubits"], qchip_size, option="random")
			ret = SABRE(args["DAG"], args["FL"], qubit_mapping, args["DM"], args["QChip"], 
					cost=args["cost"], 
					calibration_type=calibration_type, 
					write_syscode=True, 
					EM=args["EM"], 
					allow_swap=args["swap"],
					direction="forward")
			if ret: break			

	# non-calibration
	else:
		while True:
			# qubit_mapping = qubitmapping.initialize_qubit_mapping(args["algorithm_qubits"], qchip_size, option="random")
			ret = SABRE(args["DAG"], args["FL"], qubit_mapping, args["DM"], args["QChip"], 
					cost=args["cost"], 
					calibration_type=calibration_type, 
					write_syscode=True, 
					allow_swap=args["swap"], 
					direction="forward")
			if ret: break

	conn.send(qubit_mapping)


def manage_synthesize(path_QASM, path_qchip, **kwargs):
	# circuit synthesis 총괄하는 함수
	# inputs 준비 --> DAG, Distance Matrix 생성
	
	synthesis_option = kwargs.get("synthesis_option")
	if synthesis_option is None: return

	cost_function = synthesis_option.get("cost")
	if cost_function is None:  cost_function = "nnc"
	if cost_function == "lap":

		global lap_depth
		lap_depth = synthesis_option.get("lap_depth")
		if lap_depth is None: lap_depth = 5
		else: lap_depth = int(lap_depth)

	# iteration : default iteration -> 10
	try:
		iteration = int(synthesis_option.get("iteration"))
	except:
		iteration = 10

	flag_single_qubit = synthesis_option.get("single_qubit")
	if flag_single_qubit is None: flag_single_qubit = False

	flag_allow_swap = synthesis_option.get("allow_swap")
	if flag_allow_swap is None: flag_allow_swap = True

	if isinstance(path_qchip, str):
		json_qchip_data = open(path_qchip).read()
		qchip_data = json.loads(json_qchip_data)

	elif isinstance(path_qchip, dict):
		qchip_data = path_qchip
	
	else:
		raise Exception("path qchip information : {}".format(path_qchip))	

	# quantum chip data 
	qchip_data["qubit_connectivity"] = {int(k): v for k, v in qchip_data["qubit_connectivity"].items()}
	
	# gate performance가 주어지지 않으면, 모든 데이터를 1로 강제로..
	for item in ["cnot_gate_time", "cnot_error_rate"]:
		temp_data = qchip_data.get(item)
		if temp_data:
			qchip_data[item] = {int(k): {int(m): float(n) for m, n in v.items()} for k, v in temp_data.items()}
		# else:
		# 	qchip_data[item] = {k: {v: 1 for v in v_list} for k, v_list in qchip_data["qubit_connectivity"].items()}

	for item in ["measure_time", "measure_error", "error_rate", "gate_time"]:
		temp_data = qchip_data.get(item)
		if temp_data:
			qchip_data[item] = {int(k): float(v) for k, v in qchip_data[item].items()}

	flag_measurement_qubits_specified = False
	list_measure_qubits = qchip_data.get("measurable_qubits")
	if list_measure_qubits is not None and len(list_measure_qubits):
		flag_measurement_qubits_specified = True
		list_measure_qubits = [int(v) for v in list_measure_qubits]

	qchip_size = len(qchip_data["qubit_connectivity"].keys())

	# cost matrix over all qubits
	# cost : distance, CNOT time, CNOT fidelity

	reference_matrix = {}
	for criterion in ["distance", "time", "fidelity"]:
		matrix, path = DM.generateDM(qchip_data, criterion)
		reference_matrix[criterion] = {"matrix": matrix, "path": path}

	# print(" -------------------------")
	# print(" cost matrix ")
	# print(" 1. distance :")
	# print(pandas.DataFrame(reference_matrix["distance"]["matrix"]).to_string())
	# print(" 2. time :")
	# print(pandas.DataFrame(reference_matrix["time"]["matrix"]).to_string())
	# print(" 3. fidelity :")
	# print(pandas.DataFrame(reference_matrix["fidelity"]["matrix"]).to_string())
	
	# pre-analyze a qasm code
	list_qasm_commands, list_algorithm_qubits, list_algorithm_cbits, algorithm_cnot = parse_qasm.analyze_qasm(path_QASM)	
	
	# 측정 큐빗이 명시된 경우라면.. 측정 명령 앞에 측정 가능 큐빗으로 이동하라는 qasm 삽입
	if flag_measurement_qubits_specified:
		modified_list_qasm_commands = []
		for inst in list_qasm_commands:
			if inst[0] in list_measure:
				modified_list_qasm_commands.append([move, inst[1], "measurement_qubit"])
			modified_list_qasm_commands.append(inst)

		list_qasm_commands = modified_list_qasm_commands

	# directed acyclic graph from qasm
	retDAG = DirectedAcyclicGraph.createDAG(list_qasm_commands, algorithm_qubits=list_algorithm_qubits)

	# in backward graph traversal, move instructions don't need to be considered
	# therefore, we delete them for making backward graph
	for inst in list_qasm_commands:
		if inst[0] == move:
			list_qasm_commands.remove(inst)

	reverseDAG = DirectedAcyclicGraph.createDAG(reversed(list_qasm_commands), algorithm_qubits=list_algorithm_qubits)

	calibration_type = synthesis_option.get("calibration_type")
	if calibration_type is None: calibration_type = "fidelity"

	# calibration aware 합성에 사용할 데이터
	# CNOT gate time 과 CNOT gate error rate 정보중 선택
	arguments = {"QChip": qchip_data, 
				"algorithm_qubits": list_algorithm_qubits,
				"calibration_type": calibration_type,
	 		 	"DM": reference_matrix["distance"], 
	 		 	"TM": reference_matrix["time"], 
	 		 	"EM": reference_matrix["fidelity"], 
	 		 	"DAG": retDAG["DAG"], "cost": cost_function,  "swap": flag_allow_swap}
	
	try:
		arguments.update({"seed": int(synthesis_option.get("random_seed"))})
	except: 
		pass
	

	# calibration data type 을 time 으로 선정한 뒤,
	# target fidelity 를 설정해주면, target fidelity 를 만족하는 최소 시간 회로를 생성하도록 함
	try:
		target_fidelity = float(kwargs.get("target_fidelity"))
	except:
		flag_target_fidelity = False
	else:
		flag_target_fidelity = True


	best_performance_time = math.inf
	best_performance_fidelity = -1 * math.inf
	best_mapping = None
	best_syscode = None
	
	flag_must = False

	bar = Bar('Progress', max=iteration)

	collections_mapping_time = list()
	parent_conn, child_conn = multiprocessing.Pipe(duplex=False)

	# SABRE search 를 복수번 반복해서 최상의 circuit 확인 with different random initial mapping
	while True:
		for i in range(iteration):

			start_time = time.process_time()

			# front layer 는 매 SABRE search 에서 수정됨
			# iteration 마다 동일한 입력을 위해서 사전에 복제해야 함
			FL = copy.deepcopy(retDAG["roots"])
			arguments.update({"FL": FL})

			ps = multiprocessing.Process(target=manage_forward_traversal, args=(arguments, child_conn))
			ps.start()

			if not flag_must:
				ps.join(10+float(algorithm_cnot))
			else:
				ps.join()

			if ps.is_alive():
				print("시간 제한 초과로 다음 차수로 넘어감!")
				ps.join()
				ps.terminate()

			else:
				qubit_mapping = parent_conn.recv()

				# backward traversal
				while True:
					copied_reverseDAG = copy.deepcopy(reverseDAG)
					syscode = SABRE(copied_reverseDAG["DAG"], copied_reverseDAG["roots"], qubit_mapping, reference_matrix["distance"], qchip_data, 
						cost=cost_function, 
						calibration_type=calibration_type, 
						write_syscode=True, 
						TM=reference_matrix["time"], 
						EM=reference_matrix["fidelity"], 
						swap=flag_allow_swap,
						direction="backward")

					if syscode : break
					elif type(syscode) == list and len(syscode): break
				
				final_mapping = copy.deepcopy(qubit_mapping)
				
				# final forward traverse circuit
				FL = copy.deepcopy(retDAG["roots"])

				while True:
					copied_DAG = copy.deepcopy(retDAG)
					list_syscode_commands = SABRE(copied_DAG["DAG"], copied_DAG["roots"], qubit_mapping, reference_matrix["distance"], qchip_data, 
												cost=cost_function, 
												write_syscode=True, 
												calibration_type=calibration_type, 
												TM=reference_matrix["time"], 
												EM=reference_matrix["fidelity"], 
												swap=flag_allow_swap,
												direction="forward")

					if list_syscode_commands: break
					elif type(list_syscode_commands) == list and len(list_syscode_commands): break
				
				last_mapping = qubit_mapping

				# list_syscode_commands = formatconversion.cancel_redundancy(list_syscode_commands)
				time_ordered_syscode = formatconversion.transform_ordered_syscode(list_syscode_commands, qchip_data=qchip_data)

				# evaluate the synthesized circuit 
				ret = evaluate_syscode(list_syscode_commands,
									criterion=calibration_type,
									qchip_data=qchip_data,
									qubit_mapping=final_mapping,
									single_qubit=flag_single_qubit,
									fidelity_calculation="product")

				# calibration type 이 time 이면, best performance 를 평가하는 데 time 이 1 순위, fidelity 가 2 순위
				#                    fidelity 이면,                          fidelity 가 1순위, time 이 2 순위
				if calibration_type == "time":
					# 최적 회로와 최악 회로 확인 과정
					if ret <= best_performance_time:
						best_performance_time = ret
						best_mapping = final_mapping
						best_last_mapping = last_mapping
						best_syscode = copy.deepcopy(list_syscode_commands)

				elif calibration_type == "fidelity":
					# 사용자가 지정한 목적 fidelity가 있다면, 그리고 현재 합성 결과의 fidelity 가 목적 fidelity 보다 크다면 시간을 측정함
					# 해당 시간이 기존의 best perfomance time 보다 작으면, 본 합성 결과가 최적이 되는 것임
					
					if ret >= best_performance_fidelity:
						best_performance_fidelity = ret
						best_mapping = final_mapping
						best_last_mapping = last_mapping
						best_syscode = copy.deepcopy(list_syscode_commands)

			end_time = time.process_time()
			collections_mapping_time.append(end_time-start_time)

			bar.next()

		if best_mapping is None: 
			flag_must = True
			iteration = 1

		else:
			break


	# checkup the mapping result is compatible with the given qubit connectivity
	if checkup.checkup_system_code(best_syscode, final_mapping, qchip_data):
		checkup_msg = "the mapping result is compatible with the given qubit connectivity."
	else:
		checkup_msg = "the mapping result is NOT compatible with the given qubit connectivity."
		raise Exception("the mapping result is NOT compatible with the given qubit connectivity.")

	# 시스템 코드내 quantum gate 갯수 count 
	function_list = collections.defaultdict(int)
	for inst in best_syscode: function_list[inst[0]]+=1

	# the part to count the number of swap gates
	cnot_analysis = {}
	cnot_analysis["Algorithm"] = algorithm_cnot
	cnot_analysis["Circuit"] = function_list[cnot] + 3*function_list[swap]
	cnot_analysis["Overhead"] = cnot_analysis["Circuit"] - cnot_analysis["Algorithm"]

	# in the "system_code" mode,
	# the generated system code is written into a file 
	
	best_mapping = {k: v for k, v in best_mapping.items() if "dummy" not in k}
	best_last_mapping  = {k: v for k, v in best_last_mapping.items() if "dummy" not in k}

	best_circuit = formatconversion.transform_ordered_syscode(best_syscode, 
															qchip_data=qchip_data)
	
	time_ordered_circuit = {"circuit": best_circuit,
							"initial_mapping": best_mapping,
							"final_mapping" : best_last_mapping,
							"qubit": list_algorithm_qubits,
							"cbit": list_algorithm_cbits}

	ret = {"system_code": time_ordered_circuit,
		   "qchip": qchip_data,
		   "analysis": {"Function List": function_list,
		   				"CNOT Overhead": cnot_analysis,
		   				"performance" : {"fidelity": str(best_performance_fidelity), 
		   								"time": str(best_performance_time)}},
		   	"checkup": checkup_msg,
		   	"mapping_time": collections_mapping_time}

	# synthesis mode 에 따라 생산된 결과 ret 를 리턴

	return ret


if __name__ == "__main__":
	# multiprocessing.freeze_support()
	##############
	path_qchip = os.path.join("../DB-QChip", "ibmq_16_melbourne-2021-6-14-17.json")
	
	# import layout_generator
	# qchip = layout_generator.generate_regular_qchip_architecture("../DB-QChip", {"height": 4, "width": 5, "length": 1}, 
	# 					architecture=2)

	list_kisti_algorithms = [
						"Bernstein-Vazirani_5q.qasm", 
						# "Bernstein-Vazirani_5q_2.qasm",
						# "CHSH1.qasm", "CHSH2.qasm", "CHSH3.qasm", 
						# "CHSH4.qasm"
						]
	list_algorithms = [
						"1_bell-5q.qasmf", 
						# "2_qft-5q.qasmf",
						# "3_17_14.qasmf", 
						# "3_17_15.qasmf", 
						# "3_17_16.qasmf", 
						# "7_cat.qasmf",
						# "bernstein_vazirani.qasmf", 
						# "grover_01.qasmf", 
						# "hhl.qasmf", 
						]						

	synthesis_option={"random_seed": 0, "cost": "lap", "lap_depth": 5, 
					 "iteration": 10, "calibration_type" : "fidelity",
					  "initialize_qubit_mapping" : "random",
					  "allow_swap": False}

	for algorithm in list_kisti_algorithms:
		path_qasm = os.path.join("../DB-Assembly/KISTI_sample", algorithm)
		# path_qasm = os.path.join("../DB-Assembly", algorithm)
		
		transformed_code, qubit_assoc, cbit_assoc = formatconversion.transform_to_standardqasm(path_qasm)
			
		ret = manage_synthesize(transformed_code, path_qchip, synthesis_option=synthesis_option)
		pprint(ret)
		ret = formatconversion.transform_to_openqasm(ret, qubit_association=qubit_assoc,
															cbit_association=cbit_assoc)

		pprint(ret)
		file_basename = os.path.basename(path_qasm)
		algorithm_name = os.path.splitext(file_basename)[0]
		file_syscode = "{}.syscode".format(algorithm_name)

		with open(file_syscode, "w") as outfile:
			json.dump(ret, outfile, sort_keys=True, indent=4, separators=(',', ':'))


