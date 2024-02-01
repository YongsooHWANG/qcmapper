# -*-coding:utf-8-*-

"""
	working version
"""
import os
import re

import collections
import math
import copy
import time

import multiprocessing
from progress.bar import Bar
from pprint import pprint
import simplejson as json

import library.parse_qasm as parse_qasm
import library.DirectedAcyclicGraph as DirectedAcyclicGraph
import library.DistanceMatrix as DM
import qubitmapping.qubitmapping as qubitmapping

import library.formatconversion as formatconversion
import library.checkup as checkup
from library.gatelist import *
import library.depth_analysis as depth_analysis


decay = 0.1
decay_increment = 1 + decay
extended_set_weight = 0.5

lap_depth = 0
get_smaller = lambda a, b: a if a < b else b

list_measured_qubit = []


def calculate_NNC_cost(listSWAP, FL, DM, MT, **kwargs):
	'''
		cost function based on nearest neighbor cost
	'''
	cost = {k: 0 for k in listSWAP}
	for pair in listSWAP:
		MT[pair[0]], MT[pair[1]] = MT[pair[1]], MT[pair[0]]
		for node in FL:
			if node["gate"] in [swap]: 
				# 알고리즘 상에 SWAP 이 있는 경우, classically 큐빗 relabelling 하는 것으로 처리함...
				# 실제 양자 연산을 수행하지는 않음 -> 따라서 비용은 0
				continue
			
			elif node["gate"] in list_2q_gates:
				cost[pair] += DM[MT[node["ctrl"]]][MT[node["trgt"]]]
			
			elif node["gate"] in [move]:
				cost[pair] += DM[MT[node["ctrl"]]][node["trgt"]]
		
		MT[pair[0]], MT[pair[1]] = MT[pair[1]], MT[pair[0]]

		# return sum(DM[MT[node["ctrl"]]][MT[node["trgt"]]] for node in FL if node["gate"] in list_2q_gates)
	return cost


def calculate_LAP_cost(listSWAP, DAG, FL, DM, MT, listDecay, decay_factor):
	'''
		cost function based on Look-Ahead Ability and Parallelism
	'''

	cost = {}
	for pair in listSWAP:
		MT[pair[0]], MT[pair[1]] = MT[pair[1]], MT[pair[0]]

		temp_cost_E = 0
		temp_cost_F = 0

		decay = max(listDecay[pair[0]], listDecay[pair[1]]) + decay_factor
		extended_set = []

		for node in FL:
			if node["gate"] in [swap]: continue
		
			elif node["gate"] in list_2q_gates:
				associated_physical_ctrl_qubit = MT[node["ctrl"]]
				associated_physical_trgt_qubit = MT[node["trgt"]]

			elif node["gate"] == move:
				associated_physical_ctrl_qubit = MT[node["ctrl"]]
				associated_physical_trgt_qubit = node["trgt"]

			elif node["gate"] in list_barrier: continue

			temp_cost_F += DM[associated_physical_ctrl_qubit][associated_physical_trgt_qubit]
			# node 기준으로 인접해 있는 그래프 상의 노드들을 가져오는 기능
			# extended set 을 일괄적으로 다루기 위해서, 일단 먼저 병합
			extended_set.extend(DirectedAcyclicGraph.get_children_from_node(DAG, node, lap_depth))

		for node in extended_set:
			if node["gate"] in [swap]: continue
			elif node["gate"] in list_2q_gates:
				associated_physical_ctrl_qubit = MT[node["ctrl"]]
				associated_physical_trgt_qubit = MT[node["trgt"]]
			elif node["gate"] == move:
				associated_physical_ctrl_qubit = MT[node["ctrl"]]
				associated_physical_trgt_qubit = node["trgt"]

			elif node["gate"] in list_barrier: continue
			
			temp_cost_E += DM[associated_physical_ctrl_qubit][associated_physical_trgt_qubit]

		cost[pair] = float(temp_cost_F/len(FL))
		if len(extended_set):
			cost[pair] += extended_set_weight*float(temp_cost_E/len(extended_set))

		cost[pair] *= decay
		MT[pair[0]], MT[pair[1]] = MT[pair[1]], MT[pair[0]]

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

	cost_function = kwargs.get("cost")
	previous_best_SWAP = None

	physical_connectivity = qchip_data.get("qubit_connectivity")

	flag_write_syscode = kwargs.get("write_syscode")
	if flag_write_syscode is None: flag_write_syscode = False

	SABRE_direction = kwargs.get("direction")
	flag_swap = kwargs.get("allow_swap")
	lap_depth = kwargs.get("lap_depth")

	inverse_MT = {v: k for k, v in MT.items()}
	listDecay = collections.defaultdict(lambda: 0)

	list_executed_nodes = set([])

	list_waiting_operatios = collections.defaultdict(list)

	while len(FL):
		Execute_gate_list = []

		# 1. front layer 에 있는 게이트들 중 직접 수행 가능한 연산들을 Execute_gate_list 에 추가
		for node in FL:
			# 1q gates, register, swap
			if node["gate"] in list_1q_gates + list_register + [swap]:
				# 알고리즘에 포함된 swap 경우 두가지 방법이 있음
				# 실제로 큐빗의 위치를 바꾸는 경우
				# 큐빗 레이블링만 교환하는 경우..
				Execute_gate_list.append(node)

			elif node["gate"] in list_2q_gates:
				ctrl_qubit = node["ctrl"]
				trgt_qubit = node["trgt"]

				if MT[trgt_qubit] in physical_connectivity[MT[ctrl_qubit]]:
					Execute_gate_list.append(node)

			# # one-qubit gate
			# elif node["gate"] in list_1q_gates + list_register:
			# 	Execute_gate_list.append(node)

			elif node["gate"] == barrier_all:
				if all(node["gate"] == barrier_all for node in FL):
					Execute_gate_list.append(node)

			elif node["gate"] == barrier:
				if len(FL) == 1 and FL[0] == node: 
					Execute_gate_list.append(node)

				else:
					# barrier 명령으로 locked 된 큐빗
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

			else:
				raise Exception("Error Happened : {}".format(node))
				
		# 2. Execute_gate_list 가 empty 가 아니면, 
		#	1) front layer 에서 해당 게이트 삭제
		# 	2) DAG 에서 해당 게이트의 후속 게이트 확인
		#	3) 해당 후속 게이트의 큐빗 종속성 확인해서, 바로 실행 가능하면, front layer 에 추가
		#		(1) 큐빗 종속성을 확인하는 방법: 해당 게이트의 선행자가 여전히 front layer에 존재하는지 확인
		#		(2) 만약 선행자가 존재하면, 해당 게이트의 실행은 아직 불가능하고
		#		(3) 선행자가 없으면, 해당 게이트는 바로 실행 가능함
		if len(Execute_gate_list):
			# 현재 Front Layer 에 위치한 노드 (id)
			# list_FL_nodes = set([node["id"] for node in FL])
			
			# 실행 가능 게이트 처
			for gate in Execute_gate_list:
				if flag_write_syscode:
					# one qubit gate 타입에 따라 qasm instruction 의 포맷이 약간씩 다름
					if gate["gate"] in list_1q_gates:
						# measurement
						if gate["gate"] in list_measure:
							cbit = gate.get("cbit")
							if cbit is None:
								list_syscode_commands.append([gate["gate"], MT[gate["trgt"]], MT[gate["trgt"]]])
							else:
								list_syscode_commands.append([gate["gate"], MT[gate["trgt"]], gate["cbit"]])
							
							# 측정된 큐빗은 측정 큐빗 목록에 삽입: 추후 해당 큐빗 사용 금지 위해 사용함
							list_measured_qubit.append(MT[gate["trgt"]])

						# rotational gates:
						elif gate["gate"] in list_1q_rotations:
							list_syscode_commands.append([gate["gate"], gate["angle"], MT[gate["trgt"]]])

						# other H, Pauli, T, Tdag gates
						else:
							list_syscode_commands.append([gate["gate"], MT[gate["trgt"]]])

					elif gate["gate"] in list_register:
						if gate["gate"] == "Qubit": 
							list_syscode_commands.append([gate["gate"], MT[gate["trgt"]]])

						else: 
							list_syscode_commands.append([gate["gate"], gate["trgt"]])

					elif node["gate"] == barrier_all:
						list_syscode_commands.append([node["gate"]])

						list_instructions = list_waiting_operatios.get("all")
						if list_instructions is not None:
							FL.extend(list_instructions)
							del list_waiting_operatios["all"]


					elif gate["gate"] == barrier:
						list_qubits = [MT[qubit] for qubit in node["trgt"]]
						list_syscode_commands.append([node["gate"], list_qubits])
						key = gate["id"]

						list_instructions = list_waiting_operatios.get(key)
						if list_instructions is not None:
							FL.extend(list_instructions)
							del list_waiting_operatios[key]

					# two qubit gates : parameteric gates (given angle) and non-parameteric gates
					elif gate["gate"] in list_2q_gates:
						angle = gate.get("angle")
						if angle is None:
							list_syscode_commands.append([gate["gate"], MT[gate["ctrl"]], MT[gate["trgt"]]])
						
							if gate["gate"] == swap:
								MT[gate["ctrl"]], MT[gate["trgt"]] = MT[gate["trgt"]], MT[gate["ctrl"]]

						else:
							list_syscode_commands.append([gate["gate"], angle, MT[gate["ctrl"]], MT[gate["trgt"]]])

					# two qubit gates : non-parametric gates
					# elif gate["gate"] in list_2q_gates:
					# 	list_syscode_commands.append([gate["gate"], MT[gate["ctrl"]], MT[gate["trgt"]]])
						
					# 	if gate["gate"] == swap:
					# 		MT[gate["ctrl"]], MT[gate["trgt"]] = MT[gate["trgt"]], MT[gate["ctrl"]]

					else:
						raise Exception("Error Happened : {}".format(gate))

				# 실행할 수 있는 system code 명령은 FL 에서 삭제
				FL.remove(gate)
				list_executed_nodes.add(gate["id"])
				
				list_barrier_nodes_in_FL = [temp_node for temp_node in FL if temp_node["gate"] == barrier]
				list_barrier_blocked_qubits = []
				inverse_barrier_nodes = collections.defaultdict(list)

				for barrier_node in list_barrier_nodes_in_FL:
					list_barrier_blocked_qubits.extend(barrier_node["trgt"])
					for qubit in barrier_node["trgt"]:
						inverse_barrier_nodes[qubit].append(barrier_node["id"])

				list_barrier_blocked_qubits = set(list_barrier_blocked_qubits)

				# 삭제할 노드 gate["id"] 의 후속 노드들 중에서...
				for j in DAG.successors(gate["id"]):
					# 조상들과 FL 내 원소들간의 비교
					# 단순히 현재 노드의 부모 노드만을 고려하면 안된다.

					# 후속노드 중 노드 j 의 직계 부모
					ancestors = set(DAG.predecessors(j))
					
					# 1. 노드 j의 직계부모가 모두 실행된 노드라면, j 를 FL 에 추가함
					# if ancestors.issubset(list_executed_nodes):
					# 	FL.append(DAG.nodes[j])

					if ancestors.issubset(list_executed_nodes):
						# 현재 FL 에 barrier_all 이 존재하면, 노드 j를 FL에 추가하지 말고, list_waiting_operation["all"] 에 추가
						if barrier_all in [temp_node["gate"] for temp_node in FL]:
							list_waiting_operatios["all"].append(DAG.nodes[j])

						# 현재 FL 에 barrier가 존재하고, 이 barrier 에 의해서 blocked 큐빗에 node[j] 연산이 동작하면,
						# 해당 연산은 list_waiting_operation 으로...
						
						elif len(list_barrier_nodes_in_FL):
							if DAG.nodes[j]["gate"] in list_2q_gates:
								if DAG.nodes[j]["ctrl"] in list_barrier_blocked_qubits: 
									# 어느 노드의 후속으로 waiting 시킬 건지 확인
									for barrier_node_id in inverse_barrier_nodes[DAG.nodes[j]["ctrl"]]:
										list_waiting_operatios[barrier_node_id].append(DAG.nodes[j])
								
								elif DAG.nodes[j]["trgt"] in list_barrier_blocked_qubits:
									for barrier_node_id in inverse_barrier_nodes[DAG.nodes[j]["trgt"]]:
										list_waiting_operatios[barrier_node_id].append(DAG.nodes[j])

								else:
									FL.append(DAG.nodes[j])	

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
			
			for node in FL:
				if node["gate"] in list_barrier: continue

				if node["gate"] in list_2q_gates:
					ctrl_qubit = node["ctrl"]
					trgt_qubit = node["trgt"]

					physical_index_ctrl = MT[ctrl_qubit]
					physical_index_trgt = MT[trgt_qubit]

				# for control qubit	
				for j in physical_connectivity[physical_index_ctrl]:
					SWAP_candidate_list.append((ctrl_qubit, inverse_MT[j]))

				# for target qubit	
				for j in physical_connectivity[physical_index_trgt]:
					SWAP_candidate_list.append((trgt_qubit, inverse_MT[j]))


			if not len(SWAP_candidate_list): continue
			
			# evaluate SWAP_candidate_list
			if cost_function == "nnc":
				cost = calculate_NNC_cost(SWAP_candidate_list, FL, DM, MT)

			elif cost_function == "lap":
				cost = calculate_LAP_cost(SWAP_candidate_list, DAG, FL, DM, MT, listDecay, 1+decay)
			
			# minimum cost value & apply the chosen swap gate
			best_SWAP = min(cost, key=cost.get)

			if len(cost) > 1:
				# 새로운 optimal SWAP 이 이전 SWAP 과 동일하다면, 
				# 해당 SWAP을 채택하지 않고, 랜덤하게 SWAP 선택 
				while True:
					if best_SWAP != previous_best_SWAP: break
					else:
						del cost[best_SWAP]
						best_SWAP = min(cost, key=cost.get)

			listDecay[best_SWAP[0]] += decay_increment
			listDecay[best_SWAP[1]] += decay_increment

			associated_physical_qubit_ctrl = MT[best_SWAP[0]]
			associated_physical_qubit_trgt = MT[best_SWAP[1]]
			
			# swapping qubit mapping table
			MT[best_SWAP[0]], MT[best_SWAP[1]] = MT[best_SWAP[1]], MT[best_SWAP[0]]

			inverse_MT[associated_physical_qubit_ctrl], inverse_MT[associated_physical_qubit_trgt] =\
				inverse_MT[associated_physical_qubit_trgt], inverse_MT[associated_physical_qubit_ctrl]

			previous_best_SWAP = best_SWAP

			if flag_write_syscode:
				# update qubit mapping by performing the best SWAP
				# to run the system code, it is written with the qubit index not the algorithm qubit name
				# therefore, the algorithm qubit is mapped to the qubit index through the qubit mapping table, MT
				if flag_swap:
					list_syscode_commands.append([swap, MT[best_SWAP[0]], MT[best_SWAP[1]]])
				else:
					# swap a, b -> CNOT a, b / CNOT b, a / CNOT a, b
					list_syscode_commands.append([cnot, MT[best_SWAP[0]], MT[best_SWAP[1]]])
					list_syscode_commands.append([cnot, MT[best_SWAP[1]], MT[best_SWAP[0]]])
					list_syscode_commands.append([cnot, MT[best_SWAP[0]], MT[best_SWAP[1]]])
		
	if flag_write_syscode:
		return list_syscode_commands


def manage_forward_traversal(args, conn):
	'''
		first forward traversal (with random initial mapping) 관리 함수	
	'''
	qchip_size = len(args["QChip"]["qubit_connectivity"])
	seed = args.get("seed")

	initial_mapping = args.get("initial_mapping")
	initial_mapping_option = args.get("initial_mapping_option")

	lap_depth = args.get("lap_depth")

	qubit_mapping = qubitmapping.initialize_qubit_mapping(args["algorithm_qubits"], qchip_size, 
														option=initial_mapping_option, 
														seed=seed, fixed_qubits=initial_mapping)
	initial_mapping = copy.deepcopy(qubit_mapping)

	list_syscode_commands = SABRE(args["DAG"], args["FL"], qubit_mapping, args["DM"], args["QChip"], 
								cost=args["cost"], write_syscode=True, allow_swap=args["allow_swap"], lap_depth=lap_depth,
								direction="forward")

	# initial_mapping : first forward traversal 할 때의 초기 큐빗 배치
	# qubit_mapping : first forward traversal 하고 난 다음의 큐빗 배치 상태
	conn.send((initial_mapping, qubit_mapping, list_syscode_commands))


def manage_synthesize(path_QASM, path_qchip, **kwargs):
	# circuit synthesis 총괄하는 함수
	# inputs 준비 --> DAG, Distance Matrix 생성

	synthesis_option = kwargs.get("synthesis_option")

	cost_function = synthesis_option.get("cost")
	if cost_function is None: cost_function = "nnc"

	iteration = synthesis_option.get("iteration")
	if iteration is None: iteration = 100

	flag_swap = synthesis_option.get("allow_swap")
	if flag_swap is None: flag_swap = False

	flag_post_processing = synthesis_option.get("post")
	if flag_post_processing is None: flag_post_processing = False

	optimal_criterion = synthesis_option.get("optimal_criterion")
	if optimal_criterion is None:
		optimal_criterion = "circuit_depth"

	initial_mapping_option = synthesis_option.get("initial_mapping_option")
	if initial_mapping_option is None:
		initial_mapping_option = "random"

	# information about quantum chip
	if isinstance(path_qchip, str):
		json_qchip_data = open(path_qchip).read()
		qchip_data = json.loads(json_qchip_data)
	elif isinstance(path_qchip, dict):
		qchip_data = path_qchip
	else:
		raise Exception("path qchip : {}".format(path_qchip))

	qchip_data["qubit_connectivity"] = {int(k): v for k, v in qchip_data["qubit_connectivity"].items()}

	qchip_size = len(qchip_data["qubit_connectivity"].keys())
	
	# distance matrix over qubits, the information about the position of qubits is from quantum chip
	# retDM = DM.generateDM(qchip_data["qubit_connectivity"], "distance")
	retDM, _ = DM.generateDM(qchip_data, "distance")

	# file about the given qubit mapping
	# some qubits are placed at certain positions unconditionally
	file_qubit_mapping = kwargs.get("qubit_table")
	if file_qubit_mapping is not None:
		json_fixed_qubit_mapping = open(file_qubit_mapping).read()
		fixed_mapping = json.loads(json_fixed_qubit_mapping)
	else:
		fixed_mapping = None

	# initial values for evaluating the generated quantum circuits
	optimal_performance = math.inf
	best_mapping = None
	best_circuit = None
	max_circuit_depth = 0
	
	# pre-analyze a qasm code
	list_qasm_commands = []
	list_algorithm_qubits = []
	list_algorithm_cbits = []

	# analyze the given QASM
	list_qasm_commands, list_algorithm_qubits, list_algorithm_cbits, cnot_counts = parse_qasm.analyze_qasm(path_QASM)

	# 초기에 큐빗 또는 비트 선언 같은 파트(알고리즘 큐빗 목록 == none) 에서는 iteration 1회만 수행
	# 일반 non-trivial 회로 (양자 게이트 포함) 인 경우, iteration 은 사용자에 의해 주어진 횟수만큼 수

	if not len(list_algorithm_qubits): iteration = 1

	# 양자칩이 all-to-all connection 이면, iteration 1회면 충분함
	# all-to-all connection : 모든 큐빗의 degree 가 양자칩 사이즈 - 1 이고, 자신은 자신에게 연결되지 않음
	conditionA = all(len(v) == qchip_size - 1 for v in qchip_data["qubit_connectivity"].values())
	if conditionA:
		conditionB = all(k not in v for k, v in qchip_data["qubit_connectivity"].items())
		if conditionB: iteration = 1

	# directed acyclic graph from qasm
	retDAG = DirectedAcyclicGraph.createDAG(list_qasm_commands, algorithm_qubits=list_algorithm_qubits)
	reverseDAG = DirectedAcyclicGraph.createDAG(reversed(list_qasm_commands), algorithm_qubits=list_algorithm_qubits)
	
	collections_mapping_time = list()

	# arguments required for executing SABRE algorithm
	arguments = {"QChip": qchip_data, 
				 "algorithm_qubits": list_algorithm_qubits, 
				 "initial_mapping": fixed_mapping,
				 "initial_mapping_option": initial_mapping_option,
	 		 	 "DM": retDM, 
	 		 	 "DAG": retDAG["DAG"], 
	 		 	 "cost": cost_function, 
	 		 	 "allow_swap": flag_swap}

	if cost_function == "lap":
		lap_depth = synthesis_option.get("lap_depth")
		if lap_depth is None: lap_depth = 1
		arguments.update({"lap_depth": lap_depth})
	
	if "random_seed" in synthesis_option:
		arguments.update({"seed": synthesis_option["random_seed"]})

	flag_must = False
	
	bar = Bar('Progress', max=iteration)
	parent_conn, child_conn = multiprocessing.Pipe(duplex=False)

	# SABRE search 를 복수번 반복해서 최상의 circuit 확인 with different random initial mapping
	while True:
		for i in range(iteration):
			# front layer 는 매 SABRE search 에서 수정됨
			# iteration 마다 동일한 입력을 위해서 사전에 복제해야 함
			start_time = time.process_time()

			FL = copy.deepcopy(retDAG["roots"])
			arguments.update({"FL": FL})
			
			ps = multiprocessing.Process(target=manage_forward_traversal, args=(arguments, child_conn))
			ps.start()
			
			if not flag_must:
				ps.join(float(cnot_counts) + len(list_qasm_commands))
			else:
				ps.join()
			
			if ps.is_alive():
				ps.join()
				ps.terminate()

			else:
				# initial_mapping : 초기 큐빗 배치 at first forward traversal
				# qubit_mapping : first forward traversal 수행 후 큐빗 배치 상태
				initial_mapping, qubit_mapping, list_syscode_commands = parent_conn.recv()[:]

				# first forward traversal 수행 후,
				# 사전에 주어진 고정된 큐빗 배치가 결과 큐빗 배치에서 유지되고 있으면, 이후 graph traversal 수행

				if fixed_mapping is not None and not (fixed_mapping.items() <= qubit_mapping.items()):
					final_mapping = qubit_mapping
				
				else:
					# backward traversal
					SABRE(reverseDAG["DAG"], reverseDAG["roots"], qubit_mapping, retDM, qchip_data, 
						cost=cost_function, write_syscode=False, allow_swap=flag_swap,
						direction="backward")
					
					# backward traversal 수행 이후, 큐빗 배치 상태가 최종의 초기 큐빗 배치가 됨
					initial_mapping = copy.deepcopy(qubit_mapping)

					# final forward traverse circuit
					FL = copy.deepcopy(retDAG["roots"])
					list_syscode_commands = SABRE(retDAG["DAG"], FL, qubit_mapping, retDM, qchip_data, 
													cost=cost_function, write_syscode=True, allow_swap=flag_swap,
													direction="forward")
					# last forward traversal 수행 후, 큐빗 배치 상태가 최종 마지막 큐빗 배치 상태임
					final_mapping = copy.deepcopy(qubit_mapping)
				
				list_syscode_commands = formatconversion.cancel_redundancy(list_syscode_commands)
				
				if optimal_criterion == "circuit_depth":
					circuit_depth = depth_analysis.evaluate_circuit_depth(list_syscode_commands)
					if circuit_depth < optimal_performance:
						optimal_performance = circuit_depth
						best_syscode = list_syscode_commands
						best_mapping = copy.deepcopy(initial_mapping)

				elif optimal_criterion == "number_gates":
					number_instructions = len(list_syscode_commands)
					if number_instructions < optimal_performance:
						optimal_performance = number_instructions
						best_syscode = list_syscode_commands
						best_mapping = copy.deepcopy(initial_mapping)

			end_time = time.process_time()
			collections_mapping_time.append(end_time-start_time)

			bar.next()

		if not best_mapping is None: break
		else:
			flag_must = True
			iteration = 1

	# transform the system code in list to the time ordered format with dictionary	
	best_circuit = formatconversion.transform_ordered_syscode(best_syscode, qchip_data=qchip_data)	
	
	print()
	# checkup the mapping result is compatible with the given qubit connectivity
	if checkup.checkup_system_code(best_circuit, final_mapping, qchip_data):
		checkup_msg = "mapping result is compatible with the given qubit connectivity."
	else:
		checkup_msg = "mapping result is NOT compatible with the given qubit connectivity."
		raise Exception("mapping result is NOT compatible with the given qubit connectivity.")
	
	function_list = collections.defaultdict(int)
	for list_inst in list(best_circuit.values()):
		for inst in list_inst:
			token = inst.split(" ")
			function_list[token[0]]+=1
	
	best_mapping = {k: v for k, v in best_mapping.items() if "dummy" not in k}
	# best_mapping = {k: v for k, v in best_mapping.items()}
	final_mapping = {k: v for k, v in final_mapping.items() if "dummy" not in k}
	# final_mapping = {k: v for k, v in final_mapping.items()}
	
	try:
		circuit_depth = max(best_circuit.keys())+1
	except Exception as e: 
		circuit_depth = 0
	
	# in the "system_code" mode, the generated system code is written into a file 
	# if synthesis_option["goal"] == "system_code":
	time_ordered_circuit = {}
	time_ordered_circuit["circuit"] = best_circuit
	time_ordered_circuit["initial_mapping"] = best_mapping
	time_ordered_circuit["final_mapping"] = final_mapping
	time_ordered_circuit["qubit"] = list_algorithm_qubits
	time_ordered_circuit["cbit"] = list_algorithm_cbits
	
	ret = {"system_code": time_ordered_circuit,
		   "qchip": qchip_data,
		   "analysis": {"Function List": function_list,
		   				"Qubit": {"Qubit": len(best_mapping.keys()),
		   						  "Layout Size": qchip_data.get("dimension")},
		   				# "CNOT Overhead": cnot_analysis,
		   				"Circuit Depth": circuit_depth},
		   	"mapping_time": collections_mapping_time}

	# synthesis mode 에 따라 생산된 결과 ret 를 리턴
	return ret


if __name__ == "__main__":
	# multiprocessing.freeze_support()

	path_qchip = os.path.join("../DB-QChip", "ibmq_16_melbourne-2021-6-14-17.json")
	list_kisti_algorithms = [
						"Bernstein-Vazirani_5q.qasm", 
						"Bernstein-Vazirani_5q_2.qasm",
						"CHSH1.qasm", "CHSH2.qasm", "CHSH3.qasm", 
						"CHSH4.qasm"
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

	synthesis_option={"goal": "system_code", "allow_swap": True, "random_seed": 0,
					"iteration": 10, "cost": "lap", "optimal_criterion": "number_gates",
					"initial_mapping_option": "periodic_random"}

	for algorithm in list_kisti_algorithms:
		path_qasm = os.path.join("../DB-Assembly/KISTI_sample", algorithm)
		# path_qasm = os.path.join("../DB-Assembly", algorithm)
		transformed_code, qubit_assoc, cbit_assoc = formatconversion.transform_to_standardqasm(path_qasm)
		
		ret = manage_synthesize(transformed_code, path_qchip, synthesis_option=synthesis_option)
		ret = formatconversion.transform_to_openqasm(ret, qubit_association=qubit_assoc,
															cbit_association=cbit_assoc)
		
		file_basename = os.path.basename(path_qasm)
		algorithm_name = os.path.splitext(file_basename)[0]
		file_syscode = "{}.syscode".format(algorithm_name)

		with open(file_syscode, "w") as outfile:
			json.dump(ret, outfile, sort_keys=True, indent=4, separators=(',', ':'))
