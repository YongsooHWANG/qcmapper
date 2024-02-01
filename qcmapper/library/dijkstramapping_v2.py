# -*-coding:utf-8-*-
import os
import re
import sys
import collections
import copy
import time
import math
from pprint import pprint
import simplejson as json

import library.circuit as circuit
import library.graph as graph

import qubitmapping.qubitmapping as qubitmapping
import library.parse_qasm as parse_qasm

import library.formatconversion as formatconversion
from library.gatelist import *

from progress.bar import Bar


parser = re.compile("[\{a-zA-Z0-9_.*\-+/->\}]+")


def synthesize_dijkstra_manner(qasm_code, path_qchip, **kwargs):
	# over several loops
	# qubit table 생성 -> random
	# synthesize_flattened_qasm 함수 호출... 최소 time... target function

	if isinstance(path_qchip, str):
		json_qchip_data = open(path_qchip).read()
		qchip_data = json.loads(json_qchip_data)
	elif isinstance(path_qchip, dict):
		qchip_data = path_qchip
	else:
		raise Exception("qchip data is wrong : {}".format(path_qchip))

	synthesis_option = kwargs.get("synthesis_option")
	if synthesis_option is None:
		raise Exception("synthesis option is not provided.")

	calibration_type = synthesis_option.get("calibration_type")
	if calibration_type is None: calibration_type = "depth"
	
	##
	qchip_data["qubit_connectivity"] = {int(k): v for k, v in qchip_data["qubit_connectivity"].items()}
	for item in ["cnot_gate_time", "cnot_error_rate"]:
		qchip_data[item] = {int(k): {int(m): float(n) for m, n in v.items()} for k, v in qchip_data[item].items()}

	for item in ["measure_time", "measure_error", "error_rate", "gate_time"]:
		if item in qchip_data.keys():
			qchip_data[item] = {int(k): float(v) for k, v in qchip_data[item].items()}
	
	flag_measurement_qubits_specified = False
	list_measure_qubits = qchip_data.get("measurement_qubits")
	if list_measure_qubits is not None:
		flag_measurement_qubits_specified = True
		list_measure_qubits = [int(v) for v in list_measure_qubits]

	##

	flag_initial_qubit_table = False
	file_qubit_mapping = kwargs.get("qubit_table")
	if file_qubit_mapping is not None:
		flag_initial_qubit_table = True
		json_data = open(file_qubit_mapping).read()
		qubit_mapping = json.loads(json_data)

	else:
		qubit_mapping = None

	collections_mapping_time = list()

	qchip_size = len(qchip_data["qubit_connectivity"])

	list_qasm_commands, list_algorithm_qubits, list_algorithm_cbits, algorithm_cnot =\
		parse_qasm.analyze_qasm(qasm_code)	
	
	if calibration_type in ["time", "depth"]: best_performance = math.inf
	elif calibration_type in ["fidelity"]: best_performance = 0

	best_performance_data = None

	# setting the number of mapping iterations
	try:
		iteration = int(synthesis_option.get("iteration"))
	except:
		iteration = 10

	if flag_measurement_qubits_specified:
		modified_list_qasm_commands = []
		for inst in list_qasm_commands:
			if inst[0] in list_measure:
				modified_list_qasm_commands.append([move, inst[1], "measurement_qubit"])
			modified_list_qasm_commands.append(inst)

		list_qasm_commands = modified_list_qasm_commands


	if calibration_type == "time":
		# make matrix of time
		
		# swap a, b -> (1) cnot a, b * cnot b, a * cnot a, b
		#              (2) cnot b, a * cnot a, b * cnot b, a
		swap_order = collections.defaultdict(dict)
		weight_graph = collections.defaultdict(dict)
		
		qchip_performance = qchip_data.get("cnot_gate_time")
		
		for ctrl, trgt_list in qchip_performance.items():
			for trgt, cnot_time in trgt_list.items():
				if cnot_time*2 + qchip_performance[trgt][ctrl] < qchip_performance[trgt][ctrl]*2 + cnot_time:
					weight_graph[ctrl][trgt] = weight_graph[trgt][ctrl] = {"direction" : "{}>{}".format(ctrl, trgt),
																			"weight": cnot_time*2 + qchip_performance[trgt][ctrl]}
				else:
					weight_graph[ctrl][trgt] = weight_graph[trgt][ctrl] = {"direction": "{}>{}".format(trgt, ctrl),
																			"weight": qchip_performance[trgt][ctrl]*2 + cnot_time}
		
	elif calibration_type == "fidelity":
		# make matrix of fidelity
		swap_order = collections.defaultdict(dict)
		weight_graph = collections.defaultdict(dict)

		qchip_performance = qchip_data.get("cnot_error_rate")

		for ctrl, trgt_list in qchip_performance.items():
			for trgt, cnot_infidelity in trgt_list.items():
				if cnot_infidelity**2 * qchip_performance[trgt][ctrl] < qchip_performance[trgt][ctrl]**2 * cnot_infidelity:
					weight_graph[ctrl][trgt] = weight_graph[trgt][ctrl] = {"direction" : "{}>{}".format(ctrl, trgt),
																			"weight": cnot_infidelity**2 * qchip_performance[trgt][ctrl]}
				else:
					weight_graph[ctrl][trgt] = weight_graph[trgt][ctrl] = {"direction" : "{}>{}".format(trgt, ctrl),
																			"weight": qchip_performance[trgt][ctrl]**2 * cnot_infidelity}

	else:
		weight_graph = qchip_data.get("qubit_connectivity")


	random_seed = synthesis_option.get("random_seed")

	# qubit mapping table 이 주어지면, 한 번만 매핑하면 됨
	# 만약, qubit mapping table 이 주어지지 않으면, 랜덤으로 생성해서 반복적으로 매핑함
	if flag_initial_qubit_table:
		if len(qubit_mapping) < qchip_size:
			
			qubit_mapping = qubitmapping.initialize_qubit_mapping(list_algorithm_qubits, qchip_size, 
																		option=synthesis_option.get("initial_mapping_option"),
																		seed=random_seed,
																		fixed_qubits=qubit_mapping)

		synthesized_result = synthesize_flattened_qasm(list_qasm_commands,
														qchip=qchip_data,
											 	   		qubit_table=qubit_mapping,
											 	   		option=synthesis_option,
											 	   		weight_graph=weight_graph)
		
		if calibration_type == "fidelity":
			if synthesized_result.get("performance") > best_performance:
				best_performance = synthesized_result.get("performance")
				best_performance_data = synthesized_result
		else:
			if synthesized_result.get("performance") < best_performance:
				best_performance = synthesized_result.get("performance")
				best_performance_data = synthesized_result

	else:
		# perform the mapping as much as the number of **iterations**
		bar = Bar('Progress', max=iteration)

		for iter_idx in range(iteration):
			start_time = time.process_time()

			qubit_mapping = qubitmapping.initialize_qubit_mapping(list_algorithm_qubits, qchip_size, 
																		option=synthesis_option.get("initial_mapping_option"),
																		seed=random_seed,
																		fixed_qubits=qubit_mapping)
			
			synthesized_result = synthesize_flattened_qasm(list_qasm_commands,
															qchip=qchip_data,
															weight_graph=weight_graph,
												 	   		qubit_table=qubit_mapping,
												 	   		option=synthesis_option,
												 	   		algorithm_qubits=list_algorithm_qubits)

			end_time = time.process_time()
			collections_mapping_time.append(end_time-start_time)

			if calibration_type == "fidelity":
				if synthesized_result.get("performance") >= best_performance:
					best_performance = synthesized_result.get("performance")
					best_performance_data = synthesized_result
			else:
				if synthesized_result.get("performance") <= best_performance:
					best_performance = synthesized_result.get("performance")
					best_performance_data = synthesized_result

			bar.next()

	# post-operation to reflect the option "commutable_cnot_swap"
	# 현재 오류 발생해서 수정중 --------------------
	# if "commutable_cnot_swap" in synthesis_option and synthesis_option["commutable_cnot_swap"]:
	# 	best_performance_data["system_code"]["circuit"] =\
	# 		formatconversion.preanalyze_qasm(best_performance_data["system_code"]["circuit"])

	# 성능 분석도 update 필요 --> 아래 부분에서 수행함
	# syscode_cnot = best_performance_data["analysis"]["Function List"]["Logical"]["CNOT"] +\
	# 				3*best_performance_data["analysis"]["Function List"]["Logical"]["SWAP"]

	# best_performance_data["analysis"]["CNOT Overhead"] = {"Algorithm": algorithm_cnot,
	# 													  "Circuit": syscode_cnot,
	# 													  "Overhead": syscode_cnot - algorithm_cnot}

	# ic(best_performance_data["system_code"]["circuit"])
	##################################################################
	# 합성 목적에 따라 후속 작업 수행
	# - 시스템 코드: 1. time ordered format 으로 변환
	#             2. analysis data 타입 포맷 변환 (for dumping the data into json file)
	# - 성능 평가:  1. 합성 수행 결과로 성능 분석 수행 (buidlingblock file, kq file 추가 필요)
	##################################################################
	# if synthesis_option["goal"] == "system_code":
		# transform the simple list --> time ordered system code

	time_ordered_system_code = formatconversion.transform_ordered_syscode(best_performance_data["system_code"]["circuit"], qchip_data=qchip_data)
	
	function_list = collections.defaultdict(int)
	for list_instructions in list(time_ordered_system_code.values()):
		for inst in list_instructions:
			token = inst.split(" ")
			function_list[token[0]] += 1

	best_performance_data["system_code"]["circuit"] = time_ordered_system_code
	best_performance_data["analysis"]["Function List"] = function_list

	flag_time_checking = kwargs.get("time_checking")
	if flag_time_checking:
		best_performance_data["Mapping_Time"] = collections_mapping_time
	# synthesized_result = best_performance_data


	# 데이터 타입 변환: json 파일에 dump 하기 위해 serialized 되지 않는 데이터 타입(numpy64.int or numpy64.float) 을 string 으로...
	# for k, v in synthesized_result["analysis"].items():
	# 	if not isinstance(v, dict):
	# 		synthesized_result["analysis"][k] = str("{0:>}".format(v))
	# 	else:
	# 		for a, b in v.items():
	# 			synthesized_result["analysis"][k][a] = str("{0:>}".format(b))


	return best_performance_data


def synthesize_flattened_qasm(qasm_code, **kwargs):

	qchip_data = kwargs.get("qchip")
	qubit_connectivity = qchip_data.get("qubit_connectivity")
	qchip_size = len(qubit_connectivity)

	qubit_mapping_table = kwargs.get("qubit_table")
	if qubit_mapping_table is None:
		raise Exception("Qubit Mapping Table is not provided.")

	cloned_qubit_mapping = copy.deepcopy(qubit_mapping_table)
	synthesis_option = kwargs.get("option")
	if synthesis_option is None:
		raise Exception("Synthesis Option is not provided.")
	
	calibration_type = synthesis_option.get("calibration_type")
	if calibration_type is None: calibration_type = "depth"

	list_measure_qubits = qchip_data.get("measurement_qubits")

	list_algorithm_qubits = kwargs.get("algorithm_qubits")

	weight_graph = kwargs.get("weight_graph")

	if isinstance(qubit_mapping_table, str):
		json_data = open(qubit_mapping_table).read()
		qubit_mapping_table = json.loads(json_data)
	
	_circuit = circuit.QuantumCircuit(layout_size=qchip_size, 
									  qubit_table=qubit_mapping_table,
									  option=synthesis_option,
									  algorithm_qubits=list_algorithm_qubits)

	list_move_distance = collections.defaultdict(int)
	
	for inst in qasm_code:
		gate = inst[0]
		
		# case 1: one-qubit gate
		if gate in ["Qubit"]:
			qubit_idx = int(qubit_mapping_table[inst[1]])
			_circuit.register(gate, qubit_idx)
		
		elif gate in ["Cbit"]:
			_circuit.register(gate, inst[1])

		elif gate in list_1q_gates:
			if gate in [u]:
				theta, phi, lamda, qubit = inst[1:]
				qubit_idx = int(qubit_mapping_table[qubit])
				_circuit.one_qubit_gate(gate, qubit_idx, theta=theta, phi=phi, lamda=lamda)

			elif gate in list_1q_rotations:
				angle, qubit = inst[1:]
				
				qubit_idx = int(qubit_mapping_table[qubit])
				_circuit.one_qubit_gate(gate, qubit_idx, angle=angle)


			elif gate in list_measure:
				if len(inst) == 4:
					qubit, cbit = inst[1::2]
					qubit_idx = int(qubit_mapping_table[qubit])
					_circuit.one_qubit_gate(gate, qubit_idx, cbit=cbit)

				elif len(inst) == 2:
					qubit = inst[1]
					qubit_idx = int(qubit_mapping_table[qubit])
					_circuit.one_qubit_gate(gate, qubit_idx)

			else:
				qubit = inst[1]
				qubit_idx = int(qubit_mapping_table[qubit])
				_circuit.one_qubit_gate(gate, qubit_idx)

		# move 연산: movement of a qubit to a certain point on a quantum chip
		elif gate in list_move + [swap]:
			ctrl, trgt = inst[1:]

			if trgt == "measurement_qubit":
				ctrl_idx = int(qubit_mapping_table[ctrl])

				# 현재 큐빗을 기준으로.. 가장 가까운 측정 가능 큐빗 찾기
				if calibration_type in [None, "depth"]:
					connection_path = None
					shortest_path_length = math.inf

					for m_qubit in list_measure_qubits:
						path = graph.findQubitConnectionPath(weight_graph, ctrl_idx, trgt_idx, weight_type=calibration_type)
						if len(path) < shortest_path_length:
							shortest_path_length = len(path)
							connection_path = path

				else:
					connection_path = None
					shortest_path_length = math.inf

					if calibration_type == "time":
						for m_qubit in list_measure_qubits:
							path = graph.findQubitConnectionPath(weight_graph, ctrl_idx, m_qubit, weight_type=calibration_type)
							# 경로상의 이동시간 + 측정 시간
							time = sum(edge[2] for edge in path) + qchip_data.get("measure_time")[m_qubit]
							
							if time < shortest_path_length:
								shortest_path_length = time
								connection_path = path

					elif calibration_type == "fidelity":
						for m_qubit in list_measure_qubits:
							path = graph.findQubitConnectionPath(weight_graph, ctrl_idx, m_qubit, weight_type=calibration_type)
							fidelity = qchip_data.get("measure_error")[m_qubit]
							for edge in path: fidelity *= edge[2]

							if fidelity < shortest_path_length:
								shortest_path_length = fidelity
								connection_path = path

			else:
				ctrl_idx = int(qubit_mapping_table[ctrl])
				trgt_idx = int(qubit_mapping_table[trgt])
			
				connection_path = graph.findQubitConnectionPath(weight_graph, ctrl_idx, trgt_idx, weight_type=calibration_type)

			if len(connection_path):
				qubit_mapping_table = _circuit.manage_cnot([ctrl_idx, trgt_idx], path=connection_path, 
																				 qubit_table=qubit_mapping_table,
																				 weight_graph=weight_graph,
																				 calibration_type=calibration_type,
																				 final_operation=swap)
		# case 2: two-qubit gate: cnot 과 swap 고려함
		elif gate in list_2q_gates:
			*angle, ctrl, trgt = inst[1:]
			ctrl_idx = int(qubit_mapping_table[ctrl])
			trgt_idx = int(qubit_mapping_table[trgt])

			connection_path = graph.findQubitConnectionPath(weight_graph, ctrl_idx, trgt_idx, weight_type=calibration_type)
				
			# final_cnot: 원거리 이동후 마지막에 cnot 수행 여부 결정
			qubit_mapping_table = _circuit.manage_cnot([ctrl_idx, trgt_idx], path=connection_path, 
																			 qubit_table=qubit_mapping_table,
																			 weight_graph=weight_graph,
																			 calibration_type=calibration_type,
																			 final_operation=gate,
																			 angle=angle)


		elif gate in list_barrier:
			if gate == barrier_all:
				_circuit.insert_barrier_all()
			
			else:
				qubits_index = [qubit_mapping_table[k] for k in inst[1:]]
				_circuit.insert_barrier(qubits_index)
	
	list_syscode_commands = _circuit.get_system_code()["system_code"]["circuit"]
	# list_syscode_commands = formatconversion.cancel_redundancy(list_syscode_commands)

	# performance = circuit_fidelity.calculate_fidelity_product_list_form(, qchip_data)
	performance = 0	
	return {"system_code": _circuit.get_system_code()["system_code"],
		    "performance": performance,
		    "qchip": qchip_data,
		    "analysis": _circuit.get_analysis()}


if __name__ == "__main__":
	path_qchip = os.path.join("../DB-QChip", "ibmq_16_melbourne-2021-6-14-17.json")
	synthesis_option = {'allow_swap': True,
					"parallel_swap": False,
            		"calibration_type" : "fidelity",
					'cnot_depth': False,
					'goal': 'system_code',
					'iteration': 10,
					'qubit_connectivity': 'user',
					"initial_mapping_option": "random",}
					# 't_depth': True}

	list_kisti_algorithms = [
						"Bernstein-Vazirani_5q.qasm", 
						"Bernstein-Vazirani_5q_2.qasm",
						"CHSH1.qasm", "CHSH2.qasm", "CHSH3.qasm", 
						"CHSH4.qasm"
						]
	list_algorithms = [
						"dj_indep_qiskit_5.qasmf", 
						"ae_indep_qiskit_5.qasmf", 
						"grover-noancilla_indep_qiskit_5.qasmf",
						"portfolioqaoa_indep_qiskit_5.qasmf",
					 	"qaoa_indep_qiskit_5.qasmf",
					 	"qpeexact_indep_qiskit_5.qasmf",
					 	"ghz_indep_qiskit_5.qasmf", 
					 	"graphstate_indep_qiskit_5.qasmf", 
						"portfoliovqe_indep_qiskit_5.qasmf", 
						"realamprandom_indep_qiskit_5.qasmf", 
						"su2random_indep_qiskit_5.qasmf", 
						"twolocalrandom_indep_qiskit_5.qasmf",
						"vqe_indep_qiskit_5.qasmf", 
						"wstate_indep_qiskit_5.qasmf"
						]
	for algorithm in list_kisti_algorithms:
		ic(algorithm)
		# path_qasm = os.path.join("../Experiment_paper/DB-Assembly", algorithm)
		path_qasm = os.path.join("../DB-Assembly/KISTI_sample", algorithm)
		
		transformed_code, qubit_assoc, cbit_assoc = formatconversion.transform_to_standardqasm(path_qasm)

		ret = synthesize_dijkstra_manner(transformed_code, path_qchip, synthesis_option=synthesis_option)
		ret = formatconversion.transform_to_openqasm(ret, qubit_association=qubit_assoc,
															cbit_association=cbit_assoc)
		
		pprint(ret)
