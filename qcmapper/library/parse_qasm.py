
import os
import re
import collections
from math import *

from library.gatelist import *

get_bigger = lambda a, b: a if a > b else b
get_smaller = lambda a, b: a if a < b else b

parser = re.compile(r"[\{\[\]a-zA-Z0-9_.*/\->\+}]+")


def parse_qasm(list_qasm_commands):
	"""
		parsing qasm and extract the algorithm qubits
	"""

	list_algorithm_qubits = set([])
	list_algorithm_cbits = set([])
	cnot_counts = 0

	for inst in list_qasm_commands:
		if not len(inst): continue
		
		# measure 먼저
		if inst[0] in list_measure:
			# syntax : measure q -> c
			list_algorithm_qubits.add(inst[1])

		elif inst[0] in list_1q_gates: 
			# syntax : U (angle1, angle2, angle3) qubit
			# syntax : p (angle) qubit
			# syntax : H qubit

			# *angle, qubit = inst[1:]
			list_algorithm_qubits.add(inst[-1])

		# 2q gate 구분 : angle 이 주어지냐 안주어지냐.. 
		elif inst[0] in list_2q_gates: 
			*angle, ctrl, trgt = inst[1:]

			# syntax : cnot (cz, swap, .) q1, q2
			# syntax : cphase (rzz) angle q1, q2
			list_algorithm_qubits.update([ctrl, trgt])

			if inst[0] != swap: 
				cnot_counts+=1
			else: 
				cnot_counts+=3

		elif inst[0] in list_register:
			if inst[0] == "Qubit":
				# syntax : Qubit q
				list_algorithm_qubits.add(inst[1])
			else:
				# syntax : Cbit c
				list_algorithm_cbits.add(inst[1])

		elif inst[0] in list_barrier:
			# barrier 에 걸린 큐빗 이라면, 
			# 이미 앞서 선언하고 사용한 큐빗일 것이므로.. algorithm_qubit 에 추가할 필요 없음
			continue
			
		else:
			raise Exception("error happened : {}".format(inst))

	return list_qasm_commands, list(list_algorithm_qubits), list(list_algorithm_cbits), cnot_counts


def analyze_qasm(path_QASM):
	"""
		function to extract list_qasm_commands, list_algorithm_qubits from QASM
	"""
	list_qasm_commands = []
	list_algorithm_qubits = set([])
	
	cnot_counts = 0

	if isinstance(path_QASM, str):
		if os.path.isfile(path_QASM):
			with open(path_QASM, "r") as infile:
				for line in infile:
					tokens = parser.findall(line)
					if not len(tokens): continue
					list_qasm_commands.append(tokens)
		else:
			raise Exception("Given string is not valid path : {}".format(path_QASM))

	else:
		list_qasm_commands = path_QASM

	return parse_qasm(list_qasm_commands)
	


# def evaluate_syscode(system_code, **kwargs):
# 	"""

# 	"""
# 	performance_criterion = kwargs.get("criterion")
# 	qchip_data = kwargs.get("qchip_data")
# 	if qchip_data is None: 
# 		raise error.Error("qchip data is not provided.")
# 	flag_single_qubit = kwargs.get("single_qubit")

# 	if flag_single_qubit is None: flag_single_qubit = False
# 	flag_measurement = kwargs.get("measurement")
# 	if flag_measurement is None: flag_measurement = False

# 	# performance criterion 1
# 	# gates: the number of quantum gates (usually swap)
# 	if performance_criterion == "cnot":
# 		count_cnot = 0
# 		for inst in system_code:
# 			if inst[0] == g.str_gate_swap: count_cnot+=3
# 			elif inst[0] in [g.str_gate_cnot, g.str_gate_cz]: count_cnot+=1

# 		return count_cnot

# 	# performance criterion 2
# 	# circuit depth	
# 	elif performance_criterion == "depth":
# 		qubit_depth = collections.defaultdict(int)
		
# 		for inst in system_code:
# 			if inst[0] in [g.str_gate_swap, g.str_gate_cnot, g.str_gate_cz]:
# 				ctrl, trgt = inst[1:3]
# 				last_step = get_bigger(qubit_depth[ctrl], qubit_depth[trgt])
# 				qubit_depth[ctrl] = qubit_depth[trgt] = last_step + 1

# 			elif inst[0] in [g.str_gate_rx, g.str_gate_rz, g.str_gate_ry]:
# 				qubit = inst[2]
# 				qubit_depth[qubit]+=1

# 			else:
# 				qubit = inst[1]
# 				qubit_depth[qubit]+=1

# 		return max(list(qubit_depth.values()))

# 	# performance criterion 3
# 	# circuit execution time
# 	elif performance_criterion == "time":
# 		qubit_time = collections.defaultdict(float)
# 		# measurement 시간 반영
# 		if flag_measurement and qchip_data.get("measure_time") is not None:
# 			for inst in system_code:
# 				if inst[0] in [g.str_gate_cnot, g.str_gate_cz]:
# 					ctrl, trgt = inst[1:3]
# 					last_time = get_bigger(qubit_time[ctrl], qubit_time[trgt])
# 					gate_time = qchip_data["net_cnot_time"][ctrl][trgt]
# 					qubit_time[ctrl] = qubit_time[trgt] = last_time + gate_time

# 				elif inst[0] == g.str_gate_swap:
# 					ctrl, trgt = inst[1:3]
# 					last_time = get_bigger(qubit_time[ctrl], qubit_time[trgt])
# 					gate_time = qchip_data["net_cnot_time"][ctrl][trgt] * 3
# 					qubit_time[ctrl] = qubit_time[trgt] = last_time + gate_time

# 				elif inst[0] == g.str_gate_measz:
# 					qubit = inst[1]
# 					gate_time = qchip_data["measure_time"][qubit]
# 					qubit_time[qubit] += gate_time
		
# 		# measurement 시간 미반영
# 		else:
# 			for inst in system_code:
# 				if inst[0] in [g.str_gate_cnot, g.str_gate_cz]:
# 					ctrl, trgt = inst[1:3]
# 					last_time = get_bigger(qubit_time[ctrl], qubit_time[trgt])
# 					gate_time = qchip_data["net_cnot_time"][ctrl][trgt]
# 					qubit_time[ctrl] = qubit_time[trgt] = last_time + gate_time

# 				elif inst[0] == g.str_gate_swap:
# 					ctrl, trgt = inst[1:3]
# 					last_time = get_bigger(qubit_time[ctrl], qubit_time[trgt])
# 					gate_time = qchip_data["net_cnot_time"][ctrl][trgt] * 3
# 					qubit_time[ctrl] = qubit_time[trgt] = last_time + gate_time

# 		return max(list(qubit_time.values()))
		
# 	elif performance_criterion == "fidelity":
# 		qubit_fidelity = {qubit: 1 for qubit in qchip_data["qubit_connectivity"].keys()}
# 		if flag_measurement and qchip_data.get("measure_error") is not None:
# 			for inst in system_code:
# 				if inst[0] in [g.str_gate_cnot, g.str_gate_cz]:
# 					ctrl, trgt = inst[1:3]
# 					last_fidelity = get_smaller(qubit_fidelity[ctrl], qubit_fidelity[trgt])
# 					gate_fidelity = qchip_data["net_cnot_error"][ctrl][trgt]
# 					qubit_fidelity[ctrl] = qubit_fidelity[trgt] = last_fidelity * gate_fidelity

# 				elif inst[0] in [g.str_gate_swap]:
# 					ctrl, trgt = inst[1:3]
# 					last_fidelity = get_smaller(qubit_fidelity[ctrl], qubit_fidelity[trgt])
# 					gate_fidelity = qchip_data["net_cnot_error"][ctrl][trgt]**3
# 					qubit_fidelity[ctrl] = qubit_fidelity[trgt] = last_fidelity * gate_fidelity

# 				elif inst[0] in [g.str_gate_measz]:
# 					qubit = inst[1]
# 					gate_fidelity = 1 - qchip_data["measure_error"]
# 					qubit_fidelity[qubit] *= gate_fidelity

# 		else:
# 			for inst in system_code:
# 				if inst[0] in [g.str_gate_cnot, g.str_gate_cz]:
# 					ctrl, trgt = inst[1:3]
# 					last_fidelity = get_smaller(qubit_fidelity[ctrl], qubit_fidelity[trgt])
# 					gate_fidelity = qchip_data["net_cnot_error"][ctrl][trgt]
# 					qubit_fidelity[ctrl] = qubit_fidelity[trgt] = last_fidelity * gate_fidelity

# 				elif inst[0] in [g.str_gate_swap]:
# 					ctrl, trgt = inst[1:3]
# 					last_fidelity = get_smaller(qubit_fidelity[ctrl], qubit_fidelity[trgt])
# 					gate_fidelity = qchip_data["net_cnot_error"][ctrl][trgt]**3
# 					qubit_fidelity[ctrl] = qubit_fidelity[trgt] = last_fidelity * gate_fidelity

# 		return min(list(qubit_fidelity.values()))
# 					