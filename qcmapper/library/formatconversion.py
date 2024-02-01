
import collections
import re
import copy
import itertools
from math import *

from pprint import pprint
from library.gatelist import *

get_bigger = lambda a, b: a if a>b else b

parser = re.compile(r"[\{a-zA-Z0-9_.*\-+/->\}]+")


def cancel_redundancy(syscode):
	"""
		function to cancel out the redundant quantum gates in time order

		args:
			syscode in list
	"""

	table = collections.defaultdict(list)

	for idx, inst in enumerate(syscode):
		# 2-qubit gate : typeA (control, target 큐빗 명시가 중요한 게이트)
		gate = inst[0]

		if gate in list_2q_rotations:
			if len(table[inst[2]]) and len(table[inst[3]]):
				last_instA = table[inst[2]][-1]
				last_instB = table[inst[3]][-1]

				conditionA = (last_instA["gate"] == gate) and (last_instA["qubits"] == inst[2:])
				conditionB = (last_instB["gate"] == gate) and (last_instB["qubits"] == inst[2:])
				
				# 동일하면
				if conditionA and conditionB:
					# 이전 연산의 angle 조정
					new_angle = float(eval(last_instA["angle"])) + float(eval(inst[1]))
					table[inst[2]][-1]["angle"] = str(new_angle)
					table[inst[3]][-1]["angle"] = str(new_angle)

				# 다르면
				else:
					table[inst[2]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})
					table[inst[3]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})

			else:
				table[inst[2]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})
				table[inst[3]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})


		elif any(item in gate for item in list_2q_gates):
			if len(table[inst[1]]) and len(table[inst[2]]):
				# 새로운 명령과 이전 명령 비교
				last_instA = table[inst[1]][-1]
				last_instB = table[inst[2]][-1]

				conditionA = (last_instA["gate"] == gate) and (last_instA["qubits"] == inst[1:])
				conditionB = (last_instB["gate"] == gate) and (last_instB["qubits"] == inst[1:])

				# 동일하면
				if conditionA and conditionB:
					table[inst[1]].pop()
					table[inst[2]].pop()

				# 다르면
				else:
					table[inst[1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
					table[inst[2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

			else:
				table[inst[1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
				table[inst[2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

		# 2-qubit gate : typeB (control, target 큐빗 명시가 필요없는 게이트)
		elif swap in gate:
			if len(table[inst[1]]) and len(table[inst[2]]):
				# 새로운 명령과 이전 명령 비교
				last_instA = table[inst[1]][-1]
				last_instB = table[inst[2]][-1]

				conditionA = (last_instA["gate"] == gate) and (set(last_instA["qubits"]) == set(inst[1:]))
				conditionB = (last_instB["gate"] == gate) and (set(last_instB["qubits"]) == set(inst[1:]))

				# 동일하면
				if conditionA and conditionB:
					table[inst[1]].pop()
					table[inst[2]].pop()

				# 다르면
				else:
					table[inst[1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
					table[inst[2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

			else:
				table[inst[1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
				table[inst[2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

		# barrier-All
		elif gate == barrier_all:
			for qubit in table.keys():
				table[qubit].append({"gate": barrier_all, "idx": idx})

		# selective barrier : 
		elif gate == barrier:
			list_qubits = inst[1]
			for qubit in list_qubits:
				table[qubit].append({"gate": barrier, "qubits": list_qubits, "idx": idx})

		elif gate in list_measure:
			# 동일한 큐빗을 측정해서 다른 cbit 에 넣는 경우가 생길 수 있음
			# 따라서 큐빗과 cbit 이 모두 동일해야 동일한 operation 에 해당함
			if len(table[inst[1]]):
				last_inst = table[inst[1]][-1]
				if (last_inst["gate"] == gate) and (last_inst["qubits"] == inst[1]) and\
					last_inst["cbits"] == inst[2]:
					table[inst[1]].pop()

				else:
					table[inst[1]].append({"gate": gate, "qubits": [inst[1]], "cbits": [inst[2]], "idx": idx})

		# 1-qubit gate
		elif gate == u:
			if len(table[inst[2]]):
				last_inst = table[inst[2]][-1]
				if last_inst["gate"] == gate and last_inst["qubits"] == inst[2:]:
					new_angle = {}

					for axis in ["x", "z", "y"]:
						new_angle[axis] = float(eval(str(last_inst["angle"][axis]))) + float(eval(str(inst[1][axis])))
					table[inst[2]][-1]["angle"] = new_angle
				else:
					table[inst[2]].append({"gate": gate, "qubits": [inst[2]], "angle": inst[1], "idx": idx})	
			else:
				table[inst[2]].append({"gate": gate, "qubits": [inst[2]], "angle": inst[1], "idx": idx})


		elif gate in list_1q_rotations:
			if len(table[inst[2]]):
				last_inst = table[inst[2]][-1]

				# Rz 게이트가 연속되면, angle 확인 후, 앞선 게이트의 angle 값을 변경
				if (last_inst["gate"] == gate) and (last_inst["qubits"] == inst[2:]):
					new_angle = float(eval(last_inst["angle"])) + float(eval(inst[1]))
					table[inst[2]][-1]["angle"] = str(new_angle)

				else:
					table[inst[2]].append({"gate": gate, "qubits": [inst[2]], "idx": idx, "angle": inst[1]})
			else:
				table[inst[2]].append({"gate": gate, "qubits": [inst[2]], "idx": idx, "angle": inst[1]})


		elif gate in list_1q_gates:
			# 새로운 명령과 이전 명령 비교
			if len(table[inst[1]]):
				last_inst = table[inst[1]][-1]

				# 동일하면
				if (last_inst["gate"] == gate) and (last_inst["qubits"] == inst[1:]):
					table[inst[1]].pop()

				# 다르면
				else:
					table[inst[1]].append({"gate": gate, "qubits": [inst[1]], "idx": idx})	

			else:
				table[inst[1]].append({"gate": gate, "qubits": [inst[1]], "idx": idx})

		elif gate in list_register: continue


		else:
			raise Exception("Error Happened : {}".format(inst))

	temp_syscode = {}
	for v_list in list(table.values()):
		for v in v_list:
			temp_syscode[v["idx"]] = v
	
	sorted_index = sorted(temp_syscode.keys())

	post_processed_syscode = []
	for k in sorted_index:
		v = temp_syscode[k]
		gate = v["gate"]

		if swap in gate:
			post_processed_syscode.append([gate, v["qubits"][0], v["qubits"][1]])

		elif gate in list_2q_rotations:
			post_processed_syscode.append([gate, v["angle"], v["qubits"][0], v["qubits"][1]])

		elif any(item in gate for item in list_2q_gates):
			post_processed_syscode.append([gate, v["qubits"][0], v["qubits"][1]])

		elif gate in list_measure:
			post_processed_syscode.append([gate, v["qubits"][0], v["cbits"][0]])

		elif gate in list_1q_rotations:
			post_processed_syscode.append([gate, v["angle"], v["qubits"][0]])

		elif gate in list_1q_gates:
			post_processed_syscode.append([gate, v["qubits"][0]])

		# barrier-all
		elif gate == barrier_all:
			post_processed_syscode.append([gate])
			
		# selective barrier
		elif gate == barrier:
			post_processed_syscode.append([gate, v["qubits"]])

		else:
			raise Exception("Error Happened : {}".format(v))

	return post_processed_syscode	



def transform_ordered_syscode(syscode, **kwargs):
	'''
		개별 게이트의 circuit index를 분석하고, 시간순으로 정리된 회로를 생성 리턴하는 함
	'''
	
	time_index = collections.defaultdict(int)
	ordered_syscode = collections.defaultdict(list)
	
	algorithm_qubits = kwargs.get("algorithm_qubits")
	qchip_data = kwargs.get("qchip_data")
	
	# 큐빗 갯수 확인
	if qchip_data is not None:
		qchip_dimension = qchip_data.get("dimension")
		try:
			number_qubits = qchip_dimension["height"] * qchip_dimension["width"] * qchip_dimension["length"]

		except Exception as e:
			number_qubits = len(qchip_data["qubit_connectivity"].keys())

	for inst in syscode:
		flag_barrier = False

		if inst[0] in list_2q_rotations:
			ctrl, trgt = map(int, inst[2:])

			applying_index = max(time_index[ctrl], time_index[trgt])
			time_index[ctrl] = time_index[trgt] = applying_index+1
			list_command = "{}({}) {},{}".format(inst[0], inst[1], ctrl, trgt)


		elif any(item in inst[0] for item in list_2q_gates):
			ctrl, trgt = map(int, inst[1:])

			applying_index = max(time_index[ctrl], time_index[trgt])
			time_index[ctrl] = time_index[trgt] = applying_index+1
			list_command = "{} {},{}".format(inst[0], ctrl, trgt)
		
		elif inst[0] in list_register: continue
			
		else: 
			if inst[0] == u:
				*angle, qubit = inst[1:]
				qubit = int(qubit)
				if len(angle) == 1 and isinstance(angle[0], dict):
					list_command = "{}({},{},{}) {}".format(inst[0], 
						angle[0].get("x"), angle[0].get("y"), angle[0].get("z"), qubit)
				elif len(angle) == 3:
					list_command = "{}({},{},{}) {}".format(inst[0], 
						angle[0], angle[1], angle[2], qubit)
				
				applying_index = time_index[qubit]
				time_index[qubit] += 1

			elif inst[0] in list_1q_rotations:
				angle, qubit = inst[1:]
				qubit = int(qubit)
				list_command = "{}({}) {}".format(inst[0], angle, qubit)
				applying_index = time_index[qubit]
				time_index[qubit] += 1

			elif inst[0] in list_measure:
				qubit, cbit = inst[1:]
				qubit = int(qubit)
				list_command = " ".join([inst[0], str(qubit), "->", str(cbit)])
				
				applying_index = time_index[qubit]
				time_index[qubit] += 1

			elif inst[0] in list_1q_gates:
				qubit = int(inst[1])
				list_command = "{} {}".format(inst[0], qubit)
				applying_index = time_index[qubit]
				time_index[qubit] += 1
			
		
			elif inst[0] == barrier_all:
				flag_barrier = True
				list_command = g.str_barrier_all
				
				if not len(time_index): applying_index = 0
				else:
					applying_index = max(list(time_index.values()))

				if qchip_data is not None:
					for qubit in range(number_qubits):
						time_index[qubit] = applying_index
				else:			
					time_index.update({qubit: applying_index} for qubit in time_index.keys())
					# for qubit in time_index.keys():
					# 	time_index[qubit] = applying_index

			elif inst[0] == barrier:
				flag_barrier = True
				list_command = "{} {}".format(barrier, inst[1])
				applying_index = max(time_index[int(qubit)] for qubit in inst[1])

				for qubit in inst[1]: 
					time_index[qubit] = applying_index

			elif inst[0] in ["Qubit"]:
				if len(inst[1:]) == 2:
					qubit, size = inst[1:]
					list_command = "{} {} {}".format(inst[0], qubit, size)
				else:
					qubit = int(inst[1])
					list_command = "{} {}".format(inst[0], qubit)

			else:
				raise Exception("Error Happened : {}".format(inst))

		if flag_barrier: applying_index -= 1

		ordered_syscode[applying_index].append(list_command)

	return ordered_syscode


def transform_to_standardqasm(qasm, **kwargs):
	"""
		function to transform an openqasm to a standard qasm
	"""
	table_qubit_association = {}
	table_cbit_association = {}

	# unmatched_qubit_cbit : qubit <-> cbit
	unmatched_qubit_cbit = {}

	list_converted_code = []

	with open(qasm, "r") as infile:
		for line in infile:
			if any(item in line for item in ["OPENQASM", "qelib1.inc"]): continue

			tokens = parser.findall(line)[:-1]
			if not len(tokens): continue

			if tokens[0] in ["qreg", "creg"]:
				if tokens[0] == "qreg":
					qreg_name, qreg_size = tokens[1:]
					for i in range(int(qreg_size)):
						list_converted_code.append(["Qubit", "{}{}".format(qreg_name, str(i))])
						table_qubit_association[(qreg_name, i)] = "{}{}".format(qreg_name, str(i))
				else:
					creg_name, creg_size = tokens[1:]
					for i in range(int(creg_size)):
						list_converted_code.append(["Cbit", "{}{}".format(creg_name, str(i))])
						table_cbit_association[(creg_name, i)] = "{}{}".format(creg_name, str(i))

				continue

			converted_gate = gate_open2standard.get(tokens[0])
			
			if converted_gate in list_2q_gates:
				*angle, ctrl_name, ctrl_idx, trgt_name, trgt_idx = tokens[1:]
				
				ctrl_qubit = "{}{}".format(ctrl_name, ctrl_idx)
				trgt_qubit = "{}{}".format(trgt_name, trgt_idx)
				
				if len(angle):
					list_converted_code.append([converted_gate, angle[0], ctrl_qubit, trgt_qubit])
				else:
					list_converted_code.append([converted_gate, ctrl_qubit, trgt_qubit])


			elif converted_gate in list_1q_rotations:
				if converted_gate == u:
					angle = "{},{},{}".format(tokens[1], tokens[2], tokens[3])
				else:
					angle = tokens[1]
				
				qubit_argument = "{}{}".format(tokens[-2],tokens[-1])
				list_converted_code.append([converted_gate, angle, qubit_argument])


			# measure 먼저
			elif converted_gate in list_measure:
				# 만약, measure 값을 저장할 cbit register 정보가 주어지지 않으면,
				# 강제로, cbit_{qubit register index} 형태로 cbit 지정

				qubit_argument = "{}{}".format(tokens[1], tokens[2])
				if len(tokens) > 3:
					cbit_argument = "{}{}".format(tokens[4], tokens[5])
				else:
					cbit_argument = "cbit{}".format(tokens[2])

				list_converted_code.append([converted_gate, qubit_argument, "->", cbit_argument])


			# non-parametric 1qubit gates
			elif converted_gate in list_1q_gates:
				qubit_argument = "{}{}".format(tokens[1], tokens[2])
				list_converted_code.append([converted_gate, qubit_argument])

	return list_converted_code, table_qubit_association, table_cbit_association


def transform_to_openqasm(syscode, **kwargs):
	"""
		transform the system code with openqasm format
	"""	
	table_qubit_association = kwargs.get("qubit_association")
	table_cbit_association = kwargs.get("cbit_association")

	system_code = syscode.get("system_code")
	circuit = system_code.get("circuit")

	
	# update qubit_mapping based on the source code format and data
	inverse_qubit_association = {v: k for k, v in table_qubit_association.items()}
	
	for mapping in ["initial_mapping", "final_mapping"]:
		qubit_mapping = system_code.get(mapping)
		new_qubit_mapping = {}

		for k, v in qubit_mapping.items():
			reference = inverse_qubit_association.get(k)
			# only for the qubits defined in the source code
			if reference is not None:
				qubit_label = "{}[{}]".format(reference[0],reference[1])
				new_qubit_mapping[qubit_label] = v
		
		system_code[mapping] = new_qubit_mapping
	
	system_code["qubit"] = list(system_code.get("initial_mapping"))
	syscode["system_code"] = system_code
	
	# cbit_mapping
	inverse_cbit_association = {v: k for k, v in table_cbit_association.items()}
	cbit_mapping = {"{}[{}]".format(v[0], v[1]): k for k, v in inverse_cbit_association.items()}

	# circuit gate 
	new_circuit = collections.defaultdict(list)
	inverse_cbit_mapping = {v: k for k, v in cbit_mapping.items()}
	system_code["cbit"] = list(inverse_cbit_mapping.values())

	# 시스템 코드에서 명령 변환... standard -> open qasm
	for idx, list_instructions in circuit.items():
		for inst in list_instructions:
			tokens = parser.findall(inst)
			if not len(tokens): continue
			
			gate = tokens[0]
			converted_gate = gate_standard2open.get(gate)
			
			if gate == measz:
				cbit = inverse_cbit_mapping.get(tokens[-1])
				if cbit is not None:
					tokens[-1] = cbit

			tokens[0] = converted_gate
			new_circuit[idx].append(tokens)

	system_code["circuit"] = new_circuit
	syscode["system_code"] = system_code
	
	return syscode


# def transform_to_openqasm2(qasm, **kwargs):
# 	"""
# 		function to transform a quantum assembly code in standard format to an open qasm
# 		qubit labelling is also changed
# 		qasm : 
# 	"""
# 	number_qubits = kwargs.get("number_qubits")
# 	backend_name = kwargs.get("backend_name")

# 	if number_qubits is None:
# 		raise Exception("the size of quantum chip is not provided.")

# 	list_converted_code = []
	
# 	if backend_name in ["quito", "belem", "lima", "guadalupe"]:
# 		list_converted_code.append(["OPENQASM 2.0"])
# 	else:
# 		list_converted_code.append(["OPENQASM 3.0"])

# 	list_converted_code.append(["include \"qelib1.inc\""])

# 	list_converted_code.append(["qreg", "{}".format(number_qubits)])
# 	list_converted_code.append(["creg", "{}".format(number_qubits)])

# 	flag_measurement_appear = False

# 	if isinstance(qasm, dict):
# 		try:
# 			circuit_data = qasm["result"].get("system_code")
# 		except:
# 			circuit_data = qasm.get("system_code")

# 		if circuit_data is None:
# 			raise Exception("QASM : {}".format(qasm))

# 		quantum_circuit = circuit_data.get("circuit")
# 		try:
# 			circuit_depth = max(quantum_circuit.keys()) + 1
# 		except:
# 			if any(isinstance(item, str) for item in quantum_circuit.keys()):
# 				quantum_circuit = {int(k): v for k, v in quantum_circuit.items()}
		
# 			circuit_depth = max(quantum_circuit.keys()) + 1

# 		for i in range(circuit_depth):
# 			for inst in quantum_circuit[i]:
# 				token = parser.findall(inst)
# 				if token[0] in [g.str_gate_x, g.str_gate_z, g.str_gate_h, g.str_gate_y, g.str_gate_s, g.str_gate_sx]:
# 					list_converted_code.append([token[0].lower(), "{}".format(token[1])])

# 				elif token[0] in [g.str_gate_cnot, g.str_gate_cx, g.str_gate_cz, g.str_gate_swap]:
# 					list_converted_code.append([token[0].lower(), "{}".format(token[1]), "{}".format(token[2])])

# 				elif token[0] in [g.str_gate_cphase, g.str_gate_rzz]:
# 					list_converted_code.append([token[0].lower(), "{}".format(token[1]),
# 						"{}".format(token[2]), "{}".format(token[3])])

# 				elif token[0] in [g.str_gate_measz]:
# 					list_converted_code.append(["measure", "{}".format(token[1]), "{}".format(token[3])])

# 				elif token[0] in [g.str_gate_rx, g.str_gate_ry, g.str_gate_rz, g.str_gate_phase]:
# 					list_converted_code.append([token[0].lower(), "{}".format(token[1]), "{}".format(token[2])])

# 				elif token[0] in [g.str_gate_u]:
# 					list_converted_code.append([token[0].lower(), 
# 						"{}".format(token[1]), "{}".format(token[2]), "{}".format(token[3]), "{}".format(token[4])])
				
# 				elif token[0] in [g.str_barrier]:
# 					converted_inst = [token[0].lower()]
# 					for qubit in token[1:]:
# 						converted_inst.append(qubit)

# 					list_converted_code.append(converted_inst)

# 				else:
# 					raise Exception("inst : {}".format(inst))
# 	else:
# 		for inst in qasm:
# 			if inst[0] in ["PrepZ"]: 
# 				list_converted_code.append(["reset", "{}".format(inst[1])])

# 			elif inst[0] in ["X", "Z", "Y", "I", "H", "SX", "S"]:
# 				list_converted_code.append(["{}".format(inst[0].lower()), "{}".format(inst[1])])

# 			elif inst[0] in ["CNOT", "cx"]:
# 				list_converted_code.append(["cx", "{}".format(inst[1]), "{}".format(inst[2])])

# 			elif inst[0] in ["MeasZ"]:
# 				list_converted_code.append(["measure", "{}".format(inst[1]), "{}".format(inst[3])])

# 			elif inst[0] in ["Rx", "Rz", "Ry", "P"]:
# 				list_converted_code.append(["{}".format(inst[0].lower()), "{}".format(inst[1]), "{}".format(inst[2])])

# 			elif inst[0] in ["U"]:
# 				inst[0] = inst[0].lower()
# 				list_converted_code.append(inst)

# 			else:
# 				raise Exception("inst : {}".format(inst))

# 	return list_converted_code


# 필요한 함수인가 ?
# def extract_list_qubits2(syscode):
# 	list_qubits = []
# 	parser = re.compile("[\{a-zA-Z0-9_.*/\->\+}]+")

# 	if type(syscode) in [collections.defaultdict, dict]:
# 		for time_index, list_instrunctions in syscode.items():
# 			for instruction in list_instrunctions:
# 				token = parser.findall(instruction)
# 				if not len(token): continue

# 				if token[0] in ["CNOT", "cx", "cz"]: list_qubits.extend(map(int, token[1:3]))
# 				else:
# 					if token[0] in ["Rz", "rz", "Rx", "rx", "Ry", "ry", "P"]:
# 						angle, qubit = token[1:]
# 						list_qubits.append(int(qubit))
					
# 					elif token[0] in ["U"]:
# 						_, _, _, qubit = token[1:]
# 						list_qubits.append(int(qubit))
					
# 					elif token[0] not in ["Cbit"]:
# 						list_qubits.append(int(token[1]))

# 	return set(list_qubits)

# def transform_time_ordered_syscode(syscode, qubit_mapping):
	
# 	inverse_qubit_mapping = {v: k for k, v in qubit_mapping.items()}

# 	collections_qubits = collections.defaultdict(lambda: collections.defaultdict(bool))

# 	list_working_qubits = []
# 	circuit_index = 0

# 	collections_circuits = collections.defaultdict(object)
# 	circuit = collections.defaultdict(list)
# 	qubit_time_index = collections.defaultdict(int)

# 	for inst in syscode:
# 		if inst[0] in [g.str_gate_cnot, g.str_gate_cz, g.str_gate_swap]:
# 			ctrl, trgt = inst[1:]
# 			ctrl_name = inverse_qubit_mapping[ctrl]
# 			trgt_name = inverse_qubit_mapping[trgt]

# 			time_index = get_bigger(qubit_time_index[ctrl_name], qubit_time_index[trgt_name])
# 			circuit[time_index].append(inst)

# 			qubit_time_index[ctrl_name] = qubit_time_index[trgt_name] = time_index + 1

# 			if inst[0] == g.str_gate_swap:
# 				inverse_qubit_mapping[ctrl], inverse_qubit_mapping[trgt] =\
# 					inverse_qubit_mapping[trgt], inverse_qubit_mapping[ctrl]

# 		elif inst[0] in [g.str_gate_prepz, g.str_gate_prepx]:
# 			qubit_index = inst[1]
# 			qubit_name = inverse_qubit_mapping[qubit_index]

# 			qubit_type = qubit_name
# 			while qubit_type[-1].isdigit(): qubit_type = qubit_type[:-1]
# 			collections_qubits[qubit_type][qubit_name] = True

# 			circuit[qubit_time_index[qubit_name]].append(inst)
# 			qubit_time_index[qubit_name]+=1

# 		elif inst[0] in [g.str_gate_measz, g.str_gate_measx]:
# 			qubit_index = inst[1]
# 			qubit_name = inverse_qubit_mapping[qubit_index]

# 			qubit_type = qubit_name
# 			while qubit_type[-1].isdigit(): qubit_type = qubit_type[:-1]
# 			collections_qubits[qubit_type][qubit_name] = False

# 			circuit[qubit_time_index[qubit_name]].append(inst)
# 			qubit_time_index[qubit_name] += 1

# 			if not any(collections_qubits[qubit_type].values()): 
# 				collections_circuits[circuit_index] = circuit
# 				for qubit in qubit_time_index.keys(): qubit_time_index[qubit] = 0

# 				circuit = collections.defaultdict(list)
# 				circuit_index+=1
# 		else:
# 			qubit_index = inst[1]
# 			qubit_name = inverse_qubit_mapping[qubit_index]

# 			circuit[qubit_time_index[qubit_name]].append(inst)
# 			qubit_time_index[qubit_name]+=1

# 	inverse_qubit_mapping = {v: k for k, v in qubit_mapping.items()}
# 	for circuit_idx, circuit in collections_circuits.items():
# 		for time_idx, instructions in circuit.items():
# 			for inst in instructions:
# 				if inst[0] in [g.str_gate_swap, g.str_gate_cnot]:
# 					print(inst, inst[0], inverse_qubit_mapping[inst[1]], inverse_qubit_mapping[inst[2]])

# 					if inst[0] == g.str_gate_swap:
# 						inverse_qubit_mapping[inst[1]], inverse_qubit_mapping[inst[2]] =\
# 							inverse_qubit_mapping[inst[2]], inverse_qubit_mapping[inst[1]]
# 				else:
# 					print(inst, inst[0], inverse_qubit_mapping[inst[1]])
# 		print("\n")

# 	pprint(collections_circuits)




# def preanalyze_qasm(qasm):
# 	'''
# 		commutable cnot swap 함수
# 		args: qasm in list data structure
# 	'''
# 	import re
# 	p = re.compile("[\{a-zA-Z0-9_.*\->\}]+")

# 	idx = 0

# 	# list_qasm -> idx: qasm instruction
# 	list_qasm = {}

# 	# monitoring_CNOT -> [{"control": ..., "target": ..., "idx": idx}, {..}]
# 	list_monitoring_CNOT = []

# 	# qubit -> qubit : [{"qubit": .., "cnots": [.. ]}]
# 	list_monitoring_qubit = []

# 	list_memory_declaration = []
# 	# with open(file_qasm, "r") as infile:
	
# 	for command in qasm:
# 		# 일단 순서대로 list에 저장
# 		list_qasm[idx] = command

# 		if command[0] in ["Qubit", "Cbit"]: 
# 			idx += 1
# 			continue

# 		# 게이트 타입에 따라 연산 작업
# 		# case 1: CNOT
# 		if command[0] == "CNOT":
# 			ctrl, trgt = command[1:]
# 			if ctrl == trgt:
# 				raise error.Error("For CNOT, ctrl qubit and trgt qubit are the same. {}".format(command))

# 			# 새로운 CNOT의 ctrl, trgt 큐빗이 기존 CNOT에 물려 있는지 확인
# 			# check 1: ctrl 만 공유되는 상황 
# 			# check 2: ctrl 과 trgt 가 공유되는 상황
# 			# check 3: check 1과 2는 아니지만, 어떻게는 큐빗 공유는 이루어지는 상황
# 			# check 4: 새로운 cnot 과 앞선 cnot 들이 완전히 독립적인 상황

# 			flag_checking = {k: False for k in range(3)}
# 			flag_new = True
# 			for cnot in list_monitoring_CNOT:
# 				if ctrl == cnot["control"]:
# 					if trgt == cnot["target"]:
# 						print("기존 CNOT: {}, 신규 CNOT : {}".format(cnot, (ctrl, trgt)))
# 						print("중복 발생")
# 						list_monitoring_CNOT.remove(cnot)
# 						# list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})
# 						flag_new = False
					
# 					else:
# 						print("기존 CNOT: {}, 신규 CNOT : {}".format(cnot, (ctrl, trgt)))
# 						print("교환 대상")
# 						list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})
# 						list_monitoring_qubit.append({"qubit": trgt, 
# 													  "cnot_list": [cnot, {"control": ctrl, "target": trgt, "id": idx}]})
# 						flag_new = False

# 				else:
# 					if ctrl == cnot["target"] or trgt in list(cnot.values()):
# 						print("기존 CNOT: {}, 신규 CNOT : {}".format(cnot, (ctrl, trgt)))
# 						print("기존 CNOT 제거 & 신규 CNOT 추가")
# 						list_monitoring_CNOT.remove(cnot)
# 						list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})
# 						flag_new = False
# 				break

# 			if flag_new:
# 				print("독립적인 CNOT 게이트")
# 				list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})

# 			print("list of monitoring cnot: ")
# 			pprint(list_monitoring_CNOT)

# 		# case 2: one-qubit gate
# 		else:
# 			# case 2-1: rotational gate with arguments: angle and qubit
# 			if command[0] in ["Rz", "rz"]: qubit = command[2]
# 			# case 2-2: the other one-qubit gate with argument: qubit
# 			else: qubit = command[1]

# 			if not len(list_monitoring_qubit):
# 				for cnot in list_monitoring_CNOT:
# 					if cnot["control"] == qubit or cnot["target"] == qubit:
# 						list_monitoring_CNOT.remove(cnot)

# 			pprint(list_monitoring_qubit)
# 			for q in list_monitoring_qubit:
# 				# check 1: 공통 control 큐빗에 one-qubit 게이트가 인가된 경우
# 				if (q["cnot_list"][0]["control"] == qubit) and (q["cnot_list"][1]["control"] == qubit):
# 					print("case 1-1.. ")
# 					former, latter = q["cnot_list"][:]
# 					list_monitoring_CNOT.remove(former)
# 					list_monitoring_CNOT.remove(latter)
# 					list_monitoring_qubit.remove(q)

# 				# # check 2: 첫번째 cnot의 타겟 큐빗에 one-qubit 게이트가 인가된 경우
# 				# if q["cnot_list"][0]["target"] == qubit:
# 				# 	pprint(q["cnot_list"])
# 				# 	former, latter = q["cnot_list"][:]
# 				# 	list_monitoring_CNOT.remove(former)
# 				# 	list_monitoring_qubit.remove(q)
# 				# check 3
# 				if q["qubit"] == qubit:
# 					print("교환하자... {}".format(q["cnot_list"]))
# 					print("교환 이전: ")
# 					pprint(list_qasm)
# 					former, latter = q["cnot_list"][:]
# 					list_qasm[former["id"]] = ["CNOT", latter["control"], latter["target"]]
# 					list_qasm[latter["id"]] = ["CNOT", former["control"], former["target"]]
					
# 					print("교환 결과: ")
# 					pprint(list_qasm)
# 					# list_monitoring_CNOT.remove(former)
# 					list_monitoring_CNOT.remove(latter)
# 					list_monitoring_qubit.remove(q)

# 					print("monitoring cnot: ")
# 					pprint(list_monitoring_CNOT)
# 					print("monitoring qubit: ")
# 					pprint(list_monitoring_qubit)

# 		idx += 1

# 	return list(list_qasm.values())


# def preanalyze_qasm_file(file_qasm):
# 	'''
# 		commutable cnot swap 함수
# 		args: qasm in file
# 	'''

# 	import re
# 	p = re.compile("[\{a-zA-Z0-9_.*\->\}]+")

# 	idx = 0

# 	# list_qasm -> idx: qasm instruction
# 	list_qasm = {}

# 	# monitoring_CNOT -> [{"control": ..., "target": ..., "idx": idx}, {..}]
# 	list_monitoring_CNOT = []

# 	# qubit -> qubit : [{"qubit": .., "cnots": [.. ]}]
# 	list_monitoring_qubit = []

# 	list_memory_declaration = []
# 	with open(file_qasm, "r") as infile:
# 		for line in infile:
# 			tokens = p.findall(line)
# 			# print(tokens)
# 			if not len(tokens):  continue

# 			# 일단 순서대로 list에 저장
# 			list_qasm[idx] = tokens

# 			if tokens[0] in ["Qubit", "Cbit"]: 
# 				idx += 1
# 				continue

# 			# 게이트 타입에 따라 연산 작업
# 			# case 1: CNOT
# 			if tokens[0] == "CNOT":
# 				ctrl, trgt = tokens[1:]
# 				if ctrl == trgt:
# 					raise error.Error("For CNOT, ctrl qubit and trgt qubit are the same. {}".format(tokens))
# 				# 새로운 CNOT의 ctrl, trgt 큐빗이 기존 CNOT에 물려 있는지 확인
# 				# check 1: ctrl 만 공유되는 상황 
# 				# check 2: ctrl 과 trgt 가 공유되는 상황
# 				# check 3: check 1과 2는 아니지만, 어떻게는 큐빗 공유는 이루어지는 상황
# 				# check 4: 새로운 cnot 과 앞선 cnot 들이 완전히 독립적인 상황

# 				flag_checking = {k: False for k in range(3)}
# 				flag_new = True
# 				for cnot in list_monitoring_CNOT:
# 					if ctrl == cnot["control"]:
# 						if trgt == cnot["target"]:
# 							print("기존 CNOT: {}, 신규 CNOT : {}".format(cnot, (ctrl, trgt)))
# 							print("중복 발생")
# 							list_monitoring_CNOT.remove(cnot)
# 							# list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})
# 							flag_new = False
						
# 						else:
# 							print("기존 CNOT: {}, 신규 CNOT : {}".format(cnot, (ctrl, trgt)))
# 							print("교환 대상")
# 							list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})
# 							list_monitoring_qubit.append({"qubit": trgt, 
# 														  "cnot_list": [cnot, {"control": ctrl, "target": trgt, "id": idx}]})
# 							flag_new = False

# 					else:
# 						if ctrl == cnot["target"] or trgt in list(cnot.values()):
# 							print("기존 CNOT: {}, 신규 CNOT : {}".format(cnot, (ctrl, trgt)))
# 							print("기존 CNOT 제거 & 신규 CNOT 추가")
# 							list_monitoring_CNOT.remove(cnot)
# 							list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})
# 							flag_new = False
# 					break

# 				if flag_new:
# 					print("독립적인 CNOT 게이트")
# 					list_monitoring_CNOT.append({"control": ctrl, "target": trgt, "id": idx})

# 				print("list of monitoring cnot: ")
# 				pprint(list_monitoring_CNOT)

# 			# case 2: one-qubit gate
# 			else:
# 				# case 2-1: rotational gate with arguments: angle and qubit
# 				if tokens[0] in ["Rz", "rz"]: qubit = tokens[2]
# 				# case 2-2: the other one-qubit gate with argument: qubit
# 				else: qubit = tokens[1]

# 				if not len(list_monitoring_qubit):
# 					for cnot in list_monitoring_CNOT:
# 						if cnot["control"] == qubit or cnot["target"] == qubit:
# 							list_monitoring_CNOT.remove(cnot)

# 				for q in list_monitoring_qubit:
# 					# check 1: 공통 control 큐빗에 one-qubit 게이트가 인가된 경우
# 					if (q["cnot_list"][0]["control"] == qubit) and (q["cnot_list"][1]["control"] == qubit):
# 						print("case 1-1.. ")
# 						former, latter = q["cnot_list"][:]
# 						list_monitoring_CNOT.remove(former)
# 						list_monitoring_CNOT.remove(latter)
# 						list_monitoring_qubit.remove(q)

# 					# # check 2: 첫번째 cnot의 타겟 큐빗에 one-qubit 게이트가 인가된 경우
# 					# if q["cnot_list"][0]["target"] == qubit:
# 					# 	pprint(q["cnot_list"])
# 					# 	former, latter = q["cnot_list"][:]
# 					# 	list_monitoring_CNOT.remove(former)
# 					# 	list_monitoring_qubit.remove(q)
# 					# check 3
# 					if q["qubit"] == qubit:
# 						print("교환하자... {}".format(q["cnot_list"]))
# 						print("교환 이전: ")
# 						pprint(list_qasm)
# 						former, latter = q["cnot_list"][:]
# 						list_qasm[former["id"]] = ["CNOT", latter["control"], latter["target"]]
# 						list_qasm[latter["id"]] = ["CNOT", former["control"], former["target"]]
						
# 						print("교환 결과: ")
# 						pprint(list_qasm)
# 						# list_monitoring_CNOT.remove(former)
# 						list_monitoring_CNOT.remove(latter)
# 						list_monitoring_qubit.remove(q)

# 						print("monitoring cnot: ")
# 						pprint(list_monitoring_CNOT)
# 						print("monitoring qubit: ")
# 						pprint(list_monitoring_qubit)

# 			idx += 1

# 	with open(file_qasm, "w") as outfile:
# 		for idx, inst in list_qasm.items():
# 			if inst[0] in ["Qubit", "Cbit"]:
# 				str_command = "{} {}\n".format(inst[0], inst[1])

# 			elif inst[0] == "CNOT":
# 				str_command = "{} {},{}\n".format(inst[0], inst[1], inst[2])
			
# 			elif inst[0] in ["Rz", "rz"]: 
# 				str_command = "{}({}) {}\n".format(inst[0], inst[1], inst[2])
			
# 			elif inst[0] in ["MeasZ"]:
# 				str_command = "{} {} -> {}\n".format(inst[0], inst[1], inst[3])	
			
# 			else:
# 				str_command = "{} {}\n".format(inst[0], inst[1])

# 			outfile.write(str_command)

# 	return file_qasm
