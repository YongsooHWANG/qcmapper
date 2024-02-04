
import collections
import re
import copy
import itertools
from math import *

from pprint import pprint
from gatelist import *

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
			if len(table[inst[-1]]) and len(table[inst[-2]]):
				last_instA = table[inst[2]][-1]
				last_instB = table[inst[3]][-1]

				conditionA = (last_instA["gate"] == gate) and (last_instA["qubits"] == inst[2:])
				conditionB = (last_instB["gate"] == gate) and (last_instB["qubits"] == inst[2:])
				
				# 동일하면
				if conditionA and conditionB:
					# 이전 연산의 angle 조정
					new_angle = float(eval(last_instA["angle"])) + float(eval(inst[1]))
					table[inst[-1]][-1]["angle"] = str(new_angle)
					table[inst[-2]][-1]["angle"] = str(new_angle)

				# 다르면
				else:
					table[inst[-1]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})
					table[inst[-2]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})

			else:
				table[inst[-1]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})
				table[inst[-2]].append({"gate": gate, "angle": inst[1], "qubits": inst[2:], "idx": idx})


		elif any(item in gate for item in list_2q_gates + ["Relabel"]):
			if len(table[inst[1]]) and len(table[inst[2]]):
				# 새로운 명령과 이전 명령 비교
				last_instA = table[inst[-1]][-1]
				last_instB = table[inst[-2]][-1]

				conditionA = (last_instA["gate"] == gate) and (last_instA["qubits"] == inst[1:])
				conditionB = (last_instB["gate"] == gate) and (last_instB["qubits"] == inst[1:])

				# 동일하면
				if conditionA and conditionB:
					table[inst[-1]].pop()
					table[inst[-2]].pop()

				# 다르면
				else:
					table[inst[-1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
					table[inst[-2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

			else:
				table[inst[-1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
				table[inst[-2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

		# 2-qubit gate : typeB (control, target 큐빗 명시가 필요없는 게이트)
		elif swap in gate:
			if len(table[inst[-1]]) and len(table[inst[-2]]):
				# 새로운 명령과 이전 명령 비교
				last_instA = table[inst[-1]][-1]
				last_instB = table[inst[-2]][-1]

				conditionA = (last_instA["gate"] == gate) and (set(last_instA["qubits"]) == set(inst[1:]))
				conditionB = (last_instB["gate"] == gate) and (set(last_instB["qubits"]) == set(inst[1:]))

				# 동일하면
				if conditionA and conditionB:
					table[inst[-1]].pop()
					table[inst[-2]].pop()

				# 다르면
				else:
					table[inst[-1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
					table[inst[-2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

			else:
				table[inst[-1]].append({"gate": gate, "qubits": inst[1:], "idx": idx})
				table[inst[-2]].append({"gate": gate, "qubits": inst[1:], "idx": idx})

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

		elif gate in list_1q_rotations:
			if len(table[inst[-1]]):
				last_inst = table[inst[-1]][-1]

				# 동일한 회전 게이트가 연속되면, angle 확인 후, 앞선 게이트의 angle 값을 변경
				if (last_inst["gate"] == gate) and (last_inst["qubits"] == inst[-1]):
					angle_in_last_inst = last_inst["angle"]
					
					if gate in [u, u3]:
						new_angle = (float(eval(last_inst["angle"][0]))+ float(eval(inst[1][0])),
									float(eval(last_inst["angle"][1]))+ float(eval(inst[1][1])),
									float(eval(last_inst["angle"][2]))+ float(eval(inst[1][2])))
					
					elif gate in [u2]:
						new_angle = (float(eval(last_inst["angle"][0]))+ float(eval(inst[1][0])),
									float(eval(last_inst["angle"][1]))+ float(eval(inst[1][1])))

					else:
						new_angle = float(eval(last_inst["angle"])) + float(eval(inst[1]))

					table[inst[-1]][-1]["angle"] = str(new_angle)

				else:
					table[inst[-1]].append({"gate": gate, "qubits": [inst[-1]], "idx": idx, "angle": inst[1:-1]})
			else:
				table[inst[-1]].append({"gate": gate, "qubits": [inst[-1]], "idx": idx, "angle": inst[1:-1]})


		elif gate in list_1q_gates:
			# 새로운 명령과 이전 명령 비교
			if len(table[inst[-1]]):
				last_inst = table[inst[-1]][-1]

				# 동일하면
				if (last_inst["gate"] == gate) and (last_inst["qubits"] == inst[-1]):
					table[inst[-1]].pop()

				# 다르면
				else:
					table[inst[-1]].append({"gate": gate, "qubits": [inst[-1]], "idx": idx})	

			else:
				table[inst[-1]].append({"gate": gate, "qubits": [inst[-1]], "idx": idx})

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

		elif any(item in gate for item in list_2q_gates + ["Relabel"]):
			post_processed_syscode.append([gate, v["qubits"][0], v["qubits"][1]])

		elif gate in list_measure:
			post_processed_syscode.append([gate, v["qubits"][0], v["cbits"][0]])

		elif gate in list_1q_rotations:
			command = [gate] + v.get("angle") + v.get("qubits")
			post_processed_syscode.append(command)

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


		elif any(item in inst[0] for item in list_2q_gates + ["Relabel"]):
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
				*angle, qubit = inst[1:]
				list_command = "{}({}) {}".format(inst[0], ",".join(angle), qubit)

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
				list_angles = [angle for angle in tokens[1:-2]]
				qubit_argument = "{}{}".format(tokens[-2],tokens[-1])
				
				list_angles.insert(0, converted_gate)
				list_angles.append(qubit_argument)
				
				list_converted_code.append(list_angles)


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

	instruction_format = kwargs.get("format")

	system_code = syscode.get("system_code")
	circuit = system_code.get("circuit")

	# flag_gate_angle : 
	# 		True : U(theta, phi, lambda) qubit
	# 		False : U, theta, phi, lambda, qubit
	flag_gate_angle = instruction_format.get("gate_angle")
	
	# flag_str_instruction:
	# 		True : [.., "U, theta, phi, lambda, qubit", ..]
	# 		False : [.., ["U", "theta", "phi", "lambda", "qubit"]]
	flag_str_instruction = instruction_format.get("str_instruction")

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
			
			else:
				new_qubit_mapping[k] = v
		
		system_code[mapping] = new_qubit_mapping
	
	system_code["qubit"] = list(item for item in system_code.get("initial_mapping") if "dummy" not in item)

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

			if flag_gate_angle:
				if gate in list_1q_rotations: 
					angle = ",".join(tokens[1:-1])
					tokens = ["{}({})".format(converted_gate, angle), tokens[-1]]
				
				elif gate in list_2q_rotations:
					tokens = ["{}({})".format(converted_gate, tokens[1])] + tokens[-2:]
				
				else:		
					tokens[0] = converted_gate
			else:		
				tokens[0] = converted_gate
			
			new_circuit[idx].append(tokens)

	system_code["circuit"] = new_circuit
	syscode["system_code"] = system_code
	
	return syscode

