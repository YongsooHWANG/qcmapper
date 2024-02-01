# -*-coding:utf-8-*-

import collections
import itertools

from pprint import pprint

from library.gatelist import *
import library.depth_analysis as depth_analysis

get_bigger = lambda a, b: a if a > b else b


class QuantumCircuit:
	def __init__ (self, **kwargs):

		if "layout_size" in kwargs: 
			self.device_size = kwargs["layout_size"]
		
		if "qubit_table" in kwargs: 
			qubit_table = kwargs["qubit_table"]
		
		self.synthesis_option = kwargs["option"]
		
		self.flag_swap = self.synthesis_option.get("allow_swap")
		if self.flag_swap is None:
			self.flag_swap = False

		if "parallel_swap" in self.synthesis_option:
			self.parallel_swap = self.synthesis_option["parallel_swap"]
		else:
			self.parallel_swap = False

		self.logical_gate_performance = None
		if "gate" in kwargs:
			self.logical_gate_performance = kwargs["gate"]

		self.inverse_mapping_table = {v: k for k, v in qubit_table.items()}
		
		self.system_code = {"cbit": [], "circuit": [], "initial_mapping": qubit_table}
		
		algorithm_qubits = kwargs.get("algorithm_qubits")
		
		self.flag_algorithm_qubits = False
		
		if algorithm_qubits is not None:
			self.system_code["qubit"] = algorithm_qubits
			self.flag_algorithm_qubits = True
		
		self.list_algorithm_qubits = set([])
			
		
		# for i in range(self.device_size):
		# 	self.system_code["qubit"].append("qbit{}".format(i))
		# 	self.system_code["cbit"].append("bit{}".format(i))

		self.function_list = collections.defaultdict(int)
					
		# qubit performance table: [qubit_idx, gate_count, computing_cycle, computing_time]
		# plain array than numpy array
		# self.qubit_array = [[0, 0, 0, 0] for i in range(self.device_size)]


	def get_qubit_association(self):
		return self.qubit_association

	def register(self, register_type, qubit, **kwargs):
		self.system_code["circuit"].append([register_type, qubit])

	def one_qubit_gate(self, gate, qubit, **kwargs):
		if not self.flag_algorithm_qubits:
			self.list_algorithm_qubits.add(qubit)

		# system code 생성 모드이면..
		if gate in [u]:
			theta = kwargs.get("theta")
			phi = kwargs.get("phi")
			lamda = kwargs.get("lamda")

			if any(item is None for item in [theta, phi, lamda]): 
				raise Exception("Error Happened for angle ")
			else:
				_str_code = [gate, theta, phi, lamda, qubit]

		elif gate in list_1q_rotations:
			angle = kwargs.get("angle")
			_str_code = [gate, str(angle), qubit]

		elif gate in list_measure:
			if "cbit" in kwargs: 
				_cbit = kwargs["cbit"]
				self.system_code["cbit"].append(kwargs["cbit"])
				# _str_code = [gate, str(qubit), str(qubit)]
				_str_code = [gate, qubit, _cbit]
			else:
				# _str_code = [gate, str(qubit), str(qubit)]
				_str_code = [gate, qubit, qubit]
		
		else:
			# _str_code = [gate, str(qubit)]
			_str_code = [gate, qubit]
		
		self.system_code["circuit"].append(_str_code)


	def swap(self, path, **kwargs):
		_qubit1, _qubit2, _weight = path[:]
		
		calibration_type = kwargs.get("calibration_type")

		if self.flag_swap:
			# self.system_code["circuit"].append([g.str_gate_swap, str(_qubit1), str(_qubit2)])
			self.system_code["circuit"].append([swap, _qubit1, _qubit2])

		else:
			'''
				allow_swap 이 false 인 경우,
				swap a, b = cnot a, b x cnot b, a x cnot a, b 또는
							cnot b, a x cnot a, b x cnot b, a 
				의 두가지 경우로 구현 가능함 cnot a, b 와 cnot b, a 의 성능을 비교해서, 선택함
			'''
			if calibration_type in [None, "depth"]:
				if _weight == 1:
					# ctrl -> trgt
					self.cnot([_qubit1, _qubit2])
					
					self.one_qubit_gate(h, _qubit1)
					self.one_qubit_gate(h, _qubit2)
					
					self.cnot([_qubit1, _qubit2])

					self.one_qubit_gate(h, _qubit1)
					self.one_qubit_gate(h, _qubit2)
					
					self.cnot([_qubit1, _qubit2])

				elif _weight == 3:
					self.cnot([_qubit2, _qubit1])

					self.one_qubit_gate(h, _qubit1)
					self.one_qubit_gate(h, _qubit2)

					self.cnot([_qubit2, _qubit1])

					self.one_qubit_gate(h, _qubit1)
					self.one_qubit_gate(h, _qubit2)

					self.cnot([_qubit2, _qubit1])

				elif _weight == 4:
					self.cnot([_qubit1, _qubit2])
					self.cnot([_qubit2, _qubit1])
					self.cnot([_qubit1, _qubit2])

			else:
				weight_graph = kwargs.get("weight_graph")
				ctrl, trgt = map(int, weight_graph[_qubit1][_qubit2]["direction"].split(">"))

				self.cnot([ctrl, trgt])
				self.cnot([trgt, ctrl])
				self.cnot([ctrl, trgt])

		qubit_table = kwargs.get("qubit_table")
		if qubit_table is not None:
			inverse_mapping_table = {v: k for k, v in qubit_table.items()}

			inverse_mapping_table[_qubit1], inverse_mapping_table[_qubit2] =\
				inverse_mapping_table[_qubit2], inverse_mapping_table[_qubit1]

			return {v: k for k, v in inverse_mapping_table.items()}


	def manage_cnot(self, qubits, **kwargs):

		list_path = kwargs.get("path")
		qubit_mapping_table = kwargs.get("qubit_table")
		final_operation = kwargs.get("final_operation")
		calibration_type = kwargs.get("calibration_type")
		angle = kwargs.get("angle")

		if any(item is None for item in [list_path, qubit_mapping_table, final_operation]):
			raise Exception("stop Error ")

		if not self.flag_algorithm_qubits:
			self.list_algorithm_qubits.add(qubits[0])
			self.list_algorithm_qubits.add(qubits[1])

		# cnot path 의 사이즈가 1보다 크뎐, swap 이 필요함
		if len(list_path)>1:
			# cnot path 상에서 병렬적 swap 이 가능한 경우
			# 병렬적 swap : forward path: cnot path 의 처음부터 순서대로 path/2 까지
			#			: backward path: cnot path 의 마지막 부터 역순으로 path/2 까지
			if self.parallel_swap:

				length_swap_chain = len(list_path)
				# 병렬 swap 을 수행하기 위해서는 swap chain 이 길이에 따라서 다음과 같이 동작한다.
				# even 일 경우 : a - b - c - d - e
				# a<->b, e<->d 이후 b<->c 는 swap 수행, 큐빗 d 는 대기 후 마지막에 cnot c,d 수행
				# 따라서 forward path 는 list_path 의 절반, backward path 는 list_path 후반의 의 절반-1. 
				if length_swap_chain%2 == 0:
					half_length = int(length_swap_chain/2)
					
					forward_path = list_path[:half_length]
					backward_path = list(reversed(list_path[half_length+1:]))
					final_edge = list_path[half_length:half_length+1][0]
				else:
					half_length = int(length_swap_chain/2)
					forward_path = list_path[:half_length]
					backward_path = list(reversed(list_path))[:half_length]
					final_edge = list_path[half_length:half_length+1][0]

				for item in itertools.zip_longest(forward_path, backward_path):
					if item[0] is not None:
						qubit_mapping_table = self.swap(item[0], qubit_table=qubit_mapping_table,
																weight_graph=kwargs.get("weight_graph"),
																calibration_type=calibration_type)
					if item[1] is not None:
						qubit_mapping_table = self.swap(item[1], qubit_table=qubit_mapping_table,
																weight_graph=kwargs.get("weight_graph"),
																calibration_type=calibration_type)
			
			# 병렬적 swap이 허용되지 않으면, cnot path 상에서 순서대로 serial swap 수행함
			# cnot path 상의 마지막 item 은 최종 cnot 수행할 큐빗 pair
			else:
				final_edge = list_path[-1]
				# swap_path = list_path[:-1]

				# swap path 상에서 순서대로 swap 수행함
				for edge in list_path[:-1]:
					qubit_mapping_table = self.swap(edge, qubit_table=qubit_mapping_table,
															weight_graph=kwargs.get("weight_graph"),
															calibration_type=calibration_type)

				self.inverse_mapping_table = {v: k for k, v in qubit_mapping_table.items()}

		else:
			# cnot path 테이블의 길이가 1이면, 바로 cnot 수행하면 된다.
			final_edge = list_path[-1]
		
		# 최종 연산 수행 : 
		# final_edge 의 마지막 value -> weight
		# case 1: cnot
		if final_operation == cnot:  self.cnot(final_edge[:2])

		# case 2: swap
		elif final_operation == swap: 
			qubit_mapping_table = self.swap(final_edge, qubit_table=qubit_mapping_table,
														weight_graph=kwargs.get("weight_graph"),
														calibration_type=calibration_type)
		elif final_operation == cphase:
			self.cphase(final_edge[:2], angle=angle)

		# 원거리상의 2-큐빗 게이트 수행 종료 후의 큐빗 매핑 테이블
		self.inverse_mapping_table = {v: k for k, v in qubit_mapping_table.items()}
		return qubit_mapping_table


	def cphase(self, qubits, **kwargs):
		_ctrl, _trgt = qubits[:]
		angle = kwargs.get("angle")

		# self.system_code["circuit"].append([g.str_gate_cphase, angle, str(_ctrl), str(_trgt)])
		self.system_code["circuit"].append([cphase, angle, _ctrl, _trgt])

	def rzz(self, qubits, **kwargs):
		_ctrl, _trgt = qubits[:]
		angle = kwargs.get("angle")

		# self.system_code["circuit"].append([g.str_gate_rzz, angle, str(_ctrl), str(_trgt)])
		self.system_code["circuit"].append([rzz, angle, _ctrl, _trgt])


	def cnot(self, qubits):
		_ctrl, _trgt = qubits[:]
		# self.system_code["circuit"].append([g.str_gate_cnot, str(_ctrl), str(_trgt)])
		self.system_code["circuit"].append([cnot, _ctrl, _trgt])


	def insert_barrier_all(self):
		# 시스템 코드에 삽입
		# 실제로 양자 게이트를 추가하는 것은 아니기 때문에, 각 큐빗 array 의 항목들을 increment 하는 것은 아님.

		_str_code = [barrier_all]
		self.system_code["circuit"].append(_str_code)


	def insert_barrier(self, list_qubits):
		pass

	def get_analysis(self):

		# logical_circuit_depth = max([self.qubit_array[i][1] for i in range(self.device_size)])+1
		# physical_circuit_depth = max([self.qubit_array[i][3] for i in range(self.device_size)])+1
		# computing_time = max([self.qubit_array[i][2] for i in range(self.device_size)])
		# computing_time = round(computing_time, 16)
		
		gate_depth_analysis = {}
		
		flag_t_depth = self.synthesis_option.get("t_depth")
		flag_cnot_depth = self.synthesis_option.get("cnot_depth")

		if flag_t_depth:
		# t 게이트 depth 계산 함수 호출
			if "T" in self.function_list.keys():
				t_depth = depth_analysis.evaluate_t_depth(self.system_code)
			else:
				t_depth = 0
			gate_depth_analysis.update({"T-Gate": t_depth})

		if flag_cnot_depth:
			if "CNOT" in self.function_list.keys():	
				cnot_depth = depth_analysis.evaluate_cnot_depth(self.system_code)
			else:
				cnot_depth = 0
			gate_depth_analysis.update({"CNOT": cnot_depth})

		if not self.flag_algorithm_qubits:
			self.system_code["qubit"] = list(self.list_algorithm_qubits)

		return {"Algorithm Qubits": len(self.system_code["qubit"]),
			    "Gate Depth": gate_depth_analysis}


	def get_system_code(self, **kwargs):

		if "file" in kwargs: 
			_path_syscode = kwargs["file"]
			
			with open(_path_syscode, "w") as outfile:
				for i in self.system_code:
					outfile.write("{}\n".format(i))

		self.system_code["cbit"] = list(set(self.system_code["cbit"]))
		self.system_code["initial_mapping"] = {k: v for k, v in self.system_code["initial_mapping"].items() if "dummy" not in k}
		self.system_code["final_mapping"] = {v: k for k, v in self.inverse_mapping_table.items() if "dummy" not in v}

		return {"system_code": self.system_code}