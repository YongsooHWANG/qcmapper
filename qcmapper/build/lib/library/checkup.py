
import re
from pprint import pprint
from library.gatelist import *


def checkup_system_code(system_code, mapping_table, qchip, **kwargs):
	"""
		시스템 코드 checkup 함수
		합성 결과가 qubit connectivity 를 만족하는 지 확인하는 함수
	"""
	
	parser = re.compile(r"[\{a-zA-Z0-9_.*\->\+}]+")
	flag_verification = True
	
	if type(system_code) == dict:
		for instructions in list(system_code.values()):
			for inst in instructions:
				tokens = parser.findall(inst)
				if not len(tokens): continue
				
				if tokens[0] in list_2q_gates:
					*angle, ctrl, trgt = tokens[1:]
					ctrl = int(ctrl)
					trgt = int(trgt)

					if trgt not in qchip["qubit_connectivity"][ctrl]:
						flag_verification = False
						break

	elif type(system_code) == list:
		for inst in system_code:
			if inst[0] in list_2q_gates:
				*angle, ctrl, trgt = inst[1:]
				ctrl = int(ctrl)
				trgt = int(trgt)

				if trgt not in qchip["qubit_connectivity"][ctrl]:
					flag_verification = False
					break

	return flag_verification