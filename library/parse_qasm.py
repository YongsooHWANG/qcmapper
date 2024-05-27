
import os
import re
import collections
from math import *

from gatelist import *
from icecream import ic
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
    
