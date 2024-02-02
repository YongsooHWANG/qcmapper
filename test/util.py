
# -*-coding:utf-8-*-

import os
import sys
import simplejson as json
from pprint import pprint
import time

import math
import pandas
import numpy as np

from icecream import ic

def display_qubit_movements(system_code, lattice_size, **kwargs):
	
	# initial mapping
	circuit = system_code.get("circuit")
	qubit_mapping = system_code.get("initial_mapping")
	inverse_mapping = {v: k for k, v in qubit_mapping.items()}
	
	layout = [[0 for i in range(lattice_size["width"])] for j in range(lattice_size["height"])]
	
	for idx, qubit in inverse_mapping.items():
		x_coord = int(idx/lattice_size["width"])
		z_coord = int(idx%lattice_size["width"])

		layout[x_coord][z_coord] = qubit
	
	print(" =====================================================  ")
	print("Initial Mapping: ")
	print(" -----------------------------------------------------  ")
	print(pandas.DataFrame(layout).to_string())
	print(" =====================================================  ")


	circuit_depth = len(circuit)
	# circuit
	for idx in range(circuit_depth):
		instructions = circuit.get(idx)

		flag_swap = False
		print(" =====================================================  ")
		print("instructions at {}-th index : {}".format(idx, instructions))
		print(" -----------------------------------------------------  ")
		
		for inst in instructions:
			flag_swap = False	
			if inst[0] in ["swap", "cx", "cp", "rzz", "cz"]:
				*angle, ctrl, trgt = inst[1:]
				ctrl, trgt = map(int, [ctrl, trgt])
				

				print(" {} qubits ({}, {}) -> ({}, {})".format(inst[0], ctrl, trgt, inverse_mapping.get(ctrl), inverse_mapping.get(trgt)))

				if inst[0] == "swap": 
					flag_swap = True
					inverse_mapping[ctrl], inverse_mapping[trgt] =\
						inverse_mapping[trgt], inverse_mapping[ctrl]

			elif inst[0] in ["measure"]:
				qubit = int(inst[1])
				print(" {} ({}) -> {}".format(inst[0], qubit, inverse_mapping.get(qubit)))

			else:
				qubit = int(inst[-1])
				print(" {} ({}) -> {}".format(inst[0], qubit, inverse_mapping.get(qubit)))


			# 2d array 재 구성
			for idx, qubit in inverse_mapping.items():
				x_coord = int(idx/lattice_size["width"])
				z_coord = int(idx%lattice_size["width"])

				layout[x_coord][z_coord] = qubit

		print(pandas.DataFrame(layout).to_string())
		print(" =====================================================  ")


def display_qubit_mapping(qubit_mapping, layout_size):
	layout = [[0 for i in range(layout_size["width"])] for j in range(layout_size["height"])]

	for key, value in qubit_mapping.items():
		x_coord = int(value/layout_size["width"])
		z_coord = int(value%layout_size["width"])

		layout[x_coord][z_coord] = key
	
	print("===============================================")
	print(pandas.DataFrame(layout))
	print("===============================================")


def merge_qubit_layout(mapping1, mapping2, direction, layout_size):
	# function merge two blocks 

	extended_qubit_layout = {}

	if direction == "horizon":
		for key, value in mapping1.items():
			x_coord = int(key/layout_size["width"])
			z_coord = int(key%layout_size["width"])

			extended_index = x_coord * 2 * layout_size["width"] + z_coord
			extended_qubit_layout[extended_index] = value


		for key, value in mapping2.items():
			x_coord = int(key/layout_size["width"])
			z_coord = int(key%layout_size["width"])

			extended_index = x_coord * 2 * layout_size["width"] + z_coord + layout_size["width"]
			extended_qubit_layout[extended_index] = value


	elif direction == "vertical":
		extended_qubit_layout = mapping1
		for key, value in mapping2.items():
			x_coord = int(key/layout_size["width"])
			z_coord = int(key%layout_size["width"])

			index_in_extended_layout = (layout_size["height"] + x_coord) * layout_size["width"] + z_coord
			extended_qubit_layout[index_in_extended_layout] = value

	return {v: int(k) for k, v in extended_qubit_layout.items()}