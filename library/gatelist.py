
# -*-coding:utf-8-*-
# from __future__ import unicode_literals
# version 0.5

# gate list of standard qasm
pauli_x = 'X'
pauli_z = 'Z'
pauli_y = 'Y'
pauli_i = 'I'

h = 'H'
t = 'T'
tdag = 'Tdag'

s = 'S'
sdag = 'Sdag'

rz = 'Rz'
rx = 'Rx'
ry = 'Ry'
u = "U"

sx = 'SX'
phase = 'P'

cnot = 'CNOT'
cphase = 'CP'

cx = 'CX'
cz = 'CZ'
cphase = 'CP'

rzz = 'Rzz'
swap = 'SWAP'

lcnot = 'LCNOT'
lswap = 'LSWAP'

cnotv = 'CNOTV'
cnoth = 'CNOTH'

swapv = 'SWAPV'
swaph = 'SWAPH'

prepz = 'PrepZ'
measz = 'MeasZ'

prepx = 'PrepX'
measx = 'MeasX'

move = 'MOVE'
move_back = 'MoveBack'

barrier_all = "Barrier-All"
barrier = "Barrier"

# list of gates 
list_prepare = [prepz, prepx]
list_measure = [measz, measx]
list_pauli_gates = [pauli_x, pauli_y, pauli_z, pauli_i]
list_1q_rotations = [rx, ry, rz, u, phase]
list_2q_rotations = [rzz, cphase]

list_1q_gates = [pauli_x, pauli_y, pauli_z, pauli_i, h, t, tdag, s, sdag,\
					rx, ry, rz, u, phase] + list_measure + list_prepare
list_2q_gates = [cnot, cphase, cz, cx, rzz, swap]
list_barrier = [barrier, barrier_all]
list_move = [move, move_back]
list_register = ["Qubit", "Cbit"]

ibm_quantum_register = 'qreg'
ibm_classical_register = 'creg'

ibm_cnot = 'cx'
ibm_cz = 'cz'
ibm_cphase = 'cp'
ibm_cu1 = 'cu1'
ibm_rzz = 'rzz'
ibm_hadamard = 'h'
ibm_id = 'id'
ibm_x = 'x'
ibm_z = 'z'
ibm_y = 'y'
ibm_s = 's'
ibm_sdag = 'sdg'
ibm_t = 't'
ibm_tdag = 'tdg'
ibm_idle = 'id'
ibm_rx = 'rx'
ibm_ry = 'ry'
ibm_rz = 'rz'
ibm_u1 = 'u1'
ibm_u2 = 'u2'
ibm_u3 = 'u3'
ibm_u = 'u'
ibm_sx = 'sx'
ibm_h = 'h'
ibm_measure = 'measure'
ibm_prepare = 'prepare'
ibm_phase = "p"
ibm_swap = "swap"
ibm_barrier = "barrier"
ibm_reset = "reset"

gate_standard2open = {
    pauli_x: ibm_x,
    pauli_z: ibm_z,
    pauli_y: ibm_y,
    pauli_i: ibm_id,
    h: ibm_hadamard,
    t: ibm_t,
    tdag: ibm_tdag,
    s: ibm_s,
    sdag: ibm_sdag,
    rz: ibm_rz,
    rx: ibm_rx,
    ry: ibm_ry,
    rzz : ibm_rzz,
    sx: ibm_sx,
    cnot: ibm_cnot,
    cz : ibm_cz,
    measz: ibm_measure,
    phase: ibm_phase,
    u: ibm_u,
    prepz: ibm_prepare,
    cphase: ibm_cphase,
    swap: ibm_swap
}


gate_open2standard = {
    ibm_hadamard: h,
    ibm_x: pauli_x,
    ibm_z: pauli_z,
    ibm_y: pauli_y,
    ibm_id: pauli_i,
    ibm_t: t,
    ibm_tdag: tdag,
    ibm_s: s,
    ibm_sdag: sdag,
    ibm_rz: rz,
    ibm_rx: rx,
    ibm_ry: ry,
    ibm_sx: sx,
    ibm_u: u,
    ibm_u1: u,
    ibm_u2: u,
    ibm_u3: u,
    ibm_cnot: cnot,
    ibm_cphase : cphase,
    ibm_cz : cz,
    ibm_rzz : rzz,
    ibm_measure: measz,
    ibm_prepare: prepz,
    ibm_phase: phase,
    ibm_swap: swap
}
