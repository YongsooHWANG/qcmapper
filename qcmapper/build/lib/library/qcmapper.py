# -*-coding:utf-8-*-
import os
import sys
import simplejson as json
import multiprocessing
import library.formatconversion as formatconversion


def map_circuit(path_qasm, path_qchip, **kwargs):
    '''
        synthesize a quantum computing circuit written in standard qasm
    '''
    # 합성 알고리즘에 따라, 해당 함수 목록
    mapper_function = {"SABRE": call_sabre_mapper,
                       "dijkstra": call_dijkstra_mapper}

    synthesis_option = kwargs.get("option")
        
    # 선택한 mapper에 따라서 해당 함수 호출
    # 선택한 mapper가 없으면, dijkstra
    mapper = synthesis_option.get("mapper")
    if mapper is None:
        mapper = "dijkstra"

    transformed_code, qubit_assoc, cbit_assoc = formatconversion.transform_to_standardqasm(path_qasm)
    ret = mapper_function[mapper](transformed_code, path_qchip, option=synthesis_option)
    ret = formatconversion.transform_to_openqasm(ret, qubit_association=qubit_assoc,
                                                      cbit_association=cbit_assoc)
        
    return ret


def call_sabre_mapper(qasm, qchip, **kwargs):
    '''
        function to call SABRE mapper
    '''
    arguments = kwargs.get("option")
    
    calibration_aware = arguments.get("calibration")
    
    if calibration_aware is None:
        calibration_aware = False

    if calibration_aware:
        import library.SABRE_calibration_aware_v5 as SABRE
    else:
        import library.SABRE_v6 as SABRE

    return SABRE.manage_synthesize(qasm, qchip, 
                        synthesis_option=arguments, 
                        qubit_table=arguments.get("qubit_table"))


def call_dijkstra_mapper(qasm, qchip, **kwargs):
    '''
        function to call dijkstra mapper 
    '''
    import library.dijkstramapping_v2 as dijkstramapping
    arguments = kwargs.get("option")
    
    return dijkstramapping.synthesize_dijkstra_manner(qasm, qchip, 
                                                      synthesis_option=arguments)



if __name__ == "__main__":
    # multiprocessing.freeze_support()
    from pprint import pprint

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
                    "calibration" : False,
                    "iteration": 10, "cost": "lap", "optimal_criterion": "number_gates",
                    "initial_mapping_option": "random", "mapper": "SABRE"}

    for algorithm in list_kisti_algorithms:
        path_qasm = os.path.join("../DB-Assembly/KISTI_sample", algorithm)
        # path_qasm = os.path.join("../DB-Assembly", algorithm)

        ret = map_circuit(path_qasm, path_qchip, option=synthesis_option)
        
        file_basename = os.path.basename(path_qasm)
        algorithm_name = os.path.splitext(file_basename)[0]
        file_syscode = "{}.syscode".format(algorithm_name)

        with open(file_syscode, "w") as outfile:
            json.dump(ret, outfile, sort_keys=True, indent=4, separators=(',', ':'))
    