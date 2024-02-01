# -*-coding:utf-8-*-

import os
import error
import simplejson as json
import numpy as np
import itertools
import collections
from pprint import pprint

def generate_regular_qchip_architecture(parent_dir, layout_size, **kwargs):
    '''
        function to make a file of qubit architecture
    '''

    qubit_connectivity = collections.defaultdict(list)
    
    if "architecture" in kwargs:
        architecture = kwargs["architecture"]
    else:
        architecture = 2

    width = layout_size.get("width")
    height = layout_size.get("height")
    length = layout_size.get("length")

    if length is None: length = 1
    else: length = int(length)

    qubits = width * height * length

    # 0 : all-to-all connections between qubits
    if architecture == 0:
        for idx in range(qubits):
            qubit_connectivity[idx] = list(range(qubits))
            qubit_connectivity[idx].remove(idx)

    elif architecture == 23:
        # 2-dimensional triangular lattice
        # 1. 먼저 2d lattice
        for idx in itertools.product(range(height), range(width)):
            cell_idx = idx[0] * width + idx[1]

            if not idx[0]:
                if height > 1:
                    qubit_connectivity[cell_idx].append(cell_idx + width)

            elif idx[0] < height - 1:
                qubit_connectivity[cell_idx].extend([cell_idx-width, cell_idx+width])

            elif idx[0] == height - 1:
                qubit_connectivity[cell_idx].append(cell_idx-width)

            if not idx[1]:
                if width > 1:
                    qubit_connectivity[cell_idx].append(cell_idx+1)

            elif idx[1] < width - 1:
                qubit_connectivity[cell_idx].extend([cell_idx-1, cell_idx+1])

            elif idx[1] == width - 1:
                qubit_connectivity[cell_idx].append(cell_idx-1)

        # 2. 오른쪽 사선으로..
        for row_idx in range(height-1):
            for col_idx in range(width-1):
                source_idx = row_idx * width + col_idx
                target_idx = (row_idx+1) * width + (col_idx+1)

                qubit_connectivity[source_idx].append(target_idx)
                qubit_connectivity[target_idx].append(source_idx)

    elif architecture == 21:
        # 2-dimensional complex lattice 
        # 1. 먼저 2d lattice
        for idx in itertools.product(range(height), range(width)):
            cell_idx = idx[0] * width + idx[1]

            if not idx[0]:
                if height > 1:
                    qubit_connectivity[cell_idx].append(cell_idx + width)

            elif idx[0] < height - 1:
                qubit_connectivity[cell_idx].extend([cell_idx-width, cell_idx+width])

            elif idx[0] == height - 1:
                qubit_connectivity[cell_idx].append(cell_idx-width)

            if not idx[1]:
                if width > 1:
                    qubit_connectivity[cell_idx].append(cell_idx+1)

            elif idx[1] < width - 1:
                qubit_connectivity[cell_idx].extend([cell_idx-1, cell_idx+1])

            elif idx[1] == width - 1:
                qubit_connectivity[cell_idx].append(cell_idx-1)
        
        # 2. 왼쪽 사선으로...
        # 3. 오른쪽 사선으로..

        for col_idx in range(1, width, 2):
            row_idx = 0
            list_connection = []
            source_idx = row_idx*width + col_idx
            
            while True:
                if row_idx == height-1 or col_idx == 0:
                    break

                row_idx+=1
                col_idx-=1

                target_idx = row_idx*width + col_idx
                qubit_connectivity[source_idx].append(target_idx)
                qubit_connectivity[target_idx].append(source_idx)

                source_idx = target_idx

        # missing point
        if width%2 == 0: starting_idx = 2
        else: starting_idx = 1

        for row_idx in range(starting_idx, height, 2):
            col_idx = width-1

            list_connection = []
            source_idx = row_idx*width + col_idx

            while True:
                if row_idx == height-1 or col_idx == 0:
                    break
                row_idx+=1
                col_idx-=1

                target_idx = row_idx *width + col_idx
                qubit_connectivity[source_idx].append(target_idx)
                qubit_connectivity[target_idx].append(source_idx)
                source_idx = target_idx

        for col_idx in range(1, width, 2):
            row_idx = 0
            list_connection = []
            source_idx = row_idx*width + col_idx            
            
            while True:
                if row_idx == height-1 or col_idx == width-1:
                    break

                row_idx+=1
                col_idx+=1
                target_idx = row_idx * width + col_idx

                qubit_connectivity[source_idx].append(target_idx)
                qubit_connectivity[target_idx].append(source_idx)

                source_idx = target_idx
        
        # # missing point
        for row_idx in range(1, height, 2):
            col_idx = 0
            list_connection = []
            source_idx = row_idx*width + col_idx

            while True:
                if row_idx == height-1 or col_idx == width-1:
                    break
             
                row_idx+=1
                col_idx+=1
                target_idx = row_idx*width + col_idx
                qubit_connectivity[source_idx].append(target_idx)
                qubit_connectivity[target_idx].append(source_idx)
                source_idx = target_idx

    elif architecture == 2 :
        # 2-dimensional rectangular lattice
        for idx in itertools.product(range(height), range(width)):
            cell_idx = idx[0]*width + idx[1]
            
            if not idx[0]:
                if height > 1:
                    qubit_connectivity[cell_idx].append(cell_idx + width)

            elif idx[0] < height - 1:
                qubit_connectivity[cell_idx].extend([cell_idx-width, cell_idx+width])

            elif idx[0] == height - 1:
                qubit_connectivity[cell_idx].append(cell_idx-width)

            if not idx[1]:
                if width > 1:
                    qubit_connectivity[cell_idx].append(cell_idx+1)

            elif idx[1] < width - 1:
                qubit_connectivity[cell_idx].extend([cell_idx-1, cell_idx+1])

            elif idx[1] == width-1:
                qubit_connectivity[cell_idx].append(cell_idx-1)

    elif architecture == 3:
        # 3-dimensional cube lattice
        for idx in itertools.product(range(height), range(length), range(width)):
            cell_idx = idx[0] * length * width + idx[1] * width + idx[2]

            # coord : x, z, y
            list_neigbor_coord = []

            if idx[0] > 0:
                list_neigbor_coord.append((idx[0]-1, idx[1], idx[2]))
            if idx[0] < height - 1:
                list_neigbor_coord.append((idx[0]+1, idx[1], idx[2]))

            if idx[1] > 0:
                list_neigbor_coord.append((idx[0], idx[1]-1, idx[2]))
            if idx[1] < length - 1:
                list_neigbor_coord.append((idx[0], idx[1]+1, idx[2]))

            if idx[2] > 0:
                list_neigbor_coord.append((idx[0], idx[1], idx[2]-1))
            if idx[2] < width - 1:
                list_neigbor_coord.append((idx[0], idx[1], idx[2]+1))

            list_neighbor_idx = [neighbor[0]*length*width + neighbor[1]*width + neighbor[2] \
                                    for neighbor in list_neigbor_coord]

            qubit_connectivity[cell_idx].extend(list_neighbor_idx)

    else:
        raise Exception("Given architecture is wrong : {}".format(architecture))


    if length > 1:
        file_device = "".join(["file_qchip_{}x{}x{}.json".format(height, length, width)])
    else:
        file_device = "".join(["file_qchip_{}x{}.json".format(height, width)])

    qchip_architecture = {"qubit_connectivity": qubit_connectivity, 
                          "device_name": file_device,
                          "dimension": layout_size,
                          "architecture": architecture}

    full_path_device = os.path.join(parent_dir, file_device)


    with open(full_path_device, "w") as f:
        json.dump(qchip_architecture, f, sort_keys=True, indent=4, separators=(',', ':'))

    return {"result_file": full_path_device, "qchip_architecture": qchip_architecture}


if __name__ == "__main__":
    ret = generate_regular_qchip_architecture(".", {"width": 4, "height": 4}, architecture=23)

    pprint(ret)