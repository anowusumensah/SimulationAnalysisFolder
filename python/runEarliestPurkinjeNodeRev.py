#!/usr/bin/env python

import pandas as pd
from collections import defaultdict

PURKINJE_RANGE = range(547680, 548790)

def read_vtx_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()[2:] # skip first two lines
    return [int(line.strip()) for line in lines]

def read_purkinje_filelist(file_path):
    second_column = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                second_column.append(int(parts[1]))
    return second_column
                


def parse_pmjs_file(pmjs_path):
    """Return a mapping of cable_tag -> list of node indices (Purkinje + Myo)."""
    tag_map = defaultdict(list)
    with open(pmjs_path, 'r') as file:
        for idx, line in enumerate(file):
            tag = int(line.strip())
            if tag != -1:
                tag_map[tag].append(idx)
    return tag_map

def load_activation_data(dat_path):
    df = pd.read_csv(dat_path, sep=r"\s+", header=None, names=["node", "activation_time"])
    return df.groupby("node")["activation_time"].apply(list).to_dict()

def get_purkinje_and_coupled_myo(tag_map):
    """Returns a dict: purkinje_node -> list of coupled myocardial nodes."""
    purk_myo_map = {}
    for tag, nodes in tag_map.items():
        purk_nodes = [n for n in nodes if n in PURKINJE_RANGE]
        myo_nodes = [n for n in nodes if n not in PURKINJE_RANGE]
        for purk in purk_nodes:
            purk_myo_map[purk] = myo_nodes
    return purk_myo_map

def get_next_activation_after(times, ref_time):
    """Return the first activation time greater than ref_time from a list."""
    return min((t for t in times if t > ref_time), default=None)

def is_purkinje_activated(purk_node, myo_nodes, activation_dict, ectopic_time):
    """Check if >=50% of myo nodes have activation after ectopic time."""
    count = 0
    for node in myo_nodes:
        if node in activation_dict:
            next_time = get_next_activation_after(activation_dict[node], ectopic_time)
            if next_time and next_time > ectopic_time:
                count += 1
    return count >= len(myo_nodes) / 2

def find_earliest_activated_purkinje(purk_myo_map, activation_dict, ectopic_time, exclude_nodes):
    earliest_node = None
    earliest_time = float('inf')

    for purk_node, myo_nodes in purk_myo_map.items():
        if purk_node in exclude_nodes:
            continue

        if is_purkinje_activated(purk_node, myo_nodes, activation_dict, ectopic_time):
            next_purk_time = get_next_activation_after(activation_dict.get(purk_node, []), ectopic_time)
            if next_purk_time and next_purk_time < earliest_time:
                earliest_time = next_purk_time
                earliest_node = purk_node

    return earliest_node, earliest_time if earliest_node else (None, None)

def main(pmjs_path, dat_path, exclude_list, ectopic_time):
    tag_map = parse_pmjs_file(pmjs_path)
    activation_dict = load_activation_data(dat_path)
    purk_myo_map = get_purkinje_and_coupled_myo(tag_map)

    result_node, result_time = find_earliest_activated_purkinje(
        purk_myo_map, activation_dict, ectopic_time, set(exclude_list)
    )

    if result_node:
        print(f"Earliest activated Purkinje node: {result_node} at {result_time} ms")
    else:
        print("No valid Purkinje activation found after ectopic time.")

# Example call

def run():

    parentFolder="/home/anthonyowusu-mensah/Desktop/newSinusEctopy/EctopyForAll/"
     
    files=[
           "Ctrl-His-1000-GM-1.0-S1S2-230-Stim-50-Rpmj-45e3-mGNa-1.0-pGNa-1.0-vtx-LeftMultiple_two-savF-1-Mlump-1-parabSol-0-fit-1.0-New",
           "Ctrl-His-1000-GM-0.30-S1S2-250-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-LeftMultiple_two-savF-1-Mlump-1-fit-1.25",
           "LQT2-His-1000-GM-1.0-S1S2-250-Stim-50-Rpmj-45e3-mGNa-1.0-pGNa-1.0-vtx-LeftMultiple_two-savF-1-Mlump-1-parabSol-0-fit-1.0-New",
           "LQT2-His-1000-GM-0.30-S1S2-285-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-LeftMultiple_two-savF-1-Mlump-1-fit-1.25"
           
              ]
    fileNames = ['Ctrl_Sinus', 'Ctrl_Drug', 'LQT2_Sinus', 'LQT2_Drug']        
    ectopyTimes = [231, 251, 251, 286]
    idxVal = 0              
    for file in files:
        # Purkinje Nodes to Isolate
        fileName = fileNames[idxVal]
        electrodeType = "leftTerminalPurkList.dat"
        pathToElectrode="/media/anthonyowusu-mensah/TonySSD/AnalysisScripts/TerminalPurkinjeList/"
        #lstElectrode = read_vtx_file(pathToElectrode + electrodeType)
        lstElectrode = read_purkinje_filelist(pathToElectrode + electrodeType)
        #print(lstElectrode) 
        
        # Process Simulation
       
        
        #prmPath="/media/anthonyowusu-mensah/TonySSD/tryNew/noSTARTstate/lumpMatrix/testSimul2/mGNa_0.7_pGNa_0.95_GM_0.3/Ctrl/Ctrl-His-1000-GM-0.30-S1S2-250-Stim-50-Rpmj-45e3-  mGNa-0.7-pGNa-0.95-vtx-LeftMultiple_two-savF-1-Mlump-1-fit-1.25/"
        prmPath= parentFolder + '/' + file + '/'
        ectopyTime = ectopyTimes[idxVal]
        pmjPath=prmPath + "myopurk.pmjs"
        secActPath=prmPath + "secondAct-thresh.dat"
        
        print(f"Processing {fileName} ===========================")
        main(pmjPath, secActPath,lstElectrode,ectopyTime)
        print("")
        idxVal += 1;

if __name__ == '__main__':
    run()


