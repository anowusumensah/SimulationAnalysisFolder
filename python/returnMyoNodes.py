#!/usr/bin/env python

from collections import defaultdict
"""
Written by Anthony Owusu-Mensah
current Febrruary 2024

"""

PURKINJE_RANGE = range(547680, 548790)

def parse_pmjs_file(pmjs_path):
    """Return a mapping of cable_tag -> list of node indices (Purkinje + Myo)."""
    tag_map = defaultdict(list)
    with open(pmjs_path, 'r') as file:
        for idx, line in enumerate(file):
            tag = int(line.strip())
            if tag != -1:
                tag_map[tag].append(idx)
    return tag_map

def get_purkinje_and_coupled_myo(tag_map):
    """Returns a dict: purkinje_node -> list of coupled myocardial nodes."""
    purk_myo_map = {}
    for tag, nodes in tag_map.items():
        purk_nodes = [n for n in nodes if n in PURKINJE_RANGE]
        myo_nodes = [n for n in nodes if n not in PURKINJE_RANGE]
        for purk in purk_nodes:
            purk_myo_map[purk] = myo_nodes
    return purk_myo_map

def find_key_by_value(item, data_dict):
    for key in data_dict:
        value = data_dict[key]
        if (item in value):
            return key 
    
       
             
           
def run():

    purkinjeNodeNumber = 548500
    myoInterest = 437749
    pmjs_path = "/media/anthonyowusu-mensah/TonySSD/AnalysisScripts/myopurk.pmjs"
    tag_map = parse_pmjs_file(pmjs_path)
    purk_myo_map = get_purkinje_and_coupled_myo(tag_map)
    
    
    if purkinjeNodeNumber in purk_myo_map:
        for item in purk_myo_map[purkinjeNodeNumber]:
            print (item)
            #if (myoInterest == item):
            #    print (f"{myoInterest} found")  
    else:
        print(f"Terminal Purkinje {purkinjeNodeNumber} not found")
    
     
    key = find_key_by_value(myoInterest, purk_myo_map)
    find_key_by_value(myoInterest, purk_myo_map)
    print(f"Terminal Purkinje node for {myoInterest} is {key}")
    
if __name__ == '__main__':
    run()
