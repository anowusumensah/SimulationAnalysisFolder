#!/usr/bin/env python

"""
Written by Anthony Owusu-Mensah
current Febrruary 2024

"""

#import sys
import numpy as np
from datetime import date
from carputils import tools
import os
from scipy.signal import find_peaks 
from scipy.interpolate import CubicSpline
from writeElementFile import process_element_file
import time

#%%
EXAMPLE_DESCRIPTIVE_NAME = 'PEERP'
EXAMPLE_AUTHOR = 'Anthony <aowus003@odu.edu>'
EXAMPLE_DIR = os.path.dirname(__file__)
GUIinclude = False


def runIgbextractVoltage(igbFile, nodes, outputFileName):
    print(); print();
    print("Running Igbextract For Voltage")
    #myoNodes = list(map(str,np.arange(0,547680)))
    myoNodesStr = ",".join(nodes)
    cmd = 'igbextract'
    #cmd+=" " + '-list='+myoNodesStr 
    cmd+=" " + '--list=' +  myoNodesStr + ' -O' + " " + outputFileName
    cmd+=" " + igbFile
    print(cmd)
    os.system(cmd)


def run():

    parentFolder="/media/anthonyowusu-mensah/TonySSD/tryNew/noSTARTstate/lumpMatrix/testSimul2/sinus_mGNa_0.7_pGNa_0.95_GM_0.3/singleSinus/twoSinus/"
     
    files=[ "Ctrl-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-0-Mlump-1-fit-1.25",
            "Ctrl-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.25-NF-2",
            "LQT2-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-0-Mlump-1-fit-1.25",
            "LQT2-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.25-NF-2"
              ]
               
    fileNames = ["Ctrl_Drug_savF_0", "Ctrl_Drug_savF_1", "LQT2_Drug_savF_0", "LQT2_Drug_savF_1"] 
    nodes = [["548789", "349529"], ["548789", "349529"], ["548789", "349529"], ["548789", "349529"]]
    idxSel = 0
    for file in files: 
        outputFileName = fileNames[idxSel]
        igbFileName = os.path.join(parentFolder,file, file + ".igb")
        checkPath = os.path.join(os.getcwd(),outputFileName) 
        if not os.path.isfile(checkPath):   
            runIgbextractVoltage(igbFileName, nodes[idxSel], outputFileName)
        idxSel += 1

if __name__ == '__main__':
    run()

