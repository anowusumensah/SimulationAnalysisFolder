#!/usr/bin/env python

"""
Written by Anthony Owusu-Mensah
current Febrruary 2024

"""
#%%
import os
EXAMPLE_DESCRIPTIVE_NAME = 'PEERP'
EXAMPLE_AUTHOR = 'Anthony <aowus003@odu.edu>'
EXAMPLE_DIR = os.path.dirname(__file__)
GUIinclude = False

#import sys
import numpy as np
from datetime import date
from carputils import tools
#import os
from scipy.signal import find_peaks 
from scipy.interpolate import CubicSpline
from writeElementFile import process_element_file
import time

def parser():
    parser = tools.standard_parser()
    group  = parser.add_argument_group('experiment specific options')
    group.add_argument('--elemfile',type=str,default='Elemfile', help='element file')
    group.add_argument('--expType',type=str,default='LQT2_APDcheck', help='experimentType')
    group.add_argument('--startTime',type=float,default=2000, help='start Time') 
    group.add_argument('--EndTime',type=float,default=2200, help='start Time')
    group.add_argument('--tag',type=int,default=666, help='region Tag')
    group.add_argument('--IgbFileName',type=str,default='vm.igb', help='igbFile')
    group.add_argument('--scaleDvdt',type=float,default= 0.5, help='scaling for maxDvdt')
    group.add_argument('--threshold',type=float,default= -65, help='threshold voltage')
    
    return parser

def jobID(args):
    """
    Generate name of top level output directory.
    """
    today = date.today()
    out_DIR = '{}-startTime-{}-endTime-{}'.format(today.isoformat(),args.startTime,args.EndTime)
    return out_DIR 

def getIndixes(arrayDvdt, startTime, endTime):
    indexes = np.where((arrayDvdt >= startTime) & (arrayDvdt <= endTime))
    return indexes[0] 
    

def runIgbops(igbFile, dvdtFileName):
    #myoNodes = list(map(str,np.arange(0,547680)))
    #myoNodesStr = ",".join(myoNodes)
    print(); print();
    print("Igbopbs")
    cmd = '/opt/carpentry/bin/igbops'
    cmd+=" " + '--op=diff_c' # Compute the forward difference 
    cmd+=" " + '--output-file='+ dvdtFileName
    cmd+=" " + igbFile
    print(cmd)
    os.system(cmd)
    
def runIgbextract(dvdtFileName):
    print(); print();
    print("Running Igbextract For Activation Time")
    #myoNodes = list(map(str,np.arange(0,547680)))
    #myoNodesStr = ",".join(myoNodes)
    cmd = 'igbextract'
    #cmd+=" " + '-list='+myoNodesStr 
    cmd+=" " + '-O' + " " + 'out'
    cmd+=" " + dvdtFileName + '.igb'
    print(cmd)
    os.system(cmd)

def runIgbextractVoltage(igbFile):
    print(); print();
    print("Running Igbextract For Voltage")
    #myoNodes = list(map(str,np.arange(0,547680)))
    #myoNodesStr = ",".join(myoNodes)
    cmd = 'igbextract'
    #cmd+=" " + '-list='+myoNodesStr 
    cmd+=" " + '-O' + " " + 'apd.dat'
    cmd+=" " + igbFile
    print(cmd)
    os.system(cmd)

    
# Generate Tags
def tagregopt( reg, field, val ) :
    return ['-tagreg['+str(reg)+'].'+field, val ]
    
def returnTagRegion(args):
    reg = 0
    dyn_regH = ['-numtagreg', 1]
    dyn_regH += tagregopt(reg, 'type',        4)
    dyn_regH.extend(tagregopt(reg, 'elemfile',  args.elemfile))
    dyn_regH.extend( tagregopt(reg, 'tag',     args.tag))
    return dyn_regH        

# Read Elementfile
def readElemenFile(dotElemFile):
    data = np.loadtxt(dotElemFile,dtype=int,skiprows=1,usecols=(1,2,3,4),)
    return data

def intersection(lst1, lst2):
    return len(list(set(lst1) & set(lst2))), list(set(lst1) & set(lst2))

def detecPeaks(actualArray, numIterations):
    avgPeakofEachNode = []
    for i in range(numIterations):
        numpyArray = actualArray[:,i]
        peaks, _ = find_peaks(numpyArray, height=0)
        actualPeaks = np.mean(numpyArray[peaks])
        avgPeakofEachNode.append(actualPeaks)
    return avgPeakofEachNode, np.mean(np.array(avgPeakofEachNode))
        
            
def find_peaks_columns(data, val):
    #peaks_indices = []
    #for i in range(data.shape[0]):  # Iterate over columns
    column = data
    peaks, _ = find_peaks(column,height=val)
    #peaks_indices.append(peaks)
    return peaks         

@tools.carpexample(parser, jobID)
def run(args, job):
    igbFileName = args.IgbFileName
    if not os.path.isfile(os.path.join(os.getcwd(), 'apd.dat')):
        runIgbextractVoltage(igbFileName)
    
    
    # Read ElemenFile(colums 2 - 5 of .elem file)
    dotElemFile = os.path.join(os.getcwd(), 'myopurk.elem')    
    #colTwoToFive = readElemenFile(dotElemFile)
    #numRows = 3073529
    # Times to turn vertices on
    
    # Read Apd file
    apdData = np.loadtxt(os.path.join(os.getcwd(), 'apd.dat'))
    timeAxis=apdData[:,0]
    apDataMyo = apdData[:,1:547681] #
    threshVoltage = args.threshold
    #T1 = [1285, 1405, 1510]
    T1 = [4065, 4110, 4245]
    lenT1 = len(T1)
    for i in range(lenT1):
        args.startTime = T1[i]
        activatedMyoNode = [] # Activated Nodes
        for j in range(apDataMyo.shape[1]):
            apdDataAtNode = apDataMyo[:,j] # Action Potential data
            interPolator = CubicSpline(timeAxis, apdDataAtNode)
            voltageAtTime = interPolator(args.startTime) 
            if voltageAtTime >= threshVoltage:
                activatedMyoNode.append(j)                                    
        """
        #Writing activated Nodes to file 
        elemFile = args.expType + '_' + str(args.startTime) + '_threshold_' + str(threshVoltage)
        file=open(elemFile + '.elem','w')
        file.write(str(numRows) + '\n')
        for w in range(numRows):
            connectedNodes = colTwoToFive[w]
            toList = list(connectedNodes);
            intersect,_ = intersection(toList,activatedMyoNode)
            tag_value = args.tag if intersect else 100 
            file.write(f"Tt {toList[0]} {toList[1]} {toList[2]} {toList[3]} {tag_value}\n") 
        file.close()
        """
        start_time = time.time()
        print("Writing File ",  T1[i])
        print("For " + str(T1[i]) + ", the number of nodes is " + str(len(activatedMyoNode)))
        #elemFile = args.expType + '_' + str(args.startTime) + '_threshold_' + str(threshVoltage)
        elemFile = "ChatGPT_"+ str(args.startTime)
        input_file = dotElemFile
        output_file = elemFile + '.elem'
        element_type = "Tt"
        target_nodes = activatedMyoNode
        tag_tuple = (100, args.tag)
        process_element_file(input_file,
                                 output_file,
                                 element_type,
                                 target_nodes,
                                 tag_tuple)
        end_time = time.time()
        print(f"Execution Time: {end_time - start_time:.4f} secs")
                                                            
if __name__ == '__main__':
    run()










