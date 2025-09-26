# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:08:26 2024

@author: Tony
"""
import os
from carputils import format
import numpy as np
from carputils import carpio
from carputils import settings
from carputils import tools



    
def generateMeshCmd(args):
    folder_path = 'meshes'
    meshname ='Tune_' + str(args.sres) 
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created.")
    else:
        print(f"Folder '{folder_path}' already exists.")
    
    um2cm = 1/10000
    meshDim = [args.slen, args.sres * um2cm, args.sres * um2cm] #Mesh dimensions
    meshDimStr = [str(element) for element in meshDim]
    sizes = " ".join(meshDimStr)
    
    res = [args.sres, args.sres, args.sres] #Mesh resolution
    res2Str = [str(val) for val in res]
    resVals = " ".join(res2Str)

    meshfolder = os.path.join(folder_path,meshname)
    if not os.path.exists(meshfolder + '.pts'):
        print('Generating mesh: {}'.format(meshname))
        cmd = 'mesher ' 
        cmd += "-size '{" +  sizes + "}' "  
        cmd += "-resolution '{" + resVals + "}' "
        cmd += "-Elem3D " + str(args.Elem3D) + " "
        cmd += '-mesh ' + meshfolder
        print(cmd)
        os.system(cmd)
    else:
        print(meshfolder + " exists")
    return meshfolder 
    
    

# Conduction velocity observation sites
def CV_obs_sites(args, rel_x):
    """
    Extract spatial observation points for CV measurement
    """
    
    # read ub strand
    folder_path = 'meshes'
    meshname ='Tune_' + str(args.sres) 
    meshfiles = os.path.join(folder_path,meshname) + '.pts'
    pts = np.loadtxt(meshfiles, skiprows=1, usecols=(0,1,2))
   
    # figure out nodes which fall into transitions
    # between quarter 1-2 and 3-4
    lo, up = np.min(pts, axis=0), np.max(pts, axis=0)
    slen   = up[0]-lo[0]
    dx = float(args.sres)
    x  = []
    for rx in rel_x:
        x += [ lo[0] + slen*rx ]
 
    nq = []
    for xi in x: 
        # find points around xi +/- dx/2
        o_pts = np.where(np.logical_and(pts[:,0] >= (xi-dx/2), pts[:,0] < (xi+dx/2)))[0]
        nq.append(o_pts)

    # mean x coordinates of selected points
    xm = []
    for i,n in enumerate(nq):
        xs = pts[n,0]
        xm.append(np.mean(xs))

    return xm, nq

def measureCV(args,actsFile):
    """
    Measure conduction velocity along a strand

    As input the given mesh of a strand and the recorded activation times are used.
    Arrival times are computed at the interfaces of the quarters Q1/Q2 and Q3/Q4.
    """

    # observation sites x location relative to strand length
    rel_x = [0., 0.25, 0.5, 0.75, 1.0]
    obs_x, obs_nodes = CV_obs_sites(args, rel_x)
    

    # load local activation times
    if actsFile.endswith('.dat'):
        lats = np.loadtxt(actsFile)
    elif actsFile.endswith('.igb'):
        igb_file = carpio.igb.IGBFile(actsFile, 'r')
        header = igb_file.header()
        lats = igb_file.data().reshape((header['t'], header['x'])).T
    else:
        raise IOError('Error, unsupported file format, .dat or .igb file expected!')

    # arrival times at observation points
    x0   = obs_x[0]
    at0  = np.mean(lats[obs_nodes[0]])   # lat at stimulus nodes
    x12  = obs_x[1]
    at12 = np.mean(lats[obs_nodes[1]])   # lat at 25% into the strand
    x34  = obs_x[3]
    at34 = np.mean(lats[obs_nodes[3]])   # lat at 75% into the strand
    
    
# =============================================================================
#     print('x34 = {}, x12 = {}'.format(x34, x12))
#     print()
#     print('at34 = {}, at12 = {}'.format(at34, at12))
# =============================================================================


    # make sure we had capture and tissue activated
    if (at12 == -1) or (at34 == -1):
        if (at12 > 0.):
            if not settings.cli.silent:
                print('Error: Duration of simulation too short')
                #print('stimDur:', sdur)
                print('Proximal tissue activated, but not distal tissue')
                print('Estimating velocity based on activation in initial segment of the cable.')

            # estimate from initial segmentation
            CV = ((x12 - x0)/(at12 - at0))* 1e-3

        else:
            if not settings.cli.silent:
                print('No capture or decremental conduction, proximal quarter of tissue not fully activated')

            if args.conv:
                CV = 0.05
                if not settings.cli.silent:
                    print('Assuming slow conduction velocity of 0.05 m/s in next iteration')

            else:
                CV = 0.0

    else:
        CV = ((x34 - x12)/(at34 - at12))* 1e-3

    return CV  # Conduction velocity in m/s



def configStimLHFace(args):
     # read ub strand
    folder_path = 'meshes'
    meshname ='Tune_' + str(args.sres) 
    meshfiles = os.path.join(folder_path,meshname) + '.pts'
    pts = np.loadtxt(meshfiles, skiprows=1, usecols=(0,1,2))
    x_min, y_min, z_min = pts.min(axis=0)
    stimS = 250
    #sdur = (150/stimS) * 5
    #sdur = (150/stimS) * 30
    sdur = 10
    #print(f"sdur: {sdur}")
    idx  = 0
    stim = ['-stimulus[%d].duration'%(idx),  sdur,
            '-stimulus[%d].stimtype'%(idx),      0,
            '-stimulus[%d].strength'%(idx), stimS,
            '-stimulus[%d].x0'%(idx), x_min,
            '-stimulus[%d].y0'%(idx), y_min,
            '-stimulus[%d].z0'%(idx), z_min,
            '-stimulus[%d].xd'%(idx), args.sres,
            '-stimulus[%d].yd'%(idx), args.sres,
            '-stimulus[%d].zd'%(idx), args.sres]
    return stim 


    
    

def configLATs():
    """
    Set up detection of local activation times
    """
    # default, we detect only the first lat
    lidx = 0
    latstr = '-lats[{}]'.format(lidx)
    lats_first = [latstr + '.measurand',   0,
                  latstr + '.method',      1,
                  latstr + '.all',         0,
                  latstr + '.threshold', -40.,
                  latstr + '.mode',        0,
                  latstr + '.ID',        'lats']

    lidx += 1
    latstr = '-lats[{}]'.format(lidx)
    repol      = [latstr + '.measurand',   0,
                  latstr + '.method',      1,
                  latstr + '.all',         0,
                  latstr + '.threshold', -70.,
                  latstr + '.mode',        1,
                  latstr + '.ID',        'repol_time']

    

    lats = ['-num_LATs', lidx + 1 ] + lats_first + repol 

    return lats


def configNumerics(args):
    """
    Define numerical settings for tuning experiment
    """
    num_set  = [ '-dt',                 args.dt,
                 '-parab_solve',        args.ts,
                 '-cg_tol_parab',       args.stol,
                 '-mass_lumping',       int(args.lumping),
                 '-operator_splitting', int(args.opsplt) ]
    return num_set 



def buildCmd(args,parameterFile):
    #meshname = geomStrand(args, job)
    meshname = generateMeshCmd(args)
    lats = configLATs() # Activation
    runSingleCell(args) # runSinglecell
    ion_set = configCells(args) # setup ionic parameters
    stim = configStimLHFace(args) # setup stimulus
    num_set = configNumerics(args) # configure numerics
    cmd = tools.carp_cmd(parameterFile)
    cmd += ['-meshname', meshname] 
    cmd += lats 
    cmd += ion_set
    cmd += stim
    cmd += num_set
    # Simulation duration
    cableLen_mm = args.slen*10
    delay       = cableLen_mm / args.vel
    tend = delay*20
    cmd += ['-tend', tend ]
    cmd += ['-timedt', int(tend)]
    return cmd


def configCells(args):
    """
    setup ionic model we want to use for CV tuning
    """
    
        
    
    
    ion_set = ['-imp_region[0].cellSurfVolRatio', args.beta,
               '-imp_region[0].im',               args.model,
               '-imp_region[0].im_param',         args.mpar,
               '-imp_region[0].plugins',          args.plugin,
               '-imp_region[0].plug_param',       args.ppar ]
    
    if args.svinit != '' and args.setInit:
        # check whether file exists
        if not os.path.isfile(args.svinit):
            print('Given {} does not exist! Aborting...'.format(args.svinit))
            exit(1)
        ion_set += ['-imp_region[0].im_sv_init', args.svinit]
        
        
    if args.ext:
        modelLoc = args.extDir + '/' + args.model + '.so'  
        extPar = ['-num_external_imp', 1, 
                  '-external_imp[0]', modelLoc]
        ion_set += extPar
                        
    return ion_set 


def configConductivities(args):
    """
    Define myocardial conductivities
    """
    gi = args.gi
    ge = args.ge

    cset = ['-gregion[0].g_il', gi,
            '-gregion[0].g_el', ge,
            '-bidm_eqv_mono', args.beqm ]

    return cset

# iterate to converge conduction velocity
def iterateCV(args, job,parameterFile):
      # build command line
      cmd = buildCmd(args,parameterFile)
      # Run simulation
      err_flg = False
      i     =  0
      v_new = -1  # do at least one run
      vel_ = -1
      gi_ = args.gi
      ge_ = args.ge
    
      # start tuning with chosen start values
      if not settings.cli.silent:
          print('Start iterating')
          print('-'*format.terminal_width() + '\n')
    
      if not args.conv:
          args.max_it = 1 
      
      #print('Max iterations=', args.max_it)
    
    
      while abs(args.vel - vel_) > args.ctol and i < args.max_it:
    # rerun simulation
          g = configConductivities(args)
          if args.conv:
              ItJobID = '{}_It_{}'.format(job.ID, i)
          else:
              ItJobID = job.ID
          job.carp(cmd + g + ['-simID', ItJobID ])
        
          # measure macroscopic conduction velocity
          lats_file = os.path.join(ItJobID, 'init_acts_lats-thresh.dat')
          if not os.path.exists(lats_file):
              lats_file = os.path.join(ItJobID, 'init_acts_lats-thresh.igb')
          vel_ = measureCV(args,lats_file)
          #print('computed velocity = ', vel_)
        
        
          # adjust conductivities only if converge is True
          if args.max_it > 1 and  vel_ > 0:
              #print('Azonto')
              alpha  = (args.vel/vel_)**2.
              if args.beqm:
                gi_ = alpha * args.gi
                ge_ = alpha * args.ge
                hmean = (gi_* ge_)/ (gi_ +  ge_)
           
              else:
                  # beqm is off
                  gi_ = alpha * args.gi
                  ge_ = args.ge
                  hmean = gi_
          else:
                if args.beqm:                   
                    gi_ = args.gi
                    ge_ = args.ge
                    hmean = (gi_* ge_)/ (gi_ +  ge_)
                else:
                    # beqm is off
                    gi_ = args.gi
                    ge_ = args.ge
                    hmean = gi_            
          args.gi = gi_
          args.ge = ge_
          last_it = i
          i = i+1
          
      if not settings.cli.silent and args.ar == 0:
          print('\n\n')
          print('Tuning results')
          print('-'*format.terminal_width() + '\n')
          print('Conduction velocity: {:.4f} m/s [gi={:.4f}, ge={:.4f}, gm={:.4f}]'.format(vel_, gi_, ge_, hmean))
          
      if args.conv and i==args.max_it and i!=1:
          err_flg = True
          print('CV tuning not converged with given max_it number of iterations!')
      #return gi_, ge_, vel_, hmean, i, err_flg
      return gi_, ge_, hmean,vel_ 
      
          
              
              
def runSingleCell(args):
    init_file = '{}_initFile_bcl_{}.sv'.format(args.model, args.His) 
# =============================================================================
#     if os.path.exists(init_file):
#         os.system('rm ' +  init_file)  
# =============================================================================
    
    if args.ext:
        passImp = ['--load-module '] # Determines how ionic model will be passed
        modelLoc = args.extDir + '/' + args.model + '.so'
    else:
        passImp = ['--imp ']
        modelLoc = args.model
        
    
        
    # command   
    if not os.path.exists(init_file) and args.setInit:
        duration = args.His * 100
        cmd = 'bench '
        cmd +=" " + passImp[0] +  " " + modelLoc
        
        if  args.mpar:
            cmd +=" " + '--imp-par' + " " + args.mpar
            init_file = args.mpar + '{}_initFile_bcl_{}.sv'.format(args.model, args.His)
        cmd +=" " + '--bcl' + " "  + str(args.His)
        cmd +=" " + '--dt-out' + " "  + str(50)
        cmd +=" " + '--stim-curr' + " " + str(60.0)
        cmd +=" " + '--stim-dur' + " " + str(1.0)
        cmd +=" " + '--numstim' +  " " + str(100)
        cmd +=" " + '--duration' + " " + str(duration)
        cmd +=" " + '--stim-start' + " " + str(1.0)
        cmd +=" " + '--dt' + " " + str(0.01)
        cmd +=" " + '--no-trace ' 
        cmd +=" " + '-S' + " "    + str(duration)
        cmd +=" " + '-F ' +  init_file
        #cmd +=" " + '--no-trace' + " " + 'on'
      
        print(cmd)
        os.system(cmd)
        args.svinit =  init_file
    else:
        args.svinit =  init_file
         
    
    
    
        

            
             





