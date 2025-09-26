#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

To run the experiments of this tutorial do

.. code-block:: bash

    cd ${TUTORIALS}/02_EP_tissue/03A__study_prep_tuneCV


Concept
=======

This tutorial introduces the theoretical background for the relationship 
between conduction velocity and tissue conductivity and how an iterative
method to tune the conductivity values to yield a desired velocity for a 
given setup can be derived. In addition, the definition of anisotropy ratio 
and its impact on simulation results are presented. A simple strategy to enforce
desired anisotropy ratio and conduction velocities is also presented.

Tuning conduction velocity
==========================

The relationship between conduction velocity and conductivity
-------------------------------------------------------------

A minimum requirement in modeling studies which aim at making case-specific predictions on 
ventricular electrophysiology is that activation sequences are carefully matched. 
Conduction velocity in the ventricles is orthotropic and may vary in space, thus
profoundly influencing propagation patterns.

As described in Section `[Tissue scale] <https://carp.medunigraz.at/carputils/cme-manual.html#tissue-scale>`_, 
the monodomain model `[1] <https://carp.medunigraz.at/carputils/cme-manual.html#monodomain-equation>`_ 
is equivalent to the bidomain model `[2] <https://carp.medunigraz.at/carputils/cme-manual.html#equation-eq-bidm>`_
if the monodomain conductivity tensor is given by the harmonic mean tensor, :math:`{\sigma_m}`

.. math::
   :label: equ-sigma-m
   
   \mathbf{\sigma_m} = (\mathbf{\sigma_i}+\mathbf{\sigma_e})*(\mathbf{\sigma_i}*\mathbf{\sigma_e})^{-1},


:math:`\mathbf{\sigma_i}` and :math:`\mathbf{\sigma_i}` denoting the intra- and extracellular conductivity, respectively. Conduction velocity, CV, is not a parameter in the bidomain equations
`[2] <https://carp.medunigraz.at/carputils/cme-manual.html#equation-eq-bidm>`_ and as such cannot be directly
parametrized. However, assuming a continuously propagating planar wavefront along a given direction, :math:`\zeta`,
one can derive a proportionality relationship [#Costa2013]_ between :math:`CV_\zeta` and :math:`\sigma_{m\zeta}`

.. math::
   :label: equ-sigma-cv

   CV_\zeta \propto \sqrt{{\sigma_{m\zeta}}}


For simplicity, the surface-to-volume ratio, :math:`\beta`, was ommited in this tutorial, 
   as it is kept constant when tuning conductivities. [#Costa2013]_

Conduction velocity as a function of simulation parameters
----------------------------------------------------------

Experimental measurements of conductivity values are scarce and the variation in 
measured values across studies is vast. These uncertainties inevitably arise due to the significant degree of 
biological variation and the substantial errors in the measurement techniques themselves.
Modeling and technical uncertainties may also have an
impact on model predictions. Particularly, CV depends on the specific model used to describe cellular dynamics, :math:`I_\mathrm{ion}`,
the spatial discretization, :math:`\Delta_x`, and on several implementation choices including the time-stepping scheme used, which
we represent as :math:`\\xi`. Thus, CV can be represented as a function 

.. math::
   :label: equ-cv-par

   CV_\zeta = CV_\zeta({\sigma_{m\zeta}},I_\mathrm{ion},\Delta_x,\\xi)

Iterative scheme for conduction velocity tuning
-----------------------------------------------

In most practical scenarios, :math:`I_\mathrm{ion}`, :math:`\Delta_x`, and :math:`\\xi` are parameters defined by users in the course of
selecting a simulation software, an ionic model and a provided mesh to describe the geometry. Thus, only :math:`{\sigma_{m\zeta}}` is left 
which can be tuned to achieve a close match between the pre-specified conduction velocity, :math:`CV_\zeta`, and the velocity, :math:`\\overline{CV_\zeta}`,
predicted by the simulation.

Using :eq:`equ-sigma-cv`, and :eq:`equ-sigma-m`, one can find unique monodomain conductivities along all
axes :math:`\zeta`, which yield the prescribed conduction velocities, :math:`CV_\zeta`, 
by iteratively refining conductivities based on :math:`\\overline{CV_\zeta}`,
measured in simple 1D cable simulations. The iterative update scheme we propose is given as

.. _fig-algorithm:

.. figure:: /images/02_03A_tuneCV-algorithm.png
   :width: 40 %
   :align: center

   Iterative update scheme to tune conductivity values based on a prescribed velocity :math:`CV_\zeta`. [#Costa2013]_

The algorithm is implemented as follows. Given an initial conductivity "guess", :math:`{\sigma_i}^0` and :math:`{\sigma_e}^0` and all other
simulation parameters defined by :math:`I_\mathrm{ion}`, :math:`\Delta_x`, and :math:`\\xi`, a simulation is run using a
1D cable and the conduction velocity is computed, as described in the figure below:

.. _fig-setup:

.. figure:: /images/02_03A_tuneCV-Setup1D.png
   :width: 50 %
   :align: center

   Setup for convergence testing. A pseudo 1D cable is created, which is in fact a 1 element tick cable of hexahedral elements,
   where :math:`h` is the edge length of each hexahedral, which also defines the spatial resolution. 
   Propagation is initiated by applying a transmembrane stimulus current at the left corner.
   :math:`T_0`, and :math:`T_1` correspond to wavefront
   arrival times recorded at locations :math:`x_0 = -0.25` cm and :math:`x_1 = 0.25` cm,
   respectively. CV is computed as :math:`CV = (x_1 - x_0)/(T_1 - T_0)`. 
   

Using :math:`factor = (CV_{\zeta}/\\overline{CV_{\zeta}})^2`, the conductivities :math:`{\sigma_i}[i+1]` and :math:`{\sigma_e}[i+1]`
are updated. A new simulation is then run with the new conductivities. 
These steps are repeated until the error in CV is below a given tolerance, :math:`stop_\mathrm{tol}`.


Problem Setup
-------------

This tutorial will run a simulation on a 1D cable and compute the simulated CV, :math:`\\bar{CV}`, and/or run 
the iterative scheme shown Figure :numref:`fig-algorithm` to compute the conductivity :math:`{\sigma_m}` which yields
the prescribed CV for the chosen simulation setup.

Usage
-----

The experiment specific options are:

.. code-block:: bash

  --resolution RESOLUTION
                        Spatial resolution
  --velocity VELOCITY   Desired conduction velocity (m/s)
  --gi GI               Intracellular conductivity (S/m)
  --ge GE               Extracellular conductivity (S/m)
  --ts TS               Choose time stepping method. 0: explicit Euler, 1:
                        crank nicolson, 2: second-order time stepping
  --dt DT               Integration time step on micro-seconds
  --converge {0,1}      0: Measure velocity with given setup or 1: Compute
                        conductivities that yield desired velocity with given
                        setup
  --compareTuning       run tuneCV with --resolution 100 200 and 400 with and
                        without tuning and make comparison plot
  --compareTimeStepping
                        run tuneCV with --resolution 100 200 and 400 with
                        Explicit Euler and Crank Nicolson and make comparison
                        plot
  --compareMassLumping  run tuneCV with --resolution 100 200 and 400 with and
                        without mass lumping and make comparison plot
  --compareModelPar     run tuneCV with --resolution 100 200 and 400 with
                        original Ten Tusher cell model and with reduced
                        sodium conductance and make comparison plot
  --ar AR               run tuneCV with for CV_f = 0.6 and CV_s = 0.3 m/s with
                        given anisotropy ratio (ar)



The user can either run single experiments, do a comparison of the multiple
resolutions with the options ``--compareTuning``, ``--compareTimeStepping``, ``--compareMassLumping``, and ``--compareModelPar``,
or compare the change in conductivities with different anisotropy ratios


Run 1D examples
---------------

To run the experiments of this tutorial do

.. code-block:: bash

    cd ${TUTORIALS}/02_EP_tissue/03A__study_prep_tuneCV


**compareTuning:**
^^^^^^^^^^^^^^^^^^

To compute the conductivities for resolutions of 100, 200, and 400 um with and without tuning, run:
For webGUI, remove RESOLUTION value
.. code-block:: bash

    ./run.py --compareTuning --visualize
    
After running this command you should see the following plot

.. _fig-tuning:

.. figure:: /images/02_03A_tuneCV-compareTuning.png
   :width: 50 %
   :align: center

   Simulated CV for different resolutions with and without the iterative tuning scheme.



Notice that CV varies with resolution in a non-linear fashion and CV is constant when the iterative scheme is applied.


**compareTimeStepping:**
^^^^^^^^^^^^^^^^^^^^^^^^

To compare the effect using the explicit Euler or Crank–Nicolson methods for the resolutions of 100, 200, and 400 um run:
For webGUI, remove RESOLUTION and TS values
.. code-block:: bash

    ./run.py --compareTimeStepping --visualize
    
.. _fig-timestepping:

.. figure:: /images/02_03A_tuneCV-compareTimeStepping.png
   :width: 50 %
   :align: center

   Simulated CV for different resolutions with and without the iterative tuning scheme.

Notice that CV varies with resolution virtually identically with both time-stepping schemes.

**compareMassLumping:**
^^^^^^^^^^^^^^^^^^^^^^^

Mass lumping is a common numerical technique in FEM to speed up computation. The mass matrix, :math:`M`, is made diagonal,
implying that its inverse is also diagonal, and solving the system is trivial.
The lumped mass matrix, :math:`M^L`, is computed by setting its main diagonal to be

.. math:: M_{ii}^L = \Sigma_{j=1}^N M_{ij}

To compare the effect of using mass lumping on CV for the resolutions of 100, 200, and 400 um run:
For webGUI, remove RESOLUTION and LUMPING values
.. code-block:: bash

    ./run.py --compareMassLumping --visualize
    
.. _fig-lumping:

.. figure:: /images/02_03A_tuneCV-compareMassLumping.png
   :width: 50 %
   :align: center

   Simulated CV for different resolutions with and without mass lumping.

Notice that CV varies with resolution in both cases, but the decrease in CV for coarser resolutions is much steeper with mass lumping.
   

**compareModelPar:**
^^^^^^^^^^^^^^^^^^^^

To compare the effect of modifying the sodium conduction on CV for the resolutions of 100, 200, and 400 um run:
For webGUI, remove RESOLUTION, MODEL and MODELPAR values
.. code-block:: bash

    ./run.py --compareModelPar --visualize

.. _fig-modelpar:

.. figure:: /images/02_03A_tuneCV-compareModelPar.png
   :width: 50 %
   :align: center

   Simulated CV for different resolutions with the original Ten Tusher model and with reduce sodium conductance

Notice that CV varies with resolution similarly in both cases, but the CV with reduced sodium conductance is lower than with the original model.
   


Run 2D and 3D examples
----------------------

As previously mentioned, conduction velocity in the ventricles is orthotropic, that is, CV is different 
in each axis :math:`\zeta`. Estimated average CVs in the longitudinal, transverse, and sheet normal directions,
:math:`CV_f`, :math:`CV_s` and :math:`CV_n`, are about 0.67 m/s, 0.3 m/s, and 0.17 m/s, respectively. [#Caldwell2009]_

**2D:**
^^^^^^^

To run a 2D simulation where :math:`CV_f > CV_s`, that is, an anisotropic setup, **two** sets of conductivities must be defined,
one for each axis :math:`\zeta=f,s`. To compute two sets of conductivities with tuneCV, call the example without
the "compare" flags twice with two different velocities and look at the conductivities output on the screen.

For the longitudinal direction, run 

.. code-block:: bash

    ./run.py --velocity 0.6 --converge 1

You should get an output close to the following:

.. code-block:: bash

    Conduction velocity: 0.6398 m/s [gi=0.1740, ge=0.6250, gm=0.1361]
    Conduction velocity: 0.6006 m/s [gi=0.1530, ge=0.5497, gm=0.1197]
    Conduction velocity: 0.6001 m/s [gi=0.1527, ge=0.5486, gm=0.1195]


where gi, ge, and gm are :math:`{\sigma_i}`, :math:`{\sigma_e}`, and :math:`{\sigma_m}`, respectively. 
The :math:`{\sigma_i}` and :math:`{\sigma_e}` pair will yield a CV of 0.6 m/s in both monodomain and bidomain simulations.

Now, for the transverse direction, run

.. code-block:: bash

    ./run.py --velocity 0.3 --converge 1

You should get an output close to the following:

.. code-block:: bash

    Conduction velocity: 0.6398 m/s [gi=0.1740, ge=0.6250, gm=0.1361]
    Conduction velocity: 0.3070 m/s [gi=0.0383, ge=0.1374, gm=0.0299]
    Conduction velocity: 0.3001 m/s [gi=0.0365, ge=0.1312, gm=0.0286]
    Conduction velocity: 0.3000 m/s [gi=0.0365, ge=0.1311, gm=0.0285]


In this case, the :math:`{\sigma_i}` and :math:`{\sigma_e}` pair will yield a CV of 0.3 m/s.

**3D:**
^^^^^^^

Similarly, to run a 3D simulation where :math:`CV_f > CV_s > CV_n`, that is, an orthotropic setup, **three** sets of conductivities must be defined,
one for each axis :math:`\zeta=f,s,n`. To compute the set of conductivities for the sheet normal direction, :math:`\zeta=n`, call the example once
more with :math:`CV_n = 0.2`

.. code-block:: bash

    ./run.py --velocity 0.2 --converge 1

You should get an output close to the following:

.. code-block:: bash

    Conduction velocity: 0.6398 m/s [gi=0.1740, ge=0.6250, gm=0.1361]
    Conduction velocity: 0.2038 m/s [gi=0.0170, ge=0.0611, gm=0.0133]
    Conduction velocity: 0.1998 m/s [gi=0.0164, ge=0.0588, gm=0.0128]
    Conduction velocity: 0.2000 m/s [gi=0.0164, ge=0.0590, gm=0.0128]


In this case, the :math:`{\sigma_i}` and :math:`{\sigma_e}` pair will yield a CV of 0.2 m/s.

Notice that for :math:`CV_s = 0.3` and :math:`CV_n = 0.2`, four iterations were required to converge to
   the desired CV, whereas :math:`CV_f = 0.6` required three. This is because, in this example, the initial conductivity :math:`{\sigma_{m,\zeta}}^0`
   (see Figure :numref:`fig-algorithm`) yields a velocity close to 0.6 m/s and is the same in all cases.
   Thus, the initial error in CV is larger for :math:`CV_s=0.3` and :math:`CV_n=0.2` than for :math:`CV_f=0.6`.


Tuning Anisotropy Ratios
==========================

 
Defining anisotropy ratios
--------------------------

According to experimental measurements, the intracellular
domain is more anisotropic than the extracellular domain, that is, the ratio between the 
longitudinal and transverse conductivities is larger in the intracellular space
than in the extracellular space.

In the intracellular domain, the anisotropy ratio between the longitudinal, :math:`\zeta=f`, 
and transverse, :math:`\zeta=s`, directions is defined as 
   
.. math::
   :label: equ-aniso-intra
   
   \\alpha_{ifs} = {\sigma_{if}}/{\sigma_{is}}

where :math:`{\sigma_{if}}` and :math:`{\sigma_{is}}` are the intracellular conductivities in the 
longitudinal and transverse direction, respectively.

Similarly, in the extracellular domain, the anisotropy ratio between the longitudinal, :math:`\zeta=f`, 
and transverse, :math:`\zeta=s`, directions is defined as 

.. math::
   :label: equ-aniso-extra
   
   \\alpha_{efs} = {\sigma_{ef}}/{\sigma_{es}}

where :math:`{\sigma_{ef}}` and :math:`{\sigma_{es}}` are the extracellular conductivities in the 
longitudinal and transverse direction, respectively.

Now, we can express the anisotropy ratio between the two domains as in the longitudinal and transverse directions as

.. math::
   :label: equ-aniso
   
   \\alpha_{lt} = \\alpha_{ifs}/\\alpha_{efs} = ({\sigma_{if}}*{\sigma_{es}})/({\sigma_{ef}}*{\sigma_{is}})


Computing conductivities with a fixed anisotropy ratio
------------------------------------------------------

There are many ways to enforce anisotropy ratios between conductivity values. In this tutorial, we use a 
fixed anisotropy strategy, where, given an anisotropy ratio, we compute an initial transverse conductivity
for extracellular domain. The intracellular conductivities as well as the longitudinal extracellular
conductivity are kept fixed at default values. Using Equation :eq:`equ-aniso` we obtain:

.. math::
   :label: equ-aniso-ges
   
   {\sigma_{es}} = (\\alpha_{lt}*{\sigma_{ef}}*{\sigma_{is}})/{\sigma_{if}}

We then tune the longitudinal and transverse conductivities separately using the iterative scheme. In this manner,
we enforce both the desired CV and anisotropy ratio.

The effect of equal and unequal anisotropy ratios
-------------------------------------------------

If :math:`\\alpha_{ifs} = \\alpha_{efs}`, then :math:`\\alpha_{lt} = 1`. In this case, the two domains
have **equal** anisotropy ratios. On the other hand, if :math:`\\alpha_{ifs} > \\alpha_{efs}`, 
then :math:`\\alpha_{lt} > 1`. In this case, the two domains have **unequal** anisotropy ratios.
Unequal anisotropy ratios can only be represented with the bidomain model `[2] <https://carp.medunigraz.at/carputils/cme-manual.html#equation-eq-bidm>`_.

While anisotropy ratios are of rather minor relevance when simulating impulse propagation in tissue [#Costa2013]_, they play a prominent
role when the stimulation of cardiac tissue via externally applied electric fields is studied.
In this case, virtual electrodes appear around the stimulus site, 
as show in the figure below.

.. _fig-anisotropy-vep:

.. figure:: /images/02_03A_tuneCV-anisotropy-vep.png
   :align: center

   Induced virtual electrodes in response to a strong hyperpolarizing extracellular stimulus for the equal and unequal anisotropy ratios.


Run anisotropy ratio examples
-----------------------------

To compute longitudinal and transverse conductivities for the case of **equal** anisotropy ratios and :math:`CV_f = 0.6` and :math:`CV_s = 0.3`, run
For webGUI, remove VELOCITY, GI, GE and CONVERGE values

.. code-block:: bash

    ./run.py --ar=1.0

You should get conductivities similar to these:

.. code-block:: bash

    g_if: 0.152700  g_is: 0.036500  g_ef: 0.548500  g_es: 0.131100


To compute longitudinal and transverse conductivities for the case of **unequal** anisotropy ratios and :math:`CV_f = 0.6` and :math:`CV_s = 0.3`, run

.. code-block:: bash

    ./run.py --ar=3.0
    
You should get conductivities similar to these:

.. code-block:: bash

   g_if: 0.152700  g_is: 0.031200  g_ef: 0.548500  g_es: 0.336100
   

   
.. rubric:: References

.. [#Costa2013]  *Costa CM, Hoetzl E, Rocha BM, Prassl AJ, Plank G.*,
                **Automatic Parameterization Strategy for Cardiac Electrophysiology Simulations**, 
                Comput Cardiol 40:373-376, 2013.
                `[Pubmed] <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3980367>`__
           
.. [#Caldwell2009]  *B. J. Caldwell, M. L. Trew, G. B. Sands, D. A. Hooks, I. J. LeGrice, B. H. Smaill*,
                   **Three distinct directions of intramural activation reveal nonuniform side-to-side electrical coupling of ventricular myocytes.**,
                   Circ Arrhythm Electrophysiol, vol. 2, no. 4, pp. 433-440, Aug 2009. 
                   `[Pubmed] <https://www.ncbi.nlm.nih.gov/pubmed/19808500>`__

"""
import os
EXAMPLE_DESCRIPTIVE_NAME = 'Tuning velocities and anisotropy ratios'
EXAMPLE_AUTHOR = 'Caroline Mendonca Costa <caroline.mendonca-costa@kcl.ac.uk>'
EXAMPLE_DIR = os.path.dirname(__file__)
GUIinclude = False

import sys
from datetime import date
from carputils import settings
from carputils import tools
from carputils.carpio import txt
from matplotlib import pyplot
import tuning as tn


def parser():
    parser = tools.standard_parser()
    group = parser.add_argument_group('experiment specific options')
    group.add_argument('--sres',
                        type = float, default = 100,
                        help = 'Spatial resolution')
    group.add_argument('--slen',
                        type = float, default = 2,
                        help = 'Length of slab')
    group.add_argument('--conv',
                        type = int, default = 0,
                        help = '0: Measure velocity with given setup or \n'
                               '1: Compute conductivities that yield desired velocity with given setup')                     
    group.add_argument('--gi',
                        type = float, default = 0.174,
                        help = 'Intracellular conductivity (S/m)')
    group.add_argument('--ge',
                        type= float, default = 0.625,
                        help = 'Extracellular conductivity (S/m)')
    group.add_argument('--beqm',
                        type= int, default = 1,
                        help = '0: use intracellular tensor or \n'
                               '1:  Use harmonic mean tensor')
    group.add_argument('--vel',
                        type = float, default = 0.6,
                        help = 'Target conduction velocity'
                                 '1:  Use harmonic mean tensor')
    group.add_argument('--ctol',
                        type = float, default = 0.0001,
                        help = 'Stop tolerance for CV tuning')
    group.add_argument('--max_it',
                        type = int, default = 30,
                        help = 'Maximum number of iterations')
    # Solver Paramters
    group.add_argument('--ts',
                        type = int, default = 0,
                        help = 'Choose time stepping method. \n'
                               '0: explicit Euler,\n'
                               '1: crank nicolson,\n'
                               '2: second-order time stepping')
    group.add_argument('--dt',
                        type = float, default = 10.,
                        help = 'Integration time step on micro-seconds')
    group.add_argument('--stol',
                        type = float, default = 1e-09,
                        help = 'solver tolerance')    
    group.add_argument('--lumping',
                        type = bool, default = False,
                        help = 'Use mass lumping')
    group.add_argument('--opsplt',
                        type = bool, default = True,
                        help = 'Operator splitting') 
    group.add_argument('--ar',
                        type = float, default = 0.,
                        help = 'run tuneCV with for CV_f = 0.6 and CV_s = 0.3 m/s with given anisotropy ratio (ar)')
    # Ionic region parameters
    group.add_argument('--beta',
                        type = float, default = 0.14,
                        help = 'surface-to-volume-ratio used')
    group.add_argument('--model',
                        type = str, default = 'UCLA_ZhengWT_MahajRev',
                        help = 'Ionic model')
    group.add_argument('--svinit',
                        type = str, default = '',
                        help = 'starting state for single cells')
    group.add_argument('--mpar',
                        type = str, default = '',
                        help = 'Ionic model parameters')
    group.add_argument('--plugin',
                        type = str, default = '',
                        help = 'Plugin for inonic model')
    group.add_argument('--ppar',
                        type = str, default = '',
                        help = 'parameters for ionic plugin for inonic model')

    group.add_argument('--Elem3D',
                        type = str, default = 0,
                        help = 'Element type')

   
    group.add_argument('--His',
                        type = float, default = 1000.,
                        help = 'Pacing frequency')
    
    # Decide to run external release
    group.add_argument('--ext',
                        type = int, default = 1,
                        help = '0: Run original carp models \n'
                               '1: Run externally compiled model')
    group.add_argument('--extDir',
                        type = str, default = '/mnt/f/MarkovIKr/SCN5A_data/scnaMatFiles/Quentin/NewDawn/Workers/execProtocols/Models/UCLA_RAB/ZhengMahaj/RevisedModels/UCLA_WT/',
                        help = 'External model location')
    
    group.add_argument('--setInit',
                        type = int, default = 0,
                        help = '0: donot initialize \n'
                               '1: Initialize')
    return parser


def jobID(args):
    """
    Generate name of top level output directory.
    """
    today = date.today()
    if args.conv:
        return '{}_tuneCV_{}_conv'.format(today.isoformat(), args.sres)
    else:
        return '{}_tuneCV_{}'.format(today.isoformat(), args.sres)


@tools.carpexample(parser, jobID, clean_pattern='{}*|(.log)|(.dat)|imp_*|(mesh)'.format(date.today()))
def run(args, job):
    parameterFile = 'runTune.par'
    
    if args.ar:
        # user defined anisotropy ratio
        if '--vel' in sys.argv or '--gi' in sys.argv or '--ge' in sys.argv or '--conv' in sys.argv:
            raise Exception('cannot set --velocity or --gi or --ge or --converge with --ar')
    # -------------------------------------------------------------------------

        # set default values
        gil = 0.174
        gel = 0.625
        git = 0.019
        get = (args.ar * gel * git) / gil

        args.vel = 0.6
        args.gi = gil; args.ge = gel; args.conv = 1; 
        gil_, gel_,_,_ = tn.iterateCV(args, job,parameterFile)
        args.vel = 0.3
        args.gi = git; args.ge = get; args.conv = 1; 
        git_, get_,_,_ = tn.iterateCV(args, job,parameterFile)

        print('g_if: %f\tg_is: %f\tg_ef: %f\tg_es: %f\n' % (gil_, git_, gel_, get_))
    else:
        tn.iterateCV(args, job,parameterFile)

if __name__ == '__main__':
    run()
