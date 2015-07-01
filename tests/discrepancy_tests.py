# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 17:06:14 2013

@author: zirbel
"""

from unittest import TestCase
import numpy as np

from fr3d.geometry.discrepancy import discrepancy

from fr3d.data import Atom
from fr3d.data import Component

nt10_0 = Component([Atom(type='N', name='N1',x=19.679762,y=146.510834,z=106.455235),
                    Atom(type='C', name='C2',x=19.459715,y=146.692932,z=107.816836),
                    Atom(type='O', name='O2',x=19.354777,y=147.785831,z=108.333165),
                    Atom(type='N', name='N3',x=19.375964,y=145.491465,z=108.503069),
                    Atom(type='C', name='C4',x=19.482438,y=144.177244,z=108.006327),
                    Atom(type='O', name='O4',x=19.384232,y=143.219519,z=108.751357),
                    Atom(type='C', name='C6',x=19.800699,y=145.278872,z=105.856545),
                    Atom(type='C', name='C5',x=19.712414,y=144.129302,z=106.564466),
                    Atom(type='H', name='H5',x=19.808116,y=143.162711,z=106.089600),
                    Atom(type='H', name='H1',x=19.747568,y=147.360820,z=105.915883),
                    Atom(type='H', name='H3',x=19.215883,y=145.577268,z=109.499783),
                    Atom(type='H', name='H6',x=19.970458,y=145.295294,z=104.785373),
                    Atom(type='C', name='C1*',x=19.841000,y=147.761000,z=105.646000),
                    Atom(type='C', name='C2*',x=20.264000,y=147.625000,z=104.181000),
                    Atom(type='O', name='O2*',x=21.328000,y=148.530000,z=104.003000),
                    Atom(type='C', name='C3*',x=18.978000,y=148.021000,z=103.445000),
                    Atom(type='O', name='O3*',x=19.159000,y=148.571000,z=102.125000),
                    Atom(type='C', name='C4*',x=18.405000,y=149.086000,z=104.371000),
                    Atom(type='O', name='O4*',x=18.630000,y=148.501000,z=105.679000),
                    Atom(type='C', name='C5*',x=16.942000,y=149.450000,z=104.201000),
                    Atom(type='O', name='O5*',x=16.030000,y=148.355000,z=104.290000),
                    Atom(type='P', name='P',x=19.695000,y=146.523000,z=106.436000),
                    Atom(type='O', name='O1P',x=19.695000,y=146.523000,z=106.436000),
                    Atom(type='O', name='O2P',x=19.695000,y=146.523000,z=106.436000)],
                   type='rna', pdb='1S72', sequence='U', chain='0', number='10')
nt11_0 = Component([Atom(type='N', name='N9',x=24.987289,y=147.803631,z=104.303074),
                    Atom(type='C', name='C4',x=26.110439,y=147.121416,z=104.714354),
                    Atom(type='N', name='N3',x=27.347788,y=147.594850,z=104.920430),
                    Atom(type='N', name='N1',x=27.898898,y=145.306543,z=105.516171),
                    Atom(type='C', name='C6',x=26.649175,y=144.875521,z=105.298527),
                    Atom(type='N', name='N6',x=26.358055,y=143.568865,z=105.491067),
                    Atom(type='C', name='C8',x=23.966202,y=146.889301,z=104.238608),
                    Atom(type='C', name='C5',x=25.671970,y=145.801703,z=104.873618),
                    Atom(type='C', name='C2',x=28.167642,y=146.615899,z=105.317633),
                    Atom(type='N', name='N7',x=24.332543,y=145.661269,z=104.575517),
                    Atom(type='H', name='H2',x=29.197130,y=146.907332,z=105.509753),
                    Atom(type='H', name='H8',x=22.966837,y=147.174282,z=103.937596),
                    Atom(type='H', name='H9',x=24.941473,y=148.789794,z=104.093050),
                    Atom(type='1', name='1H6',x=27.087431,y=142.947086,z=105.794198),
                    Atom(type='2', name='2H6',x=25.424508,y=143.232853,z=105.331318),
                    Atom(type='C', name='C1*',x=24.920000,y=149.247000,z=104.019000),
                    Atom(type='C', name='C2*',x=25.979000,y=149.737000,z=103.029000),
                    Atom(type='O', name='O2*',x=26.289000,y=151.092000,z=103.272000),
                    Atom(type='C', name='C3*',x=25.261000,y=149.504000,z=101.712000),
                    Atom(type='O', name='O3*',x=25.713000,y=150.275000,z=100.620000),
                    Atom(type='C', name='C4*',x=23.822000,y=149.841000,z=102.032000),
                    Atom(type='O', name='O4*',x=23.656000,y=149.501000,z=103.432000),
                    Atom(type='C', name='C5*',x=22.872000,y=149.049000,z=101.168000),
                    Atom(type='O', name='O5*',x=21.538000,y=149.169000,z=101.650000),
                    Atom(type='P', name='P',x=20.413000,y=148.134000,z=101.205000),
                    Atom(type='O', name='O1P',x=20.046000,y=148.458000,z=99.805000),
                    Atom(type='O', name='O2P',x=20.863000,y=146.760000,z=101.535000)],
                   type='rna', pdb='1S72', sequence='A', chain='0', number='11')
nt12_0 = Component([Atom(type='N', name='N1',x=28.124614,y=145.500212,z=101.946248),
                    Atom(type='C', name='C2',x=28.291980,y=144.186024,z=102.371027),
                    Atom(type='O', name='O2',x=29.361265,y=143.730926,z=102.719664),
                    Atom(type='N', name='N3',x=27.105856,y=143.469055,z=102.348201),
                    Atom(type='C', name='C4',x=25.818979,y=143.895789,z=101.964649),
                    Atom(type='O', name='O4',x=24.871809,y=143.132035,z=102.003377),
                    Atom(type='C', name='C6',x=26.919410,y=146.029395,z=101.548459),
                    Atom(type='C', name='C5',x=25.784087,y=145.293563,z=101.541375),
                    Atom(type='H', name='H5',x=24.838534,y=145.713744,z=101.227795),
                    Atom(type='H', name='H1',x=28.964200,y=146.059737,z=101.946910),
                    Atom(type='H', name='H3',x=27.181710,y=142.505043,z=102.650444),
                    Atom(type='H', name='H6',x=26.944612,y=147.070012,z=101.243538),
                    Atom(type='C', name='C1*',x=29.366000,y=146.320000,z=101.853000),
                    Atom(type='C', name='C2*',x=30.207000,y=145.911000,z=100.643000),
                    Atom(type='O', name='O2*',x=31.575000,y=146.193000,z=100.914000),
                    Atom(type='C', name='C3*',x=29.619000,y=146.822000,z=99.579000),
                    Atom(type='O', name='O3*',x=30.454000,y=146.919000,z=98.434000),
                    Atom(type='C', name='C4*',x=29.495000,y=148.124000,z=100.357000),
                    Atom(type='O', name='O4*',x=29.010000,y=147.683000,z=101.655000),
                    Atom(type='C', name='C5*',x=28.560000,y=149.166000,z=99.798000),
                    Atom(type='O', name='O5*',x=27.261000,y=148.601000,z=99.560000),
                    Atom(type='P', name='P',x=26.018000,y=149.543000,z=99.231000),
                    Atom(type='O', name='O1P',x=26.451000,y=150.604000,z=98.281000),
                    Atom(type='O', name='O2P',x=24.872000,y=148.661000,z=98.883000)],type='rna', pdb='1S72', sequence='U', chain='0', number='12')
nt13_0 = Component([Atom(type='N', name='N9',x=25.799709,y=142.530131,z=98.872358),
                    Atom(type='C', name='C4',x=24.662685,y=141.783599,z=99.027532),
                    Atom(type='N', name='N3',x=24.623209,y=140.466993,z=99.376271),
                    Atom(type='N', name='N1',x=22.284330,y=140.799272,z=99.198561),
                    Atom(type='C', name='C6',x=22.261771,y=142.183818,z=98.829380),
                    Atom(type='O', name='O6',x=21.199539,y=142.742705,z=98.640139),
                    Atom(type='C', name='C8',x=25.402596,y=143.799650,z=98.519473),
                    Atom(type='C', name='C5',x=23.619314,y=142.669099,z=98.752125),
                    Atom(type='C', name='C2',x=23.391426,y=140.025691,z=99.446628),
                    Atom(type='N', name='N7',x=24.089116,y=143.924114,z=98.436204),
                    Atom(type='N', name='N2',x=23.166305,y=138.723928,z=99.784329),
                    Atom(type='H', name='H1',x=21.359967,y=140.391679,z=99.271691),
                    Atom(type='H', name='H8',x=26.117331,y=144.590102,z=98.336489),
                    Atom(type='H', name='H9',x=26.741969,y=142.190520,z=98.998684),
                    Atom(type='1', name='1H2',x=23.968829,y=138.148541,z=99.968041),
                    Atom(type='2', name='2H2',x=22.245672,y=138.332290,z=99.853357),
                    Atom(type='C', name='C1*',x=27.167000,y=142.051000,z=99.065000),
                    Atom(type='C', name='C2*',x=27.576000,y=141.058000,z=97.976000),
                    Atom(type='O', name='O2*',x=28.321000,y=140.031000,z=98.604000),
                    Atom(type='C', name='C3*',x=28.448000,y=141.936000,z=97.084000),
                    Atom(type='O', name='O3*',x=29.320000,y=141.119000,z=96.303000),
                    Atom(type='C', name='C4*',x=29.191000,y=142.727000,z=98.146000),
                    Atom(type='O', name='O4*',x=28.086000,y=143.129000,z=98.994000),
                    Atom(type='C', name='C5*',x=30.041000,y=143.915000,z=97.750000),
                    Atom(type='O', name='O5*',x=29.244000,y=144.996000,z=97.263000),
                    Atom(type='P', name='P',x=29.906000,y=146.425000,z=97.004000),
                    Atom(type='O', name='O1P',x=31.117000,y=146.240000,z=96.171000),
                    Atom(type='O', name='O2P',x=28.833000,y=147.347000,z=96.553000)],type='rna', pdb='1S72', sequence='G', chain='0', number='13')
nt14_0 = Component([Atom(type='N', name='N1',x=24.431662,y=138.047168,z=96.415061),
                    Atom(type='C', name='C2',x=23.033174,y=137.804837,z=96.331428),
                    Atom(type='O', name='O2',x=22.624966,y=136.670658,z=96.492930),
                    Atom(type='N', name='N3',x=22.246402,y=138.905955,z=96.068044),
                    Atom(type='C', name='C4',x=22.801415,y=140.089007,z=95.911224),
                    Atom(type='N', name='N4',x=21.965901,y=141.128407,z=95.654238),
                    Atom(type='C', name='C6',x=25.002181,y=139.266355,z=96.253322),
                    Atom(type='C', name='C5',x=24.217298,y=140.346612,z=95.993754),
                    Atom(type='H', name='H1',x=24.994454,y=137.230789,z=96.608487),
                    Atom(type='H', name='H6',x=26.082590,y=139.320759,z=96.341151),
                    Atom(type='H', name='H5',x=24.638876,y=141.334042,z=95.859111),
                    Atom(type='1', name='1H4',x=22.303849,y=142.063684,z=95.521178),
                    Atom(type='2', name='2H4',x=20.979785,y=140.934055,z=95.599370),
                    Atom(type='C', name='C1*',x=25.320000,y=136.868000,z=96.718000),
                    Atom(type='C', name='C2*',x=25.519000,y=135.955000,z=95.508000),
                    Atom(type='O', name='O2*',x=25.725000,y=134.638000,z=95.976000),
                    Atom(type='C', name='C3*',x=26.766000,y=136.567000,z=94.886000),
                    Atom(type='O', name='O3*',x=27.464000,y=135.669000,z=94.031000),
                    Atom(type='C', name='C4*',x=27.592000,y=136.924000,z=96.113000),
                    Atom(type='O', name='O4*',x=26.605000,y=137.350000,z=97.090000),
                    Atom(type='C', name='C5*',x=28.601000,y=138.031000,z=95.915000),
                    Atom(type='O', name='O5*',x=27.952000,y=139.206000,z=95.394000),
                    Atom(type='P', name='P',x=28.794000,y=140.473000,z=94.926000),
                    Atom(type='O', name='O1P',x=29.985000,y=139.988000,z=94.182000),
                    Atom(type='O', name='O2P',x=27.865000,y=141.432000,z=94.268000)],type='rna', pdb='1S72', sequence='C', chain='0', number='14')
nt15_0 = Component([Atom(type='N', name='N1',x=21.976931,y=136.027203,z=92.627909),
                    Atom(type='C', name='C2',x=20.646444,y=136.521737,z=92.545886),
                    Atom(type='O', name='O2',x=19.729958,y=135.726030,z=92.468542),
                    Atom(type='N', name='N3',x=20.512471,y=137.893829,z=92.561881),
                    Atom(type='C', name='C4',x=21.580262,y=138.658339,z=92.649011),
                    Atom(type='N', name='N4',x=21.373461,y=140.000576,z=92.659472),
                    Atom(type='C', name='C6',x=23.076085,y=136.815796,z=92.717629),
                    Atom(type='C', name='C5',x=22.933387,y=138.168491,z=92.732670),
                    Atom(type='H', name='H1',x=22.058497,y=135.020325,z=92.615025),
                    Atom(type='H', name='H6',x=24.038004,y=136.316341,z=92.774327),
                    Atom(type='H', name='H5',x=23.788564,y=138.827699,z=92.803312),
                    Atom(type='1', name='1H4',x=22.130371,y=140.656065,z=92.723318),
                    Atom(type='2', name='2H4',x=20.423635,y=140.328533,z=92.600456),
                    Atom(type='C', name='C1*',x=22.130000,y=134.529000,z=92.608000),
                    Atom(type='C', name='C2*',x=22.139000,y=133.962000,z=91.189000),
                    Atom(type='O', name='O2*',x=21.568000,y=132.673000,z=91.218000),
                    Atom(type='C', name='C3*',x=23.630000,y=133.916000,z=90.908000),
                    Atom(type='O', name='O3*',x=23.957000,y=133.038000,z=89.845000),
                    Atom(type='C', name='C4*',x=24.184000,y=133.459000,z=92.247000),
                    Atom(type='O', name='O4*',x=23.370000,y=134.176000,z=93.209000),
                    Atom(type='C', name='C5*',x=25.632000,y=133.779000,z=92.486000),
                    Atom(type='O', name='O5*',x=25.894000,y=135.148000,z=92.141000),
                    Atom(type='P', name='P',x=27.312000,y=135.803000,z=92.435000),
                    Atom(type='O', name='O1P',x=28.369000,y=134.948000,z=91.840000),
                    Atom(type='O', name='O2P',x=27.224000,y=137.229000,z=92.058000)],type='rna', pdb='1S72', sequence='C', chain='0', number='15')
nt16_0 = Component([Atom(type='N', name='N9',x=19.563275,y=136.400283,z=88.640605),
                    Atom(type='C', name='C4',x=18.863564,y=137.573177,z=88.815991),
                    Atom(type='N', name='N3',x=17.545204,y=137.786243,z=88.696964),
                    Atom(type='N', name='N1',x=18.074986,y=140.082599,z=89.280967),
                    Atom(type='C', name='C6',x=19.385459,y=139.825872,z=89.389166),
                    Atom(type='N', name='N6',x=20.229849,y=140.828623,z=89.721834),
                    Atom(type='C', name='C8',x=20.885249,y=136.685719,z=88.872300),
                    Atom(type='C', name='C5',x=19.845589,y=138.512394,z=89.151790),
                    Atom(type='C', name='C2',x=17.249403,y=139.066075,z=88.947581),
                    Atom(type='N', name='N7',x=21.106424,y=137.954015,z=89.185802),
                    Atom(type='H', name='H2',x=16.198310,y=139.333611,z=88.874668),
                    Atom(type='H', name='H8',x=21.651778,y=135.925774,z=88.798331),
                    Atom(type='H', name='H9',x=19.163502,y=135.508267,z=88.389180),
                    Atom(type='1', name='1H6',x=19.855344,y=141.747757,z=89.881291),
                    Atom(type='2', name='2H6',x=21.215050,y=140.649212,z=89.806238),
                    Atom(type='C', name='C1*',x=18.951000,y=135.144000,z=88.195000),
                    Atom(type='C', name='C2*',x=18.823000,y=135.099000,z=86.671000),
                    Atom(type='O', name='O2*',x=17.664000,y=134.389000,z=86.303000),
                    Atom(type='C', name='C3*',x=20.108000,y=134.382000,z=86.285000),
                    Atom(type='O', name='O3*',x=20.018000,y=133.798000,z=84.986000),
                    Atom(type='C', name='C4*',x=20.242000,y=133.353000,z=87.400000),
                    Atom(type='O', name='O4*',x=19.799000,y=134.075000,z=88.584000),
                    Atom(type='C', name='C5*',x=21.631000,y=132.816000,z=87.648000),
                    Atom(type='O', name='O5*',x=22.534000,y=133.896000,z=87.945000),
                    Atom(type='P', name='P',x=24.050000,y=133.611000,z=88.346000),
                    Atom(type='O', name='O1P',x=24.568000,y=132.504000,z=87.503000),
                    Atom(type='O', name='O2P',x=24.763000,y=134.921000,z=88.359000)],type='rna', pdb='1S72', sequence='A', chain='0', number='16')
nt17_0 = Component([Atom(type='N', name='N9',x=17.731718,y=139.432830,z=85.401299),
                    Atom(type='C', name='C4',x=17.537494,y=140.736120,z=85.772686),
                    Atom(type='N', name='N3',x=16.360699,y=141.416887,z=85.681382),
                    Atom(type='N', name='N1',x=17.654502,y=143.178040,z=86.596128),
                    Atom(type='C', name='C6',x=18.919214,y=142.514871,z=86.716466),
                    Atom(type='O', name='O6',x=19.876659,y=143.113837,z=87.164846),
                    Atom(type='C', name='C8',x=19.051108,y=139.134536,z=85.654545),
                    Atom(type='C', name='C5',x=18.783020,y=141.158744,z=86.240272),
                    Atom(type='C', name='C2',x=16.483392,y=142.647866,z=86.113384),
                    Atom(type='N', name='N7',x=19.721358,y=140.154170,z=86.163236),
                    Atom(type='N', name='N2',x=15.398836,y=143.474100,z=86.093756),
                    Atom(type='H', name='H1',x=17.677290,y=144.139432,z=86.914143),
                    Atom(type='H', name='H8',x=19.464358,y=138.157141,z=85.447020),
                    Atom(type='H', name='H9',x=17.023969,y=138.824625,z=85.016167),
                    Atom(type='1', name='1H2',x=14.536651,y=143.096982,z=85.742639),
                    Atom(type='2', name='2H2',x=15.434504,y=144.424566,z=86.411644),
                    Atom(type='C', name='C1*',x=16.735000,y=138.521000,z=84.833000),
                    Atom(type='C', name='C2*',x=16.307000,y=138.885000,z=83.407000),
                    Atom(type='O', name='O2*',x=14.916000,y=138.655000,z=83.294000),
                    Atom(type='C', name='C3*',x=17.113000,y=137.887000,z=82.583000),
                    Atom(type='O', name='O3*',x=16.493000,y=137.691000,z=81.313000),
                    Atom(type='C', name='C4*',x=17.014000,y=136.647000,z=83.460000),
                    Atom(type='O', name='O4*',x=17.251000,y=137.196000,z=84.783000),
                    Atom(type='C', name='C5*',x=17.955000,y=135.485000,z=83.207000),
                    Atom(type='O', name='O5*',x=19.307000,y=135.811000,z=83.583000),
                    Atom(type='P', name='P',x=20.423000,y=134.666000,z=83.688000),
                    Atom(type='O', name='O1P',x=20.298000,y=133.748000,z=82.525000),
                    Atom(type='O', name='O2P',x=21.723000,y=135.328000,z=83.955000)],
                   type='rna', pdb='1S72', sequence='G', chain='0', number='17')
nt18_0 = Component([Atom(type='N', name='N1',x=17.002058,y=143.265799,z=82.804934),
                    Atom(type='C', name='C2',x=17.523475,y=144.428013,z=83.436495),
                    Atom(type='O', name='O2',x=16.797122,y=145.393170,z=83.577775),
                    Atom(type='N', name='N3',x=18.839122,y=144.353528,z=83.841917),
                    Atom(type='C', name='C4',x=19.530101,y=143.252575,z=83.635289),
                    Atom(type='N', name='N4',x=20.819910,y=143.249953,z=84.060566),
                    Atom(type='C', name='C6',x=17.714887,y=142.132273,z=82.592657),
                    Atom(type='C', name='C5',x=19.012326,y=142.068689,z=82.996368),
                    Atom(type='H', name='H1',x=16.037218,y=143.335988,z=82.513764),
                    Atom(type='H', name='H6',x=17.202344,y=141.312566,z=82.099362),
                    Atom(type='H', name='H5',x=19.610545,y=141.180157,z=82.843013),
                    Atom(type='1', name='1H4',x=21.418197,y=142.453088,z=83.943453),
                    Atom(type='2', name='2H4',x=21.168182,y=144.083284,z=84.504948),
                    Atom(type='C', name='C1*',x=15.546000,y=143.310000,z=82.422000),
                    Atom(type='C', name='C2*',x=15.281000,y=144.036000,z=81.102000),
                    Atom(type='O', name='O2*',x=14.010000,y=144.651000,z=81.187000),
                    Atom(type='C', name='C3*',x=15.297000,y=142.873000,z=80.119000),
                    Atom(type='O', name='O3*',x=14.588000,y=143.173000,z=78.919000),
                    Atom(type='C', name='C4*',x=14.600000,y=141.783000,z=80.921000),
                    Atom(type='O', name='O4*',x=15.075000,y=141.982000,z=82.276000),
                    Atom(type='C', name='C5*',x=14.868000,y=140.358000,z=80.510000),
                    Atom(type='O', name='O5*',x=16.275000,y=140.088000,z=80.537000),
                    Atom(type='P', name='P',x=16.839000,y=138.673000,z=80.084000),
                    Atom(type='O', name='O1P',x=16.035000,y=138.212000,z=78.928000),
                    Atom(type='O', name='O2P',x=18.307000,y=138.788000,z=79.949000)],
                   type='rna', pdb='1S72', sequence='C', chain='0', number='18')
nt19_0 = Component([Atom(type='N', name='N1',x=18.396146,y=147.260932,z=80.217077),
                    Atom(type='C', name='C2',x=19.483911,y=147.869323,z=80.835273),
                    Atom(type='O', name='O2',x=19.557234,y=149.065306,z=81.025757),
                    Atom(type='N', name='N3',x=20.460540,y=146.956683,z=81.202339),
                    Atom(type='C', name='C4',x=20.477048,y=145.558321,z=81.030358),
                    Atom(type='O', name='O4',x=21.422511,y=144.894339,z=81.413996),
                    Atom(type='C', name='C6',x=18.299443,y=145.907425,z=79.994695),
                    Atom(type='C', name='C5',x=19.277167,y=145.050670,z=80.369505),
                    Atom(type='H', name='H5',x=19.196964,y=143.987068,z=80.192663),
                    Atom(type='H', name='H1',x=17.659798,y=147.889306,z=79.932626),
                    Atom(type='H', name='H3',x=21.272535,y=147.357256,z=81.656906),
                    Atom(type='H', name='H6',x=17.390802,y=145.579349,z=79.501508),
                    Atom(type='C', name='C1*',x=17.292000,y=148.160000,z=79.772000),
                    Atom(type='C', name='C2*',x=17.615000,y=148.777000,z=78.413000),
                    Atom(type='O', name='O2*',x=17.001000,y=150.044000,z=78.323000),
                    Atom(type='C', name='C3*',x=16.988000,y=147.751000,z=77.482000),
                    Atom(type='O', name='O3*',x=16.783000,y=148.262000,z=76.179000),
                    Atom(type='C', name='C4*',x=15.694000,y=147.440000,z=78.209000),
                    Atom(type='O', name='O4*',x=16.105000,y=147.396000,z=79.603000),
                    Atom(type='C', name='C5*',x=14.990000,y=146.159000,z=77.826000),
                    Atom(type='O', name='O5*',x=15.901000,y=145.048000,z=77.860000),
                    Atom(type='P', name='P',x=15.390000,y=143.568000,z=77.576000),
                    Atom(type='O', name='O1P',x=14.377000,y=143.628000,z=76.498000),
                    Atom(type='O', name='O2P',x=16.560000,y=142.674000,z=77.416000)],
                   type='rna', pdb='1S72', sequence='U', chain='0', number='19')

nt10_0.infer_hydrogens()
nt11_0.infer_hydrogens()
nt12_0.infer_hydrogens()
nt13_0.infer_hydrogens()
nt14_0.infer_hydrogens()
nt15_0.infer_hydrogens()
nt16_0.infer_hydrogens()
nt17_0.infer_hydrogens()
nt18_0.infer_hydrogens()
nt19_0.infer_hydrogens()

# print nt10_0.rotation_matrix


class DiscrepancyTest(TestCase):
    def test_can_compute_simple_discrepancy(self):
        val = discrepancy([nt10_0, nt11_0, nt12_0, nt13_0, nt14_0],
                          [nt15_0,nt16_0,nt17_0,nt18_0,nt19_0])
        np.testing.assert_almost_equal(0.9402, val, decimal=3)

    def test_can_compute_another_discrepancy(self):
        val = discrepancy([nt10_0, nt12_0, nt14_0, nt16_0, nt18_0],
                          [nt11_0,nt13_0,nt15_0,nt17_0,nt19_0])
        np.testing.assert_almost_equal(1.1414, val, decimal=3)


# d = discrepancy([nt10_0, nt11_0, nt12_0, nt13_0, nt14_0],[nt15_0,nt16_0,nt17_0,nt18_0,nt19_0])
# print "Discrepancy", d
# location error (squared) is 17.7790, from Matlab
# discrepancy should be about 0.9402
# discrepancy is sqrt(L^2 + A1^2 + A2^2 + ... + A5^2)/5, where L is location error and A1 ... A5 are angles between bases in radians
# angles between bases should be about     0.9014    0.3105    0.1999    0.2084    0.2959
# or else twice that, not quite sure

# rotation matrix should be about:
#   -0.9933    0.0681   -0.0929
#   -0.0533   -0.9867   -0.1537
#   -0.1022   -0.1478    0.9837

# d = discrepancy([nt10_0, nt11_0, nt12_0, nt13_0, nt14_0],[nt15_0,nt16_0,nt17_0,nt18_0,nt19_0],['P'])
# print "Discrepancy", d

# d = discrepancy([nt10_0, nt12_0, nt14_0, nt16_0, nt18_0],[nt11_0,nt13_0,nt15_0,nt17_0,nt19_0])
# print "Discrepancy", d

# location error (squared) is 30.0211, from Matlab
# discrepancy should be about 1.1414
# discrepancy is sqrt(L^2 + A1^2 + A2^2 + ... + A5^2)/5, where L is location error and A1 ... A5 are angles between bases in radians
# angles between bases should be about     0.7006    0.1180    0.2260    0.2561    0.1276

# or else twice that, not quite sure

# rotation matrix should be about:
#    0.6357    0.7690    0.0672
#   -0.7708    0.6371    0.0006
#   -0.0424   -0.0522    0.9977
