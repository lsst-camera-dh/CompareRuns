# bot_results_compare.py  11 Oct 2020
# Quick comparison of BOT EO test results for two runs, averaged over raft types

# This needs the datacat-utilities package from lsst-camera-dh
# (see https://github.com/lsst-camera-dh/datacat-utilities)

import numpy as np
from get_EO_analysis_results import get_EO_analysis_results
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import scipy.stats as stats

from astropy.table import Table
import pylab
import numpy

raft_types = {}
raft_types['e2v'] = 'R11 R12 R13 R14 R21 R22 R23 R24 R30 R31 R32 R33 R34'.split()
raft_types['ITL'] = 'R01 R02 R03 R10 R20 R41 R42 R43'.split()
raft_types['Corner'] = 'R00 R04 R40 R44'.split()

EO = get_EO_analysis_results(db='Prod')

def compare2runs( botrun1, botrun2 ):
    dev_list1, data1 = EO.get_tests(site_type='I&T-Raft', run=botrun1)
    botresults1 = EO.get_all_results(data1)
    dev_list2, data2 = EO.get_tests(site_type='I&T-Raft', run=botrun2)
    botresults2 = EO.get_all_results(data2)

    fplist1 = list(botresults1['read_noise'].keys())  # get the list of rafts [assume that at least read_noise is available
    fplist2 = list(botresults2['read_noise'].keys())  # get the list of rafts [assume that at least read_noise is available
    # Find the runs that are in common
    fplist = []
    for raft in fplist1:
        if raft in fplist2:
            fplist.append(raft)

    # Look up the e2v rafts
    included_e2v = []
    for raft in raft_types['e2v']:
        if raft in fplist:
            included_e2v.append(raft)
    # and the ITL rafts
    included_ITL = []
    for raft in raft_types['ITL']:
        if raft in fplist:
            included_ITL.append(raft)

    # Make lists of the available test results for both runs
    eotests1 = list(botresults1.keys())
    eotests2 = list(botresults2.keys())

    # Find the tests in common in both runs
    eotests = []
    for eotest in eotests1:
        if eotest in eotests2:
            eotests.append(eotest)
    eotests.sort()

    comp = {}

    for eotest in eotests:
        print(eotest)
        results1_e2v = []
        for raft in included_e2v:
            for sensor in botresults1[eotest][raft].keys():
                results1_e2v.extend(botresults1[eotest][raft][sensor])
        results1_ITL = []
        for raft in included_ITL:
            for sensor in botresults1[eotest][raft].keys():
                results1_ITL.extend(botresults1[eotest][raft][sensor])
        results2_e2v = []
        for raft in included_e2v:
            for sensor in botresults2[eotest][raft].keys():
                results2_e2v.extend(botresults2[eotest][raft][sensor])
        results2_ITL = []
        for raft in included_ITL:
            for sensor in botresults2[eotest][raft].keys():
                results2_ITL.extend(botresults2[eotest][raft][sensor])
        results1_e2v = np.array(results1_e2v)
        results1_ITL = np.array(results1_ITL)
        results2_e2v = np.array(results2_e2v)
        results2_ITL = np.array(results2_ITL)
        temp = {}
        if eotest == 'tearing_detections':
            temp['e2v'] = (np.sum(results1_e2v), 0, np.sum(results2_e2v), 0)
            temp['ITL'] = (np.sum(results1_ITL), 0, np.sum(results2_ITL), 0)
        else:
            temp['e2v'] = (np.median(results1_e2v), stats.median_absolute_deviation(results1_e2v)/np.sqrt(len(results1_e2v)), np.median(results2_e2v), stats.median_absolute_deviation(results2_e2v)/np.sqrt(len(results2_e2v)))
            temp['ITL'] = (np.median(results1_ITL), stats.median_absolute_deviation(results1_ITL)/np.sqrt(len(results1_ITL)), np.median(results2_ITL), stats.median_absolute_deviation(results2_ITL)/np.sqrt(len(results2_ITL)))
        comp[eotest] = temp


    for kind in ['e2v', 'ITL']:
        print(' ')
        print(kind)
        print('Test                      ----- ' + botrun1 + ' ------   ----- ' + botrun2 + ' ------')
        for eotest in eotests:
            m1 = comp[eotest][kind][0]
            s1 = comp[eotest][kind][1]
            m2 = comp[eotest][kind][2]
            s2 = comp[eotest][kind][3]
            flag = ''
            if np.abs(m2-m1) > 2*np.sqrt(s1**2 + s2**2):
                flag = ' <--'
            print(eotest.ljust(25, ' ') + '%9.2e %8.1e %9.2e %8.1e' % (comp[eotest][kind][0], comp[eotest][kind][1], comp[eotest][kind][2], comp[eotest][kind][3]) + flag)

    print(' ')
    outstr = ''
    for elt in included_e2v:
        outstr += elt + ' '
    print('e2v rafts:  ' + outstr)
    outstr = ''
    for elt in included_ITL:
        outstr += elt + ' '
    print('ITL rafts:  ' + outstr)
    return eotests, comp


def makeplot( botrun1, botrun2, eotests, comp ):
    tablearr = []
    for kind in ['e2v', 'ITL']:
        tablearr.append(
            ( Table( [ eotests, *[ [ comp[eotest][kind][i] for eotest in eotests ] for i in range(4) ] ] ), kind )
        )
    fig, axs = pylab.subplots(1,1,figsize=(5,8),dpi=150,sharex=True)
    #pylab.errorbar( range(len(e2v["col1"])), y=e2v["col2"], yerr=e2v["col3"] )
    x =numpy.arange(len(tablearr[0][0]["col1"]))
    for i, k  in enumerate(tablearr):
        data, label = k
        pylab.errorbar( y=x+0*(1+0.001)*i,
                   x=(data["col3"]-data["col1"])/data["col1"],
                   xerr=numpy.sqrt(data["col2"]**2+data["col4"]**2)/data["col1"],
                   fmt="o", label=label
                  )
    pylab.yticks(x,data["col0"])
    pylab.xlabel("({}-{})/{}".format(botrun1,botrun2,botrun1))
    pylab.axvline(0,color="k")
    pylab.grid()
    pylab.xlim(-1.0,0.3)
    pylab.legend()
    pylab.show()
    pylab.savefig("plot{}-{}.png".format(botrun1,botrun2))
    
def run(botrun1, botrun2):
    eotests, comp=compare2runs( botrun1, botrun2 )
    makeplot( botrun1, botrun2, eotests, comp )