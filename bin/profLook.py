#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import numpy

from astropy.io import fits

from fitPlateProfile import DuPontProfile

basePath = os.getenv("LCOQA_DATA")
basePath = "/uufs/chpc.utah.edu/common/home/u0449727/lcoqa_data.saved"
class ConorQA(object):
    def __init__(self, plateID, fscanMJD, fscanID, expMJD, expNum, xPos, yPos, signal, signalModel):
        self.plateID = plateID
        self.fscanMJD = fscanMJD
        self.fscanID = fscanID
        self.expMJD = expMJD
        self.expNum = expNum
        self.xPos = xPos
        self.yPos = yPos
        self.signal = signal
        self.signalModel = signalModel
        self.prof = DuPontProfile()
        self.prof.getProfileFromDB(plateID, fscanID, fscanMJD)
        self.zeroPoint = -2.5*(numpy.log10(self.signal/self.signalModel))
        self.zeroPointRMS = self.rms(self.zeroPoint)
        self.profErrArray = numpy.array([self.prof.getErr(x,y) for x,y in zip(self.xPos,self.yPos)])
        self.profErrRMS = self.rms(self.profErrArray)


    def rms(self, array):
        # remove nans from array
        array = array[numpy.logical_not(numpy.isnan(array))]
        return numpy.sqrt(numpy.sum(array**2)/len(array))




def getSignal(mjd, expNo):
    sigPath = os.path.join(basePath, "lco", "%i"%mjd, "signal-%i-%i.fits"%(mjd, expNo))
    if not os.path.exists(sigPath):
        print("doesn't exist", sigPath)
        return None
    sigTable = fits.open(sigPath)[-1].data
    return sigTable


if __name__ == "__main__":

    onlyGood = True
    summaryTable = fits.open(os.path.join(basePath, "lco/apogee-summary.fits"))[-1].data
    cQAList = []
    for row in summaryTable:
        if onlyGood and not row["good_conditions"]:
            print("skipping %s, bad conditions"%row["name"])
            continue
        print("scraping up %s"%row["name"])
        plateID, fscanMJD, fscanID = [int(x) for x in row["name"].split("-")]
        expNum = row["expno"]
        expMJD = row["mjd"]
        sigTable = getSignal(expMJD, expNum)
        if sigTable is None:
            print("no signal table for %s, skipping"%row["name"])
            continue
        xPos = sigTable["xFocal"]
        yPos = sigTable["yFocal"]
        signal = sigTable["signal"]
        signal_model = sigTable["signal_model"]
        try:
           cQA = ConorQA(plateID, fscanMJD, fscanID, expMJD, expNum, xPos, yPos, signal, signal_model)
        except RuntimeError:
            print("failed to get profile for %s, skipping"%row["name"])
            continue
        cQAList.append(cQA)
    cQAList.sort(key=lambda x: x.zeroPointRMS)
    nNaNs = 0
    with open("profRMS.txt", "w") as f:
        f.write("PlateID ExpMJD ExpNum FscanMJD FscanID ZeroPointRMS ProfErrRMS ProfPercent\n")
        for cQA in cQAList:
            if numpy.isnan(cQA.profErrRMS):
                nNaNs += 1
            f.write("%i %i %i %i %i %.4f %.4f %.4f\n"%(cQA.plateID, cQA.expMJD, cQA.expNum, cQA.fscanMJD, cQA.fscanID, cQA.zeroPointRMS, cQA.profErrRMS, cQA.prof.percentInSpec))
    print("percent nans: %.2f"%(nNaNs/float(len(cQAList))*100))
    #with open("allHoles.txt", "w") as f:
        #f.write("signal profErr\n")
        #for cQA in cQAList:
            #for sig, prof in zip(cQA.zeroPoint, cQA.profErrArray):
                #f.write("%.4f %.4f\n"%(sig, prof))
    
           



