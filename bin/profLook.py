#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os

from astropy.io import fits

from fitPlateProfile import DuPontProfile

basePath = os.getenv("LCOQA_DATA")

class ConorQA(object):
    def __init__(self, plateID, fscanMJD, fscanID, xPos, yPos, signal, signalModel):
        self.plateID = plateID
        self.fscanMJD = fscanMJD
        self.fscanID = fscanID
        self.xPos = xPos
        self.yPos = yPos
        self.signal = signal
        self.signalModel = signalModel
        self.prof = DuPontProfile()
        self.prof.getProfileFromDB(plateID, fscanID, fscanMJD)


def getSignal(mjd, expNo):
    sigPath = os.path.join(basePath, "lco" "%i"%mjd, "signal-%i-%i.fits"%(mjd, expNo))
    if not os.path.exists(sigPath):
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
        cQAList.append(ConorQA(plateID, fscanMJD, fscanID, xPos, yPos, signal, signal_model))
