import ROOT
import numpy as np

def myHeaviside(x, xShift=0, yCentral=1):
    return yCentral if x == xShift else type(yCentral)(x > xShift)

def fExpAndHill(tfName, xmin, xmax, inverseHorizontal=False, inverseVerticle=False):
    # myHeaviside((x - [1])*signHorizontal)*[2]*exp(([1] - x)*signHorizontal/[4]) + myHeaviside(([1] - x)*signHorizontal, yCentral=0.)*[2]*(1 + (x - [1])*signHorizontal/[3])
    strXDisplacement = "(x-[1])"
    strHeavisideTerm = f"myHeaviside({strXDisplacement})"
    strHeavisideTermInversed = f"myHeaviside(-{strXDisplacement}, yCentral=0)"
    fStrHeavisideTerm = lambda is_positive: strHeavisideTerm if is_positive else strHeavisideTermInversed 
    # strSignH = "-" if inverseHorizontal else ""
    # strSignV = "-" if inverseHorizontal else ""
    fStrSign = lambda is_negative: "-" if is_negative else ""
    return ROOT.TF1(tfName, 
               f"{fStrSign(inverseVerticle)}[2]*(exp({fStrSign(not inverseHorizontal)}{strXDisplacement}/[4])*{fStrHeavisideTerm(not inverseHorizontal)}+(1 + {fStrSign(inverseHorizontal)}{strXDisplacement}/[3]))", 
               xmin, xmax)
    
tfile = ROOT.TFile("output_proper_Mx2-150_Mx1-1_20191224.root")

tfBetatauLab = fExpAndHill("tfBetatauLab", 0, 2)
tfChi2Beta = fExpAndHill("tfChi2Beta", 0, 1, inverseHorizontal=True)
tfChi2Gamma = fExpAndHill("tfChi2Gamma", 1, 20)
tfCtauLab = fExpAndHill("tfCtauLab", 0, 5)
tfCtauProper = fExpAndHill("tfCtauLab", 0, 0.6)

