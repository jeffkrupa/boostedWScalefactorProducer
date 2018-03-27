import os
import glob
import math
import ROOT
import pdb

from ROOT import *

#ROOT.gROOT.Macro("rootlogon.C")

import FWCore.ParameterSet.Config as cms
import miniTreeFunctions
from miniTreeFunctions import *

import miniTreeProducer_WTaggingSF_TightMu
from miniTreeProducer_WTaggingSF_TightMu import *

import sys
from DataFormats.FWLite import Events, Handle

from array import *

#ROOT.gSystem.Load('libCondFormatsBTagObjects') 
#command line options 
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-f", "--pathIn", dest="inputFile",
                  help="inputFile path")

parser.add_option("-o", "--outName", dest="outName",
                  help="output file name")

parser.add_option("-i", "--min", dest="min", 
		  help="input index low end")

parser.add_option("-j", "--max", dest="max", 
		  help="input index high end")

parser.add_option("-l", "--file", dest="txt", 
		  help="input txt file")

parser.add_option("-m", "--isMC", dest="isMC", 
		  help="bool for is MC")

parser.add_option("-y", "--saveTrig", dest="saveTrig",
                  help="bool to not save triggers for background")

parser.add_option("-x", "--xsec", dest="xsec", 
		  help="cross section")

parser.add_option("-S", "--syst", dest="syst",
                  help="Systematic")

parser.add_option("-H", "--hnevents", dest="h_nevents",
                  help="High number of events")

parser.add_option("-L", "--lnevents", dest="l_nevents",
                  help="Low number of events")


(options, args) = parser.parse_args()

inputfile = options.txt 

ff_n = 1000

Lnevents = int(options.l_nevents)
Hnevents = int(options.h_nevents)

num1 = int(options.min)
num2 = int(options.max)

d1 = options.outName 
d2 = '_'
print(options.outName)
outputfilename = d1 + d2 + options.min + '.root'

print outputfilename

import copy

f_h2ddt = ROOT.TFile.Open("h3_n2ddt_26eff_36binrho5pt_Spring16_pt200.root","read")
trans_h2ddt = f_h2ddt.Get("h2ddt")
trans_h2ddt.SetDirectory(0)
f_h2ddt.Close()

f_h2ddt2 = ROOT.TFile.Open("GridOutput_v13.root","read")
trans_h2ddt2 = f_h2ddt2.Get("Rho2D")
trans_h2ddt2.SetDirectory(0)
f_h2ddt2.Close()

f_h2ddt_photon = ROOT.TFile.Open("DDTmaps.root","read")
trans_h2ddt_photon = f_h2ddt_photon.Get("SmooDDT")
trans_h2ddt_photon.SetDirectory(0)
f_h2ddt_photon.Close()

#f_trig = ROOT.TFile.Open("RUNTriggerEfficiencies_SingleMuon_Run2016_V2p1_v03.root","read")
#trig_denom = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
#trig_numer = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
#trig_denom.SetDirectory(0)
#trig_numer.SetDirectory(0)
#trig_denom.RebinX(2)
#trig_numer.RebinX(2)
#trig_denom.RebinY(5)
#trig_numer.RebinY(5)        
#trig_eff = ROOT.TEfficiency()
#if (ROOT.TEfficiency.CheckConsistency(trig_numer,trig_denom)):
#  trig_eff = ROOT.TEfficiency(trig_numer,trig_denom)
#  trig_eff.SetDirectory(0)
#f_trig.Close()

f_pu= ROOT.TFile.Open("puWeights_All.root","read")
puw      = f_pu.Get("puw")
puw.SetDirectory(0)
puw_up   = f_pu.Get("puw_p")
puw_up.SetDirectory(0)
puw_down = f_pu.Get("puw_m")
puw_down.SetDirectory(0)
f_pu.Close()

lumi_GH = 16.146
lumi_BCDEF = 19.721
lumi_total = lumi_GH + lumi_BCDEF
        
f_mutrig_GH = ROOT.TFile.Open("EfficienciesAndSF_Period4.root","read")
mutrig_eff_GH = f_mutrig_GH.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
mutrig_eff_GH.Sumw2()
mutrig_eff_GH.SetDirectory(0)
f_mutrig_GH.Close()
        
f_mutrig_BCDEF = ROOT.TFile.Open("EfficienciesAndSF_RunBtoF.root","read")
mutrig_eff_BCDEF = f_mutrig_BCDEF.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
mutrig_eff_BCDEF.Sumw2()
mutrig_eff_BCDEF.SetDirectory(0)
f_mutrig_BCDEF.Close()

mutrig_eff = mutrig_eff_GH.Clone('pt_abseta_DATA_mutrig_ave')
mutrig_eff.Scale(lumi_GH/lumi_total)        
mutrig_eff.Add(mutrig_eff_BCDEF,lumi_BCDEF/lumi_total)

f_muid_GH = ROOT.TFile.Open("EfficienciesAndSF_GH.root","read")
muid_eff_GH = f_muid_GH.Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
muid_eff_GH.Sumw2()
muid_eff_GH.SetDirectory(0)
f_muid_GH.Close()
       
f_muid_BCDEF = ROOT.TFile.Open("EfficienciesAndSF_BCDEF.root","read")
muid_eff_BCDEF = f_muid_BCDEF.Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
muid_eff_BCDEF.Sumw2()
muid_eff_BCDEF.SetDirectory(0)
f_muid_BCDEF.Close()
        
muid_eff = muid_eff_GH.Clone('pt_abseta_DATA_muid_ave')
muid_eff.Scale(lumi_GH/lumi_total)        
muid_eff.Add(muid_eff_BCDEF,lumi_BCDEF/lumi_total)



f_muiso_GH = ROOT.TFile.Open("EfficienciesAndSF_ISO_GH.root","read")
muiso_eff_GH = f_muiso_GH.Get("TightISO_TightID_pt_eta/efficienciesDATA/pt_abseta_DATA")
muiso_eff_GH.Sumw2()
muiso_eff_GH.SetDirectory(0)
f_muiso_GH.Close()
        
f_muiso_BCDEF = ROOT.TFile.Open("EfficienciesAndSF_ISO_BCDEF.root","read")
muiso_eff_BCDEF = f_muiso_BCDEF.Get("TightISO_TightID_pt_eta/efficienciesDATA/pt_abseta_DATA")
muiso_eff_BCDEF.Sumw2()
muiso_eff_BCDEF.SetDirectory(0)
f_muiso_BCDEF.Close()
        
muiso_eff = muiso_eff_GH.Clone('pt_abseta_DATA_muiso_ave')
muiso_eff.Scale(lumi_GH/lumi_total)        
muiso_eff.Add(muiso_eff_BCDEF,lumi_BCDEF/lumi_total)

print sys.argv[1]

f =  ROOT.TFile(outputfilename, 'recreate')

f.cd()

myTree =  ROOT.TTree('myTree', 'myTree')

test = miniTreeProducer(options.isMC, options.saveTrig, options.syst, myTree, options.xsec, Lnevents, Hnevents)
test.runProducer(options.inputFile, inputfile, num1, num2, trans_h2ddt,trans_h2ddt2,trans_h2ddt_photon,puw,puw_up,puw_down,mutrig_eff,muid_eff,muiso_eff)

f.cd()
f.Write()
f.Close()



