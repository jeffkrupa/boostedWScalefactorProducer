import os
import glob
import math
import ROOT
from ROOT import *


#ROOT.gROOT.Macro("rootlogon.C")

import FWCore.ParameterSet.Config as cms

import sys
from DataFormats.FWLite import Events, Handle

from array import *

import copy

#defining functions
def div_except(a, b):
    if b>0:
        return float(a)/b
    else:
        return 1 
def btagging_efficiency_medium(pt):
   result = 0.898 + 0.0002254*pt -1.74e-6*pt*pt +2.71e-9*pt*pt*pt -1.39e-12*pt*pt*pt*pt
   return result	

def trigger_function(histo_efficiency,htJet30=700):
    result = histo_efficiency.GetBinContent(htJet30)
    return result

def ClosestJet(jets, fourvec): #returns the index of the jet (from a collection "jets") closest to the given four-vector
	DR = 9999.
	index = -1
	for j in range(0,len(jets)):
	    if jets[j].Pt() > 0 :
		dR = fourvec.DeltaR(jets[j])
		if dR < DR :
			DR = fourvec.DeltaR(jets[j])
			index = j
	return index

def MatchCollection(Col, jet): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
		C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def MatchCollection2(Col, jet, index): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
            if i != index:
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def MatchCollection3(Col, jet, index1, index2): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
            if i != index1 and i != index2:
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def MatchCollection4(Col, jet, index1, index2, index3): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
            if i != index1 and i != index2 and i != index3:
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def open_files(file_name,path) : #opens files to run on
    ff_n = 1000
    g = open(file_name)
    list_file = []
    final_list = []
    for i in range(ff_n):  # this is the length of the file
        list_file.append(g.readline().split())
    s = path

    for i in range(len(list_file)):
        for j in range(len(list_file[i])) :
            final_list.append(s + list_file[i][j])
    print final_list
    return final_list


def deltaR( particle, jet ) : #gives deltaR between two particles
    DeltaPhiHere = math.fabs( particle.phi() - jet.phi() )
    if DeltaPhiHere > math.pi :
        DeltaPhiHere = 2*math.pi - DeltaPhiHere

    deltaRHere = math.math.sqrt( (particle.eta() - jet.eta() )**2 + ( DeltaPhiHere )**2  )
    return deltaRHere

def getPUPPIweight( puppipt, puppieta ):
    PuppiWeightFile = ROOT.TFile.Open("puppiCorr.root","R")
    
    puppisd_corrGEN      = PuppiWeightFile.Get("puppiJECcorr_gen")
    puppisd_corrRECO_cen = PuppiWeightFile.Get("puppiJECcorr_reco_0eta1v3")
    puppisd_corrRECO_for = PuppiWeightFile.Get("puppiJECcorr_reco_1v3eta2v5")

    genCorr  = 1.
    recoCorr = 1.
    totalWeight = 1.

    genCorr =  puppisd_corrGEN.Eval( puppipt )

    if math.fabs(puppieta) <= 1.3:
      recoCorr = puppisd_corrRECO_cen.Eval( puppipt )
    else:
      recoCorr = puppisd_corrRECO_for.Eval( puppipt )

    totalWeight = genCorr * recoCorr
    PuppiWeightFile.Close()

    return totalWeight


class miniTreeProducer:
    def __init__(self, isMC, saveTrig, syst, tree, xsec, lnevents, hnevents):
        self.isMC = isMC
        self.saveTrig = saveTrig
        self.theTree = tree
        self.syst = syst
        self.xsecs = xsec
	self.lnevents = lnevents
	self.hnevents = hnevents

    def runProducer(self,location,inputfile, num1, num2,trans_h2ddt,trans_h2ddt2,trans_h2ddt_photon,puw,puw_up,puw_down,mutrig_eff,muid_eff,muiso_eff):

        self.AK8Puppijet0_tau21 = array('f', [-100.0])
        self.AK8Puppijet0_N2DDT_5Per = array('f', [-100.0])
        self.AK8Puppijet0_N2DDT_26Per = array('f', [-100.0])
        self.AK8Puppijet0_N2DDT_Photon = array('f', [-100.0])
        self.AK8Puppijet0_tau21DDT = array('f', [-100.0])
        self.AK8Puppijet0_N2 = array('f', [-100.0])
        self.AK8Puppijet0_doublecsv = array('f', [-100.0])
	self.AK8Puppijet0_pt = array('f', [-100.0])
	self.AK8CHSjet0_pt = array('f', [-100.0])
        self.AK8CHSjet0_eta = array('f', [-100.0])
        self.AK8CHSjet0_phi = array('f', [-100.0])
	self.AK8Puppijet0_msd = array('f', [-100.0])
        self.AK8Puppijet0_TheaCorr = array('f', [-100.0])
        self.AK8Puppijet0_msd_TheaCorr = array('f', [-100.0])
	self.AK8Puppijet0_vMatching = array('f', [-100.0])
	self.AK8Puppijet0_isHadronicV = array('f', [-100.0])
	self.pfmet = array('f', [-100.0])
	self.AK4Puppijet0_pt = array('f', [-100.0])
        self.AK4Puppijet0_eta = array('f', [-100.0])
        self.AK4Puppijet0_phi = array('f', [-100.0])
	self.AK4Puppijet0_csv = array('f', [-100.0])
	self.AK4Puppijet0_mass = array('f', [-100.0])
        self.AK4Puppijet1_pt = array('f', [-100.0])
        self.AK4Puppijet1_eta = array('f', [-100.0])
        self.AK4Puppijet1_phi = array('f', [-100.0])
        self.AK4Puppijet1_csv = array('f', [-100.0])
        self.AK4Puppijet1_mass = array('f', [-100.0])
        self.AK4Puppijet2_pt = array('f', [-100.0])
        self.AK4Puppijet2_eta = array('f', [-100.0])
        self.AK4Puppijet2_phi = array('f', [-100.0])
        self.AK4Puppijet2_csv = array('f', [-100.0])
        self.AK4Puppijet2_mass = array('f', [-100.0])
        self.AK4Puppijet3_pt = array('f', [-100.0])
        self.AK4Puppijet3_eta = array('f', [-100.0])
        self.AK4Puppijet3_phi = array('f', [-100.0])
        self.AK4Puppijet3_csv = array('f', [-100.0])
        self.AK4Puppijet3_mass = array('f', [-100.0])
	self.vmuoLoose0_pt = array('f', [-100.0])
        self.vmuoLoose0_eta = array('f', [-100.0])
        self.vmuoLoose0_phi = array('f', [-100.0])
        self.veleLoose0_pt = array('f', [-100.0])
        self.veleLoose0_eta = array('f', [-100.0])
        self.veleLoose0_phi = array('f', [-100.0])
	self.triggerBits = array('f', [-100.0])
        self.nLooseMu = array('f', [-100.0])
        self.nTightMu = array('f', [-100.0])
        self.nLooseEl = array('f', [-100.0])
        self.nTightEl = array('f', [-100.0])
	self.nEvents = array('f', [-100.0])
        self.xsec = array('f', [-100.0])
        self.LeadingAK8Jet_MatchedHadW = array('f', [-100.0])
	self.puWeight = array('f', [-100.0])
        self.triggerpassbb = array('f', [-100.0])
	self.scale1fb = array('f',[-100.0])
	self.weight = array('f',[-100.0])
	self.mutrigweight = array('f',[-100.0])
        self.mutrigweightDown = array('f',[-100.0])
        self.mutrigweightUp = array('f',[-100.0])
	self.puweight = array('f',[-100.0])
        self.puweight_up = array('f',[-100.0])
        self.puweight_down = array('f',[-100.0])
	self.muidweight = array('f',[-100.0])
        self.muidweightDown = array('f',[-100.0])
        self.muidweightUp = array('f',[-100.0])
        self.muisoweight = array('f',[-100.0])
        self.muisoweightDown = array('f',[-100.0])
        self.muisoweightUp = array('f',[-100.0])
	self.topPtWeight = array('f',[-100.0])

        #creating the tree branches we need
        self.theTree.Branch('AK8Puppijet0_tau21', self.AK8Puppijet0_tau21, 'AK8Puppijet0_tau21/F')
        self.theTree.Branch('AK8Puppijet0_N2DDT_5Per', self.AK8Puppijet0_N2DDT_5Per, 'AK8Puppijet0_N2DDT_5Per/F')
        self.theTree.Branch('AK8Puppijet0_N2DDT_26Per', self.AK8Puppijet0_N2DDT_26Per, 'AK8Puppijet0_N2DDT_26Per/F')
        self.theTree.Branch('AK8Puppijet0_N2DDT_Photon', self.AK8Puppijet0_N2DDT_Photon, 'AK8Puppijet0_N2DDT_Photon/F')
        self.theTree.Branch('AK8Puppijet0_tau21DDT', self.AK8Puppijet0_tau21DDT, 'AK8Puppijet0_tau21DDT/F')
        self.theTree.Branch('AK8Puppijet0_N2', self.AK8Puppijet0_N2, 'AK8Puppijet0_N2/F')
        self.theTree.Branch('AK8Puppijet0_doublecsv', self.AK8Puppijet0_doublecsv, 'AK8Puppijet0_doublecsv/F')
        self.theTree.Branch('AK8Puppijet0_pt', self.AK8Puppijet0_pt, 'AK8Puppijet0_pt/F')
        self.theTree.Branch('AK8CHSjet0_pt', self.AK8CHSjet0_pt, 'AK8CHSjet0_pt/F')
        self.theTree.Branch('AK8CHSjet0_eta', self.AK8CHSjet0_eta, 'AK8CHSjet0_eta/F')
        self.theTree.Branch('AK8CHSjet0_phi', self.AK8CHSjet0_phi, 'AK8CHSjet0_phi/F')
        self.theTree.Branch('AK8Puppijet0_msd', self.AK8Puppijet0_msd, 'AK8Puppijet0_msd/F')
        self.theTree.Branch('AK8Puppijet0_TheaCorr', self.AK8Puppijet0_TheaCorr, 'AK8Puppijet0_TheaCorr/F')
        self.theTree.Branch('AK8Puppijet0_msd_TheaCorr', self.AK8Puppijet0_msd_TheaCorr, 'AK8Puppijet0_msd_TheaCorr/F')
        self.theTree.Branch('AK8Puppijet0_vMatching', self.AK8Puppijet0_vMatching, 'AK8Puppijet0_vMatching/F')
        self.theTree.Branch('AK8Puppijet0_isHadronicV', self.AK8Puppijet0_isHadronicV, 'AK8Puppijet0_isHadronicV/F')
        self.theTree.Branch('nLooseMu', self.nLooseMu, 'nLooseMu/F')
        self.theTree.Branch('nTightMu', self.nTightMu, 'nTightMu/F')
        self.theTree.Branch('nLooseEl', self.nLooseEl, 'nLooseEl/F')
        self.theTree.Branch('nTightEl', self.nTightEl, 'nTightEl/F')
        self.theTree.Branch('nEvents', self.nEvents, 'nEvents/F')
        self.theTree.Branch('xsec', self.xsec, 'xsec/F')
        self.theTree.Branch('LeadingAK8Jet_MatchedHadW', self.LeadingAK8Jet_MatchedHadW, 'LeadingAK8Jet_MatchedHadW/F')
        self.theTree.Branch('puWeight', self.puWeight, 'puWeight/F')
        self.theTree.Branch('scale1fb', self.scale1fb, 'scale1fb/F')
	self.theTree.Branch('triggerBits', self.triggerBits, 'triggerBits/F')
        self.theTree.Branch('weight', self.weight, 'weight/F')
        self.theTree.Branch('mutrigweight', self.mutrigweight, 'mutrigweight/F')
        self.theTree.Branch('mutrigweightUp', self.mutrigweightUp, 'mutrigweightUp/F')
        self.theTree.Branch('mutrigweightDown', self.mutrigweightDown, 'mutrigweightDown/F')
        self.theTree.Branch('puweight', self.puweight, 'puweight/F')
        self.theTree.Branch('puweight_up', self.puweight_up, 'puweight_up/F')
        self.theTree.Branch('puweight_down', self.puweight_down, 'puweight_down/F')
        self.theTree.Branch('muidweight', self.muidweight, 'muidweight/F')
        self.theTree.Branch('muidweightUp', self.muidweightUp, 'muidweightUp/F')
        self.theTree.Branch('muidweightDown', self.muidweightDown, 'muidweightDown/F')
        self.theTree.Branch('muisoweight', self.muisoweight, 'muisoweight/F')
        self.theTree.Branch('muisoweightUp', self.muisoweightUp, 'muisoweightUp/F')
        self.theTree.Branch('muisoweightDown', self.muisoweightDown, 'muisoweightDown/F')
        self.theTree.Branch('topPtWeight', self.topPtWeight, 'topPtWeight/F')

        self.bb = ROOT.TH1F("bb", "No Cuts", 3, -0.5, 1.5)
        self.bb0 = ROOT.TH1F("bb0", "After MET", 3, -0.5, 1.5)
        self.bb1 = ROOT.TH1F("bb1", "After Tight Muon", 3, -0.5, 1.5)
        self.bb2 = ROOT.TH1F("bb2", "After Loose Muon", 3, -0.5, 1.5)
        self.bb3 = ROOT.TH1F("bb3", "After Hadronic AK8 Jet", 3, -0.5, 1.5)
        self.bb4 = ROOT.TH1F("bb4", "After b-tagged AK4 Jet", 3, -0.5, 1.5)
        self.bb5 = ROOT.TH1F("bb5", "After Leptonic AK8  Jet", 3, -0.5, 1.5)
        self.bb6 = ROOT.TH1F("bb6", "After Leptonic AK8  Jet", 3, -0.5, 1.5)
        self.bb7 = ROOT.TH1F("bb7", "After Leptonic AK8  Jet", 3, -0.5, 1.5)
        self.bb8 = ROOT.TH1F("bb8", "After Leptonic AK8  Jet", 3, -0.5, 1.5)

        self.Files_list	= open_files( inputfile, location )

        self.count = 0
        #loop over files
        for i in range(num1, num2):
            self.files = self.Files_list[i]
            print self.files
            self.f1 = ROOT.TFile.Open(self.files, "READ")
            self.treeMine  = self.f1.Get('otree')
            self.nevent = self.treeMine.GetEntries();
            self.nFills = 0
            self.nFills2 = 0

	    print "low # of events ", self.lnevents
	    print "high # of events ", self.hnevents
            #loop over events in file
            print "Start looping"
	    self.count = self.lnevents
            for j in range(self.lnevents,self.nevent):
#            for j in range(0,20):
                self.treeMine.GetEntry(j)
                self.count = self.count + 1
                if self.count % 1000 == 0 :
                    print "processing events", self.count
		if self.count == self.hnevents:
		   break;
	
                #variables we need from the baconbits ntuple
     		self.AK8Puppijet0_tau21[0] = self.treeMine.AK8Puppijet0_tau21
		self.AK8Puppijet0_N2[0] = self.treeMine.AK8Puppijet0_N2sdb1
                self.AK8Puppijet0_doublecsv[0] = self.treeMine.AK8Puppijet0_doublecsv
		self.AK8Puppijet0_pt[0] = self.treeMine.AK8Puppijet0_pt
		self.AK8CHSjet0_pt[0] = self.treeMine.AK8CHSjet0_pt
                self.AK8CHSjet0_eta[0] = self.treeMine.AK8CHSjet0_eta
                self.AK8CHSjet0_phi[0] = self.treeMine.AK8CHSjet0_phi
		self.pfmet[0] = self.treeMine.pfmet
		self.AK8Puppijet0_vMatching[0] = self.treeMine.AK8Puppijet0_vMatching
		self.AK8Puppijet0_isHadronicV[0] = self.treeMine.AK8Puppijet0_isHadronicV
		self.AK4Puppijet0_pt[0] = self.treeMine.AK4Puppijet0_pt
                self.AK4Puppijet0_eta[0] = self.treeMine.AK4Puppijet0_eta
                self.AK4Puppijet0_phi[0] = self.treeMine.AK4Puppijet0_phi
                self.AK4Puppijet0_csv[0] = self.treeMine.AK4Puppijet0_csv
                self.AK4Puppijet1_pt[0] = self.treeMine.AK4Puppijet1_pt
                self.AK4Puppijet1_eta[0] = self.treeMine.AK4Puppijet1_eta
                self.AK4Puppijet1_phi[0] = self.treeMine.AK4Puppijet1_phi
                self.AK4Puppijet1_csv[0] = self.treeMine.AK4Puppijet1_csv
                self.AK4Puppijet2_pt[0] = self.treeMine.AK4Puppijet2_pt
                self.AK4Puppijet2_eta[0] = self.treeMine.AK4Puppijet2_eta
                self.AK4Puppijet2_phi[0] = self.treeMine.AK4Puppijet2_phi
                self.AK4Puppijet2_csv[0] = self.treeMine.AK4Puppijet2_csv
                self.AK4Puppijet3_pt[0] = self.treeMine.AK4Puppijet3_pt
                self.AK4Puppijet3_eta[0] = self.treeMine.AK4Puppijet3_eta
                self.AK4Puppijet3_phi[0] = self.treeMine.AK4Puppijet3_phi
                self.AK4Puppijet3_csv[0] = self.treeMine.AK4Puppijet3_csv
		self.triggerBits[0] = self.treeMine.triggerBits&4
		self.AK8Puppijet0_msd[0] = self.treeMine.AK8Puppijet0_msd
		self.AK8Puppijet0_TheaCorr[0] =  getPUPPIweight(self.treeMine.AK8Puppijet0_pt, self.treeMine.AK8Puppijet0_eta)
		self.AK8Puppijet0_msd_TheaCorr[0] = self.AK8Puppijet0_msd[0]*self.AK8Puppijet0_TheaCorr[0]
		self.nEvents[0] = self.f1.NEvents.GetBinContent(1)
		self.xsec[0] = float(self.xsecs)
		self.puWeight[0] = self.treeMine.puWeight
		self.scale1fb[0] = self.treeMine.scale1fb


		if self.AK8Puppijet0_msd_TheaCorr[0] > 0:
		  self.AK8Puppijet0_tau21DDT[0] = self.treeMine.AK8Puppijet0_tau21+0.063*math.log(self.AK8Puppijet0_msd_TheaCorr[0]*self.AK8Puppijet0_msd_TheaCorr[0]/self.treeMine.AK8Puppijet0_pt)
		else:
		  self.AK8Puppijet0_tau21DDT[0] = -9999

		#Muon Trigger Weight
		mutrigweight = 1
                mutrigweightDown = 1
                mutrigweightUp = 1
                if self.treeMine.nmuLoose> 0:
                  muPtForTrig = max(52.,min(self.treeMine.vmuoLoose0_pt,700.))
                  muEtaForTrig = min(abs(self.treeMine.vmuoLoose0_eta),2.3)
                  mutrigweight = mutrig_eff.GetBinContent(mutrig_eff.FindBin(muPtForTrig,muEtaForTrig))
                  mutrigweightUp = mutrigweight + mutrig_eff.GetBinError(mutrig_eff.FindBin(muPtForTrig,muEtaForTrig))
                  mutrigweightDown = mutrigweight - mutrig_eff.GetBinError(mutrig_eff.FindBin(muPtForTrig,muEtaForTrig))
                  if mutrigweight<=0 or mutrigweightDown<=0 or mutrigweightUp<=0:
                      print 'mutrigweights are %f, %f, %f, setting all to 1'%(mutrigweight,mutrigweightUp,mutrigweightDown)
                      mutrigweight = 1
                      mutrigweightDown = 1
                      mutrigweightUp = 1
		self.mutrigweight[0] = mutrigweight
                self.mutrigweightDown[0] = mutrigweightDown
                self.mutrigweightUp[0] = mutrigweightUp

		#puweight
                nPuForWeight = min(self.treeMine.npu,49.5)
                puweight = puw.GetBinContent(puw.FindBin(nPuForWeight))
                puweight_up = puw_up.GetBinContent(puw_up.FindBin(nPuForWeight))
		puweight_down = puw_down.GetBinContent(puw_down.FindBin(nPuForWeight))
		self.puweight[0] = puweight
                self.puweight_up[0] = puweight_up
                self.puweight_down[0] = puweight_down

		#muIsoWeight
                muidweight = 1
                muidweightDown = 1
                muidweightUp = 1
                if self.treeMine.nmuLoose > 0:
                  muPtForId = max(20.,min(self.treeMine.vmuoLoose0_pt,100.))
                  muEtaForId = min(abs(self.treeMine.vmuoLoose0_eta),2.3)
                  muidweight = muid_eff.GetBinContent(muid_eff.FindBin(muPtForId,muEtaForId))
                  muidweightUp = muidweight + muid_eff.GetBinError(muid_eff.FindBin(muPtForId,muEtaForId))
                  muidweightDown = muidweight - muid_eff.GetBinError(muid_eff.FindBin(muPtForId,muEtaForId))
                  if muidweight<=0 or muidweightDown<=0 or muidweightUp<=0:
                     print 'muidweights are %f, %f, %f, setting all to 1'%(muidweight,muidweightUp,muidweightDown)
                     muidweight = 1
                     muidweightDown = 1
                     muidweightUp = 1

		self.muidweight[0] = muidweight
                self.muidweightDown[0] = muidweightDown
                self.muidweightUp[0] = muidweightUp

		#muIsoWeight
                muisoweight = 1
                muisoweightDown = 1
                muisoweightUp = 1
                if self.treeMine.nmuLoose > 0:
                  muPtForIso = max(20.,min(self.treeMine.vmuoLoose0_pt,100.))
                  muEtaForIso = min(abs(self.treeMine.vmuoLoose0_eta),2.3)
                  muisoweight = muiso_eff.GetBinContent(muiso_eff.FindBin(muPtForIso,muEtaForIso))
                  muisoweightUp = muisoweight + muiso_eff.GetBinError(muiso_eff.FindBin(muPtForIso,muEtaForIso))
                  muisoweightDown = muisoweight - muiso_eff.GetBinError(muiso_eff.FindBin(muPtForIso,muEtaForIso))
                  if muisoweight<=0 or muisoweightDown<=0 or muisoweightUp<=0:
                     print 'muisoweights are %f, %f, %f, setting all to 1'%(muisoweight,muisoweightUp,muisoweightDown)
                     muisoweight = 1
                     muisoweightDown = 1
		     muisoweightUp = 1
	
                self.muisoweight[0] = muisoweight
                self.muisoweightDown[0] = muisoweightDown
                self.muisoweightUp[0] = muisoweightUp

                if self.isMC == 'False':
		  self.weight[0] = 1
                if self.isMC == 'True':		
                  self.weight[0] = puweight*self.treeMine.scale1fb*self.treeMine.kfactor*mutrigweight


		#26% eff N2DDT
		jpt_8 = self.treeMine.AK8Puppijet0_pt
		jmsd_8 = self.AK8Puppijet0_msd_TheaCorr[0]
		if jmsd_8 <= 0: jmsd_8 = 0.01
		rh_8 = math.log(jmsd_8*jmsd_8/jpt_8/jpt_8)
		jtN2b1sd_8 = self.treeMine.AK8Puppijet0_N2sdb1

		# define N2ddt
		cur_rho_index = trans_h2ddt.GetXaxis().FindBin(rh_8)
		cur_pt_index = trans_h2ddt.GetYaxis().FindBin(jpt_8)
		if rh_8 > trans_h2ddt.GetXaxis().GetBinUpEdge(trans_h2ddt.GetXaxis().GetNbins()):    
		    cur_rho_index = trans_h2ddt.GetXaxis().GetNbins()
		if rh_8 < trans_h2ddt.GetXaxis().GetBinLowEdge(1):    
		    cur_rho_index = 1
		if jpt_8 > trans_h2ddt.GetYaxis().GetBinUpEdge(trans_h2ddt.GetYaxis().GetNbins()):    
		    cur_pt_index = trans_h2ddt.GetYaxis().GetNbins()
		if jpt_8 < trans_h2ddt.GetYaxis().GetBinLowEdge(1):    
		    cur_pt_index = 1
		jtN2b1sddt_8 = jtN2b1sd_8 - trans_h2ddt.GetBinContent(cur_rho_index, cur_pt_index)
		self.AK8Puppijet0_N2DDT_26Per[0] = jtN2b1sddt_8

                #5% eff N2DDT
                jpt_2 = self.treeMine.AK8Puppijet0_pt
                jmsd_2 = self.AK8Puppijet0_msd_TheaCorr[0]
                if jmsd_2 <= 0: jmsd_2 = 0.01
                rh_2 = math.log(jmsd_2*jmsd_2/jpt_2/jpt_2)
                jtN2b1sd_2 = self.treeMine.AK8Puppijet0_N2sdb1

                # define N2ddt
                cur_rho_index2 = trans_h2ddt2.GetXaxis().FindBin(rh_2)
                cur_pt_index2 = trans_h2ddt2.GetYaxis().FindBin(jpt_2)
                if rh_2 > trans_h2ddt2.GetXaxis().GetBinUpEdge(trans_h2ddt2.GetXaxis().GetNbins()):
                    cur_rho_index2 = trans_h2ddt2.GetXaxis().GetNbins()
                if rh_2 < trans_h2ddt2.GetXaxis().GetBinLowEdge(1):
                    cur_rho_index2 = 1
                if jpt_2 > trans_h2ddt2.GetYaxis().GetBinUpEdge(trans_h2ddt2.GetYaxis().GetNbins()):
                    cur_pt_index2 = trans_h2ddt2.GetYaxis().GetNbins()
                if jpt_2 < trans_h2ddt2.GetYaxis().GetBinLowEdge(1):
                    cur_pt_index2 = 1
                jtN2b1sddt_2 = jtN2b1sd_2 - trans_h2ddt2.GetBinContent(cur_rho_index2, cur_pt_index2)
                self.AK8Puppijet0_N2DDT_5Per[0] = jtN2b1sddt_2

                # define Photon N2ddt
                cur_rho_index_photon = trans_h2ddt_photon.GetXaxis().FindBin(rh_2)
                cur_pt_index_photon = trans_h2ddt_photon.GetYaxis().FindBin(jpt_2)
                if rh_2 > trans_h2ddt_photon.GetXaxis().GetBinUpEdge(trans_h2ddt_photon.GetXaxis().GetNbins()):
                    cur_rho_index_photon = trans_h2ddt_photon.GetXaxis().GetNbins()
                if rh_2 < trans_h2ddt_photon.GetXaxis().GetBinLowEdge(1):
                    cur_rho_index_photon = 1
                if jpt_2 > trans_h2ddt_photon.GetYaxis().GetBinUpEdge(trans_h2ddt_photon.GetYaxis().GetNbins()):
                    cur_pt_index_photon = trans_h2ddt_photon.GetYaxis().GetNbins()
                if jpt_2 < trans_h2ddt_photon.GetYaxis().GetBinLowEdge(1):
                    cur_pt_index_photon = 1
                jtN2b1sddt_photon = jtN2b1sd_2 - trans_h2ddt_photon.GetBinContent(cur_rho_index_photon, cur_pt_index_photon)
                self.AK8Puppijet0_N2DDT_Photon[0] = jtN2b1sddt_photon

                #Determining the number of loose Muons for WTagger SF
                self.nLooseMuons = 0

		if self.treeMine.nmuTight == 1:
		  if self.treeMine.vmuoLoose0_pt > 20 and abs(self.treeMine.vmuoLoose0_eta) < 2.4:
		    self.nLooseMuons+= 1

		if self.treeMine.nmuTight == 2:
		  if self.treeMine.vmuoLoose0_pt > 20 and abs(self.treeMine.vmuoLoose0_eta) < 2.4:
		    self.nLooseMuons+= 1
		  if self.treeMine.vmuoLoose1_pt > 20 and abs(self.treeMine.vmuoLoose1_eta) < 2.4:
                    self.nLooseMuons+= 1


		self.nLooseMu[0] = self.nLooseMuons

		#Determining the number of Tight Muons for WTaggerSF
		self.nTightMuons = 0

                if self.treeMine.nmuTight == 1:
		  if self.treeMine.vmuoLoose0_pt > 53 and abs(self.treeMine.vmuoLoose0_eta) < 2.1:
		    self.nTightMuons+= 1
                    self.vmuoLoose0_pt[0] = self.treeMine.vmuoLoose0_pt
                    self.vmuoLoose0_eta[0] = self.treeMine.vmuoLoose0_eta
                    self.vmuoLoose0_phi[0] = self.treeMine.vmuoLoose0_phi

		if self.treeMine.nmuTight == 2:
		  if self.treeMine.vmuoLoose0_pt > 53 and abs(self.treeMine.vmuoLoose0_eta) < 2.1:
                    self.nTightMuons+= 1
                    self.vmuoLoose0_pt[0] = self.treeMine.vmuoLoose0_pt
                    self.vmuoLoose0_eta[0] = self.treeMine.vmuoLoose0_eta
                    self.vmuoLoose0_phi[0] = self.treeMine.vmuoLoose0_phi
		  if self.treeMine.vmuoLoose1_pt > 53 and abs(self.treeMine.vmuoLoose1_eta) < 2.1:
                    self.nTightMuons+= 1
                    self.vmuoLoose0_pt[0] = self.treeMine.vmuoLoose1_pt
                    self.vmuoLoose0_eta[0] = self.treeMine.vmuoLoose1_eta
                    self.vmuoLoose0_phi[0] = self.treeMine.vmuoLoose1_phi

		
                self.nTightMu[0] = self.nTightMuons

		#Determining the number of Loose Electrons - only use Muon decays in final calculation, unnecessary
		self.nLooseElectrons = 0


		if self.treeMine.neleHEEP == 1:
		  if self.treeMine.veleLoose0_pt > 35 and ((abs(self.treeMine.veleLoose0_eta) < 1.442) or (abs(self.treeMine.veleLoose0_eta) > 1.56 and abs(self.treeMine.veleLoose0_eta) < 2.5)):
		    self.nLooseElectrons+= 1

                if self.treeMine.neleHEEP == 2:
                  if self.treeMine.veleLoose0_pt > 35 and ((abs(self.treeMine.veleLoose0_eta) < 1.442) or (abs(self.treeMine.veleLoose0_eta) > 1.56 and abs(self.treeMine.veleLoose0_eta) < 2.5)):
                    self.nLooseElectrons+= 1
		if self.treeMine.veleLoose1_pt > 35 and ((abs(self.treeMine.veleLoose1_eta) < 1.442) or (abs(self.treeMine.veleLoose1_eta) > 1.56 and abs(self.treeMine.veleLoose1_eta) < 2.5)):
                    self.nLooseElectrons+= 1		

                self.nLooseEl[0] = self.nLooseElectrons

		#Determining the number of tight Electrons - only use Muon decays in final calculation, unnecessary
		self.nTightElectrons = 0


                if self.treeMine.neleHEEP == 1:
		  if self.treeMine.veleLoose0_pt > 120 and ((abs(self.treeMine.veleLoose0_eta) < 1.442) or (abs(self.treeMine.veleLoose0_eta) > 1.56 and abs(self.treeMine.veleLoose0_eta) < 2.5)):
                    self.nTightElectrons+= 1
                    self.veleLoose0_pt[0] = self.treeMine.veleLoose0_pt
                    self.veleLoose0_eta[0] = self.treeMine.veleLoose0_eta
                    self.veleLoose0_phi[0] = self.treeMine.veleLoose0_phi

                if self.treeMine.neleHEEP == 2:
                  if self.treeMine.veleLoose0_pt > 120 and ((abs(self.treeMine.veleLoose0_eta) < 1.442) or (abs(self.treeMine.veleLoose0_eta) > 1.56 and abs(self.treeMine.veleLoose0_eta) < 2.5)):
                    self.nTightElectrons+= 1
                    self.veleLoose0_pt[0] = self.treeMine.veleLoose0_pt
                    self.veleLoose0_eta[0] = self.treeMine.veleLoose0_eta
                    self.veleLoose0_phi[0] = self.treeMine.veleLoose0_phi
                if self.treeMine.veleLoose1_pt > 120 and ((abs(self.treeMine.veleLoose1_eta) < 1.442) or (abs(self.treeMine.veleLoose1_eta) > 1.56 and abs(self.treeMine.veleLoose1_eta) < 2.5)):
                    self.nTightElectrons+= 1
                    self.veleLoose0_pt[0] = self.treeMine.veleLoose1_pt
                    self.veleLoose0_eta[0] = self.treeMine.veleLoose1_eta
                    self.veleLoose0_phi[0] = self.treeMine.veleLoose1_phi


                self.nTightEl[0] = self.nTightElectrons
	

		#Checking whether the AK8 jet passes selection
		HadronicAK8Jet = 0.
		DPhi_AK8Jet_TightMuon = self.AK8CHSjet0_phi[0] - self.vmuoLoose0_phi[0]
		if DPhi_AK8Jet_TightMuon >= math.pi:
		  DPhi_AK8Jet_TightMuon -= 2.*math.pi
		elif DPhi_AK8Jet_TightMuon < -math.pi:
		  DPhi_AK8Jet_TightMuon += 2.*math.pi

		DR_AK8Jet_TightMuon = math.sqrt((self.AK8CHSjet0_eta[0] - self.vmuoLoose0_eta[0])*(self.AK8CHSjet0_eta[0] - self.vmuoLoose0_eta[0])+DPhi_AK8Jet_TightMuon*DPhi_AK8Jet_TightMuon)

                DPhi_AK8Jet_TightElectron = self.AK8CHSjet0_phi[0] - self.veleLoose0_phi[0]
                if DPhi_AK8Jet_TightElectron >= math.pi:
                  DPhi_AK8Jet_TightElectron -= 2.*math.pi
                elif DPhi_AK8Jet_TightElectron < -math.pi:
                  DPhi_AK8Jet_TightElectron += 2.*math.pi

                DR_AK8Jet_TightElectron = math.sqrt((self.AK8CHSjet0_eta[0] - self.veleLoose0_eta[0])*(self.AK8CHSjet0_eta[0] - self.veleLoose0_eta[0])+DPhi_AK8Jet_TightElectron*DPhi_AK8Jet_TightElectron)

		if abs(self.AK8CHSjet0_eta[0]) < 2.4 and DR_AK8Jet_TightMuon > 1.0:
		  HadronicAK8Jet = 1

		#Checking whether there is at least one b-tagged AK4 jet. Has to be isolated from the AK8 jet and the Muon. Update to new medium CSV WP???
		nBTaggedAK4Jet = 0
		DR_TightMuon_Matched = 0
		DR_AK8Jet_Matched = 0
                DPhi_AK4Jet_TightMuon = self.AK4Puppijet0_phi[0] - self.vmuoLoose0_phi[0]
                if DPhi_AK4Jet_TightMuon >= math.pi:
                  DPhi_AK4Jet_TightMuon -= 2.*math.pi
                elif DPhi_AK4Jet_TightMuon < -math.pi:
                  DPhi_AK4Jet_TightMuon += 2.*math.pi

                DR_AK4Jet_TightMuon = math.sqrt((self.AK4Puppijet0_eta[0] - self.vmuoLoose0_eta[0])*(self.AK4Puppijet0_eta[0] - self.vmuoLoose0_eta[0])+DPhi_AK4Jet_TightMuon*DPhi_AK4Jet_TightMuon)
		if DR_AK4Jet_TightMuon < 0.3:
		  DR_TightMuon_Matched = 1

                DPhi_AK4Jet_AK8Jet = self.AK4Puppijet0_phi[0] - self.AK8CHSjet0_phi[0]
                if DPhi_AK4Jet_AK8Jet >= math.pi:
                  DPhi_AK4Jet_AK8Jet -= 2.*math.pi
                elif DPhi_AK4Jet_AK8Jet < -math.pi:
                  DPhi_AK4Jet_AK8Jet += 2.*math.pi

                DR_AK4Jet_AK8Jet = math.sqrt((self.AK4Puppijet0_eta[0] - self.AK8CHSjet0_eta[0])*(self.AK4Puppijet0_eta[0] - self.AK8CHSjet0_eta[0])+DPhi_AK4Jet_AK8Jet*DPhi_AK4Jet_AK8Jet)
                if DR_AK4Jet_AK8Jet < 0.8:
                  DR_AK8Jet_Matched = 1


                if abs(self.AK4Puppijet0_eta[0]) < 2.4 and self.AK4Puppijet0_csv[0] > 0.8838 and DR_TightMuon_Matched == 0 and DR_AK8Jet_Matched == 0 and self.AK4Puppijet0_pt[0] > 30.:
		  nBTaggedAK4Jet += 1

                DR_TightMuon_Matched = 0
                DR_AK8Jet_Matched = 0
                DPhi_AK4Jet_TightMuon = self.AK4Puppijet1_phi[0] - self.vmuoLoose0_phi[0]
                if DPhi_AK4Jet_TightMuon >= math.pi:
                  DPhi_AK4Jet_TightMuon -= 2.*math.pi
                elif DPhi_AK4Jet_TightMuon < -math.pi:
                  DPhi_AK4Jet_TightMuon += 2.*math.pi

                DR_AK4Jet_TightMuon = math.sqrt((self.AK4Puppijet1_eta[0] - self.vmuoLoose0_eta[0])*(self.AK4Puppijet1_eta[0] - self.vmuoLoose0_eta[0])+DPhi_AK4Jet_TightMuon*DPhi_AK4Jet_TightMuon)
                if DR_AK4Jet_TightMuon < 0.3:
                  DR_TightMuon_Matched = 1

                DPhi_AK4Jet_AK8Jet = self.AK4Puppijet1_phi[0] - self.AK8CHSjet0_phi[0]
                if DPhi_AK4Jet_AK8Jet >= math.pi:
                  DPhi_AK4Jet_AK8Jet -= 2.*math.pi
                elif DPhi_AK4Jet_AK8Jet < -math.pi:
                  DPhi_AK4Jet_AK8Jet += 2.*math.pi

                DR_AK4Jet_AK8Jet = math.sqrt((self.AK4Puppijet1_eta[0] - self.AK8CHSjet0_eta[0])*(self.AK4Puppijet1_eta[0] - self.AK8CHSjet0_eta[0])+DPhi_AK4Jet_AK8Jet*DPhi_AK4Jet_AK8Jet)
                if DR_AK4Jet_AK8Jet < 0.8:
                  DR_AK8Jet_Matched = 1


                if abs(self.AK4Puppijet1_eta[0]) < 2.4 and self.AK4Puppijet1_csv[0] > 0.8838 and DR_TightMuon_Matched == 0 and DR_AK8Jet_Matched == 0 and self.AK4Puppijet1_pt[0] > 30:
                  nBTaggedAK4Jet += 1

                DR_TightMuon_Matched = 0
                DR_AK8Jet_Matched = 0
                DPhi_AK4Jet_TightMuon = self.AK4Puppijet2_phi[0] - self.vmuoLoose0_phi[0]
                if DPhi_AK4Jet_TightMuon >= math.pi:
                  DPhi_AK4Jet_TightMuon -= 2.*math.pi
                elif DPhi_AK4Jet_TightMuon < -math.pi:
                  DPhi_AK4Jet_TightMuon += 2.*math.pi

                DR_AK4Jet_TightMuon = math.sqrt((self.AK4Puppijet2_eta[0] - self.vmuoLoose0_eta[0])*(self.AK4Puppijet2_eta[0] - self.vmuoLoose0_eta[0])+DPhi_AK4Jet_TightMuon*DPhi_AK4Jet_TightMuon)
                if DR_AK4Jet_TightMuon < 0.3:
                  DR_TightMuon_Matched = 1

                DPhi_AK4Jet_AK8Jet = self.AK4Puppijet2_phi[0] - self.AK8CHSjet0_phi[0]
                if DPhi_AK4Jet_AK8Jet >= math.pi:
                  DPhi_AK4Jet_AK8Jet -= 2.*math.pi
                elif DPhi_AK4Jet_AK8Jet < -math.pi:
                  DPhi_AK4Jet_AK8Jet += 2.*math.pi

                DR_AK4Jet_AK8Jet = math.sqrt((self.AK4Puppijet2_eta[0] - self.AK8CHSjet0_eta[0])*(self.AK4Puppijet2_eta[0] - self.AK8CHSjet0_eta[0])+DPhi_AK4Jet_AK8Jet*DPhi_AK4Jet_AK8Jet)
                if DR_AK4Jet_AK8Jet < 0.8:
                  DR_AK8Jet_Matched = 1


                if abs(self.AK4Puppijet2_eta[0]) < 2.4 and self.AK4Puppijet2_csv[0] > 0.8838 and DR_TightMuon_Matched == 0 and DR_AK8Jet_Matched == 0 and self.AK4Puppijet2_pt[0] > 30.:
                  nBTaggedAK4Jet += 1

                DR_TightMuon_Matched = 0
                DR_AK8Jet_Matched = 0
                DPhi_AK4Jet_TightMuon = self.AK4Puppijet3_phi[0] - self.vmuoLoose0_phi[0]
                if DPhi_AK4Jet_TightMuon >= math.pi:
                  DPhi_AK4Jet_TightMuon -= 2.*math.pi
                elif DPhi_AK4Jet_TightMuon < -math.pi:
                  DPhi_AK4Jet_TightMuon += 2.*math.pi

                DR_AK4Jet_TightMuon = math.sqrt((self.AK4Puppijet3_eta[0] - self.vmuoLoose0_eta[0])*(self.AK4Puppijet3_eta[0] - self.vmuoLoose0_eta[0])+DPhi_AK4Jet_TightMuon*DPhi_AK4Jet_TightMuon)
                if DR_AK4Jet_TightMuon < 0.3:
                  DR_TightMuon_Matched = 1

                DPhi_AK4Jet_AK8Jet = self.AK4Puppijet3_phi[0] - self.AK8CHSjet0_phi[0]
                if DPhi_AK4Jet_AK8Jet >= math.pi:
                  DPhi_AK4Jet_AK8Jet -= 2.*math.pi
                elif DPhi_AK4Jet_AK8Jet < -math.pi:
                  DPhi_AK4Jet_AK8Jet += 2.*math.pi

                DR_AK4Jet_AK8Jet = math.sqrt((self.AK4Puppijet3_eta[0] - self.AK8CHSjet0_eta[0])*(self.AK4Puppijet3_eta[0] - self.AK8CHSjet0_eta[0])+DPhi_AK4Jet_AK8Jet*DPhi_AK4Jet_AK8Jet)
                if DR_AK4Jet_AK8Jet < 0.8:
                  DR_AK8Jet_Matched = 1


                if abs(self.AK4Puppijet3_eta[0]) < 2.4 and self.AK4Puppijet3_csv[0] > 0.8838 and DR_TightMuon_Matched == 0 and DR_AK8Jet_Matched == 0 and self.AK4Puppijet3_pt[0] > 30:
                  nBTaggedAK4Jet += 1

		#Leptonic W pT
		PassingLeptonicW = 0.
		self.pxMet = self.pfmet[0]*math.cos(self.treeMine.pfmetphi)
		self.pyMet = self.pfmet[0]*math.sin(self.treeMine.pfmetphi)
		self.pxMuon = self.vmuoLoose0_pt[0]*math.cos(self.vmuoLoose0_phi[0])
                self.pyMuon = self.vmuoLoose0_pt[0]*math.sin(self.vmuoLoose0_phi[0])
		if math.sqrt((self.pxMet+self.pxMuon)*(self.pxMet+self.pxMuon) + (self.pyMet+self.pyMuon)*(self.pyMet+self.pyMuon)) > 200.:
		  PassingLeptonicW = 1.

		#Hadronic W matching
                if self.treeMine.AK8Puppijet0_isHadronicV == 1 and self.treeMine.AK8Puppijet0_vMatching < 0.8 and self.treeMine.AK8Puppijet0_vMatching > 0.0:
		  self.LeadingAK8Jet_MatchedHadW[0] = 1.
		else:
		  self.LeadingAK8Jet_MatchedHadW[0] = 0.

                self.triggerpassbb[0] = 1.
                #Determining the importance of each cut
		self.bb.Fill(self.triggerpassbb[0])
                if self.pfmet[0] > 40:
                  self.bb0.Fill(self.triggerpassbb[0])
                  if self.nTightMu[0] == 1:
                    self.bb1.Fill(self.triggerpassbb[0])
                    if self.nLooseMu[0] == 1:
                      self.bb2.Fill(self.triggerpassbb[0])
                      if HadronicAK8Jet == 1:
                        self.bb3.Fill(self.triggerpassbb[0])
                        if nBTaggedAK4Jet > 0:
                          self.bb4.Fill(self.triggerpassbb[0])		
			  if PassingLeptonicW > 0:
			    self.bb5.Fill(self.triggerpassbb[0])

		if self.treeMine.AK8Puppijet0_isHadronicV == 1:
                  self.bb6.Fill(self.triggerpassbb[0])
  	  	if self.treeMine.AK8Puppijet0_vMatching < 0.8:
		  self.bb7.Fill(self.triggerpassbb[0])
		if self.treeMine.AK8Puppijet0_vSize < 0.8:
		  self.bb8.Fill(self.triggerpassbb[0])

		#Determining whether to the event passes the selection
		if self.isMC == 'False':
                  if self.pfmet[0] > 40 and self.nTightMu[0] == 1 and self.nLooseMu[0] == 1 and HadronicAK8Jet == 1 and nBTaggedAK4Jet > 0 and PassingLeptonicW > 0 and self.treeMine.passJson == 1 and self.treeMine.triggerBits&4==4:
                     self.theTree.Fill()
                if self.isMC == 'True':
                  if self.pfmet[0] > 40 and self.nTightMu[0] == 1 and self.nLooseMu[0] == 1 and HadronicAK8Jet == 1 and nBTaggedAK4Jet > 0 and PassingLeptonicW > 0:
                     self.theTree.Fill()

	
            self.f1.Close()




