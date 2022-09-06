#!/usr/bin/env python3

"""
Draw invariant mass distribution in different pt bins.
Arguments: input file, output file without extension, HF task name
Use: ./compare_ratio.py AnalysisResults_O2.root Ratios flow
"""

import argparse

from ROOT import TH1, TF1, TH2, TCanvas, TFile, gROOT
from ROOT import kFullCircle, kBlack, kGreen, kBlue, kRed
from ROOT import TMath

d0task = "hf-task-d0"
flowtask = "hf-task-flow"

pt_d0_to_pik = [1.61, 2.12] # cuts used during skimming
pt_range_an = [1.7, 2.05] # inv mass x-range in the AN

ptbins = [1, 2, 4, 6, 8, 12, 24]
sig_range = [1.84, 1.91]
nptbins = len(ptbins)

def background(x, params):
    if x[0] > sig_range[0] and x[0] < sig_range[1]:
        TF1.RejectPoint()
        return 0
    return TMath.Exp(params[0] + params[1] * x[0])

def signal(x, params):
    if x[0] <= sig_range[0] or x[0] >= sig_range[1]:
        TF1.RejectPoint()
        return 0
    if params[2] == 0.:
        return 0
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2]))

def sigbkg(x, params):
    if params[2] == 0.:
        return 0
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2])) + TMath.Exp(params[3] + params[4] * x[0])

def main(file, outfile, task):
    gROOT.SetBatch(True)
    print(f"Processing file {file}")
    f = TFile(file)
    fout = TFile(f"{outfile}.root", "RECREATE")

    hMass = f.Get(f"{task}/hMass")

    for ind in range(nptbins - 1):
        binmin = ptbins[ind]
        binmax = ptbins[ind + 1]

        # Draw invariant mass distribution
        hname = f"hMass_{binmin}-{binmax}"
        hMassProj = hMass.ProjectionX(hname, ind + 1, ind + 2)
        #hMassProj.GetXaxis().SetRangeUser(pt_d0_to_pik[0], pt_d0_to_pik[1])
        hMassProj.GetXaxis().SetRangeUser(pt_range_an[0], pt_range_an[1])
        hMassProj.SetTitle("%s < #it{p}_{T} < %s" % (binmin, binmax))

        c = TCanvas(hname, hname)
        c.cd()

        hMassProj.SetMarkerStyle(kFullCircle)
        hMassProj.SetLineColor(kBlack)
        hMassProj.SetMarkerColor(kBlack)
        #hMassProj.SetMarkerSize(2)
        hMassProj.Draw("PE")

        # Fit signal and background (background and signal)
        fitBkg = TF1("fitBkg", background, pt_range_an[0], pt_range_an[1], 2)
        hMassProj.Fit(fitBkg, "NQ+", "", pt_range_an[0], pt_range_an[1])

        hMassSig = hMassProj.Clone("hMassSig")
        hMassSig.Add(fitBkg, -1.)
        nbins = hMassSig.GetNbinsX()
        for i in range(int(nbins)):
            if hMassSig.GetBinContent(i + 1) < 0.:
                hMassSig.SetBinContent(i + 1, 0.)
        estSigma = (sig_range[1] - sig_range[0]) / 6.
        estMean = sig_range[0] + (sig_range[1] - sig_range[0]) / 2
        fitSig = TF1("fitSig", signal, pt_range_an[0], pt_range_an[1], 3)
        fitSig.SetParameters(15000, estMean, estSigma)
        hMassSig.Fit(fitSig, "NQ", "", pt_range_an[0], pt_range_an[1])

        # Draw fits
        fBkg = TF1("fBkg", "expo", pt_range_an[0], pt_range_an[1])
        offsetBkg = fitBkg.GetParameter(0)
        scaleBkg = fitBkg.GetParameter(1)
        fBkg.SetParameters(offsetBkg, scaleBkg)
        print(f"Background parameters: offset: {offsetBkg} scale: {scaleBkg}")
        fBkg.SetLineColor(kRed)
        fBkg.Draw("same")

        scaleSig = fitSig.GetParameter(0)
        meanSig = fitSig.GetParameter(1)
        sigmaSig = fitSig.GetParameter(2)
        print(f"Signal parameters: scale: {scaleSig} mean: {meanSig} sigma: {sigmaSig}")
        #fSigBkg = TF1("fSigBkg", "gaus(0)+expo(3)", pt_range_an[0], pt_range_an[1])
        #fSigBkg.SetParameters(scaleSig, meanSig, sigmaSig, offsetBkg, scaleBkg)
        #fSigBkg.SetLineColor(kGreen)
        #fSigBkg.Draw("same")

        fSigBkgSig = TF1("fSigBkgSig", "gaus(0)+expo(3)", sig_range[0], sig_range[1])
        fSigBkgSig.SetParameters(scaleSig, meanSig, sigmaSig, offsetBkg, scaleBkg)
        fSigBkgSig.SetLineColor(kBlue)
        fSigBkgSig.Draw("same")

        # Calculate yield
        intAll = fSigBkgSig.Integral(sig_range[0], sig_range[1])
        intBkg = fBkg.Integral(sig_range[0], sig_range[1])
        yieldSig = intAll - intBkg
        print(f"Full integral: {intAll} background: {intBkg} yield: {yieldSig}")

        c.Write(hname)
        c.SaveAs(f"{outfile}_{hname}.png")

        # Draw signal-only histogram
        c2 = TCanvas(f"{hname}-sigonly", f"{hname}-sigonly")
        c2.cd()
        hMassSig.SetMarkerStyle(kFullCircle)
        hMassSig.SetLineColor(kBlack)
        hMassSig.SetMarkerColor(kBlack)
        hMassSig.Draw("PE")
        fSig = TF1("fSig", "gaus", pt_range_an[0], pt_range_an[1])
        fSig.SetParameters(scaleSig, meanSig, sigmaSig)
        fSig.SetLineColor(kBlue)
        fSig.Draw("same")
        c2.Write(f"{hname}-sigonly")
        c2.SaveAs(f"{outfile}_{hname}_sigonly.png")

        yieldSigSub = fSig.Integral(sig_range[0], sig_range[1])
        print(f"Yield from signal hist: {yieldSigSub}")

    fout.Close()

if __name__ == "__main__":
    pass
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=str, help="Input file")
    parser.add_argument("outfile", type=str, help="Output file")
    parser.add_argument("task", type=str, help="HF final task name")
    args = parser.parse_args()
    task = d0task
    if args.task == "flow":
        task = flowtask
    main(file=args.file, outfile=args.outfile, task=task)
