#!/usr/bin/env python3

"""
Draw invariant mass distribution in different pt bins.
Arguments: input file, output file without extension, HF task name
Use: ./compare_ratio.py AnalysisResults_O2.root Ratios flow
"""

import argparse

from ROOT import TH1, TF1, TH2, TCanvas, TFile, gROOT
from ROOT import kFullCircle, kBlack, kGreen, kBlue, kRed, kMagenta
from ROOT import TMath

d0task = "hf-task-d0"
flowtask = "hf-task-flow"

pt_d0_to_pik = [1.61, 2.12] # cuts used during skimming
pt_range_an = [1.7, 2.05] # inv mass x-range in the AN

ptbins = [1, 2, 4, 6, 8, 12, 24]
sig_init_range = [1.85, 1.9]
sig_range = [1.86, 1.9]
nptbins = len(ptbins)

def background(x, params):
    if x[0] > sig_init_range[0] and x[0] < sig_init_range[1]:
        TF1.RejectPoint()
        return 0
    return TMath.Exp(params[0] + params[1] * x[0])

def signal(x, params):
    if x[0] <= sig_init_range[0] or x[0] >= sig_init_range[1]:
        TF1.RejectPoint()
        return 0
    if params[2] == 0.:
        return 0
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2]))

def gauss_shifted(x, params):
    if x[0] <= sig_init_range[0] or x[0] >= sig_init_range[1]:
        TF1.RejectPoint()
        return 0
    if params[2] == 0.:
        return 0
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2])) + params[3]

def exp_gauss_shifted(x, params):
    if params[2] == 0.:
        return 0
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2])) + params[3] + TMath.Exp(params[4] + params[5] * x[0])

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
        print(f"Processing pT bin: [{binmin}, {binmax})")

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

        # Fit signal and background
        estSigma = (sig_init_range[1] - sig_init_range[0]) / 6.
        estMean = sig_init_range[0] + (sig_init_range[1] - sig_init_range[0]) / 2
        fitSig = TF1("fitSig", gauss_shifted, sig_init_range[0], sig_init_range[1], 4)
        fitSig.SetParameters(20000, estMean, estSigma, 5000)
        hMassProj.Fit(fitSig, "NQ", "", sig_init_range[0], sig_init_range[1])
        scaleSig = fitSig.GetParameter(0)
        meanSig = fitSig.GetParameter(1)
        sigmaSig = fitSig.GetParameter(2)
        shiftSig = fitSig.GetParameter(3)
        sig_range[:] = [meanSig - 3. * sigmaSig, meanSig + 3. * sigmaSig]
        print(f"Signal parameters: scale: {scaleSig:.3f} mean: {meanSig:.3f} sigma: {sigmaSig:.3f} shift: {shiftSig:.3f}\n"
                f"Initial signal region: {sig_init_range[0]:.2f}, {sig_init_range[1]:.2f}, final: {sig_range[0]:.2f}, {sig_range[1]:.2f}")

        fitBkg = TF1("fitBkg", background, pt_range_an[0], pt_range_an[1], 2)
        hMassProj.Fit(fitBkg, "NQ+", "", pt_range_an[0], pt_range_an[1])
        offsetBkg = fitBkg.GetParameter(0)
        scaleBkg = fitBkg.GetParameter(1)
        print(f"Background parameters: offset: {offsetBkg:.3f} scale: {scaleBkg:.3f}")

        fitSigBkg = TF1("fitSigBkg", exp_gauss_shifted, pt_range_an[0], pt_range_an[1], 6)
        fitSigBkg.SetParameters(20000, meanSig, sigmaSig, 5000, offsetBkg, scaleBkg)
        hMassProj.Fit(fitSigBkg, "NQ+", "", pt_range_an[0], pt_range_an[1])
        scaleSigBkg = fitSigBkg.GetParameter(0)
        meanSigBkg = fitSigBkg.GetParameter(1)
        sigmaSigBkg = fitSigBkg.GetParameter(2)
        shiftSigBkg = fitSigBkg.GetParameter(3)
        offsetSigBkg = fitSigBkg.GetParameter(4)
        expScaleSigBkg = fitSigBkg.GetParameter(5)
        print(f"Joint fit parameters: scale: {scaleSigBkg:.3f} mean: {meanSigBkg:.3f} sigma: {sigmaSigBkg:.3f} shift: {shiftSigBkg:.3f} "
                f"offset: {offsetSigBkg:.3f} exp scale: {expScaleSigBkg:.3f}")

        # Draw fits
        fBkg = TF1("fBkg", "expo", pt_range_an[0], pt_range_an[1])
        fBkg.SetParameters(offsetBkg, scaleBkg)
        fBkg.SetLineColor(kRed)
        fBkg.Draw("same")

        fSig = TF1("fSig", gauss_shifted, sig_init_range[0], sig_init_range[1], 4)
        fSig.SetParameters(scaleSig, meanSig, sigmaSig, shiftSig)
        fSig.SetLineColor(kGreen)
        fSig.Draw("same")

        fSigBkg = TF1("fSigBkg", exp_gauss_shifted, pt_range_an[0], pt_range_an[1], 6)
        fSigBkg.SetParameters(scaleSigBkg, meanSigBkg, sigmaSigBkg, shiftSigBkg, offsetSigBkg, expScaleSigBkg)
        fSigBkg.SetLineColor(kBlue)
        fSigBkg.Draw("same")

        # Calculate yield
        intAll = fSigBkg.Integral(sig_range[0], sig_range[1])
        intSig = fSig.Integral(sig_range[0], sig_range[1])
        intBkg = fBkg.Integral(sig_range[0], sig_range[1])
        yieldSig = intAll - intBkg
        print(f"Full integral: {intAll:.3f} sig: {intSig:.3f} background: {intBkg:.3f} yield: {yieldSig:.3f}")

        c.Write(hname)
        c.SaveAs(f"{outfile}_{hname}.png")
        print("*" * 80)

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
