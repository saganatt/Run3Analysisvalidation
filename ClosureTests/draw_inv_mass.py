#!/usr/bin/env python3

"""
Draw invariant mass distribution in different pt bins.
Arguments: input file, output file without extension, HF task name, optional draw
Use: ./compare_ratio.py AnalysisResults_O2.root Ratios flow --draw
"""

import argparse

from ROOT import TH1, TF1, TH2, TLine, TCanvas, TFile, gROOT
from ROOT import kFullCircle, kDashed, kBlack, kGreen, kBlue, kRed, kMagenta, kGray
from ROOT import TMath, TGraph

d0task = "hf-task-d0"
flowtask = "hf-task-flow"

pt_d0_to_pik = [1.61, 2.12] # cuts used during skimming
pt_range_an = [1.7, 2.05] # inv mass x-range in the AN

ptbins = [1, 2, 4, 6, 8, 12, 24]
sig_init_range = [1.85, 1.9]
sig_range = [1.85, 1.9]
nptbins = len(ptbins)

def background(x, params):
    if x[0] > sig_range[0] and x[0] < sig_range[1]:
        TF1.RejectPoint()
        return 0
    return TMath.Exp(params[0] + params[1] * x[0])

def signal(x, params):
    if x[0] <= sig_init_range[0] or x[0] >= sig_init_range[1]:
        TF1.RejectPoint()
        return 0
    #if params[2] == 0.:
    #    return 0
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2])) +\
            params[3] * TMath.Exp(-0.5 * ((x[0] - params[4]) / params[5]) * ((x[0] - params[4]) / params[5]))

def sigbkg(x, params):
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2]) * ((x[0] - params[1]) / params[2])) +\
            params[3] * TMath.Exp(-0.5 * ((x[0] - params[4]) / params[5]) * ((x[0] - params[4]) / params[5])) +\
           TMath.Exp(params[6] + params[7] * x[0])

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

def fit_sig_bkg(hist, doMC, mcSigmaSignal = -1, mcSigmaReflBkg = -1):
    params = {}
    estSigma1 = (sig_init_range[1] - sig_init_range[0]) / 6.
    estSigma2 = (sig_init_range[1] - sig_init_range[0]) / 6.
    if doMC:
        estSigma1 = mcSigmaSignal
        estSigma2 = mcSigmaReflBkg
        print(f"MC sigma from signal: {mcSigmaSignal:.3f} reflected background: {mcSigmaReflBkg:.3f}")

    estMean = sig_init_range[0] + (sig_init_range[1] - sig_init_range[0]) / 2
    fitSig = TF1("fitSig", "gaus(0)+gaus(3)", sig_init_range[0], sig_init_range[1], 6)
    fitSig.SetParameters(20000, estMean, estSigma1, 20000, estMean, estSigma2)
    hist.Fit(fitSig, "NQ", "", sig_init_range[0], sig_init_range[1])
    params["scaleSig1"] = fitSig.GetParameter(0)
    params["meanSig1"] = fitSig.GetParameter(1)
    params["sigmaSig1"] = fitSig.GetParameter(2)
    params["scaleSig2"] = fitSig.GetParameter(3)
    params["meanSig2"] = fitSig.GetParameter(4)
    params["sigmaSig2"] = fitSig.GetParameter(5)
    sig_range[:] = [params["meanSig1"] - 3. * params["sigmaSig1"], params["meanSig1"] + 3. * params["sigmaSig1"]]
    print(f'Signal parameters: scale: {params["scaleSig1"]:.3f} mean: {params["meanSig1"]:.3f} '
            f'sigma: {params["sigmaSig1"]:.3f}\nscale: {params["scaleSig2"]:.3f} '
            f'mean: {params["meanSig2"]:.3f} sigma: {params["sigmaSig2"]:.3f}\n'
            f'Initial signal region: {sig_init_range[0]:.2f}, {sig_init_range[1]:.2f}, '
            f'final: {sig_range[0]:.2f}, {sig_range[1]:.2f}')

    fitBkg = TF1("fitBkg", background, pt_range_an[0], pt_range_an[1], 2)
    hist.Fit(fitBkg, "NQ+", "", pt_range_an[0], pt_range_an[1])
    params["offsetBkg"] = fitBkg.GetParameter(0)
    params["scaleBkg"] = fitBkg.GetParameter(1)
    print(f'Background parameters: offset: {params["offsetBkg"]:.3f} scale: {params["scaleBkg"]:.3f}')

    fitSigBkg = TF1("fitSigBkg", "gaus(0)+gaus(3)+expo(6)", pt_range_an[0], pt_range_an[1], 8)
    fitSigBkg.SetParameters(params["scaleSig1"], params["meanSig1"], params["sigmaSig1"], params["scaleSig2"], params["meanSig2"], params["sigmaSig2"], params["offsetBkg"], params["scaleBkg"])
    hist.Fit(fitSigBkg, "NQ+", "", pt_range_an[0], pt_range_an[1])
    params["scaleSigBkg1"] = fitSigBkg.GetParameter(0)
    params["meanSigBkg1"] = fitSigBkg.GetParameter(1)
    params["sigmaSigBkg1"] = fitSigBkg.GetParameter(2)
    params["scaleSigBkg2"] = fitSigBkg.GetParameter(0)
    params["meanSigBkg2"] = fitSigBkg.GetParameter(1)
    params["sigmaSigBkg2"] = fitSigBkg.GetParameter(2)
    params["offsetSigBkg"] = fitSigBkg.GetParameter(4)
    params["expScaleSigBkg"] = fitSigBkg.GetParameter(5)
    print(f'Joint fit parameters: scale: {params["scaleSigBkg1"]:.3f} mean: {params["meanSigBkg1"]:.3f} '
            f'sigma: {params["sigmaSigBkg1"]:.3f}\nscale: {params["scaleSigBkg2"]:.3f} '
            f'mean: {params["meanSigBkg2"]:.3f} sigma: {params["sigmaSigBkg2"]:.3f}\n'
            f'offset: {params["offsetSigBkg"]:.3f} exp scale: {params["expScaleSigBkg"]:.3f}')

    return params

def draw_inv_mass(hist, binmin, binmax):
    #hist.GetXaxis().SetRangeUser(pt_d0_to_pik[0], pt_d0_to_pik[1])
    hist.GetXaxis().SetRangeUser(pt_range_an[0], pt_range_an[1])
    hist.SetTitle("%s #leq #it{p}_{T} < %s" % (binmin, binmax))
    hist.SetMarkerStyle(kFullCircle)
    hist.SetLineColor(kBlack)
    hist.SetMarkerColor(kBlack)
    #hist.SetMarkerSize(2)
    hist.Draw("PE")

def get_vertical_line(x, c):
    c.Update()
    yMin = c.GetUymin()
    yMax = c.GetUymax()
    lVert = TLine(x, yMin, x, yMax)
    lVert.SetLineColor(kBlack)
    lVert.SetLineStyle(kDashed)
    return lVert

def draw_fits(params):
    fits = {}
    fBkg = TF1("fBkg", "expo", pt_range_an[0], pt_range_an[1])
    fBkg.SetParameters(params["offsetBkg"], params["scaleBkg"])
    fBkg.SetLineColor(kRed)
    fBkg.Draw("same")
    fits["fBkg"] = fBkg

    fSig = TF1("fSig", "gaus(0)+gaus(3)", sig_init_range[0], sig_init_range[1], 6)
    fSig.SetParameters(params["scaleSig1"], params["meanSig1"], params["sigmaSig1"], params["scaleSig2"], params["meanSig2"], params["sigmaSig2"])
    fSig.SetLineColor(kGreen)
    fSig.Draw("same")
    fits["fSig"] = fSig

    fSigBkg = TF1("fSigBkg", "gaus(0)+gaus(3)+expo(6)", pt_range_an[0], pt_range_an[1], 8)
    fSigBkg.SetParameters(params["scaleSigBkg1"], params["meanSigBkg1"], params["sigmaSigBkg1"],
                          params["scaleSigBkg2"], params["meanSigBkg2"], params["sigmaSigBkg2"],
                          params["offsetSigBkg"], params["expScaleSigBkg"])
    fSigBkg.SetLineColor(kBlue)
    fSigBkg.Draw("same")
    fits["fSigBkg"] = fSigBkg

    return fits

def shade_signal(c, fits):
    gr = TGraph()
    gr.SetFillColor(kGray)
    gr.SetFillStyle(3013)
    c.Update()
    ymin = c.GetUymin()
    ymax = c.GetUymax()
    xmin = sig_range[0]
    xmax = sig_range[1]

    npx = fits["fSig"].GetNpx()
    npoints = 0
    dx = (xmax - xmin) / npx

    x = xmin + 0.5 * dx
    while x <= xmax:
        y = fits["fSigBkg"].Eval(x)
        if y < ymin:
            y = ymin
        if y > ymax:
            y = ymax
        gr.SetPoint(npoints, x, y)
        npoints += 1
        x += dx

    x = xmax - 0.5 * dx
    while x >= xmin:
        y = fits["fBkg"].Eval(x)
        if y < ymin:
            y = ymin
        if y > ymax:
            y = ymax
        gr.SetPoint(npoints, x, y)
        npoints += 1
        x -= dx

    return gr

def fit_mc(hMassMC, histname, binmin, binmax):
    ind1 = hMassMC.FindBin(binmin)
    ind2 = hMassMC.FindBin(binmax)

    hname = f"{histname}_{binmin}-{binmax}"
    c = TCanvas(hname, hname)
    c.cd()

    # Draw invariant mass distribution
    hMassProjMC = hMassMC.ProjectionX(ind1, ind2)
    draw_inv_mass(hMassProjMC, binmin, binmax)

    # Fit signal
    params = {}

    # Draw histogram and fits

    c.Write(hname)
    c.SaveAs(f"{outfile}_{hname}.png")

    return params

def main(file, outfile, task, drawMore, doMC):
    gROOT.SetBatch(True)
    print(f"Processing file {file}")
    f = TFile(file)
    fout = TFile(f"{outfile}.root", "RECREATE")

    hMass = f.Get(f"{task}/hMass")
    if doMC:
        hMassSigD0 = f.Get(f"{task}/hMassSigD0")
        hMassReflBkgD0 = f.Get(f"{task}/hMassReflBkgD0")

    for ind in range(nptbins - 1):
        binmin = ptbins[ind]
        binmax = ptbins[ind + 1]
        print(f"Processing pT bin: [{binmin}, {binmax})")

        if doMC:
            print("-" * 40)
            print("Fitting MC signal")
            mc_sig_params = fit_mc(hMassSigD0, "hMassSigD0", binmin, binmax)
            #mcSigmaSignal = mc_sig_params["sigma"]
            print("-" * 40)
            print("Fitting MC reflected background")
            mc_bkg_params = fit_mc(hMassReflBkgD0, "hMassReflBkgD0", binmin, binmax)
            #mcSigmaReflBkg = mc_bkg_params["sigma"]
            print("-" * 40)
        mcSigmaSignal = -1
        mcSigmaReflBkg = -1

        hname = f"hMass_{binmin}-{binmax}"
        c = TCanvas(hname, hname)
        c.cd()

        # Draw invariant mass distribution
        hMassProj = hMass.ProjectionX(hname, ind + 1, ind + 2)
        draw_inv_mass(hMassProj, binmin, binmax)

        # Fit signal and background
        params = fit_sig_bkg(hMassProj, doMC, mcSigmaSignal, mcSigmaReflBkg)

        # Draw histogram and fits
        fits = draw_fits(params)

        # Optionally: draw vertical lines and shade signal area
        if drawMore:
            lMin = get_vertical_line(sig_range[0], c)
            lMin.Draw("same")
            lMax = get_vertical_line(sig_range[1], c)
            lMax.Draw("same")
            gr = shade_signal(c, fits)
            gr.Draw("f")

        # Calculate yield
        intAll = fits["fSigBkg"].Integral(sig_range[0], sig_range[1])
        intSig = fits["fSig"].Integral(sig_range[0], sig_range[1])
        intBkg = fits["fBkg"].Integral(sig_range[0], sig_range[1])
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
    parser.add_argument("--draw", default=False, action="store_true", help="Draw vertical lines and shade signal area")
    parser.add_argument("--doMC", default=False, action="store_true", help="Process MC histograms")
    args = parser.parse_args()
    task = d0task
    if args.task == "flow":
        task = flowtask
    main(file=args.file, outfile=args.outfile, task=task, drawMore=args.draw, doMC=args.doMC)
