#!/usr/bin/env python3

"""
Draw invariant mass distribution, fit signal and background, calculate yield
Arguments: input file, input MC file, output file without extension,
           HF task name, optional draw, optional doMC
Use: ./draw_inv_mass.py AnalysisResults.root MC_AnalysisResults.root Results d0 --doMC
"""

# pylint: disable=missing-function-docstring, invalid-name, too-many-locals, too-many-arguments

import argparse

from ROOT import TCanvas, TFile, gROOT # pylint: disable=import-error

import draw_utils
import fit_utils
import ranges

d0task = "hf-task-d0"
flowtask = "hf-task-flow"

ptbins = [1, 2, 4, 6, 8, 12, 24]
nptbins = len(ptbins)

def process_mc(hMassMC, histname, outfile, binmin, binmax, init_range, fin_range):
    ind1 = hMassMC.FindBin(binmin)
    ind2 = hMassMC.FindBin(binmax)

    hname = f"{histname}_{binmin}-{binmax}"
    draw_utils.draw_mc_raw(hMassMC, hname, outfile)

    c = TCanvas(hname, hname)
    c.cd()

    # Draw invariant mass distribution
    binsZ = hMassMC.GetNbinsZ()
    hMassProjMC = hMassMC.ProjectionX(hname, ind1, ind2, 1, binsZ)
    draw_utils.draw_inv_mass(hMassProjMC, (binmin, binmax), ranges.pt_range_an)

    # Fit
    params = fit_utils.fit_mc(hMassProjMC, init_range, fin_range)

    # Draw histogram and fit
    draw_utils.draw_mc_fit(params, init_range)

    c.Write(hname)
    c.SaveAs(f"{outfile}_{hname}.png")

    return params

def main(file, mcfile, outfile, task, drawMore, doMC):
    gROOT.SetBatch(True)
    print(f"Processing file {file}")
    f = TFile(file)
    if doMC:
        print(f"Processing MC file {mcfile}")
        fMC = TFile(mcfile)
    fout = TFile(f"{outfile}.root", "RECREATE")

    hMass = f.Get(f"{task}/hMass")
    if doMC:
        hMassSigD0 = fMC.Get(f"{task}/hMassSigD0")
        hMassReflBkgD0 = fMC.Get(f"{task}/hMassReflBkgD0")

    for ind in range(nptbins - 1):
        binmin = ptbins[ind]
        binmax = ptbins[ind + 1]
        print(f"Processing pT bin: [{binmin}, {binmax})")

        if doMC:
            print("-" * 40)
            print("Fitting MC signal")
            mc_sig_params = process_mc(hMassSigD0, "hMassSigD0", outfile, binmin, binmax,
                                       ranges.sig_mc_init_range, ranges.sig_mc_range)
            mcSigmaSignal = mc_sig_params["sigma"]
            print("-" * 40)
            print("Fitting MC reflected background")
            mc_bkg_params = process_mc(hMassReflBkgD0, "hMassReflBkgD0", outfile, binmin, binmax,
                                       ranges.bkg_mc_init_range, ranges.bkg_mc_range)
            mcSigmaReflBkg = mc_bkg_params["sigma"]
            print("-" * 40)
        else:
            mcSigmaSignal = -1
            mcSigmaReflBkg = -1

        hname = f"hMass_{binmin}-{binmax}"
        c = TCanvas(hname, hname)
        c.cd()

        # Draw invariant mass distribution
        hMassProj = hMass.ProjectionX(hname, ind + 1, ind + 2)
        draw_utils.draw_inv_mass(hMassProj, (binmin, binmax), ranges.pt_range_an)

        # Fit signal and background
        params = fit_utils.fit_sig_bkg(hMassProj, doMC, mcSigmaSignal, mcSigmaReflBkg)

        # Draw histogram and fits
        fits = draw_utils.draw_fits(params, ranges.pt_range_an, ranges.sig_init_range)

        # Optionally: draw vertical lines and shade signal area
        if drawMore:
            draw_utils.draw_more(ranges.sig_range, fits, c)

        # Calculate yield
        intAll = fits["fSigBkg"].Integral(ranges.sig_range[0], ranges.sig_range[1])
        intSig = fits["fSig"].Integral(ranges.sig_range[0], ranges.sig_range[1])
        intBkg = fits["fBkg"].Integral(ranges.sig_range[0], ranges.sig_range[1])
        yieldSig = intAll - intBkg
        print(f"Full integral: {intAll:.3f} sig: {intSig:.3f} "
              f"background: {intBkg:.3f} yield: {yieldSig:.3f}")

        c.Write(hname)
        c.SaveAs(f"{outfile}_{hname}.png")
        print("*" * 80)

    fout.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=str, help="Input file")
    parser.add_argument("mcfile", type=str, help="Input MC file")
    parser.add_argument("outfile", type=str, help="Output file")
    parser.add_argument("task", type=str, help="HF final task name")
    parser.add_argument("--draw", default=False, action="store_true",
                        help="Draw vertical lines and shade signal area")
    parser.add_argument("--doMC", default=False, action="store_true",
                        help="Process MC histograms")
    args = parser.parse_args()
    seltask = d0task
    if args.task == "flow":
        seltask = flowtask
    main(file=args.file, mcfile=args.mcfile, outfile=args.outfile,
         task=seltask, drawMore=args.draw, doMC=args.doMC)
