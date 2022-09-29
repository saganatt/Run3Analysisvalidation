#!/usr/bin/env python3

"""
Draw invariant mass distribution, fit signal and background, calculate yield
Arguments: input file, input MC file, output file without extension,
           HF task name, optional draw, optional doMC
Use: ./draw_inv_mass.py AnalysisResults.root MC_AnalysisResults.root Results d0 --doMC
"""

# pylint: disable=missing-function-docstring, invalid-name, too-many-arguments, too-many-locals

import argparse

from ROOT import TCanvas, TFile, gROOT # pylint: disable=import-error

import draw_utils
import fit_utils
import ranges

d0task = "hf-task-d0"
flowtask = "hf-task-flow"

ptbins = [1, 2, 4, 6, 8, 12, 24]
nptbins = len(ptbins)

def calculate_yield(fits, int_range):
    yields = {}
    yields["all"] = fits["fSigBkg"].Integral(int_range[0], int_range[1])
    yields["signal"] = fits["fSig"].Integral(int_range[0], int_range[1])
    yields["refl bkg"] = fits["fReflBkg"].Integral(int_range[0], int_range[1])
    if yields["refl bkg"] == 0.:
        yields["ratio"] = 0.
    else:
        yields["ratio"] = yields["signal"] / yields["refl bkg"]
    yields["bkg"] = fits["fBkg"].Integral(int_range[0], int_range[1])
    print(f'Full integral: {yields["all"]:.3f} '
          f'yield=signal: {yields["signal"]:.3f} '
          f'refl bkg: {yields["refl bkg"]:.3f} '
          f'ratio signal / refl bkg: {yields["ratio"]:.3f} '
          f'background: {yields["bkg"]:.3f}')
    return yields

def process_mc_single(hMassMC, histname, outfile, binmin, binmax,
                      init_range, fin_range, init_scale):
    hname = f"{histname}_{binmin}-{binmax}"
    #draw_utils.draw_mc_raw(hMassMC, hname, outfile)

    c = TCanvas(hname, hname)
    c.cd()

    # Draw invariant mass distribution
    ind1 = hMassMC.GetYaxis().FindBin(binmin)
    ind2 = hMassMC.GetYaxis().FindBin(binmax) - 1 # open upper interval
    print(f"Projecting MC, pt range: {binmin}, {binmax}, found bins: {ind1}, {ind2}")
    hMassProjMC = hMassMC.ProjectionX(hname, ind1, ind2, 1, hMassMC.GetNbinsZ())
    draw_utils.draw_inv_mass(hMassProjMC, (binmin, binmax), ranges.pt_range_an)

    # Fit
    params = fit_utils.fit_mc(hMassProjMC, init_range, fin_range, init_scale)

    # Draw histogram and fit
    fSig = draw_utils.draw_mc_fit(params, init_range)

    c.Write(hname)
    c.SaveAs(f"{outfile}_{hname}.png")

    return params, fSig

def process_mc(hMassSigD0, hMassReflBkgD0, outfile, binmin, binmax):
    print("-" * 40)
    print("Fitting MC signal")
    mc_sig_params, mc_sig_f = process_mc_single(hMassSigD0, "hMassSigD0", outfile,
                                                binmin, binmax,
                                                ranges.sig_mc_init_range, ranges.sig_mc_range,
                                                3000)
    # Global range needed for background fitting on real data
    ranges.sig_range[0] = ranges.sig_mc_range[0]
    ranges.sig_range[1] = ranges.sig_mc_range[1]
    print("-" * 40)
    print("Fitting MC reflected background")
    mc_bkg_params, mc_bkg_f = process_mc_single(hMassReflBkgD0, "hMassReflBkgD0", outfile,
                                                binmin, binmax,
                                                ranges.bkg_mc_init_range, ranges.bkg_mc_range,
                                                300)
    #ranges.bkg_range[0] = ranges.bkg_mc_range[0]
    #ranges.bkg_range[1] = ranges.bkg_mc_range[1]
    print("-" * 40)

    # Ratio of MC sig / refl bkg
    mc_yields = {}
    mc_yields["signal"] = mc_sig_f.Integral(ranges.sig_mc_range[0], ranges.sig_mc_range[1])
    mc_yields["bkg"] = mc_bkg_f.Integral(ranges.sig_mc_range[0], ranges.sig_mc_range[1])
    mc_yields["ratio"] = mc_yields["signal"] / mc_yields["bkg"]
    print(f'MC signal integral: {mc_yields["signal"]:.3f} '
          f'reflection background: {mc_yields["bkg"]:.3f} '
          f'ratio: {mc_yields["ratio"]:.3f}')

    return mc_sig_params, mc_bkg_params, mc_yields

def main(file, mcfile, outfile, task, drawMore, doMC):
    gROOT.SetBatch(True)

    print(f"Processing file {file}")
    f = TFile(file)
    hMass = f.Get(f"{task}/hMass")

    if doMC:
        print(f"Processing MC file {mcfile}")
        fMC = TFile(mcfile)
        hMassSigD0 = fMC.Get(f"{task}/hMassSigD0")
        hMassReflBkgD0 = fMC.Get(f"{task}/hMassReflBkgD0")

    fout = TFile(f"{outfile}.root", "RECREATE")

    for ind in range(nptbins - 1):
        binmin = ptbins[ind]
        binmax = ptbins[ind + 1]
        print(f"Processing pT bin: [{binmin}, {binmax})")

        if doMC:
            mc_sig_params, mc_bkg_params, mc_yields = process_mc(hMassSigD0, hMassReflBkgD0,
                                                                 outfile,
                                                                 binmin, binmax)
        else:
            mc_sig_params = {}
            mc_bkg_params = {}
            mc_yields = {100, 200}

        hname = f"hMass_{binmin}-{binmax}"
        c = TCanvas(hname, hname)
        c.cd()

        # Draw invariant mass distribution
        hMassProj = hMass.ProjectionX(hname, ind + 1, ind + 2)
        hMassProj.GetYaxis().SetRangeUser(0, 30000)
        draw_utils.draw_inv_mass(hMassProj, (binmin, binmax), ranges.pt_range_an)

        # Fit signal and background
        params = fit_utils.fit_sig_bkg(hMassProj, ranges.sig_range, ranges.pt_range_an,
                                       doMC, mc_sig_params, mc_bkg_params, mc_yields)

        # Draw histogram and fits
        fits = draw_utils.draw_fits(params, ranges.pt_range_an, ranges.sig_range)

        # Optionally: draw vertical lines and shade signal area
        if drawMore:
            draw_utils.draw_more(ranges.sig_range, fits, c)

        # Calculate yield
        yields = calculate_yield(fits, ranges.sig_range)

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
