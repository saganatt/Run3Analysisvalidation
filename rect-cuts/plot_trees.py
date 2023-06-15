#!/usr/bin/env python
"""
file: plot_tree.py
brief: Plotting most important distributions for rectangular cuts
usage: ./plot_tree.py AnalysisResults_tree_signal.root AnalysisResults_tree_bkg.root
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
"""

import argparse
from array import array

# pylint: disable=import-error,no-name-in-module
from ROOT import (
    TCanvas,
    TFile,
    TLegend,
    TH1F,
    gROOT,
    gStyle,
)


def save_canvas(canvas, title):
    """
    Save canvas in png, pdf, root.
    """
    format_list = [".png", ".pdf", ".root"]
    for file_format in format_list:
        canvas.SaveAs(title + file_format)


def prepare_hists(settings):
    """
    Initialize histograms.
    """
    hists_sig = { var: [] for var in settings["var_list"] }
    hists_bkg = { var: [] for var in settings["var_list"] }
    gStyle.SetOptStat(0)
    gStyle.SetFrameLineWidth(2)
    for i in range(len(settings["pt_ranges"]) - 1):
        for var, varu, var_range in zip(settings["var_list"], settings["var_list_u"],
                                        settings["var_ranges"]):
            histname = f'{varu}_pt_{settings["pt_ranges"][i]}-{settings["pt_ranges"][i + 1]}'
            hist_title = f'Normalized {var} for {settings["pt_ranges"][i]} #leq #it{{p}}_{{T}} <' \
                         f'{settings["pt_ranges"][i + 1]}'
            hist_sig = TH1F(f"h_sig_{histname}", hist_title, *var_range)
            hist_bkg = TH1F(f"h_bkg_{histname}", hist_title, *var_range)
            hist_sig.SetLineColor(1)
            hist_bkg.SetLineColor(2)
            hists_sig[var].append(hist_sig)
            hists_bkg[var].append(hist_bkg)
    return hists_sig, hists_bkg


def fill_hists(settings, infile, hists, sig_or_bkg, outfile): # pylint: disable=too-many-locals
    """
    Fill histograms from trees.
    """
    tree_name = "O2hfcand3pfull"
    for key in infile.GetListOfKeys(): # pylint: disable=too-many-nested-blocks
        key_name = key.GetName()
        if key_name.startswith("DF_"): # is the dataframe directory
            tree = infile.Get(f"{key_name}/{tree_name}")
            pt_val = array("f", [ 0. ])
            tree.SetBranchAddress("fPt", pt_val)
            mc_flag = array("b", [ 0 ])
            tree.SetBranchAddress("fMCflag", mc_flag)
            var_leaves = { var: array("f", [0.]) for var in settings["var_list"] }
            for var, leaf in zip(settings["var_list"], settings["leaf_list"]):
                tree.SetBranchAddress(leaf, var_leaves[var])
            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                # mc_flag == 1 << DecayType::LcToPKPi, decay type == 1
                if (sig_or_bkg == "sig" and (mc_flag[0] == 2 or mc_flag[0] == -2)) \
                        or (sig_or_bkg == "bkg" and mc_flag[0] == 0):
                    for j in range(len(settings["pt_ranges"]) - 1):
                        if settings["pt_ranges"][j] <= pt_val[0] < settings["pt_ranges"][j + 1]:
                            for var in var_leaves:
                                hists[var][j].Fill(var_leaves[var][0])
                            break
    for var in settings["var_list"]:
        for i in range(len(settings["pt_ranges"]) - 1):
            outfile.WriteObject(hists[var][i], hists[var][i].GetName())


def plot_single(h_sig, h_bkg, filename):
    """
    Plot a single signal vs background comparison.
    """
    canv = TCanvas(f"c_{filename}", filename, 800, 600)
    canv.cd()
    canv.SetGridx()
    canv.SetGridy()

    int_sig = h_sig.Integral()
    int_bkg = h_bkg.Integral()
    if int_sig != 0.0:
        h_sig.Scale(1. / int_sig) # probability distribution, sum of content = 1.0
    if int_bkg != 0.0:
        h_bkg.Scale(1. / int_bkg) # probability distribution, sum of content = 1.0
    h_sig.Draw("hist")
    h_bkg.Draw("hist;same")
    y_min = min(h_sig.GetMinimum(), h_bkg.GetMinimum())
    y_max = max(h_sig.GetMaximum(), h_bkg.GetMaximum())
    margin = 0.05
    k = 1.0 - 2 * margin
    y_range = y_max - y_min
    h_sig.GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)
    h_bkg.GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)

    legend = TLegend(0.50, 0.72, 0.70, 0.90)
    legend.AddEntry(h_sig, "signal", "L")
    legend.AddEntry(h_bkg, "background", "L")
    legend.Draw()

    save_canvas(canv, filename)


def plot(hists_sig, hists_bkg, settings):
    """
    Plot distributions for each variable.
    """
    for var, varu in zip(settings["var_list"], settings["var_list_u"]):
        for i, (h_sig, h_bkg) in enumerate(zip(hists_sig[var], hists_bkg[var])):
            filename = f'{varu}_pt_{settings["pt_ranges"][i]}-{settings["pt_ranges"][i + 1]}'
            plot_single(h_sig, h_bkg, filename)


def plot_file(outfile, settings):
    """
    Plot distributions from file.
    """
    for varu in settings["var_list_u"]:
        for i in range(len(settings["pt_ranges"]) - 1):
            histname = f'{varu}_pt_{settings["pt_ranges"][i]}-{settings["pt_ranges"][i + 1]}'
            h_sig = outfile.Get(f"h_sig_{histname}")
            h_bkg = outfile.Get(f"h_bkg_{histname}")
            plot_single(h_sig, h_bkg, histname)


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)
    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument("sig_input_file", help="input signal tree AnalysisResults_tree.root file")
    parser.add_argument("bkg_input_file", help="input bkg tree AnalysisResults_tree.root file")
    parser.add_argument("outfile", help="output file with all histograms saved for later plotting")
    parser.add_argument("--plot_file", default=False, action="store_true",
                        help="if set, plot histograms from file instead of calculating from trees")

    args = parser.parse_args()

    settings = { "var_list": ["decay length", "decay length XY", "CPA", "CPA XY", "Chi2PCA"],
                 "var_list_u": ["decay_length", "decay_length_XY", "CPA", "CPA_XY", "Chi2PCA"],
                 "leaf_list": ["fDecayLength", "fDecayLengthXY", "fCPA", "fCPAXY", "fChi2PCA"],
                 "var_ranges": [[100, 0.0, 0.1], [100, 0.0, 0.1], [100, 0.9, 1.], [100, 0.9, 1.],
                                [200, 0., 0.01]],
                 "pt_ranges": [0, 1, 2, 4, 6, 8, 12, 24]
                }

    if args.plot_file:
        outfile = TFile(args.outfile)
        plot_file(outfile, settings)
    else:
        outfile = TFile(args.outfile, "RECREATE")
        infile_sig = TFile(args.sig_input_file)
        infile_bkg = TFile(args.bkg_input_file)

        hists_sig, hists_bkg = prepare_hists(settings)

        print("Filling histos")
        fill_hists(settings, infile_sig, hists_sig, "sig", outfile)
        print("Filled histos for signal")
        fill_hists(settings, infile_bkg, hists_bkg, "bkg", outfile)
        print("Filled histos for background")

        plot(hists_sig, hists_bkg, settings)

        infile_sig.Close()
        infile_bkg.Close()

    outfile.Close()


if __name__ == "__main__":
    main()
