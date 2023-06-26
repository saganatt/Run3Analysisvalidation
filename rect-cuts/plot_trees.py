#!/usr/bin/env python
"""
file: plot_tree.py
brief: Plotting most important distributions for rectangular cuts.
       Input: trees from tree creator homogenized with homogenize_output.py and merged with hadd.
usage: ./plot_tree.py tree_signal.root tree_bkg.root
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
"""

from array import array
import argparse

# pylint: disable=import-error,no-name-in-module
from ROOT import (
    TCanvas,
    TFile,
    TLegend,
    TH1F,
    TH2F,
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


def plot_tree_single(tree, settings_var, hist_title, sig_or_bkg, htype):
    """
    Get a single histogram from a tree.
    """
    hist_color = { "sig": 1, "bkg": 2}
    # mc_flag == 1 << DecayType::LcToPKPi, decay type == 1
    tree_cond = { "sig": "(fMCflag == 2 || fMCflag == -2)", "bkg": "fMCflag == 0" }

    histname = f"h_{sig_or_bkg}_{settings_var[0]}"
    if htype == "TH1F":
        hist = TH1F(histname, f"{hist_title}", *settings_var[2])
        hist.SetLineColor(hist_color[sig_or_bkg])
        tree.Draw(f"{settings_var[1]}>>{histname}", tree_cond[sig_or_bkg])
    else:
        hist = TH2F(histname, f"{hist_title}", *settings_var[2], 7,
                    array("d", [0, 1, 2, 4, 6, 8, 12, 24]))
        hist.SetLineColor(hist_color[sig_or_bkg])
        tree.Draw(f"fPt:{settings_var[1]}>>{histname}", tree_cond[sig_or_bkg])

    int_hist = hist.Integral()
    if int_hist != 0.0:
        hist.Scale(1. / int_hist) # probability distribution, sum of content = 1.0

    return hist


def plot_single(hists, c_name):
    """
    Plot a single signal vs bkg comparison.
    """
    margin = 0.05
    k = 1.0 - 2 * margin

    canv = TCanvas(f"c_{c_name}", c_name, 800, 600)
    canv.cd()
    canv.SetGridx()
    canv.SetGridy()

    sig_entries = hists["sig"].GetEntries()
    bkg_entries = hists["bkg"].GetEntries()
    print(f"Sig histogram for {c_name} counts: {sig_entries} bkg: {bkg_entries}")
    hists["sig"].Draw("hist")
    hists["bkg"].Draw("hist;same")
    y_min = min(hists["sig"].GetMinimum(), hists["bkg"].GetMinimum())
    y_max = max(hists["sig"].GetMaximum(), hists["bkg"].GetMaximum())
    y_range = y_max - y_min
    hists["sig"].GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)
    hists["bkg"].GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)

    legend = TLegend(0.50, 0.72, 0.70, 0.90)
    legend.AddEntry(hists["sig"], "signal", "L")
    legend.AddEntry(hists["bkg"], "background", "L")
    legend.Draw()

    hists["sig"].Write()
    hists["bkg"].Write()
    save_canvas(canv, c_name)


def plot_hists(settings, pt_ranges, trees):
    """
    Plot histograms from trees.
    """

    for var, settings_var in settings.items():
        hists = {}
        htype = "TH1F" if var.startswith("pt") else "TH2F"
        for sig_or_bkg in ("sig", "bkg"):
            hists[sig_or_bkg] = plot_tree_single(trees[sig_or_bkg], settings_var,
                                                 f"Normalized {var}", sig_or_bkg,
                                                 htype)

        if var.startswith("pt"):
            plot_single(hists, settings_var[0])
        else:
            for i in range(len(pt_ranges) - 1):
                histname = f"{settings_var[0]}_pt_{pt_ranges[i]}-{pt_ranges[i + 1]}"
                hist_title = f"Normalized {var} for {pt_ranges[i]}" \
                             f" #leq #it{{p}}_{{T}} < {pt_ranges[i + 1]}"
                projs = {}
                for sig_or_bkg in ("sig", "bkg"):
                    ind = hists[sig_or_bkg].GetYaxis().FindBin(pt_ranges[i])
                    projs[sig_or_bkg] = hists[sig_or_bkg].ProjectionX(f"h_{sig_or_bkg}_{histname}",
                                                                      ind, ind + 1)
                    projs[sig_or_bkg].SetTitle(hist_title)

                plot_single(projs, histname)


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetFrameLineWidth(2)

    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument("sig_input_file", help="input signal tree AnalysisResults_tree.root file")
    parser.add_argument("bkg_input_file", help="input bkg tree AnalysisResults_tree.root file")
    parser.add_argument("output_file", help="output histograms file")

    args = parser.parse_args()

    settings = { "decay length": ("decay_length", "fDecayLength", [100, 0.0, 0.1]),
                 "decay length XY": ("decay_length_XY", "fDecayLengthXY", [100, 0.0, 0.1]),
                 "CPA": ("CPA", "fCPA", [100, 0.9, 1.0]),
                 "CPA XY": ("CPA_XY", "fCPAXY", [100, 0.9, 1.0]),
                 "Chi2PCA": ("Chi2PCA", "fChi2PCA", [200, 0.0, 0.01]),
                 "mass": ("mass", "fM", [600, 1.98, 2.58]),
                 "impact parameter 0": ("impact_parameter_0", "fImpactParameter0",
                                        [100, -0.02, 0.02]),
                 "impact parameter 1": ("impact_parameter_1", "fImpactParameter1",
                                        [100, -0.02, 0.02]),
                 "impact parameter 2": ("impact_parameter_2", "fImpactParameter2",
                                        [100, -0.02, 0.02]),
                 "pt": ("pt", "fPt", [200, 0, 24]),
                 "pt prong0": ("pt_prong0", "fPtProng0", [200, 0, 24]),
                 "pt prong1": ("pt_prong1", "fPtProng1", [200, 0, 24]),
                 "pt prong2": ("pt_prong2", "fPtProng2", [200, 0, 24])
                }
    pt_ranges = [0, 1, 2, 4, 6, 8, 12, 24]

    infile_sig = TFile(args.sig_input_file)
    infile_bkg = TFile(args.bkg_input_file)
    outfile = TFile(args.output_file, "RECREATE") # pylint: disable=unused-variable

    trees = {}
    trees["sig"] = infile_sig.Get("DF_0/O2hfcand3pfull")
    trees["bkg"] = infile_bkg.Get("DF_0/O2hfcand3pfull")

    plot_hists(settings, pt_ranges, trees)


if __name__ == "__main__":
    main()
