#!/usr/bin/env python
"""
file: plot_histos.py
brief: Plotting most important distributions for rectangular cuts
usage: ./plot_histos.py AnalysisResults_signal.root AnalysisResults_bkg.root Lc
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
"""

import argparse

# pylint: disable=import-error,no-name-in-module
from ROOT import (
    TCanvas,
    TFile,
    TLegend,
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


def prepare_canvas(var):
    """
    Initialize canvas.
    """
    canv = TCanvas(f"c_{var}", var)
    canv.SetCanvasSize(800, 600)
    canv.cd()
    canv.SetGridy()
    canv.SetGridx()
    return canv


def get_proj(hist, var, var_range, var_file, pt_range, sig_or_bkg):
    """
    Get the projected histogram.
    """
    bin1 = hist.GetYaxis().FindBin(pt_range[0])
    bin2 = hist.GetYaxis().FindBin(pt_range[1])
    h_proj = hist.ProjectionX(f"hproj_{sig_or_bkg}_{var_file}_"
                              f"pt_{pt_range[0]}-{pt_range[1]}", bin1, bin2)
    h_proj.Scale(1. / h_proj.Integral()) # probability distribution, sum of content = 1.0

    h_proj.SetTitle(f"Normalized {var} for {pt_range[0]} #leq #it{{p}}_{{T}} < {pt_range[1]}")
    h_proj.GetXaxis().SetRangeUser(var_range[0], var_range[1])
    if sig_or_bkg == "sig":
        h_proj.SetLineColor(1)
    else:
        h_proj.SetLineColor(2)
    #hist.GetYaxis().CenterTitle()
    #hist.GetXaxis().CenterTitle()
    #hist.GetXaxis().SetNoExponent()
    gStyle.SetOptStat(0)
    gStyle.SetFrameLineWidth(2)
    #gStyle.SetTitleSize(0.045, "x")
    #gStyle.SetTitleSize(0.045, "y")
    #gStyle.SetMarkerSize(1)
    #gStyle.SetLabelOffset(0.015, "x")
    #gStyle.SetLabelOffset(0.02, "y")
    #gStyle.SetTitleOffset(1.1, "x")
    #gStyle.SetTitleOffset(1.0, "y")
    return h_proj


def plot(hist_sig, hist_bkg, var, var_range, pt_range):
    """
    Plot a 1D histogram for a given variable and pT range.
    """
    canv = prepare_canvas(f"{var}_pt_{pt_range[0]}-{pt_range[1]}")
    var_file = var.replace(" ", "-")

    sig_h = get_proj(hist_sig, var, var_range, var_file, pt_range, "sig")
    bkg_h = get_proj(hist_bkg, var, var_range, var_file, pt_range, "bkg")
    sig_h.Draw("hist")
    bkg_h.Draw("hist;same")
    legend = TLegend(0.50, 0.72, 0.70, 0.90)
    legend.AddEntry(sig_h, "signal", "L")
    legend.AddEntry(bkg_h, "background", "L")
    legend.Draw()

    save_canvas(canv, f"{var_file}_pt_{pt_range[0]}-{pt_range[1]}")


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)
    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument("sig_input_file", help="input signal AnalysisResults.root")
    parser.add_argument("bkg_input_file", help="input background AnalysisResults.root")
    parser.add_argument("part", help="particle kind (Lc or D0)")

    args = parser.parse_args()

    all_settings = [{ "hist_path": "hf-task-lc/Data",
                  "hist_list": ["hDecLengthVsPt", "hDecLengthxyVsPt", "hCPAVsPt", "hCPAxyVsPt"],
                  "var_list": ["decay length", "decay length XY", "CPA", "CPA XY"],
                  "var_ranges": [[-0.05, 0.45], [-0.05, 0.45], [0.45, 1.05], [0.45, 1.05]],
                  "pt_ranges": [0, 1, 2, 4, 6, 8, 12, 24]
                },
                { "hist_path": "hf-task-d0",
                  "hist_list": ["hDecLengthFinerBinning", "hDecLengthxyFinerBinning",
                                "hCPAFinerBinning"],
                  "var_list": ["decay length", "decay length XY", "CPA"],
                  "var_ranges": [[-0.05, 0.45], [-0.05, 0.45], [0.45, 1.05], [0.45, 1.05]],
                  "pt_ranges": [0, 1, 2, 4, 6, 8, 12, 24]
                }]
    if args.part == "Lc":
        settings = all_settings[0]
    else:
        settings = all_settings[1]

    infile_sig = TFile(args.sig_input_file)
    infile_bkg = TFile(args.bkg_input_file)

    for histname, var, var_range in \
            zip(settings["hist_list"], settings["var_list"], settings["var_ranges"]):
        print(f'Getting histogram: {settings["hist_path"]}/{histname}')
        hist_sig = infile_sig.Get(f'{settings["hist_path"]}/{histname}')
        hist_bkg = infile_bkg.Get(f'{settings["hist_path"]}/{histname}')
        for ind in range(len(settings["pt_ranges"]) - 1):
            pt_range = (settings["pt_ranges"][ind], settings["pt_ranges"][ind + 1])
            plot(hist_sig, hist_bkg, var, var_range, pt_range)


if __name__ == "__main__":
    main()
