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
import json

# pylint: disable=import-error,no-name-in-module
from ROOT import (
    TCanvas,
    TFile,
    TLegend,
    TH1F,
    TH2F,
    gROOT,
    gStyle,
    kBlack,
    kRed,
    kBlue,
    kGreen,
)


def save_canvas(canvas, title):
    """
    Save canvas in png, pdf, root.
    """
    format_list = [".png", ".pdf", ".root"]
    for file_format in format_list:
        canvas.SaveAs(title + file_format)


def plot_tree_single(tree, settings_var, hist_title, selection, hcolor, htype):
    """
    Get a single histogram from a tree.
    """
    # mc_flag == 1 << DecayType::LcToPKPi, decay type == 1
    tree_cond = { "sig": "fM > 2.226 && fM < 2.346" if settings_var[1] != "fM" else "", "bkg": "fFlagMc == 0",
                  "sig_mc": "(fFlagMc == 2 || fFlagMc == -2)"}

    histname = f"h_{sig_or_bkg}_{settings_var[0]}"
    if htype == "TH1F":
        hist = TH1F(histname, f"{hist_title}", *settings_var[2])
        hist.SetLineColor(hcolor)
        tree.Draw(f"{settings_var[1]}>>{histname}", tree_cond[selection])
    else:
        hist = TH2F(histname, f"{hist_title}", *settings_var[2], 7,
                    array("d", [0, 1, 2, 4, 6, 8, 12, 24]))
        hist.SetLineColor(hcolor)
        tree.Draw(f"fPt:{settings_var[1]}>>{histname}", tree_cond[selection])

    int_hist = hist.Integral()
    if int_hist != 0.0 and settings_var[1] != "fM":
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

    counts_str = f"Sig histogram for {c_name} counts:"
    option = "hist"
    for hist in hists:
        entries = hist.GetEntries()
        counts_str = f"{counts_str} {entries}"
        hist.Draw(option)
        option = "hist;same"
    print(counts_str)

    legend = TLegend(0.50, 0.72, 0.70, 0.90)

    y_min = min((hist.GetMinimum() for hist in hists))
    y_max = min((hist.GetMaximum() for hist in hists))
    y_range = y_max - y_min
    for label, hist in hists.items():
        hist.GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)
        hist.GetYaxis().SetTitle("Normalized counts")
        legend.AddEntry(hist, label, "L")
    legend.Draw()

    for hist in hists:
        hist.Write()
    save_canvas(canv, c_name)


def plot_total_mass(hists_mass, c_name, hist_title):
    """
    Plot total mass from signal MC for mass fitter.
    """
    canv = TCanvas(f"c_{c_name}", c_name, 800, 600)
    canv.cd()
    canv.SetGridx()
    canv.SetGridy()

    for i in sorted(hists_mass)[:1]:
        mass_total = hists_mass[i].Clone(f"h_{c_name}")
    for i in sorted(hists_mass)[2:]:
        mass_total.Add(hists_mass[i])
    mass_total.SetTitle(hist_title)
    mass_total.Draw()
    mass_total.Write()
    save_canvas(canv, c_name)


def plot_hists(settings, pt_ranges, trees, selections):
    """
    Plot histograms from trees.
    """

    hist_colors = [kBlack, kRed+1, kBlue-3, kGreen+3]
    for var, settings_var in settings.items():
        hists = {}
        htype = "TH1F" if "#it{p}_{T}" in var else "TH2F"
        for ind, label, tree, selection in zip(trees.items(), selections):
            hists[label] = plot_tree_single(tree, settings_var, f"Normalized {var}",
                                            selection, hist_colors[ind], htype)

        if htype == "TH1F":
            plot_single(hists, settings_var[0])
        else:
            for i in range(len(pt_ranges) - 1):
                histname = f"{settings_var[0]}_pt_{pt_ranges[i]}-{pt_ranges[i + 1]}"
                hist_title = f"{var} for {pt_ranges[i]}" \
                             f" #leq #it{{p}}_{{T}} < {pt_ranges[i + 1]}"
                projs = {}
                for label in trees:
                    ind = hists[label].GetYaxis().FindBin(pt_ranges[i])
                    ind2 = hists[label].GetYaxis().FindBin(pt_ranges[i + 1] - 0.05)
                    projs[label] = hists[label].ProjectionX(f"h_{label}_{histname}",
                                                                      ind, ind2)
                    projs[label].SetTitle(hist_title)

                plot_single(projs, histname)

                if var == "mass":
                    hist_title = f"Total mass for {pt_ranges[i]}" \
                                 f" #leq #it{{p}}_{{T}} < {pt_ranges[i + 1]}"
                    plot_total_mass(projs, f"total_{histname}", hist_title)


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetFrameLineWidth(2)

    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument("config_file", type=str, help="JSON config file")
    args = parser.parse_args()

    with open(args.config_file, "r", encoding="utf-8") as config_f:
        config_text = config_f.read()
    config = json.loads(config_text)

    settings = { "mass": ("mass", "fM", [600, 2.18, 2.38]),
                 "decay length": ("decay_length", "fDecayLength", [100, 0.0, 0.1]),

                 "decay length XY": ("decay_length_XY", "fDecayLengthXY", [100, 0.0, 0.1]),
                 "CPA": ("CPA", "fCpa", [100, 0.9, 1.0]),
                 "CPA XY": ("CPA_XY", "fCpaXY", [100, 0.9, 1.0]),
                 "Chi2PCA": ("Chi2PCA", "fChi2PCA", [400, 0.0, 2.0]),
                 "impact parameter 0": ("impact_parameter_0", "fImpactParameter0",
                                        [100, -0.02, 0.02]),
                 "impact parameter 1": ("impact_parameter_1", "fImpactParameter1",
                                        [100, -0.02, 0.02]),
                 "impact parameter 2": ("impact_parameter_2", "fImpactParameter2",
                                        [100, -0.02, 0.02]),
                 "#Lambda_{c} #it{p}_{T}": ("pt", "fPt", [200, 0, 24]),
                 "#it{p}_{T} prong_{0}": ("pt_prong0", "fPtProng0", [200, 0, 24]),
                 "#it{p}_{T} prong_{1}": ("pt_prong1", "fPtProng1", [200, 0, 24]),
                 "#it{p}_{T} prong_{2}": ("pt_prong2", "fPtProng2", [200, 0, 24])
                }
    pt_ranges = [1, 2, 4, 6, 8, 12, 24]

    files = []
    for file in config["input_files"]:
        files.append(TFile(file))
    outfile = TFile(config["output_file"], "RECREATE") # pylint: disable=unused-variable

    trees = {}
    for file, label in zip(files, config["labels"]):
        trees[label] = file.Get("DF_0/O2hfcandlcfull")

    plot_hists(settings, pt_ranges, trees, config["selections"])


if __name__ == "__main__":
    main()
