#!/usr/bin/env python
"""
file: multiprojection.py
brief: A projection macro that generates overlapping 1D invariant mass histograms
       from several input AnalysisResults.root files
usage: ./multiprojection.py config_multiprojection.json
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
"""

import argparse
import json

# pylint: disable=import-error,no-name-in-module
from ROOT import (
    TCanvas,
    TFile,
    TLatex,
    TLegend,
    gROOT,
    gStyle,
    kBlack,
    kRed,
    kBlue,
    kGreen,
)

TASK_BIN_EDGES = [-1, 0, 1, 2, 4, 6, 8, 12, 24] # -1 is dummy to mark the underflow bin
NBINS = len(TASK_BIN_EDGES)
HIST_COLORS = [kBlack, kRed+1, kBlue-3, kGreen+3]
MARGIN = 0.05
SCALE = 1.0 - 2 * MARGIN

def save_canvas(canvas, title):
    """
    Save canvas in png, pdf, root.
    """
    format_list = [".png"]
    for file_format in format_list:
        canvas.SaveAs(title + file_format)


def get_bin_indices(config):
    """
    Calculate bin indices of histograms for pT limits provided in JSON.
    """
    ind = 0
    indices = []
    for binmin, binmax in zip(config["pt_min"], config["pt_max"]):
        while ind < NBINS and TASK_BIN_EDGES[ind] < binmin + 0.001:
            ind += 1
        minind = ind - 1
        while ind < NBINS and TASK_BIN_EDGES[ind] < binmax - 0.001:
            ind += 1
        maxind = ind - 1
        if maxind > NBINS - 2 or minind > NBINS - 2:
            raise ValueError(f'Requested range [{config["pt_min"]}, {config["pt_max"]}) '
                             f'outside task maximum {TASK_BIN_EDGES[-1]}')
        indices.append((minind, maxind))

    return indices


def plot(indices, config, files):
    """
    Get and plot invariant mass projections.
    """
    for (minind, maxind), binmin, binmax in zip(indices, config["pt_min"], config["pt_max"]):
        print(f"Processing pT bin: [{binmin}, {binmax})")
        hname = f"proj_{binmin}-{binmax}"
        canv = TCanvas(f"c_{hname}", f"c_{hname}")
        legend = TLegend(0.75, 0.70, 0.90, 0.90)
        legend.SetTextSize(0.03)

        latexa = TLatex()
        latexa.SetTextSize(0.08)
        latexa.SetTextFont(62) # bold helvetica
        latexa.SetTextAlign(2) # centered

        projs = []

        for ind, (file, label) in enumerate(zip(files, config["labels"])):
            print(f"Processing file {label}")
            hist = file.Get(config["hist_path"])
            h_proj = hist.ProjectionX(f"{hname}_{label}", minind, maxind)
            h_proj = h_proj.Rebin(8)
            h_proj.SetLineColor(HIST_COLORS[ind])
            if ind == 0:
                h_proj.Draw()
            else:
                h_proj.Draw("same")
            legend.AddEntry(h_proj, label, "L")
            projs.append(h_proj)

        y_min = min((proj.GetMinimum() for proj in projs))
        y_max = max((proj.GetMaximum() for proj in projs))
        y_range = y_max - y_min

        for h_proj in projs:
            h_proj.GetYaxis().SetRangeUser(0.0, y_max + MARGIN / SCALE * y_range)
            h_proj.Write()

        legend.Draw()
        latexa.DrawLatexNDC(0.72, 0.95, config["output_label"])
        save_canvas(canv, f'{config["output_dir"]}/{hname}-{config["output_label"]}')


def main():
    """
    Main
    """
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", type=str, help="JSON configuration file")
    args = parser.parse_args()

    with open(args.config_file, "r", encoding="utf-8") as config_f:
        config_text = config_f.read()
    config = json.loads(config_text)

    files = []
    for file in config["input_files"]:
        fin = TFile(file)
        files.append(fin)

    fout = TFile(f'{config["output_dir"]}/{config["output_file"]}', "RECREATE")

    indices = get_bin_indices(config)

    plot(indices, config, files)

    fout.Close()

if __name__ == "__main__":
    main()
