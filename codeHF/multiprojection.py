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

from ROOT import TFile, gROOT, TCanvas # pylint: disable=import-error

TASK_BIN_EDGES = [-1, 0, 1, 2, 4, 6, 8, 12, 24] # -1 is dummy to mark the underflow bin
NBINS = len(TASK_BIN_EDGES)

def save_canvas(canvas, title):
    """
    Save canvas in png, pdf, root.
    """
    format_list = [".png"]
    for file_format in format_list:
        canvas.SaveAs(title + file_format)


def main():
    """
    Main
    """
    gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file", type=str, help="JSON configuration file")
    args = parser.parse_args()

    with open(args.config_file, "r", encoding="utf-8") as config_f:
        config_text = config_f.read()
    config = json.loads(config_text)

    print(f"Processing file {args.infile}")
    fin = TFile(config["input_file"])
    hist = fin.Get(config["histpath"])
    fout = TFile(f'{config["output_dir"]}/{config["output_file"]}', "RECREATE")

    ind = 0
    for binmin, binmax in zip(config["pt_min"], config["pt_max"]):
        print(f"Processing pT bin: [{binmin}, {binmax})")
        hname = f"proj_{binmin}-{binmax}"
        canv = TCanvas(f"c_{hname}", f"c_{hname}")

        while ind < NBINS and TASK_BIN_EDGES[ind] < binmin + 0.001:
            ind += 1
        minind = ind - 1
        while ind < NBINS and TASK_BIN_EDGES[ind] < binmax - 0.001:
            ind += 1
        maxind = ind - 1
        if maxind > NBINS - 2 or minind > NBINS - 2:
            raise ValueError(f'Requested range [{config["pt_min"]}, {config["pt_max"]}) '
                             f'outside task maximum {TASK_BIN_EDGES[-1]}')

        h_proj = hist.ProjectionX(hname, minind, maxind)
        h_proj = h_proj.Rebin(8)
        h_proj.Draw()
        #hMassProj.GetYaxis().SetRangeUser(0, 30000)

        h_proj.Write()
        save_canvas(canv, f'{config["output_dir"]}/{hname}')

    fout.Close()

if __name__ == "__main__":
    main()
