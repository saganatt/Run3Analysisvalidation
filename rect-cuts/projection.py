#!/usr/bin/env python
"""
file: projection.py
brief: Simple projection macro that prepares AnalysisResults.root for the PWGHF mass fitter.
usage: ./projection.py AnalysisResults.root hf-task-lc/Data/hMassVsPt
                       ~/alice/O2Physics/PWGHF/D2H/Macros/config_massfitter.json invmass.root
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
"""

import argparse
import json

from ROOT import TFile, gROOT  # pylint: disable=import-error

task_bin_edges = [-1, 0, 1, 2, 4, 6, 8, 12, 24]  # -1 is dummy to mark the underflow bin


def main():
    """
    Main
    """

    gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", type=str, help="Input file")
    parser.add_argument("histpath", type=str, help="Path to histogram to project")
    parser.add_argument(
        "fitter_config_file", type=str, help="Mass fitter JSON config file"
    )
    parser.add_argument("outfile", type=str, help="Output file")
    args = parser.parse_args()

    with open(args.fitter_config_file, "r", encoding="utf-8") as fitter_config_f:
        fitter_config_text = fitter_config_f.read()
    fitter_config = json.loads(fitter_config_text)
    pt_min = fitter_config["PtMin"]
    pt_max = fitter_config["PtMax"]

    print(f"Processing file {args.infile}")
    fin = TFile(args.infile)
    hist = fin.Get(f"{args.histpath}")

    fout = TFile(args.outfile, "RECREATE")

    ind = 0
    nbins = len(task_bin_edges)
    for binmin, binmax in zip(pt_min, pt_max):
        print(f"Processing pT bin: [{binmin}, {binmax})")
        hname = f"proj_{binmin}-{binmax}"
        while ind < nbins and task_bin_edges[ind] < binmin + 0.001:
            ind += 1
        minind = ind - 1
        while ind < nbins and task_bin_edges[ind] < binmax - 0.001:
            ind += 1
        maxind = ind - 1
        if maxind > nbins - 2 or minind > nbins - 2:
            raise ValueError(
                f"Requested range [{pt_min}, {pt_max}) outside task maximum {task_bin_edges[-1]}"
            )

        h_proj = hist.ProjectionX(hname, minind, maxind)
        # hMassProj.GetYaxis().SetRangeUser(0, 30000)

        h_proj.Write()

    fout.Close()


if __name__ == "__main__":
    main()
