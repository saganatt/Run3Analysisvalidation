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

    for ind, (binmin, binmax) in enumerate(zip(pt_min, pt_max)):
        print(f"Processing pT bin: [{binmin}, {binmax})")
        hname = f"proj_{binmin}-{binmax}"

        h_proj = hist.ProjectionX(hname, ind + 1, ind + 2)
        # hMassProj.GetYaxis().SetRangeUser(0, 30000)

        h_proj.Write()

    fout.Close()


if __name__ == "__main__":
    main()
