#!/usr/bin/env python3

"""
Draw invariant mass distribution in different pt bins.
Arguments: input file, output file without extension, HF task name
Use: ./compare_ratio.py AnalysisResults_O2.root Ratios flow
"""

import argparse

from ROOT import TH1, TH2, TCanvas, TFile, gROOT

d0task = "hf-task-d0"
flowtask = "hf-task-flow"

ptbins = [1, 2, 4, 6, 8, 12, 24]
nptbins = len(ptbins)

def main(file, outfile, task):
    gROOT.SetBatch(True)
    print(f"Processing file {file}")
    f = TFile(file)
    fout = TFile(f"{outfile}.root", "RECREATE")

    hMass = f.Get(f"{task}/hMass")

    for ind in range(nptbins - 1):
        binmin = ptbins[ind]
        binmax = ptbins[ind + 1]

        hname = f"hMass_{binmin}-{binmax}"
        hMassProj = hMass.ProjectionX(hname, ind + 1, ind + 2)

        c = TCanvas(hname, hname)
        c.cd()
        hMassProj.Draw()
        c.Write(hname)
        c.SaveAs(f"{outfile}_{hname}.png")

    fout.Close()

if __name__ == "__main__":
    pass
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=str, help="Input file")
    parser.add_argument("outfile", type=str, help="Output file")
    parser.add_argument("task", type=str, help="HF final task name")
    args = parser.parse_args()
    task = d0task
    if args.task == "flow":
        task = flowtask
    main(file=args.file, outfile=args.outfile, task=task)
