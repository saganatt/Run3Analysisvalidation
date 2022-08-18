#! /usr/bin/env python3

"""
Compare histograms from a single file by producing their ratios.
Arguments: input file, output file without extension
Use: ./compare_ratio.py AnalysisResults_O2.root Ratios
"""

import argparse

from ROOT import TH1, TCanvas, TColor, TFile, TLegend, gPad, gROOT

histnames = ["hPtCand_befsel", "hPtCand_wrongdecay", "hPtCand_wrongy", "hPtCand"]
mc_histnames = ["hPtRecSig", "hPtRecBg", "hPtRecSigPrompt", "hPtRecSigNonPrompt"]
d0task = "hf-task-d0"
corrtask = "hf-task-flow"

def main(file, outfile):
    gROOT.SetBatch(True)
    print(f"Processing file {file}")
    f = TFile(file)
    fout = TFile(f"{outfile}.root", "RECREATE")
    ratios = []

    for hn in histnames:
        print(f"Drawing ratio for {hn}")
        hd0 = f.Get(f"{d0task}/{hn}")
        hcorr = f.Get(f"{corrtask}/{hn}")
        ratio = hcorr.Clone("ratio")
        ratio.Divide(hd0)

        c = TCanvas(hn, hn)
        c.cd()
        ratio.Draw()

        ratio.Write(ratio.GetName().replace("/", "_folder_"))

        for i in range(ratio.GetNbinsX()):
            binc = ratio.GetBinContent(i)
            binhd0 = hd0.GetBinContent(i)
            binhcorr = hcorr.GetBinContent(i)
            if binc != 1.0 and binc != 0.0:
                print(f"Bin: {i} contains: {binc}, d0: {binhd0}, hfcorr: {binhcorr}")

        ol = f"both_{hn}"
        c2 = TCanvas(ol, ol)
        c2.cd()
        hd0.SetLineColor(TColor.GetColor("#e41a1c"))
        hd0.DrawClone()
        hcorr.SetLineColor(TColor.GetColor("#377eb8"))
        hcorr.DrawClone("same")

        d0_entries = hd0.GetEntries()
        corr_entries = hcorr.GetEntries()
        if d0_entries != corr_entries:
            raise ValueError(f"Different number of entries! d0: {d0_entries} hcorr: {corr_entries} hist: {hn}")

        c2.Write(ol)
        c2.SaveAs(f"{outfile}_{ol}.png")

    fout.Close()

if __name__ == "__main__":
    pass
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=str, help="Input file")
    parser.add_argument("outfile", type=str, help="Output file")
    args = parser.parse_args()
    main(file=args.file, outfile=args.outfile)
