#! /usr/bin/env python3

"""
Compare histograms from a single file by producing their ratios.
Use: ./compare_ratio.py AnalysisResults_O2.root
"""

import argparse

from ROOT import TH1, TCanvas, TColor, TFile, TLegend, gPad

histnames = ["hptcand_befsel", "hptcand_wrongdecay", "hptcand_wrongy", "hptcand"]
mc_histnames = ["hPtRecSig", "hPtRecBg", "hPtRecSigPrompt", "hPtRecSigNonPrompt"]
d0task = "hf-task-d0"
corrtask = "task-hf-correlations"

def save_ratios(ratios):
    for obj in ratios:
        obj.SaveAs("Ratios.pdf")

    fout = TFile("Ratios.root", "RECREATE")
    for obj in ratios:
        print("Writing", obj.GetName())
        obj.Write(obj.GetName().replace("/", "_folder_"))
    fout.Close()


def main(file):
    f = TFile(file)
    ratios = []
    drawn = {}

    for hn in histnames:
        print(f"Drawing ratio for {hn}")
        hd0 = f.Get(f"{d0task}/{hn}")
        hcorr = f.Get(f"{corrtask}/{hn}")
        ratio = hcorr.Clone("ratio")
        ratio.Divide(hd0)

        c = TCanvas(hn, hn)
        c.cd()
        ratio.Draw()
        ratios.append(ratio)
        drawn[hn] = [c, ratio]

        for i in range(ratio.GetNbinsX()):
            binc = ratio.GetBinContent(i)
            binhd0 = hd0.GetBinContent(i)
            binhcorr = hcorr.GetBinContent(i)
            if binc != 1.0 and binc != 0.0:
                print(f"Bin: {i} contains: {binc}, d0: {binhd0}, hcorr: {binhcorr}")

    save_ratios(ratios)


if __name__ == "__main__":
    pass
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=str, help="Input file")
    args = parser.parse_args()
    main(file=args.file)
