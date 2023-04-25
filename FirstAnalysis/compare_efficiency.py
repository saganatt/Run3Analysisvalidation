#!/usr/bin/env python
"""
file: compare_efficiency.py
brief: Plotting macro to compare O2 and AliPhysics efficiencies.
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
usage: ./compare_efficiency.py <O2 *.root file pattern> <Ali input *.root file>
       ./compare_efficiency.py efficiency_tracking_ITS-TPC Merged_ALI.root

For O2, efficiency_tracking_*.root files are produced by efficiency_studies macro and contain
a canvas with all particles respective efficiencies.

The AliPhysics file is produced by the Run3Analysisvalidation workflow using
AliAnalysisTaskTrackingEfficiencyPID.cxx and ProcessTrackingEffPID.C
"""

import argparse
from array import array

# pylint: disable=import-error,no-name-in-module
from ROOT import TH1F, TCanvas, TFile, TLatex, TLegend, gROOT, gStyle, kWhite


def save_canvas(canvas, title):
    """
    Save canvas in png, pdf, root.
    """
    format_list = [".png", ".pdf", ".root"]
    for file_format in format_list:
        canvas.SaveAs(title + file_format)


def get_titles(var, sign, had, det, ytitle):
    """
    Compose titles and names for histograms and canvas.
    """
    hname = f"hempty_{sign}_{had}_{det}_{var}"
    cname = f"c_{sign}_{had}_{det}_{var}"
    ctitle = f"{ytitle} for {sign} {had}, {det}"
    if var in ("Pt", "Eta"):
        ctitle = f"{ctitle} primaries"
    return hname, cname, ctitle, ytitle


def prepare_canvas(var, titles, single):
    """
    Initialize canvas, axes, legend.
    `hempty` must be captured at return, otherwise ROOT crashes.
    """

    hname, cname, ctitle, ytitle = titles

    def get_pt_hist():
        hempty = TH1F(
            hname, f"{ctitle};Transverse Momentum (GeV/c);{ytitle}", 43, array('d', [0, 0.05, 0.075, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.225, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 16, 20, 30, 35]))
        # gPad.SetLogx()
        return hempty

    def get_eta_hist():
        return TH1F(hname, f"{ctitle};Pseudorapidity;{ytitle}", 12, -0.8 - 0.16, 0.8 + 0.16)

    def get_phi_hist():
        return TH1F(
            hname,
            f"{ctitle};Azimuthal angle (rad);{ytitle}",
            20,
            0 - 6.28319 / 18,
            6.28319 + 6.28319 / 18
        )

    hists = {"Pt": get_pt_hist, "Eta": get_eta_hist, "Phi": get_phi_hist}

    canv = TCanvas(cname, "Efficiency")
    canv.SetCanvasSize(800, 600)
    canv.cd()
    canv.SetGridy()
    canv.SetGridx()

    hempty = hists[var]()
    hempty.GetYaxis().CenterTitle()
    hempty.GetXaxis().CenterTitle()
    hempty.GetXaxis().SetNoExponent()
    hempty.GetXaxis().SetMoreLogLabels(1)
    hempty.GetYaxis().SetNdivisions(11)
    hempty.Draw()

    latexa = TLatex()
    latexa.SetTextSize(0.04)
    latexa.SetTextFont(42)
    latexa.SetTextAlign(3)
    latexa.DrawLatexNDC(0.15, 0.3, "-0.9 #geq #eta #geq 0.9")
    latexa.DrawLatexNDC(0.15, 0.25, "-2#pi #geq #varphi #geq 2#pi")

    leg = TLegend(0.55, 0.15, 0.89, 0.25 if single else 0.35, "P")
    leg.SetNColumns(2)
    leg.SetHeader("Minimum bias pp #sqrt{s} = 5.02TeV", "C")
    leg.SetFillColor(0)
    if not single:
        hempty.SetMarkerColor(kWhite)
        leg.AddEntry(hempty, "Run 2", "p")
        leg.AddEntry(hempty, "Run 3", "p")

    return canv, leg, hempty


def retrieve_points(canvas, iso2):
    """
    Retrieve efficiency points from input canvas.
    AliPhysics efficiency is in TH1D, O2 -- in a TGraph.
    """
    clp = canvas.GetListOfPrimitives()
    eff = []
    leg_labels = []

    def o2_cond(elem):
        return iso2 and elem.GetName() == "eff_graph"

    def ali_cond(elem):
        return not iso2 and type(elem).__name__ == "TH1D"

    for elem in clp:
        if o2_cond(elem) or ali_cond(elem):
            eff.append(elem)
        elif type(elem).__name__ == "TLegend":
            count = 0
            for pelem in elem.GetListOfPrimitives():
                if iso2 and count > 0:
                    lab = pelem.GetLabel()
                    leg_labels.append(lab)
                count += 1

    typestr = "O2" if iso2 else "AliPhysics"
    if iso2 and (len(eff) != len(leg_labels)):
        raise RuntimeError(
            f"Different number of plots and legend entries in {typestr} file"
        )

    return eff, leg_labels


def compare_efficiency(
    alifile, canvo2, var, sign
):  # pylint: disable=too-many-locals, too-many-statements
    """
    Compare O2 vs AliPhysics efficiency vs pT, eta, phi for all hadron species.
    """
    ali_marker_list = [47, 21, 22, 34, 20]  # full figures
    o2_marker_list = [46, 25, 26, 28, 24]  # open figures

    o2_to_ali_signs = { "Positive": "pos", "Negative": "neg", "All": "posneg"}
    o2_to_ali_names = {"Electron": "e", "Pion": "pi", "Proton": "p", "Kaon": "K"}

    gStyle.SetOptStat(0)
    gStyle.SetFrameLineWidth(2)
    gStyle.SetTitleSize(0.045, "x")
    gStyle.SetTitleSize(0.045, "y")
    gStyle.SetMarkerSize(1)
    gStyle.SetLabelOffset(0.015, "x")
    gStyle.SetLabelOffset(0.02, "y")
    gStyle.SetTitleOffset(1.1, "x")
    gStyle.SetTitleOffset(1.0, "y")
    gStyle.SetErrorX(0)

    results = prepare_canvas(
        var, get_titles(var, sign, "all", "ITS-TPC", "Efficiency"), False
    )
    c_all = results[0]
    leg_all = results[1]

    o2_eff, o2_leg_labels = retrieve_points(canvo2, True)

    for ind, (o2_elem, o2_lab) in enumerate(zip(o2_eff, o2_leg_labels)
    ):
        if o2_lab == "All":
            continue
        ali_elem = alifile.Get(f"hEff{var}_{o2_to_ali_names[o2_lab]}_{o2_to_ali_signs[sign]}")
        ali_elem.SetMarkerStyle(ali_marker_list[ind])
        o2_elem.SetMarkerStyle(o2_marker_list[ind])

        c_all.cd()
        ali_elem.Draw(" same p")
        leg_all.AddEntry(ali_elem, o2_lab, "p")
        o2_elem.Draw(" same p")
        leg_all.AddEntry(o2_elem, o2_lab, "p")

        results_single = prepare_canvas(
            var, get_titles(var, sign, o2_lab, "ITS-TPC", "Efficiency"), True
        )
        c_single = results_single[0]
        leg_single = results_single[1]

        c_single.cd()
        ali_elem.Draw(" same p")
        leg_single.AddEntry(ali_elem, "Run 2", "p")
        o2_elem.Draw(" same p")
        leg_single.AddEntry(o2_elem, "Run 3", "p")

        leg_single.Draw()
        save_canvas(c_single, f"comparison_efficiency_ITS-TPC_{sign}_{var}_{o2_lab}")

    c_all.cd()
    leg_all.Draw()
    save_canvas(c_all, f"comparison_efficiency_ITS-TPC_{sign}_{var}")


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument(
        "o2_input_file", help="input O2 efficiency plot ROOT file pattern"
    )
    parser.add_argument(
        "ali_input_file", help="input AliPhysics efficiency plot ROOT file"
    )
    args = parser.parse_args()

    alifile = TFile(args.ali_input_file)
    #for var in ("Pt", "Eta"):
    var = "Pt"
    for sign in ("Positive", "Negative", "All"):
        o2file = TFile(f"{args.o2_input_file}_{sign}_{var}.root")
        canvo2 = o2file.Get(f"c_{sign}_all_ITS-TPC_{var}")
        compare_efficiency(alifile, canvo2, var, sign)


if __name__ == "__main__":
    main()
