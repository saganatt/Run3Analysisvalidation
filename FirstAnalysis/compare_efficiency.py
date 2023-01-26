#!/usr/bin/env python
"""
file: compare_efficiency.py
brief: Plotting macro to compare O2 and AliPhysics efficiencies.
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
usage: ./compare_efficiency.py <O2 *.root file> <AliPhysics *.root file> <variable> <hadron>
       ./compare_efficiency.py efficiency_tracking_ITS-TPC_Positive_Pt.root
                               STE_RecPPRoverGenAcc_species_pt.root Pt All

For O2, efficiency_tracking_*.root files are produced by efficiency_studies macro and contain
a canvas with all particles respective efficiencies.

AliPhysics files were extracted separately from AliPhysics' correction framework.
"""

import argparse

from ROOT import TH1F, TCanvas, TLegend, TLatex # pylint: disable=import-error,no-name-in-module
#from ROOT import gPad, gStyle, gROOT # pylint: disable=import-error,no-name-in-module
from ROOT import TFile, gStyle, gROOT
from ROOT import kGreen, kOrange # pylint: disable=import-error,no-name-in-module


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
        hempty = TH1F(hname, f"{ctitle};Transverse Momentum (GeV/c);{ytitle}", 16, 0.00, 16)
        #gPad.SetLogx()
        return hempty
    def get_eta_hist():
        return TH1F(hname, f"{ctitle};Pseudorapidity;{ytitle}", 16, -1.5, 1.5)
    def get_phi_hist():
        return TH1F(hname, f"{ctitle};Azimuthal angle (rad);{ytitle}",
                    16, -2 * 3.1416 - 0.5, 2 * 3.1416 + 0.5)

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
                if (iso2 and count > 0) or not iso2:
                    lab = pelem.GetLabel()
                    if lab == "Charged part.":
                        lab = "All"
                    leg_labels.append(lab)
                count += 1

    typestr = "O2" if iso2 else "AliPhysics"
    if len(eff) != len(leg_labels):
        raise RuntimeError(f"Different number of plots and legend entries in {typestr} file")

    return eff, leg_labels


def compare_efficiency(canvali, canvo2, var): # pylint: disable=too-many-locals, too-many-statements
    """
    Compare O2 vs AliPhysics efficiency vs pT, eta, phi for all hadron species.
    """
    color_list = [1, 2, 4, kGreen + 2, kOrange - 3]
    ali_marker_list = [47, 21, 22, 34, 20] # full figures
    o2_marker_list = [46, 25, 26, 28, 24] # open figures

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

    results = prepare_canvas(var, get_titles(var, "Positive", "all", "ITS-TPC", "Efficiency"),
                             False)
    c_all = results[0]
    leg_all = results[1]
    #results_all_ratio = prepare_canvas(var, get_titles(var, "Positive", "all_ratio",
    #                                   "ITS-TPC", "Efficiency ratio"), False)
    #c_all_ratio = results_all_ratio[0]
    #leg_all_ratio = results_all_ratio[1]

    o2_eff, o2_leg_labels = retrieve_points(canvo2, True)
    ali_eff, ali_leg_labels = retrieve_points(canvali, False)

    if len(o2_eff) != len(ali_eff):
        raise RuntimeError("Different number of plots in AliPhysics and O2 files")

    for ind, (ali_elem, ali_lab, o2_elem, o2_lab) in \
            enumerate(zip(ali_eff, ali_leg_labels, o2_eff, o2_leg_labels)):
        ali_elem.SetMarkerStyle(ali_marker_list[ind])
        o2_elem.SetMarkerStyle(o2_marker_list[ind])

        c_all.cd()
        ali_elem.Draw(" same p")
        leg_all.AddEntry(ali_elem, ali_lab, "p")
        o2_elem.Draw(" same p")
        leg_all.AddEntry(o2_elem, o2_lab, "p")

        #if ind == 0:
            #c_all_ratio.cd()
            #max_val = ali_elem.GetXaxis().GetXmax()
            #hratio = ali_elem.Clone(f"ratio {o2_lab}")
            #for j in range(o2_elem.GetN() + 1):
            #    x_o2 = o2_elem.GetPointX(j)
            #    if x_o2 > max_val:
            #        break
            #    y_o2 = o2_elem.GetPointY(j)
            #    y_ali = ali_elem.GetBinContent(j + 1)
            #    hratio.SetBinContent(j + 1, y_o2 / y_ali)
            #hratio_max = hratio.GetYaxis().GetXmax() + 0.1
            #hratio.Draw(" same p")
            #hratio.GetYaxis().SetRangeUser(0., 1.5)
            #print(f"hratio y min: {hratio.GetYaxis().GetXmin()} "
                  #f"max: {hratio.GetYaxis().GetXmax()}")
            #leg_all_ratio.AddEntry(hratio, ali_lab, "p")

        results_single = prepare_canvas(var, get_titles(var, "Positive", o2_lab,
                                                        "ITS-TPC", "Efficiency"), True)
        c_single = results_single[0]
        leg_single = results_single[1]

        c_single.cd()
        ali_elem.Draw(" same p")
        leg_single.AddEntry(ali_elem, "Run 2", "p")
        o2_elem.Draw(" same p")
        leg_single.AddEntry(o2_elem, "Run 3", "p")

        leg_single.Draw()
        save_canvas(c_single, f"comparison_efficiency_ITS-TPC_Positive_{var}_{o2_lab}")

    c_all.cd()
    leg_all.Draw()
    save_canvas(c_all, f"comparison_efficiency_ITS-TPC_Positive_{var}")

    #c_all_ratio.cd()
    #leg_all_ratio.Draw()
    #save_canvas(c_all_ratio, f"comparison_efficiency_ratio_ITS-TPC_Positive_{var}")


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument("o2_input_file", help="input O2 efficiency plot ROOT file")
    parser.add_argument("ali_input_file", help="input AliPhysics efficiency plot ROOT file")
    parser.add_argument("var", default="Pt",
                        help="Capitalized name of efficiency variable (Pt/Eta)")
    args = parser.parse_args()

    o2file = TFile(args.o2_input_file)
    canvo2 = o2file.Get(f"c_all_ITS-TPC_Positive_{args.var}")
    alifile = TFile(args.ali_input_file)
    canvali = alifile.Get("c1")

    compare_efficiency(canvali, canvo2, args.var)

if __name__ == "__main__":
    main()
