#!/usr/bin/env python
"""
file: plot_tree.py
brief: Plotting most important distributions for rectangular cuts
usage: ./plot_tree.py AnalysisResults_tree_signal.root AnalysisResults_tree_bkg.root
author: Maja Kabus <mkabus@cern.ch>, CERN / Warsaw University of Technology
"""

import argparse
from array import array

# pylint: disable=import-error,no-name-in-module
from ROOT import (
    TCanvas,
    TFile,
    #TTree,
    #TList,
    TLegend,
    TH1F,
    gROOT,
    gStyle,
)


def save_canvas(canvas, title):
    """
    Save canvas in png, pdf, root.
    """
    format_list = [".png", ".pdf", ".root"]
    for file_format in format_list:
        canvas.SaveAs(title + file_format)


def prepare_hists(settings):
    """
    Initialize histograms.
    """
    hists_sig = { var: [] for var in settings["var_list"] }
    hists_bkg = { var: [] for var in settings["var_list"] }
    gStyle.SetOptStat(0)
    gStyle.SetFrameLineWidth(2)

    def prepare_hist_single(histname, hist_title):
        hist_sig = TH1F(f"h_sig_{histname}", f"{hist_title}", *var_range)
        hist_bkg = TH1F(f"h_bkg_{histname}", f"{hist_title}", *var_range)
        hist_sig.SetLineColor(1)
        hist_bkg.SetLineColor(2)
        hists_sig[var].append(hist_sig)
        hists_bkg[var].append(hist_bkg)

    for var, varu, var_range in zip(settings["var_list"], settings["var_list_u"],
                                    settings["var_ranges"]):
        if "pt_ranges" in settings:
            for i in range(len(settings["pt_ranges"]) - 1):
                histname = f'{varu}_pt_{settings["pt_ranges"][i]}-{settings["pt_ranges"][i + 1]}'
                hist_title = f'Normalized {var} for {settings["pt_ranges"][i]}' \
                             f' #leq #it{{p}}_{{T}} < {settings["pt_ranges"][i + 1]}'
                prepare_hist_single(histname, hist_title)
        else:
            prepare_hist_single(varu, f"Normalized {var}")

    return hists_sig, hists_bkg


def fill_hists(settings, settings_pt, infile, hists, hists_pt, sig_or_bkg, outfile): # pylint: disable=too-many-locals, too-many-arguments, too-many-branches
    """
    Fill histograms from trees.
    """
    tree_name = "O2hfcand3pfull"
    #tree_list = TList()
    #count = 0
    for key in infile.GetListOfKeys(): # pylint: disable=too-many-nested-blocks
        key_name = key.GetName()
        if key_name.startswith("DF_"): # is the dataframe directory
            tree = infile.Get(f"{key_name}/{tree_name}")
            #tree.SetName(f"merged_{sig_or_bkg}")
            #tree_list.Add(tree)
            #count = count + 1
            pt_val = array("f", [ 0. ])
            tree.SetBranchAddress("fPt", pt_val)
            mc_flag = array("b", [ 0 ])
            tree.SetBranchAddress("fMCflag", mc_flag)
            var_leaves = { var: array("f", [0.]) for var in settings["var_list"] }
            var_leaves_pt = { var: array("f", [0.]) for var in settings_pt["var_list"] }
            for var, leaf in zip(settings["var_list"], settings["leaf_list"]):
                tree.SetBranchAddress(leaf, var_leaves[var])
            for var, leaf in zip(settings_pt["var_list"], settings_pt["leaf_list"]):
                tree.SetBranchAddress(leaf, var_leaves_pt[var])
            for i in range(tree.GetEntries()):
                tree.GetEntry(i)
                # mc_flag == 1 << DecayType::LcToPKPi, decay type == 1
                if (sig_or_bkg == "sig" and (mc_flag[0] == 2 or mc_flag[0] == -2)) \
                        or (sig_or_bkg == "bkg" and mc_flag[0] == 0):
                    for j in range(len(settings["pt_ranges"]) - 1):
                        if settings["pt_ranges"][j] <= pt_val[0] < settings["pt_ranges"][j + 1]:
                            for var in var_leaves:
                                hists[var][j].Fill(var_leaves[var][0])
                            break
                    for var in var_leaves_pt:
                        hists_pt[var][0].Fill(var_leaves_pt[var][0])
            #if count == 5:
            #    break
    #tree_count = tree_list.GetEntries()
    #print(f"Merging {tree_count} trees")
    #out_tree = TTree.MergeTrees(tree_list)
    #outfile.cd()
    #print("Writing the merged tree")
    #out_tree.Write()

    #def cond(i):
    #    #return f'(("{sig_or_bkg}" == "sig"&& (fMCflag == 2 || fMCflag == -2)) ' \
    #    #            f'|| ("{sig_or_bkg}" == "bkg" && fMCflag == 0)) ' \
    #    return f'({settings["pt_ranges"][i]} <= fPt < {settings["pt_ranges"][i+1]})'
    for var, leaf in zip(settings["var_list"], settings["leaf_list"]):
        for i in range(len(settings["pt_ranges"]) - 1):
            #print(f'Plotting tree for {var} {settings["pt_ranges"][i]}')
            #cond_i = cond(i)
            #print(f"cond: {cond_i}")
            #cur_hist = hists[var][i].GetName()
            #tree_line = f"{leaf}>>+{cur_hist}"
            #print(f"Plotting {tree_line}")
            #out_tree.Draw(tree_line, cond(i))
            outfile.WriteObject(hists[var][i], hists[var][i].GetName())
    for var, leaf in zip(settings_pt["var_list"], settings_pt["leaf_list"]):
        outfile.WriteObject(hists_pt[var][0], hists_pt[var][0].GetName())


def plot_single(h_sig, h_bkg, filename):
    """
    Plot a single signal vs background comparison.
    """
    canv = TCanvas(f"c_{filename}", filename, 800, 600)
    canv.cd()
    canv.SetGridx()
    canv.SetGridy()

    int_sig = h_sig.Integral()
    int_bkg = h_bkg.Integral()
    if int_sig != 0.0:
        h_sig.Scale(1. / int_sig) # probability distribution, sum of content = 1.0
    if int_bkg != 0.0:
        h_bkg.Scale(1. / int_bkg) # probability distribution, sum of content = 1.0
    h_sig.Draw("hist")
    h_bkg.Draw("hist;same")
    y_min = min(h_sig.GetMinimum(), h_bkg.GetMinimum())
    y_max = max(h_sig.GetMaximum(), h_bkg.GetMaximum())
    margin = 0.05
    k = 1.0 - 2 * margin
    y_range = y_max - y_min
    h_sig.GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)
    h_bkg.GetYaxis().SetRangeUser(0.0, y_max + margin / k * y_range)

    legend = TLegend(0.50, 0.72, 0.70, 0.90)
    legend.AddEntry(h_sig, "signal", "L")
    legend.AddEntry(h_bkg, "background", "L")
    legend.Draw()

    save_canvas(canv, filename)


def plot(hists_sig, hists_bkg, settings):
    """
    Plot distributions for each variable.
    """
    for var, varu in zip(settings["var_list"], settings["var_list_u"]):
        for i, (h_sig, h_bkg) in enumerate(zip(hists_sig[var], hists_bkg[var])):
            filename = f'{varu}_pt_{settings["pt_ranges"][i]}-{settings["pt_ranges"][i + 1]}'
            plot_single(h_sig, h_bkg, filename)


def plot_file(outfile, settings):
    """
    Plot distributions from file.
    """
    for varu in settings["var_list_u"]:
        for i in range(len(settings["pt_ranges"]) - 1):
            histname = f'{varu}_pt_{settings["pt_ranges"][i]}-{settings["pt_ranges"][i + 1]}'
            h_sig = outfile.Get(f"h_sig_{histname}")
            h_bkg = outfile.Get(f"h_bkg_{histname}")
            plot_single(h_sig, h_bkg, histname)


def plot_pt(hists_sig, hists_bkg, settings):
    """
    Plot momenta distributions.
    """
    for var, varu in zip(settings["var_list"], settings["var_list_u"]):
        h_sig = hists_sig[var][0]
        h_bkg = hists_bkg[var][0]
        filename = varu
        plot_single(h_sig, h_bkg, filename)


def plot_pt_file(outfile, settings):
    """
    Plot momenta distributions from file.
    """
    for varu in settings["var_list_u"]:
        histname = varu
        h_sig = outfile.Get(f"h_sig_{histname}")
        h_bkg = outfile.Get(f"h_bkg_{histname}")
        plot_single(h_sig, h_bkg, histname)


def main():
    """
    Main function.
    """
    gROOT.SetBatch(True)
    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument("sig_input_file", help="input signal tree AnalysisResults_tree.root file")
    parser.add_argument("bkg_input_file", help="input bkg tree AnalysisResults_tree.root file")
    parser.add_argument("outfile", help="output file with all histograms saved for later plotting")
    parser.add_argument("--plot_file", default=False, action="store_true",
                        help="if set, plot histograms from file instead of calculating from trees")

    args = parser.parse_args()

    settings = { "var_list": ["decay length", "decay length XY", "CPA", "CPA XY",
                              "Chi2PCA", "mass",
                              "impact parameter 0", "impact parameter 1", "impact parameter 2"],
                 "var_list_u": ["decay_length", "decay_length_XY", "CPA", "CPA_XY",
                                "Chi2PCA", "mass",
                                "impact_parameter_0", "impact_parameter_1", "impact_parameter_2"],
                 "leaf_list": ["fDecayLength", "fDecayLengthXY", "fCPA", "fCPAXY",
                               "fChi2PCA", "fM",
                               "fImpactParameter0", "fImpactParameter1", "fImpactParameter2"],
                 "var_ranges": [[100, 0.0, 0.1], [100, 0.0, 0.1], [100, 0.9, 1.], [100, 0.9, 1.],
                                [200, 0., 0.01], [600, 1.98, 2.58],
                                [100, -0.1, 0.1], [100, -0.1, 0.1], [100, -0.1, 0.1]],
                 "pt_ranges": [0, 1, 2, 4, 6, 8, 12, 24]
                }
    settings_pt = { "var_list": ["pt", "pt prong0", "pt prong1", "pt prong2"],
                    "var_list_u": ["pt", "pt_prong0", "pt_prong1", "pt_prong2"],
                    "leaf_list": ["fPt", "fPtProng0", "fPtProng1", "fPtProng2"],
                    "var_ranges": [[200, 0, 24], [200, 0, 24], [200, 0, 24], [200, 0, 24]],
                  }

    if args.plot_file:
        outfile = TFile(args.outfile)
        plot_file(outfile, settings)
        plot_pt_file(outfile, settings_pt)
    else:
        infile_sig = TFile(args.sig_input_file)
        infile_bkg = TFile(args.bkg_input_file)

        hists_sig, hists_bkg = prepare_hists(settings)
        hists_sig_pt, hists_bkg_pt = prepare_hists(settings_pt)

        print("Filling histos")
        outfile = TFile(args.outfile, "RECREATE")
        fill_hists(settings, settings_pt, infile_sig, hists_sig, hists_sig_pt, "sig", outfile)
        print("Filled histos for signal")
        fill_hists(settings, settings_pt, infile_bkg, hists_bkg, hists_bkg_pt, "bkg", outfile)
        print("Filled histos for background")

        plot(hists_sig, hists_bkg, settings)
        plot_pt(hists_sig_pt, hists_bkg_pt, settings_pt)

        infile_sig.Close()
        infile_bkg.Close()

    outfile.Close()


if __name__ == "__main__":
    main()
