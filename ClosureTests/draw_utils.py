"""
Drawing utils for the invariant mass script
"""

# pylint: disable=missing-function-docstring, invalid-name

from ROOT import TF1, TLine, TCanvas, TGraph # pylint: disable=import-error
from ROOT import kFullCircle, kDashed, kBlack, kGreen, kBlue, kRed, kGray, kMagenta # pylint: disable=import-error

from fit_utils import gausn_gausn_expo, gausn_gausn_single_scale_expo, gausn_expo, gauss_root, gausn_root, expo

def draw_inv_mass(hist, pt_range, hist_range):
    #hist.GetXaxis().SetRangeUser(pt_d0_to_pik[0], pt_d0_to_pik[1])
    hist.GetXaxis().SetRangeUser(hist_range[0], hist_range[1])
    hist.SetTitle("%s #leq #it{p}_{T} < %s" % pt_range)
    hist.SetMarkerStyle(kFullCircle)
    hist.SetLineColor(kBlack)
    hist.SetMarkerColor(kBlack)
    #hist.SetMarkerSize(2)
    hist.Draw("PE")

def get_vertical_line(x, c):
    c.Update()
    yMin = c.GetUymin()
    yMax = c.GetUymax()
    lVert = TLine(x, yMin, x, yMax)
    lVert.SetLineColor(kBlack)
    lVert.SetLineStyle(kDashed)
    return lVert

def draw_fits(params, hist_range, sig_range):
    fits = {}

    fSigBkg = TF1("fSigBkg", gausn_gausn_single_scale_expo, hist_range[0], hist_range[1], 7)
    fSigBkg.SetParameters(params["scaleSigBkg1"], params["meanSigBkg1"], params["sigmaSigBkg1"],
                          params["meanSigBkg2"], params["sigmaSigBkg2"],
                          params["offsetSigBkg"], params["expScaleSigBkg"])
    fSigBkg.SetLineColor(kBlue)
    fSigBkg.Draw("same")
    fits["fSigBkg"] = fSigBkg

    fBkg = TF1("fBkg", expo, hist_range[0], hist_range[1], 2)
    fBkg.SetParameters(params["offsetBkg"], params["scaleBkg"])
    fBkg.SetLineColor(kRed)
    fBkg.Draw("same")
    fits["fBkg"] = fBkg

    real_sig_range = [params["meanSigBkg1"] - 3 * params["sigmaSigBkg1"],
                      params["meanSigBkg1"] + 3 * params["sigmaSigBkg1"]]
    fSig = TF1("fSig", gausn_root, real_sig_range[0], real_sig_range[1], 3)
    fSig.SetParameters(params["scaleSigBkg1"], params["meanSigBkg1"], params["sigmaSigBkg1"])
    fSig.SetLineColor(kGreen)
    fSig.Draw("same")
    fits["fSig"] = fSig

    refl_bkg_range = [params["meanSigBkg2"] - 3 * params["sigmaSigBkg2"],
                      params["meanSigBkg2"] + 3 * params["sigmaSigBkg2"]]
    fReflBkg = TF1("fReflBkg", gausn_root, refl_bkg_range[0], refl_bkg_range[1], 3)
    fReflBkg.SetParameters(params["scaleSigBkg1"], params["meanSigBkg2"], params["sigmaSigBkg2"])
    fReflBkg.SetLineColor(kMagenta)
    fReflBkg.Draw("same")
    fits["fReflBkg"] = fReflBkg

    return fits

def draw_mc_fit(params, sig_range):
    fSig = TF1("fSig", gausn_root, sig_range[0], sig_range[1], 3)
    fSig.SetParameters(params["scale"], params["mean"], params["sigma"])
    fSig.SetLineColor(kGreen)
    fSig.Draw("same")
    return fSig

def draw_mc_raw(mchist, hname, outfile):
    cdeb = TCanvas(hname, hname)
    cdeb.cd()
    mchist.Draw("colz")
    cdeb.Write(hname)
    cdeb.SaveAs(f"{outfile}_{hname}_raw.png")

def shade_signal(c, fits, sig_range):
    gr = TGraph()
    gr.SetFillColor(kGray)
    gr.SetFillStyle(3013)
    c.Update()
    ymin = c.GetUymin()
    ymax = c.GetUymax()
    xmin = sig_range[0]
    xmax = sig_range[1]

    npx = fits["fSigBkg"].GetNpx()
    npoints = 0
    dx = (xmax - xmin) / npx

    x = xmin + 0.5 * dx
    while x <= xmax:
        y = fits["fSigBkg"].Eval(x)
        if y < ymin:
            y = ymin
        if y > ymax:
            y = ymax
        gr.SetPoint(npoints, x, y)
        npoints += 1
        x += dx

    x = xmax - 0.5 * dx
    while x >= xmin:
        y = fits["fBkg"].Eval(x)
        if y < ymin:
            y = ymin
        if y > ymax:
            y = ymax
        gr.SetPoint(npoints, x, y)
        npoints += 1
        x -= dx

    return gr

def draw_more(sig_range, fits, c):
    lMin = get_vertical_line(sig_range[0], c)
    lMin.Draw("same")
    lMax = get_vertical_line(sig_range[1], c)
    lMax.Draw("same")
    gr = shade_signal(c, fits, sig_range)
    gr.Draw("f")
