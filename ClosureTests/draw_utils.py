"""
Drawing utils for the invariant mass script
"""

# pylint: disable=missing-function-docstring, invalid-name

from ROOT import TF1, TLine, TCanvas, TGraph # pylint: disable=import-error
from ROOT import kFullCircle, kDashed, kBlack, kGreen, kBlue, kRed, kGray # pylint: disable=import-error

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
    fBkg = TF1("fBkg", "expo", hist_range[0], hist_range[1])
    fBkg.SetParameters(params["offsetBkg"], params["scaleBkg"])
    fBkg.SetLineColor(kRed)
    fBkg.Draw("same")
    fits["fBkg"] = fBkg

    fSig = TF1("fSig", "gaus(0)+gaus(3)", sig_range[0], sig_range[1], 6)
    fSig.SetParameters(params["scaleSig1"], params["meanSig1"], params["sigmaSig1"],
                       params["scaleSig2"], params["meanSig2"], params["sigmaSig2"])
    fSig.SetLineColor(kGreen)
    fSig.Draw("same")
    fits["fSig"] = fSig

    fSigBkg = TF1("fSigBkg", "gaus(0)+gaus(3)+expo(6)", hist_range[0], hist_range[1], 8)
    fSigBkg.SetParameters(params["scaleSigBkg1"], params["meanSigBkg1"], params["sigmaSigBkg1"],
                          params["scaleSigBkg2"], params["meanSigBkg2"], params["sigmaSigBkg2"],
                          params["offsetSigBkg"], params["expScaleSigBkg"])
    fSigBkg.SetLineColor(kBlue)
    fSigBkg.Draw("same")
    fits["fSigBkg"] = fSigBkg

    return fits

def draw_mc_fit(params, sig_range):
    fSig = TF1("fSig", "gaus", sig_range[0], sig_range[1], 6)
    fSig.SetParameters(params["scale"], params["mean"], params["sigma"])
    fSig.SetLineColor(kGreen)
    fSig.Draw("same")

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

    npx = fits["fSig"].GetNpx()
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
    lMin = draw_utils.get_vertical_line(sig_range[0], c)
    lMin.Draw("same")
    lMax = draw_utils.get_vertical_line(sig_range[1], c)
    lMax.Draw("same")
    gr = draw_utils.shade_signal(c, fits, sig_range)
    gr.Draw("f")
