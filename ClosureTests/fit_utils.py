"""
Fitting utils for the invariant mass script
"""

# pylint: disable=missing-function-docstring, invalid-name

from ROOT import TF1, TMath # pylint: disable=import-error
import ranges

pdg_mass = 1.865

def expo(x, params):
    return TMath.Exp(params[0] + params[1] * x[0])

def gauss_root(x, params):
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2])**2)

def gausn_root(x, params):
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2])**2) /\
            (params[2] * TMath.Sqrt(TMath.TwoPi()))

def gausn_gausn(x, params):
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2])**2) /\
            (params[2] * TMath.Sqrt(TMath.TwoPi())) +\
            params[3] * TMath.Exp(-0.5 * ((x[0] - params[4]) / params[5])**2) /\
            (params[5] * TMath.Sqrt(TMath.TwoPi()))

def gausn_expo(x, params):
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2])**2) /\
            (params[2] * TMath.Sqrt(TMath.TwoPi())) +\
            TMath.Exp(params[3] + params[4] * x[0])

def gausn_gausn_expo(x, params):
    return params[0] * TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2])**2) /\
            (params[2] * TMath.Sqrt(TMath.TwoPi())) +\
            params[3] * TMath.Exp(-0.5 * ((x[0] - params[4]) / params[5])**2) /\
            (params[5] * TMath.Sqrt(TMath.TwoPi())) +\
            TMath.Exp(params[6] + params[7] * x[0])

def gausn_gausn_single_scale_expo(x, params):
    return params[0] * (TMath.Exp(-0.5 * ((x[0] - params[1]) / params[2])**2) /\
            (params[2] * TMath.Sqrt(TMath.TwoPi())) +\
            TMath.Exp(-0.5 * ((x[0] - params[3]) / params[4])**2) /\
            (params[4] * TMath.Sqrt(TMath.TwoPi()))) +\
            TMath.Exp(params[5] + params[6] * x[0])

def gausn_gausn_more_scales_expo(x, params):
    return params[0] * (params[1] * TMath.Exp(-0.5 * ((x[0] - params[2]) / params[3])**2) /\
            (params[3] * TMath.Sqrt(TMath.TwoPi())) +\
            params[4] * TMath.Exp(-0.5 * ((x[0] - params[5]) / params[6])**2) /\
            (params[6] * TMath.Sqrt(TMath.TwoPi()))) +\
            TMath.Exp(params[7] + params[8] * x[0])

def background(x, params):
    if x[0] >= ranges.sig_range[0] and x[0] <= ranges.sig_range[1]:
        TF1.RejectPoint()
        return 0
    return TMath.Exp(params[0] + params[1] * x[0])

def fit_sig_bkg(hist, sig_range, full_range, doMC, mc_sig_params, mc_bkg_params, mc_yields):
    params = {}

    for i in range(hist.GetNbinsX()):
        if hist.GetBinCenter(i) > ranges.pt_range_an[0] and hist.GetBinCenter(i) < ranges.pt_range_an[1]:
            unc = hist.GetBinError(i)
            pt = hist.GetBinCenter(i)
            print(f"Uncertainty for bin {i}, pt {pt:.3f}: {unc:.3f}")

    fitBkg = TF1("fitBkg", background, full_range[0], full_range[1], 2)
    hist.Fit(fitBkg, "NQ", "", full_range[0], full_range[1])
    params["offsetBkg"] = fitBkg.GetParameter(0)
    params["scaleBkg"] = fitBkg.GetParameter(1)
    print(f'Background parameters: offset: {params["offsetBkg"]:.3f} '
          f'scale: {params["scaleBkg"]:.3f}\n'
          f'range: {full_range[0]:.2f} {full_range[1]:.2f}')

    if doMC:
        estSigma1 = mc_sig_params["sigma"]
        estSigma2 = mc_bkg_params["sigma"]
        estMean1 = pdg_mass
        estMean2 = pdg_mass
    else:
        estSigma1 = (sig_range[1] - sig_range[0]) / 6.
        estSigma2 = (sig_range[1] - sig_range[0]) / 6.
        estMean1 = sig_range[0] + (sig_range[1] - sig_range[0]) / 2.
        estMean2 = sig_range[0] + (sig_range[1] - sig_range[0]) / 2.

    fitSigBkg = TF1("fitSigBkg", gausn_gausn_single_scale_expo,
                    full_range[0], full_range[1], 7)
    fitSigBkg.FixParameter(0, mc_sig_params["scale"] / mc_bkg_params["scale"])
    fitSigBkg.SetParameter(1, estMean1)
    fitSigBkg.SetParameter(2, estSigma1)
    fitSigBkg.SetParameter(3, estMean2)
    fitSigBkg.FixParameter(4, estSigma2)
    fitSigBkg.SetParameter(5, params["offsetBkg"])
    fitSigBkg.SetParameter(6, params["scaleBkg"])
    # Default: chi-square method, L: log likelihood method, when histogram represents counts
    hist.Fit(fitSigBkg, "NQ+", "", full_range[0], full_range[1])
    params["scaleSigBkg1"] = fitSigBkg.GetParameter(0)
    params["meanSigBkg1"] = fitSigBkg.GetParameter(1)
    params["sigmaSigBkg1"] = fitSigBkg.GetParameter(2)
    params["meanSigBkg2"] = fitSigBkg.GetParameter(3)
    params["sigmaSigBkg2"] = fitSigBkg.GetParameter(4)
    params["offsetSigBkg"] = fitSigBkg.GetParameter(5)
    params["expScaleSigBkg"] = fitSigBkg.GetParameter(6)
    print(f'Joint fit parameters: '
          f'scale: {params["scaleSigBkg1"]:.3f} '
          f'mean: {params["meanSigBkg1"]:.3f} '
          f'sigma: {params["sigmaSigBkg1"]:.3f}'
          f'\nmean: {params["meanSigBkg2"]:.3f} sigma: {params["sigmaSigBkg2"]:.3f}'
          f'\noffset: {params["offsetSigBkg"]:.3f} exp scale: {params["expScaleSigBkg"]:.3f}'
          f'\nrange: {full_range[0]:.2f} {full_range[1]:.2f}')

    return params

def fit_mc(hMassProjMC, init_range, fin_range, init_scale):
    params = {}
    estSigma = (init_range[1] - init_range[0]) / 6.

    estMean = init_range[0] + (init_range[1] - init_range[0]) / 2
    fitSig = TF1("fitSig", gausn_root, init_range[0], init_range[1], 4)
    fitSig.SetParameters(init_scale, estMean, estSigma) # setting scale does not affect the results for pT [6, 8)
    #fitSig.SetParameter(1, estMean)
    #fitSig.SetParameter(2, estSigma)
    hMassProjMC.Fit(fitSig, "NQ", "", init_range[0], init_range[1])
    params["scale"] = fitSig.GetParameter(0)
    params["mean"] = fitSig.GetParameter(1)
    params["sigma"] = fitSig.GetParameter(2)
    fin_range[:] = [params["mean"] - 3. * params["sigma"], params["mean"] + 3. * params["sigma"]]
    print(f'MC parameters: scale: {params["scale"]:.3f} mean: {params["mean"]:.3f} '
          f'sigma: {params["sigma"]:.3f}\n'
          f'Initial region: {init_range[0]:.2f}, {init_range[1]:.2f}, '
          f'final: {fin_range[0]:.2f}, {fin_range[1]:.2f}')
    return params
