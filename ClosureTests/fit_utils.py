"""
Fitting utils for the invariant mass script
"""

# pylint: disable=missing-function-docstring, invalid-name

from ROOT import TF1, TMath # pylint: disable=import-error
import ranges

def background(x, params):
    if x[0] > ranges.sig_range[0] and x[0] < ranges.sig_range[1]:
        TF1.RejectPoint()
        return 0
    return TMath.Exp(params[0] + params[1] * x[0])

def fit_sig_bkg(hist, doMC, mcSigmaSignal=-1, mcSigmaReflBkg=-1):
    params = {}
    estSigma1 = (ranges.sig_init_range[1] - ranges.sig_init_range[0]) / 6.
    estSigma2 = (ranges.sig_init_range[1] - ranges.sig_init_range[0]) / 6.
    if doMC:
        estSigma1 = mcSigmaSignal
        estSigma2 = mcSigmaReflBkg
        print(f"MC sigma from signal: {mcSigmaSignal:.3f} "
              f"reflected background: {mcSigmaReflBkg:.3f}")

    estMean = ranges.sig_init_range[0] + (ranges.sig_init_range[1] - ranges.sig_init_range[0]) / 2
    fitSig = TF1("fitSig", "gaus(0)+gaus(3)", ranges.sig_init_range[0], ranges.sig_init_range[1], 6)
    fitSig.SetParameters(20000, estMean, estSigma1, 20000, estMean, estSigma2)
    hist.Fit(fitSig, "NQ", "", ranges.sig_init_range[0], ranges.sig_init_range[1])
    params["scaleSig1"] = fitSig.GetParameter(0)
    params["meanSig1"] = fitSig.GetParameter(1)
    params["sigmaSig1"] = fitSig.GetParameter(2)
    params["scaleSig2"] = fitSig.GetParameter(3)
    params["meanSig2"] = fitSig.GetParameter(4)
    params["sigmaSig2"] = fitSig.GetParameter(5)
    ranges.sig_range[:] = [params["meanSig1"] - 3. * params["sigmaSig1"],
                           params["meanSig1"] + 3. * params["sigmaSig1"]]
    print(f'Signal parameters: scale: {params["scaleSig1"]:.3f} mean: {params["meanSig1"]:.3f} '
          f'sigma: {params["sigmaSig1"]:.3f}\nscale: {params["scaleSig2"]:.3f} '
          f'mean: {params["meanSig2"]:.3f} sigma: {params["sigmaSig2"]:.3f}\n'
          f'Initial ranges.signal region: {ranges.sig_init_range[0]:.2f}, '
          f'{ranges.sig_init_range[1]:.2f}, '
          f'final: {ranges.sig_range[0]:.2f}, {ranges.sig_range[1]:.2f}')

    fitBkg = TF1("fitBkg", background, ranges.pt_range_an[0], ranges.pt_range_an[1], 2)
    hist.Fit(fitBkg, "NQ+", "", ranges.pt_range_an[0], ranges.pt_range_an[1])
    params["offsetBkg"] = fitBkg.GetParameter(0)
    params["scaleBkg"] = fitBkg.GetParameter(1)
    print(f'Background parameters: offset: {params["offsetBkg"]:.3f} '
          f'scale: {params["scaleBkg"]:.3f}')

    fitSigBkg = TF1("fitSigBkg", "gaus(0)+gaus(3)+expo(6)",
                    ranges.pt_range_an[0], ranges.pt_range_an[1], 8)
    fitSigBkg.SetParameters(params["scaleSig1"], params["meanSig1"], params["sigmaSig1"],
                            params["scaleSig2"], params["meanSig2"], params["sigmaSig2"],
                            params["offsetBkg"], params["scaleBkg"])
    hist.Fit(fitSigBkg, "NQ+", "", ranges.pt_range_an[0], ranges.pt_range_an[1])
    params["scaleSigBkg1"] = fitSigBkg.GetParameter(0)
    params["meanSigBkg1"] = fitSigBkg.GetParameter(1)
    params["sigmaSigBkg1"] = fitSigBkg.GetParameter(2)
    params["scaleSigBkg2"] = fitSigBkg.GetParameter(0)
    params["meanSigBkg2"] = fitSigBkg.GetParameter(1)
    params["sigmaSigBkg2"] = fitSigBkg.GetParameter(2)
    params["offsetSigBkg"] = fitSigBkg.GetParameter(4)
    params["expScaleSigBkg"] = fitSigBkg.GetParameter(5)
    print(f'Joint fit parameters: scale: {params["scaleSigBkg1"]:.3f} '
          f'mean: {params["meanSigBkg1"]:.3f} '
          f'sigma: {params["sigmaSigBkg1"]:.3f}\nscale: {params["scaleSigBkg2"]:.3f} '
          f'mean: {params["meanSigBkg2"]:.3f} sigma: {params["sigmaSigBkg2"]:.3f}\n'
          f'offset: {params["offsetSigBkg"]:.3f} exp scale: {params["expScaleSigBkg"]:.3f}')

    return params

def fit_mc(hMassProjMC, init_range, fin_range):
    params = {}
    estSigma = (init_range[1] - init_range[0]) / 6.

    estMean = init_range[0] + (init_range[1] - init_range[0]) / 2
    fitSig = TF1("fitSig", "gaus", init_range[0], init_range[1], 6)
    fitSig.SetParameters(20000, estMean, estSigma)
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
