import sys
import argparse
import ROOT
sys.path.append('..')
from utils.analysis_utils import get_resolution
from utils.style_utils import SetObjectStyle, SetFrameStyle, kRed, kFullCircle, SetGlobalStyle, LatLabel
SetGlobalStyle(padleftmargin=0.15, padlebottommargin=0.25,
               padrightmargin=0.15, titleoffsety=1.1, maxdigits=3, titlesizex=0.03,
               labelsizey=0.04, setoptstat=0, setopttitle=0)

ROOT.gROOT.SetBatch(False)

def compute_reso(an_res_file, doEP, doAbs, outputdir, suffix):

    detA = ['FT0c', 'FT0c', 'FT0c', 'FT0c']
    detB = ['FT0a', 'FT0a', 'FT0a', 'TPCpos']
    detC = ['FV0a', 'TPCpos', 'TPCneg', 'TPCneg']

    outfile_name = f'{outputdir}resoSP{suffix}.root' if not doEP else f'{outputdir}resoEP{suffix}.root'
    outfile_name = outfile_name.replace('.root', '_abs.root') if doAbs else outfile_name
    outfile = ROOT.TFile(outfile_name, 'RECREATE')
    for i, (det_a, det_b, det_c) in enumerate(zip(detA, detB, detC)):
        hist_reso, histos_det,\
        histo_means = get_resolution(an_res_file,
                                     [det_a, det_b, det_c],
                                     0,
                                     100,
                                     doEP,
                                     doAbs)
        hist_reso.SetDirectory(0)
        hist_reso.SetName(f'hist_reso_{det_a}_{det_b}_{det_c}')
        hist_reso.GetXaxis().SetTitle('centrality (%)')
        hist_reso.GetYaxis().SetTitle(f'#it{{R}}_{{2}}{{SP}} ({det_a}, {det_b}, {det_c})') if not doEP else hist_reso.GetYaxis().SetTitle(f'#it{{R}}_{{2}}{{EP}} ({det_a}, {det_b}, {det_c})')

        for i, (hist_det, hist_mean) in enumerate(zip(histos_det, histo_means)):
            outfile.cd()
            outfile.mkdir(f'{det_a}_{det_b}_{det_c}')
            outfile.cd(f'{det_a}_{det_b}_{det_c}')
            SetObjectStyle(hist_mean, color=kRed, markerstyle=kFullCircle,
                   markersize=1, fillstyle=0, linewidth=2)
            canv_name = hist_det.GetName().replace('h', 'canv')
            canvas = ROOT.TCanvas(canv_name, canv_name, 800, 800)
            canvas.SetLogz()
            hFrame = canvas.cd().DrawFrame(0, -2, 100, 2)
            SetFrameStyle(hFrame,
                          xtitle='Cent. FT0c (%)',
                          ytitle='cos(2(#Psi^{A}-#Psi^{B}))' if doEP else 'Q^{A} Q^{B}',
                          ytitleoffset=1.15,
                          ytitlesize=0.05,
                          ylabelsize=0.04,
                          ylabeloffset=0.01,
                          xticklength=0.04,
                          yticklength=0.03,
                          xtitlesize=0.05,
                          xlabelsize=0.04,
                          xtitleoffset=0.9,
                          xlabeloffset=0.020,
                          ydivisions=406,
                          xmoreloglabels=True,
                          ycentertitle=True,
                          ymaxdigits=5,)
            hist_det.Draw('same colz')
            hist_mean.Draw('same pl')
            if i == 0:
                LatLabel(f'A: {det_a}, B: {det_b}', 0.2, 0.85, 0.05)
            elif i == 1:
                LatLabel(f'A: {det_a}, B: {det_c}', 0.2, 0.85, 0.05)
            elif i == 2:
                LatLabel(f'A: {det_b}, B: {det_c}', 0.2, 0.85, 0.05)
            canvas.Write()
            hist_mean.Write()
            hist_det.Write()
        hist_reso.SetDirectory(outfile)
        hist_reso.Write()
        outfile.cd('..')

    input(f'Press enter to continue')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--doEP",  action="store_true", default=False,
                        help="do EP resolution")
    parser.add_argument("--doAbs",  action="store_true", default=False,
                        help="do absolute value of resolution")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    compute_reso(
        an_res_file=args.an_res_file,
        doEP=args.doEP,
        doAbs=args.doAbs,
        outputdir=args.outputdir,
        suffix=args.suffix
        )
