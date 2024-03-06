'''
Analysis utilities for flow analysis
'''
import ROOT
import sys
import numpy as np

# Function to get project vn versus pt
# Input: thnsparse, pt_axis, vn_axis
# Output: th1d 
def get_vn_versus_pt(thnSparse, pt_axis, vn_axis):
    hist_vn_proj = thnSparse.Projection(vn_axis, pt_axis)
    hist_pt_proj = thnSparse.Projection(pt_axis)
    hist_pt_proj.Reset()
    for i in range(hist_pt_proj.GetNbinsX()):
        hist_sp_pt_proj = hist_vn_proj.ProjectionY(f'hist_sp_pt_proj_{i}', i+1, i+1)
        mean_sp = hist_sp_pt_proj.GetMean()
        mean_sp_err = hist_sp_pt_proj.GetMeanError()

        hist_pt_proj.SetBinContent(i+1, mean_sp)
        hist_pt_proj.SetBinError(i+1, mean_sp_err)
    return hist_pt_proj


def get_vn_versus_mass(thnSparse, inv_mass_bins, mass_axis, vn_axis, debug=False):
    '''
    Function to get project vn versus mass

    Input:
        - thnSparse:
            THnSparse, input THnSparse obeject (already projected in centrality and pt)
        - inv_mass_bins:
            list of floats, bin edges for the mass axis
        - mass_axis:
            int, axis number for mass
        - vn_axis:
            int, axis number for vn
        - debug:
            bool, if True, create a debug file with the projections (default: False)

    Output:
        - hist_mass_proj:
            TH1D, histogram with vn as a function of mass
    '''
    hist_vn_proj = thnSparse.Projection(vn_axis, mass_axis)
    hist_mass_proj = thnSparse.Projection(mass_axis)
    hist_mass_proj.Reset()
    invmass_bins = np.array(inv_mass_bins)
    hist_mass_proj = ROOT.TH1D('hist_mass_proj', 'hist_mass_proj', len(invmass_bins)-1, invmass_bins)
    if debug:
        outfile = ROOT.TFile('debug.root', 'RECREATE')
    for i in range(hist_mass_proj.GetNbinsX()):
        bin_low = hist_vn_proj.GetXaxis().FindBin(invmass_bins[i])
        bin_high = hist_vn_proj.GetXaxis().FindBin(invmass_bins[i+1])
        profile = hist_vn_proj.ProfileY(f'profile_{i}', bin_low, bin_high)
        mean_sp = profile.GetMean()
        mean_sp_err = profile.GetMeanError()
        #hist_sp_pt_proj = hist_vn_proj.ProjectionY(f'hist_sp_mass_proj_{bin_low}_{bin_high}', bin_low, bin_high)
        #if debug: hist_sp_pt_proj.Write()
        #mean_sp = hist_sp_pt_proj.GetMean()
        #mean_sp_err = hist_sp_pt_proj.GetMeanError()
        hist_mass_proj.SetBinContent(i+1, mean_sp)
        hist_mass_proj.SetBinError(i+1, mean_sp_err)
    
    if debug:
        hist_mass_proj.Write()
        outfile.Close()
    return hist_mass_proj

def get_resolution(resolution_file_name, dets, centMin, centMax, doEP=False, doAbs=False):
    '''
    Compute resolution for SP or EP method

    Input:
        - resolution_file_name: 
            str, path to the resolution file
        - dets: 
            list of strings, subsystems to compute the resolution
            if 2 subsystems are given, the resolution is computed as the product of the two
            if 3 subsystems are given, the resolution is computed as the product of the first two divided by the third
        - centMin:
            int, minimum centrality bin
        - centMax:
            int, maximum centrality bin
        - doEP:
            bool, if True, compute EP resolution
            if False, compute SP resolution (default)
        - doAbs:
            bool, if True, first apply absoulte value to the projections and then compute the resolution
            if False, do not apply absolute value to the projections (default)

    Output:
        - histo_reso:
            TH1D, histogram with the resolution value as a function of centrality
        - histo_dets:
            list of TH1D, list of projections used to compute the resolution
        - histo_means:
            list of TH1D, list of histograms with the mean value of the projections as a function of centrality
    '''
    if doEP:
        path = 'hf-task-flow-charm-hadrons/epReso/'
        prefix = 'EpReso'
    else:
        path = 'hf-task-flow-charm-hadrons/spReso/'
        prefix = 'SpReso'
    infile = ROOT.TFile(resolution_file_name, 'READ')

    histo_projs, histo_dets, histo_means = [], [], []
    detA = dets[0]
    detB = dets[1]
    if len(dets) == 3:
        detC = dets[2]
    dets = [f'{detA}{detB}', f'{detA}{detC}', f'{detB}{detC}'] if len(dets) == 3 else [f'{detA}{detB}']
   
    # collect the qvecs and the prepare histo for mean and resolution
    for det in dets:
        histo_dets.append(infile.Get(f'{path}h{prefix}{det}'))
        histo_dets[-1].SetDirectory(0)
        histo_dets[-1].SetName(f'h{prefix}{det}')
        histo_means.append(histo_dets[-1].ProjectionX(f'proj_{histo_dets[-1].GetName()}_mean'))
        histo_means[-1].SetDirectory(0)
        histo_means[-1].Reset()
        histo_projs.append([])

        # collect projections
        for cent in range(centMin, centMax):
            cent_min = cent
            cent_max = cent
            bin_cent_low = histo_dets[-1].GetXaxis().FindBin(cent_min) # common binning
            bin_cent_high = histo_dets[-1].GetXaxis().FindBin(cent_max)
            histo_projs[-1].append(histo_dets[-1].ProjectionY(f'proj_{histo_dets[-1].GetName()}_{cent_min}_{cent_max}', bin_cent_low, bin_cent_high))
            histo_projs[-1][-1].SetDirectory(0)

        # Appllying absolute value to the projections
        for ihist, histo in enumerate(histo_projs[-1]):
            if doAbs:
                histo_dum = histo.Clone()
                histo_projs[-1][ihist] = get_abs_histo(histo)
            histo_means[-1].SetBinContent(ihist+1, histo_projs[-1][ihist].GetMean())
    infile.Close()

    histo_reso = ROOT.TH1F('', '', 100, 0, 100)
    histo_reso.SetDirectory(0)
    for icent in range(centMin, centMax):
        histo_reso.SetBinContent(icent+1, compute_resolution([histo_means[i].GetBinContent(icent+1) for i in range(len(dets))]))

    if doAbs:
        for ihist, histo in enumerate(histo_dets):
            histo_dets[ihist] = get_absy_histo2(histo)

    return histo_reso, histo_dets, histo_means


def compute_resolution(subMean):
    '''
    Compute resolution for SP or EP method

    Input:
        - subMean:
            list of floats, list of mean values of the projections

    Output:
        - resolution:
            float, resolution value
    '''
    if len(subMean) == 1:
        resolution =  subMean[0]
        if resolution <= 0:
            return 0
        else:
            return np.sqrt(resolution)
    elif len(subMean) == 3:
        resolution = (subMean[2] * subMean[1]) / subMean[0] if subMean[0] != 0 else 0
        if resolution <= 0:
            return 0
        else:
            return np.sqrt(resolution)
    else:
        print('ERROR: dets must be a list of 2 or 3 subsystems')
        sys.exit(1)

def get_abs_histo(histo):
    '''
    Apply absolute value to the histogram

    Input:
        - histo:
            TH1D, input histogram
    Output:
        - histo:
            TH1D, histogram with absolute value of the bin content
    '''
    
    for ibin in range(histo.GetNbinsX()):
        if histo.GetBinCenter(ibin+1) < 0:
            histo.SetBinContent(histo.GetNbinsX() - ibin, histo.GetBinContent(ibin+1) + histo.GetBinContent(histo.GetNbinsX() - ibin))
            histo.SetBinContent(ibin+1, 0)

    return histo

def get_absy_histo2(histo):
    '''
    Apply absolute value to the histogram

    Input:
        - histo:
            TH1D, input histogram
    Output:
        - histo:
            TH1D, histogram with absolute value of the bin content
    '''
    
    for ibin in range(histo.GetNbinsX()):
        for jbin in range(histo.GetNbinsY()):
            if histo.GetYaxis().GetBinCenter(jbin+1) < 0:
                bincontent = histo.GetBinContent(ibin+1, jbin+1) + histo.GetBinContent(ibin+1, histo.GetNbinsY() - jbin)
                histo.SetBinContent(ibin+1, histo.GetNbinsY() - jbin, bincontent)
                histo.SetBinContent(ibin+1, jbin+1, 0)
    return histo

    

'''
    #hist_detA_detB = infile.Get(f'{path}h{prefix}{detA}{detB}')
    #hist_detA_detC = infile.Get(f'{path}h{prefix}{detA}{detC}')
    #hist_detB_detC = infile.Get(f'{path}h{prefix}{detB}{detC}')

    #bin_cent_low = hist_detA_detB.GetXaxis().FindBin(centMin) # common binning
    #bin_cent_high = hist_detA_detB.GetXaxis().FindBin(centMax)
    #proj_detA_detB = hist_detA_detB.ProjectionY(f'proj_detA_detB_{centMin}_{centMax}', bin_cent_low, bin_cent_high)
    #proj_detA_detC = hist_detA_detC.ProjectionY(f'proj_detA_detC_{centMin}_{centMax}', bin_cent_low, bin_cent_high)
    #proj_detB_detC = hist_detB_detC.ProjectionY(f'proj_detB_detC_{centMin}_{centMax}', bin_cent_low, bin_cent_high)
    resolution = (proj_detA_detB.GetMean() * proj_detA_detC.GetMean())/ proj_detB_detC.GetMean() if proj_detB_detC.GetMean() != 0 else 0

    if doAbs:
        proj_detA_detB_abs = proj_detA_detB.Clone()
        proj_detA_detC_abs = proj_detA_detC.Clone()
        proj_detB_detC_abs = proj_detB_detC.Clone()

        for i in range(proj_detA_detB_abs.GetNbinsX()):
            if proj_detA_detB_abs.GetBinCenter(i+1) < 0:
                proj_detA_detB_abs.SetBinContent(proj_detA_detB_abs.GetNbinsX() - i, proj_detA_detB_abs.GetBinContent(i+1) + proj_detA_detB_abs.GetBinContent(proj_detA_detB_abs.GetNbinsX() - i))
                proj_detA_detB_abs.SetBinContent(i+1, 0)
                proj_detA_detB_abs.SetBinError(i+1, 0)
        for i in range(proj_detA_detC_abs.GetNbinsX()):
            if proj_detA_detC_abs.GetBinCenter(i+1) < 0:
                proj_detA_detC_abs.SetBinContent(proj_detA_detC_abs.GetNbinsX() - i, proj_detA_detC_abs.GetBinContent(i+1) + proj_detA_detC_abs.GetBinContent(proj_detA_detC_abs.GetNbinsX() - i))
                proj_detA_detC_abs.SetBinContent(i+1, 0)
                proj_detA_detC_abs.SetBinError(i+1, 0)
        for i in range(proj_detB_detC_abs.GetNbinsX()):
            if proj_detB_detC_abs.GetBinCenter(i+1) < 0:
                proj_detB_detC_abs.SetBinContent(proj_detB_detC_abs.GetNbinsX() - i, proj_detB_detC_abs.GetBinContent(i+1) + proj_detB_detC_abs.GetBinContent(proj_detB_detC_abs.GetNbinsX() - i))
                proj_detB_detC_abs.SetBinContent(i+1, 0)
                proj_detB_detC_abs.SetBinError(i+1, 0)

        resolution = (proj_detA_detB_abs.GetMean() * proj_detA_detC_abs.GetMean())/ proj_detB_detC_abs.GetMean() if proj_detB_detC_abs.GetMean() != 0 else 0
        
            
    if resolution < 0:
        return 0
    resolution = np.sqrt(resolution)

    return resolution
'''

'''
def get_resolution_2sub(resolution_file_name, detA, detB, centMin, centMax, doEP=False):
    if doEP:
        path = 'hf-task-flow-charm-hadrons/epReso/'
        prefix = 'EpReso'
    else:
        path = 'hf-task-flow-charm-hadrons/spReso/'
        prefix = 'SpReso'
    infile = ROOT.TFile(resolution_file_name, 'READ')
    hist_detA_detB = infile.Get(f'{path}h{prefix}{detA}{detB}')
    bin_cent_low = hist_detA_detB.GetXaxis().FindBin(centMin) # common binning
    bin_cent_high = hist_detA_detB.GetXaxis().FindBin(centMax)
    proj_detA_detB = hist_detA_detB.ProjectionY(f'proj_detA_detB_{centMin}_{centMax}', bin_cent_low, bin_cent_high)
    resolution = np.sqrt(proj_detA_detB.GetMean()) if proj_detA_detB.GetMean() >= 0 else 0

    return resolution
'''

def compute_v2_oldstyle(thnSparse, deltaphiaxis):
    hist_cosDeltaPhi = thnSparse.Projection(deltaphiaxis)
    hist_cosDeltaPhi.GetXaxis().SetRangeUser(np.cos(-np.pi/4), np.cos(+np.pi/4))
    #hist_cosDeltaPhi.Draw()
    #input()
    Nin = hist_cosDeltaPhi.Integral()
    #hist_cosDeltaPhi.GetXaxis().SetRangeUser(-1, -0.924)
    #Nout_neg = hist_cosDeltaPhi.Integral()
    #print(f'Nin_pos = {Nin_pos}, Nout_neg = {Nout_neg}')
    #print(f'Nout = {Nout}')
    hist_cosDeltaPhi.GetXaxis().SetRangeUser(np.cos(np.pi/4), 1)
    Nout_pos = hist_cosDeltaPhi.Integral()
    hist_cosDeltaPhi.GetXaxis().SetRangeUser(-1, np.cos(-np.pi/4))
    Nout_neg = hist_cosDeltaPhi.Integral()
    Nout = Nout_pos + Nout_neg
    #print(f'Nin = {Nin}, Nout = {Nout}')
    vn = (np.pi * (Nin - Nout)) / (4 * (Nin + Nout)) if Nin+Nout != 0 else 0
    #print(f'vn = {vn}')
    return vn

def get_vn_versus_pt_oldstyle(thnSparse, pt_axis, deltaphi_axis):
    hist_pt_proj = thnSparse.Projection(pt_axis)
    hist_pt_proj.Reset()
    for i in range(hist_pt_proj.GetNbinsX()):
        pt_min = hist_pt_proj.GetXaxis().GetBinLowEdge(i+1)
        pt_max = hist_pt_proj.GetXaxis().GetBinUpEdge(i+1)
        thnSparse.GetAxis(pt_axis).SetRangeUser(pt_min, pt_max)
        vn = compute_v2_oldstyle(thnSparse, deltaphi_axis)
        hist_pt_proj.SetBinContent(i+1, vn)
    hist_pt_proj.Draw()
    input() 
    return hist_pt_proj
    

'''
def build_simfit_pdf(mass_min, mass_max, mass_axis, sp_axis, pt_axis, cent_axis):
    x = ROOT.RooRealVar("inv mass", "inv mass", mass_min, mass_max)
 
    # Construct signal pdf
    mean = ROOT.RooRealVar("mean", "mean", 1.86, 1.84, 1.88)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.01, 0.001, 0.1)
    gx = ROOT.RooGaussian("gx", "gx", x, mean, sigma)
 
    # Construct background pdf
    a0 = ROOT.RooRealVar("a0", "a0", -0.1, -1, 1)
    a1 = ROOT.RooRealVar("a1", "a1", 0.004, -1, 1)
    px = ROOT.RooChebychev("px", "px", x, [a0, a1])
 
    # Construct composite pdf
    f = ROOT.RooRealVar("f", "f", 0.2, 0.0, 1.0)
    model = ROOT.RooAddPdf("model", "model", [gx, px], [f])
 

    # Construct signal pdf.
    # NOTE that sigma is shared with the signal sample model
    mean_ctl = ROOT.RooRealVar("mean_ctl", "mean_ctl", -3, -8, 8)
    gx_ctl = ROOT.RooGaussian("gx_ctl", "gx_ctl", x, mean_ctl, sigma)
 
    # Construct the background pdf
    a0_ctl = ROOT.RooRealVar("a0_ctl", "a0_ctl", -0.1, -1, 1)
    a1_ctl = ROOT.RooRealVar("a1_ctl", "a1_ctl", 0.5, -0.1, 1)
    px_ctl = ROOT.RooChebychev("px_ctl", "px_ctl", x, [a0_ctl, a1_ctl])
 
    # Construct the composite model
    f_ctl = ROOT.RooRealVar("f_ctl", "f_ctl", 0.5, 0.0, 1.0)
    model_ctl = ROOT.RooAddPdf("model_ctl", "model_ctl", [gx_ctl, px_ctl], [f_ctl])
 
 
    # Construct combined dataset in (x,sample)
    combData = ROOT.RooDataSet(
        "combData",
        "combined data",
        {x},
        Index=sample,
        Import={"physics": data, "control": data_ctl},
    )
 
    # Construct a simultaneous pdf in (x, sample)
    # -----------------------------------------------------------------------------------
    
    # Construct a simultaneous pdf using category sample as index: associate model
    # with the physics state and model_ctl with the control state
    simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", {"physics": model, "control": model_ctl}, sample)
    
    # Perform a simultaneous fit
    # ---------------------------------------------------
    
    # Perform simultaneous fit of model to data and model_ctl to data_ctl
    fitResult = simPdf.fitTo(combData, PrintLevel=-1, Save=True)
    fitResult.Print()
    
    # Plot model slices on data slices
    # ----------------------------------------------------------------
    
    # Make a frame for the physics sample
'''