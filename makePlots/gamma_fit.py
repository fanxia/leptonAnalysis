import ROOT

def makeFit(varname, varmin, varmax, signalHist, backgroundHist, dataHist, plotName):

    # RooFit variables
    var = ROOT.RooRealVar(varname, varname, varmin, varmax)
    argList = ROOT.RooArgList()
    argList.add(var)
    argSet = ROOT.RooArgSet()
    argSet.add(var)

    # create PDFs
    signalDataHist = ROOT.RooDataHist('signalDataHist', 'signal RooDataHist', argList, signalHist)
    signalPdf = ROOT.RooHistPdf('signalPdf', varname+' of signal', argSet, signalDataHist)

    backgroundDataHist = ROOT.RooDataHist('backgroundDataHist', 'background RooDataHist', argList, backgroundHist)
    backgroundPdf = ROOT.RooHistPdf('backgroundPdf', varname+' of background', argSet, backgroundDataHist)

    # data
    dataDataHist = ROOT.RooDataHist('data '+varname, varname+' in Data', argList, dataHist)

    # signal fraction parameter
    signalFractionVar = ROOT.RooRealVar('signal fraction', 'signal fraction', 0.5, 0.0, 1.0)
    sumPdf = ROOT.RooAddPdf('totalPdf', 'signal and background', signalPdf, backgroundPdf, signalFractionVar)

    # fit
    sumPdf.fitTo(dataDataHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1))

    print 'fit returned value ', signalFractionVar.getVal(), ' +/- ', signalFractionVar.getError()
    return (signalFractionVar.getVal(), signalFractionVar.getError())

openfiles = {}

def get1DHist(filename, histname):
    if filename not in openfiles:
        openfiles[filename] = ROOT.TFile(filename,'READ')
    file = openfiles[filename]

    hist = ROOT.TH1D()
    hist = file.Get(histname)
    hist.SetDirectory(0)
    hist.SetFillColor(0)
    #hist.Sumw2()
    return hist

def doFit(typeToFit):

    varToFit = 'w_mT'
    inFile = 'mcPlots_'+typeToFit+'_bjj.root'

    ttgammaHist = get1DHist(inFile, varToFit+'_ttA_2to5_'+typeToFit+'_bjj')
    
    ttbarHist = get1DHist(inFile, varToFit+'_ttJetsHadronic_'+typeToFit+'_bjj')
    ttbarHist.Add(get1DHist(inFile, varToFit+'_ttJetsSemiLep_'+typeToFit+'_bjj'))
    ttbarHist.Add(get1DHist(inFile, varToFit+'_ttJetsFullLep_'+typeToFit+'_bjj'))

    DataHist = get1DHist(inFile, varToFit+'_gg_'+typeToFit+'_bjj')
    DataHist.Add(get1DHist(inFile, varToFit+'_TBar_s_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_TBar_t_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_TBar_tW_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_T_s_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_T_t_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_T_tW_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_dyJetsToLL_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_WJetsToLNu_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_WW_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_WZ_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_ZZ_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_TTWJets_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_TTZJets_'+typeToFit+'_bjj'), -1.0)
    DataHist.Add(get1DHist(inFile, varToFit+'_qcd_'+typeToFit+'_bjj'), -1.0)

    (ttgammaFrac, ttgammaFracError) = makeFit(varToFit, 70.0, 500.0, ttgammaHist, ttbarHist, DataHist, varToFit+'_fit.png')

    dataInt = DataHist.Integral()
    ttbarInt = ttbarHist.Integral()
    ttgammaInt = ttgammaHist.Integral()
    ttgammaSF = ttgammaFrac * dataInt / ttgammaInt
    ttgammaSFerror = ttgammaFracError * dataInt / ttgammaInt
    print '#'*80
    print 'Correction to the ttgamma scale factor: ', ttgammaSF, ' +-', ttgammaSFerror, '(fit error only)'
    ttbarSF = (1.0-ttgammaFrac) * dataInt / ttbarInt
    ttbarSFerror = ttgammaFracError * dataInt / ttbarInt
    print 'Correction to ttbar scale factor: ', ttbarSF, ' +-', ttbarSFerror, '(fit error only)'
    print '#'*80

    DataHist.Draw('e1')
    ttgammaHist.SetLineColor(8)
    ttgammaHist.Draw('hist same')
    ttbarHist.Draw('hist same')

doFit('ele')

print ' '*80
print '-'*80
print ' '*80

doFit('muon')
