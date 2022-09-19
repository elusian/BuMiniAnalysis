

void roo_pulls(TCanvas& c, RooPlot* frame, const char* data, const char* fit, int chi2ndf = -1) {
    c.Divide(1, 2);
    c.GetPad(1)->SetLogy(c.GetLogy());
    c.GetPad(1)->SetPad(0.0, 0.25, 1, 1);
    c.GetPad(1)->SetBottomMargin(0.015);
    c.GetPad(1)->SetRightMargin(0.05);
//    c.GetPad(1)->SetTicks(1, 1);
    c.GetPad(2)->SetPad(0.0, 0.0, 1, 0.25);
    c.GetPad(2)->SetBottomMargin(0.32);
    c.GetPad(2)->SetTopMargin(0.0);
    c.GetPad(2)->SetRightMargin(0.05);
//    c.GetPad(2)->SetTicks(1, 1);
    
    c.GetPad(1)->SetBorderMode(0);
    c.GetPad(2)->SetBorderMode(0);
    
    auto pulls = frame->pullHist(data, fit, true); // true for useAverage
    pulls->SetMarkerStyle(frame->getAttMarker(data)->GetMarkerStyle());
    
    auto pullFrame = frame->emptyClone("pulls");
    auto line1 = RooFit::RooConst(3);
    auto linen1 = RooFit::RooConst(-3);
    auto line0 = RooFit::RooConst(0);
    line0.plotOn(pullFrame, RooFit::Name("line0"), RooFit::LineWidth(1), RooFit::LineColor(kRed));
    line1.plotOn(pullFrame, RooFit::Name("line1"), RooFit::LineWidth(1), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    linen1.plotOn(pullFrame, RooFit::Name("linen1"), RooFit::LineWidth(1), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    pullFrame->SetTitle("");
    pullFrame->addPlotable(pulls, "PE1");
    
    auto label_scale = 1.;
    
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->CenterTitle();
    pullFrame->GetXaxis()->SetTitleSize(0.18);
    pullFrame->GetYaxis()->SetTitleSize(0.18);
    pullFrame->GetYaxis()->SetTitleOffset(0.28);
    pullFrame->GetXaxis()->SetTitleOffset(.82);
    pullFrame->GetXaxis()->SetLabelSize(0.12 * label_scale);
    pullFrame->GetYaxis()->SetLabelSize(0.12 * label_scale);
    pullFrame->GetYaxis()->SetLabelOffset(0.006);
    
    frame->GetXaxis()->SetLabelOffset(0.006);
    frame->GetXaxis()->SetLabelSize(0.1);
    frame->GetYaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetTitleSize(0.06);
    
    if (pulls->GetMaximum() > 3.5 or pulls->GetMinimum() < -3.5)
    {
        pullFrame->SetMinimum(-5.5);
        pullFrame->SetMaximum(5.5);
    }
    else
    {
        pullFrame->SetMinimum(-3.5);
        pullFrame->SetMaximum(3.5);
    }
    
    c.GetPad(1)->cd();
    frame->Draw();
    
    if (chi2ndf >= 0) {
        Double_t chisquare = frame->chiSquare(fit, data, chi2ndf);
        
        int nParTest = 0;
        if (chi2ndf == 0) {
            nParTest = 1;
        }
        Double_t chisquareTest = frame->chiSquare(fit, data, nParTest);
        
        auto realChi2 = chisquare*chisquareTest*(chi2ndf-nParTest)/(chisquare - chisquareTest);
        auto nBins = (chisquare*chi2ndf - chisquareTest*nParTest)/(chisquare - chisquareTest);
        
        clog << "CHI2 TEST: " << realChi2 << "/" << (nBins - chi2ndf) << " = " << realChi2/(nBins - chi2ndf) << " vs " << chisquare << endl;
        clog << "NBINS: " << nBins << endl;
          
        auto& var = *frame->getPlotVar();

        auto range = var.getMax() - var.getMin();
        auto chi2 = new TLatex(var.getMin() + range/1000, frame->GetMaximum()*1.04, Form("#chi^{2}/NDF = %f", chisquare));
        chi2->Draw();
    }
    
    
    c.GetPad(2)->cd();
    pullFrame->Draw();
}
