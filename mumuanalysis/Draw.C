#include <iostream>
#include <iomanip>

#include "TSystem.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"

using namespace std;

void Draw() {

  TString fileHistos, massrange, filename_prepend;
  Bool_t dualPlot, extractScalingFactor, drawInvmBox, drawLegend, drawInvm;
  Double_t integrated_lumi, integrated_lumiA, integrated_lumiB;

  integrated_lumiA = 2376.;
  integrated_lumiB = 2864.;

  integrated_lumi = 0.;
  integrated_lumi+= integrated_lumiA;
  integrated_lumi+= integrated_lumiB;

  // General plotting parameters
  Style_t mainFont = 43;
  Float_t legendSize = 16.;
  Float_t labelSize = 22.;
  Float_t axisTitleSize = 28;
  Float_t axisTitleOffset = 1.1;
  Float_t infoSize = 21.;

  drawInvm = false;

  drawInvmBox = drawInvm;
  drawLegend = drawInvm;

  //massrange = "lowm";
  //massrange = "Zpeak";
  //massrange = "highm";
  massrange = "full";

  filename_prepend = "ElasticSelection_";
  //filename_prepend = "ElasticSelection_runA_";
  //filename_prepend = "ElasticSelection_runB_";
  //filename_prepend = "InelasticSelection_";
  //filename_prepend = "NoSelection_";

  fileHistos = "XXX_EDITME_XXX.root";

  gSystem->Load("includes/GlobalFunctions_C.so");
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame 
  gStyle->SetPadTickY(1); 
  gStyle->SetEndErrorSize(2);
  

  dualPlot = true;
  extractScalingFactor = false;

  Double_t weight_elel = integrated_lumi/1e5*0.40892114; // August 2012 generation with 0<M<500 and pT>15
  Double_t weight_inelel = integrated_lumi/1e5*0.41048574*2; // August 2012 generation with pT>15
  Double_t weight_inelinel = integrated_lumi/1e5*0.61697446; // August 2012 generation with 0<M<500 and pT>15
  //Double_t weight_dymumu = integrated_lumi/27452498*1626.0; // FIXME no k factor
  Double_t weight_dymumu = integrated_lumi/25477423*1626.0;
  Double_t weight_dytautau = integrated_lumi/19937479*1627;
  Double_t weight_tt = integrated_lumi/10263328*149.6;
  Double_t weight_qcd_pt20to30 = integrated_lumi/10076800*2.363e8*0.00568;
  Double_t weight_ww_ct10 = integrated_lumi/539746*43.7;
  Double_t weight_ww_z2 = integrated_lumi/210667*2.927;
  Double_t weight_ww_jets = integrated_lumi/1197558*3.783;
  Double_t weight_qcd_muenriched = integrated_lumi/8797418*2.966e8*0.00118;
  Double_t weight_signal = integrated_lumi/5000*0.0407*(10.57e-2*10.57e-2)*2; //FIXME Jonathan's cross section value

  Double_t num_data(0);
  Double_t num_mc(0);
  Double_t num_elel(0);
  Double_t num_inelel(0);
  Double_t num_inelinel(0);
  Double_t num_others(0);
  Double_t num_data_evts(0);
  Double_t num_mc_evts(0);
  Double_t num_mc_evts_ds(0);

  Int_t numBins;
  Bool_t hasUnit;
  std::stringstream ss;
  TString xlabel, ylabel, hist_name, out_file;
  TH1D *hist_mc;
  Double_t error_mc, error_elel, error_inelel, error_inelinel;
  Double_t max_data, max_mc, maximum, minValue, maxValue, binSize;
  Double_t min, max, dist;
  Double_t num_data, num_mc, error_data, error_mc;
  
  ///////////////// QUANTITITES TO PLOT
    
  AddQuantity("acopl", "", "1-|#Delta#phi(#mu#mu)/#pi|", "");
  AddQuantity("Cuts", "", "", "");
  AddQuantity("numVtxAfterCuts", "numVtx_2", "Number of vertices in the event", "");
  AddQuantity("numVtx", "numVtx", "Number of vertices in the event", "");
  AddQuantity("vtxT", "vtxT", "", "");
  AddQuantity("vtxZ", "vtxZ", "", "");
  AddQuantity("pTpairZoom3", "pTpair-0to5", "p_{T}(#mu#mu)", "GeV");
  AddQuantity("invariantMass", "invariantMass", "m(#mu#mu)", "GeV");
  AddQuantity("etaSingleM", "etasinglem", "#eta(#mu^{-})", "");
  AddQuantity("etaSingleP", "etasinglep", "#eta(#mu^{+})", "");
  AddQuantity("etaPair", "etapair", "#eta(#mu#mu)", "");
  AddQuantity("numExtraTracks", "numExtraTracks", "Number of extra tracks on dimuon vertex", "", true);
  AddQuantity("dpt", "dpt-0to1", "#Delta p_{T}", "GeV");
  AddQuantity("pTSingleM", "ptsinglem", "p_{T}(#mu^{-})", "GeV");
  AddQuantity("pTSingleP", "ptsinglep", "p_{T}(#mu^{+})", "GeV");
  AddQuantity("acoplZoom2", "acopl-0to1", "1-|#Delta#phi(#mu#mu)/#pi|", "", true);
  AddQuantity("acoplZoom4", "acopl-0to0p1", "1-|#Delta#phi(#mu#mu)/#pi|", "", true);
  AddQuantity("pTpair", "pTpair", "p_{T}(#mu#mu)", "GeV", true);
  AddQuantity("pTpairZoom1", "pTpair-0to50", "p_{T}(#mu#mu)", "GeV");
  /*AddQuantity("etPtAfter", "", "p_{T} (extra tracks)", "GeV");
  AddQuantity("etEtaAfter", "", "#eta (extra tracks)", "");
  AddQuantity("etPt", "", "p_{T} (extra tracks)", "GeV");
  AddQuantity("etEta", "", "#eta (extra tracks)", "");
  AddQuantity("etPhi", "", "#phi (extra tracks)", "");
  AddQuantity("etSumPt", "", "#Sigma p_{T} (extra tracks)", "GeV");
  //AddQuantity("etDvtx", "", "Distance to primary vertex (extra tracks)", "cm");
  /*AddQuantity("pTSingleP", "", "p_{T}(#mu#mu)", "GeV");
    AddQuantity("etaSingleP", "", "#eta(#mu^{+})", "");*/

  /*AddQuantity("invariantMass", "", "m(#mu#mu)", "GeV");
  AddQuantity("pTSingleM", "", "p_{T}(#mu^{-})", "GeV");
  AddQuantity("pTSingleP", "", "p_{T}(#mu^{+})", "GeV");
  /*AddQuantity("pTpair", "", "p_{T}(#mu#mu)", "GeV");
  AddQuantity("dpt", "", "#Delta p_{T}", "GeV");
  AddQuantity("numExtraTracks", "", "Number of extra tracks on dimuon vertex", "");
  AddQuantity("numTrackerHitsM", "", "Number of hits in the tracker (#mu^{-})", "");
  AddQuantity("numTrackerHitsP", "", "Number of hits in the tracker (#mu^{+})", "");*/

  TFile *file = new TFile(fileHistos);
  TH1D *hist_tmp;

  //////////////// DATASETS TO PLOT

  AddDistribution("data", "Data", 1., 0, 0);
  AddDistribution("inelinel", "LPAIR #gamma#gamma#rightarrow#mu^{+}#mu^{-} (double-diss.)", weight_inelinel, 419);
  AddDistribution("inelel", "LPAIR #gamma#gamma#rightarrow#mu^{+}#mu^{-} (single-diss.)", weight_inelel, 30);
  AddDistribution("elel", "LPAIR #gamma#gamma#rightarrow#mu^{+}#mu^{-} (elastic)", weight_elel, 800);
  AddDistribution("ww-z2", "Inclusive WW (TuneZ2)", weight_ww_z2, 3);
  AddDistribution("tt", "t#bar{t}", weight_tt, 4);
  //AddDistribution("qcd", "QCD (#mu enriched)", weight_qcd_muenriched, 7);
  //AddDistribution("qcd-pt20to30", "QCD (#mu enriched, 30 GeV/c<p_{T}<50 GeV/c)", weight_qcd_pt20to30, 5);
  AddDistribution("ww-jets", "Madgraph W^{+}W^{-}", weight_ww_jets, 5);
  AddDistribution("dytautau", "PYTHIA Z2 DY #tau^{+}#tau^{-}", weight_dytautau, 6);
  AddDistribution("dymumu", "PYTHIA Z2 DY #mu^{+}#mu^{-} (Tune CT10)", weight_dymumu, 2);

  for (int q=0; q<quantities.size(); q++) {

    TLegend *leg = new TLegend(0.45,0.55,0.8,0.85);
    leg->SetTextFont(mainFont);
    leg->SetTextSize(legendSize);
    THStack *hs_q = new THStack(quantities.at(q), quantities.at(q));
    num_mc_evts = 0;
    num_data_evts = 0;
    num_elel = 0;
    num_inelel = 0;
    num_inelinel = 0;
    num_others = 0;
    Bool_t has_signal(false);
    Bool_t has_signal2(false);

    for (int d=0; d<dists.size(); d++) {
      hist_name = quantities.at(q)+"_"+massrange+"_"+dists.at(d);
      if (quantities.at(q)=="Cuts")
	hist_name = quantities.at(q)+"_"+dists.at(d);
      TH1D *hist = (TH1D*)(file->Get(hist_name));

      if (d==0) {
	hist_mc = (TH1D*)(hist->Clone());
	hist_mc->Scale(0.);
      }
      TString leg_style;


      if (distsIsMC.at(d)==2) {

	////////////////////////////////////////////////////////////////////////
	//                               SIGNAL                               //
	////////////////////////////////////////////////////////////////////////

	if (!extractScalingFactor) {
	  TH1D *hist_signal = (TH1D*)(hist->Clone());
	  hist_signal->Scale(distsWeights.at(d));
	  hist_signal->SetLineWidth(2);
	  hist_signal->SetLineColor(3);
	  hist = (TH1D*)(hist_signal->Clone());
	  
	  Double_t error;
	  num_mc_evts_ds = hist->IntegralAndError(0,hist->GetNbinsX()+1, error);
	  cout << "[Signal] Number of events for " << dists.at(d) << " : " << num_mc_evts_ds << " +/- " << error << endl;
	}
	leg_style="lpf";
      }
      else if (distsIsMC.at(d)!=1) {

        ////////////////////////////////////////////////////////////////////////
	//                                DATA                                //
	////////////////////////////////////////////////////////////////////////

	TH1D *hist_data = (TH1D*)(hist->Clone());	

	hist_data->SetMarkerStyle(20);
	hist_data->SetMarkerSize(.8);
	hist_data->SetLineColor(1);
	
	num_data_evts=hist_data->Integral(0, hist_data->GetNbinsX()+1);
	hist_data->Sumw2();
	leg_style="lpf";
	hist = (TH1D*)(hist_data->Clone());
	cout << "[DATA] Number of events for " << dists.at(d) << " : " << num_data_evts << endl;
      }
      else {

	////////////////////////////////////////////////////////////////////////
	//                                 MC                                 //
	////////////////////////////////////////////////////////////////////////

	hist->SetFillColor((Int_t)(distsColours.at(d)));
	hist->Scale(distsWeights.at(d));
	hist->SetLineStyle(1);
	hist->SetLineColor(1);

	num_mc_evts_ds = hist->IntegralAndError(0,hist->GetNbinsX()+1, error_mc);
	cout << "[ MC ] Number of events for " << dists.at(d) << " : " << num_mc_evts_ds << " +/- " << error_mc << endl;
	num_mc_evts+=num_mc_evts_ds;
	hs_q->Add(hist);
	hist_mc->Add(hist);
	leg_style="f";

	if (dists.at(d)=="elel") {
	  num_elel = num_mc_evts_ds;
	  error_elel = error_mc;
	  if (extractScalingFactor) {
	    TH1D *hist_signal = (TH1D*)(hist->Clone());
	    hist_signal->Scale(2.5);
	    hist_signal->SetFillColor(0);
	    hist_signal->SetFillStyle(0);
	    hist_signal->SetLineWidth(2);
	    hist_signal->SetLineColor(4);
	    if (hist_signal->Integral() != 0)
	      has_signal = true;
	    TH1D *hist_signal2 = (TH1D*)(hist->Clone());
	    hist_signal2->Scale(2.);
	    hist_signal2->SetFillColor(0);
	    hist_signal2->SetLineStyle(2);
	    hist_signal2->SetLineWidth(2);
	    hist_signal2->SetLineColor(3);
	    if (hist_signal2->Integral() != 0)
	      has_signal2 = true;
	  }
	}
	else if (dists.at(d)=="inelel") {
	  num_inelel = num_mc_evts_ds;
	  error_inelel = error_mc;
	}
	else if (dists.at(d)=="inelinel") {
	  num_inelinel = num_mc_evts_ds;
	  error_inelinel = error_mc;
	}
	else if (distsIsMC.at(d)==1) {
	  num_others += num_mc_evts_ds;
	}
      }


      if (hist->Integral() != 0)
	leg->AddEntry(hist, distsTitles.at(d), leg_style);
    }
    cout << "number of [data    ] events : " << num_data_evts << endl;
    cout << "number of [data    ] events : " << hist_data->GetSumOfWeights() << endl;
    cout << "number of [expected] events : " << num_mc_evts << endl;
    cout << "number of [expected] events : " << hist_mc->GetSumOfWeights() << endl;
    cout << "number of [data - background - elastic] events : " << num_data_evts-num_others-num_elel << endl;


    if (extractScalingFactor) {
      if (has_signal)
	leg->AddEntry(hist_signal, "LPAIR (elastic) * 2.5");
      if (has_signal2)
	leg->AddEntry(hist_signal2, "LPAIR (elastic) * 2");
    }      
    Double_t extra_factor = num_data_evts/num_mc_evts;
    if (num_elel!=0.)
      cout << "===========================\n"
	   << "Quantity -> " << quantities.at(q) << "\n"
	   << "===========================\n"
	   << "   data : " << num_data_evts << "\n"
	   << "     MC : " << num_mc_evts << "\n"
	   << "data/MC : " << extra_factor << "\n"
	   << "-----------------" << "\n"
	   << "  Elastic : " << num_elel << "\n"
	   << "Inelastic : " << (num_inelel+num_inelinel) << "\n"
	   << " ( Single : " << num_inelel << "\n"
	   << "   double : " << num_inelinel << " )\n"
	   << "--- Ratio : " << (num_inelel+num_inelinel+num_elel)/num_elel << " +/- " << sqrt(pow((num_inelel+num_inelinel)/num_elel,2)*pow(error_elel,2)+pow(error_inelel,2)+pow(error_inelinel,2))/num_elel << "\n"
	   << endl;


    ////////////////////////////////////////////////////////////////////////////
    //                              COMBINATION                               //
    ////////////////////////////////////////////////////////////////////////////

    leg->SetFillColor(0);
    leg->SetLineColor(0);

    TCanvas *c1 = new TCanvas(quantities.at(q), quantities.at(q), 800, 650);

    if (dualPlot) {
      c1->Divide(1,2);
      TPad *c1_1 = (TPad*)(c1->GetPad(1));
      c1_1->SetPad(0.,.250,1.,1.);
      c1_1->SetRightMargin(0.03);
      c1_1->SetBottomMargin(0.);
      TPad *c1_2 = (TPad*)(c1->GetPad(2));
      c1_2->SetPad(0.,0.,1.,.250);
      c1_2->SetBottomMargin(0.3);
      c1_2->SetRightMargin(0.03);
      c1_2->SetTopMargin(0.);
      c1->cd(1);
    }

    hist_mc->Draw("E2");

    hist_mc->SetTitle("");
    hist_mc->SetTitle("");
    xlabel = quantitiesNames.at(q);
    if (quantitiesUnits.at(q) != "")
      xlabel += " ["+quantitiesUnits.at(q)+"]";
    hist_mc->GetXaxis()->SetTitle(xlabel);
    minValue = hist_mc->GetXaxis()->GetBinLowEdge(hist_mc->GetXaxis()->GetFirst());
    maxValue = hist_mc->GetXaxis()->GetBinUpEdge(hist_mc->GetXaxis()->GetLast());

    numBins = hist_mc->GetNbinsX();
    binSize = (maxValue-minValue)/numBins;
    hasUnit = (quantitiesUnits.at(q) != "");

    ss.str(""); ss << setprecision(4) << binSize;
    ylabel = "Events";
    if (binSize!=1. || quantitiesUnits.at(q)!="") {
      ylabel += " / ";
      //if (hasUnit) ylabel += "(";
      ylabel += ss.str();
      if (hasUnit) ylabel += " "+quantitiesUnits.at(q)/*+")"*/;
    }

    hist_mc->GetYaxis()->SetTitle(ylabel);
    hist_mc->GetXaxis()->SetTitleOffset(axisTitleOffset*1.03);
    hist_mc->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hist_mc->GetXaxis()->SetLabelFont(mainFont);
    hist_mc->GetXaxis()->SetLabelSize(labelSize);
    hist_mc->GetXaxis()->SetTitleFont(mainFont);
    hist_mc->GetXaxis()->SetTitleSize(axisTitleSize);
    hist_mc->GetYaxis()->SetLabelFont(mainFont);
    hist_mc->GetYaxis()->SetLabelSize(labelSize);
    hist_mc->GetYaxis()->SetTitleFont(mainFont);
    hist_mc->GetYaxis()->SetTitleSize(axisTitleSize);

    max_data = hist_data->GetBinContent(hist_data->GetMaximumBin());
    max_mc = hist_mc->GetBinContent(hist_mc->GetMaximumBin());
    maximum = TMath::Max(max_data, max_mc);

    if (drawInvmBox) {
      TBox *box = new TBox(70., 0.1, 106., maximum*1.25);
      box->SetFillColor(kBlue);
      box->SetFillStyle(3003);
      box->SetLineColor(kBlue);
      box->SetLineWidth(2);
      box->Draw();
    }

    hs_q->Draw("SAME");

    hist_mc->SetMarkerStyle(0);
    hist_mc->SetFillColor(1);
    hist_mc->SetFillStyle(3004);
    hist_mc->Draw("E2SAME");

    if (extractScalingFactor) {
      if (has_signal)
	hist_signal->Draw("SAME");
      if (has_signal2)
	hist_signal2->Draw("SAME");
    }

    // Poisson error bars
    const double alpha = 1-0.6827;
    TGraphAsymmErrors * g = new TGraphAsymmErrors(hist_data);

    for (int i = 0; i < g->GetN(); ++i) {
      Int_t N = g->GetY()[i];
      Double_t L = (N==0) 
	? 0.
	: (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      Double_t U = (N==0) 
	? (ROOT::Math::gamma_quantile_c(alpha,N+1,1))
	: (ROOT::Math::gamma_quantile_c(alpha/2,N+1,1));
      if(N==0) U=0;
      g->SetPointEXlow(i, 0);
      g->SetPointEXhigh(i, 0);
      g->SetPointEYlow(i, N-L);
      g->SetPointEYhigh(i, U-N);
    }
    g->Draw("P SAME");

    max_data = hist_data->GetBinContent(hist_data->GetMaximumBin());
    max_mc = hist_mc->GetBinContent(hist_mc->GetMaximumBin());
    maximum = TMath::Max(max_data, max_mc);
    hist_mc->GetYaxis()->SetRangeUser(0.01,maximum*1.25);


    ss.str(""); ss << "CMS Preliminary 2011, #sqrt{s}=7 TeV, L="<< integrated_lumi/1000 << " fb^{-1}";
    TPaveText *plotlabel = new TPaveText(0.45,0.91,0.99,0.97,"NDC");
    plotlabel->SetTextColor(kBlack);
    plotlabel->SetFillColor(kWhite);
    plotlabel->SetBorderSize(0);
    plotlabel->SetTextAlign(32);
    plotlabel->SetTextSize(infoSize);
    plotlabel->SetTextFont(mainFont);
    plotlabel->AddText(ss.str().c_str());
    plotlabel->Draw();

    hist_mc->Draw("SAMEAXIS");

    if (drawLegend) leg->Draw();

    if (dualPlot) {
      c1->cd(2);
      TH1D *rap = (TH1D*)(hist_data->Clone());
      rap->Divide(hist_mc);

      TH1D *herror = (TH1D*)(hist_mc->Clone());
      TH1D *hmcerror = (TH1D*)(hist_mc->Clone());
      hmcerror->Sumw2();
      hmcerror->Divide(herror);

      num_data = hist_data->IntegralAndError(0, hist_data->GetNbinsX()+1, error_data);
      num_mc = hist_mc->IntegralAndError(0, hist_mc->GetNbinsX()+1, error_mc);
      cout << "Total number of events : \n"
	   << " -> Data : " << num_data << " +/- " << error_data << "\n"
	   << " ->   MC : " << num_mc << " +/- " << error_mc << "\n"
	   << " -> D/MC : " << num_data/num_mc << " +/- " << num_data/num_mc*sqrt(pow(error_data/num_data, 2)+pow(error_mc/num_mc, 2)) << endl;



      TLine *line = new TLine(rap->GetXaxis()->GetXmin(),1.,rap->GetXaxis()->GetXmax(),1.);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      line->SetLineStyle(2);
      
      rap->Draw();
      rap->GetXaxis()->SetTitle(xlabel);
      rap->GetYaxis()->SetTitle("Data/MC");
      rap->GetYaxis()->SetRangeUser(-2.,4.);
      rap->GetYaxis()->SetTickLength(.04);
      max = rap->GetBinContent(rap->GetMaximumBin())*1.15-1.;
      min = rap->GetBinContent(rap->GetMinimumBin())*1.15-1.;
      dist = TMath::Max(fabs(min),fabs(max));
      dist = TMath::Min(dist, 8.);
      rap->GetYaxis()->SetRangeUser(1.-dist,1.+dist);
      
      rap->GetXaxis()->SetTitleOffset(3.);
      rap->GetYaxis()->SetTitleOffset(axisTitleOffset);
      rap->GetXaxis()->SetLabelFont(mainFont);
      rap->GetXaxis()->SetLabelSize(labelSize);
      rap->GetXaxis()->SetTitleFont(mainFont);
      rap->GetXaxis()->SetTitleSize(axisTitleSize);
      rap->GetYaxis()->SetLabelFont(mainFont);
      rap->GetYaxis()->SetLabelSize(labelSize);
      rap->GetYaxis()->SetTitleFont(mainFont);
      rap->GetYaxis()->SetTitleSize(axisTitleSize);

      rap->SetTitle("");
      rap->SetStats(false);

      line->Draw();
      hmcerror->Draw("E2SAME");
      rap->Draw("SAME");

    }
    
    //canv->SetGrayscale();
    out_file = filename_prepend+quantitiesTitles.at(q)+"_"+massrange+"_1009";
    out_file = "tmp/"+out_file;
    c1->SaveAs(out_file+".pdf");
  }

}
