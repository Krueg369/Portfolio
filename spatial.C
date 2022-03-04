
// This program creates a plot of the angular difference between double muon tracks

// Define Libraries
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include "TMath.h"
#include <cmath>
#define _USE_MATH_DEFINES
#include "TAttAxis.h"

using namespace std;

void spatial() {

  // use TChain just to get list of files in directory
  TChain chain;
  std::string infilename;
  infilename = "muons-*.root";
  chain.Add(infilename.c_str());

  TObjArray* filelist = chain.GetListOfFiles();
  TIter iter(filelist);


  // make histogram
  TH1D *angleHisto = new TH1D("Legend","Angular Difference Between Double Muon Events Over Root 2; Degrees; Frequency",140,0.0,5.0);
  TH1D *angle2Histo = new TH1D("Legend","68% of Angular Difference Histogram; Degrees; Frequency",100,0.0,0.875);
  TH1D *chiHisto = new TH1D("Legend","chi2",100,-1.0,70.0);
  TH1D *ndigitHisto = new TH1D("Legend","ndigit",100,-1.0,300);
  TH1D *lengthHisto = new TH1D("lengthHisto","length per track",60,0.0,60.0);
  TH1D *distHisto = new TH1D("distHisto","ditance between ntracks > 2",60,0.0,30.0);
  TH1D *dist2Histo = new TH1D("dist2Histo","ditance between ntracks == 2",60,0.0,30.0);

  // Loop over all files in chain element list until done selecting those
  TChainElement* element = 0;
  while ((element = dynamic_cast<TChainElement*>(iter.Next())) )
    {
      cout<<" new input file "<<element->GetTitle()<<endl;

      TFile* CondenseFile = new TFile(element->GetTitle());
      if ( !CondenseFile->IsOpen() || CondenseFile->IsZombie() )
        {
          cerr << "Warning! Failed to open file " << CondenseFile->GetName() << endl;
          delete CondenseFile;
          continue;
        }

      // make a tree reader instance: get the "data" tree out of CondenseFile
      // This is the root6 magic.  It's gross in root5
      TTreeReader muonReader("data", CondenseFile);

      // Let's look at number of tracks.  Things where there are one value
      // per event are accessed like this

      TTreeReaderValue<Int_t> ntrack(muonReader, "ntrack");
      TTreeReaderArray<UShort_t> ndigit(muonReader, "ndigit");
      TTreeReaderArray<Double_t> dcosz(muonReader, "dcosz");
      TTreeReaderArray<Double_t> dcosy(muonReader, "dcosy");
      TTreeReaderArray<Double_t> dcosx(muonReader, "dcosx");
      TTreeReaderArray<Double_t> chi2(muonReader, "chi2");
      TTreeReaderArray<Double_t> length(muonReader, "length");
      TTreeReaderArray<Double_t> vtxx(muonReader, "vtxx");
      TTreeReaderArray<Double_t> vtxy(muonReader, "vtxy");
      TTreeReaderArray<Double_t> vtxz(muonReader, "vtxz");
      // Loop over all entries of the TTree

      Double_t dist;
      Double_t dist2;


      while (muonReader.Next()) {
	        // TTreeReaderValue's are accessed with the "*" dereference
          UChar_t i;
          UChar_t j;
          UChar_t h;
          if(*ntrack > 2){

            for (i = 0; i < *ntrack; i++ ) {

              TVector3 a1(dcosx[i],dcosy[i],dcosz[i]);
              TVector3 q1(vtxx[i], vtxy[i],vtxz[i]);

              for(j = i; j < *ntrack; j++){
                TVector3 b2(dcosx[j],dcosy[j],dcosz[j]);
                TVector3 q2(vtxx[j], vtxy[j],vtxz[j]);

                if (i != j){
                  TVector3 l;
                  TVector3 diff(q1-q2);
                  l = a1.Cross(b2);
                  Double_t k = l.Mag();

                  dist = abs((l.Dot(diff))/k);


                  distHisto->Fill(dist);
                }

              }
            }

          }

          if (*ntrack == 2){

            TVector3 v1(dcosx[0],dcosy[0],dcosz[0]);
            TVector3 v2(dcosx[1],dcosy[1],dcosz[1]);
            TVector3 p1(vtxx[0], vtxy[0],vtxz[0]);
            TVector3 p2(vtxx[1], vtxy[1],vtxz[1]);
            TVector3 n;
            TVector3 dif(p1-p2);
            n = v1.Cross(v2);
            Double_t m = n.Mag();

            dist2 = abs((n.Dot(dif))/m);

            dist2Histo->Fill(dist2);


            Double_t angle = v1.Angle(v2);
            Double_t thet = angle * 180.0 * M_1_PI ;
            Double_t theta = thet * M_SQRT1_2 ;
            angleHisto->Fill(theta);
            angle2Histo->Fill(theta);


          }

      }

      // free input tree stuff, close input file
      if (CondenseFile) { delete CondenseFile; CondenseFile = 0; } // also destructs input file tree


    }

    //Histogram style
    TStyle* minosStyle = new  TStyle("minosStyle", "MINOS Style");

    //set the background color to white
    minosStyle->SetFillColor(10);
    minosStyle->SetFrameFillColor(10);
    minosStyle->SetCanvasColor(10);
    minosStyle->SetPadColor(10);
    minosStyle->SetTitleFillColor(0);
    minosStyle->SetStatColor(10);

    //dont put a colored frame around the plots
    minosStyle->SetFrameBorderMode(0);
    minosStyle->SetCanvasBorderMode(0);
    minosStyle->SetPadBorderMode(0);

    //use the primary color palette
    minosStyle->SetPalette(1,0);

    //set the default line color for a histogram to be black
    minosStyle->SetHistLineColor(kBlack);

    //set the default line color for a fit function to be red
    minosStyle->SetFuncColor(kRed);

    //make the axis labels black
    minosStyle->SetLabelColor(kBlack,"xyz");

    //set the default title color to be black
    minosStyle->SetTitleColor(kBlack);

    //set the margins
    minosStyle->SetPadBottomMargin(0.2);
    minosStyle->SetPadTopMargin(0.075);
    minosStyle->SetPadLeftMargin(0.15);

    //set axis label and title text sizes
    minosStyle->SetLabelSize(0.07,"xyz");
    minosStyle->SetTitleSize(0.08,"xyz");
    minosStyle->SetTitleOffset(0.9,"x");
    minosStyle->SetTitleOffset(0.8,"yz");
    minosStyle->SetStatFontSize(0.07);
    minosStyle->SetTextSize(0.08);
    minosStyle->SetTitleBorderSize(0);
    minosStyle->SetStatBorderSize(0);

    //set line widths
    minosStyle->SetHistLineWidth(2);
    minosStyle->SetFrameLineWidth(2);
    minosStyle->SetFuncWidth(2);

    //set the number of divisions to show
    minosStyle->SetNdivisions(506, "xy");

    //turn off xy grids
    minosStyle->SetPadGridX(0);
    minosStyle->SetPadGridY(0);

    //set the tick mark style
    minosStyle->SetPadTickX(1);
    minosStyle->SetPadTickY(1);

    //show the fit parameters in a box
    minosStyle->SetOptFit(1111);

    //turn off all other stats
    //minosStyle->SetOptStat(0000000);

    //marker settings
    minosStyle->SetMarkerStyle(8);
    minosStyle->SetMarkerSize(0.9);

    // Fonts
    const int kMinosFont = 42;

    minosStyle->SetStatFont(kMinosFont);
    minosStyle->SetLabelFont(kMinosFont,"xyz");
    minosStyle->SetTitleFont(kMinosFont,"xyz");
    minosStyle->SetTextFont(kMinosFont);

    //done
    minosStyle->cd();

    gROOT->ForceStyle();
    gStyle->ls();

    Int_t entries = angleHisto->GetEntries();
    Int_t cutentries = entries * 0.68;
    Int_t max = angleHisto->GetNbinsX();
    Double_t count = 0;
    Double_t total = 0;
    Int_t mybin;

    Double_t angle68 = angleHisto->GetBinCenter(25);

    TCanvas* distCanvas = new TCanvas("distCanvas","distance Histogram",640,480);
    distHisto->Draw("E");
    TCanvas* dist2Canvas = new TCanvas("dist2Canvas","distance Histogram",640,480);
    dist2Histo->Draw("E");

}
