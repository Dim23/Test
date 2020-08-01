//double* test(){double *a=new double[2];a[0]=1;a[1]=3;return a;}
#include <iostream>
#include <fstream>
double start(double N0){TFile *f = TFile::Open("/home/dim2/FLOW5/All.root");
TH1I *hNpart, *hNcoll;hNpart = (TH1I *)f->Get("hNpart");double c;c=hNpart->Integral(N0,600)/hNpart->Integral();//hNpart->Draw();
return c;}

void start2(){TFile *f = TFile::Open("/home/dim2/FLOW5/All.root");
TH1I *hNpart;
TH1D *hBimp,*hETA,*hdPHI,*hPt;


gROOT->SetStyle("Pub");
gStyle->SetOptStat(1111);
TCanvas *c1=new TCanvas("c1","Graph Draw Options1",650,500);
hNpart = (TH1I *)f->Get("hNpart");hNpart->GetXaxis()->SetTitle("N_{part}");hNpart->Draw();

TCanvas *c2=new TCanvas("c2","Graph Draw Options2",650,500);
hBimp = (TH1D *)f->Get("hBimp");hBimp->GetXaxis()->SetTitle("b,fm");hBimp->Draw();

TCanvas *c3=new TCanvas("c3","Graph Draw Options3",650,500);
hETA = (TH1D *)f->Get("hETA");hETA->GetXaxis()->SetTitle("#eta");hETA->Draw();

TCanvas *c4=new TCanvas("c4","Graph Draw Options4",650,500);
hdPHI = (TH1D *)f->Get("hdPHI");hdPHI->GetXaxis()->SetTitle("#varphi-#psi_{R},rad");hdPHI->Draw();

TCanvas *c5=new TCanvas("c5","Graph Draw Options5",650,500);
hPt = (TH1D *)f->Get("hPt");hPt->GetXaxis()->SetTitle("Pt,GeV/c");hPt->Draw();
}

 
TCanvas *hlabels2()
{gROOT->SetStyle("Pub");
gStyle->SetOptStat(1111);
TFile *fDim = TFile::Open("/home/dim2/FLOW5/PLOT/urqmd30_40_eta_m1p1.root");
                
TGraphErrors *DgrMC, *DgrV2,*DgrV4;
DgrMC=(TGraphErrors *)fDim->Get("grIntMC");
DgrV2=(TGraphErrors *)fDim->Get("grIntV2");
DgrV4=(TGraphErrors *)fDim->Get("grIntV4");

   const Int_t nx = 3;
   const char *month[nx]  = {"V_{2}{MC}","V_{2}{2}","V_{2}{4}"};

   TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,600,600);
   c1->SetLeftMargin(0.15);c1->SetRightMargin(0.02);
   c1->SetBottomMargin(0.05);c1->SetTopMargin(0.02);
   TH1F *h = new TH1F("h","Reference flow V_{2}",3,0.,3.);
   
   h->SetLineColorAlpha(kBlack, 1);
   h->SetCanExtend(TH1::kAllAxes);
   h->SetStats(0);
   for (Int_t i=0;i<3;i++) {
      h->Fill(month[i],1);
   }
   h->LabelsDeflate("X");
   
   double *vMC=DgrMC->GetY(); 
   double *v4=DgrV4->GetY();
   double *ev4=DgrV4->GetEY();
   double *v2=DgrV2->GetY();
   double *ev2=DgrV2->GetEY();
   TLine *line=new TLine(0.,vMC[0],3.,vMC[0]);
   line->SetLineWidth(1);
   line->SetLineStyle(5);
double ymin=min(v4[0]-ev4[0],v2[0]-ev4[0]);
double ymax=max(v4[0]+ev4[0],v2[0]+ev4[0]);
   h->GetYaxis()->SetRangeUser(0.9*ymin,1.1*ymax);
   h->Draw();
   DgrV2->Draw("P SAME");DgrV4->Draw("P SAME");
   DgrMC->Draw("P SAME");line->Draw("SAME");
   return c1;
}
