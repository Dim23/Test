//double* test(){double *a=new double[2];a[0]=1;a[1]=3;return a;}

double start(double N0){TFile *f = TFile::Open("/home/dim2/FLOW5/All.root");
TH1I *hNpart, *hNcoll;hNpart = (TH1I *)f->Get("hNpart");double c;c=hNpart->Integral(N0,600)/hNpart->Integral();//hNpart->Draw();
return c;}

void start2(){TFile *fDim = TFile::Open("/home/dim2/FLOW5/PLOT/urqmd30_40_eta_m1p1.root");
                
TGraphErrors *DgrMC, *DgrV2,*DgrV4;
TGraphErrors *VgrMC, *VgrV2,*VgrV4;
DgrMC=(TGraphErrors *)fDim->Get("grMc");
DgrV2=(TGraphErrors *)fDim->Get("grV2");
DgrV4=(TGraphErrors *)fDim->Get("grV4");


TFile *fVin = TFile::Open("/home/dim2/FLOW5/TGraph.root");
VgrMC=(TGraphErrors *)fVin->Get("gr_0");
VgrV2=(TGraphErrors *)fVin->Get("gr_1");
VgrV4=(TGraphErrors *)fVin->Get("gr_2");
DgrMC->SetTitle("V_{2}{MC} Dim");
DgrV2->SetTitle("V_{2}{2} Dim");
DgrV4->SetTitle("V_{2}{4} Dim");
VgrMC->SetTitle("V_{2}{MC} Vin");
VgrV2->SetTitle("V_{2}{2} Vin");
VgrV4->SetTitle("V_{2}{4} Vin");


gROOT->SetStyle("Pub");
gStyle->SetOptStat(1111);
TCanvas *c1=new TCanvas("c1","Graph Draw Options1",650,500);
DgrMC->SetMarkerStyle(22);
DgrMC->SetMarkerSize(1);
DgrMC->SetMarkerColorAlpha(kGreen, 1);
DgrMC->SetLineColorAlpha(kGreen, 1);
DgrMC->SetLineWidth(1);
VgrMC->Draw("AP");DgrMC->Draw("P SAME");

TLegend *leg = new TLegend(0.2,.7,0.5,.8);
   leg->AddEntry(VgrMC,"V_{2}{MC} Vin","pe");
   leg->AddEntry(DgrMC,"V_{2}{MC} Dim","pe");
  leg -> SetFillColor(0);
  leg -> SetTextSize(0.04);
  leg -> SetTextFont(62);
  leg -> SetBorderSize(0);
  leg -> Draw();
c1->Print("/home/dim2/FLOW5/grvMC.png");

TCanvas *c2=new TCanvas("c2","Graph Draw Options2",650,500);
DgrV2->SetMarkerStyle(20);
DgrV2->SetMarkerSize(1);
DgrV2->SetMarkerColorAlpha(kBlue, 1);
DgrV2->SetLineColorAlpha(kBlue, 1);
DgrV2->GetYaxis()->SetTitle("V_{2}");
DgrV2->GetXaxis()->SetTitle("Pt,Gev/c");
DgrV2->SetLineWidth(2);
VgrV2->Draw("AP");DgrV2->Draw("P SAME");
   auto *leg1 = new TLegend(0.2,0.7,0.5,0.8);
   leg1->SetHeader("V_{2} from 2 particle Qum");
   leg1->AddEntry(VgrV2,"V_{2}{2} Vin","pe");
   leg1->AddEntry(DgrV2,"V_{2}{2} Dim","pe");
   leg1->Draw("P SAME");
c2->Print("/home/dim2/FLOW5/grv2.png");

TCanvas *c3=new TCanvas("c3","Graph Draw Options3",650,500);
DgrV4->SetMarkerStyle(21);
DgrV4->SetMarkerSize(1);
DgrV4->SetMarkerColorAlpha(kBlack, 1);
DgrV4->SetLineColorAlpha(kBlack, 1);
DgrV4->SetLineWidth(2);

   auto *leg2 = new TLegend(0.2,0.7,0.5,0.8);
   leg2->SetHeader("V_{4} fom r4 particle Qum");
   leg2->AddEntry(VgrV4,"V_{2}{4} Vin","pe");
   leg2->AddEntry(DgrV4,"V_{2}{4} Dim","pe");
VgrV4->Draw("AP");DgrV4->Draw("P SAME");
   leg2->Draw("P SAME");c3->Print("/home/dim2/FLOW5/grv4.png");

}
