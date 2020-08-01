#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>

double sigmaX(TH1F *X,double weight){
double RMS=X->GetStdDev();
double neff=X->GetEffectiveEntries();
return weight*RMS*RMS/(neff-1);
}

double sigmaXY(TH1F *X,TH1F *Y,TH1F *XY,double weight){
double RMS=XY->GetMean()-(X->GetMean())*(Y->GetMean());
double Xsum=X->GetSum();
double Ysum=Y->GetSum();
double XYsum=XY->GetSum();
double neff=Xsum*Ysum/XYsum;
return weight*RMS*RMS/(neff-1);
}
void read(const char *outfile,const char *savefile="~/FLOW5/PLOT/urqmd30_40_eta_m1p1.root"){
//double pt_bin[11]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};int NN=10;
//double pt_bin[15]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2};int NN=14;//бины пот импульсу для дифференциального потока
double b_bin[9]={0.0,4.18,6.01,7.37,8.52,9.57,10.55,11.46,12.31};static const int Nb=8;
double pt_bin[13]={0.,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.7,2.0}; static const int NN=12;
//double pt_bin[11]={0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0}; static const int NN=10;
double RMSSV2[Nb],RMSSV4[Nb],RMSVMC[Nb];//RMS для SV2 и SV4

double RMSdn2[Nb][NN];double RMSv2[Nb][NN];
double RMSdn4[Nb][NN];double RMSv4[Nb][NN];
double Rv2[Nb][NN];double Rv4[Nb][NN];
double mc[Nb][NN];double RMSmc[Nb][NN];
double RMSbinPt[Nb][NN];

double dsv4[Nb][NN],dsv2[Nb][NN],sv4[Nb],sv2[Nb],v4[Nb],v2[Nb],vMC[Nb],cn2[Nb],cn4[Nb];


TH1F *SV2[Nb];//2 частичная референсная корреляция
TH1F *SV4[Nb];//4 частичная референсная корреляция
TH1F *COSf1[Nb];
TH1F *SINf1[Nb];
TH1F *COSf1f2[Nb];
TH1F *SINf1f2[Nb];
TH1F *COSf1f2f3[Nb];
TH1F *SINf1f2f3[Nb];

TH1F *SV2_SV4[Nb];
TH1F *HMC[Nb];
TH1F *DiffSV2[Nb][NN];TH1F *DiffSV4[Nb][NN];TH1F *DiffMC[Nb][NN];
TH1F *SV2_DiffSV2[Nb][NN];TH1F *SV2_DiffSV4[Nb][NN];
TH1F *SV4_DiffSV2[Nb][NN];TH1F *SV4_DiffSV4[Nb][NN];
TH1F *DiffSV2_DiffSV4[Nb][NN];
TH1F *hBin_Pt[Nb][NN];
TH1F *COSp1[Nb][NN];TH1F *SINp1[Nb][NN];
TH1F *COSp1f2[Nb][NN];TH1F *SINp1f2[Nb][NN];
TH1F *COSp1f2mf3[Nb][NN];TH1F *SINp1f2mf3[Nb][NN];
TH1F *COSp1mf2mf3[Nb][NN];TH1F *SINp1mf2mf3[Nb][NN];


double binPt[Nb][NN];
double QMC[Nb][NN];
double dn2[Nb][NN];double dn4[Nb][NN];
double Diffv2[Nb][NN]; double Diffv4[Nb][NN];
double Diffsv2[Nb][NN]; double Diffsv4[Nb][NN];

char strCOVsv2sv2[200];char strCOVsv2sv4[200];
char strCOVsv4sv2[200];char strCOVsv4sv4[200];
char DiffCOVsv2sv4[200];
char strE[200];char strW[200];char strMC[200];

char strCOSp1[200];char strSINp1[200];
char strCOSp1f2[200];char strSINp1f2[200];
char strCOSp1f2mf3[200];char strSINp1f2mf3[200];
char strCOSp1mf2mf3[200];char strSINp1mf2mf3[200];char hBinPt[200];
char str[200];
const char *STR;
cout <<" v2{2} is "<<endl;

TFile *f=new TFile(outfile);
for(int k=0;k<Nb;k++){


sprintf (str, "%sCENT%d_%d", "sv2", k*10,(k+1)*10);
STR=(char*)str;
SV2[k]=(TH1F*)f->Get(STR);

sprintf (str, "%sCENT%d_%d", "sv4", k*10,(k+1)*10);
STR=(char*)str;
SV4[k]=(TH1F*)f->Get(STR);

sprintf (str, "%sCENT%d_%d", "sv2_sv4", k*10,(k+1)*10);
STR=(char*)str;
SV2_SV4[k]=(TH1F*)f->Get(STR);

sprintf (str, "%sCENT%d_%d", "mc", k*10,(k+1)*10);
STR=(char*)str;
HMC[k]=(TH1F*)f->Get(STR);

for(int m=0;m<NN;m++){
sprintf(strE,"%s%d_CENT%d_%d","sv2Diff",m,k*10,(k+1)*10);
sprintf(strW,"%s%d_CENT%d_%d","sv4Diff",m,k*10,(k+1)*10);
sprintf(strMC,"%s%d_CENT%d_%d","DiffMC",m,k*10,(k+1)*10);
sprintf(strCOVsv2sv2,"%s%d_CENT%d_%d","SV2_DiffSV2",m,k*10,(k+1)*10);
sprintf(strCOVsv2sv4,"%s%d_CENT%d_%d","SV2_DiffSV4",m,k*10,(k+1)*10);
sprintf(strCOVsv4sv2,"%s%d_CENT%d_%d","SV4_DiffSV2",m,k*10,(k+1)*10);
sprintf(strCOVsv4sv4,"%s%d_CENT%d_%d","SV4_DiffSV4",m,k*10,(k+1)*10);
sprintf(DiffCOVsv2sv4,"%s%d_CENT%d_%d","DiffSV2_DiffSV4",m,k*10,(k+1)*10);

sprintf(strCOSp1,"%s%d_CENT%d_%d","COSp1",m,k*10,(k+1)*10);
sprintf(strCOSp1f2,"%s%d_CENT%d_%d","COSp1f2",m,k*10,(k+1)*10);
sprintf(strCOSp1f2mf3,"%s%d_CENT%d_%d","COSp1f2mf3",m,k*10,(k+1)*10);
sprintf(strCOSp1mf2mf3,"%s%d_CENT%d_%d","COSp1mf2mf3",m,k*10,(k+1)*10);
sprintf(strSINp1,"%s%d_CENT%d_%d","SINp1",m,k*10,(k+1)*10);
sprintf(strSINp1f2,"%s%d_CENT%d_%d","SINp1f2",m,k*10,(k+1)*10);
sprintf(strSINp1f2mf3,"%s%d_CENT%d_%d","SINp1f2mf3",m,k*10,(k+1)*10);
sprintf(strSINp1mf2mf3,"%s%d_CENT%d_%d","SINp1mf2mf3",m,k*10,(k+1)*10);
sprintf(hBinPt,"%s%d_CENT%d_%d","BIN_pt_",m,k*10,(k+1)*10);
const char *EastH=(char*)strE;
const char *WastH=(char*)strW;
const char *SMC=(char*)strMC;
const char *stCOVsv2sv2=(char*)strCOVsv2sv2;
const char *stCOVsv2sv4=(char*)strCOVsv2sv4;
const char *stCOVsv4sv2=(char*)strCOVsv4sv2;
const char *stCOVsv4sv4=(char*)strCOVsv4sv4;
const char *DifCOVsv2sv4=(char*)DiffCOVsv2sv4;
const char *HBinPt=(char*)hBinPt;
DiffSV2[k][m]=(TH1F*)f->Get(EastH);
DiffSV4[k][m]=(TH1F*)f->Get(WastH);
DiffMC[k][m]=(TH1F*)f->Get(SMC);
SV2_DiffSV2[k][m]=(TH1F*)f->Get(stCOVsv2sv2);
SV2_DiffSV4[k][m]=(TH1F*)f->Get(stCOVsv2sv4);
SV4_DiffSV2[k][m]=(TH1F*)f->Get(stCOVsv4sv2);
SV4_DiffSV4[k][m]=(TH1F*)f->Get(stCOVsv4sv4);
DiffSV2_DiffSV4[k][m]=(TH1F*)f->Get(DifCOVsv2sv4);
hBin_Pt[k][m]=(TH1F*)f->Get(HBinPt);
RMSbinPt[k][m]=0;
}}


//Референсный поток
for(int k=0;k<Nb;k++){
sv2[k]=SV2[k]->GetMean();v2[k]=pow(fabs(sv2[k]),0.5);
sv4[k]=SV4[k]->GetMean();cn4[k]=sv4[k]-2*(sv2[k]*sv2[k]);v4[k]=pow(fabs(cn4[k]),0.25);
vMC[k]=HMC[k]->GetMean();
cout <<" v2{2} is "<<v2[k]<<" v2{4} is "<<v4[k]<<" v2{MC} is "<<vMC[k] <<endl; 

//Дифференциальный поток
cout <<"Different. flow for cent "<<10*(k+0.5)<<endl;
for(int m=0;m<NN;m++){
//DiffSV2[m]->Sumw2();DiffSV4[m]->Sumw2();
dn2[k][m]=DiffSV2[k][m]->GetMean();Diffv2[k][m]=dn2[k][m]/pow(fabs(sv2[k]),0.5);
dsv4[k][m]=DiffSV4[k][m]->GetMean();dn4[k][m]=dsv4[k][m]-2*dn2[k][m]*sv2[k];Diffv4[k][m]=-dn4[k][m]/pow(fabs(cn4[k]),0.75);
binPt[k][m]=hBin_Pt[k][m]->GetMean();

//cout <<"Different. flow "<<" v2{2} "<<Diffv2[k][m]<"+- "<<Diffv2[k][m]<<", v2{4} "<<Diffv4[k][m]<<" Pt "<<binPt[m]<<endl ;}
}
//погрешность для дифференциального потока

for(int m=0;m<NN;m++){
mc[k][m]=DiffMC[k][m]->GetMean();
RMSmc[k][m]=pow(sigmaX(DiffMC[k][m],1),0.5);RMSdn2[k][m]=sigmaX(DiffSV2[k][m],1);RMSdn4[k][m]=DiffSV4[k][m]->GetStdDev();

Rv2[k][m]=0.25*pow(sv2[k],-3)*(pow(dn2[k][m],2)*sigmaX(SV2[k],1)+4*sv2[k]*sv2[k]*sigmaX(DiffSV2[k][m],1)-4*dn2[k][m]*sv2[k]*sigmaXY(SV2[k],DiffSV2[k][m],SV2_DiffSV2[k][m],1));
RMSv2[k][m]=pow(Rv2[k][m],0.5);

Rv4[k][m]=(pow(fabs(cn4[k]),-3.5))*(pow((2*sv2[k]*sv2[k]*dn2[k][m]-3*sv2[k]*dsv4[k][m]+2*sv4[k]*dn2[k][m]),2)*sigmaX(SV2[k],1)+

9/16*pow(dn4[k][m],2)*sigmaX(SV4[k],1)+ 4*sv2[k]*sv2[k]*pow(fabs(cn4[k]),2)*sigmaX(DiffSV2[k][m],1)+

pow(fabs(cn4[k]),2)*sigmaX(DiffSV4[k][m],1)-

1.5*fabs(dn4[k][m])*(2*sv2[k]*sv2[k]*dn2[k][m]-3*sv2[k]*dsv4[k][m]+2*sv4[k]*dn2[k][m])*sigmaXY(SV2[k],SV4[k],SV2_SV4[k],1)-

4*sv2[k]*fabs(cn4[k])*(2*sv2[k]*sv2[k]*dn2[k][m]-3*sv2[k]*dsv4[k][m]+2*sv4[k]*dn2[k][m])*sigmaXY(SV2[k],DiffSV2[k][m],SV2_DiffSV2[k][m],1)+

2*fabs(cn4[k])*(2*sv2[k]*sv2[k]*dn2[k][m]-3*sv2[k]*dsv4[k][m]+2*sv4[k]*dn2[k][m])*sigmaXY(SV2[k],DiffSV4[k][m],SV2_DiffSV4[k][m],1)+

3*sv2[k]*(fabs(cn4[k]))*dn4[k][m]*sigmaXY(SV4[k],DiffSV2[k][m],SV4_DiffSV2[k][m],1)-

1.5*fabs(cn4[k])*sigmaXY(SV4[k],DiffSV4[k][m],SV4_DiffSV4[k][m],1)-

4*sv2[k]*pow(fabs(cn4[k]),2)*sigmaXY(DiffSV2[k][m],DiffSV4[k][m],DiffSV2_DiffSV4[k][m],1));
/*
cout <<"8 "<< 4*sv2[k]*pow(fabs(cn4[k]),2)*sigmaXY(DiffSV2[k][m],DiffSV4[k][m],DiffSV2_DiffSV4[k][m],1)<<endl;
cout <<"7 "<< 1.5*fabs(cn4[k])*sigmaXY(SV4[k],DiffSV4[k][m],SV4_DiffSV4[k][m],1)<<endl;
cout <<"6 "<< sigmaXY(SV4[k],DiffSV2[k][m],SV4_DiffSV2[k][m],1)<<endl;
cout <<"5 "<< sigmaXY(SV2[k],DiffSV4[k][m],SV2_DiffSV4[k][m],1) <<endl;
cout <<"4 "<< sigmaXY(SV2[k],DiffSV2[k][m],SV2_DiffSV2[k][m],1) <<endl;
cout <<"3 "<< sigmaXY(SV2[k],SV4[k],SV2_SV4[k],1) <<endl;
cout <<"2 "<< pow(fabs(cn4[k]),2)*sigmaX(DiffSV4[k][m],1) <<endl;
cout <<"1 "<< 9/16*pow(dn4[k][m],2)*sigmaX(SV4[k],1)+ 4*sv2[k]*sv2[k]*pow(fabs(cn4[k]),2)*sigmaX(DiffSV2[k][m],1) <<endl;
cout<<"0 "<< pow((2*sv2[k]*sv2[k]*dn2[k][m]-3*sv2[k]*dsv4[k][m]+2*sv4[k]*dn2[k][m]),2)*sigmaX(SV2[k],1)<< endl;*/
RMSv4[k][m]=pow(Rv4[k][m],0.5);
cout <<"Different. flow "<<" v2{2} "<<Diffv2[k][m]<<"+- "<<RMSv2[k][m]<<", v2{4} "<<Diffv4[k][m]<<"+- "<<RMSv4[k][m]<<" Pt "<<binPt[k][m]<<endl ;
}}

TCanvas *c1=new TCanvas("c1","Graph Draw Options1",650,500);

//TF1 *momentdist = new TF1("momentdist","[1]*x+[0]", 0., 2.1);
//TF1 *momentdist = new TF1("momentdist","[1]*(exp(2*x)-1)/(exp(2*x)+1)+[0]", 0., 2.3);
/*
momentdist->FixParameter(0, 0.02);
momentdist->FixParameter(1, 0.1);
momentdist->SetLineColorAlpha(kBlack, 0.8);
momentdist->SetLineWidth(2);momentdist->SetLineStyle(5);
momentdist->GetXaxis()->SetTitle("Pt,Gev/c");
momentdist->GetYaxis()->SetTitle("V_{2}");
momentdist->SetTitle("V_{2}(Pt)");
momentdist->Draw();*/

const int kk=3;
TGraphErrors *grv4 = new TGraphErrors(NN,binPt[kk],Diffv4[kk],RMSbinPt[kk],RMSv4[kk]);
grv4->SetName("diff_v4");
grv4->GetYaxis()->SetRangeUser(0,0.15);
grv4->SetMarkerStyle(21);
grv4->SetMarkerSize(1);
grv4->SetMarkerColorAlpha(kBlack, 1);
grv4->SetLineColorAlpha(kBlack, 1);
grv4->SetLineWidth(2);
grv4->GetXaxis()->SetTitle("Pt,Gev/c");
grv4->GetYaxis()->SetTitle("V_{2}");grv4->Write();
grv4->Draw("AP");
grv4->SetTitle("V_{2}{4}");

TGraphErrors *grv2 = new TGraphErrors(NN,binPt[kk],Diffv2[kk],RMSbinPt[kk],RMSv2[kk]);
grv2->SetName("diff_v2");
grv2->SetMarkerStyle(20);
grv2->SetMarkerSize(1);
grv2->SetMarkerColorAlpha(kBlue, 1);
grv2->SetLineColorAlpha(kBlue, 1);
grv2->GetYaxis()->SetTitle("V_{2}");
grv2->GetXaxis()->SetTitle("Pt,Gev/c");
grv2->SetLineWidth(2);
grv2->Draw("SAME P");
grv2->SetTitle("V_{2}{2}");

TGraphErrors *grv = new TGraphErrors(NN,binPt[kk],mc[kk],RMSbinPt[kk],RMSmc[kk]);
grv->SetName("diff_vMC");
grv->SetMarkerStyle(22);
grv->SetMarkerSize(1);
grv->SetMarkerColorAlpha(kGreen, 1);
grv->SetLineColorAlpha(kGreen, 1);
grv->SetLineWidth(1);
grv->Draw("SAME P");
grv->SetTitle("V_{2}{MC}");


gPad->BuildLegend();
gROOT->SetStyle("Pub");
gStyle->SetOptStat(1111);

//Рефренсный поток с погрешностями
TCanvas *c2=new TCanvas("c2","Graph Draw Options",650,500);
double x[1]={1},y[1]={vMC[kk]},ex[1]={0.},ey[1]={pow(sigmaX(HMC[kk],1),0.5)};
TGraphErrors *grMC = new TGraphErrors(1,x,y,ex,ey);
grMC->SetName("vMC");
grMC->SetMarkerStyle(20);
grMC->GetYaxis()->SetTitle("V_{2}");

double x2[1]={1.4},y2[1]={v2[kk]},ex2[1]={0.},ey2[1]={pow(sigmaX(SV2[kk],0.25/sv2[kk]),0.5)};
TGraphErrors *grsv2 = new TGraphErrors(1,x2,y2,ex2,ey2);
grsv2->SetName("v2");
grsv2->SetMarkerStyle(20);

double x4[1]={1.8},y4[1]={v4[kk]},ex4[1]={0.},ey4[1]={0};
double s4=pow(fabs(cn4[kk]),-1.5)*(sigmaX(SV2[kk],sv2[kk]*sv2[kk])+sigmaX(SV4[kk],0.125)-0.5*sv2[kk]*sigmaXY(SV2[kk],SV4[kk],SV2_SV4[kk],1));
ey4[0]=pow(s4,0.5);
TGraphErrors *grsv4 = new TGraphErrors(1,x4,y4,ex4,ey4);
grsv4->SetName("v4");
grsv4->SetMarkerStyle(20);
cout <<"MC eror "<<ey[0]<<" SV2 eror "<<ey2[0]<<" SV4 eror "<<ey4[0]<<endl;
grMC->GetXaxis()->SetRangeUser(0,1.9);
grMC->GetYaxis()->SetRangeUser(y4[0]-0.003,y4[0]+0.003);


TLine *line=new TLine(0.9,vMC[kk],1.9,vMC[kk]);
line->SetLineWidth(1);
line->SetLineStyle(5);
grMC->Draw("AP");
grsv4->Draw("SAME P");
grsv2->Draw("SAME P");
line->Draw("SAME");


TFile *d_outfile = new TFile(savefile,"recreate");
d_outfile->cd();

grv4->Write("grV4");
grv2->Write("grV2");
grv->Write("grMc");
grMC->Write("grIntMC");
grsv4->Write("grIntV2");
grsv2->Write("grIntV4");
d_outfile->Close();
}
