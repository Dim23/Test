
void Trees(const char *files,const char *outfile){
TChain *chain = new TChain("mсtree");
chain->Add(files);
TFile *f0=new TFile (outfile,"recreate");
f0->cd();
chain->Write();
f0->Close();
}
