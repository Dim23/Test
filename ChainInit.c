#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
using namespace std;

void MakeChain(const char* fileList, const char *outFile, const int line1, const int lineLast){
    TChain *chain = new TChain("mctree");
    ifstream ifile(fileList);
    char fileName[200];
    int nFiles = 0;
    while(ifile.getline(fileName,200)) {
        ++nFiles;
        if (nFiles > lineLast) { break; }
        if (nFiles < line1)  { continue; }
        cout << "Adding "<< fileName <<" to chain "<< outFile << endl;
        chain -> Add(fileName);
    }
    TFile *f = new TFile (outFile,"recreate");
    f -> cd();
    chain -> Write();
    f -> Close();
    cout << "Done. " << nFiles << " files are added to chain" << outFile <<endl;
}

void ChainInit(){
    char list[500];
    sprintf(list,"/home/dim2/FLOW5/file.list");
    
        
        MakeChain(list,"/home/dim2/FLOW5/chain/chainBIG.root",1,132);
    
    cout << "Tree chains were created successfully!" << endl;
}

