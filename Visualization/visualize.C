#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPad.h"
#include "TLine.h"
#include "TGraph.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
 
void visualization(int eventNumber, std::string fileName, std::string tag = "", float maxZ = 10, float maxY = 10){

  TFile * f = TFile::Open(fileName.c_str(),"read");
  TTree * t = (TTree*)f->Get("vertextree_ana/vertexTree");

  std::vector< float > * zVtx_gen = 0;   
  unsigned int nVertices = 0;
  float zVtx[300] = {0};
  unsigned int nTracks[300] = {0};
  float zTrk[300][300] = {0};
  float etaTrk[300][300] = {0};
  int trackGenVertex[300][300] = {0};
  bool trackHPVtx[300][300] = {0};
  float contaminationFracByWeight[300] = {0};
  int bestMatchGenByWeight[300] = {0};

  t->SetBranchAddress("zVtx_gen",&zVtx_gen); 
  t->SetBranchAddress("nVertices",&nVertices);
  t->SetBranchAddress("zVtx",zVtx);
  t->SetBranchAddress("nTracks",nTracks);
  t->SetBranchAddress("trackZVtx",zTrk);
  t->SetBranchAddress("trackEtaVtx",etaTrk);
  t->SetBranchAddress("trackGenVertex",trackGenVertex);
  t->SetBranchAddress("trackHPVtx",trackHPVtx);
  t->SetBranchAddress("contaminationFracByWeight",&contaminationFracByWeight);
  t->SetBranchAddress("bestMatchGenByWeight",&bestMatchGenByWeight);

  t->GetEntry(eventNumber);
  
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas * c1 =new TCanvas("c1","c1",1200,800);
  c1->SetBorderSize(0);
  TPad * p1 = new TPad("p1","p1",0,0.2,1,1,0);
  TPad * p2 = new TPad("p2","p2",0,0,1,0.2,0);
  c1->SetLineWidth(0);
  p1->SetBottomMargin(0);
  p1->SetLeftMargin(0.1);
  p1->SetRightMargin(0.05);
  p1->SetTopMargin(0.07);
  p1->SetBorderSize(0);
  p1->Draw();
  p2->SetTopMargin(0);
  p2->SetLeftMargin(0.1);
  p2->SetRightMargin(0.05);
  p2->SetBottomMargin(0.42);
  p2->SetBorderSize(0);
  p2->Draw();
  p1->cd();

  TH1D * topDummy = new TH1D("topDummy",";z (cm); GenTop",1,-maxZ,maxZ);
  topDummy->Fill(0.0010);
  topDummy->GetYaxis()->SetRangeUser(-0.25,maxY);
  topDummy->GetYaxis()->SetTitle("Reco Radius (cm)");
  topDummy->GetYaxis()->CenterTitle();
  topDummy->Draw("p");
  
  p2->cd();
  TH1D * botDummy = new TH1D("botDummy",";z (cm); Gen",1,-maxZ,maxZ);
  botDummy->Fill(0.001);
  botDummy->GetYaxis()->SetRangeUser(-0.1,0.1);
  botDummy->GetYaxis()->CenterTitle();
  botDummy->GetXaxis()->CenterTitle();
  botDummy->GetYaxis()->SetTitleSize(0.13);
  botDummy->GetYaxis()->SetTitleOffset(0.06);
  botDummy->GetYaxis()->SetLabelSize(0);
  botDummy->GetXaxis()->SetTitleSize(0.2);
  botDummy->GetXaxis()->SetLabelSize(0.2);
  botDummy->Draw("p");
  

  TMarker * m[300];
  TLine * l[300];
  std::sort(zVtx_gen->begin(), zVtx_gen->end());
  std::vector< int > colorID;
  int c = 0;
  for(unsigned int i = 0; i<zVtx_gen->size(); i++){
    m[i] = new TMarker(zVtx_gen->at(i), 0, 8);
    m[i]->SetMarkerColor(c%2+1);
    m[i]->SetMarkerSize(1.4);
    m[i]->Draw("same");
    l[i] = new TLine(zVtx_gen->at(i),-0.1,zVtx_gen->at(i),0.1);
    l[i]->SetLineColor(c%2 +1);
    l[i]->SetLineStyle(2);
    l[i]->Draw();
    colorID.push_back(c);
    c++;
  }

  p1->cd();
  TLine * lR[300];
  std::vector< float > z;
  std::vector< float > matchedIndexes;
  for(unsigned int i = 0; i<nVertices; i++){
    z.push_back(zVtx[i]);
    matchedIndexes.push_back(bestMatchGenByWeight[i]);
  }

  int nSplit = 0;
  for(unsigned int i = 0; i<zVtx_gen->size(); i++){
    int n = std::count(matchedIndexes.begin(), matchedIndexes.end(),i);
    if(n>1) nSplit++;
  }
  std::sort( z.begin(),z.end());

  TLine * l1[10000];
  int lIndex = 0;
  for(unsigned int i = 0; i<z.size(); i++){
    int index = -1;
    for(unsigned int j = 0; j<z.size(); j++){
      if( TMath::Abs( zVtx[j] - z.at(i) ) < 0.001) index = j;
    }
    for(unsigned int j = 0; j<nTracks[index]; j++){
      //if(j==0) std::cout << zTrk[index][j] << " " << etaTrk[index][j] << std::endl;
      float theta = 2*TMath::ATan(TMath::Exp(-etaTrk[index][j]));
      float xStop = 0;
      float yStop = 3;
      if( (theta<TMath::Pi()/2.0) && (TMath::Tan(theta) < maxY/(maxZ-zTrk[index][j])) ){
        xStop = maxZ;
        yStop = TMath::Tan(theta) * (maxZ-zTrk[index][j]);
      }
      else if( theta<TMath::Pi()/2.0 ){
        xStop = zTrk[index][j] + 1.0/TMath::Tan(theta)*(maxY);
        yStop = maxY;
      } 
      else if( (theta>TMath::Pi()/2.0) && (TMath::Tan(TMath::Pi()-theta) < maxY/(TMath::Abs(zTrk[index][j] + maxZ))) ){
        xStop = -maxZ;
        yStop = TMath::Tan(TMath::Pi()-theta) * (TMath::Abs(zTrk[index][j] + maxZ));
      }
      else if( theta>TMath::Pi()/2.0 ){
        xStop = zTrk[index][j] - 1.0/TMath::Tan(TMath::Pi()-theta)*(maxY) ;
        yStop = maxY;
      }
      else{
        std::cout << "uh oh" << std::endl;
      }
      
      l1[lIndex] = new TLine(zTrk[index][j] ,0,xStop ,yStop);
      l1[lIndex]->SetLineColor(i%2+1);
      if(trackGenVertex[index][j] != bestMatchGenByWeight[index] ){
        l1[lIndex]->SetLineColor(4);
        l1[lIndex]->SetLineWidth(2);
        l1[lIndex]->SetLineStyle(2);
      }
      //if(trackHPVtx[index][j] == false) l1[lIndex]->SetLineStyle(2);
      l1[lIndex]->Draw("same");
      lIndex++;      
    }
  }
  TMarker * mR[300];
  TMarker * mR_noContam[300];
  TMarker * mR_Split[300];
  int uncontam = 0;
  int nRecoSplit = 0;
  int nUncontamNoSplit = 0;
  for(unsigned int i = 0; i<z.size(); i++){
    bool isUncontam = false;
    bool isSplit = false;
    int index = -1;
    for(unsigned int j = 0; j<z.size(); j++){
      if( TMath::Abs( zVtx[j] - z.at(i) ) < 0.001) index = j;
    }
    mR[i] = new TMarker(zVtx[index], 0, 8);
    mR[i]->SetMarkerColor(i%2+1);
    mR[i]->SetMarkerSize(1.5);
    mR[i]->Draw("same");
    if(contaminationFracByWeight[index] < 0.0001){
      mR_noContam[i] = new TMarker(zVtx[index], -0.1, 8);
      mR_noContam[i]->SetMarkerColor(kGreen+1);
      mR_noContam[i]->SetMarkerSize(0.9);
      mR_noContam[i]->Draw("same");   
      uncontam++; 
      isUncontam = true;
    }
    
    int splitNum = std::count(matchedIndexes.begin(), matchedIndexes.end(), bestMatchGenByWeight[index]);
    if(splitNum>1){
      mR_Split[i] = new TMarker(zVtx[index], -0.17, 8);
      mR_Split[i]->SetMarkerColor(kMagenta);
      mR_Split[i]->SetMarkerSize(0.9);
      mR_Split[i]->Draw("same");   
      nRecoSplit++;
      isSplit = true;
    }
    if(!isSplit && isUncontam) nUncontamNoSplit++;
  } 

  TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry((TObject*)0,Form("%d Gen Vtxs",(int)zVtx_gen->size()),"");
  leg->AddEntry((TObject*)0,Form("%d Split Gen Vtxs",(int)nSplit),"");
  leg->AddEntry((TObject*)0,Form("%d Reco Vtxs",(int)nVertices),"");
  leg->AddEntry((TObject*)0,Form("%d Uncontaminated Vtxs",(int)uncontam),"");
  leg->AddEntry((TObject*)0,Form("%d Split Reco Vtxs",(int)nRecoSplit),"");
  leg->AddEntry((TObject*)0,Form("%d Uncontam+Unsplit Reco Vtxs",(int)nUncontamNoSplit),"");
  leg->Draw("same");

  c1->SaveAs(Form("visualization_Evt%d_Tag%s.png",eventNumber, tag.c_str()));
  c1->SaveAs(Form("visualization_Evt%d_Tag%s.pdf",eventNumber, tag.c_str()));
  c1->SaveAs(Form("visualization_Evt%d_Tag%s.C",eventNumber, tag.c_str()));

}

void loopVisualization(int nEvents, std::string fileName, std::string tag = "", float maxZ = 10, float maxY = 10){
 for(int i = 0; i<nEvents; i++){
   visualization(i, fileName, tag+"_"+std::to_string(i), maxZ, maxY);
 }
}
