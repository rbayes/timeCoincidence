#include "HistContainer.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

bool sortByTime(const std::vector<int> &lhs, const std::vector<int> &rhs) 
{return (lhs[5] < rhs[5]) || ((lhs[5] == rhs[5]) && (lhs[4] < lhs[4]))
    || ((lhs[5] == rhs[5]) && (lhs[4] == lhs[4]) && (lhs[3] < rhs[3]));} 

void cycleTree(std::vector<std::string> infile, std::string outfile, TVector3 spos, int forceCal){
  
  TChain* tree = new TChain("output");
  
  for(int itr=0; itr<infile.size(); itr++){
    tree->Add(infile.at(itr).c_str());
    std::cout<<"Added "<<infile.at(itr).c_str()<<" to chain with "
	     <<tree->GetEntries()<<" entries"<<std::endl;
  }

  double energy, beta14, thetaij, itr,
    posx, posy, posz, posr,
    dirx, diry, dirz;
  double biPolike212, biPolike214, biPoCumul, alphaBeta212, alphaBeta214;
  double meanTime=0;
  int nhit=0, unstime=0, ustime=0, udtime=0;
  double timingPeaks=0.0;
  ULong64_t dcApplied, dcFlagged;
  bool fitValid=false, isCal=false;
  double utime=0; // in microseconds
  tree->SetBranchAddress("time",&meanTime);
  tree->SetBranchAddress("isCal", &isCal);
  tree->SetBranchAddress("dcApplied", &dcApplied);
  tree->SetBranchAddress("dcFlagged", &dcFlagged);
  tree->SetBranchAddress("timingPeaks",&timingPeaks);
  tree->SetBranchAddress("fitValid",&fitValid);
  tree->SetBranchAddress("nhits", &nhit);
  tree->SetBranchAddress("uTNSecs", &unstime);
  tree->SetBranchAddress("uTSecs", &ustime);
  tree->SetBranchAddress("uTDays", &udtime);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchAddress("beta14", &beta14);
  tree->SetBranchAddress("thetaij", &thetaij);
  tree->SetBranchAddress("itr", &itr);
  tree->SetBranchAddress("posx", &posx);
  tree->SetBranchAddress("posy", &posy);
  tree->SetBranchAddress("posz", &posz);
  tree->SetBranchAddress("posr", &posr);
  tree->SetBranchAddress("dirx", &dirx);
  tree->SetBranchAddress("diry", &diry);
  tree->SetBranchAddress("dirz", &dirz);
  tree->SetBranchAddress("biPoLikelihood212", &biPolike212);
  tree->SetBranchAddress("biPoLikelihood214", &biPolike214);
  tree->SetBranchAddress("biPoCumul",&biPoCumul);
  tree->SetBranchAddress("alphaBeta212",&alphaBeta212);
  tree->SetBranchAddress("alphaBeta214",&alphaBeta214);
  HistContainer* histos = new HistContainer(outfile, spos);
  
  int entries = tree->GetEntries();
  int i=0;
  /*
  double Iutime, Ienergy, Ibeta14, Ithetaij, 
    Iposx, Iposy, Iposz, Iposr;
  double Futime, Fenergy, Fbeta14, Fthetaij, 
    Fposx, Fposy, Fposz, Fposr;
  */
  int mintime=99999999;
  int Pnhit=0, Ptime=0, Pcal=0, Ptmean=0, Pvalid=0,
    Pustime=0, Punstime=0, Pudtime=0;
  std::cout<<"Tree has "<<entries<<" entries"<<std::endl;
  std::vector<std::vector<int> > indices0;
  std::vector<std::vector<int> > indices1;
  std::vector<int> iprompt;
  std::vector<int> idelay;
  tree->GetEntry(0);
  double firsttime = double(ustime);
  double firstday  = double(udtime);
  for (int j=0; j < entries; j++){
    tree->GetEntry(j);
    double ltime = ((double(udtime) - firstday)*86400 + double(ustime) + double(unstime)*1e-9)/1e-6;
    if ( ( ( dcApplied & 0xFB0000017FFE) & dcFlagged) == (dcApplied & 0xFB0000017FFE)){
      if(isCal || forceCal==1){
	if(j%1000==0) std::cout<<j<<std::endl;
	histos->SetPromptData(fitValid, meanTime,
			      ltime, nhit, energy, beta14, 
			      thetaij, itr, TVector3(posx, posy, posz),
			      TVector3(dirx, diry, dirz));
	histos->FillAntiHistograms();
      }
      else {
	if(j%1000==0) std::cout<<j<<std::endl;
	std::vector<int> tempt;
	utime = ((double(udtime) - firstday)*86400 + (double(ustime) + double(unstime)*1e-9) - firsttime)/ 1e-4;
	tempt.push_back(j);
	tempt.push_back(utime);
	tempt.push_back(nhit);
	tempt.push_back(ustime);
	tempt.push_back(unstime);
	tempt.push_back(udtime);
	tempt.push_back(isCal || forceCal==1 ? 1:0);
	tempt.push_back(fitValid);
	tempt.push_back(int(meanTime));
	tempt.push_back(posx);
	tempt.push_back(posy);
	tempt.push_back(posz);
	indices0.push_back(tempt);
	// indices1.push_back(tempt);

      }
    }
  }
  if (indices0.size() > 0){
    std::sort(indices0.begin(), indices0.end(), sortByTime);
    // std::sort(indices1.begin(), indices1.end(), sortByTime);
    // for (int l=0; l<entries; l++){
    indices1 = indices0;
    auto l = indices0.begin();
    int ltest=0;
    do {
      
      Ptime = (*l)[1];
      Pnhit = (*l)[2];
      Pustime = (*l)[3];
      Punstime = (*l)[4];
      Pudtime = (*l)[5];
      Pcal  = (*l)[6];
      Pvalid = (*l)[7];
      Ptmean = (*l)[8];
      TVector3 promptPos(double((*l)[9]) - posx,
			 double((*l)[10]) - posy,
			 double((*l)[11]) - posz);
      // Pposx  = (*l)[9];
      // Pposy  = (*l)[10];
      // Pposz  = (*l)[11];

      
      int minj = -1;
      int minnhit = 9999;
      int itrmin = -1;
      int k=0;
      if((Pnhit > 12 && indices1.size() > 0) && Pvalid==1) {
	auto pev = find(indices1.begin(), indices1.end(), (*l));
	if(!Pcal){ // If this is not a calibration event run time coincidence
	  for(auto j=indices1.begin(); j < pev+3; j++){
	    
	    // std::cout<<indices0[l][0]<<std::endl;
	    int loctime = (*j)[1];
	    int locnhit = (*j)[2];
	    int locdtime = (*j)[5];
	    int locvalid = (*j)[7];
	    int loctmean = (*j)[8];
	    TVector3 locpos(double((*j)[9])-posx,
			    double((*j)[10])-posy,
			    double((*j)[11])-posz);
	    double tdiff = (loctime - Ptime)*100;
	    if(locdtime != Pudtime){
	      // Assume we have gained a day
	      double tdiff = (loctime + 86400 - Ptime)*100;
	    }
	    double tDdelta = (locpos.Mag()/3);
	    double tPdelta = (promptPos.Mag()/3);
	    // std::cout<<"inner loop"<<j[0]<<std::endl;
	    if(tdiff >= 0 && tdiff < mintime && tdiff < 50000){
	      if(locnhit > 50){
		// Check the relative number of hits in the event.
		if(locvalid==1 && tDdelta > tPdelta){
		  mintime = (loctime - Ptime)*100;
		  minnhit = locnhit;
		  minj = (*j)[0];
		  itrmin = k;
		}
	      }
	    }
	    k++;
	  }
	  if(ltest % 100 == 0){
	    
	    // std::cout<<"Length of index vector is "<<indices1.size()<<std::endl;
	    // 
	    
	    // std::cout<<"Length of index vector is "<<indices1.size()<<std::endl;
	    std::cout<<ltest<<" Matched events "<<(*l)[0]<<" and "<<minj<<" with tdiff = "
		     <<mintime<<", nhit1 = "<<Pnhit<<", nhit2 = "<<minnhit<<std::endl;
	  }
	}
	bool classifierdat = false;
	if(fitValid>0){
	  classifierdat = true;
	}
	tree->GetEntry((*l)[0]);
	
	double ltime = (double(ustime) + double(unstime)*1e-9)/1e-6;
	histos->SetPromptData(fitValid, meanTime,
			      ltime, nhit, energy, beta14, 
			      thetaij, itr, TVector3(posx, posy, posz),
			      TVector3(dirx, diry, dirz));
	
	if(classifierdat){
	  histos->SetClassifierData(biPolike212, biPolike214, biPoCumul,
				    alphaBeta212, alphaBeta214);
	}
	/* histos->SetPromptSelectedData(ltime, nhit, energy, beta14, 
	   thetaij, TVector3(posx, posy, posz),
	   TVector3(dirx, diry, dirz));    
	*/
	
	if (minj > 0){
	  tree->GetEntry(minj);
	  double jtime = (double(ustime) + double(unstime)*1e-9)/1e-6;
	  histos->SetDelayedData(fitValid, meanTime,
				 jtime, nhit, energy, beta14, 
				 thetaij, itr, TVector3(posx, posy, posz), 
				 TVector3(dirx, diry, dirz));
	  /*
	    if(IsEventGood())
	    histos->SetDelaySelectedData(jtime, nhit, energy, beta14, 
	    thetaij, TVector3(posx, posy, posz),
	    TVector3(dirx, diry, dirz));    
	  */
	  histos->FillHistograms(classifierdat);
	  
	  // drop the two events of interest from the search list.
	  // if(ltest % 100 == 0 ) std::cout<<pdiff<<std::endl;
	  /*
	  auto dev = find(indices0.begin(), indices0.end(), 
			  *(indices1.begin() + itrmin));
	  indices0.erase(dev);
	  indices1.erase(indices1.begin() + itrmin);
	  auto pev = find(indices1.begin(), indices1.end(), (*l));
	  // int pdiff = pev - indices1.begin();
	  indices1.erase(pev);
	  */
	  // Don't want a delayed event used as a prompt event
	  // indices0.erase(indices0.begin() + ltest);
	}
	else if(classifierdat && Pcal==0){
	  histos->FillAntiHistograms();
	}
	mintime=999999999;
	
      }
      ++l; 
      ltest++;
    }
    while(l < indices0.end()-10 && indices1.size() > 0);
  }
  histos->WriteHistograms();
  
}

int main(int argc, char** argv){
  if(argc > 5){
    std::string ofilename = argv[1];
    double xpos = atof(argv[2]);
    double ypos = atof(argv[3]);
    double zpos = atof(argv[4]);
    int forceCal = atoi(argv[5]);
    std::vector<std::string> targetruns;
    for (int i=6;i<argc;i++){
      targetruns.push_back(argv[i]);
    }
    cycleTree(targetruns, ofilename, TVector3(xpos, ypos, zpos), forceCal);
    
  }
  else 
    std::cout<<"Insufficent arguments to run program."<<std::endl;
  return 0;
}
