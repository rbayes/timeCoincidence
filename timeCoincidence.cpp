#include "HistContainer.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

bool sortByTime(const std::vector<double> &lhs, const std::vector<double> &rhs) 
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
  int nhit=0, unstime=0, ustime=0, udtime=0, pdg1=0, pdg2=0;
  double timingPeaks=0.0, mcke1=0, mcke2=0, mcEdep=0;
  ULong64_t dcApplied, dcFlagged;
  bool fitValid=false, isCal=false;
  double utime=0; // in microseconds
  tree->SetBranchAddress("pdg1",&pdg1);
  tree->SetBranchAddress("pdg2",&pdg2);
  tree->SetBranchAddress("mcke1",&mcke1);
  tree->SetBranchAddress("mcke2",&mcke2);
  tree->SetBranchAddress("mcEdep",&mcEdep);
  tree->SetBranchAddress("time",&meanTime);
  tree->SetBranchAddress("isCal", &isCal);
  tree->SetBranchAddress("dcApplied", &dcApplied);
  tree->SetBranchAddress("dcFlagged", &dcFlagged);
  tree->SetBranchAddress("timingPeaks",&timingPeaks);
  tree->SetBranchAddress("fitValid",&fitValid);
  tree->SetBranchAddress("nhitsCleaned", &nhit);
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
  double mintime=99999999;
  int Pnhit=0, Pvalid=0, Ppdg1=0, Ppdg2=0;
  double Ptime=0, Pcal=0, Ptmean=0, 
    Pustime=0, Punstime=0, Pudtime=0;
  double Penergy=0, Pitr=0, Pbeta14=0, Pthetaij=0,
    Pbplk212=0, Pbplk214=0, Pbpcml=0, Pab212=0, Pab214=0;
  ULong64_t PdcApplied, PdcFlagged;
  double Pmcke1, Pmcke2, PmcEdep;
  std::cout<<"Tree has "<<entries<<" entries"<<std::endl;
  std::vector<std::vector<double> > indices0;
  std::vector<std::vector<double> > indices1;
  // std::vector<int> iprompt;
  // std::vector<int> idelay;
  tree->GetEntry(0);
  double firsttime = double(ustime);
  double firstday  = double(udtime);
  for (int j=0; j < entries; j++){
    tree->GetEntry(j);
    double ltime = ((double(udtime) - firstday)*86400 + double(ustime) + double(unstime)*1e-9)/1e-6;
    // 
    if(nhit>180){
      if(isCal || forceCal==1){
	if(j%1000==0) std::cout<<j<<std::endl;
	histos->SetPromptData(fitValid, meanTime,
			      ltime, nhit, energy, beta14, 
			      thetaij, itr, TVector3(posx, posy, posz),
			      TVector3(dirx, diry, dirz),  dcApplied,
			      dcFlagged, mcke1, mcke2,
			      mcEdep, pdg1, pdg2);
	histos->FillAntiHistograms();
      }
      else {
	if(j%1000==0) std::cout<<j<<std::endl;
	std::vector<double> tempt;
	utime = ((double(udtime) - firstday)*86400 + (double(ustime) + double(unstime)*1e-9) - firsttime)/ 1e-6;
	tempt.push_back(double(j));
	tempt.push_back(utime);
	tempt.push_back(double(nhit));
	tempt.push_back(ustime);
	tempt.push_back(unstime);
	tempt.push_back(udtime);
	tempt.push_back(isCal || forceCal==1 ? 1:0);
	tempt.push_back(double(fitValid));
	tempt.push_back(meanTime);
	tempt.push_back(posx);
	tempt.push_back(posy);
	tempt.push_back(posz);
	tempt.push_back(energy);
	tempt.push_back(itr);
	tempt.push_back(beta14);
	tempt.push_back(thetaij);
	tempt.push_back(biPolike212);
	tempt.push_back(biPolike214);
	tempt.push_back(biPoCumul);
	tempt.push_back(alphaBeta212);
	tempt.push_back(alphaBeta214);
	tempt.push_back(double(dcApplied));
	tempt.push_back(double(dcFlagged));
	tempt.push_back(mcke1);
	tempt.push_back(mcke2);
	tempt.push_back(mcEdep);
	tempt.push_back(pdg1);
	tempt.push_back(pdg2);
	indices0.push_back(tempt);
	// indices1.push_back(tempt);

      }
    }
  }
  if (indices0.size() > 0){
    std::sort(indices0.begin(), indices0.end(), sortByTime);
    // std::sort(indices1.begin(), indices1.end(), sortByTime);
    // for (int l=0; l<entries; l++){
    std::cout<<"indices list sorted"<<std::endl;
    indices1 = indices0;
    auto l = indices0.begin();
    int ltest=0;
    do {

      bool classifierdat = false;
      if(int((*l)[7]) >0){
	classifierdat = true;
      }
      
      
      Ptime = (*l)[1];
      Pnhit = int((*l)[2]);
      Pustime = (*l)[3];
      Punstime = (*l)[4];
      Pudtime = (*l)[5];
      Pcal  = (*l)[6];
      Pvalid = int((*l)[7]);
      Ptmean = (*l)[8];
      TVector3 promptPos((double((*l)[9])),
			 (double((*l)[10])),
			 (double((*l)[11])));
      Penergy = (*l)[12];
      Pitr    = (*l)[13];
      Pbeta14 = (*l)[14];
      Pthetaij = (*l)[15];
      Pbplk212 = (*l)[16];
      Pbplk214 = (*l)[17];
      Pbpcml   = (*l)[18];
      Pab212   = (*l)[19];
      Pab214   = (*l)[20];
      PdcApplied = ULong64_t((*l)[21]);
      PdcFlagged = ULong64_t((*l)[22]);
      Pmcke1 = (*l)[23];
      Pmcke2 = (*l)[24];
      PmcEdep = (*l)[25];
      Ppdg1 = int((*l)[26]);
      Ppdg2 = int((*l)[27]);
      // Pposx  = (*l)[9];
      // Pposy  = (*l)[10];
      // Pposz  = (*l)[11];
      // std::cout<<"Copied variables"<<std::endl;
      
      int minj = -1;
      int minnhit = 9999;
      int itrmin = -1;
      int k=0;
      
      if(Pvalid>0 && Pnhit > 180 && Pnhit < 2000) {
	histos->SetPromptData(int(Pvalid), Ptmean, Ptime,
			      int(Pnhit), Penergy, Pbeta14, 
			      Pthetaij, Pitr, promptPos,
			      TVector3(0, 0, 1), PdcApplied,
			      PdcFlagged, Pmcke1, Pmcke2, PmcEdep,
			      Ppdg1, Ppdg2);
	
	if(classifierdat){
	  histos->SetPromptClassifierData(Pbplk212, Pbplk214, Pbpcml,
					  Pab212, Pab214);
	}
	
	double tPdelta = -1.;
	if((Pnhit > 180 && indices1.size() > 0) && !Pcal){
	  // If this is not a calibration event run time coincidence
	   
	  // auto pev = find(indices1.begin(), indices1.end(), (*l));
	  // std::cout<<"Evaluating event "<<int((*pev)[0])<<std::endl;
   
	  auto j = indices1.begin();
	  auto jend = indices1.end();
	  bool sequenceDone = false, anchorFound=false;
	  // for(auto j=indices1.begin(); j==std::next(pev,3); j++){
	  while(j!=jend){ 
	    // std::cout<<(*j)[0]<<std::endl;
	    if(!anchorFound){
	      if((*j)[0]==(*l)[0]){
		anchorFound=true;
		jend=std::next(j,2);
		// j-=1;
	      } else {
		++j;
		continue;
	      }
	    }
	    double loctime = (*j)[1];
	    double locnhit = (*j)[2];
	    double locdtime = (*j)[5];
	    double locvalid = (*j)[7];
	    double loctmean = (*j)[8];
	    TVector3 locpos((double((*j)[9])),
			    (double((*j)[10])),
			    (double((*j)[11])));
	    
	    ULong64_t locdcApplied = ULong64_t((*l)[21]);
	    ULong64_t locdcFlagged = ULong64_t((*l)[22]);
	    double locmcke1 = (*l)[23];
	    double locmcke2 = (*l)[24];
	    double locmcEdep = (*l)[25];
	    double tdiff = (loctime - Ptime);
	    
	    // std::cout<<"Evaluating event "<<int((*l)[0])<<" with "
	    // <<(*j)[0]<<" tdiff = "<<tdiff<<std::endl;
	    /*
	      if(locdtime != Pudtime){
	      // Assume we have gained a day
	      double tdiff = (loctime  - Ptime);
	      }
	    */
	    double tDdelta = 1000.;
	    if (locvalid==1 ){
	      tPdelta = (promptPos - locpos).Mag();
	      if(tdiff >= 0 && tdiff < mintime && tdiff < 50000){
		if(locnhit != Pnhit){
		  // Check the relative number of hits in the event.
		  if(tDdelta > tPdelta){
		    mintime = (loctime - Ptime);
		    minnhit = locnhit;
		    minj = int((*j)[0]);
		    itrmin = k;
		  }
		}
	      }
	    }
	    k++;
	    ++j;
	  }
	  if(ltest % 1000 == 0 || minj > 0){
	    
	    // std::cout<<"Length of index vector is "<<indices1.size()<<std::endl;
	    // 
	    
	    //std::cout<<"Length of index vector is "<<indices1.size()<<std::endl;
	    std::cout<<ltest<<" Matched events "<<int((*l)[0])<<" and "<<minj<<" with tdiff = "
		     <<mintime<<" and posdiff = "<<tPdelta<<", nhit1 = "<<Pnhit<<", nhit2 = "<<minnhit<<std::endl;
	  }
	}
	// tree->GetEntry((*l)[0]);
	
	// double ltime = (double(ustime) + double(unstime)*1e-9)/1e-6;
	/* histos->SetPromptSelectedData(ltime, nhit, energy, beta14, 
	   thetaij, TVector3(posx, posy, posz),
	   TVector3(dirx, diry, dirz));    
	*/
	
	if (minj > 0 && minj < tree->GetEntries()){
	  tree->GetEntry(minj);
	  double jtime = ((double(udtime) - firstday)*86400 + (double(ustime) + double(unstime)*1e-9) - firsttime)/ 1e-6;
	  // (double(ustime) + double(unstime)*1e-9)/1e-6;
	  histos->SetDelayedData(fitValid, meanTime,
				 jtime, nhit, energy, beta14, 
				 thetaij, itr, TVector3(posx, posy, posz), 
				 TVector3(dirx, diry, dirz), dcApplied,
				 dcFlagged, mcke1, mcke2, mcEdep, pdg1, pdg2);
	  
	  if(classifierdat){
	    histos->SetDelayedClassifierData(biPolike212, biPolike214, biPoCumul,
					    alphaBeta212, alphaBeta214);
	  }
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
	  
	  while((*l)[0]<minj) ++l; 
	}
	else if(classifierdat && Pcal==0){
	  histos->FillAntiHistograms();
	}
	mintime=999999999;
	
      }
      ++l;
      ltest++;
    }
    while(l < indices0.end()-3 && indices1.size() > 0);
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
