#include "HistContainer.h"


HistContainer::HistContainer(std::string outfile, TVector3 sourcepos)
{
  ofile = new TFile(outfile.c_str(), "RECREATE");

  promptT   = new TTree("promptT", "Prompt Time Coincidence Events");
  promptT->Branch("sourcepos",&sourcepos, "x:y:z");
  promptT->Branch("posx", &_Pposx, "posx/D");
  promptT->Branch("posy", &_Pposy, "posy/D");
  promptT->Branch("posz", &_Pposz, "posz/D");
  promptT->Branch("dir", &_Pdir, "x:y:z");
  promptT->Branch("nhit", &_Pnhit, "nhit/I");
  promptT->Branch("energy", &_Penergy, "energy/D");
  promptT->Branch("utime", &_Putime, "utime/D");
  promptT->Branch("beta14", &_Pbeta14, "beta14/D");
  promptT->Branch("thetaij", &_Pthetaij, "thetaij/D");
  promptT->Branch("itr", &_Pitr, "itr/D");
  promptT->Branch("bipolike212", &_bipolike212, "bipolike212/D");
  promptT->Branch("bipolike214", &_bipolike214, "bipolike214/D");
  promptT->Branch("bipocumulat", &_bipocumulat, "bipocumulat/D");
  promptT->Branch("alphabeta212",&_alphabeta212, "alphabeta212/D");
  promptT->Branch("alphabeta214",&_alphabeta214, "alphabeta214/D");
  promptT->Branch("timediff",&_TimeDiff, "timediff/D");
  delayedT  = new TTree("delayedT","Delayed Time Coincidence Events");
  delayedT->Branch("sourcepos",&sourcepos, "x:y:z");
  delayedT->Branch("posx", &_Dposx, "posx/D");
  delayedT->Branch("posy", &_Dposy, "posy/D");
  delayedT->Branch("posz", &_Dposz, "posz/D");
  delayedT->Branch("dir", &_Ddir, "x:y:z");
  delayedT->Branch("nhit", &_Dnhit, "nhit/I");
  delayedT->Branch("energy", &_Denergy, "energy/D");
  delayedT->Branch("utime", &_Dutime, "utime/D");
  delayedT->Branch("beta14", &_Dbeta14, "beta14/D");
  delayedT->Branch("thetaij", &_Dthetaij, "thetaij/D");
  delayedT->Branch("itr", &_Ditr, "itr/D");
  delayedT->Branch("bipolike212", &_bipolike212, "bipolike212/D");
  delayedT->Branch("bipolike214", &_bipolike214, "bipolike214/D");
  delayedT->Branch("bipocumulat", &_bipocumulat, "bipocumulat/D");
  delayedT->Branch("alphabeta212",&_alphabeta212, "alphabeta212/D");
  delayedT->Branch("alphabeta214",&_alphabeta214, "alphabeta214/D");
  delayedT->Branch("timediff",&_TimeDiff, "timediff/D");
  
  unmatchT  = new TTree("unmatchT","Unmatched Calibration Events");
  unmatchT->Branch("sourcepos",&sourcepos, 64000);
  unmatchT->Branch("posx", &_posx, "posx/D");
  unmatchT->Branch("posy", &_posy, "posy/D");
  unmatchT->Branch("posz", &_posz, "posz/D");
  unmatchT->Branch("dir", &_dir, 64000);
  unmatchT->Branch("nhit", &_nhit, "nhit/I");
  unmatchT->Branch("energy", &_energy, "energy/D");
  unmatchT->Branch("utime", &_utime, "utime/D");
  unmatchT->Branch("beta14", &_beta14, "beta14/D");
  unmatchT->Branch("thetaij", &_thetaij, "thetaij/D");
  unmatchT->Branch("itr", &_itr, "itr/D");
  unmatchT->Branch("bipolike212", &_bipolike212, "bipolike212/D");
  unmatchT->Branch("bipolike214", &_bipolike214, "bipolike214/D");
  unmatchT->Branch("bipocumulat", &_bipocumulat, "bipocumulat/D");
  unmatchT->Branch("alphabeta212",&_alphabeta212, "alphabeta212/D");
  unmatchT->Branch("alphabeta214",&_alphabeta214, "alphabeta214/D");
  
  TimeDiff = new TH1D("TimeDiff",";#Delta t (#mu s)",100,0,5000);
  TimeDiffZD = new TH2D("TimeDiffZD",";#Delta t (#mu s); Z_{Delayed}",
			100,0,5000,200,-800,800);
  STimeDiff = new TH1D("STimeDiff",";#Delta t (#mu s)",100,0,5000);
  STimeDiffZD = new TH2D("STimeDiffZD",";#Delta t (#mu s); Z_{Delayed}",
			 100,0,5000,200,-800,800);
  PromptNhit = new TH1D("nhit_prompt",";T_{e} (MeV)",100,0,100);
  DelayedNhit = new TH1D("nhit_delay",";T_{e} (MeV)",100,0,100);
  PromptEnergy = new TH1D("Erec_prompt",";T_{e} (MeV)",100,0,10);
  DelayedEnergy = new TH1D("Erec_delay",";T_{e} (MeV)",100,0,10);
  PromptBeta14 = new TH1D("bt14_prompt",";#beta_{14}; Events per 0.04",100,-1,3);
  DelayedBeta14 = new TH1D("bt14_delay",";#beta_{14}; Events per 0.04",100,-1,3);
  PromptThetaij = new TH1D("thij_prompt",";#theta_{ij}; Events per 0.01",
			   628,-3.14,3.14);
  DelayedThetaij = new TH1D("thij_delay",";#theta_{ij}; Events per 0.01",
			    628,-3.14,3.14);
  PromptUdotR = new TH1D("PromptUdotR",";#vec{U} * #vec{R}",400,-8,8);
  DelayedUdotR = new TH1D("DelayedUdotR",";#vec{U} * #vec{R}",400,-1000,1000);

  SPromptNhit = new TH1D("nhit_Sprompt",";T_{e} (MeV)",100,0,100);
  SDelayedNhit = new TH1D("nhit_Sdelay",";T_{e} (MeV)",100,0,100);
  SPromptEnergy = new TH1D("Erec_Sprompt",";T_{e} (MeV)",100,0,10);
  SDelayedEnergy = new TH1D("Erec_Sdelay",";T_{e} (MeV)",100,0,10);
  SPromptBeta14 = new TH1D("bt14_Sprompt",";#beta_{14}; Events per 0.04",100,-1,3);
  SDelayedBeta14 = new TH1D("bt14_Sdelay",";#beta_{14}; Events per 0.04",100,-1,3);
  SPromptThetaij = new TH1D("thij_Sprompt",";#theta_{ij}; Events per 0.01",
			   628,-3.14,3.14);
  SDelayedThetaij = new TH1D("thij_Sdelay",";#theta_{ij}; Events per 0.01",
			    628,-3.14,3.14);
  SPromptUdotR = new TH1D("SPromptUdotR",";#vec{U} * #vec{R}",400,-8,8);
  SDelayedUdotR = new TH1D("SDelayedUdotR",";#vec{U} * #vec{R}",400,-1000,1000);

  Nhit = new TH1D("nhit",";T_{e} (MeV)",100,0,100);
  Energy = new TH1D("Erec",";T_{e} (MeV)",100,0,10);
  Beta14 = new TH1D("bt14",";#beta_{14}; Events per 0.04",100,-1,3);
  Thetaij = new TH1D("thij",";#theta_{ij}; Events per 0.01",
			   628,-3.14,3.14);
  HUdotR = new TH1D("UdotR",";#vec{U} * #vec{R}",400,-8,8);

  BiPoLike212 = new TH1D("bipolike212_TC",";Likelihood",200,-150,150);
  EBiPoLike212 = new TH2D("Ebipolike212_TC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  BiPoLike214 = new TH1D("bipolike214_TC",";Likelihood",200,-150,150);
  EBiPoLike214 = new TH2D("Ebipolike214_TC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  BiPoCumulat = new TH1D("bipocumul_TC",";Cumulative Probability",1000,0,0.01);
  EBiPoCumulat = new TH2D("Ebipocumul_TC",";Cumulative Probability;Energy (MeV)",1000,0,0.01,100,0,2.5);
  AlphaBeta212 = new TH1D("bipoAlBe212_TC",";Likelihood;",200,-150,150);
  EAlphaBeta212 = new TH2D("EbipoAlBe212_TC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  AlphaBeta214 = new TH1D("bipoAlBe214_TC",";Likelihood;",200,-150,150);
  EAlphaBeta214 = new TH2D("EbipoAlBe214_TC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  
  BiPoLike212_N = new TH1D("bipolike212_NC",";Likelihood",200,-150,150);
  EBiPoLike212_N = new TH2D("Ebipolike212_NC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  BiPoLike214_N = new TH1D("bipolike214_NC",";Likelihood",200,-150,150);
  EBiPoLike214_N = new TH2D("Ebipolike214_NC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  BiPoCumulat_N = new TH1D("bipocumul_NC",";Cumulative Probability",1000,0,0.01);
  EBiPoCumulat_N = new TH2D("Ebipocumul_NC",";Cumulative Probability;Energy (MeV)",1000,0,0.01,100,0,2.5);
  AlphaBeta212_N = new TH1D("bipoAlBe212_NC",";Likelihood;",200,-150,150);
  EAlphaBeta212_N = new TH2D("EbipoAlBe212_NC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);
  AlphaBeta214_N = new TH1D("bipoAlBe214_NC",";Likelihood;",200,-150,150);
  EAlphaBeta214_N = new TH2D("EbipoAlBe214_NC",";Likelihood;Energy(MeV)",200,-150,150,100,0,5);

  
}

void HistContainer::SetClassifierData(double bipolike212, double bipolike214,
				      double bipocumulat, double alphabeta212,
				      double alphabeta214){
  _bipolike212 = bipolike212;
  _bipolike214 = bipolike214;
  _bipocumulat = bipocumulat;
  _alphabeta212 = alphabeta212;
  _alphabeta214 = alphabeta214;
}

void HistContainer::SetPromptData(int fitValid, double meanTime,
				  double utime, int nhit, double energy, 
				  double beta14, double thetaij, double itr,
				  TVector3 pos, TVector3 dir){
  _Putime = utime;
  _Pnhit  = nhit;
  _Penergy = energy;
  _Pbeta14 = beta14;
  _Pthetaij = thetaij;
  _Pitr = itr;
  _Ppos = pos;
  _Pposx = pos[0];
  _Pposy = pos[1];
  _Pposz = pos[2];
  _Pdir = dir;
  _nhit = nhit;
  _energy = energy;
  _beta14 = beta14;
  _thetaij = thetaij;
  _utime = utime;
  _itr   = itr;
  _mtime = meanTime;
  _GoodPEvent = IsEventGood(fitValid);
}


bool HistContainer::IsEventGood(int fitValid){
  bool isGood=false;
  bool unshaded=false;
  if(fitValid==1){
    TVector3 vertpos(_pos - sourcepos);
    TVector3 vertdir(_dir); 
    if ( vertpos.Mag() < 120.0 ){
      double _costh1 = vertpos.Dot(vertdir)/(vertpos.Mag() * vertdir.Mag());
      if ( _costh1 < 0. )
	unshaded = false;
      else 
	unshaded = true;
    } else
      unshaded = true;
    isGood=unshaded;
    // Select early time
    // if(_mtime < -5 || _mtime > 10) isGood=false;
    // Check classifiers
    if ( _itr < 0.55 ) isGood = false;
    if ( _beta14 > 0.9 ) isGood = false;
    
  }

  return isGood;
}


void HistContainer::SetDelayedData(int fitValid, double meanTime,
				   double utime, int nhit, double energy, 
				   double beta14, double thetaij, double itr, 
				   TVector3 pos, TVector3 dir){
  _Dutime = utime;
  _Dnhit  = nhit;
  _Denergy = energy;
  _Dbeta14 = beta14;
  _Dthetaij = thetaij;
  _Ditr = itr;
  _Dpos = pos;
  _Dposx = pos[0];
  _Dposy = pos[1];
  _Dposz = pos[2];
  _Ddir = dir;
  _pos = pos;
  _dir = dir;
  _beta14 = beta14;
  _mtime = meanTime;
  _GoodDEvent = IsEventGood(fitValid);
}

void HistContainer::FillHistograms(bool fillclassifier){
  _TimeDiff = _Dutime - _Putime;
  promptT->Fill();
  delayedT->Fill();
  TimeDiff->Fill(_Dutime - _Putime);
  PromptNhit->Fill(_Pnhit);
  DelayedNhit->Fill(_Dnhit);
  PromptEnergy->Fill(_Penergy);
  DelayedEnergy->Fill(_Denergy);
  PromptBeta14->Fill(_Pbeta14);
  DelayedBeta14->Fill(_Dbeta14);
  PromptThetaij->Fill(_Pthetaij);
  DelayedThetaij->Fill(_Dthetaij);
  PromptUdotR->Fill(UdotR(true));
  DelayedUdotR->Fill(UdotR(false));
  TimeDiffZD->Fill(_Dutime - _Putime, _Dpos[2]);
  
  if(_GoodDEvent){
    SDelayedNhit->Fill(_Dnhit);
    SDelayedEnergy->Fill(_Denergy);
    SDelayedBeta14->Fill(_Dbeta14);
    SDelayedThetaij->Fill(_Dthetaij);
    SDelayedUdotR->Fill(UdotR(false));
  }
  if(_GoodPEvent){
    SPromptNhit->Fill(_Pnhit);
    SPromptEnergy->Fill(_Penergy);
    SPromptBeta14->Fill(_Pbeta14);
    SPromptThetaij->Fill(_Pthetaij);
    SPromptUdotR->Fill(UdotR(true));

  }
  if(_GoodPEvent && _GoodDEvent){
     STimeDiff->Fill(_Dutime - _Putime);
     STimeDiffZD->Fill(_Dutime - _Putime, _Dpos[2]);
  }

  if(fillclassifier){
    BiPoLike212->Fill(_bipolike212);
    EBiPoLike212->Fill(_bipolike212, _Penergy);
    BiPoLike214->Fill(_bipolike214);
    EBiPoLike214->Fill(_bipolike214, _Penergy);
    BiPoCumulat->Fill(_bipocumulat);
    EBiPoCumulat->Fill(_bipocumulat, _Penergy);
    AlphaBeta212->Fill(_alphabeta212);
    EAlphaBeta212->Fill(_alphabeta212, _Penergy);
    AlphaBeta214->Fill(_alphabeta214);
    EAlphaBeta214->Fill(_alphabeta214, _Penergy);
  }

}


void HistContainer::FillAntiHistograms(){
  unmatchT->Fill();
  if(_GoodPEvent){
    Nhit->Fill(_Pnhit);
    Energy->Fill(_Penergy);
    Beta14->Fill(_Pbeta14);
    Thetaij->Fill(_Pthetaij);    
    HUdotR->Fill(UdotR(true));
  }

  BiPoLike212_N->Fill(_bipolike212);
  EBiPoLike212_N->Fill(_bipolike212, _Penergy);
  BiPoLike214_N->Fill(_bipolike214);
  EBiPoLike214_N->Fill(_bipolike214, _Penergy);
  BiPoCumulat_N->Fill(_bipocumulat);
  EBiPoCumulat_N->Fill(_bipocumulat, _Penergy);
  AlphaBeta212_N->Fill(_alphabeta212);
  EAlphaBeta212_N->Fill(_alphabeta212, _Penergy);
  AlphaBeta214_N->Fill(_alphabeta214);
  EAlphaBeta214_N->Fill(_alphabeta214, _Penergy);
  
}

void HistContainer::WriteHistograms(){
  ofile->cd();
  promptT->Write();
  delayedT->Write();
  unmatchT->Write();
  
  TimeDiff->Write();
  PromptNhit->Write();
  DelayedNhit->Write();
  PromptEnergy->Write();
  DelayedEnergy->Write();
  PromptBeta14->Write();
  DelayedBeta14->Write();
  PromptThetaij->Write();
  DelayedThetaij->Write();
  PromptUdotR->Write();
  DelayedUdotR->Write();
  TimeDiffZD->Write();

  STimeDiff->Write();
  SPromptNhit->Write();
  SDelayedNhit->Write();
  SPromptEnergy->Write();
  SDelayedEnergy->Write();
  SPromptBeta14->Write();
  SDelayedBeta14->Write();
  SPromptThetaij->Write();
  SDelayedThetaij->Write();
  SPromptUdotR->Write();
  SDelayedUdotR->Write();
  STimeDiffZD->Write();
  
  Nhit->Write();
  Energy->Write();
  Beta14->Write();
  Thetaij->Write();
  HUdotR->Write();

  BiPoLike212->Write();
  EBiPoLike212->Write();
  BiPoLike214->Write();
  EBiPoLike214->Write();
  BiPoCumulat->Write();
  EBiPoCumulat->Write();
  AlphaBeta212->Write();
  EAlphaBeta212->Write();
  AlphaBeta214->Write();
  EAlphaBeta214->Write();

  BiPoLike212_N->Write();
  EBiPoLike212_N->Write();
  BiPoLike214_N->Write();
  EBiPoLike214_N->Write();
  BiPoCumulat_N->Write();
  EBiPoCumulat_N->Write();
  AlphaBeta212_N->Write();
  EAlphaBeta212_N->Write();
  AlphaBeta214_N->Write();
  EAlphaBeta214_N->Write();

  ofile->Close();
}
