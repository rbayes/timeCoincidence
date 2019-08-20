#include <vector>
#include <string>

#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"



class HistContainer{
 public:
  HistContainer(std::string outfile, TVector3 sourcepos);
  ~HistContainer();

  void SetPromptClassifierData(double bipolike212, double bipolike214,
			 double bipocumulat, double alphabeta212,
			 double alphabeta214);
  void SetDelayedClassifierData(double bipolike212, double bipolike214,
			 double bipocumulat, double alphabeta212,
			 double alphabeta214);
  bool IsEventGood(int fitValid);
  void SetPromptData(int fitValid, double meanTime,
		     double utime, int nhit, double energy, double beta14, 
		     double thetaij, double itr, TVector3 pos, TVector3 dir);
  void SetDelayedData(int fitValid, double meanTime,
		      double utime, int nhit, double energy, double beta14, 
		      double thetaij, double itr, TVector3 pos, TVector3 dir);
  
  void FillHistograms(bool fillclassifier);
  void FillAntiHistograms();
  void WriteHistograms();

 private:
  TVector3 sourcepos;
  TVector3 _pos;
  TVector3 _dir;
  int _fitValid;
  double _posx;
  double _posy;
  double _posz;
  double _itr;
  double _mtime;
  int _nhit;
  double _utime;
  double _energy;
  double _beta14;
  double _thetaij;
  bool _GoodPEvent;
  bool _GoodDEvent;
 
  int    _Pnhit;
  double    _Putime;
  double _Penergy;
  double _Pbeta14;
  double _Pthetaij;
  double _Pitr;
  TVector3 _Ppos;
  double _Pposx;
  double _Pposy;
  double _Pposz;
  TVector3 _Pdir;

  
  int _DfitValid;
  int    _Dnhit;
  double _Dutime;
  double _Denergy;
  double _Dbeta14;
  double _Dthetaij;
  double _Ditr;
  TVector3 _Dpos;
  double _Dposx;
  double _Dposy;
  double _Dposz;
  TVector3 _Ddir;
  double _TimeDiff;
  
  double _bipolike212;
  double _bipolike214;
  double _bipocumulat;
  double _alphabeta212;
  double _alphabeta214;
    
  double d_bipolike212;
  double d_bipolike214;
  double d_bipocumulat;
  double d_alphabeta212;
  double d_alphabeta214;
  
  TFile* ofile;

  TTree* promptT;
  TTree* delayedT;
  TTree* unmatchT;
  
  TH1D* TimeDiff;
  TH2D* TimeDiffZD;
  TH1D* PromptNhit;
  TH1D* DelayedNhit;
  TH1D* PromptEnergy;
  TH1D* DelayedEnergy;
  TH1D* PromptBeta14;
  TH1D* DelayedBeta14;
  TH1D* PromptThetaij;
  TH1D* DelayedThetaij;
  TH1D* PromptUdotR;
  TH1D* DelayedUdotR;

  TH1D* STimeDiff;
  TH2D* STimeDiffZD;
  TH1D* SPromptNhit;
  TH1D* SDelayedNhit;
  TH1D* SPromptEnergy;
  TH1D* SDelayedEnergy;
  TH1D* SPromptBeta14;
  TH1D* SDelayedBeta14;
  TH1D* SPromptThetaij;
  TH1D* SDelayedThetaij;
  TH1D* SPromptUdotR;
  TH1D* SDelayedUdotR;

  TH1D* NPromptNhit;
  TH1D* NPromptEnergy;
  TH1D* NPromptBeta14;
  TH1D* NPromptThetaij;
  TH1D* NPromptUdotR;

  TH1D* Nhit;
  TH1D* Energy;
  TH1D* Beta14;
  TH1D* Thetaij;
  TH1D* HUdotR;
  
  TH1D* BiPoLike212;
  TH2D* EBiPoLike212;
  TH1D* BiPoLike214;
  TH2D* EBiPoLike214;
  TH1D* BiPoCumulat;
  TH2D* EBiPoCumulat;
  TH1D* AlphaBeta212;
  TH2D* EAlphaBeta212;
  TH1D* AlphaBeta214;
  TH2D* EAlphaBeta214;

  TH1D* BiPoLike212_N;
  TH2D* EBiPoLike212_N;
  TH1D* BiPoLike214_N;
  TH2D* EBiPoLike214_N;
  TH1D* BiPoCumulat_N;
  TH2D* EBiPoCumulat_N;
  TH1D* AlphaBeta212_N;
  TH2D* EAlphaBeta212_N;
  TH1D* AlphaBeta214_N;
  TH2D* EAlphaBeta214_N;

  double UdotR(bool isPrompt){ 
    return isPrompt ? _Pdir.Angle(_Ppos - sourcepos): _Ddir.Dot(_Dpos - sourcepos);}

};
