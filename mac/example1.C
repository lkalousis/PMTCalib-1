
R__LOAD_LIBRARY(libPMTCalib)

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

#include "PMTModel.h"
#include "SPEFitter.h"

Int_t example1()
{
  time_t start;  
  
  time( &start );
  
  cout << "" << endl;
  
  cout << " The macro starts ( example1.C ) ... " << endl;

  cout << "" << endl;

  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();

  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();

  Double_t Q0 = 0.0;
  Double_t s0 = 2.0;
  Pedestal ped( Q0, s0 );
  
  Double_t Q = 40.0;
  Double_t s = 13.0;
  Double_t alpha = 1.0/8.0;
  Double_t w = 0.2;
  Double_t p[4] = { Q, s, alpha, w };
  SPEResponse gaus( PMType::GAUSS, p );

  Int_t nbins = 250;
  Double_t xmin = -20.0;
  Double_t xmax = 480.0;
    
  PMT specimen( nbins, xmin, xmax, ped, gaus );
  Double_t mu = 1.2;
  Int_t ntot = 1.0e+6;
  specimen.GenSpectrum( ntot, mu );
  specimen.GetSpectrum()->SetStats(0);
  specimen.DrawSpectrum();
  
  
  SPEFitter fit;
  
  PMTModel mod( 2.0*nbins, xmin, xmax, PMType::SIMPLEGAUSS1 );
  mod.wbin = specimen.GetSpectrum()->GetBinWidth(1);
  
  Double_t mu_test = fit.FindMu( specimen.GetSpectrum(), Q0, s0 );
  Double_t g_test = fit.FindG( specimen.GetSpectrum(), Q0, mu_test );
  
  //Double_t p_test[8] = { 1.0*ntot*1.0, Q0, s0, mu, Q, s, alpha, w };
  Double_t p_test[8] = { 1.0*ntot*1.0, Q0, s0, mu_test, g_test, 0.3*g_test, 1.0/(0.5*g_test), 0.2 };
  mod.SetParams( p_test );
  
  fit.SetPMTModel( mod );
  fit.FitwPMTModel( specimen.GetSpectrum() );
  
  Double_t p_bf[8] = { fit.vals[0], fit.vals[1], fit.vals[2], fit.vals[3], fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
  mod.SetParams( p_bf );
   
  TGraph *grBF = mod.GetGraph();
  grBF->Draw( "SAME,L" );
  
  Double_t Gtrue = ( w/alpha+(1.0-w)*Q );
  Double_t Gfit = ( fit.vals[7]/fit.vals[6]+(1.0-fit.vals[7])*fit.vals[4] ); 
  
  cout << " True Gain : " << Gtrue << endl;
  cout << " BF Gain   : " << Gfit  << endl;
  cout << " Deviation : " << ( Gfit/Gtrue - 1.0 )*100.0 << endl;
  
  cout << "" << endl;
  cout << "" << endl;
  
  
  
  
  cout << " ... the macro ends ! " << endl;
	  
  cout << "" << endl;
  
  time_t end;      
  
  time( &end );
      
  Int_t dura = difftime( end, start );      
   
  Int_t min = dura / 60; Int_t sec = dura % 60;
    
  cout << " ---> "<< Form( "%02d:%02d", min, sec ) << endl;  
    
  cout << "" << endl;

  cout << "" << endl;
    
  return 0;
  
}
