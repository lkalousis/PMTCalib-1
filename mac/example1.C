
R__LOAD_LIBRARY(libPMTCalib)

#include "PMTStyle.h"
#include "PMType.h"

#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"

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
  Double_t s = 8.0;
  Double_t lambda = 20.0;
  Double_t w = 0.1;
  Double_t p[4] = { Q, s, lambda, w };
  SPEResponse gaus( PMType::GAUSS, p );
  cout << gaus.spefunc->Integral( 0.0, 1000.0 ) << endl; 
  
  PMT specimen( 200, -50.0, 350.0, ped, gaus );
  Double_t mu = 0.8;
  specimen.GenSpectrum( 1.0e+6, mu );
  specimen.DrawSpectrum();
  cout << Q0 + mu*( w*lambda+(1.0-w)*Q ) << endl;
  c1->Update();
  c1->WaitPrimitive();
  
  
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
