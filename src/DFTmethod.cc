
#include "DFTmethod.h"

using namespace std;

int N0;
int M0;

double xx0[2000];
double yy0[2000];

double wbin0;

SPEResponse spef0;

double fit_func( const double *x )
{
  double result = 0.0;
  
  double Norm = x[0]; 
  
  double params0[20];
  for( Int_t i=0; i<spef0.nparams; i++ ) params0[i] = x[i+1];
  spef0.SetParams( params0 );

  
  fftw_plan FWfft;
  Double_t wfin[N0]; fftw_complex wfout[M0];
    
  for ( Int_t i=0; i<N0; i++ )
    {
      Double_t xx = xx0[i];
      wfin[i] = Norm*wbin0*spef0.GetValue( xx );
            
    }
  
  FWfft = fftw_plan_dft_r2c_1d( N0, wfin, wfout, FFTW_ESTIMATE );
  fftw_execute( FWfft );
  fftw_destroy_plan( FWfft );
  
  
  Double_t fftout[N0];
  fftw_plan BWfft;
  BWfft = fftw_plan_dft_c2r_1d( N0, wfout, fftout, FFTW_ESTIMATE );
  fftw_execute( BWfft );
  fftw_destroy_plan( BWfft );
  
  for ( Int_t i=0; i<N0; i++ )
    {
      fftout[i] = fftout[i]/( 1.0*N0*1.0 );
      
    }
    
  for ( Int_t i=0; i<N0; i++ )
    {
      //if ( fftout[i]>1.0e-9 ) result += pow( fftout[i]-yy0[i], 2.0 )/( fftout[i] );
      if ( yy0[i]>0.0 ) result += 2.0*( yy0[i] - fftout[i] + fftout[i]*TMath::Log( fftout[i]/yy0[i] ) );
      else result += 2.0*( yy0[i] - TMath::Abs( fftout[i] ) ); 
            
    }
  
  return result;
  
}


ClassImp( DFTmethod )

DFTmethod::DFTmethod()
{}

DFTmethod::~DFTmethod()
{}

DFTmethod::DFTmethod( SPEResponse _spef )
{
  spef = _spef;
  spef0 = _spef;
    
}

void DFTmethod::SetHisto( TH1D* _h )
{
  hdist = _h;

  wbin = hdist->GetXaxis()->GetBinWidth(1);

  N = hdist->GetXaxis()->GetNbins();
  M = N/2+1;
    
  //cout << " N : " << N << endl;
  //cout << " M : " << M << endl;
  
  return;
  
};

void DFTmethod::CalculateValues()
{
  fftw_plan FWfftBG;
  fftw_plan FWfftSG;
  
  Double_t wfinBG[N]; fftw_complex wfoutBG[M];
  Double_t wfinSG[N]; fftw_complex wfoutSG[M];
    
  for ( UInt_t i=0; i<N; i++ )
    {
      Double_t xx = hdist->GetXaxis()->GetBinCenter( i+1 );

      Double_t arg = 0.0; 
      if ( s0!=0.0 ) arg = ( xx - Q0 )/s0;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t yy = 1.0/( sqrt( 2.0 * TMath::Pi() ) * s0 ) * TMath::Exp( -0.5*arg*arg );
      wfinBG[i] = yy;
      
      wfinSG[i] = spef.GetValue( xx );

    }
  
  FWfftBG = fftw_plan_dft_r2c_1d( N, wfinBG, wfoutBG, FFTW_ESTIMATE );
  fftw_execute( FWfftBG );
  fftw_destroy_plan( FWfftBG );

  FWfftSG = fftw_plan_dft_r2c_1d( N, wfinSG, wfoutSG, FFTW_ESTIMATE );
  fftw_execute( FWfftSG );
  fftw_destroy_plan( FWfftSG );
  

  fftw_complex wfout[M];
  Double_t fftout[N];
  
  for ( UInt_t i=0; i<M; i++ )
    {
      Double_t amp_BG = sqrt( pow( wfoutBG[i][0], 2.0 )+pow( wfoutBG[i][1], 2.0 ) );
      Double_t amp_SG = sqrt( pow( wfoutSG[i][0], 2.0 )+pow( wfoutSG[i][1], 2.0 ) );
      
      Double_t ph_BG = fftPhase( wfoutBG[i][1], wfoutBG[i][0] );
      Double_t ph_SG = fftPhase( wfoutSG[i][1], wfoutSG[i][0] );
            
      //wfout[i][0] = amp_BG*( TMath::Cos( ph_BG )  + mu*amp_SG*TMath::Cos( ph_BG+ph_SG ) ); 
      //wfout[i][1] = amp_BG*( TMath::Sin( ph_BG )  + mu*amp_SG*TMath::Sin( ph_BG+ph_SG ) );
      
      double ph = ( ph_BG+ph_SG );
      double ph0 = fmod( ph, 2.0*TMath::Pi() );
      if ( ph0>TMath::Pi() ) { ph0=ph0-2.0*TMath::Pi(); }
      
      wfout[i][0] = mu*amp_BG*amp_SG*TMath::Cos( ph0 );
      wfout[i][1] = mu*amp_BG*amp_SG*TMath::Sin( ph0 );
      
      //wfout[i][0] = TMath::Exp( mu*wfoutSG[i][0] )*( wfoutBG[i][0]*TMath::Cos( mu*wfoutSG[i][1] ) - wfoutBG[i][1]*TMath::Sin( mu*wfoutSG[i][1] ) );
      //wfout[i][1] = TMath::Exp( mu*wfoutSG[i][0] )*( wfoutBG[i][0]*TMath::Sin( mu*wfoutSG[i][1] ) + wfoutBG[i][1]*TMath::Cos( mu*wfoutSG[i][1] ) );

    }
  
  
  fftw_plan BWfft;
  BWfft = fftw_plan_dft_c2r_1d( N, wfout, fftout, FFTW_ESTIMATE );
  fftw_execute( BWfft );
  fftw_destroy_plan( BWfft );
  
  values.clear();
  for ( UInt_t i=0; i<N; i++ )
    {
      values.push_back( Norm * wbin * TMath::Exp( -1.0*mu ) * fftout[i]/( 1.0*N*1.0 ) );
      
    }
  
  return;

}

TGraph* DFTmethod::GetGraph()
{
  CalculateValues();
  
  const UInt_t nsize = hdist->GetXaxis()->GetNbins();
  cout << nsize << endl;
  
  Double_t x[nsize];
  Double_t y[nsize];

  for ( UInt_t i=0; i<nsize; i++ )
    {
      x[i] = hdist->GetXaxis()->GetBinCenter( i+1 );
      
      if ( values.at(i)<1.0e-10 ) y[i] = 1.e-3;
      else y[i] = values.at( i );
      //cout << i << ", " << x[i] << ", " << y[i] << endl;

    }
  
  TGraph *_gr = new TGraph( nsize, x, y );
  
  return _gr;
  
}

void DFTmethod::Fit()
{
  N0 = N;
  M0 = M;

  wbin0 = wbin;
  
  for ( UInt_t i=0; i<N; i++ )
    {
      xx0[i] = hdist->GetXaxis()->GetBinCenter(i+1);
      yy0[i] = hdist->GetBinContent(i+1);
      
    }
  
  mFFT = new ROOT::Minuit2::Minuit2Minimizer();
  
  ROOT::Math::Functor FCA;
  FCA = ROOT::Math::Functor( &fit_func, spef.nparams+1 );
  
  mFFT->SetFunction(FCA);
    
  mFFT->SetLimitedVariable( 0, "Norm", 1.0e+6, 1.0e+6*0.01, 1.0e+6*0.1, 1.0e+6*10.0 );
  mFFT->SetLimitedVariable( 1, "Q", 50.0, 0.1, 10.0, 80.0 );
  mFFT->SetLimitedVariable( 2, "s", 8.0, 0.1, 2.0, 40.0 );
  mFFT->SetLimitedVariable( 3, "lambda", 20.0, 0.1, 5.0, 40.0 );
  mFFT->SetLimitedVariable( 4, "w", 0.15, 0.01, 0.001, 0.8 );
  
  mFFT->SetMaxFunctionCalls(1.E9);
  mFFT->SetMaxIterations(1.E9);
  mFFT->SetTolerance(0.01);
  mFFT->SetStrategy(2);
  mFFT->SetErrorDef(1.0);
  mFFT->Minimize();
  mFFT->Hesse();
  
  Int_t ifits = 0;
  while( mFFT->Status() != 0 && ifits < 5 )
    { 
      mFFT->Minimize();
      mFFT->Hesse();
      ifits++;
      
    }
  
  if( mFFT->Status() != 0 )
    {
      cout << "" << endl;
      cout << " Fit has failed ! " << endl;
      return;

    }
  

  cout << " * " << endl;
  cout << " * Minimization results, " << mFFT->NCalls() << endl;
  cout << " * " << endl;
  
  Int_t ndim = mFFT->NDim();
  const double *pars = mFFT->X();  
  const double *erpars = mFFT->Errors();
    
  for ( int i=0; i<ndim; i++ )
    {
      cout << " * " << setw(10)  << mFFT->VariableName(i) << " : " << Form( "%.3f", pars[i] ) << " +/- " << Form( "%.3f", erpars[i] ) << endl; 
      cout << " * " << endl;

      //vals[i]=pars[i];
      //errs[i]=erpars[i];
            
    }
  
  cout << " * " << setw(10) << "chi2/NDOF" << " : " << Form( "%.2f", mFFT->MinValue()/( N-spef.nparams-1 ) ) << endl;
  cout << " * " << endl;
  
  delete pars;
  delete erpars;
  
  cout << "" << endl;
  

  return;
  
}

Double_t DFTmethod::fftPhase( Double_t vy, Double_t vz )
{
  Double_t thetayz = -999.0;

  if ( vz>0 && vy>0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); }

  else if ( vz<0 && vy>0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); thetayz=3.14159-thetayz; }

  else if ( vz<0 && vy<0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); thetayz=thetayz+3.14159; }

  else if ( vz>0 && vy<0 ) { Double_t ratio=TMath::Abs( vy/vz ); thetayz=TMath::ATan( ratio ); thetayz=2.0*3.14159-thetayz; }

  else if ( vz==0 && vy>0 ) { thetayz=3.14159/2.0; }

  else if ( vz==0 && vy<0 ) { thetayz=3.0*3.14159/2.0; }

  else if ( vz>0 && vy==0 ) { thetayz=0.0; }

  else if ( vz<0 && vy==0 ) { thetayz=3.14159; }

  if ( thetayz>3.14159 ) { thetayz=thetayz-2.0*3.14159; }
  
  return thetayz;

}
