
#include "PMTModel.h"

using namespace std;

ClassImp( PMTModel )

PMTModel::PMTModel()
{}

PMTModel::~PMTModel()
{}

PMTModel::PMTModel( Int_t _nbins, Double_t _xmin, Double_t _xmax, PMType::Model _modtype )
{
  nbins = _nbins;

  xmin = _xmin;
  xmax = _xmax;

  step = ( xmax-xmin )/( 1.0*nbins*1.0 );
  
  modtype = _modtype;
  
  if ( _modtype==PMType::SIMPLEGAUSS1 ) nparams = 8;
  if ( _modtype==PMType::SIMPLEGAUSS2 ) nparams = 8;
  
}

void PMTModel::SetParams( Double_t _params[] )
{
  for ( Int_t i=0; i<nparams; i++ )
    {
      params[i] = _params[i];
      
    }
  
  return;
  
}

Double_t PMTModel::GetValue( Double_t xx )
{
  Double_t result = -666;
  
  if ( modtype==PMType::SIMPLEGAUSS1 ) result = F1( xx );
  if ( modtype==PMType::SIMPLEGAUSS2 ) result = F2( xx );
  
  return result;
  
}

Double_t PMTModel::F1( Double_t xx )
{
  Double_t result = 0.0; 
  
  Double_t Norm = params[0];
  
  Double_t Q0 = params[1];
  Double_t s0 = params[2];
  
  Double_t mu = params[3];
    
  Double_t Q1 = params[4];
  Double_t s1 = params[5];
  
  Double_t alpha = params[6];
  Double_t w = params[7];
  
  Double_t arg0 = 0.0; 
  if ( s0!=0.0 ) arg0 = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero, 0" << endl;
  
  result += TMath::Exp( -mu )/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg0*arg0 );
  
  
  Double_t sp = sqrt( pow( s0, 2.0 ) + pow( s1, 2.0 ) );
  Double_t argp = 0.0; 
  if ( sp!=0.0 ) argp = ( xx - Q0 - Q1 )/sp;    
  else cout << "Error: The code tries to divide by zero, p" << endl;
  
  Double_t S1 = w/2.0*alpha*TMath::Exp( ( 2.0*( Q0-xx )+s0*s0*alpha )/2.0*alpha )*( 1.0-TMath::Erf( ( Q0 - xx + s0*s0*alpha )/( sqrt(2.0)*s0 ) ) );
  S1 += (1.0-w)/( sqrt( 2.0*TMath::Pi() )*sp )*TMath::Exp( -0.5*argp*argp );
  
  result += TMath::Exp( -mu )*pow( mu, 1.0 )/TMath::Factorial( 1.0 )*S1;
  
  for ( Int_t n=2; n<20; n++ )
    {
      Double_t Qn = 1.0*n*( w/alpha+(1.0-w)*Q1 );
      //Double_t Qn = 1.0*n*Q1 + w/alpha;
      Double_t sn = sqrt( pow( s0, 2.0 ) + 1.0*n*( (1.0-w)*pow(s1,2.0)+w*pow(1.0/alpha,2.0)+(1.0-w)*w*pow( Q1-1.0/alpha, 2.0 )  )   );
      //Double_t sn = sqrt( pow( s0, 2.0 ) + 1.0*n*pow( s1, 2.0 ) );
            
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx - Q0 - Qn )/sn;    
      else cout << "Error: The code tries to divide by zero, 1 " << endl;
      
      result += TMath::Exp( -mu )*pow( mu, n )/TMath::Factorial( n )/( sqrt( 2.0*TMath::Pi() )*sn ) * TMath::Exp( -0.5*argn*argn );
            
    }

  result *= Norm*wbin;
  
  return result;

}

Double_t PMTModel::F2( Double_t xx )
{
  Double_t result = 0.0; 

  Double_t Norm = params[0];
  
  Double_t Q0 = params[1];
  Double_t s0 = params[2];
  
  Double_t mu = params[3];
    
  Double_t Q1 = params[4];
  Double_t s1 = params[5];
  
  Double_t alpha = params[6];
  Double_t w = params[7];

  Double_t arg0 = 0.0; 
  if ( s0!=0.0 ) arg0 = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero, 0" << endl;
  
  result += TMath::Exp( -mu )/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg0*arg0 );
  
  
  Double_t sp = sqrt( pow( s0, 2.0 ) + pow( s1, 2.0 ) );
  Double_t argp = 0.0; 
  if ( sp!=0.0 ) argp = ( xx - Q0 - Q1 )/sp;    
  else cout << "Error: The code tries to divide by zero, p" << endl;
  
  Double_t S1 = w/2.0*alpha*TMath::Exp( ( 2.0*( Q0-xx )+s0*s0*alpha )/2.0*alpha )*( 1.0-TMath::Erf( ( Q0 - xx + s0*s0*alpha )/( sqrt(2.0)*s0 ) ) );
  S1 += (1.0-w)/( sqrt( 2.0*TMath::Pi() )*sp )*TMath::Exp( -0.5*argp*argp );

  result += TMath::Exp( -mu )*pow( mu, 1.0 )/TMath::Factorial( 1.0 )*S1;
  
  
  Double_t S2 = 0.0;//pow( w/lambda, 2.0 )*xx*TMath::Exp( -xx/lambda );

  Double_t spp = sqrt( pow( s0, 2.0 ) + 2.0*pow( s1, 2.0 ) );
  Double_t argpp = 0.0; 
  if ( spp!=0.0 ) argpp = ( xx - Q0 - 2.0*Q1 )/spp;    
  else cout << "Error: The code tries to divide by zero, p" << endl;
    
  S2 += pow( 1.0-w, 2.0 )/( sqrt( 2.0*TMath::Pi() )*spp )*TMath::Exp( -0.5*argpp*argpp );
  S2 += 2.0*w*(1.0-w)/2.0*alpha*TMath::Exp( ( 2.0*( Q1+Q0-xx )+(s0*s0+s1*s1)*alpha )/2.0*alpha )*( 1.0-TMath::Erf( ( Q1+Q0-xx + (s0*s0+s1*s1)*alpha )/( sqrt(2.0)*sqrt(s0*s0+s1*s1) ) ) );
  
  result += TMath::Exp( -mu )*pow( mu, 2.0 )/TMath::Factorial( 2.0 )*S2;
  
  for ( Int_t n=3; n<20; n++ )
    {
      Double_t Qn = 1.0*n*( w/alpha+(1.0-w)*Q1 );
      Double_t sn = sqrt( pow( s0, 2.0 ) + 1.0*n*( (1.0-w)*pow(s1,2.0)+w*pow(1.0/alpha,2.0)+(1.0-w)*w*pow( Q1-1.0/alpha, 2.0 )  )   );
      
      
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx - Q0 - Qn )/sn;    
      else cout << "Error: The code tries to divide by zero, 1 " << endl;
      
      result += TMath::Exp( -mu )*pow( mu, n )/TMath::Factorial( n )/( sqrt( 2.0*TMath::Pi() )*sn ) * TMath::Exp( -0.5*argn*argn );
      
    }

  result *= Norm*wbin;
  
  return result;

}

TGraph* PMTModel::GetGraph()
{
  Double_t x[nbins];
  Double_t y[nbins];
  
  for ( Int_t i=0; i<nbins; i++ )
    {
      x[i] = xmin + step/2 + 1.0*i*step;

      Double_t y_ = GetValue( x[i] );
      
      if ( y_<1.0e-10 ) y[i] = 1.e-3;
      else y[i] = y_;
      
    }
  
  TGraph *_gr = new TGraph( nbins, x, y );

  _gr->SetLineWidth( 2 );
  _gr->SetLineColor( kBlue );
  _gr->SetMarkerColor( kBlue );
  _gr->SetMarkerSize( 0.1 );
  
  return _gr;

}
