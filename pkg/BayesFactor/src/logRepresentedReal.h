#ifndef _RcppBayesFactor_LOGREAL_HPP
#define _RcppBayesFactor_LOGREAL_HPP

#include <Rcpp.h>
#include "bfcommon.h"

class logRepresentedReal {
  int s;
  double m;

  public:
  logRepresentedReal() {
    s = NA_INTEGER;
    m = NA_REAL;
	}

  logRepresentedReal(double mod, int sign ) {
		if ( std::abs(sign) > 1 )  // invalid sign
			Rcpp::stop("ERROR: sign must be -1, 0, or 1.");
    if( !R_FINITE(mod) && sgn(mod)==-1 ){
      sign = 0;
    }
    if( sign == 0 ){
      mod = R_NegInf;
    }
    s = sign;
    m = mod;
		}
	

public: // ==== USER INTERFACE =====
	const int sign() const { return s; }
  
  const double modulo() const { return m; }
  
  const double log() const {
    if( s == 1 ){
      return m;
    }else if(s == 0){
      return R_NegInf;
    }else if(s == -1){
      return NA_REAL;
    }
    
    Rcpp::stop("ERROR: Invalid sign in logRepresentedReal.");  
    return(NA_REAL);
  };
  
  logRepresentedReal abs() const {
    return logRepresentedReal(m, 1 );  
  }
  
  logRepresentedReal negative() const {
    return logRepresentedReal(m, -s );  
  } 

  logRepresentedReal reciprocal() const {
    return logRepresentedReal(-m, s);
  }

  logRepresentedReal pow( int e ) const {
    if( e == 0 ){ logRepresentedReal( 0, 1 ); } 
    if( !e%2 ){ logRepresentedReal( e * m, 1 ); }
    
    return logRepresentedReal( e * m, s );
  }

  logRepresentedReal pow( double e ) const {
    if( e == 0 ){ logRepresentedReal( 0, 1 ); } 

    return logRepresentedReal( e * m, s );
  }

  bool isZero() const {
    if( ( !R_FINITE(m) && sgn(m)==-1 ) | (s == 0) ){
      return true; 
    }else{
      return false;
    } 
  }

  operator double() { return s * exp(m); }
  
  logRepresentedReal operator+(const logRepresentedReal& right) const;
  logRepresentedReal operator-(const logRepresentedReal& right) const;
  logRepresentedReal operator*(const logRepresentedReal& right) const;
  logRepresentedReal operator/(const logRepresentedReal& right) const;

  logRepresentedReal operator+(const double right) const;
  logRepresentedReal operator-(const double right) const;
  logRepresentedReal operator*(const double right) const;
  logRepresentedReal operator/(const double right) const;

  bool operator==(const logRepresentedReal& right) const;
  bool operator>(const logRepresentedReal& right) const;
  bool operator<(const logRepresentedReal& right) const;
  bool operator<=(const logRepresentedReal& right) const;
  bool operator>=(const logRepresentedReal& right) const;


};

#endif
