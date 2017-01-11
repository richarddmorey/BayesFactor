#include "bfcommon.h"
#include "logRepresentedReal.h"
using namespace Rcpp;


// double sum overload
logRepresentedReal logRepresentedReal::operator+(const double right) const 
{
  if( right == 0.0) return *this;
  return *this + logRepresentedReal( std::log(std::abs(right)), sgn(right) );
}

// double minus overload
logRepresentedReal logRepresentedReal::operator-(const double right) const 
{
  if( right == 0.0) return *this;
  return *this + logRepresentedReal( std::log(std::abs(right)), -sgn(right) );
}

// Sum overload
logRepresentedReal logRepresentedReal::operator+(const logRepresentedReal& right) const
{
    logRepresentedReal result;
    
    if( isZero() ) return right;
    if( right.isZero() ) return *this;
    
    double x1m = modulo();
    int x1s = sign();
    double x2m = right.modulo();
    int x2s = right.sign();
       
    if( x1s == -1 && x2s == -1 ){
      result = this->negative() + right.negative();
      return result.negative();
    }
    
    if( x1s == 1 && x2s == -1 ){
      return *this - right.negative();
    }

    if( x1s == -1 && x2s == 1 ){
      return right - this->negative();
    }
    
    // Both are positive
    result = logRepresentedReal( logExpXplusExpY( x1m, x2m ), 1 );
    return result;

}

logRepresentedReal logRepresentedReal::operator-(const logRepresentedReal& right) const
{
    logRepresentedReal result;
    
    if( isZero() ) return right.negative();
    if( right.isZero() ) return *this;
    
    double x1m = modulo();
    int x1s = sign();
    double x2m = right.modulo();
    int x2s = right.sign();
        
    if( x1s == -1 && x2s == -1 ){ 
      return this->abs().negative() + right.abs();
    }
    
    if( x1s == 1 && x2s == -1 ){
      return *this + right.negative();
    }

    if( x1s == -1 && x2s == 1 ){
      result = right + this->negative();
      return result.negative();
    }

    // Both are positive
    
    if( *this > right ){
      return logRepresentedReal( logExpXminusExpY(x1m, x2m), 1 );  
    }else if( *this < right ){
      return logRepresentedReal( logExpXminusExpY(x2m, x1m) , -1);
    }
    
    // They must be equal; return 0.
    return logRepresentedReal( 0, 0 );
}

logRepresentedReal logRepresentedReal::operator*(const logRepresentedReal& right) const
{
  return logRepresentedReal( this->modulo() + right.modulo(), this->sign() * right.sign() );
}

logRepresentedReal logRepresentedReal::operator/(const logRepresentedReal& right) const
{
  return logRepresentedReal( this->modulo() - right.modulo(), this->sign() * right.sign() );
}

logRepresentedReal logRepresentedReal::operator*(const double right) const
{
  return logRepresentedReal( this->modulo() + std::log(std::abs(right)), this->sign() * sgn(right) );
}

logRepresentedReal logRepresentedReal::operator/(const double right) const
{
  return logRepresentedReal( this->modulo() - std::log(std::abs(right)), this->sign() * sgn(right) );
}

bool logRepresentedReal::operator==(const logRepresentedReal& right) const
{
  
    // Both are zero
    if( this->isZero() && right.isZero() ) { return true; }
    
    // One is zero (but not both)
    if( this->isZero() || right.isZero() ) { return false; }

    // Signs are different
    if( this->sign() != right.sign() ) { return false; }
    
    // Modulos are different
    if( this->modulo() != right.modulo() ) { return false; }

    return true;

}

bool logRepresentedReal::operator>(const logRepresentedReal& right) const
{
  
    // Both are equal
    if( *this == right ) { return false; }
    
    // signs are different
    if( this->sign() > right.sign() ) { return true; }
    if( this->sign() < right.sign() ) { return false; }
    
    // Both are positive
    if( this->sign() > 0){
      return this->modulo() > right.modulo();
    }
    
    // If we're here, both must be negative.
    return this->modulo() < right.modulo();  
}

bool logRepresentedReal::operator<(const logRepresentedReal& right) const
{
    return right > *this;
}

bool logRepresentedReal::operator<=(const logRepresentedReal& right) const
{
    return !(*this > right);
}

bool logRepresentedReal::operator>=(const logRepresentedReal& right) const
{
    return !(*this < right);
}



