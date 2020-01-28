#include "TGenPhaseSpace.h"
#include "TRandom.h"

#include <dlfcn.h>

namespace k3pi_binning {
  /// wrapper around function pointers/dlsym to give a somewhat friendlier interface 
  template <class RETURN_TYPE, class ...IN_TYPES> class DynamicFCN;
  template <class RETURN_TYPE, class ...IN_TYPES>
    class DynamicFCN<RETURN_TYPE( IN_TYPES... )>
    {
      private:
        void* m_handle                        = {nullptr};
        RETURN_TYPE ( *m_fcn )( IN_TYPES... ) = {nullptr};

      public:
        DynamicFCN() = default; 
        DynamicFCN( const std::string& lib, const std::string& name ) :
          m_handle(dlopen( lib.c_str(), RTLD_NOW )) 
      {
        m_fcn = (RETURN_TYPE( * )( IN_TYPES... ))dlsym( m_handle, name.c_str() );
        if ( m_fcn == nullptr ) std::cout << dlerror() << std::endl;
      }
        ~DynamicFCN() = default;
        RETURN_TYPE operator()( IN_TYPES... input ) const { return ( *m_fcn )( input... ); }
    };
  
  /// generate an "event" from a set of four-vectors. The "event" format is
  /// { particle1_px, particle1_py, particle1_pz, particle1_E, particle2_px ..., particleN_E},
  std::vector<double> eventFromVectors(const std::vector<TLorentzVector>& p); 
  
  /// Generate an unweighted event from TGenPhaseSpace
  std::vector<TLorentzVector> makeUnweighted( TGenPhaseSpace& phsp );

  struct binning {
    /// The amplitude code takes as argument an "event" as described above, 
    /// and the charge of the first particle in the decay, i.e. the charge of the kaon. 
    /// The amplitudes are by convention defined for D0->K+,pi-,pi-,pi+ and Dbar0->K+,pi-,pi-,pi+ for dcs 
    /// and cf, respectively. So for Dbar0->K-,pi+,pi+,pi- and D0->K-,pi+,pi+,pi- decays, the second argument should be (-1). 
    DynamicFCN<std::complex<double>(double const*, const int&) > dcs;
    DynamicFCN<std::complex<double>(double const*, const int&) >  cf;
    std::complex<double> dcsOffset; 
    std::vector<double> binLimits;

    binning( const std::string& dcsModel, 
        const std::string& cfModel, 
        const std::complex<double>& dcsOffset,  
        const std::vector<double>& binLimits ) : 
      dcs( dcsModel, "AMP"),
      cf( cfModel, "AMP"),
      dcsOffset( dcsOffset),
      binLimits(binLimits) {}

    bool ksVeto( const double& mass ) const {
      return abs( mass - 0.497614 ) < 0.010; /// veto anything within 10 MeV of the nominal KS mass.
    } 
    int binFromRelativePhase( const double& x) const {
      for( size_t i = 0 ; i < binLimits.size() ; ++i )
      {
        if( x < binLimits[i] ) return i;
      }
      std::cout << "ERROR: I Should never get here!" << std::endl; 
      return 999;
    }

    double relPhase(const double* event, const int& q) const {
      auto d1 = dcs(event, q ) * dcsOffset ;
      auto d2 = cf(event , q );
      auto c = std::conj(d1) * d2 ; 
      return 180 * std::arg(c) / M_PI;
    }
    int bin( const std::vector<TLorentzVector>& p , const int& charge ) const {
      auto event = eventFromVectors(p);
      if( ksVeto( (p[1]+p[3]).Mag() ) || ksVeto((p[2]+p[3]).Mag() ) )return binLimits.size(); 
      return binFromRelativePhase( relPhase( event.data(), charge) );
    }
  };
}

std::vector<double> k3pi_binning::eventFromVectors( const std::vector<TLorentzVector>& p )
{
  std::vector<double> event( 4 * p.size() );
  for( int i = 0 ; i < 4; ++i ){
    for( size_t j = 0 ; j < p.size(); ++j ){
      event[ 4 *j + i] = p[j][i];
    }
  }
  return event; 
}

std::vector<TLorentzVector> k3pi_binning::makeUnweighted( TGenPhaseSpace& phsp )
{
  double wt = 1; 
  do {
    wt = phsp.Generate();
  } while( wt < gRandom->Rndm() );
  return { *phsp.GetDecay(0), *phsp.GetDecay(1), *phsp.GetDecay(2), *phsp.GetDecay(3) };
}
