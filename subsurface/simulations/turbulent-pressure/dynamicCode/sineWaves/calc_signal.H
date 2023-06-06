// Import arrays and global variables
#include "fft_pairs.H"

// Function declararion
scalar signalBuild(scalar Xm, scalar t);

// Function definition
scalar signalBuild(scalar Xm, scalar t){

  // Parameter to modify:
  scalar LX     = 0.15;
  double delf   = 0.00628677737438413;
  
  // Initializations (do not change)
  double freq   = 0.0;
  scalar signal = 0.0;

  // Mapper for X-coordinate
  int xi = round((nXcoor-1) * Xm / LX);

  // Calculate pressure field
  for ( int f = 0; f < nFreqs; f++ ){
    freq = delf * f;
    signal += RealValues[xi][f] * Foam::cos(2*pi*t*freq) 
            - ImagValues[xi][f] * Foam::sin(2*pi*t*freq);    
    
  }
  
  return signal;
}