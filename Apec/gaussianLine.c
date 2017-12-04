// Add the input gaussian line to the fluxArray. Converted from Fortran
// in xsgaul.f.

//#include <xsTypes.h>
//#include <XSstreams.h>
#include <math.h>

double energyFraction( double linecentre,  double linewidth,  double energy);
int binarySearch( double *energyArray, int nbins, double energy);

void calcGaussianLine(double *energyArray, int nbins, double ecenter, double width, double lineflux, double crtsig, double *fluxArray) {

  // to calculate gaussian line shapes. If qspeedy is true then uses a tabulation
  // otherwise does a proper calculation of the erf function to +/- crtsig.

  // first trap case of zero width line


  double efirst, elast, emin, emax;
  double winv;
  int ielow = 0;
  int ie;
  double alow, ahi;

  if ( width <= 0.0 ) {
    ie = binarySearch(energyArray, nbins, ecenter);
    if ( ie >= 0 && ie < nbins ) fluxArray[ie] += lineflux;
    return;
  }

  efirst = energyArray[0];
  elast = energyArray[nbins];

  // Set the line energy, width and the energy range over which the
  // line will be calculated (essentially +/- crtsig sigma)

  emin = ecenter - crtsig*width;
  emax = ecenter + crtsig*width;

  // If the line calculation energies lie outside the energy range
  // then give up

  if ( emax < efirst || emin > elast ) return;


  // and this is the slower method using erf calculation.

  winv = 1./width/sqrt(2.0);

  // If there is part of the gaussian below the response energy range then
  // store that in alow.

  ielow = 0;
  ie = 1;

  if ( emin <= efirst ) {
    alow = 0.5*erf(winv*(efirst-ecenter));
    ielow = -1;
  } else {
    alow = -0.5;
  }

  // find the start of the line

  if ( ielow == 0 ) ielow = binarySearch(energyArray, nbins,  emin);
  if ( ielow != -1 ) ie = ielow;

  // integrate the gaussian for the energy ranges covering the line

  while ( energyArray[ie] < emax ) {
    ahi = 0.5*erf(winv*(energyArray[ie]-ecenter));
    fluxArray[ie] += lineflux*(ahi - alow);
    alow = ahi;
    ie++;
    if ( ie > nbins ) return;
  }

  // if all of the model was in a single bin

  if ( ie == ielow ) {
    fluxArray[ie] += lineflux;
  } else {
    fluxArray[ie] += lineflux*(0.5 - alow);
  }

  return;
}

double energyFraction( double linecentre,  double linewidth,  double energy)
{

  // Function to return amount of energy emitted between line centre and energy

  int GaussMaxStep = 120 ;
  int GaussMax = 6;
  double GaussStep = (double)GaussMax / GaussMaxStep;

  double intGaus[GaussMaxStep+1];
  int first = 1;
  double sigma, remainder;
  int i, minindex;

  // if the first time through then calculate the erf grid on which we
  // will interpolate

  if ( first ) {
    for (i=0; i<=GaussMaxStep; i++) {
      intGaus[i] = erf(i*GaussStep/sqrt((double)2.0));
    }
    first = 0;
  }
  
  // Compute how many sigmas away we are
  sigma = fabs( (energy-linecentre)/linewidth );

  // Now interpolate from the table
  minindex = (int)(sigma/GaussStep);

  // If we're past the edge of the tabulated data return -1.
  if ( minindex >= GaussMaxStep ) return -1;

    remainder = (sigma - (double)minindex*GaussStep) * (1.0/GaussStep);

  // Do the interpolation
  return (1.0-remainder)*intGaus[minindex] + remainder*intGaus[minindex+1];

}

int binarySearch( double *energyArray, int nbins, double energy) {

  // Function to do a bisection search for which element of the energy array
  // the input energy lies in.

  if ( energy < energyArray[0] || energy > energyArray[nbins] ) {
    return -1;
  }

  int low = 0;
  int high = nbins;

  int bisearch;
  while ( (high-low) > 1 ) {

    bisearch = (low + high) / 2;
    if ( energy > energyArray[bisearch-1] ) {
      low = bisearch;
    } else {
      high = bisearch;
    }

  }

  if ( energy > energyArray[bisearch] ) {
    bisearch = high;
  } else {
    bisearch = low;
  }

  return bisearch-1;
}


