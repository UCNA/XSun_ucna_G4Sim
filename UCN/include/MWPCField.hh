#ifndef MWPCField_h
#define MWPCField_h 1

#include <vector>
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;

using namespace std;

class MWPCField: public G4MagneticField
{
public:
  MWPCField();

  void LoadFieldMap();
  void GetFieldValue( const G4double Point[4], G4double *Bfield ) const;

protected:
  G4double fE0;		// apparently a field scaling constant

private:
  void AddPoint(G4double zPositions, G4double BValues);
  vector<G4double> Bpoints; ///< field profile B values
  vector<G4double> Zpoints; ///< field profile z positions

  G4double fSqOfMaxRadius;			// These two needed for making Mag field component of global field
  double fFieldScale;				// dimensionless scaling factor

};

#endif
