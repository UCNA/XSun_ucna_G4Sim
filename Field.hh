#ifndef FIELD_HH
#define FIELD_HH

#include <vector>

#include <globals.hh>
#include <G4MagneticField.hh>
#include <G4ThreeVector.hh>

using namespace std;

class Field: public G4MagneticField {
  public:
    Field(const G4String& filename = "");	// constructor
    void GetFieldValue( const  G4double Point[3], G4double *Bfield ) const;	// get field at given point
    void SetFieldScale(G4double val) { fieldScale = val; }	// set fieldmap scaling factor
    void SetAFPDipole(G4double val) { afp_m = val; }	// set AFP dipole fringe
    void LoadFieldMap(const G4String& filename);	// load fieldmap from file

  private:
    void addPoint(G4double z, G4double B) { Zpoints.push_back(z); Bpoints.push_back(B); }	// add a point to field profile
    vector<G4double> Bpoints;		///< field profile B values
    vector<G4double> Zpoints;		///< field profile z positions
    G4double rmax2;					///< max radius squared (position in world volume) to apply field

    void addAFPFringeField(const G4double Point[3], G4double *Bfield) const;	// add AFP fringe field construction

    G4double fieldScale;			///< scaling factor for field strength
    G4double afp_m;					///< magnitude of AFP fringe dipole contribution, in A * m^2 (~14600 A*m^2)
};

#endif

