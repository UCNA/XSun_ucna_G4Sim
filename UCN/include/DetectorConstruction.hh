#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <G4Material.hh>		// stole from Michael Mendenhall's code.
#include <G4Element.hh>
#include <G4SystemOfUnits.hh>


const G4double inch = 2.54*cm;
const G4double torr = atmosphere/760.;

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    G4Material* Be; 		///< Beryllium for trap windows
    G4Material* Al; 		///< Aluminum
    G4Material* Si; 		///< Silicon
    G4Material* Cu; 		///< Copper for decay trap
    G4Material* Wu; 		///< Tungsten for anode wires
    G4Material* Au; 		///< Gold for cathode wires coating

    G4Material* Vacuum; ///< our slightly crappy vacuum
    G4Material* Brass; 		///< brass for source holder
    G4Material* SS304; 		///< 304 Stainless Steel
    G4Material* Kevlar; 	///< kevlar for wirechamber window support strings
    G4Material* Mylar; 		///< mylar for windows
    G4Material* Polyethylene; ///< poly for collimator
    G4Material* WCPentane; ///< Wirechamber fill: (neo)pentane @ 100torr
    G4Material* WCNitrogen; ///< Wirechamber fill: Nitrogen @ 100torr
    G4Material* Sci; ///< scintillator material

    void setVacuumPressure(G4double pressure);


  protected:
    G4LogicalVolume*  fScoringVolume;

  private:
    void DefineMaterials();

    float fScintStepLimit;

    G4LogicalVolume* experimentalHall_log;
    G4VPhysicalVolume* experimentalHall_phys;

};

#endif

