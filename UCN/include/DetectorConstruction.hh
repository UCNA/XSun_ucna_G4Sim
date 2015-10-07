#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <G4Material.hh>		// stole from Michael Mendenhall's code.
#include <G4Element.hh>
#include <G4SystemOfUnits.hh>

#include <string>
#include <sstream>

using 	namespace	std;

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

    G4Material* Vacuum; 	///< our slightly crappy vacuum
    G4Material* Brass; 		///< brass for source holder
    G4Material* SS304; 		///< 304 Stainless Steel
    G4Material* Kevlar; 	///< kevlar for wirechamber window support strings
    G4Material* Mylar; 		///< mylar for windows
    G4Material* Polyethylene; 	///< poly for collimator
    G4Material* WCPentane; 	///< Wirechamber fill: (neo)pentane @ 100torr
    G4Material* WCNitrogen; 	///< Wirechamber fill: Nitrogen @ 100torr
    G4Material* Sci; 		///< scintillator material

    void setVacuumPressure(G4double pressure);

// ---- Below are public variables from Source Holder class
    /// get thickness
    G4double getHolderThick() const { return fSourceHolderThickness; }

    G4double fWindowThick; ///< source foil window single-side thickness
    G4double fCoatingThick; ///< source foil coating thickness
    G4Material* fWindowMat; ///< source foil window material
    G4Material* fCoatingMat; ///< source foil coating material

    G4LogicalVolume* container_log;
    G4LogicalVolume* window_log;
    G4LogicalVolume* coating_log[2];



  protected:
    G4LogicalVolume*  fScoringVolume;	// from B1 example

// ---- Below are private variables from Source Holder class
    G4VPhysicalVolume* window_phys;
    G4VPhysicalVolume* coating_phys[2];
    G4VPhysicalVolume* holder_phys;
    G4VPhysicalVolume* ring_phys;

    G4double fSourceHolderThickness;


  private:
    void DefineMaterials();
    string Append(int i, string str);

    float fScintStepLimit;

    G4LogicalVolume* experimentalHall_log;
    G4VPhysicalVolume* experimentalHall_phys;

};

#endif

