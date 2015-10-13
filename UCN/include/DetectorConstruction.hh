#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <G4Material.hh>		// stole from Michael Mendenhall's code.
#include <G4Element.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

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
    G4double getHolderThick() const { return fSourceHolderThickness; }

    G4double fSourceWindowThick; 			///< source foil window single-side thickness
    G4double fSourceCoatingThick; 			///< source foil coating thickness
    G4Material* fSourceWindowMat; 			///< source foil window material
    G4Material* fSourceCoatingMat; 			///< source foil coating material

    G4LogicalVolume* container_log;
    G4LogicalVolume* window_log;
    G4LogicalVolume* coating_log[2];

    G4VPhysicalVolume* source_phys;

//  ---- Below are the public variables from Decay Trap
    G4double fTrapWindowThick;
    G4double fTrapCoatingThick;
    G4double fTrapIRtrap; 				///< decay trap IR
    G4double fTrapDecayTube_Wall; 			///< decay trap wall thickness
    G4double fTrapIRcollimator; 			///< collimator IR

    G4Material* fTrapTubeMat; 			///< decay tube material
    G4Material* fTrapCollimatorMat; 		///< collimator material
    G4Material* fTrapWindowMat; 			///< decay tube window material
    G4Material* fTrapCoatingMat; 			///< decay tube coating material

    G4LogicalVolume* decayTube_log; 		///< decay trap tube
    G4LogicalVolume* trap_win_log[2]; 		///< trap window volume
    G4LogicalVolume* mylar_win_log[2]; 		///< mylar layer of window
    G4LogicalVolume* be_win_log[2]; 		///< berillium layer of window
    G4LogicalVolume* trap_monitor_log[2]; 	///< extra event monitoring region
    G4LogicalVolume* collimator_log[2]; 	///< collimator
    G4LogicalVolume* collimatorBack_log[2]; 	///< bracket behind collimator
    G4LogicalVolume* plug_log;

// ---- Below are public variables from Scintillator Construction

    G4double getScintFacePos() const { return fScintFacePos; }
    G4double GetScintWidth() const { return fN2_volume_Z; }

    G4double fScint_Radius; 			///< scintillator disc radius
    G4double fBacking_Radius; 			///< backing veto (and overall volume) radius
    G4double fScint_thick; 			///< scintillator disc thickness
    G4double fDead_thick; 			///< dead scintillator thickness, 3 um according to Junhua's thesis
    G4double fBacking_thick; 			///< backing vecto thickness (guess)
    G4double fLightguide_thick; 		///< light guide thickness at scintillator edge (guess), sets scintillator to backing distance

    G4LogicalVolume* N2_container_log; 		///< overall container (nitrogen volume)
    G4LogicalVolume* Dscint_log; 		///< scintillator dead layer logical volume
    G4LogicalVolume* scint_log; 		///< scintillator logical volume
    G4LogicalVolume* backing_log; 		///< backing veto logical volume
    G4LogicalVolume* lightguide_log; 		///< lightguide material logical volume

// ----- Below are public variables from Wire Volume Construction
    G4double GetWireVolWidth() const { return 2*cm; }

    G4Material* fMWPCGas; 			///< MWPC fill gas
    G4double fAnode_R; 				///< anode wire radius
    G4double fCathode_R; 			///< cathode wire radius
    G4double fPlating_thick; 			///< thickness of gold plating on wires
    G4double fSpacing; 				///< wire spacing
    G4int fNbOfWires; 				///< number of wires
    G4double fPlaneSpacing; 			///< spacing between wireplanes

    G4LogicalVolume* gas_log; 			///< constructed logical volume containing wireplanes
    G4LogicalVolume* cathSeg_log; 		///< cathode "segment" containing one wire in gas
    G4LogicalVolume* anodeSeg_log; 		///< anode "segment" containing one wire in gas
    G4LogicalVolume* cathode_wire_log; 		///< cathode wires logical volume
    G4LogicalVolume* cath_plate_log; 		///< gold plating on cathode segment
    G4LogicalVolume* anode_wire_log; 		///< anode wires logical volume

  protected:
    G4LogicalVolume*  fScoringVolume;	// from B1 example

// ---- Below are private variables from Source Holder class
    G4VPhysicalVolume* window_phys;
    G4VPhysicalVolume* coating_phys[2];
    G4VPhysicalVolume* holder_phys;
    G4VPhysicalVolume* ring_phys;

    G4double fSourceHolderThickness;

// ---- Below are protected variables from Scintillator Construction
    G4double fScintFacePos; 			///< position of scintillator face in container
    G4double fN2_volume_Z; 			///< z width of N2 volume

    G4VPhysicalVolume* Dscint_phys; 		///< dead layer physical volume
    G4VPhysicalVolume* scint_phys; 		///< scintillator physical volume
    G4VPhysicalVolume* backing_phys; 		///< backing veto physical volume
    G4VPhysicalVolume* lightguide_phys; 	///< lightguide material physical volume


  private:
    void DefineMaterials();
    string Append(int i, string str);

    float fScintStepLimit;

    G4LogicalVolume* experimentalHall_log;
    G4VPhysicalVolume* experimentalHall_phys;

    G4ThreeVector fSourceHolderPos;	// here and below is returning to Mendenhall's DetectorConstruction class

//    float fCrinkleAngle;		// Decay trap foil crinkle angle. NOT USING

};

#endif

