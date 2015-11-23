#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

//#include "Field.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <G4Material.hh>		// stole from Michael Mendenhall's code.
#include <G4Element.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

#include <G4ElectroMagneticField.hh>	// Taken from WirechamberConstruction.
#include <G4MagneticField.hh>
#include <G4RotationMatrix.hh>

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
    DetectorConstruction();		// Constructor/destructors
    virtual ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct();

    void setVacuumPressure(G4double pressure);

    G4LogicalVolume* GetScoringVolume1() const { return fScoreVol1; }	// SteppingAction needs access
    G4LogicalVolume* GetScoringVolume2() const { return fScoreVol2; }	// to scoring volumes.
    G4LogicalVolume* GetScoringVolume3() const { return fScoreVol3; }
    G4LogicalVolume* GetScoringVolume4() const { return fScoreVol4; }

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

    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* source_container_log;
    G4LogicalVolume* source_window_log;
    G4LogicalVolume* source_coating_log[2];
    G4LogicalVolume* decayTrap_tube_log;
    G4LogicalVolume* decayTrap_window_log[2];
    G4LogicalVolume* decayTrap_mylarWindow_log[2];
    G4LogicalVolume* decayTrap_beWindow_log[2];
    G4LogicalVolume* decayTrap_collimator_log[2];
    G4LogicalVolume* decayTrap_collimatorBack_log[2];
    G4LogicalVolume* decayTrap_innerMonitors_log[2];
    G4LogicalVolume* scint_container_log[2];
    G4LogicalVolume* scint_deadLayer_log[2];
    G4LogicalVolume* scint_scintillator_log[2];
    G4LogicalVolume* scint_lightGuide_log[2];
    G4LogicalVolume* scint_backing_log[2];
    G4LogicalVolume* wireVol_gas_log[2];
    G4LogicalVolume* wireVol_cathSeg_log[2];
    G4LogicalVolume* wireVol_anodeSeg_log[2];
    G4LogicalVolume* wireVol_cathodeWire_log[2];
    G4LogicalVolume* wireVol_cathPlate_log[2];
    G4LogicalVolume* wireVol_anodeWire_log[2];
    G4LogicalVolume* mwpc_container_log[2];
    G4LogicalVolume* mwpc_kevContainer_log[2];
    G4LogicalVolume* mwpc_kevSeg_log[2];
    G4LogicalVolume* mwpc_kevStrip_log[2];
    G4LogicalVolume* mwpc_winIn_log[2];
    G4LogicalVolume* mwpc_winOut_log[2];
    G4LogicalVolume* frame_mwpcEntrance_log[2];
    G4LogicalVolume* frame_entranceFront_log[2];
    G4LogicalVolume* frame_entranceMid_log[2];
    G4LogicalVolume* frame_entranceBack_log[2];
    G4LogicalVolume* frame_container_log[2];
    G4LogicalVolume* frame_mwpcExit_log[2];
    G4LogicalVolume* frame_mwpcExitGasN2_log[2];
    G4LogicalVolume* frame_backStuff_log[2];

  protected:
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* source_holder_phys;
    G4VPhysicalVolume* source_window_phys;
    G4VPhysicalVolume* source_coating_phys[2];
    G4VPhysicalVolume* source_ring_phys;
    G4VPhysicalVolume* source_phys;
    G4VPhysicalVolume* scint_deadLayer_phys[2];
    G4VPhysicalVolume* scint_scintillator_phys[2];
    G4VPhysicalVolume* scint_lightGuide_phys[2];
    G4VPhysicalVolume* scint_backing_phys[2];
    G4VPhysicalVolume* scint_container_phys[2];
    G4VPhysicalVolume* mwpc_container_phys[2];
    G4VPhysicalVolume* frame_entranceFront_phys[2];
    G4VPhysicalVolume* frame_entranceMid_phys[2];
    G4VPhysicalVolume* frame_entranceBack_phys[2];
    G4VPhysicalVolume* frame_mwpcExit_phys[2];
    G4VPhysicalVolume* frame_mwpcExitGasN2_phys[2];
    G4VPhysicalVolume* frame_backStuff_phys[2];
    G4VPhysicalVolume* frame_mwpcEntrance_phys[2];
    G4VPhysicalVolume* frame_container_phys[2];

  private:
//    void ConstructGlobalField();
//    void ConstructLocalField(int flag);

//    Field* fpMagField;
//    Field* mwpc_internalBfield_EAST;
//    Field* mwpc_internalBfield_WEST;

    void DefineMaterials();
    string Append(int i, string str);

    G4double fScintStepLimit;

    G4LogicalVolume* fScoreVol1;
    G4LogicalVolume* fScoreVol2;
    G4LogicalVolume* fScoreVol3;
    G4LogicalVolume* fScoreVol4;



  public:
// ----- Below are public variables from Wirechamber construction
	/// whether the field changes particle energy
    virtual G4bool DoesMWPCFieldChangeEnergy() const { return fE0 != 0; }
	/// set anode voltage
//    void setMWPCPotential(G4double Vanode);

  protected:
// ---- Below are protected variables from Wirechamber construction
    G4double fE0; 				///< field scaling constant
};

#endif

