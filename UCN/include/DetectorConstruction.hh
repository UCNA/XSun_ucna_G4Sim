#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "Field.hh"

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

    void SetScoringVolumes(G4LogicalVolume* vol, int type, int location);

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
    G4LogicalVolume* scint_container_log;
    G4LogicalVolume* scint_deadLayer_log;
    G4LogicalVolume* scint_scintillator_log;
    G4LogicalVolume* scint_lightGuide_log;
    G4LogicalVolume* scint_backing_log;
    G4LogicalVolume* wireVol_gas_log;
    G4LogicalVolume* wireVol_cathSeg_log;
    G4LogicalVolume* wireVol_anodeSeg_log;
    G4LogicalVolume* wireVol_cathodeWire_log;
    G4LogicalVolume* wireVol_cathPlate_log;
    G4LogicalVolume* wireVol_anodeWire_log;
    G4LogicalVolume* mwpc_container_log;
    G4LogicalVolume* mwpc_kevContainer_log;
    G4LogicalVolume* mwpc_kevSeg_log;
    G4LogicalVolume* mwpc_kevStrip_log;
    G4LogicalVolume* mwpc_winIn_log;
    G4LogicalVolume* mwpc_winOut_log;
    G4LogicalVolume* frame_mwpcEntrance_log;
    G4LogicalVolume* frame_entranceFront_log;
    G4LogicalVolume* frame_entranceMid_log;
    G4LogicalVolume* frame_entranceBack_log;
    G4LogicalVolume* frame_container_log;
    G4LogicalVolume* frame_mwpcExit_log;
    G4LogicalVolume* frame_mwpcExitGasN2_log;
    G4LogicalVolume* frame_backStuff_log;

  protected:
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* source_holder_phys;
    G4VPhysicalVolume* source_window_phys;
    G4VPhysicalVolume* source_coating_phys[2];
    G4VPhysicalVolume* source_ring_phys;
    G4VPhysicalVolume* source_phys;
    G4VPhysicalVolume* scint_deadLayer_phys;
    G4VPhysicalVolume* scint_scintillator_phys;
    G4VPhysicalVolume* scint_lightGuide_phys;
    G4VPhysicalVolume* scint_backing_phys;
    G4VPhysicalVolume* scint_container_phys;
    G4VPhysicalVolume* mwpc_container_phys;
    G4VPhysicalVolume* frame_entranceFront_phys;
    G4VPhysicalVolume* frame_entranceMid_phys;
    G4VPhysicalVolume* frame_entranceBack_phys;
    G4VPhysicalVolume* frame_mwpcExit_phys;
    G4VPhysicalVolume* frame_mwpcExitGasN2_phys;
    G4VPhysicalVolume* frame_backStuff_phys;

  private:




  public:
    G4VPhysicalVolume* detPackage_phys[2];	// Will be an array later

// ---- Below are public variables from Source Holder class
    G4double getHolderThick() const { return fSourceHolderThickness; }

// ----- Below are public variables from Wirechamber construction

	/// electromagnetic field
//    void GetFieldValue(G4double Point[4], G4double* Bfield);
	/// whether the field changes particle energy
    virtual G4bool DoesMWPCFieldChangeEnergy() const { return fE0 != 0; }
	/// set up tracking in field
    void ConstructMWPCField();
	/// set anode voltage
//    void setMWPCPotential(G4double Vanode);

    G4MagneticField* fMyBField; 			///< Magnetic field pointer
    G4RotationMatrix* fMyRotation; 		///< rotation from global frame to local coordinates

    G4ThreeVector fMyTranslation; 		///< translation from global coordinates to center of anode plane

// ---- Below are public variables from the Detector Package Construction
    G4double fDPC_detPackageRadius;
    G4double fDPC_mwpc_entrance_thickness; 	///< MWPC entrance tube wall thickness
    G4double fDPC_mwpc_entrance_r; 		///< MWPC entrance tube radius
    G4double fDPC_mwpc_entrance_depth; 		///< MWPC entrance tube depth
    G4double fDPC_frontwin_frame_thick; 	///< MWPC front window frame thickness
    G4double fDPC_backwin_frame_thick; 		///< MWPC exit window frame thickness

    G4LogicalVolume* fDPC_container_log; 	///< overall positioning container
    G4LogicalVolume* fDPC_mwpc_entrance_log; 	///< entrance port container
    G4LogicalVolume* fDPC_entrance_front_log;	///< entrance port front plate
    G4LogicalVolume* fDPC_entrance_mid_log; 	///< entrance port tube
    G4LogicalVolume* fDPC_entrance_back_log; 	///< entrance port back plate (MWPC box cover)
    G4LogicalVolume* fDPC_mwpc_exit_log; 	///< aluminum exit window from wirechamber
    G4LogicalVolume* fDPC_mwpc_exit_N2_log; 	///< N2 between exit window and scintillator
    G4LogicalVolume* fDPC_backstuff_log; 	///< miscellaneous mass behind detectors

    G4double fDPC_entrance_face_pos; 		///< entrance window port entrance relative to scint face
    G4double fDPC_entrance_win_pos; 		///< MWPC entrance window position relative to scint face
    G4double fDPC_exit_frame_pos; 		///< exit window frame pos


  protected:
// ---- Below are private variables from Source Holder class
    G4double fSourceHolderThickness;

// ---- Below are protected variables from Wirechamber construction
    G4double fE0; 				///< field scaling constant

// ---- Below are protected variables from DetectorPackageConstruction
    G4VPhysicalVolume* fDPC_scint_phys;
    G4VPhysicalVolume* fDPC_mwpc_phys;
    G4VPhysicalVolume* fDPC_mwpc_entrance_phys;
    G4VPhysicalVolume* fDPC_entrance_front_phys;
    G4VPhysicalVolume* fDPC_entrance_mid_phys;
    G4VPhysicalVolume* fDPC_entrance_back_phys;
    G4VPhysicalVolume* fDPC_mwpc_exit_phys;
    G4VPhysicalVolume* fDPC_mwpc_exit_N2_phys;
    G4VPhysicalVolume* fDPC_backstuff_phys;

  private:
    void ConstructField();

    Field* fpMagField;
    bool fCallFieldConstructor;

    void DefineMaterials();
    string Append(int i, string str);

    G4double fScintStepLimit;

    G4double fDetRot;
    G4ThreeVector fDetOffset;

    G4LogicalVolume* fScoreVol1;
    G4LogicalVolume* fScoreVol2;
    G4LogicalVolume* fScoreVol3;
    G4LogicalVolume* fScoreVol4;



};

#endif

