#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

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

    G4LogicalVolume* experimentalHall_log;
    G4VPhysicalVolume* experimentalHall_phys;

    G4VPhysicalVolume* detPackage_phys[2];	// Will be an array later

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

//    G4double getScintFacePos() const { return fScintFacePos; }
//    G4double GetScintWidth() const { return fN2_volume_Z; }

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

// ----- Below are public variables from Wirechamber construction
//    G4double GetWChamWidth() const { return 2*fmwpcContainer_halfZ; }

    G4double fWChamWindowThick; 		///< mylar window thickness
    G4double fmwpc_entrance_R; 			///< entrance window radius
    G4double fmwpc_exit_R; 			///< exit window radius

//    G4Material* fMWPCGas; 			///< MWPC fill gas 	// LIKELY SAME GAS AS USED BEFORE
    G4double fEntranceToCathodes; 		///< entrance-window-to-cathode distance
    G4double fExitToCathodes; 			///< exit-window-to-cathode distance

    G4LogicalVolume* WCham_container_log; 	///< overall gas box
    G4LogicalVolume* winIn_log; 		///< inner window
    G4LogicalVolume* winOut_log; 		///< outer window
    G4LogicalVolume* kevContainer_log; 		///< container volume for kevlar strip array
    G4LogicalVolume* kevSeg_log; 		///< one segment of kevlar strip array
    G4LogicalVolume* kevStrip_log; 		///< kevlar strip in one segment

	/// electromagnetic field	- this is causing issues. Not sure why.
    void GetFieldValue(G4double Point[4], G4double* Bfield) const;
	/// whether the field changes particle energy
//    virtual G4bool DoesFieldChangeEnergy() const { return fE0 != 0; }
	/// set up tracking in field
//    void ConstructField();
	/// set anode voltage
//    void setPotential(G4double Vanode);

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

// ---- Below are protected variables from Wirechamber construction
    G4double fmwpcContainer_halfZ; 		///< half-width of wirechamber
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
    void DefineMaterials();
    string Append(int i, string str);

    float fScintStepLimit;

    G4ThreeVector fSourceHolderPos;	// here and below is returning to Mendenhall's DetectorConstruction class
    G4double fMWPCBowing;
    G4double fDetRot;
    G4ThreeVector fDetOffset;

//    float fCrinkleAngle;		// Decay trap foil crinkle angle. NOT USING

};

#endif

