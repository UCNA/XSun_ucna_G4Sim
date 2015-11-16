#include "DetectorConstruction.hh"
#include "Field.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include <G4UserLimits.hh>		// stole from Michael Mendenhall's code.

#include <G4SubtractionSolid.hh>	// taken from Source holder class

#include <globals.hh>			// taken from Detector construction utils class
#include <G4Material.hh>
#include <G4Element.hh>
#include <G4Tubs.hh>			// nothing used by decay trap construction
#include <G4VPhysicalVolume.hh>		// not using wiggle sheet
#include <G4LogicalVolume.hh>		// or silicon detector construction
#include <G4ThreeVector.hh>
#include <G4PVReplica.hh>
#include <G4RotationMatrix.hh>
#include <G4VisAttributes.hh>
#include <G4SystemOfUnits.hh>

#include <cassert>			// scintillator construction classes
#include <G4Polycone.hh>

#include <math.h>			// Used in WirechamberConstruction
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4EqMagElectricField.hh>
#include <G4ClassicalRK4.hh>

#include <Randomize.hh>			// Stolen from Analysis Manager
#include <G4ios.hh>			// Pretty sure needed for TrackerSD
#include <G4Run.hh>			// Leave them here since we use registerSD in DetectorConstruction
#include <G4Event.hh>			// And the registerSD is totally not working without it
#include <G4Track.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4IonTable.hh>
#include <G4SDManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4SDManager.hh>
#include <G4EventManager.hh>

#include <G4MagneticField.hh>		// Bottom half of detector construction
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4PropagatorInField.hh>
#include <G4TransportationManager.hh>
#include <G4UserLimits.hh>
#include <G4PVParameterised.hh>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), fpMagField(NULL),
  fScintStepLimit(1.0*mm), fDetRot(0.)	// note: fScintStepLimit initialized here
{
  fCallFieldConstructor = false;
}


DetectorConstruction::~DetectorConstruction()
{ }

void DetectorConstruction::DefineMaterials()
{
  Vacuum = NULL;		//This value is set later using the setVacuumPressure method.
  string name,symbol;
  int z;
  G4double a;
  G4int nAtoms;
  G4double massFrac;

  new G4Element(name="H", symbol="H", z=1, a=1.0079*g/mole);
  new G4Element(name="C", symbol="C", z=6, a=12.0107*g/mole);
  new G4Element(name="N", symbol="N", z=7, a=14.0067*g/mole);
  new G4Element(name="O", symbol="O", z=8, a=15.9994*g/mole);
  new G4Element(name="Al", symbol="Al",z=13, a=26.9815*g/mole);
  new G4Element(name="Cr", symbol="Cr",z=24, a=51.9961*g/mole);
  new G4Element(name="Fe", symbol="Fe",z=26, a=55.845*g/mole);
  new G4Element(name="Ni", symbol="Ni",z=28, a=58.6934*g/mole);
  new G4Element(name="Cu", symbol="Cu",z=29, a=63.55*g/mole);
  new G4Element(name="Zn", symbol="Zn",z=30, a=65.39*g/mole);

  Be = new G4Material("Beryllium",4.,9.01*g/mole,1.848*g/cm3);
  Al = new G4Material("Aluminum",13.,26.98*g/mole,2.7*g/cm3);
  Si = new G4Material("Silicon",14.,28.09*g/mole,2.33*g/cm3);
  Cu = new G4Material("Copper", 29., 63.55*g/mole, 8.96*g/cm3);
  Wu = new G4Material("Tungsten",74.,183.84*g/mole,19.3*g/cm3);
  Au = new G4Material("Gold",79.,196.97*g/mole,19.3*g/cm3);

  Brass = new G4Material("Brass",8.5*g/cm3,2);
  Brass->AddElement(G4Element::GetElement("Cu"),massFrac=0.70);
  Brass->AddElement(G4Element::GetElement("Zn"),massFrac=0.30);

  SS304 = new G4Material("Stainless304",8.03*g/cm3,3);
  SS304->AddElement(G4Element::GetElement("Fe"),massFrac=0.70);
  SS304->AddElement(G4Element::GetElement("Cr"),massFrac=0.20);
  SS304->AddElement(G4Element::GetElement("Ni"),massFrac=0.10);

  Kevlar = new G4Material("Kevlar",1.44*g/cm3,4);
  Kevlar->AddElement(G4Element::GetElement("N"),nAtoms=2);
  Kevlar->AddElement(G4Element::GetElement("C"),nAtoms=14);
  Kevlar->AddElement(G4Element::GetElement("H"),nAtoms=10);
  Kevlar->AddElement(G4Element::GetElement("O"),nAtoms=2);

  Mylar = new G4Material("Mylar",1.4*g/cm3,3);
  Mylar->AddElement(G4Element::GetElement("C"),nAtoms=5);
  Mylar->AddElement(G4Element::GetElement("H"),nAtoms=4);
  Mylar->AddElement(G4Element::GetElement("O"),nAtoms=2);

  Polyethylene = new G4Material("Polyethylene",0.95*g/cm3,2);
  Polyethylene->AddElement(G4Element::GetElement("C"),nAtoms=2);
  Polyethylene->AddElement(G4Element::GetElement("H"),nAtoms=4);

  // Wirechamber fill: pentane @ 100torr
  double P_MWPC = 100*torr;
  double T_MWPC = 298*kelvin;
  WCPentane = new G4Material("Pentane",(72.17*mg)/(22.4*cm3)*P_MWPC/(760*torr)*(273.15*kelvin)/T_MWPC,2,kStateGas,T_MWPC,P_MWPC);
  WCPentane->AddElement(G4Element::GetElement("C"),nAtoms=5);
  WCPentane->AddElement(G4Element::GetElement("H"),nAtoms=12);

  // Wirechamber fill: N2 @ 95torr
  double P_N2 = P_MWPC - 5*torr;
  WCNitrogen = new G4Material("MWPC_N2",(28*mg)/(22.4*cm3)*P_N2/(760*torr)*(273.15*kelvin)/T_MWPC,1,kStateGas,T_MWPC,P_N2);
  WCNitrogen->AddElement(G4Element::GetElement("N"),nAtoms=2);

  // Scintillator, per Eljen EJ-204 datasheet
  Sci=new G4Material("Scintillator",1.032*g/cm3,2);
  Sci->AddElement(G4Element::GetElement("C"),nAtoms=4.68);
  Sci->AddElement(G4Element::GetElement("H"),nAtoms=5.15);
}

void DetectorConstruction::setVacuumPressure(G4double pressure)
{
  // our slightly crappy vacuum: low-pressure air (density @20c; 1.290*mg/cm3 @STP)
  G4cout<<"------------- Detector vacuum is set at "<<pressure/torr<<" Torr"<<G4endl;
  Vacuum = new G4Material("Vacuum",1.2048*mg/cm3*pressure/atmosphere,2,kStateGas,293*kelvin,pressure);
  Vacuum->AddElement(G4Element::GetElement("N"),0.78);
  Vacuum->AddElement(G4Element::GetElement("O"),0.22);
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();	// immediate call to define all materials used as class properties (so ~global access)
  setVacuumPressure(0);	// this is the set vacuum pressure that was warned about in DefineMaterials()

  // user step limits
  G4UserLimits* UserCoarseLimits = new G4UserLimits();
  UserCoarseLimits->SetMaxAllowedStep(10*m);
  G4UserLimits* UserGasLimits = new G4UserLimits();
  UserGasLimits->SetMaxAllowedStep(1*cm);
  G4UserLimits* UserSolidLimits = new G4UserLimits();
  UserSolidLimits->SetMaxAllowedStep(fScintStepLimit);	// default value from Messenger class.

  // Experimental Hall. World volume.
  G4double expHall_x = 2.0*m;
  G4double expHall_y = 2.0*m;
  G4double expHall_z = 8.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box", expHall_x/2, expHall_y/2, expHall_z/2);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box, Vacuum, "World_log");
  experimentalHall_log -> SetVisAttributes(G4VisAttributes::Invisible);
  experimentalHall_log -> SetUserLimits(UserCoarseLimits);
  experimentalHall_phys = new G4PVPlacement(NULL, G4ThreeVector(), "World_phys", experimentalHall_log, 0, false, 0);

  //----- Source holder object. Used if it is a calibration source.
  G4double source_windowThick = 4.7*um;
  G4double source_coatingThick = 0.1*um;
  G4Material* source_windowMaterial = Mylar;
  G4Material* source_coatingMaterial = Al;
  G4double source_holderThick = (3./16.)*inch;
  G4ThreeVector source_holderPos(0,0,0);
  G4double source_ringRadius = 0.5*inch;
  G4double source_windowRadius = source_ringRadius-3.0*mm;
  G4double source_ringThickness = 3.2*mm;
  G4double source_holderHeight = 1.5*inch;
  G4double source_holderWidth = 1.5*inch;

  // source holder container
  G4Box* source_holderBox = new G4Box("source_holder_box", 0.5*source_holderWidth, 0.5*source_holderHeight, 0.5*source_holderThick);
  source_container_log = new G4LogicalVolume(source_holderBox, Vacuum, "source_container_log");

  // source holder paddle
  G4Tubs* source_holderHole = new G4Tubs("source_holder_hole", 0., source_ringRadius, source_holderThick, 0., 2*M_PI);
  G4SubtractionSolid* source_holder = new G4SubtractionSolid("source holder", source_holderBox, source_holderHole);
  G4LogicalVolume* source_holder_log = new G4LogicalVolume(source_holder, Brass, "source_holder_log");
  source_holder_log -> SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.7,0,0.5)));
  source_holder_phys = new G4PVPlacement(NULL, G4ThreeVector(), source_holder_log, "source_holder_phys", source_container_log, false, 0);

  // sealed source foil
  G4Tubs* source_windowTube = new G4Tubs("window_tube", 0., source_windowRadius, source_windowThick, 0., 2*M_PI);
  source_window_log = new G4LogicalVolume(source_windowTube, source_windowMaterial, "source_window_log");
  G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
  source_window_log->SetVisAttributes(visWindow);
  source_window_phys = new G4PVPlacement(NULL, G4ThreeVector(), source_window_log, "source_window_phys", source_container_log, false, 0);

  // source foil coating
  G4Tubs* source_coating_tube = new G4Tubs("source_coating_tube", 0., source_windowRadius, source_coatingThick*0.5, 0., 2*M_PI);
  for(int i = 0; i <= 1; i++)	// 0 = EAST, 1 = WEST
  {
    source_coating_log[i] = new G4LogicalVolume(source_coating_tube, source_coatingMaterial, Append(i, "source_coating_log"));
    source_coating_log[i] -> SetVisAttributes(new G4VisAttributes(G4Colour(0,1,0,0.5)));
  }

  source_coating_phys[0] = new G4PVPlacement(NULL, G4ThreeVector(0,0, (-1)*(source_windowThick + source_coatingThick*0.5)),
					source_coating_log[0], "source_coating_phys_0", source_container_log, false, 0);
  source_coating_phys[1] = new G4PVPlacement(NULL, G4ThreeVector(0,0, source_windowThick + source_coatingThick*0.5),
					source_coating_log[1], "source_coating_phys_1", source_container_log, false, 0);

  // source retaining ring
  G4Tubs* source_ringTube = new G4Tubs("source_ring_tube", source_windowRadius, source_ringRadius, source_ringThickness/2., 0., 2*M_PI);
  G4LogicalVolume* source_ring_log = new G4LogicalVolume(source_ringTube, Al, "source_ring_log");
  source_ring_log -> SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.5)));
  source_ring_phys = new G4PVPlacement(NULL, G4ThreeVector(), source_ring_log, "source_ring_phys", source_container_log, false, 0);

  // place entire source holder object
  source_phys = new G4PVPlacement(NULL, source_holderPos, source_container_log,"source_container_phys", experimentalHall_log, false, 0, true);


  //----- Decay Trap object (length 3m, main tube)
  G4double decayTrap_windowThick = 0.180*um;
  G4double decayTrap_coatingThick = 0.150*um;
  G4double decayTrap_innerRadiusOfTrap = 2.45*inch;
  G4double decayTrap_tubeWallThick = 2*mm;
  G4double decayTrap_innerRadiusCollimator = 2.3*inch;
  G4Material* decayTrap_tubeMaterial = Cu;
  G4Material* decayTrap_collimatorMaterial = Polyethylene;
  G4Material* decayTrap_windowMaterial = Mylar;
  G4Material* decayTrap_coatingMaterial = Be;

  // decay tube construction
  G4double decayTrap_tube_outerRadius = decayTrap_innerRadiusOfTrap + decayTrap_tubeWallThick;
  G4double decayTrap_tube_length = 3.0*m;

  G4Tubs* decayTrap_tube = new G4Tubs("decayTrap_tube", decayTrap_innerRadiusOfTrap, decayTrap_tube_outerRadius,
					decayTrap_tube_length/2., 0., 2*M_PI);
  decayTrap_tube_log = new G4LogicalVolume(decayTrap_tube, decayTrap_tubeMaterial, "decayTrap_tube_log");
  decayTrap_tube_log -> SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
  new G4PVPlacement(NULL, G4ThreeVector(), decayTrap_tube_log, "decayTrap_tube", experimentalHall_log, false, 0, true);

  // decay trap windows, collimator, monitors
  G4double decayTrap_totalWindowThickness = decayTrap_windowThick + decayTrap_coatingThick;
  G4double decayTrap_collimatorThick = 0.8*inch;
  G4double decayTrap_beWindow_PosZ = -decayTrap_totalWindowThickness/2. + decayTrap_coatingThick/2.;
  G4double decayTrap_mylarWindow_PosZ = decayTrap_totalWindowThickness/2. - decayTrap_windowThick/2.;
  G4double decayTrap_window_PosZ = (decayTrap_tube_length + decayTrap_totalWindowThickness)/2.;
  G4double decayTrap_monitorThickness = 1.0*mm;
  G4double decayTrap_monitor_PosZ = 0.5*m;

  G4Tubs* decayTrap_trapWindowTube = new G4Tubs("trap_win_tube", 0., decayTrap_tube_outerRadius, decayTrap_totalWindowThickness/2., 0, 2*M_PI);
  G4Tubs* decayTrap_mylarTube = new G4Tubs("mylarTube", 0., decayTrap_tube_outerRadius, decayTrap_windowThick/2., 0., 2*M_PI);
  G4Tubs* decayTrap_beTube = new G4Tubs("beTube", 0., decayTrap_tube_outerRadius, decayTrap_coatingThick/2., 0., 2*M_PI);
  G4Tubs* decayTrap_collimatorTube = new G4Tubs("decayTrap_collimatorTube", decayTrap_innerRadiusCollimator,
				decayTrap_innerRadiusCollimator + decayTrap_collimatorThick, decayTrap_collimatorThick/2.,0., 2*M_PI);
  G4Tubs* decayTrap_collimatorBackTube = new G4Tubs("decayTrap_collimatorBackTube", decayTrap_tube_outerRadius + 1.*mm,
				decayTrap_innerRadiusCollimator + decayTrap_collimatorThick, decayTrap_collimatorThick/2., 0., 2*M_PI);
  G4Tubs* decayTrap_monitorTube = new G4Tubs("trap_monitor_tube", 0., decayTrap_innerRadiusOfTrap, decayTrap_monitorThickness/2.0,
					0., 2*M_PI);

  for(int i = 0; i <= 1; i++)
  {
    decayTrap_window_log[i] = new G4LogicalVolume(decayTrap_trapWindowTube, Vacuum, Append(i, "trap_win_log_"));
    decayTrap_window_log[i] -> SetVisAttributes(visWindow);
    decayTrap_mylarWindow_log[i] = new G4LogicalVolume(decayTrap_mylarTube, decayTrap_windowMaterial, Append(i, "mylar_win_log_"));
    decayTrap_beWindow_log[i] = new G4LogicalVolume(decayTrap_beTube, decayTrap_coatingMaterial, Append(i, "be_win_log"));
  }
  new G4PVPlacement(NULL, G4ThreeVector(0.,0.,(-1)*decayTrap_window_PosZ), decayTrap_window_log[0], "trap_win_0",
			experimentalHall_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, decayTrap_window_PosZ), decayTrap_window_log[1], "trap_win_1",
			experimentalHall_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., (-1)*decayTrap_mylarWindow_PosZ), decayTrap_mylarWindow_log[0],
			"mylar_win_0", decayTrap_window_log[0], false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., (-1)*decayTrap_beWindow_PosZ), decayTrap_beWindow_log[0],
			"be_win_0", decayTrap_window_log[0], false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., decayTrap_mylarWindow_PosZ), decayTrap_mylarWindow_log[1],
                        "mylar_win_1", decayTrap_window_log[1], false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., decayTrap_beWindow_PosZ), decayTrap_beWindow_log[1],
                        "be_win_1", decayTrap_window_log[1], false, 0);

  G4double decayTrap_collimator_PosZ = (decayTrap_tube_length + decayTrap_collimatorThick)/2.;
  decayTrap_collimator_PosZ += (decayTrap_totalWindowThickness)/2.;
  G4double decayTrap_collimatorBack_PosZ = decayTrap_tube_length/2. - decayTrap_collimatorThick;
  for(int i = 0; i <= 1; i++)
  {
    decayTrap_collimator_log[i] = new G4LogicalVolume(decayTrap_collimatorTube, decayTrap_collimatorMaterial,
							Append(i, "collimator_log_"));
    decayTrap_collimatorBack_log[i] = new G4LogicalVolume(decayTrap_collimatorBackTube, decayTrap_collimatorMaterial,
							Append(i, "collimator_back_log_"));
    decayTrap_innerMonitors_log[i] = new G4LogicalVolume(decayTrap_monitorTube, Vacuum, Append(i, "trap_monitor_log_"));
  }
  // place everything at -z i.e. EAST.
  new G4PVPlacement(NULL, G4ThreeVector(0, 0, (-1)*decayTrap_collimator_PosZ), decayTrap_collimator_log[0],
			"collimator_0", experimentalHall_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, (-1)*decayTrap_collimatorBack_PosZ), decayTrap_collimatorBack_log[0],
			"collimator_back_0", experimentalHall_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, (-1)*decayTrap_monitor_PosZ), decayTrap_innerMonitors_log[0],
			"trap_monitor_0", experimentalHall_log, false, 0);
  // copy but place at +z i.e. WEST
  new G4PVPlacement(NULL, G4ThreeVector(0, 0, decayTrap_collimator_PosZ), decayTrap_collimator_log[1],
                        "collimator_1", experimentalHall_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, decayTrap_collimatorBack_PosZ), decayTrap_collimatorBack_log[1],
                        "collimator_back_1", experimentalHall_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, decayTrap_monitor_PosZ), decayTrap_innerMonitors_log[1],
                        "trap_monitor_1", experimentalHall_log, false, 0);

  decayTrap_window_log[0] -> SetUserLimits(UserSolidLimits);
  decayTrap_window_log[1] -> SetUserLimits(UserSolidLimits);

  G4cout << "PAST DECAY TRAP CODE" << G4endl;

  //----- Scintillator construction. Used as Sensitive Volume
  G4double scint_scintRadius = 7.5*cm;
  G4double scint_scintBackingRadius = 10*cm;
  G4double scint_scintThick = 3.5*mm;
  G4double scint_deadLayerThick = 3.0*um;
  G4double scint_scintBackingThick = 1.*inch;
  G4double scint_lightGuideThick = 1.0*cm;
  G4double scint_N2Volume_Z = scint_lightGuideThick + scint_scintBackingThick;
  G4double scint_face_PosZ = -scint_N2Volume_Z/2.;

  if((scint_scintBackingRadius < scint_scintRadius) || (scint_lightGuideThick < scint_scintThick))
	G4cout << "\n\nMajor geometry error! Scintillator measurements don't make sense! \n \n" << G4endl;

//  ConstructField(); // since you need access to fpMagField

  // Here, in your old code, you'd loop over sd to create EAST and WEST objects.
  // Now, you'll build the East first (for debugging) scint, mwpc and frame.

  // Overall container layer for the scintillator
  G4Tubs* scint_N2VolTube = new G4Tubs("N2_vol_tube", 0., scint_scintBackingRadius, scint_N2Volume_Z/2., 0., 2*M_PI);
  scint_container_log = new G4LogicalVolume(scint_N2VolTube, WCNitrogen, "N2_Vol_log_EAST");
  scint_container_log -> SetVisAttributes(G4VisAttributes::Invisible);

  // dead layer in scint
  G4Tubs* scint_deadLayerTube = new G4Tubs("Dead_scint_tube", 0, scint_scintRadius, scint_deadLayerThick/2., 0., 2*M_PI);
  scint_deadLayer_log = new G4LogicalVolume(scint_deadLayerTube, Sci, "Dead_scint_log_EAST");
  G4VisAttributes* visDScint= new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));
  scint_deadLayer_log -> SetVisAttributes(visDScint);
  scint_deadLayer_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -(scint_N2Volume_Z - scint_deadLayerThick)/2.),
				scint_deadLayer_log, "Dead_scint_phys_EAST", scint_container_log, false, 0);

  // scintillator
  G4Tubs* scint_scintTube = new G4Tubs("scint_tube", 0, scint_scintRadius, (scint_scintThick - scint_deadLayerThick)/2., 0., 2*M_PI);
  scint_scintillator_log = new G4LogicalVolume(scint_scintTube, Sci, "scint_log_EAST");
  G4VisAttributes* visScint= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.2));
  scint_scintillator_log -> SetVisAttributes(visScint);
  scint_scintillator_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -scint_N2Volume_Z/2. + scint_deadLayerThick + (scint_scintThick - scint_deadLayerThick)/2.),
				scint_scintillator_log, "scint_phys_EAST", scint_container_log, false, 0);

  // light guides around and behind detector
  G4double scint_zPlane[] = {0., scint_scintThick, scint_scintThick, scint_lightGuideThick};
  G4double scint_lightGuideRadius = scint_scintRadius - (scint_lightGuideThick - scint_scintThick);
  G4double scint_innerRad[] = {scint_scintRadius, scint_scintRadius, scint_lightGuideRadius, scint_lightGuideRadius};
  G4double scint_outerRad[] = {scint_scintBackingRadius, scint_scintBackingRadius, scint_scintBackingRadius, scint_scintBackingRadius};

  G4Polycone* scint_lightGuidePoly = new G4Polycone("lightguide_polycone", 0., 2*M_PI, 4, scint_zPlane, scint_innerRad, scint_outerRad);
  scint_lightGuide_log = new G4LogicalVolume(scint_lightGuidePoly, Sci, "light_guide_log_EAST");
  G4VisAttributes* visLG = new G4VisAttributes(G4Colour(0.0,1.0,0.5,0.2));
  scint_lightGuide_log -> SetVisAttributes(visLG);
  scint_lightGuide_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, scint_N2Volume_Z/2.), scint_lightGuide_log, "light_guide_phys_EAST",
				scint_container_log, false, 0);

  // backing veto
  G4Tubs* scint_backingTube = new G4Tubs("backing_tube", 0., scint_scintBackingRadius, scint_scintBackingThick/2., 0., 2*M_PI);
  scint_backing_log = new G4LogicalVolume(scint_backingTube, Sci, "backing_log_EAST");
  G4VisAttributes* visBacking= new G4VisAttributes(G4Colour(0.0,0.0,1,0.2));
  scint_backing_log->SetVisAttributes(visBacking);
  scint_backing_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, (scint_N2Volume_Z - scint_scintBackingThick)/2.),
				scint_backing_log, "backing_phys_EAST", scint_container_log, false, 0);

  G4ThreeVector sideTransScint = G4ThreeVector(0., 0., (-1)*(2.2*m - scint_face_PosZ));
  G4RotationMatrix* sideRot = new G4RotationMatrix();
  sideRot -> rotateY(M_PI*rad);

  // This places correctly on the EAST side!
  scint_container_phys = new G4PVPlacement(sideRot, sideTransScint, scint_container_log, "N2_vol_phys_EAST", experimentalHall_log, false, 0, true);

  G4cout << "PAST THE SCINTILLATOR CODE" << G4endl;

  //----- Begin Wire volume construction. Active region inside wire chamber.
  G4double wireVol_anodeRadius = 5*um;
  G4double wireVol_cathodeRadius = 25*um;
  G4double wireVol_platingThick = 0.2*um;
  G4double wireVol_wireSpacing = 2.54*mm;
  G4double wireVol_NbOfWires = 64;
  G4double wireVol_planeSpacing = 1*cm;

  G4Material* wireVol_cathodeWireMat = Al;
  G4Material* wireVol_anodeWireMat = Wu;
  G4Material* wireVol_cathodePlateMat = Au;
  G4Material* wireVol_activeGas = WCPentane;

  G4double wireVol_wirePlaneWidth = wireVol_NbOfWires*wireVol_wireSpacing;

  // effective mwpc gas volume containing the cathodes and anodes
  G4Box* wireVol_mwpcGasBox = new G4Box("mpwc_gas_box", wireVol_wirePlaneWidth/2., wireVol_wirePlaneWidth/2., wireVol_planeSpacing);
  wireVol_gas_log = new G4LogicalVolume(wireVol_mwpcGasBox, wireVol_activeGas, "mwpc_gas_log_EAST");
  wireVol_gas_log -> SetVisAttributes(G4VisAttributes::Invisible);

  // anode, cathode wire containers. These are "on their side" to allow wireplane parametrization. Rotated later.
  G4Box* wireVol_cathContainerBox = new G4Box("cathContainer_Box", wireVol_wirePlaneWidth/2., wireVol_cathodeRadius, wireVol_wirePlaneWidth/2.);
  G4Box* wireVol_anodeContainerBox = new G4Box("anodeContainer_Box", wireVol_wirePlaneWidth/2., wireVol_anodeRadius, wireVol_wirePlaneWidth/2.);

  // anode, cathode wires and surrouding gas
  G4Tubs* wireVol_cathPlateTube = new G4Tubs("cathplate_tube", wireVol_cathodeRadius - wireVol_platingThick,
						wireVol_cathodeRadius, wireVol_wirePlaneWidth/2., 0., 2*M_PI);
  G4Tubs* wireVol_cathodeTube = new G4Tubs("cathode_tube", 0, wireVol_cathodeRadius- wireVol_platingThick,
						wireVol_wirePlaneWidth/2., 0., 2*M_PI);
  G4Tubs* wireVol_anodeTube = new G4Tubs("anode_tube", 0, wireVol_anodeRadius, wireVol_wirePlaneWidth/2., 0, 2*M_PI);
  G4Box* wireVol_cathSegBox = new G4Box("cathodeSegmentBox", wireVol_wireSpacing/2., wireVol_cathodeRadius, wireVol_wirePlaneWidth/2.);
  G4Box* wireVol_anodeSegBox = new G4Box("anodeSegmentBox", wireVol_wireSpacing/2., wireVol_anodeRadius, wireVol_wirePlaneWidth/2.);

  // Rotate 90 degrees around X axis
  G4RotationMatrix* xRot90 = new G4RotationMatrix;
  xRot90->rotateX(M_PI/2.*rad);
  // Rotate 90 degrees around X then Z axis (Y axis in object coordinates)
  G4RotationMatrix* xzRot90 = new G4RotationMatrix;
  xzRot90->rotateX(M_PI/2.*rad);
  xzRot90->rotateY(M_PI/2.*rad);

  G4VisAttributes* visCathWires = new G4VisAttributes(G4Colour(1,0.7,0,0.8));
  G4VisAttributes* visAnodeWires = new G4VisAttributes(G4Colour(1,0.3,0,0.8));

  // make anode, cathode segments
  wireVol_cathSeg_log = new G4LogicalVolume(wireVol_cathSegBox, wireVol_activeGas, "cathSeg_log_EAST");
  wireVol_anodeSeg_log = new G4LogicalVolume(wireVol_anodeSegBox, wireVol_activeGas, "anodeSeg_log_EAST");
  wireVol_cathodeWire_log = new G4LogicalVolume(wireVol_cathodeTube, wireVol_cathodeWireMat, "cathode_log_EAST");
  wireVol_cathPlate_log = new G4LogicalVolume(wireVol_cathPlateTube, wireVol_cathodePlateMat, "cathode_plate_log_EAST");
  wireVol_anodeWire_log = new G4LogicalVolume(wireVol_anodeTube, wireVol_anodeWireMat, "anode_log_EAST");
  wireVol_cathSeg_log -> SetVisAttributes(G4VisAttributes::Invisible);
  wireVol_anodeSeg_log -> SetVisAttributes(G4VisAttributes::Invisible);
  wireVol_cathodeWire_log -> SetVisAttributes(visCathWires);
  wireVol_cathPlate_log -> SetVisAttributes(visCathWires);
  wireVol_anodeWire_log -> SetVisAttributes(visAnodeWires);

  new G4PVPlacement(NULL, G4ThreeVector(), wireVol_cathodeWire_log, "cathode_wire_phys_EAST", wireVol_cathSeg_log, true, 0);
  new G4PVPlacement(NULL, G4ThreeVector(), wireVol_cathPlate_log, "cathode_plate_phys_EAST", wireVol_cathSeg_log, true, 0);
  new G4PVPlacement(NULL, G4ThreeVector(), wireVol_anodeWire_log, "anode_wire_phys_EAST", wireVol_anodeSeg_log, true, 0);

  // make and place anode and cathode plane container volumes
  G4LogicalVolume* wireVol_cathContainer1_log = new G4LogicalVolume(wireVol_cathContainerBox, wireVol_activeGas, "cathContainer1_log_EAST");
  G4LogicalVolume* wireVol_cathContainer2_log = new G4LogicalVolume(wireVol_cathContainerBox, wireVol_activeGas, "cathContainer2_log_EAST");
  G4LogicalVolume* wireVol_anodeContainer_log = new G4LogicalVolume(wireVol_anodeContainerBox, wireVol_activeGas, "anodeContainer_log_EAST");

  new G4PVPlacement(xRot90, G4ThreeVector(0., 0., wireVol_cathodeRadius - wireVol_planeSpacing),
	wireVol_cathContainer1_log, "cathContainer1_phys_EAST", wireVol_gas_log, false, 0);
  new G4PVPlacement(xzRot90, G4ThreeVector(0., 0., (-1)*(wireVol_cathodeRadius - wireVol_planeSpacing)),
	wireVol_cathContainer2_log, "cathContainer2_phys_EAST", wireVol_gas_log, false, 0);
  new G4PVPlacement(xRot90, G4ThreeVector(0,0,0), wireVol_anodeContainer_log, "anodeContainer_phys_EAST", wireVol_gas_log, false,0);

  // replicate the segments defined above into cathode, anode arrays
  new G4PVReplica("CathodeArray1_EAST", wireVol_cathSeg_log, wireVol_cathContainer1_log, kXAxis, wireVol_NbOfWires, wireVol_wireSpacing);
  new G4PVReplica("CathodeArray2_EAST", wireVol_cathSeg_log, wireVol_cathContainer2_log, kXAxis, wireVol_NbOfWires, wireVol_wireSpacing);
  new G4PVReplica("AnodeArray_EAST", wireVol_anodeSeg_log, wireVol_anodeContainer_log, kXAxis, wireVol_NbOfWires, wireVol_wireSpacing);

  G4cout << "PAST THE WIRE CHAMBER ACTIVE VOLUME CODE" << G4endl;

  //----- Begin wirechamber construction. MWPC used in front of Scintillator.
  G4double mwpc_windowThick = 6*um;
  G4double mwpc_entranceRadius = 7.0*cm;
  G4double mwpc_exitRadius = 7.5*cm;
  G4double mwpc_entranceToCathodes = 5.0*mm;
  G4double mwpc_exitToCathodes = 5.0*mm;
  G4MagneticField* mwpc_internalBfield = NULL;	// this is tricky. Need to figure out.
  G4double mwpc_fieldE0 = 0;	// initialize field potential to 0. Also tricky.
  G4Material* mwpc_fillGas = wireVol_activeGas;	// want it to be WCPentane

  G4double mwpc_containerHalf_Z = 0.5*(mwpc_entranceToCathodes + mwpc_exitToCathodes + 2*cm);
  G4double mwpc_gasVolumeWidth = 8.0*inch;	// MWPC gas box width

  // container volume for all MWPC
  G4Box* mwpc_containerBox = new G4Box("mwpc_container_box", mwpc_gasVolumeWidth/2., mwpc_gasVolumeWidth/2., mwpc_containerHalf_Z);
  mwpc_container_log = new G4LogicalVolume(mwpc_containerBox, mwpc_fillGas, "mwpc_container_log_EAST");
  mwpc_container_log -> SetVisAttributes(G4VisAttributes::Invisible);

  // MWPC active gas volume placement with  wireplane, relative to MWPC container volume
  G4ThreeVector mwpc_activeRegionTrans(0, 0, (mwpc_entranceToCathodes - mwpc_exitToCathodes)/2.);
  // Note: this line places wireVol inside mwpc.
  new G4PVPlacement(NULL, mwpc_activeRegionTrans, wireVol_gas_log, "mwpc_phys_EAST", mwpc_container_log, false, 0);

  // construct kevlay string. Rectangular cross section strings with equal volume to nominal 140um cylinders.
  G4double mwpc_kevRadius = 0.07*mm;
  G4double mwpc_kevSpacing = 5.0*mm;
  G4int mwpc_NbKevWires = 32;
  G4double mwpc_kevLength = 15.0*cm;
  double mwpc_kevAspectRatio = 16.0;	// aspect ratio, width:depth.
  G4double mwpc_kevArea = M_PI*mwpc_kevRadius*mwpc_kevRadius;
  G4double mwpc_kevEffWidth = sqrt(mwpc_kevArea*mwpc_kevAspectRatio);
  G4double mwpc_kevEffThick = sqrt(mwpc_kevArea/mwpc_kevAspectRatio);
  G4double mwpc_kev_PosZ = -mwpc_containerHalf_Z + mwpc_kevEffThick/2.;

  G4Box* mwpc_kevContainerBox = new G4Box("kevContainer_box", mwpc_NbKevWires*mwpc_kevSpacing/2., mwpc_kevLength/2., mwpc_kevEffThick/2.);
  G4Box* mwpc_kevSegBox = new G4Box("kevSeg_box", mwpc_kevSpacing/2., mwpc_kevLength/2., mwpc_kevEffThick/2.);
  G4Box* mwpc_kevStripBox = new G4Box("kevStrip_box", mwpc_kevEffWidth/2., mwpc_kevLength/2., mwpc_kevEffThick/2.);
  mwpc_kevContainer_log = new G4LogicalVolume(mwpc_kevContainerBox, Vacuum, "kevContainer_log_EAST");
  mwpc_kevSeg_log = new G4LogicalVolume(mwpc_kevSegBox, Vacuum, "kevSeg_log_EAST");
  mwpc_kevStrip_log = new G4LogicalVolume(mwpc_kevStripBox, Kevlar, "kevStrip_log_EAST");

  // place components of kevlar strings. Replicate the array
  new G4PVPlacement(NULL, G4ThreeVector(0,0, mwpc_kev_PosZ), mwpc_kevContainer_log, "kevContainer_phys_EAST", mwpc_container_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0,0), mwpc_kevStrip_log, "kevStrip_phys_EAST", mwpc_kevSeg_log, false, 0);
  new G4PVReplica("kevlar_place_EAST", mwpc_kevSeg_log, mwpc_kevContainer_log, kXAxis, mwpc_NbKevWires, mwpc_kevSpacing);

  // Make the Mylar windows in the MWPC.
  G4Tubs* mwpc_winInnerTube = new G4Tubs("winInnerTube", 0., mwpc_entranceRadius, mwpc_windowThick/2., 0., 2*M_PI);
  G4Tubs* mwpc_winOuterTube = new G4Tubs("winOuterTube", 0., mwpc_exitRadius, mwpc_windowThick/2., 0., 2*M_PI);
  mwpc_winIn_log = new G4LogicalVolume(mwpc_winInnerTube, Mylar, "winIn_log_EAST");
  mwpc_winIn_log -> SetVisAttributes(visWindow);
  mwpc_winOut_log = new G4LogicalVolume(mwpc_winOuterTube, Mylar, "winOut_log_EAST");
  mwpc_winOut_log -> SetVisAttributes(visWindow);

  new G4PVPlacement(NULL, G4ThreeVector(0,0, -mwpc_containerHalf_Z + mwpc_kevEffThick + mwpc_windowThick/2.),
	mwpc_winIn_log, "winIn_phys_EAST", mwpc_container_log, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, mwpc_containerHalf_Z - mwpc_windowThick/2.),
	mwpc_winOut_log," winOut_phys_EAST", mwpc_container_log, false, 0);

  G4double frame_backWinFrameThick = 0.5*inch;	// originally placed further down but needed here
  G4double mwpc_PosZ = -mwpc_containerHalf_Z - frame_backWinFrameThick - (scint_N2Volume_Z/2. + scint_face_PosZ);
  G4ThreeVector sideTransMWPC = G4ThreeVector(0,0, (-1)*(2.2*m + mwpc_PosZ));

  mwpc_container_phys = new G4PVPlacement(sideRot, sideTransMWPC, mwpc_container_log, "mwpc_container_phys_EAST", experimentalHall_log, false, 0, true);


  G4cout << "PAST THE WIRECHAMBER CONSTRUCTION CODE." << G4endl;



  //----- Begin DetectorPackageConstruction. This is the frame that holds the scintillator and MWPC.

  G4double frame_packageRadius = 6.0*inch;
  G4double frame_mwpcEntranceThick = 0.375*inch;
  G4double frame_mwpcEntranceRadius = 3.0*inch;	// not the same value as mwpc_entranceRadius
  G4double frame_mwpcEntranceDepth = 5.0*inch;
  G4double frame_frontWinFrameThick = 1.0*inch;
//  G4double frame_backWinFrameThick = 0.5*inch;	// moved a few lines up since needed for mwpc placement

  // aluminum entrance collimator to detector package
  G4double frame_entranceSectionLength = frame_mwpcEntranceDepth + frame_frontWinFrameThick;
  G4Tubs* frame_mwpcEntranceTube = new G4Tubs("mwpc_entrance_tube", 0., frame_packageRadius, 0.5*frame_entranceSectionLength, 0., 2*M_PI);
  G4Tubs* frame_entranceFrontTube = new G4Tubs("entrance_front_tube", frame_mwpcEntranceRadius + frame_mwpcEntranceThick,
					frame_packageRadius, 0.5*frame_mwpcEntranceThick, 0., 2*M_PI);
  G4Tubs* frame_entranceMidTube = new G4Tubs("entrance_mid_tube", frame_mwpcEntranceRadius, frame_mwpcEntranceRadius + frame_mwpcEntranceThick,
					0.5*frame_mwpcEntranceDepth, 0., 2*M_PI);
  G4Tubs* frame_entranceBackTube = new G4Tubs("entrance_back_tube", mwpc_entranceRadius, frame_packageRadius,
					0.5*frame_frontWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCEntrance = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.8));

  frame_mwpcEntrance_log = new G4LogicalVolume(frame_mwpcEntranceTube, Vacuum, "mwpc_entrance_log_EAST");
  frame_entranceFront_log = new G4LogicalVolume(frame_entranceFrontTube, Al, "entrance_front_log_EAST");
  frame_entranceMid_log = new G4LogicalVolume(frame_entranceMidTube, Al, "entrance_mid_log_EAST");
  frame_entranceBack_log = new G4LogicalVolume(frame_entranceBackTube, Al, "entrance_back_log_EAST");
  frame_mwpcEntrance_log -> SetVisAttributes(G4VisAttributes::Invisible);
  frame_entranceFront_log -> SetVisAttributes(visMWPCEntrance);
  frame_entranceMid_log -> SetVisAttributes(visMWPCEntrance);
  frame_entranceBack_log -> SetVisAttributes(visMWPCEntrance);

//  frame_entranceFront_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*(frame_entranceSectionLength - frame_mwpcEntranceThick)),
//				frame_entranceFront_log, "entrance_front_phys_EAST", frame_mwpcEntrance_log, false, 0);
//  frame_entranceMid_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*frame_frontWinFrameThick),
//				frame_entranceMid_log, "entrance_mid_phys_EAST", frame_mwpcEntrance_log, false, 0);
//  frame_entranceBack_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, 0.5*(frame_entranceSectionLength - frame_frontWinFrameThick)),
//				frame_entranceBack_log, "entrance_back_phys_EAST", frame_mwpcEntrance_log, false, 0);

  // create overall detector package frame.
  G4double frame_detFrameHalf_Z = frame_mwpcEntranceDepth + 2*mwpc_containerHalf_Z + 1.0*inch;
  G4Tubs* frame_framePackageTube = new G4Tubs("detPackage_tube_EAST", 0, frame_packageRadius, frame_detFrameHalf_Z, 0, 2*M_PI);
  frame_container_log = new G4LogicalVolume(frame_framePackageTube, Vacuum, "frame_container_log_EAST");
  frame_container_log -> SetVisAttributes(G4VisAttributes::Invisible);

  // HERE BE CODE THAT PLACES THE COMPONENTS i.e. scint and mwpc INSIDE THE DETECTOR PACKAGE.
  // Ignore this shit and place it yourself in experimental hall. Use these measurements to find centre!

  // aluminum exit window and N2 volume at back of gas box
  G4Tubs* frame_mwpcExitTube = new G4Tubs("mwpc_exit_tube", mwpc_exitRadius, frame_packageRadius, 0.5*frame_backWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCExit = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.8));
  frame_mwpcExit_log = new G4LogicalVolume(frame_mwpcExitTube, Al, "mwpc_exit_log_EAST");
  frame_mwpcExit_log -> SetVisAttributes(visMWPCExit);

  //****** You need to fix this z position! It is defined relative to the internal placement of mwpc in detector package. *****//
  G4double frame_exitWin_PosZ = 0;	// can't fill this variable w/o dealing with relative placement of mwpc
					// ignore for now. Will do TOMORROW.
//  frame_mwpcExit_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, frame_exitWin_PosZ), frame_mwpcExit_log,
//				"mwpc_exit_EAST", frame_container_log, false, 0);

  G4Tubs* frame_mwpcExitGasN2Tube = new G4Tubs("mwpc_exit_N2_tube", 0, mwpc_exitRadius, 0.5*frame_backWinFrameThick, 0., 2*M_PI);
  frame_mwpcExitGasN2_log = new G4LogicalVolume(frame_mwpcExitGasN2Tube, WCNitrogen, "mwpc_exit_N2_log_EAST");
  frame_mwpcExitGasN2_log -> SetVisAttributes(G4VisAttributes::Invisible);
//  frame_mwpcExitGasN2_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, frame_exitWin_PosZ), frame_mwpcExitGasN2_log,
//				"mwpc_exit_N2_phys_EAST", frame_container_log, false, 0);

  // material behind the detector. Misc stuff that can cause back scattering events.
  G4double frame_backStuffThick = 1.0*inch;
  G4Tubs* frame_backStuffTube = new G4Tubs("backstuff_tube_EAST", 0, 0.5*frame_packageRadius, frame_backStuffThick, 0., 2*M_PI);
  frame_backStuff_log = new G4LogicalVolume(frame_backStuffTube, SS304, "backStuff_log_EAST");
//  frame_backStuff_phys = new G4PVPlacement(NULL, G4ThreeVector(0, 0, frame_detFrameHalf_Z - 0.5*frame_backStuffThick),
//				frame_backStuff_log, "backStuff_phys_EAST", frame_container_log, false, 0);

                                                                        
  for(int sd = 0; sd <=1; sd++)
  {

  G4double fN2_volume_Z = scint_N2Volume_Z;
  G4double fScintFacePos = scint_face_PosZ;
//  G4Material* fMWPCGas = WCPentane;
  G4double fmwpc_entrance_R = mwpc_entranceRadius;
  G4double fmwpc_exit_R = mwpc_exitRadius;
  G4double fmwpcContainer_halfZ = mwpc_containerHalf_Z;

  G4cout << "Begin detector package construction" << G4endl;
  //----- Begin DetectorPackageConstruction class code -----//
  fDPC_detPackageRadius = 6.0*inch;		// initialize the DetectorPackageConstruction variables
  fDPC_mwpc_entrance_thickness = 0.375*inch;
  fDPC_mwpc_entrance_r = 3.0*inch;
  fDPC_mwpc_entrance_depth = 5.0*inch;
  fDPC_frontwin_frame_thick = 1.0*inch;
  fDPC_backwin_frame_thick = 0.5*inch;
  //Note: here in Michael's code an instance of mwpc and scint are declared given an input Side sd
  //These are instances of WirechamberConstruction and ScintillatorConstruction.

  // aluminum entrance collimator to detector package
  const G4double entrance_section_length = fDPC_mwpc_entrance_depth+fDPC_frontwin_frame_thick;
  G4Tubs* mwpc_entrance_tube = new G4Tubs("mwpc_entrance_tube",0,fDPC_detPackageRadius,0.5*entrance_section_length,0.,2*M_PI);
  G4Tubs* entrance_front_tube = new G4Tubs("entrance_front_tube",fDPC_mwpc_entrance_r+fDPC_mwpc_entrance_thickness,
						fDPC_detPackageRadius,0.5*fDPC_mwpc_entrance_thickness,0.,2*M_PI);
  G4Tubs* entrance_mid_tube = new G4Tubs("entrance_mid_tube",fDPC_mwpc_entrance_r,
						fDPC_mwpc_entrance_r+fDPC_mwpc_entrance_thickness,0.5*fDPC_mwpc_entrance_depth,0.,2*M_PI);
  G4Tubs* entrance_back_tube = new G4Tubs("entrance_back_tube",fmwpc_entrance_R,
						fDPC_detPackageRadius,0.5*fDPC_frontwin_frame_thick,0.,2*M_PI);

  fDPC_mwpc_entrance_log = new G4LogicalVolume(mwpc_entrance_tube, Vacuum, Append(sd,"mwpc_entrance_log_"));
  fDPC_mwpc_entrance_log->SetVisAttributes(G4VisAttributes::Invisible);
  fDPC_entrance_front_log = new G4LogicalVolume(entrance_front_tube, Al, Append(sd,"entrance_front_log_"));
  fDPC_entrance_mid_log = new G4LogicalVolume(entrance_mid_tube, Al, Append(sd,"entrance_mid_log_"));
  fDPC_entrance_back_log = new G4LogicalVolume(entrance_back_tube, Al, Append(sd,"entrance_back_log_"));
  fDPC_entrance_front_log->SetVisAttributes(visMWPCEntrance);
  fDPC_entrance_mid_log->SetVisAttributes(visMWPCEntrance);
  fDPC_entrance_back_log->SetVisAttributes(visMWPCEntrance);

//  fDPC_entrance_front_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-0.5*(entrance_section_length-fDPC_mwpc_entrance_thickness)),
//						fDPC_entrance_front_log,Append(sd,"entrance_front_"),fDPC_mwpc_entrance_log,false,0);
//  fDPC_entrance_mid_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-0.5*fDPC_frontwin_frame_thick),
//						fDPC_entrance_mid_log,Append(sd,"entrance_mid_"),fDPC_mwpc_entrance_log,false,0);
//  fDPC_entrance_back_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,0.5*(entrance_section_length-fDPC_frontwin_frame_thick)),
//						fDPC_entrance_back_log,Append(sd,"entrance_back_"),fDPC_mwpc_entrance_log,false,0);

  // overall detector package
  const G4double detPackageHalfZ = fDPC_mwpc_entrance_depth+(2*fmwpcContainer_halfZ)+1.0*inch;
  G4Tubs* detPackage_tube = new G4Tubs(Append(sd,"detPackage_tube_"),0,fDPC_detPackageRadius,detPackageHalfZ,0.,2*M_PI);
  fDPC_container_log = new G4LogicalVolume(detPackage_tube,Vacuum,Append(sd,"container_log_"));
  fDPC_container_log->SetVisAttributes(G4VisAttributes::Invisible);

  // place components relative to scintillator face at 0
//  fDPC_scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-fScintFacePos),
//				N2_container_log,Append(sd,"N2_vol_phys_"),fDPC_container_log,false,0);

  const G4double mwpc_pos = -(2*fmwpcContainer_halfZ)/2.-fDPC_backwin_frame_thick-(fN2_volume_Z/2.+fScintFacePos);
//  fDPC_mwpc_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpc_pos),
//  					WCham_container_log,Append(sd,"mwpcContainer_"),fDPC_container_log,false,0);
  fMyTranslation[2] += mwpc_pos;	// I literally don't understand this syntax at all.
					// But Michael had it and causes no errors so it must work.
  const G4double entrance_pos = mwpc_pos-((2*fmwpcContainer_halfZ)+entrance_section_length)/2;
  fDPC_entrance_face_pos = entrance_pos - 0.5*entrance_section_length;
  fDPC_entrance_win_pos = entrance_pos + 0.5*entrance_section_length;
//  fDPC_mwpc_entrance_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,entrance_pos),
 // 						fDPC_mwpc_entrance_log,Append(sd,"mwpc_entrance"),fDPC_container_log,false,0);

  // aluminum exit window and N2 volume at back of gas box
  G4Tubs* mwpc_exit_tube = new G4Tubs("mwpc_exit_tube",fmwpc_exit_R,fDPC_detPackageRadius,0.5*fDPC_backwin_frame_thick,0.,2*M_PI);
  fDPC_mwpc_exit_log = new G4LogicalVolume(mwpc_exit_tube, Al, Append(sd,"mwpc_exit_log_"));
  fDPC_mwpc_exit_log->SetVisAttributes(visMWPCExit);
  fDPC_exit_frame_pos = mwpc_pos+((2*fmwpcContainer_halfZ)+fDPC_backwin_frame_thick)/2;
//  fDPC_mwpc_exit_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,fDPC_exit_frame_pos),
//					fDPC_mwpc_exit_log,Append(sd,"mwpc_exit_"),fDPC_container_log,false,0);

  G4Tubs* mwpc_exit_N2_tube = new G4Tubs("mwpc_exit_N2_tube",0,fmwpc_exit_R,0.5*fDPC_backwin_frame_thick,0.,2*M_PI);
  fDPC_mwpc_exit_N2_log = new G4LogicalVolume(mwpc_exit_N2_tube, WCNitrogen, Append(sd,"mwpc_exit_N2_log_"));
  fDPC_mwpc_exit_N2_log->SetVisAttributes(G4VisAttributes::Invisible);
//  fDPC_mwpc_exit_N2_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,fDPC_exit_frame_pos),
//					fDPC_mwpc_exit_N2_log,Append(sd,"mwpc_exit_"),fDPC_container_log,false,0);
  // material behind detector
  const G4double backstuff_thick = 1.*inch;
  G4Tubs* backstuff_tube = new G4Tubs(Append(sd,"backstuff_tube"),0,fDPC_detPackageRadius,backstuff_thick*0.5,0.,2*M_PI);
  fDPC_backstuff_log = new G4LogicalVolume(backstuff_tube,SS304,Append(sd,"backstuff_log_"));
//  fDPC_backstuff_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,detPackageHalfZ-0.5*backstuff_thick),
  //						fDPC_backstuff_log,Append(sd,"backstuff_"),fDPC_container_log,false,0);
  //----- End of Detector Package Construction -----//
//----- Return to Detector Construction. Port the remainder of the code -----//
//  fEntranceToCathodes += fMWPCBowing;	// nothing changes since defaults to 0

  G4cout << "Complete construction of components. Begin placement of final DPC" << G4endl;

  G4RotationMatrix* sideFlip = new G4RotationMatrix();
  fDetOffset = G4ThreeVector();
  G4ThreeVector sideTrans;
  if(sd == 0)				// also doesn't seem to matter since fDetRot defaults to 0
  {
    sideFlip->rotateZ(fDetRot*(-1)*rad);
    sideTrans = G4ThreeVector(0.,0.,(-1)*(2.2*m-0))+fDetOffset*(-1);
  }							// there are two getScintFacePos() methods
  else if(sd == 1)					// One in DetectorPackageConstruction
  {							// One in ScintillatorConstruction
    sideFlip->rotateZ(fDetRot*(+1)*rad);
    sideTrans = G4ThreeVector(0.,0.,(+1)*(2.2*m-0))+fDetOffset*(+1);
  }
  if(sd == 0)
  {
    sideFlip->rotateY(M_PI*rad);
  }

//  detPackage_phys[sd] = new G4PVPlacement(sideFlip,sideTrans, fDPC_container_log,
//						Append(sd,"detPackage_phys_"),experimentalHall_log,false,0);

	// these lines below redefine these rotations/translations for the second iteration of sd
  fMyRotation = sideFlip;	// equivalent to dets[sd].mwpc.myRotation = sideFlip;
  fMyTranslation = (*sideFlip)(fMyTranslation);	// don't understand this line at all
  fMyTranslation += sideTrans;

//  WCham_container_log -> SetUserLimits(UserGasLimits);
//  winIn_log -> SetUserLimits(UserSolidLimits);
//  winOut_log -> SetUserLimits(UserSolidLimits);
//  kevStrip_log -> SetUserLimits(UserSolidLimits);

//  SetScoringVolumes(scint_log, 0, sd);
//  SetScoringVolumes(gas_log, 1, sd);	// for now, this is the MWPC SD. Will be modified later.

/*  if(fCallFieldConstructor)
  {
    fMyBField = fpMagField;
    ConstructMWPCField();
  } */
  }	// end of sd for loop which makes multiple detector packages
                                                                


//  G4cout << "Print check that we registered correctly: fScoreVol1 = " << fScoreVol1 -> GetName()
//         << ", fScoreVol2 = " << fScoreVol2 -> GetName()
//         << ", fScoreVol3 = " << fScoreVol3 -> GetName()
//         << ", fScoreVol4 = " << fScoreVol4 -> GetName() << G4endl;

  return experimentalHall_phys;
}

string DetectorConstruction::Append(int i, string str)
{
  stringstream newString;
  newString << str << i;
  return newString.str();
}

#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
//#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixMixedStepper.hh"

void DetectorConstruction::ConstructField()
{
  G4double BfieldArray[3] = {0, 0, 0};

  if(!fpMagField)
  {
    cout << "##### Constructing Detector Field #####" << endl;
    fpMagField = new Field();

    ofstream fieldprint;
    fieldprint.open("fieldValues.txt", ios::app);

    G4double location[3] = {1*cm, 1*cm, 0*cm};
    for(int i = 0; i < 1000; i++)
    {
      G4double zValue = (i-500)*cm;
      location[2] = zValue;
      fpMagField -> GetFieldValue(location, BfieldArray);
      fieldprint << "B(" << location[0]/cm << ", " << location[1]/cm << ", " << location[2]/cm << ") cm  = ("
	<< BfieldArray[0]/tesla << ", " << BfieldArray[1]/tesla << ", " << BfieldArray[2]/tesla << ") tesla \n";
    }
    fieldprint.close();

    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr -> SetDetectorField(fpMagField);
    fieldMgr -> CreateChordFinder(fpMagField);

    G4MagIntegratorStepper* pStepper;
    G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(fpMagField); // equation of motion in magnetic field
    //pStepper = new G4ClassicalRK4 (fEquation); // general case for "smooth" EM fields
    //pStepper = new G4SimpleHeum( fEquation ); // for slightly less smooth EM fields
    //pStepper = new G4HelixHeum( fEquation ); // for "smooth" pure-B fields
    //pStepper = new G4HelixImplicitEuler( fEquation ); // for less smooth pure-B fields; appears ~50% faster than above
    //pStepper = new G4HelixSimpleRunge( fEquation ); // similar speed to above
    //pStepper = new G4HelixExplicitEuler( fEquation ); // about twice as fast as above
    pStepper = new G4HelixMixedStepper(fEquation,6); // avoids "Stepsize underflow in Stepper" errors
    fieldMgr->GetChordFinder()->GetIntegrationDriver()->RenewStepperAndAdjust(pStepper);

    fieldMgr->GetChordFinder()->SetDeltaChord(100.0*um);
    fieldMgr->SetMinimumEpsilonStep(1e-6);
    fieldMgr->SetMaximumEpsilonStep(1e-5);
    fieldMgr->SetDeltaOneStep(0.1*um);
    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetMaxLoopCount(INT_MAX);

    fCallFieldConstructor = true;
  }
  else
  {
    fCallFieldConstructor = false;
  }
}

void DetectorConstruction::ConstructMWPCField()
{
  G4cout << "Setting up wirechamber electromagnetic field...";

  // local field manager
  G4FieldManager* localFieldMgr = new G4FieldManager();
  localFieldMgr->SetDetectorField(fMyBField);

  // equation of motion, stepper for field
  G4EqMagElectricField* pEquation = new G4EqMagElectricField(fMyBField);
  G4ClassicalRK4* pStepper = new G4ClassicalRK4(pEquation,8);
  G4MagInt_Driver* pIntgrDriver = new G4MagInt_Driver(0.01*um,pStepper,pStepper->GetNumberOfVariables());
  G4ChordFinder* pChordFinder = new G4ChordFinder(pIntgrDriver);
  localFieldMgr->SetChordFinder(pChordFinder);

  // accuracy settings
  localFieldMgr->GetChordFinder()->SetDeltaChord(10*um);
  localFieldMgr->SetMinimumEpsilonStep(1e-6);
  localFieldMgr->SetMaximumEpsilonStep(1e-5);
  localFieldMgr->SetDeltaOneStep(0.1*um);

  // apply field manager to wirechamber and all daughter volumes
//  WCham_container_log->SetFieldManager(localFieldMgr,true);

  G4cout << " Done." << G4endl;
}

/*void DetectorConstruction::GetFieldValue(G4double Point[4], G4double* Bfield)
{
  // set magnetic field
  if(fMyBField) fMyBField->GetFieldValue(Point,Bfield);
  else Bfield[0]=Bfield[1]=Bfield[2]=0;

  if(!fE0) { Bfield[3]=Bfield[4]=Bfield[5]=0; return; }

  // local position
  G4ThreeVector localPos = G4ThreeVector(Point[0],Point[1],Point[2])-fMyTranslation;
  if(fMyRotation) localPos = (*fMyRotation)(localPos);

  // electric field components
  G4ThreeVector E(0,0,0);
  double l = localPos[2];
  if(fabs(l)<fL)
  {
    double a = localPos[0]/fd;
    a = (a-floor(a)-0.5)*fd;
    if(a*a+l*l > fr*fr)
    {
      double denom = cosh(2*M_PI*l/fd)-cos(2*M_PI*a/fd);
      E[2] = fE0*sinh(2*M_PI*l/fd)/denom;
      E[0] = fE0*sin(2*M_PI*a/fd)/denom;
    }
  }
  // return to global coordinates
  if(fMyRotation) E = fMyRotation->inverse()(E);
  for(unsigned int i=0; i<3; i++) Bfield[3+i] = E[i];
} */

/*void DetectorConstruction::setMWPCPotential(G4double Vanode)
{
  fE0 = M_PI*Vanode/fd/log(sinh(M_PI*fL/fd)/sinh(M_PI*fr/fd));
  G4cout << "Wirechamber voltage set to " << Vanode/volt <<" V => fE0 = " << fE0/(volt/cm) << " V/cm" << G4endl;
}*/

void DetectorConstruction::SetScoringVolumes(G4LogicalVolume* vol, int type, int location)
{
  G4cout << "Setting scoring volume using " << vol->GetName()
	 << " i.e. of type " << type << " and location " << location << G4endl;

  if(location == 0)
  {
    if(type == 0)
    {
      fScoreVol1 = vol;
    }
    else if(type == 1)
    {
      fScoreVol2 = vol;
    }
  }
  else if(location == 1)
  {
    if(type == 0)
    {
      fScoreVol3 = vol;
    }
    else if(type == 1)
    {
      fScoreVol4 = vol;
    }
  }
  else if((location != 0) || (location != 1) || (type != 0) || (type != 1))
  {
    G4cout << "Yabai. Input flags make no sense" << G4endl;
  }



}
