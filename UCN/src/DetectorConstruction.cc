#include "DetectorConstruction.hh"

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
#include <G4Box.hh>
#include <G4Tubs.hh>			// nothing used by decay trap construction
#include <G4VPhysicalVolume.hh>		// not using wiggle sheet
#include <G4LogicalVolume.hh>		// or silicon detector construction
#include <G4ThreeVector.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4RotationMatrix.hh>
#include <G4VisAttributes.hh>
#include <G4SystemOfUnits.hh>

#include <cassert>			// scintillator construction classes
#include <G4Polycone.hh>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0), fScintStepLimit(1)/*,	// note: fScintStepLimit initialized here
  experimentalHall_log(0), experimentalHall_phys(0),
  container_log(0), holder_phys(0), ring_phys(0),
  window_log(0), window_phys(0),
  source_phys(0), decayTube_log(0), plug_log(0) */
{
  coating_log[0] = NULL;	// source holder member variables should go here or in line above
  coating_log[1] = NULL;	// but then ordering is wrong with the materials
  coating_phys[0] = NULL;	// so it's done right before source holder is made.
  coating_phys[1] = NULL;
}


DetectorConstruction::~DetectorConstruction()
{ }

void DetectorConstruction::DefineMaterials()
{

  G4cout << "Initial Vacuum has been set to null pointer." << G4endl;
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

  cout << "\n \n \n \n Constructing Materials! \n \n \n \n" << endl;

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

  G4cout << "\n Completed construction of materials. Exiting ConstructMaterials(). \n" << G4endl;

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
  UserSolidLimits->SetMaxAllowedStep(fScintStepLimit*mm);	// hard coded scint step limit to mm units

  G4bool checkOverlaps = true;					// Just a flag to check overlaps of geometry

  //----- experimental Hall -----//
  const G4double expHall_halfx=1.0*m;
  const G4double expHall_halfy=1.0*m;
  const G4double expHall_halfz=4.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_halfx,expHall_halfy,expHall_halfz);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"World_Log");
  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  experimentalHall_log->SetUserLimits(UserCoarseLimits);
  experimentalHall_phys = new G4PVPlacement(NULL,G4ThreeVector(),"World_Phys", experimentalHall_log,0,false,0,checkOverlaps);


  //----- source holder object -----//
  fSourceWindowThick = 4.7*um;		// class initialization variables.
  fSourceCoatingThick = 0.1*um;		// really should go in constructor, but
  fSourceWindowMat = Mylar;		// materials might be defined in wrong order
  fSourceCoatingMat = Al;
  fSourceHolderThickness = (3./16.)*inch;

  fSourceHolderPos = G4ThreeVector(0,0,0);	// initialize the position. Should be done in Constructor but w/e

  const G4double SourceRingRadius = 0.5*inch;
  const G4double SourceWindowRadius = SourceRingRadius-3.*mm;
  const G4double SourceRingThickness = 3.2*mm; // suspiciously close to 1/8 in
  const G4double SourceHolderHeight = 1.5*inch;
  const G4double SourceHolderWidth = 1.5*inch;

  // source holder container
  G4Box* holder_box = new G4Box("source_holder_box",0.5*SourceHolderWidth,0.5*SourceHolderHeight,0.5*fSourceHolderThickness);
  container_log = new G4LogicalVolume(holder_box,Vacuum,"source_container_log");
  container_log->SetVisAttributes(G4VisAttributes::Invisible);	// keep this guy!

  // source holder paddle
  G4Tubs* holder_hole = new G4Tubs("source_holder_hole",0.,SourceRingRadius,fSourceHolderThickness,0.,2*M_PI);
  G4SubtractionSolid* holder = new G4SubtractionSolid("source_holder", holder_box, holder_hole);
  G4LogicalVolume* holder_log = new G4LogicalVolume(holder,Brass,"source_holder_log");
//  holder_log->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.7,0,0.5)));
  holder_log->SetVisAttributes(G4VisAttributes::Invisible); // use above
  holder_phys = new G4PVPlacement(NULL,G4ThreeVector(),holder_log,"source_holder_phys",container_log,false,0);

  // sealed source foil
  G4Tubs* window_tube = new G4Tubs("window_tube",0.,SourceWindowRadius,fSourceWindowThick,0.,2*M_PI);
  window_log = new G4LogicalVolume(window_tube,fSourceWindowMat,"source_window_log");
  G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
//  window_log->SetVisAttributes(visWindow);
  window_log->SetVisAttributes(G4VisAttributes::Invisible);	//use above
  window_phys = new G4PVPlacement(NULL,G4ThreeVector(),window_log,"source_window_phys",container_log,false,0);

  // source foil coating
  G4Tubs* coating_tube = new G4Tubs("source_coating_tube", 0., SourceWindowRadius, fSourceCoatingThick*0.5, 0., 2*M_PI);
	//	flag: 0 means EAST
	//            1 means WEST
	//	and int sd is an instance of SIDE sd (in Mendenhall code) which is just a renaming/alias
  for(int sd = 0; sd <= 1; sd++)
  {
    coating_log[sd] = new G4LogicalVolume(coating_tube, fSourceCoatingMat, Append(sd, "source_coating_log_"));
//    coating_log[sd]->SetVisAttributes(new G4VisAttributes(G4Colour(0,1,0,0.5)));
    coating_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);	// use above

    if(sd == 0)
    {
      coating_phys[sd] = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*(fSourceWindowThick+fSourceCoatingThick*0.5)),
				coating_log[sd],Append(sd,"source_coating_phys_"),container_log,false,0);
    }
    if(sd == 1)
    {
      coating_phys[sd] = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*(fSourceWindowThick+fSourceCoatingThick*0.5)),
                                coating_log[sd],Append(sd,"source_coating_phys_"),container_log,false,0);
    }
  }

  // source retaining ring
  G4Tubs* ring_tube = new G4Tubs("source_ring_tube",SourceWindowRadius,SourceRingRadius,SourceRingThickness/2,0.,2*M_PI);
  G4LogicalVolume* ring_log = new G4LogicalVolume(ring_tube,Al,"source_ring_log");
//  ring_log->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.5)));
  ring_log->SetVisAttributes(G4VisAttributes::Invisible);	// use above
  ring_phys = new G4PVPlacement(NULL,G4ThreeVector(),ring_log,"source_ring_phys",container_log,false,0);

  // ----- end of source holder class code.

  source_phys = new G4PVPlacement(NULL,fSourceHolderPos,container_log,"source_container_phys",experimentalHall_log,false,0);


  //----- geometry-dependent settings -----//
  // Note: we only need geometry settings "thinFoil"


  //----- Detector Decay Trap construction -----//

  fTrapWindowThick = 0.180*um;
  fTrapCoatingThick = 0.150*um;
  fTrapIRtrap = 2.45*inch;
  fTrapDecayTube_Wall = 2*mm;
  fTrapIRcollimator = 2.3*inch;
  fTrapTubeMat = Cu;
  fTrapCollimatorMat = Polyethylene;
  fTrapWindowMat = Mylar;
  fTrapCoatingMat = Be;

//  fCrinkleAngle = 0.;	// can be set independently (in Mendenhall code) using user input
			// For now, we are not using this at all.

  // decay tube construction
  G4double decayTube_OR = fTrapIRtrap+fTrapDecayTube_Wall;	// OR = outer radius
  G4double decayTube_Length = 3.0*m;
  G4Tubs* decayTube_tube = new G4Tubs("decayTube_tube",fTrapIRtrap, decayTube_OR, decayTube_Length/2.,0.,2*M_PI);
  decayTube_log = new G4LogicalVolume(decayTube_tube,Cu,"decayTube_log");
//  decayTube_log->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
  decayTube_log->SetVisAttributes(G4VisAttributes::Invisible);	// use above
  new G4PVPlacement(NULL,G4ThreeVector(),decayTube_log,"decayTube",experimentalHall_log,false,0);

  // plug - not sure if this is even used or what it is
  if(false)
  {
    G4Tubs* plug_tube = new G4Tubs("plug_tube",fTrapIRtrap-1*mm, fTrapIRtrap, 2.*cm, -2*cm/fTrapIRtrap, 2*cm/fTrapIRtrap);
    plug_log = new G4LogicalVolume(plug_tube,Cu,"plug_log");
    plug_log->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
    new G4PVPlacement(NULL,G4ThreeVector(),plug_log,"plug",experimentalHall_log,false,0);
  }

  // trap windows, collimator, monitors
  G4double thicknessOfTrapWindow = fTrapWindowThick + fTrapCoatingThick;
  G4double collimator_thick = 0.8*inch;
  G4Tubs* trap_win_tube = new G4Tubs("trap_win_tube",0.,decayTube_OR,thicknessOfTrapWindow/2.,0.,2*M_PI);
  G4Tubs* mylarTube = new G4Tubs("mylarTube",0.,decayTube_OR,fTrapWindowThick/2.,0.,2*M_PI);
  G4Tubs* beTube = new G4Tubs("beTube",0.,decayTube_OR,fTrapCoatingThick/2.,0.,2*M_PI);
  G4Tubs* collimatorTube = new G4Tubs("collimatorTube",fTrapIRcollimator,fTrapIRcollimator+collimator_thick,collimator_thick/2.,0.,2*M_PI);
  G4Tubs* collimatorBackTube = new G4Tubs("collimatorBackTube",decayTube_OR+1.*mm,fTrapIRcollimator+collimator_thick,collimator_thick,0.,2*M_PI);

  G4double beWinPosZ = -thicknessOfTrapWindow/2.+fTrapCoatingThick/2.;
  G4double mylarWinPosZ = thicknessOfTrapWindow/2.-fTrapWindowThick/2.;

  G4double trap_winPosZ=(decayTube_Length+thicknessOfTrapWindow)/2.;
  G4double trap_monitor_thickness = 1.0*mm;
  G4double trap_monitor_posZ = 0.5*m;
  G4Tubs* trap_monitor_tube = new G4Tubs("trap_monitor_tube",0.,fTrapIRtrap,trap_monitor_thickness/2,0.,2*M_PI);

  G4RotationMatrix* wigRot = new G4RotationMatrix;
  wigRot->rotateY(M_PI/2.*rad);
  wigRot->rotateX(M_PI/2.*rad);

  for(int sd = 0; sd <= 1; sd++)	// exact same "legend" as source holder loop
  {
    trap_win_log[sd] = new G4LogicalVolume(trap_win_tube,Vacuum,Append(sd,"trap_win_log_"));
//    trap_win_log[sd]->SetVisAttributes(visWindow);
    trap_win_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);	// get rid of, use above
	// note the weird indexing because there is an if-else statement dependent on fCrinkleAngle
      if(sd == 0)
      {
        new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*trap_winPosZ),trap_win_log[sd],
								Append(sd,"trap_win_"),experimentalHall_log,false,0);
      }
      if(sd == 1)
      {
        new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*trap_winPosZ),trap_win_log[sd],
                                                                Append(sd,"trap_win_"),experimentalHall_log,false,0);
      }

    mylar_win_log[sd] = new G4LogicalVolume(mylarTube, fTrapWindowMat, Append(sd,"mylar_win_log_"));
    be_win_log[sd] = new G4LogicalVolume(beTube, fTrapCoatingMat,Append(sd,"be_win_log_"));
	// add some visualizations to these logical volumes
    mylar_win_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);	// get rid of and replace these visualizations
    be_win_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);

    if(sd == 0)
    {
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*mylarWinPosZ),mylar_win_log[sd],
                                        Append(sd,"mylar_win_"),trap_win_log[sd],false,0);
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*beWinPosZ),be_win_log[sd],
                                        Append(sd,"be_win_"),trap_win_log[sd],false,0);
    }
    if(sd == 1)
    {
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*mylarWinPosZ),mylar_win_log[sd],
                                        Append(sd,"mylar_win_"),trap_win_log[sd],false,0);
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*beWinPosZ),be_win_log[sd],
                                        Append(sd,"be_win_"),trap_win_log[sd],false,0);
    }

    G4double collimatorPosZ = (decayTube_Length+collimator_thick)/2.;
    collimatorPosZ += (thicknessOfTrapWindow)/2.;
    collimator_log[sd] = new G4LogicalVolume(collimatorTube, fTrapCollimatorMat, Append(sd,"collimator_log_"));
    collimator_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);	// get rid of

    G4double collimatorBackZ = decayTube_Length/2.-collimator_thick;
    collimatorBack_log[sd] = new G4LogicalVolume(collimatorBackTube, fTrapCollimatorMat, Append(sd,"collimatorBack_log_"));
    collimatorBack_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);	// get rid of

    trap_monitor_log[sd] = new G4LogicalVolume(trap_monitor_tube,Vacuum,Append(sd,"trap_monitor_log_"));
    trap_monitor_log[sd]->SetVisAttributes(G4VisAttributes::Invisible);	// get rid of

    if(sd == 0)
    {
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*collimatorPosZ),collimator_log[sd],
						Append(sd,"collimator_"),experimentalHall_log,false,0);
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*collimatorBackZ),collimatorBack_log[sd],
						Append(sd,"collimatorBack_"),experimentalHall_log,false,0);
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(-1)*trap_monitor_posZ),
						trap_monitor_log[sd],Append(sd,"trap_monitor_"),experimentalHall_log,false,0);
    }
    if(sd == 1)
    {
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*collimatorPosZ),collimator_log[sd],
                                                Append(sd,"collimator_"),experimentalHall_log,false,0);
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*collimatorBackZ),collimatorBack_log[sd],
                                                Append(sd,"collimatorBack_"),experimentalHall_log,false,0);
      new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(1)*trap_monitor_posZ),
                                                trap_monitor_log[sd],Append(sd,"trap_monitor_"),experimentalHall_log,false,0);
    }
  }

  //----- End of decay trap construction code. Return to detector construction -----//

  //********** To do remaining geometry, need to import the ingredients for Detector Package Construction ****/////

  //----- Begin scintillator construction -----//
  fScint_Radius = 7.5*cm;	// initialize some constructor variables
  fBacking_Radius = 10*cm;
  fScint_thick = 3.5*mm;
  fDead_thick = 3.0*um;
  fBacking_thick = 1.*inch;
  fLightguide_thick = 1.0*cm;

  assert(fBacking_Radius > fScint_Radius);
  assert(fLightguide_thick >= fScint_thick);
  fN2_volume_Z = fLightguide_thick+fBacking_thick; // length of N2 volume
  fScintFacePos = -fN2_volume_Z/2;


  int sd = 0;	// this will become a for-loop to take values 0 and 1.

	// overall container
  G4Tubs* N2_vol_tube = new G4Tubs(Append(sd,"N2_vol_tube_"),0.,fBacking_Radius,fN2_volume_Z/2,0.,2*M_PI);
  N2_container_log = new G4LogicalVolume(N2_vol_tube,WCNitrogen,Append(sd,"N2_vol_log_"));
  N2_container_log->SetVisAttributes(G4VisAttributes::Invisible);	// keep this guy!

	// dead layer
  G4Tubs* Dscint_tube = new G4Tubs(Append(sd,"Dscint_tube_"),0.,fScint_Radius,fDead_thick/2,0.,2*M_PI);
  Dscint_log = new G4LogicalVolume(Dscint_tube,Sci,Append(sd,"Dscint_log_"));
  G4VisAttributes* visDScint= new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));
  Dscint_log->SetVisAttributes(visDScint);	// can't see anything cause didn't place it
  Dscint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-(fN2_volume_Z-fDead_thick)/2.),
  					Dscint_log,Append(sd,"Dscint_phys_"),N2_container_log,false,0);

	// scintillator
  G4Tubs* scint_tube = new G4Tubs("scint_tube",0.,fScint_Radius,(fScint_thick-fDead_thick)/2.,0.,2*M_PI);
  scint_log = new G4LogicalVolume(scint_tube,Sci,Append(sd,"scint_log_"));
  G4VisAttributes* visScint= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.2));
  scint_log->SetVisAttributes(visScint);
  scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-fN2_volume_Z/2+fDead_thick+(fScint_thick-fDead_thick)/2.),
  					scint_log,Append(sd,"scint_phys_"),N2_container_log,false,0);

	// light guides around/behind scintillator
  const G4double zPlane[] = {0., fScint_thick, fScint_thick, fLightguide_thick};
  G4double lg_rad = fScint_Radius-(fLightguide_thick-fScint_thick);
  const G4double rInner[] = {fScint_Radius, fScint_Radius, lg_rad, lg_rad};
  const G4double rOuter[] = {fBacking_Radius, fBacking_Radius, fBacking_Radius, fBacking_Radius};

  G4Polycone* lightguide_polycone = new G4Polycone("lightguide_polycone",0,2*M_PI,4,zPlane,rInner,rOuter);
  lightguide_log = new G4LogicalVolume(lightguide_polycone,Sci,Append(sd,"lightguide_log_"));
  G4VisAttributes* visLG = new G4VisAttributes(G4Colour(0.0,1.0,0.5,0.2));
  lightguide_log->SetVisAttributes(visLG);
  lightguide_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-fN2_volume_Z/2.),
  						lightguide_log,Append(sd,"scint_phys_"),N2_container_log,false,0);

	// backing veto
  G4Tubs* backing_tube = new G4Tubs("backing_tube",0.,fBacking_Radius,fBacking_thick/2.,0.,2*M_PI);
  backing_log = new G4LogicalVolume(backing_tube,Sci,Append(sd,"backing_log_"));
  G4VisAttributes* visBacking= new G4VisAttributes(G4Colour(0.0,0.0,1,0.2));
  backing_log->SetVisAttributes(visBacking);
  backing_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(fN2_volume_Z-fBacking_thick)/2.),
  					backing_log,Append(sd,"backing_phys_"),N2_container_log,false,0);
  //----- End scintillator construction -----//

  // The sensible way to do this is to port as much code as possible before the detector package construction for loop

  //----- Begin Wire Volume Construction -----//

  fAnode_R = 5*um;	// initialize wire volume construction constructor variables
  fCathode_R = 25*um;
  fPlating_thick = 0.2*um;
  fSpacing = 2.54*mm;
  fNbOfWires = 64;
  fPlaneSpacing = 1*cm;

  fMWPCGas = WCPentane;	// NOTE: this is declared in WirechamberConstruction, not the volume constr.

  G4Material* fCathodeWireMaterial = Al;
  G4Material* fAnodeWireMaterial = Wu;
  G4Material* fCathodePlateMaterial = Au;

  	// wireplane box size
  G4double wireplane_half_width = fNbOfWires*fSpacing/2;

	//effective mwpc gas volume containing the cathodes and anode
  G4Box* mwpc_box = new G4Box("mwpc_box",wireplane_half_width,wireplane_half_width,fPlaneSpacing);
  gas_log = new G4LogicalVolume(mwpc_box,fMWPCGas,Append(sd,"mwpc_log_"));
  gas_log->SetVisAttributes(G4VisAttributes::Invisible);

	// anode, cathode wire containers... note these are "on their side" to allow wireplane parametrization,
	// and will be rotated later
  G4Box* cathContainer_box = new G4Box("cathContainer_box",wireplane_half_width,fCathode_R,wireplane_half_width);
  G4Box* anodeContainer_box = new G4Box("anodeContainer_box",wireplane_half_width,fAnode_R,wireplane_half_width);

	// anode, cathode wires and surrounding gas
  G4Tubs* cathplate_tube = new G4Tubs("cathplate_tube",fCathode_R-fPlating_thick,fCathode_R,wireplane_half_width,0,2*M_PI);
  G4Tubs* cathode_tube = new G4Tubs("cathode_tube",0,fCathode_R-fPlating_thick,wireplane_half_width,0.,2*M_PI);
  G4Tubs* anode_tube = new G4Tubs("anode_tube",0,fAnode_R,wireplane_half_width,0.,2*M_PI);
  G4Box* cathodeSegmentBox = new G4Box("cathodeSegmentBox",fSpacing/2,fCathode_R,wireplane_half_width);
  G4Box* anodeSegmentBox = new G4Box("anodeSegmentBox",fSpacing/2,fAnode_R,wireplane_half_width);

	// Rotate 90 degrees around X axis
  G4RotationMatrix* xRot90 = new G4RotationMatrix;
  xRot90->rotateX(M_PI/2.*rad);
	// Rotate 90 degrees around X then Z axis (Y axis in object coordinates)
  G4RotationMatrix* xzRot90 = new G4RotationMatrix;
  xzRot90->rotateX(M_PI/2.*rad);
  xzRot90->rotateY(M_PI/2.*rad);

  G4VisAttributes* visCathWires = new G4VisAttributes(G4Colour(1,0.7,0,0.8));
  G4VisAttributes* visAnodeWires = new G4VisAttributes(G4Colour(1,0.3,0,0.8));

	// anode, cathode segments
  cathSeg_log = new G4LogicalVolume(cathodeSegmentBox,fMWPCGas,Append(sd,"cathSeg_log_"));
  cathSeg_log->SetVisAttributes(G4VisAttributes::Invisible);
  anodeSeg_log = new G4LogicalVolume(anodeSegmentBox,fMWPCGas,Append(sd,"anodeSeg_log_"));
  anodeSeg_log->SetVisAttributes(G4VisAttributes::Invisible);
  cathode_wire_log = new G4LogicalVolume(cathode_tube,fCathodeWireMaterial,Append(sd,"cathode_log_"));
  cathode_wire_log->SetVisAttributes(visCathWires);
  cath_plate_log = new G4LogicalVolume(cathplate_tube,fCathodePlateMaterial,Append(sd,"cathplate_log_"));
  cath_plate_log->SetVisAttributes(visCathWires);
  anode_wire_log = new G4LogicalVolume(anode_tube,fAnodeWireMaterial,Append(sd,"anode_log_"));
  anode_wire_log->SetVisAttributes(visAnodeWires);

  new G4PVPlacement(NULL,G4ThreeVector(),cathode_wire_log, Append(sd,"cathode_wire_phys_"),cathSeg_log,true,0);
  new G4PVPlacement(NULL,G4ThreeVector(),cath_plate_log,Append(sd,"cath_plate_phys_"),cathSeg_log,true,0);
  new G4PVPlacement(NULL,G4ThreeVector(),anode_wire_log,Append(sd,"anode_wire_phys_"),anodeSeg_log,true,0);

	// anode, cathode plane container volumes
  G4LogicalVolume* cathContainer1_log = new G4LogicalVolume(cathContainer_box,fMWPCGas,Append(sd,"cathContainer1_log_"));
  G4LogicalVolume* cathContainer2_log = new G4LogicalVolume(cathContainer_box,fMWPCGas,Append(sd,"cathContainer2_log_"));
  G4LogicalVolume* anodContainer_log = new G4LogicalVolume(anodeContainer_box,fMWPCGas,Append(sd,"anodContainer_log_"));

  new G4PVPlacement(xRot90,G4ThreeVector(0.,0.,fCathode_R-fPlaneSpacing),
				cathContainer1_log,Append(sd,"cathContainer1_phys_"),gas_log,false,0);
  new G4PVPlacement(xzRot90,G4ThreeVector(0.,0.,fPlaneSpacing-fCathode_R),
				cathContainer2_log,Append(sd,"cathContainer2_phys_"),gas_log,false,0);
  new G4PVPlacement(xRot90,G4ThreeVector(0.,0.,0.),
				anodContainer_log,Append(sd,"anodContainer_phys_"),gas_log,false,0);

	// replicate segments into cathode, anode arrays
  new G4PVReplica(Append(sd,"cathode_array_1_"), cathSeg_log, cathContainer1_log, kXAxis, fNbOfWires, fSpacing);
  new G4PVReplica(Append(sd,"cathode_array_2_"), cathSeg_log, cathContainer2_log, kXAxis, fNbOfWires, fSpacing);
  new G4PVReplica(Append(sd,"anode_array_"), anodeSeg_log, anodContainer_log, kXAxis, fNbOfWires, fSpacing);



















  // Get nist material manager
  //G4NistManager* nist = G4NistManager::Instance();


  // Set Shape2 as scoring volume
  //
  fScoringVolume = experimentalHall_log;

  //
  //always return the physical World
  //
  return experimentalHall_phys;



}

string DetectorConstruction::Append(int i, string str)
{
  stringstream newString;
  newString << str << i;
  return newString.str();
}
