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

using	namespace	std;

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }


DetectorConstruction::~DetectorConstruction()
{ }

void DetectorConstruction::DefineMaterials()
{

  G4cout << "Initial Vacuum has been set to null pointer." << G4endl;
  Vacuum = NULL;

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

  setVacuumPressure(0);

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //
  // Envelope
  //
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // Shape 1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Cons* solidShape1 =
    new G4Cons("Shape1",
    shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
    shape1_phimin, shape1_phimax);

  G4LogicalVolume* logicShape1 =
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;
  G4Trd* solidShape2 =
    new G4Trd("Shape2",                      //its name
              0.5*shape2_dxa, 0.5*shape2_dxb,
              0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size

  G4LogicalVolume* logicShape2 =
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name

  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicShape2;

  //
  //always return the physical World
  //
  return physWorld;
}