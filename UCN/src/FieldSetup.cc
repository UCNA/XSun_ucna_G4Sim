#include "FieldSetup.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

FieldSetup::FieldSetup()
 : fFieldManager(0), fEastMWPCFieldManager(0), fWestMWPCFieldManager(0),
   fChordFinder(0), fEastMWPCChordFinder(0), fWestMWPCChordFinder(0),
   fEquation(0), fEastMWPCEquation(0), fWestMWPCEquation(0),
   fMagneticField(0), fEastMWPCMagneticField(0), fWestMWPCMagneticField(0),
   fStepper(0), fEastMWPCStepper(0), fWestMWPCStepper(0)
{
  G4cout << "Entering field constructor" << G4endl;

  fMagneticField = new G4UniformMagField(G4ThreeVector(0.0*tesla,
                                                       0.0*tesla, // 0.5*tesla,
                                                       1.0*tesla));
  fEastMWPCMagneticField = new G4UniformMagField(G4ThreeVector(0.0*tesla,
                                                               0.0*tesla, // 0.5*tesla,
                                                               1.0*tesla));
  fWestMWPCMagneticField = new G4UniformMagField(G4ThreeVector(0.0,
							       0.0,
							       1.0*tesla));

  fEquation = new G4Mag_UsualEqRhs(fMagneticField);
  fEastMWPCEquation = new G4Mag_UsualEqRhs(fEastMWPCMagneticField);
  fWestMWPCEquation = new G4Mag_UsualEqRhs(fWestMWPCMagneticField);

  fMinStep     = 0.25*mm ; // minimal step of 1 mm is default

  fFieldManager = GetGlobalFieldManager();
  fEastMWPCFieldManager = new G4FieldManager();
  fWestMWPCFieldManager = new G4FieldManager();

  UpdateField();
}

FieldSetup::~FieldSetup()
{
  delete fMagneticField;
  delete fChordFinder;
  delete fStepper;
}

void FieldSetup::UpdateField()
{
  if (fStepper) delete fStepper;
  // These are the steppers (not really) used in Michael's code. Go with it.
  fStepper = new G4HelixSimpleRunge(fEquation);	 	// this needs to be HelixMixedStepper, apparently second argument 6
  fEastMWPCStepper = new G4ClassicalRK4(fEastMWPCEquation);	// I don't understand this local stepper.
  fWestMWPCStepper = new G4ClassicalRK4(fWestMWPCEquation);	// In Michael's code, it's object type G4ClassicalRK4*
							// Which doesn't match. And I don't know if EoM works the same.

  fFieldManager -> SetDetectorField(fMagneticField);	// set some field managers
  fEastMWPCFieldManager -> SetDetectorField(fEastMWPCMagneticField);
  fWestMWPCFieldManager -> SetDetectorField(fWestMWPCMagneticField);

  if (fChordFinder) delete fChordFinder;	// check to see if these objects exist on update
  if (fEastMWPCChordFinder) delete fEastMWPCChordFinder;
  if (fWestMWPCChordFinder) delete fWestMWPCChordFinder;

  fChordFinder = new G4ChordFinder(fMagneticField, fMinStep, fStepper);	// Instead of chord finder, need to use create chord
  fEastMWPCChordFinder = new G4ChordFinder(fEastMWPCMagneticField, fMinStep, fEastMWPCStepper);
  fWestMWPCChordFinder = new G4ChordFinder(fWestMWPCMagneticField, fMinStep, fWestMWPCStepper);

  fFieldManager->SetChordFinder(fChordFinder);
  fEastMWPCFieldManager->SetChordFinder(fEastMWPCChordFinder);
  fWestMWPCFieldManager->SetChordFinder(fWestMWPCChordFinder);

  G4cout << "Field has been updated." << G4endl;
}

void FieldSetup::SetFieldValue(G4double fieldStrength)
{
  if (fMagneticField) delete fMagneticField;

  G4ThreeVector fieldVector(0., 0., fieldStrength);

  if(fieldVector != G4ThreeVector(0., 0., 0.))
  {
    fMagneticField = new G4UniformMagField(fieldVector);
  }
  else
  {
    fMagneticField = 0;
  }

  UpdateField();	// simplest way to reset everything.

  // or we update the field manager and equation of motion with pointer to new field
  // Note: new field comes from setting new field value.
//  GetGlobalFieldManager()->SetDetectorField(fMagneticField);
//  fEquation -> SetFieldObj(fMagneticField);

}

G4FieldManager* FieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}






