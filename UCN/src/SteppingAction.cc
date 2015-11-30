#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fSV1(0),fSV2(0),fSV3(0),fSV4(0)
{}


SteppingAction::~SteppingAction()
{}


void SteppingAction::UserSteppingAction(const G4Step* step)
{

/*  G4String preStepName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  G4ThreeVector preStepPosition = step->GetPreStepPoint()->GetPosition();

  ofstream outfile;
  outfile.open("DebuggingGeometry.txt", ios::app);
  outfile << "Pre step volume is: " << preStepName <<"\t\t at position " << preStepPosition/m << "m \n";
  outfile.close();
*/

  const DetectorConstruction* detectorConstruction =
        static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  fSV1 = detectorConstruction -> GetScoringVolume1();
  fSV2 = detectorConstruction -> GetScoringVolume2();
  fSV3 = detectorConstruction -> GetScoringVolume3();
  fSV4 = detectorConstruction -> GetScoringVolume4();

  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  if((volume != fSV1) && (volume != fSV2) && (volume != fSV3) && (volume != fSV4))
  {
    return;	// if not scoring volume, GTFO
  }

  G4double edepStep = step->GetTotalEnergyDeposit();

  if(volume == fSV1)
  {
    fEventAction -> AddEdep(edepStep, 0, 0);
  }
  else if(volume == fSV2)
  {
    fEventAction -> AddEdep(edepStep, 1, 0);
  }
  else if(volume == fSV3)
  {
    fEventAction -> AddEdep(edepStep, 0, 1);
  }
  else if(volume == fSV4)
  {
    fEventAction -> AddEdep(edepStep, 1, 1);
  }

}

