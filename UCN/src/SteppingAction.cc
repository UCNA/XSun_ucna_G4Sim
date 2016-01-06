#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include <time.h>

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{ }


SteppingAction::~SteppingAction()
{ }


void SteppingAction::UserSteppingAction(const G4Step* step)
{
/*
  // The below is used for debugging the geometry in a step-by-step fashion.
  // Not needed but saved in case will be used in future.
  G4String preStepName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  G4ThreeVector preStepPosition = step->GetPreStepPoint()->GetPosition();

  ofstream outfile;
  outfile.open("DebuggingGeometry.txt", ios::app);
  outfile << "Pre step volume is: " << preStepName <<"\t\t at position " << preStepPosition/m << "m \n";
  outfile.close();
*/

  // kill switch in case simulation runs too long
  G4int stepNo = step -> GetTrack() -> GetCurrentStepNumber();
  clock_t timeSpentSoFar = clock() - ((EventAction*)G4EventManager::GetEventManager()->GetUserEventAction())->GetStartTime();

  double time = ((double)timeSpentSoFar)/CLOCKS_PER_SEC;

  if(stepNo >= 2000000 || time > 60)
  {
    G4cout << "----> Tracking killed by computation time limit." << G4endl;
    step -> GetTrack() -> SetTrackStatus(fStopAndKill);
    fEventAction -> SetTrappedTrue();	// sets the event action flag as true.
  }

  // stores energy deposition. i.e. important information tracker
  // THIS WILL BE MOVED OVER TO TrackerSD AND TrackerHit!
  const DetectorConstruction* detectorConstruction =
        static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  G4double edepStep = step->GetTotalEnergyDeposit();

  // check if the volume we are in is one of the logical volumes we're interested in
  if(volume == (*detectorConstruction).scint_scintillator_log[0])
  {
    fEventAction -> AddEdep(edepStep, 0, 0);
  }
  if(volume == (*detectorConstruction).mwpc_container_log[0])
  {
    fEventAction -> AddEdep(edepStep, 1, 0);
  }
  if(volume == (*detectorConstruction).scint_scintillator_log[1])
  {
    fEventAction -> AddEdep(edepStep, 0, 1);
  }
  if(volume == (*detectorConstruction).mwpc_container_log[1])
  {
    fEventAction -> AddEdep(edepStep, 1, 1);
  }

}

