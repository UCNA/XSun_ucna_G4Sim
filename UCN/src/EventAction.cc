#include "EventAction.hh"
#include "Run.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

EventAction::EventAction()
: G4UserEventAction(),
  fEdep(0.)
{}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fEdep = 0.;	// reset value for B1Example.

  if((evt->GetEventID())%100 == 0)
  {
    G4cout << "/n -------------- Begin of event: " << evt->GetEventID() << G4endl;
  }
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  // Here Michael uses his analysis manager class
//  gAnalysisManager -> FillTrackerData(evt);
//  gAnalysisManager -> FillEventTree();
  // Need to examine Analysis Manager class and check what these do/ recreate it.




  // accumulate statistics in B1Run
  Run* run	// All this stuff related to B1Example.
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddEdep(fEdep);
}
