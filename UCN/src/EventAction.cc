#include "EventAction.hh"
#include "Run.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

EventAction::EventAction()
: G4UserEventAction(),
  fEdep(0.)
{}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
}


void EventAction::EndOfEventAction(const G4Event*)
{
  // accumulate statistics in B1Run
  Run* run
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddEdep(fEdep);
}
