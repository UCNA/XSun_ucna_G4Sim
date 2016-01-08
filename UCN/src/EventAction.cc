#include "EventAction.hh"
#include "RunAction.hh"
#include "TrackerHit.hh"
#include "TrackerSD.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
using   namespace       std;

#define	OUTPUT_FILE	"TrackerTestOutput.txt"

EventAction::EventAction()
: G4UserEventAction(), fStartTime(0),
  fEastScintHCID(-1), fWestScintHCID(-1),
  fEastwireVolHCID(-1), fWestwireVolHCID(-1)
{}


EventAction::~EventAction()
{}

TrackerHitsCollection* EventAction::GetHitsCollection(int hcID, const G4Event* event) const
{
  TrackerHitsCollection* hitsCollection = static_cast<TrackerHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  return hitsCollection;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fTrapped = false;
  fEdep_East_Scint = 0;	// Ensuring these values are reset.
  fEdep_West_Scint = 0;
  fEdep_East_MWPC = 0;
  fEdep_West_MWPC = 0;
  fStartTime = 0;

  // Sets the start of the C++ 'clock' used for tracking trapped ptcl's
  fStartTime = clock();

  if((evt->GetEventID())%1000 == 0)
  {
    G4cout << "\n -------------- Begin of event: " << evt->GetEventID() << "\n" << G4endl;
  }
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  if(fTrapped == true)
  {
    G4cout << "Event " << evt->GetEventID() << " trapped. Incrementing kill count." << G4endl;
    // some code that M. Mendenhall came up with that allows me to access my own RunAction class
    ((RunAction*)G4RunManager::GetRunManager()->GetUserRunAction()) -> IncrementKillCount();
  }

  // store how much time the event ran for
  clock_t timeOfEvent = clock() - fStartTime;
  double compTime = ((double)timeOfEvent)/CLOCKS_PER_SEC;

  // get hits collection IDs only once in an overall run though.
  if(fEastScintHCID == -1)
  {
    fEastScintHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HC_scint_0");
    fWestScintHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HC_scint_1");
    fEastwireVolHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HC_wireVol_0");
    fWestwireVolHCID = G4SDManager::GetSDMpointer()->GetCollectionID("HC_wireVol_1");
  }

  TrackerHitsCollection* EScint_totHC = GetHitsCollection(fEastScintHCID, evt);
  TrackerHitsCollection* WScint_totHC = GetHitsCollection(fWestScintHCID, evt);
  TrackerHitsCollection* EwireVol_totHC = GetHitsCollection(fEastwireVolHCID, evt);
  TrackerHitsCollection* WwireVol_totHC = GetHitsCollection(fWestwireVolHCID, evt);

  // since my HC was made with 2 hits (1 for current event, 1 for total accumulation)
  // accessing total accumulation  should be # of entries -1, or equivalently, index 1 in the array
  // I've gotten rid of the hitTotal TrackerHit so now the only entry for each HC should be 0
  TrackerHit* EScint_hit = (*EScint_totHC)[0];
  TrackerHit* WScint_hit = (*WScint_totHC)[WScint_totHC->entries()-1];
  TrackerHit* EwireVol_hit = (*EwireVol_totHC)[EwireVol_totHC->entries()-1];
  TrackerHit* WwireVol_hit = (*WwireVol_totHC)[WwireVol_totHC->entries()-1];

  G4double edep_e_scint, edep_w_scint, edep_e_wireVol, edep_w_wireVol;
  edep_e_scint = EScint_hit->GetEdep();	// variables for printing out.
  edep_w_scint = WScint_hit->GetEdep();
  edep_e_wireVol = EwireVol_hit->GetEdep();
  edep_w_wireVol = WwireVol_hit->GetEdep();

  // print out of class member variables due to stepping action
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << fTrapped << "\t"
	  << compTime << "\n"
	  << fEdep_East_Scint/keV << "\t" << edep_e_scint/keV << "\n"
	  << fEdep_East_MWPC/keV << "\t" << edep_e_wireVol/keV << "\n"
	  << fEdep_West_Scint/keV << "\t" << edep_w_scint/keV << "\n"
	  << fEdep_West_MWPC/keV << "\t" << edep_w_wireVol/keV << "\n";
  outfile.close();
}
// locFlag = 0 -> EAST
// 	     1 -> WEST
// typeFlag = 0 -> Scint
//	      1 -> MWPC
void EventAction::AddEdep(G4double edep, int typeFlag, int locFlag)
{
  if(locFlag == 0)
  {
    if(typeFlag == 0)
    {
      fEdep_East_Scint = fEdep_East_Scint + edep;
    }
    else if(typeFlag == 1)
    {
      fEdep_East_MWPC = fEdep_East_MWPC + edep;
    }
  }
  else if(locFlag == 1)
  {
    if(typeFlag == 0)
    {
      fEdep_West_Scint = fEdep_West_Scint + edep;
    }
    else if(typeFlag == 1)
    {
      fEdep_West_MWPC = fEdep_West_MWPC + edep;
    }
  }
}
