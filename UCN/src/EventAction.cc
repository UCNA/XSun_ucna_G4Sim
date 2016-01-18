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
: G4UserEventAction(), fStartTime(0)
{
  fMyDetectorConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  // initialize all values to -1 since that doesn't correspond to any actual SD

  for(int i = 0; i < fNbSDs; i++)	// see comment in EventAction.hh on fNbSDs
  {
    fHitsCollectionIDs[i] = -1;
  }

  fScintEast_index = -1;
  fScintWest_index = -1;
  fActiveWireVolEast_index = -1;
  fActiveWireVolWest_index = -1;
}


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
  if(fHitsCollectionIDs[0] == -1)
  {
    for(int i = 0; i < fNbSDs; i++)
    {
      fHitsCollectionIDs[i] = G4SDManager::GetSDMpointer()->GetCollectionID((*fMyDetectorConstruction).fHCNamesArray[i]);
      G4cout << "HC ID: " << fHitsCollectionIDs[i] << ", belongs to HC name: " << (*fMyDetectorConstruction).fHCNamesArray[i] << G4endl;
      // these are the special ones we want to record extra info from.
      if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_scint_EAST")
      {
        fScintEast_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_scint_WEST")
      {
        fScintWest_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_wireVol_EAST")
      {
        fActiveWireVolEast_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_wireVol_WEST")
      {
        fActiveWireVolWest_index = i;
      }
    }
  }

  TrackerHitsCollection* SD_totalHC[fNbSDs];
  TrackerHit* SD_hitsFullStep[fNbSDs];
  TrackerHit* SD_hitsTrack[fNbSDs];
  G4double SD_edepFullStep[fNbSDs];
  G4double SD_edepTrack[fNbSDs];

  for(int j = 0; j < fNbSDs; j++)
  {
    SD_edepFullStep[j] = 0;
    SD_edepTrack[j] = 0;
  }

  for(int i = 0; i < 2; i++)
  {
    SD_totalHC[i] = GetHitsCollection(fHitsCollectionIDs[i], evt);
    SD_hitsFullStep[i] = (*SD_totalHC[i])[0];
    SD_edepFullStep[i] = SD_hitsFullStep[i] -> GetEdep();

    SD_hitsTrack[i] = (*SD_totalHC[i])[1];
    SD_edepTrack[i] += SD_hitsTrack[i] -> GetEdepTrack();


    G4cout << "Index " << i << " has HC of size = " << SD_totalHC[i] -> GetSize() << G4endl;

  }




/*  for(int i = 0; i < fNbSDs; i++)
  {
    if( (i != 2) && (i != 3)) {
    SD_totalHC[i] = GetHitsCollection(fHitsCollectionIDs[i], evt);
*/

  /* In the old code, you made a TrackerHit item for each SD and only incremented the ones
     that we went inside of. This was because .insert(new TrackerHit()) happened in initialize.
     Now, we only .insert a new TrackerHit during the ProcessHits method. Hence, if we do not
     enter some of the other SD's, there is no default TrackerHit with zeros everywhere.
     Thus we need some check statement to make sure we don't try to access a null pointer.

     Furthermore, you are currently saving a TrackerHit object in the fHitsCollection for each
     track made inside the SD's. This means you have some arbitrary number of TrackerHits in each fHitCollection.
     You want to access the TrackerHit object in the LAST entry of fHitsCollection. This is the one
     that has the final accumulation variables! However, this is non-trivial to do. Also, perhaps
     storing so many TrackerHits is a time-intensive operation. We could just use a single TrackerHit
     as an event-by-event accumulator, which is what the Basic Example does.

     As you wanted to test and unfortunately confirmed, there's no guarantee the last entry stores all
     the info. It seems like each TrackerHit stores it's own energy and all of it would get summed up later.
     That isn't how Michael's code seems to be laid out. But it is the reality.

     NOTE: you got the ordering mixed up. A HitsCollection item doesn't exist for i = 2, 3.
     So the whole if-statement needs to encompass everything.
     NOTE: OK so you don't understand what is happening with this Seg Fault. Instead of trying to fix it,
     you may choose to redo everything anyway so don't worry too much about it.

     Soon, you need to choose which way we'll record info. Is it easier to do it like Basic Example?
     What information do we lose if that is how we choose to do it?

  */
/*
    G4cout << "Index " << i << " has HC of size = " << SD_totalHC[i] -> GetSize() << G4endl;

    SD_hitsFullStep[i] = (*SD_totalHC[i])[0];   // 0th entry is the hit item in HC we are using to record everything
    SD_edepFullStep[i] = SD_hitsFullStep[i] -> GetEdep();

    for(int t = 1; t < SD_totalHC[i] -> GetSize(); t++)
    {
    SD_hitsTrack[i] = (*SD_totalHC[i])[t];

    SD_edepTrack[i] = SD_edepTrack[i] + (SD_hitsTrack[i] -> GetEdepTrack());
    }

    }
  }

*/


  // print out of class member variables due to stepping action
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << fTrapped << "\t"
	  << compTime << "\n"
	  << "\t" << SD_edepFullStep[fScintEast_index]/keV << "\t"
	  << "\t" << SD_edepFullStep[fActiveWireVolEast_index]/keV << "\t"
	  << "\t" << SD_edepFullStep[fScintWest_index]/keV << "\t"
	  << "\t" << SD_edepFullStep[fActiveWireVolWest_index]/keV << "\n"
          << "\t" << SD_edepTrack[fScintEast_index]/keV << "\t"
          << "\t" << SD_edepTrack[fActiveWireVolEast_index]/keV << "\t"
          << "\t" << SD_edepTrack[fScintWest_index]/keV << "\t"
          << "\t" << SD_edepTrack[fActiveWireVolWest_index]/keV << "\n";
  outfile.close();
}

// locFlag = 0 -> EAST
// 	     1 -> WEST
// typeFlag = 0 -> Scint
//	      1 -> MWPC

