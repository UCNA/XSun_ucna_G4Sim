#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <G4Event.hh>

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void EndOfEventAction(const G4Event* evt);

    void AddEdep(G4double edep, int typeFlag, int locFlag);

  private:
    G4double  fEdep_East_Scint;
    G4double  fEdep_East_MWPC;
    G4double  fEdep_West_Scint;
    G4double  fEdep_West_MWPC;

};

#endif


