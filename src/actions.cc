#include <sstream>

#include "actions.hh"

//#include "nain4.hh"
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-volumes.hh"
#include "generator.hh"

#include "G4AnalysisManager.hh"

//bool fS1 = false;

n4::actions* create_actions(G4ThreeVector& init_pos, params params) {
  G4AnalysisManager* ana_man = G4AnalysisManager::Instance();
  double el_rad = params.el_rad;

  auto start_run_action = [&, ana_man] (const G4Run* run){
    ana_man->OpenFile();
    ana_man->CreateNtuple("Hits","Hits");
    ana_man->CreateNtupleIColumn("Event");
    ana_man->CreateNtupleDColumn("Xi");
    ana_man->CreateNtupleDColumn("Yi");
    ana_man->CreateNtupleDColumn("Zi");
    ana_man->CreateNtupleDColumn("Xf");
    ana_man->CreateNtupleDColumn("Yf");
    ana_man->CreateNtupleDColumn("Zf");
    ana_man->CreateNtupleIColumn("SensorID");
    ana_man->FinishNtuple(0);
    };
  auto end_run_action = [&] (const G4Run* run){
    G4AnalysisManager* ana_man = G4AnalysisManager::Instance();
    ana_man->Write();
	ana_man->CloseFile();
  };
    if (params.fCylinder){el_rad = sqrt((params.cone_rad*params.cone_rad + params.el_rad*params.el_rad + params.cone_rad*params.el_rad)/3);}
    double emission_rad = el_rad;
    double emission_h   = -(params.cone_h/2 + params.buffer_gap + params.neck_gap);
    if(params.fCylinder){emission_h = params.cone_h/2 - params.sensor_gap - params.buffer_gap;}

    double angle = atan(params.cone_h/(params.cone_rad-el_rad));
    double h_aux = tan(angle) * el_rad;
    
    if (params.fS1){
    std::cout << "h = " << emission_h << std::endl;
    emission_h   = emission_h + params.buffer_gap + params.neck_gap + params.source_h;
    std::cout << "h = " << emission_h << std::endl;
    emission_rad = (params.source_h + h_aux)/tan(angle);
    std::cout << "r = " << emission_rad << std::endl;
    std::cout << "el r = " << el_rad << std::endl;
    //getchar();
    }    

    auto generate = [&init_pos, emission_rad, emission_h](auto event){el_generator(event, init_pos, emission_rad, emission_h); };
  return (new n4::        actions{          generate})
  -> set((new n4::     run_action{                  }) -> begin(start_run_action) -> end(end_run_action));
}