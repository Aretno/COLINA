#include <sstream>

#include "actions.hh"

//#include "nain4.hh"
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-volumes.hh"
#include "generator.hh"

#include "G4AnalysisManager.hh"

int  fMaxCasesS1 = 6;


bool fCylinder     = false;
bool fS1 = false;
int  pos_S1 = 6;
double sphere_rad   = 24.44*cm;
double form_factor  = 1;
double sphere_h     = sphere_rad * form_factor;
double el_rad       = 32. *mm;
double buffer_gap   = 5  *mm;
double neck_gap     = 10 *mm;
double sensor_gap   = 0.5  *cm;

n4::actions* create_actions(G4ThreeVector& init_pos) {
  auto start_run_action = [&] (const G4Run* run){
    G4AnalysisManager* ana_man = G4AnalysisManager::Instance();

    G4String filename =  fCylinder ? "Cylinder_1Sensor_SiPM_PTFE_R" : "Cone_1Sensor_SiPM_PTFE_R";
    //G4String filename =  ;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << sphere_rad/10;
    filename += stream.str() + "cm_FF";
    stream.str(std::string());

    stream << std::fixed << std::setprecision(2) << form_factor;
    filename += stream.str() + "_EL";
    stream.str(std::string());

    stream << std::fixed << std::setprecision(0) << el_rad;
    filename += stream.str() + "mm_dist";
    stream.str(std::string()); 

    stream << std::fixed << std::setprecision(0) << sensor_gap;
    filename += stream.str() + "mm_S";
    stream.str(std::string()); 

    filename += fS1 ? std::to_string(1) : std::to_string(2);

    if(fS1){filename += "_Case" + std::to_string(pos_S1) + ".root";}
    //else {filename += "neck0.root";}
    else {filename += ".root";}

    filename = "test2.root";
    
    ana_man->SetFileName(filename);

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
//    auto generate = [&init_pos](auto event){el_generator(event, init_pos, 40.5*mm, -10.5*cm); };
    if (fCylinder){el_rad = sqrt((sphere_rad*sphere_rad + el_rad*el_rad + sphere_rad*el_rad)/3);}
    double emission_rad = el_rad;
    double emission_h   = -(sphere_h/2 + buffer_gap+ neck_gap);
    if(fCylinder){emission_h = sphere_h/2 - sensor_gap - buffer_gap;}

    //emission_rad = sphere_h/2;
    //emission_h   = 0;

    // For cylinder test
    //el_rad   = 18.94*cm;//15.11 * cm;
    //sphere_h = 31.08*cm;//24.48*cm;
//
    //emission_rad = el_rad;
    //emission_h   = -(sphere_h/2 + buffer_gap+ neck_gap); //S1 cyl

    //emission_h   = sphere_h/2 - 5*mm - buffer_gap; // S2 cyl
    //bool fS1 = false;
    double angle = atan(sphere_h/(sphere_rad-el_rad));
    double h_aux = tan(angle) * el_rad;
    double jump_h = sphere_h/fMaxCasesS1;
    double displacement = (pos_S1 == 0) ? 1 : pos_S1*jump_h;
    
    if (fS1){
      emission_h   = sphere_h/2 - displacement;
      emission_rad = (sphere_h - displacement + h_aux)/tan(angle);
    }    

    auto generate = [&init_pos, emission_rad, emission_h](auto event){el_generator(event, init_pos, emission_rad, emission_h); };
    //auto generate = [&init_pos](auto event){sphere_generator(event, init_pos); };
  return (new n4::        actions{          generate})
  -> set((new n4::     run_action{                  }) -> begin(start_run_action) -> end(end_run_action));
}