#include <n4-all.hh>
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-volumes.hh"

#include "geometry.hh"
#include "generator.hh"
#include "actions.hh"

#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include <FTFP_BERT.hh>         // our choice of physics list

#include <G4EmStandardPhysics_option4.hh>
#include <G4OpticalPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>

#include <G4ThreeVector.hh>
#include <cstdlib>


int main(int argc, char** argv)
{
	G4ThreeVector init_pos = G4ThreeVector(0, 0, 0);
	
	G4int verbosity=0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list ->  ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    physics_list -> RegisterPhysics(new G4RadioactiveDecayPhysics);

	n4::run_manager::create()	
	.ui("COLINA", argc, argv)
    .macro_path("macs")
    //.apply_command("/my/straw_radius 0.5 m")
    //.apply_early_macro("early-hard-wired.mac")
	.apply_cli_early()
	.physics(physics_list)
	.geometry([&] {return colina_detector(init_pos);})
  	.actions(create_actions(init_pos))
    //.apply_command("/my/particle e-")
    //.apply_late_macro("late-hard-wired.mac")
	.apply_cli_late() // CLI?
	.run();
  	//n4::ui(argc, argv);
}
