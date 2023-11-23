#include <n4-all.hh>
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-volumes.hh"

#include "geometry.hh"
#include "generator.hh"
#include "actions.hh"
#include "colina.hh"

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
	params params;
	load_messenger(params);
	G4int verbosity=0;
    auto physics_list = new FTFP_BERT{verbosity};
    physics_list ->  ReplacePhysics(new G4EmStandardPhysics_option4());
    physics_list -> RegisterPhysics(new G4OpticalPhysics{});
    physics_list -> RegisterPhysics(new G4RadioactiveDecayPhysics);

	n4::run_manager::create()	
	.ui("COLINA", argc, argv)
    .macro_path("macs")
	.apply_cli_early()
	.physics(physics_list)
	.geometry([&] {return colina_detector(init_pos, params);})
  	.actions(create_actions(init_pos, params))
	.apply_cli_late()
	.run();
  	//n4::ui(argc, argv);
}

void load_messenger(params& params)
{
	//std::cout << "HEHE" << getchar() << std::endl;
	auto geom_messenger = new G4GenericMessenger{nullptr, "/geom/", "Geometry messenger"};
	geom_messenger -> DeclareProperty        ("cylinder"    , params.fCylinder       );
	geom_messenger -> DeclareProperty        ("neck_refl"   , params.fR_Holder       );
	geom_messenger -> DeclareProperty        ("inner_refl"  , params.fR_Inner        );
	geom_messenger -> DeclareProperty        ("cap_refl"    , params.fR_Cap          );
	geom_messenger -> DeclareProperty        ("cap_specular", params.fSpecular_Mirror);

	geom_messenger -> DeclareProperty        ("cathode_ring_1", params.fCathodeRing1  );
	geom_messenger -> DeclareProperty        ("cathode_ring_2", params.fCathodeRing2  );
	geom_messenger -> DeclareProperty        ("n_sensors_1"   , params.n_sensors_ring1);
	geom_messenger -> DeclareProperty        ("n_sensors_2"   , params.n_sensors_ring2);

	geom_messenger -> DeclarePropertyWithUnit("cone_r"     ,"m", params.cone_rad  );
	geom_messenger -> DeclarePropertyWithUnit("cone_h"     ,"m", params.cone_h    );
	geom_messenger -> DeclarePropertyWithUnit("el_r"       ,"m", params.el_rad    );
	geom_messenger -> DeclarePropertyWithUnit("neck_gap"   ,"m", params.neck_gap  );
	geom_messenger -> DeclarePropertyWithUnit("el_gap"     ,"m", params.el_gap    );
	geom_messenger -> DeclarePropertyWithUnit("buffer_gap" ,"m", params.buffer_gap);
	geom_messenger -> DeclarePropertyWithUnit("sensor_gap" ,"m", params.sensor_gap);
	geom_messenger -> DeclarePropertyWithUnit("sensor_side","m", params.sensor_gap);

	auto gen_messenger = new G4GenericMessenger{nullptr, "/gen/", "Generator messenger"};
	gen_messenger -> DeclareProperty        ("fS1"                       , params.fS1);
	gen_messenger -> DeclarePropertyWithUnit("source_h", "m", params.source_h);

}