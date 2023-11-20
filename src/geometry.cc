#include "geometry.hh"

#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include "G4AnalysisManager.hh"
#include <G4ThreeVector.hh>

#include "material_props.hh"

G4PVPlacement* colina_detector(G4ThreeVector& init_pos) {

  bool fCylinder     = false;

  bool fR_Holder     = true;
  bool fR_Inner      = true;
  bool fR_Cap        = true;
  bool fCathodeRing1 = false;
  int n_sensors_ring1 = 1;
  bool fCathodeRing2 = false;
  int n_sensors_ring2 = 6;//8 and 3 for 12, 10 and 5 for 16

  bool fSpecular_Mirror = false;
//19.17, 38.342
  double form_factor  = 1;
  double sphere_rad   = 24.44 *cm;//15.11*cm;
  double sphere_h     = form_factor * sphere_rad;//24.48*cm;
  double el_rad       = 32.*mm;//32. *mm; 
//  double el_rad       = fCylinder ? sphere_rad : 32.*mm;//32. *mm; 
  if (fCylinder){sphere_rad = sqrt((sphere_rad*sphere_rad + el_rad*el_rad + sphere_rad*el_rad)/3); el_rad = sphere_rad;}

  double vessel_thick = 5  *mm;
  double ptfe_thick   = 5  *mm;
  double grid_thickness    = 0.1*mm;
  double grid_transparency = 0.956;
  double neck_gap     = fCylinder ? 0*mm : 10 *mm;
  double el_gap       = 5  *mm + neck_gap + grid_thickness;
  double buffer_gap   = 5    *mm;
  double sensor_gap   = 0.5  *cm;

  bool f12699        = false;
  bool fCenterSquare = fCylinder ? false : true;
  //bool 

  double sensor_width = 14.8*mm;//ptfe_thick;
  double window_rad   = fCylinder ? sphere_rad : el_rad;
  double window_side  = f12699 ? 56*  mm : 64.5*mm;
  double sensor_side  = f12699 ? 48.5*mm : window_side;
  double sensor_rad   = fCenterSquare ? sqrt(2*(window_side/2)*(window_side/2)) : window_rad ; //32 mm for PMT, 8 for APD, 32 for SiPM

  auto liquid = lxe_with_props();
  auto grid   = fakegrid_with_props(liquid, grid_transparency, grid_thickness);
  auto quartz = quartz_with_props();
  auto air    = n4::material("G4_AIR");
  auto steel  = steel_with_props();//n4::material("G4_STAINLESS-STEEL");
  auto ptfe   = n4::material("G4_TEFLON");
  auto world  = n4::box("world").cube(0.5*m).volume(air);
  //double refl_ptfe = 0.8;
  double refl_ptfe = 0.98;//0.98;
  G4OpticalSurface* sphere_reflector =
    new G4OpticalSurface("sphere_reflector", unified, polished, dielectric_metal);
  G4OpticalSurface* ptfe_reflector =
    new G4OpticalSurface("ptfe_reflector", unified, polished, dielectric_metal);
  ptfe_reflector->SetMaterialPropertiesTable(PTFE_props(refl_ptfe));
  G4OpticalSurface* sensor_reflector =
    new G4OpticalSurface("sensor_reflector", unified, polished, dielectric_metal);
  double refl = f12699 ? 0 : 0.2;
  //refl = 0.2;// CYLINDER
  double trans = 1-refl;
  sensor_reflector->SetMaterialPropertiesTable(Sensor_mirror_props(refl, trans));

  if(fSpecular_Mirror) {sphere_reflector->SetMaterialPropertiesTable(Al2O3_mirror_props(0.825));} //0.825 for EXO Cu-Al
  else {sphere_reflector->SetMaterialPropertiesTable(PTFE_props(refl_ptfe));}

  auto vessel = n4::cons("vessel")  .r1(el_rad + vessel_thick)
                                    .r2(sphere_rad + vessel_thick)
                                    .z(sphere_h)
                                    .place(steel)
                                    .in(world)
                                    //.rotate_z(180*deg)
                                    .at(0, 0, 0)
                                    .now();


  double holder_w = el_gap + buffer_gap + 2 * grid_thickness;
  double holder_y = -holder_w/2 - sphere_h/2;
  double holder_rad = el_rad;
  auto   holder   = n4::tubs("holder").r(holder_rad + ptfe_thick)
                                      .r_inner(holder_rad)
                                      .z(holder_w)
                                      .place(ptfe)
                                      .in(world)
                                      //.rotate_x(90*deg)
                                      .at_z(holder_y)
                                      .check_overlaps()
                                      .now();



  double sensor_holder_w = sensor_gap + ptfe_thick;   
  double sensor_holder_pos = holder_y - holder_w/2 - sensor_holder_w/2;
  auto   sensor_holder   = n4::cons("sensor_holder").r1(sensor_rad+ptfe_thick)//.r1(sensor_rad + 2*ptfe_thick)
                                             .r2(sensor_rad+ptfe_thick)
                                             //.r_delta(ptfe_thick)
                                             .z(sensor_holder_w)
                                             .place(ptfe)
                                             .in(world)
                                             //.rotate_x(90*deg)
                                             .at_z(sensor_holder_pos)
                                             .check_overlaps()
                                             .now();

  double cap_w = vessel_thick;
  double cap_y = cap_w/2 + sphere_h/2;
  auto cap = n4::tubs("cap").r(sphere_rad + vessel_thick)
                            .z(cap_w)
                            .place(steel)
                            //.place(ptfe)
                            .in(world)
                            //.rotate_x(90*deg)
                            .at_z(cap_y)
                            .check_overlaps()
                            .now();

  auto target = n4::cons("target")  .r1(el_rad)
                                    .r2(sphere_rad)
                                    .z(sphere_h)
                                    .place(liquid)
                                    .in(vessel)
                                    //.rotate_z(180*deg)
                                    .at(0, 0, 0)
                                    .check_overlaps()
                                    .now();
                    


  double el_pos = - el_gap/2 - sphere_h/2 ;
  auto el_phys = n4::tubs("el_gap").r(el_rad)
                                   .z(el_gap)
                                   .place(liquid)
                                   .in(world)
                                   //.rotate_x(90*deg)
                                   .at_z(el_pos)
                                   .check_overlaps()
                                   .now();
                                   
  double el_grid_pos = -el_gap/2 + buffer_gap;
  auto el_grid = n4::tubs("el_grid").r(el_rad)
                                   .z(grid_thickness)
                                   .place(grid)
                                   .in(fCylinder ? target : el_phys)
                                   //.rotate_x(90*deg)
                                   .at_z(fCylinder ? sphere_h/2 - sensor_gap : el_grid_pos)
                                   .check_overlaps()
                                   .now();

  double anode_pos = el_pos - el_gap/2 - grid_thickness/2;
  auto anode_grid = n4::tubs("anode_grid").r(el_rad)
                                      .z(grid_thickness)
                                      .place(grid)
                                      .in(world)
                                      //.rotate_x(90*deg)
                                      .at_z(anode_pos)
                                      .check_overlaps()
                                      .now();


  double buffer_pos = anode_pos - grid_thickness/2 - buffer_gap/2;
  auto buffer_phys = n4::tubs("buffer_gap").r(el_rad)
                                           .z(buffer_gap)
                                           .place(liquid)
                                           .in(world)
                                           //.rotate_x(90*deg)
                                           .at_z(buffer_pos)
                                           .check_overlaps()
                                           .now();

  double shield_pos = buffer_pos - buffer_gap/2 - grid_thickness/2;
  auto shield_grid  = n4::tubs("shield_grid").r(el_rad)
                                        .z(grid_thickness)
                                        .place(grid)
                                        .in(world)
                                        //.rotate_x(90*deg)
                                        .at_z(shield_pos)
                                        .check_overlaps()
                                        .now();

  
  double sensor_gap_pos = sensor_holder_w/2 - sensor_gap/2;
  auto sg_phys = n4::cons("sensor_gap").r1(holder_rad)//r1(sensor_rad)
                                       .r2(holder_rad)
                                       .z(sensor_gap)
                                       .place(liquid)
                                       .in(sensor_holder)
                                       .at_z(sensor_gap_pos)
                                       .check_overlaps();

 auto sensor_gap_phys = sg_phys.now();
 double window_width = f12699 ? 2.5*mm : 0.0001*mm;
 double window_pos   = sensor_gap_pos - sensor_gap/2 - window_width/2; 
 auto sd = sensitive(init_pos, 0, "Photosensor_0");
 double sensor_pos = sensor_gap_pos - sensor_gap/2 - window_width/2 * 2 - sensor_width/2 ;

 if (fCenterSquare){
 new G4LogicalBorderSurface("reflectorSensor_0"  , sensor_gap_phys, n4::box("Window_0").xy(window_side)
                                                                                       .z(window_width)
                                                                                       .place(liquid)
                                                                                       .in(sensor_holder)
                                                                                       .at_z(window_pos)
                                                                                       .check_overlaps()
                                                                                       .now(), 
                                                                                       sensor_reflector);
  n4::box("Photosensor_0").xy(sensor_side)
                          .z(sensor_width)
                          .sensitive(sd)
                          .place(liquid)
                          .in(sensor_holder)
                          .at_z(sensor_pos)
                          .check_overlaps()
                          .now();
  }
  else{
  new G4LogicalBorderSurface("reflectorSensor_0"  , sensor_gap_phys, n4::tubs("Window_0").r(window_rad)
                                                                                       .z(window_width)
                                                                                       .place(quartz)
                                                                                       .in(sensor_holder)
                                                                                       .at_z(window_pos)
                                                                                       .check_overlaps()
                                                                                       .now(), 
                                                                                       fCylinder ? ptfe_reflector : sensor_reflector);//sensor_reflector);
  n4::tubs("Photosensor_0").r(window_rad)
                          .z(sensor_width)
                          .sensitive(sd)
                          .place(liquid)
                          .in(sensor_holder)
                          .at_z(sensor_pos)
                          .check_overlaps()
                          .now();
  }

  if(fCathodeRing1){
    double rad1         = 0*cm; //7 for 5, 6 for 3
    double angle1       = 360/n_sensors_ring1 *deg;
    double x;
    double y;
    for(int i=1; i<=n_sensors_ring1;i++)
    { 
      window_pos = -cap_w/2 + window_width/2;
      sensor_pos = window_pos + window_width/2 + sensor_width/2;
      x = rad1*cos(angle1 * (i-1));
      y = rad1*sin(angle1 * (i-1));
      new G4LogicalBorderSurface("reflectorSensor_"+ std::to_string(i), target, 
                                  n4::box("Window_"+ std::to_string(i)).xy(window_side)
                                                       .z(window_width)
                                                       .place(liquid)
                                                       .in(cap)
                                                       .at(x, y, window_pos)
                                                       //.rotate_z(angle1 * (i-1))
                                                       .check_overlaps()
                                                       .now(), 
                                                       sensor_reflector);
      n4::box("Photosensor_"+ std::to_string(i)).xy(sensor_side)
                                                .z(sensor_width)
                                                .sensitive(sensitive(init_pos, i, "Photosensor_"+ std::to_string(i)))
                                                .place(liquid)
                                                .in(cap)
                                                .at(x, y, sensor_pos)
                                                //.rotate_z(angle1 * (i-1))
                                                .check_overlaps()
                                                .now();
    }
  }

  if(fCathodeRing2){
  double rad2         = 10*cm;//20*cm; //14 for 8, 15 for 10
  double angle2       = 360/n_sensors_ring2 *deg;
  double angle0       = 0;angle2/2;
  double x;
  double y;
  for(int i=1+n_sensors_ring1; i<=n_sensors_ring1+n_sensors_ring2;i++)
  { 
      window_pos = -cap_w/2 + window_width/2;
      sensor_pos = window_pos + window_width/2 + sensor_width/2;
      x = rad2*cos(angle0 + angle2 * (i-1));
      y = rad2*sin(angle0 + angle2 * (i-1));
      new G4LogicalBorderSurface("reflectorSensor_"+ std::to_string(i), target, 
                                  n4::box("Window_"+ std::to_string(i)).xy(window_side)
                                                       .z(window_width)
                                                       .place(liquid)
                                                       .in(cap)
                                                       .at(x, y, window_pos)
                                                       //.rotate_z(angle2 * (i-1))
                                                       .check_overlaps()
                                                       .now(), 
                                                       sensor_reflector);
      n4::box("Photosensor_"+ std::to_string(i)).xy(sensor_side)
                                                .z(sensor_width)
                                                .sensitive(sensitive(init_pos, i, "Photosensor_"+ std::to_string(i)))
                                                .place(liquid)
                                                .in(cap)
                                                .at(x, y, sensor_pos)
                                                //.rotate_z(angle2 * (i-1))
                                                .check_overlaps()
                                                .now();
    }
  }
  if (fR_Inner ) {new G4LogicalBorderSurface("reflectorInner1"  , target, vessel, ptfe_reflector);}
  if (fR_Cap   ) {new G4LogicalBorderSurface("reflectorCap1"    , target, cap   , sphere_reflector);}
//  if (fR_Cap   ) {new G4LogicalBorderSurface("reflectorCap1"    , target, cap   , sphere_reflector);}
  if (fR_Holder) {new G4LogicalBorderSurface("reflectorHoldA"  , sensor_gap_phys, sensor_holder, ptfe_reflector);
                  new G4LogicalBorderSurface("reflectorHoldB"  , target         , holder       , ptfe_reflector);
                  new G4LogicalBorderSurface("reflectorHoldC"  , el_phys        , holder       , ptfe_reflector);
                  new G4LogicalBorderSurface("reflectorHoldD"  , buffer_phys    , holder       , ptfe_reflector);}

  return n4::place(world).now();
}

n4::sensitive_detector* sensitive(G4ThreeVector& init_pos, unsigned int sensorID, G4String sensor_name) {

    n4::sensitive_detector::process_hits_fn process_hits = [&init_pos, sensorID] (auto step) {
	step->GetTrack()->SetTrackStatus(fStopAndKill);
	G4StepPoint *preStepPoint  = step->GetPreStepPoint();  //ENTER DETECTOR
	G4ThreeVector posPhoton    = preStepPoint->GetPosition();		
	G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4AnalysisManager * ana_man = G4AnalysisManager::Instance();
	ana_man->FillNtupleIColumn(0, evt         );
	ana_man->FillNtupleDColumn(1, init_pos [0]);
	ana_man->FillNtupleDColumn(2, init_pos [1]);
	ana_man->FillNtupleDColumn(3, init_pos [2]);
	ana_man->FillNtupleDColumn(4, posPhoton[0]);
	ana_man->FillNtupleDColumn(5, posPhoton[1]);
	ana_man->FillNtupleDColumn(6, posPhoton[2]);
	ana_man->FillNtupleIColumn(7,     sensorID);
	ana_man->AddNtupleRow(0);	
    return true;
  };

  return (new n4::sensitive_detector{sensor_name, process_hits});
}