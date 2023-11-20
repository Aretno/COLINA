#include "material_props.hh"

//#include "nain4.hh"
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-utils.hh"
//#include "n4-volumes.hh"

G4MaterialPropertiesTable* PTFE_props()
{
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
  // REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
//  std::vector<G4double> REFLECTIVITY = {0.8, 0.8};
  std::vector<G4double> REFLECTIVITY = {0.98, 0.98};
  //std::vector<G4double> REFLECTIVITY = {1, 1};
  // REFLEXION BEHAVIOR
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  std::vector<G4double> specularlobe  = {0., 0.};
  // specular reflection about the average normal
  std::vector<G4double> specularspike = {0., 0.};
  // 180 degrees reflection.
  std::vector<G4double> backscatter   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  return n4::material_properties()
         .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
         .add("SPECULARLOBECONSTANT", ENERGIES, specularlobe)
         .add("SPECULARSPIKECONSTANT",ENERGIES, specularspike)
         .add("BACKSCATTERCONSTANT",  ENERGIES, backscatter)
         //.add("RINDEX", ENERGIES, {1.5, 1.5})
         .done();
}

G4MaterialPropertiesTable* PTFE_props(double refl)
{
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
  // REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
//  std::vector<G4double> REFLECTIVITY = {0.8, 0.8};
  std::vector<G4double> REFLECTIVITY = {refl, refl};
  //std::vector<G4double> REFLECTIVITY = {1, 1};
  // REFLEXION BEHAVIOR
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  std::vector<G4double> specularlobe  = {0., 0.};
  // specular reflection about the average normal
  std::vector<G4double> specularspike = {0., 0.};
  // 180 degrees reflection.
  std::vector<G4double> backscatter   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  return n4::material_properties()
         .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
         .add("SPECULARLOBECONSTANT", ENERGIES, specularlobe)
         .add("SPECULARSPIKECONSTANT",ENERGIES, specularspike)
         .add("BACKSCATTERCONSTANT",  ENERGIES, backscatter)
         //.add("RINDEX", ENERGIES, {1.5, 1.5})
         .done();
}

G4MaterialPropertiesTable* Al2O3_mirror_props(double refl)
{
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
  // REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
//  std::vector<G4double> REFLECTIVITY = {0.8, 0.8};
  std::vector<G4double> REFLECTIVITY = {refl, refl}; // 0.825 for EXO Al-Cu
  // REFLEXION BEHAVIOR
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  std::vector<G4double> specularlobe  = {0., 0.};
  // specular reflection about the average normal
  std::vector<G4double> specularspike = {1., 1.};
  // 180 degrees reflection.
  std::vector<G4double> backscatter   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  return n4::material_properties()
         .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
         .add("SPECULARLOBECONSTANT", ENERGIES, specularlobe)
         .add("SPECULARSPIKECONSTANT",ENERGIES, specularspike)
         .add("BACKSCATTERCONSTANT",  ENERGIES, backscatter)
         .done();
}

G4MaterialPropertiesTable* Sensor_mirror_props(double refl, double trans)
{
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};
  // REFLECTIVITY IN LXE (from https://link.springer.com/content/pdf/10.1140/epjc/s10052-020-7800-6.pdf)
//  std::vector<G4double> REFLECTIVITY = {0.8, 0.8};
  std::vector<G4double> REFLECTIVITY = {refl, refl}; // 0.825 for EXO Al-Cu
  std::vector<G4double> TRANSMITTANCE = {trans, trans}; // 0.825 for EXO Al-Cu
  // REFLEXION BEHAVIOR
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  std::vector<G4double> specularlobe  = {0., 0.};
  // specular reflection about the average normal
  std::vector<G4double> specularspike = {1., 1.};
  // 180 degrees reflection.
  std::vector<G4double> backscatter   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  return n4::material_properties()
         .add("REFLECTIVITY", ENERGIES, REFLECTIVITY)
         .add("TRANSMITTANCE", ENERGIES, TRANSMITTANCE)
         .add("SPECULARLOBECONSTANT", ENERGIES, specularlobe)
         .add("SPECULARSPIKECONSTANT",ENERGIES, specularspike)
         .add("BACKSCATTERCONSTANT",  ENERGIES, backscatter)
         .done();
}


G4MaterialPropertiesTable* LXe_props()
{
    /// The time constants are taken from E. Hogenbirk et al 2018 JINST 13 P10031
    const G4int ri_entries = 200;
    // This is the range where Cerenkov photons are produced by G4.
    // Since the recommendation is to use the PDE range of the photosensors used in the simulation,
    // we choose that of Hamamatsu S15779(ES1) arrays as a typical range.
    const G4double minE_n = 1.38 * eV; // corresponds to 900 nm, where Hamamatsu's pde go to zero
    const G4double maxE_n = 8.21 * eV; // corresponds to 151 nm, where Hamamatsu's pde go to zero
    // + above this value A(omega) starts to diverge
    G4double eWidth = (maxE_n - minE_n) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
      ri_energy.push_back(minE_n + i * eWidth);
    }

    G4double density = 2.85 * g/cm3;;

    std::vector<G4double> ri_index;
    for (G4int i=0; i<ri_entries; i++) {
      //G4cout << ri_energy[i] << "," << XenonRefractiveIndex(ri_energy[i], density) << G4endl;
      ri_index.push_back(XenonRefractiveIndex(ri_energy[i], density));
    }
    const G4int sc_entries = 500;
    const G4double minE_sc = 6.20625*eV;
    const G4double maxE_sc = 8.21*eV;
    eWidth = (maxE_sc - minE_sc) / sc_entries;

    std::vector<G4double> sc_energy;
    for (G4int j=0; j<sc_entries; j++){
      sc_energy.push_back(minE_sc + j * eWidth);
    }
    std::vector<G4double> intensity;
    for (G4int i=0; i<sc_entries; i++) {
      intensity.push_back(LXeScintillation(sc_energy[i]));
    }

    std::vector<G4double> abs_energy =
      {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> abs_length =
      {noAbsLength_, noAbsLength_};

    std::vector<G4double> rayleigh_energy =
      {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> rayleigh_length = {36.*cm, 36.*cm};

    return n4::material_properties()
           .add("RINDEX", ri_energy, ri_index)
           .add("SCINTILLATIONCOMPONENT1", sc_energy, intensity)
           .add("SCINTILLATIONCOMPONENT2", sc_energy, intensity)
           .add("SCINTILLATIONYIELD", 58708./MeV)
           .add("RESOLUTIONSCALE", 1)
           .add("SCINTILLATIONTIMECONSTANT1", 2.*ns)
           .add("SCINTILLATIONTIMECONSTANT2", 43.5*ns)
           .add("SCINTILLATIONYIELD1", .03)
           .add("SCINTILLATIONYIELD2", .97)
           //.add("ATTACHMENT", 1000.*ms)
           .add("ABSLENGTH", abs_energy, abs_length)
           .add("RAYLEIGH", rayleigh_energy, rayleigh_length)
           .done();
  }

G4MaterialPropertiesTable* steel_props(){
  std::vector<G4double> ENERGIES = {optPhotMinE_, optPhotMaxE_};

  return n4::material_properties()
            .copy_from(LXe_props(), {
            "RINDEX"})
           .add("ABSLENGTH", ENERGIES, {10*m, 10*m})
           .done();
}

//////////////////////////////////////////////////////////////////////// OK
G4MaterialPropertiesTable* quartz_props(){

  // REFRACTIVE INDEX
  // The range is chosen to be up to ~10.7 eV because Sellmeier's equation
  // for fused silica is valid only in that range
  const size_t ri_entries = 200;
  auto ri_energy = n4::linspace(optPhotMinE_, optPhotMaxE_, ri_entries);
  

  // The following values for the refractive index have been calculated
  // using Sellmeier's equation:
  //    n^2 - 1 = B_1 * \lambda^2 / (\lambda^2 - C_1) +
  //            + B_2 * \lambda^2 / (\lambda^2 - C_2) +
  //            + B_3 * \lambda^2 / (\lambda^2 - C_3),
  // with wavelength \lambda in micrometers and
  //    B_1 = 4.73E-1, B_2 = 6.31E-1, B_3 = 9.06E-1
  //    C_1 = 1.30E-2, C_2 = 4.13E-3, C_3 = 9.88E+1.

  auto ref_index_sellmeier = [] (auto e) {
    auto B_1 = 4.73e-1;
    auto B_2 = 6.31e-1;
    auto B_3 = 9.06e-1;
    auto C_1 = 1.30e-2;
    auto C_2 = 4.13e-3;
    auto C_3 = 9.88e+1;
    auto lambda  = CLHEP::h_Planck*CLHEP::c_light/e*1000; // in micron  
    auto lambda2 = std::pow(lambda, 2);
    auto n2 = 1 + lambda2 * ( B_1 / (lambda2 - C_1)
                            + B_2 / (lambda2 - C_2)
                            + B_3 / (lambda2 - C_3));
    G4double result;
    
    if (n2 >0){result = std::sqrt(n2);}
    else {result = 1e+11;}
        
    return result;
  };
  
  auto rIndex = n4::map<G4double>(ref_index_sellmeier, ri_energy);
  // ABSORPTION LENGTH
  auto abs_energy = n4::scale_by(eV, {
    optPhotMinE_ / eV,
              6.46499, 6.54000, 6.59490, 6.64000, 6.72714, 6.73828, 6.75000,
              6.82104, 6.86000, 6.88000, 6.89000, 7.00000, 7.01000, 7.01797,
              7.05000, 7.08000, 7.08482, 7.30000, 7.36000, 7.40000, 7.48000,
              7.52000, 7.58000, 7.67440, 7.76000, 7.89000, 7.93000, 8.00000,
    optPhotMaxE_ / eV
    });

  auto absLength = n4::scale_by(cm, {
    noAbsLength_ / cm,
    noAbsLength_ / cm,  200.0 ,  200.0 ,   90.0 ,   45.0 ,  45.0  ,   30.0 ,
                 24.0,   21.0 ,   20.0 ,   19.0 ,   16.0 ,  14.0  ,   13.0 ,
                  8.5,    8.0 ,    6.0 ,    1.5 ,    1.2 ,   1.0  ,     .65,
                   .4,     .37,     .32,     .28,     .22,    .215,   5e-5 ,
    5e-5
    });

  return n4::material_properties()
    .add("RINDEX"   ,  ri_energy, rIndex)
    .add("ABSLENGTH", abs_energy, absLength)
    .done();
}  

  /// Fake Grid ///
  G4MaterialPropertiesTable* fakegrid_props(G4double transparency,
                                      G4double thickness)
  {
    // ABSORPTION LENGTH
    G4double abs_length   = -thickness/log(transparency);
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {abs_length, abs_length};
    return n4::material_properties()
           .copy_from(LXe_props(), {
            "SCINTILLATIONCOMPONENT1",
            "SCINTILLATIONCOMPONENT2",
            "RINDEX"})
//           .copy_NEW_from(xenon_pt, {
           .copy_from(LXe_props(), {
            "SCINTILLATIONTIMECONSTANT1",
            "SCINTILLATIONTIMECONSTANT2",
            "SCINTILLATIONYIELD",
            "SCINTILLATIONYIELD1",
            "SCINTILLATIONYIELD2",
            "RESOLUTIONSCALE"})
           .add("ABSLENGTH", abs_energy, absLength)
           .done();
  }

  G4MaterialPropertiesTable* tpb_props()
  {
    // Data from https://doi.org/10.1140/epjc/s10052-018-5807-z
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    std::vector<G4double> rIndex_energies = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> TPB_rIndex      = {1.67    , 1.67};
    mpt->AddProperty("RINDEX", rIndex_energies, TPB_rIndex);

    // ABSORPTION LENGTH
    // Assuming no absorption except WLS
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // WLS ABSORPTION LENGTH (Version NoSecWLS)
    // The NoSecWLS is forced by setting the WLS_absLength to noAbsLength_
    // for wavelengths higher than 380 nm where the WLS emission spectrum starts.
    std::vector<G4double> WLS_abs_energy = {
      optPhotMinE_,
      CLHEP::h_Planck * CLHEP::c_light / (380. * nm),  CLHEP::h_Planck * CLHEP::c_light / (370. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (360. * nm),  CLHEP::h_Planck * CLHEP::c_light / (330. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (320. * nm),  CLHEP::h_Planck * CLHEP::c_light / (310. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (300. * nm),  CLHEP::h_Planck * CLHEP::c_light / (270. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (250. * nm),  CLHEP::h_Planck * CLHEP::c_light / (230. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (210. * nm),  CLHEP::h_Planck * CLHEP::c_light / (190. * nm),
      CLHEP::h_Planck * CLHEP::c_light / (170. * nm),  CLHEP::h_Planck * CLHEP::c_light / (150. * nm),
      optPhotMaxE_
    };

    std::vector<G4double> WLS_absLength = {
      noAbsLength_,                 // ~6200 nm
      noAbsLength_,   50. * nm,     // 380 , 370 nm
      30. * nm,      30. * nm,      // 360 , 330 nm
      50. * nm,      80. * nm,      // 320 , 310 nm
      100. * nm,     100. * nm,     // 300 , 270 nm
      400. * nm,     400. * nm,     // 250 , 230 nm
      350. * nm,     250. * nm,     // 210 , 190 nm
      350. * nm,     400. * nm,     // 170 , 150 nm
      400. * nm                     // ~108 nm
    };

    //for (int i=0; i<WLS_abs_energy.size(); i++)
    //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
    //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
    //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
    mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy, WLS_absLength);

    // WLS EMISSION SPECTRUM
    // Implemented with formula (7), with parameter values in table (3)
    // Sampling from ~380 nm to 600 nm <--> from 2.06 to 3.26 eV
    const G4int WLS_emi_entries = 120;
    std::vector<G4double> WLS_emi_energy;
    for (int i=0; i<WLS_emi_entries; i++)
      WLS_emi_energy.push_back(2.06 * eV + 0.01 * i * eV);

    std::vector<G4double> WLS_emiSpectrum;
    G4double A      = 0.782;
    G4double alpha  = 3.7e-2;
    G4double sigma1 = 15.43;
    G4double mu1    = 418.10;
    G4double sigma2 = 9.72;
    G4double mu2    = 411.2;

    for (int i=0; i<WLS_emi_entries; i++) {
      G4double wl = (CLHEP::h_Planck * CLHEP::c_light / WLS_emi_energy[i]) / nm;
      WLS_emiSpectrum.push_back(A * (alpha/2.) * exp((alpha/2.) *
                          (2*mu1 + alpha*pow(sigma1,2) - 2*wl)) *
                          erfc((mu1 + alpha*pow(sigma1,2) - wl) / (sqrt(2)*sigma1)) +
                          (1-A) * (1 / sqrt(2*pow(sigma2,2)*3.1416)) *
                                exp((-pow(wl-mu2,2)) / (2*pow(sigma2,2))));
      // G4cout << "* TPB WLSemi:  " << std::setw(4)
      //        << wl << " nm -> " << WLS_emiSpectrum[i] << G4endl;
    };
    mpt->AddProperty("WLSCOMPONENT", WLS_emi_energy, WLS_emiSpectrum);

    // WLS Delay
    mpt->AddConstProperty("WLSTIMECONSTANT", 1.2 * ns);

    // WLS Quantum Efficiency
    // According to the paper, the QE of TPB depends on the incident wavelength.
    // As Geant4 doesn't allow this possibility, it is set to the value corresponding
    // to Xe scintillation spectrum peak.
    mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.65);

    return mpt;
  }