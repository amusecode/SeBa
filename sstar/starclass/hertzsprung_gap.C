#include "hertzsprung_gap.h"
#include "main_sequence.h"

// This constructor 
// copies the old main_sequence star into the newly created 
// hertzsprung_gap star and destroys the old main_sequence star 
// and finally evolves the hertzsprung_gap in order to determine its
// appearence.
//
// ANSI C++ first creates the base class before the derived classes are
// created. 

     hertzsprung_gap::hertzsprung_gap(main_sequence & m) : single_star(m) {

      delete &m; 	

      //      real m_tot    = get_total_mass();
      //      core_mass     = min(TAMS_helium_core_mass(), m_tot);
      //      envelope_mass = m_tot - core_mass;
      //      core_radius   = helium_core_radius();
                  
         if (relative_mass != get_total_mass()){
             cerr<<"error constructor HG: relative_mass != get_total_mass()"<<endl;   
         }
         
      // (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      adjust_next_update_age();
 
      // (GN Oct 25 2010) try to fix transition MS -> HG in merger MS + *
//      if (core_mass > 0) add_mass_to_accretor(1e-5, false);


      instantaneous_element();
      evolve_core_mass();
      small_envelope_perturbation();   
         
      update();
      post_constructor();
   }


#if 0
void hertzsprung_gap::adjust_initial_star() {

  if(relative_age<=0)
    relative_age = max(main_sequence_time(), relative_age);
}
#endif

star* hertzsprung_gap::reduce_mass(const real mdot) {

    if (envelope_mass<=mdot) {
        envelope_mass = 0;

//        (SPZ+GN: 27 Jul 2000)
//        // non degenerate core < helium_dwarf_mass_limit always(!) become
//        // white dwarfs
//        if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
//            relative_mass < cnsts.parameters(
//                    upper_ZAMS_mass_for_degenerate_core)) {
//            star_transformation_story(Helium_Dwarf);
//            return dynamic_cast(star*, new white_dwarf(*this));
//        } 
//        else {
//            star_transformation_story(Helium_Star);
//            return dynamic_cast(star*, new helium_star(*this));
//        }
        
        real m_HeF = helium_flash_mass(metalicity);
        if (relative_mass < m_HeF){
            star_transformation_story(Helium_Dwarf);
            return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
        }
        else {
            star_transformation_story(Helium_Star);
            return dynamic_cast(star*, new helium_star(*this));
        }
    }
    real mc_bgb = terminal_hertzsprung_gap_core_mass(get_total_mass()-mdot, metalicity);
    real m_FGB = helium_ignition_mass(metalicity);
    
    if ( core_mass > mc_bgb || relative_mass > m_FGB){
        //the evolution of m_core, luminosity, timescales decouples from M
        // relative mass is no longer kept at same value as total mass
        envelope_mass -= mdot;
    }
    else{
        adjust_age_after_mass_loss(mdot, true);
        envelope_mass -= mdot;
        if (relative_mass > get_total_mass()){
            update_relative_mass(get_total_mass());
        }
    }
    return this;
}

star* hertzsprung_gap::subtrac_mass_from_donor(const real dt, real& mdot) {

    mdot = mdot_limit(dt, mdot);
      
      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;

        // (SPZ+GN: 27 Jul 2000)
        // non degenerate core < helium_dwarf_mass_limit always(!) become
        // white dwarfs
        //if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
        //    relative_mass < cnsts.parameters(
        //            upper_ZAMS_mass_for_degenerate_core)) {
        //  star_transformation_story(Helium_Dwarf);
        //  return dynamic_cast(star*, new white_dwarf(*this));
        //} 
        //else {
        //  star_transformation_story(Helium_Star);
        //  return dynamic_cast(star*, new helium_star(*this));
        //}
        
          real m_HeF = helium_flash_mass(metalicity);
          if (relative_mass < m_HeF){
              star_transformation_story(Helium_Dwarf);
              return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
          }
          else {
              star_transformation_story(Helium_Star);
              return dynamic_cast(star*, new helium_star(*this));
              
          }
      }
    
    
    real mc_bgb = terminal_hertzsprung_gap_core_mass(get_total_mass()-mdot, metalicity);
    real m_FGB = helium_ignition_mass(metalicity);

    if ( core_mass > mc_bgb || relative_mass > m_FGB){
        //the evolution of m_core, luminosity, timescales decouples from M
        // relative mass is no longer kept at same value as total mass
        envelope_mass -= mdot;
    }
    else{
        adjust_age_after_mass_loss(mdot, true);
        envelope_mass -= mdot;
        if (relative_mass > get_total_mass()){
            update_relative_mass(get_total_mass());
        }
    }
    
    adjust_donor_radius(mdot);
    return this;    
}



// add mass to accretor
// is a separate function (see single_star.C) because rejuvenation
real hertzsprung_gap::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {

    if (mdot<0) {
        cerr << "hertzsprung_gap::add_mass_to_accretor(mdot=" << mdot 
        << ", dt=" << dt << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        cerr << "Action: put mdot to zero!" << endl;
        return 0;
    }
    
    bool update_age = false;

    if(hydrogen){
        //hydrogen accretion
        mdot = accretion_limit(mdot, dt);

	envelope_mass += mdot;
	accreted_mass += mdot;

	if (relative_mass < get_total_mass()){
	  update_age = true;
	  update_relative_mass(get_total_mass());
	}
	
        
        adjust_accretor_radius(mdot, dt);
        
    }
    else{
        //for the moment assume helium accretion
        // for the moment no adjust_accretor_radius
        
        mdot = accretion_limit_eddington(mdot, dt);
        core_mass += mdot;
        update_relative_mass(relative_mass + mdot);
        
	update_age = true;
    }


    if (update_age) {

        //adjust age part
        real mc_ehg = terminal_hertzsprung_gap_core_mass(relative_mass, metalicity);
        real m5_25 = pow(relative_mass, 5.25);            
        real rho = (1.586 + m5_25) / (2.434 + 1.02*m5_25);
        real tau = (core_mass / mc_ehg - rho ) / (1. - rho);
        real t_ms = main_sequence_time();
        real t_bgb = base_giant_branch_time(relative_mass, metalicity);
        relative_age = t_ms + tau * (t_bgb - t_ms);
        last_update_age = t_ms;

        if (tau < 0.){
            real (single_star::*fptr)(const real, real) = &single_star::initial_hertzsprung_gap_core_mass;                   
            real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity);     
            update_relative_mass(m_rel);
            last_update_age = main_sequence_time(relative_mass, metalicity);
            relative_age = last_update_age;     
            evolve_core_mass();
        }
        if (tau > 1.){
            real (single_star::*fptr)(const real, real) = &single_star::terminal_hertzsprung_gap_core_mass;        
            real m_rel = linear_function_inversion(fptr, relative_mass, core_mass, metalicity);     
            update_relative_mass(m_rel);
            last_update_age = main_sequence_time(relative_mass, metalicity);
            relative_age = next_update_age;            
            evolve_core_mass();
        }
    }
    set_spec_type(Accreting);
    return mdot;
}

#if 0
// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void hertzsprung_gap::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {

      real m_rel_new;
      real m_tot_new = get_total_mass() + mdot;
      if (m_tot_new>relative_mass)
         m_rel_new = m_tot_new;
      else m_rel_new = relative_mass;

      real t_ms_old = main_sequence_time();
      real t_hg_old = hertzsprung_gap_time() - t_ms_old;

      real z_new = get_metalicity();
      real t_ms_new = main_sequence_time(m_rel_new, z_new);

      //For now, we keep metalicity constant (SPZ: 29 May 2001)
      real t_bgb = hertzsprung_gap_time(m_rel_new, z_new); 
      real t_hg_new = t_bgb - t_ms_new;

//      real dtime = relative_age - t_ms_old;

      // (GN+SPZ May  4 1999) update last_update_age
      last_update_age = t_ms_new; 
   
//      relative_age = t_ms_new 
//                   + dtime*(t_hg_new/t_hg_old);
//
//      if (rejuvenate)
//           relative_age *= rejuvenation_fraction(mdot/m_tot_new); 
//


    // rejuvenation based on keeping mc constant and finding relative age for new mass
    real mc_ehg = terminal_hertzsprung_gap_core_mass(relative_mass, metalicity);
    real m5_25 = pow(relative_mass, 5.25);    
    real rho = (1.586 + m5_25) / (2.434 + 1.02*m5_25);

    real mc_ehg_new = terminal_hertzsprung_gap_core_mass(m_rel_new, z_new);
    real m5_25_new = pow(m_rel_new, 5.25);
    real rho_new = (1.586 + m5_25_new) / (2.434 + 1.02*m5_25_new);

    
    relative_age = t_ms_new + t_hg_new * ((1-rho)/(1-rho_new) * mc_ehg/mc_ehg_new * (relative_age - t_ms_old) / t_hg_old + 
                (rho*mc_ehg - rho_new * mc_ehg_new)/ (1-rho_new) / mc_ehg_new);
                
                
      relative_age = max(relative_age, 
			 last_update_age + cnsts.safety(minimum_timestep)); 
    
      relative_age = min(relative_age, t_bgb);
      
      
        // next_update_age should not be reset here
        // next_update_age = t_ms_new + t_hg_new;
}
#endif

// Age adjustment especially for (wind) mass loss
// It is part of the single star evolution, 
// so it can include information from tracks
void hertzsprung_gap::adjust_age_after_mass_loss(const real mdot, const bool rejuvenate=true) {

    real m_rel_new;
    real m_tot_new = get_total_mass() - mdot;
    if (m_tot_new<relative_mass)
        m_rel_new = m_tot_new;
    else m_rel_new = relative_mass;
    
    real t_ms_old = main_sequence_time();
    real t_hg_old = hertzsprung_gap_time() - t_ms_old;
    
    real z_new = get_metalicity();
    real t_ms_new = main_sequence_time(m_rel_new, z_new);
    //For now, we keep metalicity constant (SPZ: 29 May 2001)
    real t_bgb = hertzsprung_gap_time(m_rel_new, z_new); 
    real t_hg_new = t_bgb - t_ms_new;
    
    real dtime = relative_age - t_ms_old;
    
    // (GN+SPZ May  4 1999) update last_update_age
    last_update_age = t_ms_new;
    //following HPT tracks only update the relative_age relative to the length of the phase

    relative_age = t_ms_new + dtime*(t_hg_new/t_hg_old);
//    if (rejuvenate){
//        real mdot_fr = -1.*mdot/m_tot_new;
//        real rejuvenation = (1-pow(mdot_fr,
//                                   cnsts.parameters(rejuvenation_exponent)));
//        relative_age *= rejuvenation; 
//    
//    }
//    

    relative_age = max(relative_age, 
                       last_update_age + cnsts.safety(minimum_timestep)); 
    relative_age = min(relative_age, t_bgb);

    // next_update_age should not be reset here
    // next_update_age = t_ms_new + t_hg_new;
}


// Adiabatic response function for hertzsprung_gap star.
// Polynomial fit to Hjellming and Webbink 1987 ApJ, 318, 804
real hertzsprung_gap::zeta_adiabatic() {

#if 0
      real z;

      real x = core_mass/get_total_mass();
      real A = -0.220823;
      real B = -2.84699;
      real C = 32.0344;
      real D = -75.6863;
      real E = 57.8109;

      if (get_relative_mass()<=0.4)
         z = -cnsts.mathematics(one_third);
      else if (low_mass_star())
         z = A + x*(B + x*(C + x*(D + x*E)));
      else if (medium_mass_star())
         z = 2.25;                 // 15 according to Pols & Marinus 1994
      else                         // We do it differently.
         z = 2.25;                 // lekker puh.
#endif
      real z = 4; // this is neede to prevent Thermal in extreme mass
                     // ratio systems ... 
// (GN+SPZ Apr 29 1999) Pols & Marinus 1994 were maybe right: not!

      return z;
   }

// Thermal response function for hertzsprung_gap star.
real hertzsprung_gap::zeta_thermal() {

      real z;

      if (get_relative_mass()<=0.4)
         z = 0;         // no better estimate present.
      else if (low_mass_star())
         z = -2; 	// -10 according to Pols & Marinus 1994
      else              // Changed to current values
         z = -2;	// by (SPZ+GN: 1 Oct 1998)

      return z;
   }
 
void hertzsprung_gap::adjust_next_update_age() {
   
    real t_ms = main_sequence_time();
    if (relative_age<t_ms) {
        cerr << "WARNING: relative_age < t_ms in Hertzsprung_gap"<<endl;
        relative_age = t_ms;
    }
    
    real t_bhg = hertzsprung_gap_time();
    if(t_ms<t_bhg) 
        next_update_age = hertzsprung_gap_time();
    else {
        cerr << "WARNING:hertzsprung_gap::adjust_next_update_age()" << endl;
        cerr << "main_sequence time exceeds hertzprung_gap time"<<endl;
        PRC(t_ms);PRL(t_bhg);
        dump(cerr, false);
        cerr << flush;
        exit(-1);
    }


}

void hertzsprung_gap::detect_spectral_features() {

// 		Use standard spectral feature detection.
      single_star::detect_spectral_features();

      if (accreted_mass>=cnsts.parameters(B_emission_star_mass_limit))
	spec_type[Emission]=Emission;
   }


// Stellar Gyration radii squared for detmination of
// angular momentum.
// Implemented by (SPZ+GN: 1 Oct 1998)
real hertzsprung_gap::gyration_radius_sq() {

//  return cnsts.parameters(radiative_star_gyration_radius_sq); 


  
  
// (SilT 13 Feb 22) 
// based on 

    real t_ms = main_sequence_time();
    real t_bgb = base_giant_branch_time(relative_mass, metalicity);
//    PRL((relative_age - t_ms)/(t_bgb-t_ms));
    real tau = (relative_age - t_ms)/(t_bgb-t_ms);
    tau = max(min(tau,1.),0.);



    const int size_mbins = 19;
	const int size_tbins = 4;
	real mass_array [size_mbins] = {0.6, 0.8, 0.85, 0.92, 1., 1.2, 1.4, 1.6,1.8, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10.}; 

    int im_l, im_u; 
    if (relative_mass <= mass_array[0]) im_l = im_u = 0;
    else if (relative_mass >= mass_array[size_mbins-1]) im_l = im_u = size_mbins-1;
    else {
        for (im_u = 0; im_u < size_mbins; im_u++){
                if (relative_mass < mass_array[im_u]) break; }
        im_l = im_u-1;
    }
    real m_u = mass_array[im_u];
    real m_l = mass_array[im_l];
//    PRC(m_u);PRC(m_l);PRL(relative_mass);



    real timestamp_array[size_mbins][size_tbins] =  {{ 75741300000.0 , 77257131625.6 , 78969864419.1 , 79733800000.0 }, { 26518300000.0 , 27114270661.8 , 27657265899.2 , 27915800000.0 }, {20847300000.0 , 21321153499.0 , 21744367312.3 , 21945800000.0 }, {15229700000.0 , 15591061040.5 , 15883890570.0 , 16032100000.0 }, {11003100000.0 , 11257154503.0 , 11466963383.2 , 11582700000.0 }, {5615540000.0 , 5783000053.66 , 5849741001.99 , 5911260000.0 }, {3371530000.0 , 3406562666.61 , 3443368437.44 , 3502810000.0 }, {2245330000.0 , 2252296758.5 , 2264108571.6 , 2289900000.0 }, {1582920000.0 , 1584938154.99 , 1588669835.1 , 1600400000.0 }, {1164160000.0 , 1165089484.55 , 1167178336.27 , 1173580000.0 }, {619370000.0 , 620547256.574 , 621598319.144 , 623600000.0 }, {377537000.0 , 378326688.348 , 378752266.532 , 379788000.0 }, {179236000.0 , 179926477.152 , 180029442.673 , 180102000.0 }, {104016000.0 , 104326729.76 , 104358147.878 , 104444000.0 }, {68374100.0 , 68487683.8658 , 68533114.8176 , 68621800.0 }, {48900700.0 , 49001667.0953 , 49043677.7384 , 49060300.0 }, {37183500.0 , 37215242.4738 , 37240648.6147 , 37295000.0 }, {29555300.0 , 29586259.3519 , 29607085.8442 , 29637500.0 }, {24305900.0 , 24333846.5383 , 24346878.3625 , 24369000.0 }};



    real beta_array[size_mbins][size_tbins] =  {{
        0.299535 , 0.304388 , 0.330665 , 0.325938 }, {
        0.277057 , 0.287441 , 0.348146 , 0.323157 }, {
        0.269815 , 0.28061 , 0.351148 , 0.328652 }, {
        0.258695 , 0.272566 , 0.355065 , 0.335492 }, {
        0.244822 , 0.258965 , 0.35817 , 0.341594 }, {
        0.20762 , 0.219628 , 0.355149 , 0.354192 }, {
        0.177592 , 0.206151 , 0.366193 , 0.362943 }, {
        0.161864 , 0.174966 , 0.369244 , 0.36979 }, {
        0.155747 , 0.167015 , 0.37857 , 0.373385 }, {
        0.153653 , 0.167129 , 0.378847 , 0.378327 }, {
        0.152961 , 0.157364 , 0.18149 , 0.30257 }, {
        0.153845 , 0.155735 , 0.170567 , 0.289095 }, {
        0.1549 , 0.1577 , 0.1624 , 0.1723 }, {
        0.153 , 0.1588 , 0.1666 , 0.2267 }, {
        0.1461 , 0.1475 , 0.1798 , 0.3322 }, {
        0.1433 , 0.1439 , 0.1457 , 0.1477 }, {
        0.1403 , 0.1434 , 0.1541 , 0.2164 }, {
        0.1366 , 0.1376 , 0.1441 , 0.2201 }, {
        0.1309 , 0.133 , 0.1406 , 0.2262 }};

    

    //find time for m_u
    real time_yrs_u = timestamp_array[im_u][0] + tau*(timestamp_array[im_u][size_tbins-1] - timestamp_array[im_u][0]);
    int it_uu, it_ul;
    if (time_yrs_u <=  timestamp_array[im_u][0]) it_uu = it_ul = 0 ;
    else if (time_yrs_u >=  timestamp_array[im_u][size_tbins-1]) it_uu = it_ul = size_tbins-1 ;
    else{
        for (it_uu = 0; it_uu < size_tbins; it_uu++){
                if (time_yrs_u < timestamp_array[im_u][it_uu]) break; }
        it_ul = it_uu-1;        
    } 
    real t_uu = timestamp_array[im_u][it_uu];
    real t_ul = timestamp_array[im_u][it_ul];
//    PRC(t_uu);PRC(t_ul);PRC(relative_age);PRL(time_yrs_u);


    //find time for m_l
    real time_yrs_l = timestamp_array[im_l][0] + tau*(timestamp_array[im_l][size_tbins-1] - timestamp_array[im_l][0]);
    int it_lu, it_ll;
    if (time_yrs_l <=  timestamp_array[im_l][0]) it_lu = it_ll = 0 ;
    else if (time_yrs_l >=  timestamp_array[im_l][size_tbins-1]) it_lu = it_ll = size_tbins-1 ;
    else{
        for (it_lu = 0; it_lu < size_tbins; it_lu++){
                if (time_yrs_l < timestamp_array[im_l][it_lu]) break; }
        it_ll = it_lu-1;        
    } 
    real t_lu = timestamp_array[im_l][it_lu];
    real t_ll = timestamp_array[im_l][it_ll];
//    PRC(t_lu);PRC(t_ll);PRC(relative_age);PRL(time_yrs_l);
    

    // for readibility keep this separate from previous two loops            
    //step three, interpolation in mass and time...
    real beta_m_l_at_t, beta_m_u_at_t, beta;

    if (it_uu == it_ul){
        beta_m_u_at_t = beta_array[im_u][it_uu];
    }
    else {
        beta_m_u_at_t = beta_array[im_u][it_ul] + (beta_array[im_u][it_uu] - beta_array[im_u][it_ul]) * (time_yrs_u - t_ul) / (t_uu - t_ul);
    }
    if (it_lu == it_ll){
        beta_m_l_at_t = beta_array[im_l][it_lu];
    }
    else {
        beta_m_l_at_t = beta_array[im_l][it_ll] + (beta_array[im_l][it_lu] - beta_array[im_l][it_ll]) * (time_yrs_l - t_ll) / (t_lu - t_ll);
    }


    if (im_u==im_l)
        beta = beta_m_u_at_t;
    else    
        beta = beta_m_l_at_t + (beta_m_u_at_t - beta_m_l_at_t) * (relative_mass - m_l) / (m_u - m_l);
        
        
//    PRC( beta_array[im_l][it_ll]);PRC(beta_array[im_l][it_lu]);       
//    PRC( beta_array[im_u][it_ul]);PRL(beta_array[im_u][it_uu]);       
//    PRC(beta_m_l_at_t);PRC(beta_m_u_at_t);PRL(beta); 
//    PRL(cnsts.parameters(radiative_star_gyration_radius_sq));
    return beta*beta;           

//  return cnsts.parameters(radiative_star_gyration_radius_sq); 

  
}


// Section 7.2 in Hurley, Pols & Tout 2000
real hertzsprung_gap::convective_envelope_mass(){

    real t_ms = main_sequence_time(); // function of relative mass and metalicity
    real t_bgb = base_giant_branch_time(); // function of relative mass and metalicity
    real tau = (relative_age - t_ms)/ (t_bgb - t_ms);
    return tau * envelope_mass;
}

// Eq. 39-40 in Hurley, Tout & Pols 2002
real hertzsprung_gap::convective_envelope_radius(){

    real t_ms = main_sequence_time(); // function of relative mass and metalicity
    real t_bgb = base_giant_branch_time(); // function of relative mass and metalicity
    real tau = (relative_age - t_ms)/ (t_bgb - t_ms);
    return max(0., sqrt(tau) * (radius-helium_core_radius())); //function of relative_mass, core_mass and metalicity
}


void hertzsprung_gap::update_wind_constant() {
#if 0
// (GN+SPZ Apr 28 1999) new fits to Maeder, de Koter and common sense
// wind_constant is fraction of envelope lost in nuclear lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 1 Oct 1998)
    
  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    
    if (relative_mass < 85)
      wind_constant = meader_fit_dm;
    else {// constant
      real final_mass = 30;
      wind_constant = relative_mass - final_mass;
    }

  } 
  else { // (GN+SPZ Apr 29 1999) 1% loss on hg

    wind_constant = 0.01*relative_mass;
  }

  wind_constant = max(wind_constant, 0.0);
#endif
    
    // wind_constant is in solar masses per year
    // Should be updated after mass accretion
    // (ST: 17 Sep 2009)
    
    real dm_dj_v = 0;
    
    // Nieuwenhuijzen & de Jager 1990
    // Massive stars
    real dm_dj = 0;
    if (luminosity > 4000.) {
        real x = min(1.0, (luminosity - 4000.0)/500.0);
        dm_dj = x * 9.6310E-15 * pow(radius, 0.81) * pow(luminosity, 1.24) * 
        pow(get_total_mass(), 0.16)*pow(metalicity/cnsts.parameters(solar_metalicity), 0.85);
    }
    
    // Vink 2000, 2001
    // Massive stars, including multi scattering effects
    real dm_v = 0;
    if (luminosity > 4000 && metalicity > cnsts.parameters(solar_metalicity)/30. && metalicity < 3*cnsts.parameters(solar_metalicity)){
        real temp = temperature();
        real sigma;//electron scattering cross section
        if (temp >= 35000){
            sigma = 0.34;
        }
        else if(temp < 30000){
            sigma = 0.31;
        }
        else {
            sigma = 0.32;
        }
        real rad_acc = 7.66E-5 * sigma * luminosity / get_total_mass();
        real log_density = -14.94 + 0.85 * log10(metalicity/cnsts.parameters(solar_metalicity)) +3.1857*rad_acc; //Eq.23 Vink 2001
        real Tjump = (61.2 + 2.59*log_density)*1000; //Eq.15 Vink 2001
        real Tjump_low = (100. + 6. * log_density)*1000; //Eq.6 Vink 2000
        real arg_dm_v;
        real T_smooth = 1500.;
        real T_smooth_below = 1000.;
        real cnsts_dm_v[9];
        real cnsts_dm_v_above[] = {2.6, -6.697, 2.194, -1.313, -1.226, 0.933, 0, -10.92, 0.85};
        real cnsts_dm_v_below[] = {1.3, -6.688, 2.210, -1.339, -1.601, 1.07, 1.07, 0, 0.85};
        
        
        if (rad_acc >0.5 || temp > 50000) {
            // vink approaches LBV, stop? transition needed? possible for low metallicities
            dm_v = 0;
            dm_dj_v = dm_dj;
        }
        else {
            if (temp <= Tjump-T_smooth){
                // smooth out second instability jump
                for (int i_t=0; i_t< 9;i_t++){
                    cnsts_dm_v[i_t] = cnsts_dm_v_below[i_t];
                }
                
                if(temp <= Tjump_low - T_smooth_below){
                    cnsts_dm_v[0] = 0.7;
                    cnsts_dm_v[1] = -5.990;
                }
                else if (temp < Tjump_low + T_smooth_below){
                    real scale_T = (Tjump_low+T_smooth_below - temp)/(2*T_smooth_below);
                    cnsts_dm_v[0] = (1.-scale_T) * cnsts_dm_v_below[0] + scale_T * 0.7;
                    cnsts_dm_v[1] = (1.-scale_T) * cnsts_dm_v_below[1] + scale_T * -5.990;
                }
            }
            else if(temp > Tjump + T_smooth){
                for (int i_t=0; i_t< 9;i_t++){
                    cnsts_dm_v[i_t] = cnsts_dm_v_above[i_t];
                }
            }
            else {
                //smooth out first instability jump
                real scale_T = (Tjump+T_smooth - temp)/(2*T_smooth);
                for (int i_t=0; i_t< 9;i_t++){
                    cnsts_dm_v[i_t] = scale_T * cnsts_dm_v_below[i_t] 
                    + (1.-scale_T) * cnsts_dm_v_above[i_t];
                }
            }
            
            arg_dm_v  = cnsts_dm_v[1] 
            +cnsts_dm_v[2]*log10(luminosity/1.E5) 
            +cnsts_dm_v[3]*log10(get_total_mass()/ 30) 
            +cnsts_dm_v[4]*log10(cnsts_dm_v[0]/2.0)
            +cnsts_dm_v[5]*log10(temp/40000)
            +cnsts_dm_v[6]*log10(2.) // (T/20000) below Tjump
            +cnsts_dm_v[7]*pow(log10(temp/40000),2)
            +cnsts_dm_v[8]*log10(metalicity/cnsts.parameters(solar_metalicity));
            
            dm_v = pow(10, arg_dm_v);
	    // decrease Vink winds by a factor 3 (e.g. see Bjorklund et al. 2020)
            // (FK 7 Oct 2021)
	    dm_v = dm_v / 3.;
            dm_dj_v = dm_v;
            
            
            if (temp < 8000){
                // line driven winds no longer efficient
                // see Achmad et al 1997
                dm_v = dm_v * 200. / (8200.-temp);
                dm_dj_v = max(max(dm_v, dm_dj), 0.);
            }
        }
    }   
    else
        dm_dj_v = dm_dj;
    
    
    // Reimers 1975
    // GB like stars
    real neta = 0.5; 
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();
    
//    //Schroder & Cuntz
//    // cool GB like stars
//    real neta_sc = 8.E-14; 
//    real surface_gravity = get_total_mass()/ pow(radius, 2);
//    real dm_sc = neta_sc * 4.E-13 * radius * luminosity / get_total_mass() 
//           * pow(temperature()/4000, 3.5) * (1 + 1./(4300*surface_gravity));
    
    
    //based on Nugis & Lamers
    // eq 8.4 in Gijs' thesis Chapter 8
    //Reduced WR-like mass loss for small H-envelope mass
    //real mu = (get_total_mass()-core_mass)/get_total_mass() * min(5.0,max(1.2, pow(luminosity/7.E4,-0.5)));
    real dm_wr = 0;
    //if ( mu < 1.){
    //    //factor (1.-mu) should be checked e.g. with resulting # BH in binaries
    //    dm_wr = 1.38E-08 * pow(get_total_mass(), 2.87) * (1.-mu);
    //}
    
    //LBV
    real dm_lbv = 0;
    real x_lbv = 1.0E-5*radius*sqrt(luminosity);
    if(luminosity > 6.0E5 && x_lbv > 1.0) {
        dm_lbv = 0.1 * pow(x_lbv-1.0, 3)*(luminosity/6.0E5-1.0);
    }
    
    wind_constant = max(max(max(dm_wr, dm_dj_v), dm_r), 0.0) + dm_lbv;
        
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Metalicity dependency from HPT 200

void hertzsprung_gap::instantaneous_element() {

    luminosity       = hertzsprung_gap_luminosity(relative_age, 
                                                  relative_mass, metalicity);
    radius           = hertzsprung_gap_radius(relative_age, relative_mass, 
                                              get_total_mass(), metalicity);
}

// Evolve a main_sequence star upto time argument according to
// the new 2000 models.
void hertzsprung_gap::evolve_element(const real end_time) {
      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {
          instantaneous_element();
          evolve_core_mass();
          small_envelope_perturbation();
          
          // if no envelope make transition to remnants
          // just as a procedure: reduce_mass with 1
          if (envelope_mass <= 0){
              reduce_mass(1.);
              return;
          }
      }
      else {
        if (relative_mass <= helium_ignition_mass(metalicity) ){
	        star_transformation_story(Sub_Giant);
            new sub_giant(*this);
	    return;
	    }
	    else {
          star_transformation_story(Horizontal_Branch);
          new horizontal_branch(*this);
          return;
	    }
      }

      update();
      stellar_wind(dt);
}


// (SilT Dec 10 2020) 
real hertzsprung_gap::get_evolve_timestep() {

  real timestep = min((next_update_age - last_update_age )/ cnsts.safety(number_of_steps), 
                      next_update_age - relative_age - 0.5 * cnsts.safety(minimum_timestep));   
   
   
//     (SilT Dec 10 2020) extra safety measure
//     For large L, large Mdot, that can fluctuate due to bi-instability jumps in the line driven winds of Vink
    real dt_mdot = timestep * pow(50000/luminosity, 0.75) / pow(metalicity/cnsts.parameters(solar_metalicity), 0.85) ;    

   return max(min(timestep, dt_mdot), cnsts.safety(minimum_timestep));
   
}


//Eq.7+
real hertzsprung_gap::terminal_hertzsprung_gap_luminosity(const real mass, 
						      const real z) {
  
  real l_ehg;
  if (mass<helium_ignition_mass(z))
    l_ehg = base_giant_branch_luminosity(mass, z);
  else
    l_ehg = helium_ignition_luminosity(mass, z);

    return  l_ehg;
}

//Eq.7+
real hertzsprung_gap::terminal_hertzsprung_gap_radius(const real mass, 
						      const real mass_tot, const real z) {
  
    real r_ehg;
    if (mass<helium_ignition_mass(z)) {
        real l_bgb = base_giant_branch_luminosity(mass, z);
        r_ehg = giant_branch_radius(l_bgb, mass_tot, z);
    }
    else {
        r_ehg = helium_ignition_radius(mass, mass_tot, z);
        //safety check
        //these lines are not in the HPT2000 article, but they are in the HPT2000 code
        //in case a massive star skips the blue loop phase, 
        // the stellar radius should continue smoothly 
        real l_HeI = helium_ignition_luminosity(mass, z);
        real r_agb =  AGB_radius(l_HeI, mass, mass_tot, z);
        if (r_ehg > r_agb){
            if (mass >= 12.0){
                r_ehg = r_agb;
                cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: R_ehg is set to r_agb"<<endl;
            }
            else {
                cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: R_AGB(L_HeI) < R_mHe, skipping blue loop?"<<endl;
            }
        }
//        if (blue_phase_timescale(mass, mass_tot, z) < cnsts.safety(tiny) && abs(r_ehg-r_agb)> cnsts.safety(tiny)) {
//            cerr<<"WARNING in hertzsprung_gap:: terminal_hertzsprung_gap_radius: t_bl <0, but r_ehg != r_agb)"<<endl;; 
//        }
        
    }  
    return r_ehg;
}

real hertzsprung_gap::terminal_hertzsprung_gap_radius() {
    cerr<<"terminal_hertzsprung_gap_radius() without parameters is used. "<<endl;
  return terminal_hertzsprung_gap_radius(relative_mass, get_total_mass(), metalicity);
}


real hertzsprung_gap::terminal_hertzsprung_gap_luminosity() {

  return terminal_hertzsprung_gap_luminosity(relative_mass, metalicity);
}

//real hertzsprung_gap::base_giant_branch_luminosity() {
//
//  return terminal_hertzsprung_gap_luminosity(relative_mass, metalicity);
//}

//Eq.26
real hertzsprung_gap::hertzsprung_gap_luminosity(const real time,
					     const real mass, 
					     const real z) {

  real t_ms = main_sequence_time(mass, z);
  real t_bgb = base_giant_branch_time(mass, z);
  real tau = (time - t_ms)/(t_bgb - t_ms);

  real l_ehg = terminal_hertzsprung_gap_luminosity(mass, z);
  real l_tms = terminal_main_sequence_luminosity(mass, z);

  real l_hg = l_tms * pow(l_ehg/l_tms, tau);
  return l_hg;
}

real hertzsprung_gap::hertzsprung_gap_luminosity(const real time) {

  return hertzsprung_gap_luminosity(time, relative_mass, metalicity);
}

real hertzsprung_gap::hertzsprung_gap_luminosity() {

  return hertzsprung_gap_luminosity(relative_age, relative_mass, metalicity);
}


//Eq.27
real hertzsprung_gap::hertzsprung_gap_radius(const real time,
					     const real mass, 
                         const real mass_tot, 
					     const real z) {

  real t_ms = main_sequence_time(mass, z);
  real t_bgb = base_giant_branch_time(mass, z);
  real tau = (time - t_ms)/(t_bgb - t_ms);
    
  real r_tms = terminal_main_sequence_radius(mass_tot, z);
  real r_thg = terminal_hertzsprung_gap_radius(mass, mass_tot, z);
 
  real r_hg = r_tms * pow(r_thg/r_tms, tau);
  return r_hg;
}

real hertzsprung_gap::hertzsprung_gap_radius(const real time) {
 
  return hertzsprung_gap_radius(time, relative_mass, get_total_mass(), metalicity);
}

real hertzsprung_gap::hertzsprung_gap_radius() {

  return hertzsprung_gap_radius(relative_age, relative_mass, get_total_mass(), metalicity);
}

void hertzsprung_gap::evolve_core_mass(const real time,
				       const real mass,
				       const real z, const real m_core_old) {

  real mc_Hg = hertzsprung_gap_core_mass(time, mass, z, m_core_old);

    if(!update_core_and_envelope_mass(mc_Hg)) {
    cerr << "Update core mass failed in hertzsprung_gap()"<<endl;
  }
}


void hertzsprung_gap::evolve_core_mass() {
    
    evolve_core_mass(relative_age, relative_mass, metalicity, core_mass);
}


// Eq.4 (base_giant_branch_time(mass, z);
// Supersedes ::hertzsprung_gap_time(mass, t_ms) in single star
real hertzsprung_gap::hertzsprung_gap_time(const real mass, const real z) {
    
    real t_Hg = base_giant_branch_time(mass, z);
    return t_Hg;
}

real hertzsprung_gap::hertzsprung_gap_time() {
    
    return hertzsprung_gap_time(get_relative_mass(), get_metalicity());
}


//Eq.30
real hertzsprung_gap::hertzsprung_gap_core_mass(const real time, 
						const real mass,
						const real z, const real m_core_old) {

    real t_ms = main_sequence_time(mass, z);
    real t_bgb = base_giant_branch_time(mass, z);
    real tau = (time - t_ms)/(t_bgb - t_ms);

    real mc_ehg = terminal_hertzsprung_gap_core_mass(mass, z);
    real mc_ihg = initial_hertzsprung_gap_core_mass(mass, z);  
    real m_core = mc_ihg + (mc_ehg - mc_ihg) * tau;
        
    // according to HPT this is important in case of mass loss 
    m_core = max(m_core, m_core_old);   
  
  return m_core;
}

real hertzsprung_gap::helium_core_radius(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real r_c;
    if(mass > m_HeF){
        r_c = helium_star_radius_for_solar_metalicity(m_core);
    }
    else{
        // due to small nucleair burning layer 
        // r_c > white_dwarf_radius
        r_c = 5.* white_dwarf_radius(m_core, 10000.);
    }
    return r_c;
}

real hertzsprung_gap::helium_core_radius(){
    return helium_core_radius(relative_mass, core_mass, metalicity);
}    

real hertzsprung_gap::small_envelope_core_radius(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real r_c;
    if(mass > m_HeF){
        r_c = helium_star_radius_for_solar_metalicity(m_core);
    }
    else{
        r_c = white_dwarf_radius(m_core, 10000.);
    }
    return r_c;
}
  
real hertzsprung_gap::small_envelope_core_radius(){
    return small_envelope_core_radius(relative_mass, core_mass, metalicity);
}    



real hertzsprung_gap::small_envelope_core_luminosity(const real mass, const real m_core, const real z){
    real m_HeF = helium_flash_mass(z);
    real l_c;
    if(mass > m_HeF){
    l_c = helium_star_luminosity_for_solar_metalicity(m_core);
    }
    else{
        l_c = 40.;
    }
    return l_c;
}
  
real hertzsprung_gap::small_envelope_core_luminosity(){
    return small_envelope_core_luminosity(relative_mass, core_mass, metalicity);
}    



