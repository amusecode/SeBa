//
// helium_giant.C
//

#include "super_giant.h"
#include "helium_star.h"
#include "helium_giant.h"

// (GN+SPZ May  3 1999) only for stars more massive than 
// super_giant2neutron_star, below become white_dwarf
// to be done (SPZ+GN: 27 Jul 2000)
// (ST: 10 Sep 2000)  it's different now with the HPT tracks 
helium_giant::helium_giant(super_giant & g) : single_star(g) {

    delete &g;   
    accreted_mass = 0;
    update_relative_helium_mass(get_total_mass());
    
    relative_age = helium_giant_age_core_mass_relation(core_mass, relative_helium_mass);
    last_update_age = helium_main_sequence_time_for_solar_metalicity(relative_helium_mass); 
    
    if (relative_age < last_update_age){     
        
        real (single_star::*fptr)(const real, real) = &single_star::helium_giant_initial_core_mass;        
        real m_he_rel = linear_function_inversion(fptr, relative_mass, core_mass);             
        update_relative_helium_mass(m_he_rel);
        last_update_age = helium_main_sequence_time_for_solar_metalicity(relative_helium_mass);
        relative_age = last_update_age;

        //for orderliness core_mass should be updated,
        //though it also happens at end of constructor
        evolve_core_mass();
    }
    
    adjust_next_update_age();  
    // (SPZ+GN: 27 Jul 2000)
    // core mass has been set in super_giant before constructor is called.
    
    if (relative_age >= next_update_age) {
        // this can happen when mc_max = 1.45*get_total_mass()-0.31        
        // or mc_max = 0.773*m_he_rel - 0.35 --> in this case should m_he_rel be updated?
                
        // happens in evolve_element
        // create_remnant(relative_helium_mass, get_total_mass(), core_mass);  

        //reset radius to prevent continues CE 
        radius = core_radius; 

        return;
    }

    instantaneous_element();
    evolve_core_mass();
    small_envelope_perturbation();   
    effective_radius = radius;    
    update();
    post_constructor();
}


helium_giant::helium_giant(helium_star & h) : single_star(h) {
    delete &h;
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    last_update_age = next_update_age;
    update_relative_helium_mass(get_total_mass());

    adjust_next_update_age();
    //(GN Feb 16 2011) dirty trick to avoid crashes after mergers
    if (relative_age > next_update_age) relative_age = next_update_age;

    instantaneous_element();
    evolve_core_mass();
    small_envelope_perturbation();   
    
    update();
    post_constructor();

}

// Adjust radius & luminosity at relative_age
void helium_giant::instantaneous_element() {

    luminosity       = helium_giant_luminosity_core_mass_relation(relative_age, relative_helium_mass, metalicity);
    radius           = helium_giant_radius(luminosity, relative_helium_mass, get_total_mass(), metalicity);
    
    // don't do here:
    //effective_radius = max(effective_radius, radius)
    //because of small_envelope_perturbation
}

// Evolve a helium_giant upto time argument according to
// the new 2000 models.
void helium_giant::evolve_element(const real end_time) {
    real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;
    
    if (relative_age<=next_update_age) {
        instantaneous_element();
        evolve_core_mass();
        small_envelope_perturbation();
        
        // if no envelope make transition to remnants
        // just as a procedure: reduce_mass with 1
        if(envelope_mass <= 0){
            reduce_mass(1.);
            return;    
        }       
    }
    else {
        create_remnant(relative_helium_mass, get_total_mass(), core_mass);
        return;
    }
    
    update();
    stellar_wind(dt);
}


// (SilT June 15 2012) needed with new stellar wind based on HPT2000
real helium_giant::get_evolve_timestep() {

  real timestep = min((next_update_age - last_update_age )/ cnsts.safety(number_of_steps), 
                      next_update_age - relative_age - 0.5 * cnsts.safety(minimum_timestep));   
   
    //extra safety measure
    // when L and R increase rapidly, so will mdot
    // useful for stars 0.5-2.5 Msun. 
    real dt_mdot = timestep;
    real A_He = AGB_A_He_estimator();
    real R_hezams = helium_star_radius_for_solar_metalicity(relative_mass);
    
    real new_luminosity = helium_giant_luminosity_core_mass_relation(relative_age+timestep, relative_helium_mass, metalicity);
    real new_radius = helium_giant_radius(new_luminosity, relative_helium_mass, get_total_mass(), metalicity);
    dt_mdot = core_mass / ( new_luminosity * A_He) * max(R_hezams/new_radius, 0.01);
   return max(min(timestep, dt_mdot), cnsts.safety(minimum_timestep));
}

void helium_giant::update_relative_helium_mass(const real new_relative_helium_mass) {
    relative_helium_mass = new_relative_helium_mass;
    adjust_next_update_age();
    update_wind_constant();
    
}



real helium_giant::nucleair_evolution_timescale() {        
    return helium_giant_end_time(relative_helium_mass, get_total_mass());
}


void helium_giant::update() {

//  real m_tot = get_total_mass();
//  core_mass = COcore_mass = CO_core_mass();
//  envelope_mass = m_tot - core_mass;

  // removed (SPZ+GN:10 Nov 1998)
  //core_radius = helium_core_radius();
  //added (ST: 15 Sept 2009)
    core_radius = co_core_radius();
    
  detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
  effective_radius = max(radius, effective_radius);
}

// should only be used in constructors.
// (SPZ+GN:30 Sep 1998)
// ST:8 Sep 2009 can also be used in other functions.
void helium_giant::adjust_next_update_age() {
    real t_Heg = helium_giant_end_time(relative_helium_mass, get_total_mass());
    next_update_age = t_Heg;    
}

real helium_giant::helium_giant_end_time(const real mass, const real mass_tot) {
    
    //stop evolving star if the core mass reaches
    // max(m_Ch = cnsts.parameters(Chandrasekar_mass),
    // mc_ignite_CO = 0.773*relative_helium_mass - 0.35) or
    // min(get_total_mass(), 1.45*get_total_mass()-0.31)
    real mc_max = min(mass_tot, 1.45*mass_tot-0.31);
    mc_max = min(maximum_helium_giant_core_mass(mass), mc_max);
    
    //core should minimally grow 5% as a helium giant
    real t_hems = helium_main_sequence_time_for_solar_metalicity(mass);   
    real mco_hems = helium_giant_core_mass(t_hems, mass);  //co core mass
    mc_max = max(mc_max, 1.05*mco_hems);
    
    real t_Heg = helium_giant_age_core_mass_relation(mc_max, mass);
    real t_Hems  = helium_main_sequence_time_for_solar_metalicity(mass);
    if (t_Hems > t_Heg){
        cerr<<"In helium_giant_end_time t_Hems > t_Heg"<<endl;
        exit(-1);
    }
    return t_Heg;
}


void helium_giant::create_remnant(const real mass, const real mass_tot, const real mc_core) {

    stellar_type type;
    //real mc_SN = maximum_helium_giant_core_mass(mass);
    if (mass_tot < cnsts.parameters(Chandrasekar_mass)){
        // if mc_core equals 1.45*mass-0.31
        // shell burning stops before whole envelope is converted into C and O
        if (core_mass < mass_tot){
            if(!update_core_and_envelope_mass(get_total_mass())) {
                cerr << "Update core mass failed in helium_giant()"<<endl;
            }
        }
        
//      if (1.45*mass-0.31 < mass){
//            cerr<<"Warning: not homogeneous WD"<<endl;
//            type = Carbon_Dwarf;
//            
//            if(!update_core_and_envelope_mass(get_total_mass())) {
//                cerr << "Update core mass failed in helium_giant()"<<endl;
//            }
//        }


        // mc_core equals total mass
        // core mass reaches outside of star, no envelope anymore
        if (mass < 1.6 || mass_tot < 1.1)
            type = Carbon_Dwarf;
        else if (mass <= 2.25)
            type = Oxygen_Dwarf;
        else {
             type = Oxygen_Dwarf;
        }
    }
    else {    
        if (mass < 1.6) 
            type = Disintegrated;
        else {
    	  // (GN Oct  5 2016) Hurley uses theoretical Base AGB core mass to determine
    	  //     real mc_SN = maximum_helium_giant_core_mass(mass);
    	  //            if (mc_SN <= 7.)
    	  // in SeBa we use the current real properties
    	  // note core_mass = COcore_mass
            if (core_mass  <  cnsts.parameters(COcore2black_hole))
                type = Neutron_Star;
            else
                type = Black_Hole;
        }
    }
 
	switch (type) {
        case Black_Hole : star_transformation_story(Black_Hole);
                new black_hole(*this); 
                break;
        case Neutron_Star : star_transformation_story(Neutron_Star);
                new neutron_star(*this);
                break;
        case Disintegrated : star_transformation_story(Disintegrated);
			    new disintegrated(*this);
			    break;
        case Carbon_Dwarf : star_transformation_story(Carbon_Dwarf);
	            new white_dwarf(*this, Carbon_Dwarf);
                break;
        case Oxygen_Dwarf : star_transformation_story(Oxygen_Dwarf);
                new white_dwarf(*this, Oxygen_Dwarf);
                break;
            
        default :   cerr << "helium_giant::create_remnant()" <<endl;
                    cerr << "star_type not recognized." << endl;
                    exit(-1);
       }

}

//used by subtrac_mass_from_donor and double_star::perform_mass_transfer
real helium_giant::mdot_limit(const real dt, real mdot){
    //real mdot = relative_mass*dt/get_binary()->get_donor_timescale();
    if (is_binary_component()) {
	mdot = get_total_mass()*dt/get_binary()->get_donor_timescale();
    }
    return mass_ratio_mdot_limit(mdot);
    
}


star* helium_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

    mdot = mdot_limit(dt, mdot);
          
    if (mdot<envelope_mass){
        envelope_mass -= mdot;
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass
        adjust_next_update_age();
    }
      else {
        mdot = envelope_mass;
        envelope_mass = 0.;	
        // if (core_mass <= cnsts.parameters(minimum_helium_star)) {
        //    star_transformation_story(Carbon_Dwarf);
        //    return dynamic_cast(star*, new white_dwarf(*this));
        // }
          
        if (core_mass >= cnsts.parameters(Chandrasekar_mass)){
            real mc_SN = maximum_helium_giant_core_mass(relative_helium_mass);
            if (mc_SN < cnsts.parameters(COcore2black_hole)){
                star_transformation_story(Neutron_Star);
                return dynamic_cast(star*, new neutron_star(*this));
            }
            else{
                star_transformation_story(Black_Hole);
                return dynamic_cast(star*, new black_hole(*this));
            }
        }         
        else if(relative_helium_mass >= 1.6) {      
              star_transformation_story(Oxygen_Dwarf);
              return dynamic_cast(star*, new white_dwarf(*this, Oxygen_Dwarf));
        }
        else {
              star_transformation_story(Carbon_Dwarf);	   
              return dynamic_cast(star*, new white_dwarf(*this, Carbon_Dwarf));
        }          
      }
      return this;
}

star* helium_giant::reduce_mass(const real mdot) {   
    if (mdot < envelope_mass){
        envelope_mass -= mdot;
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass
        adjust_next_update_age();        
    }
    else {
 //       real rest_mdot = mdot;
//        rest_mdot -= envelope_mass;
//        envelope_mass = 0.;
//	   
//        if (core_mass>rest_mdot) {
//            if (core_mass-rest_mdot<=
//                cnsts.parameters(minimum_helium_star)) {
//                core_mass -= rest_mdot;
//                COcore_mass = core_mass;
//            }
//            else {
//                core_mass -= rest_mdot;
//                COcore_mass = core_mass;
//            }
//        }
//        else {
//            cerr << "ERROR!:"<<endl;
//            cerr << "void helium_giant::reduce_mass(mdot="
//                << rest_mdot<<")"<<endl;
//            cerr << "mdot exceeds helium core mass ("<<core_mass
//                << ")"<<endl;
//            cerr << "Decision: Disintegrate helium star!"<<endl;
//	  
//            star_transformation_story(Disintegrated);
//            return dynamic_cast(star*, new disintegrated(*this));
//        }
        envelope_mass = 0.;
        
        if (core_mass >= cnsts.parameters(Chandrasekar_mass)){
            real mc_SN = maximum_helium_giant_core_mass(relative_helium_mass);
            if (mc_SN < cnsts.parameters(COcore2black_hole)){
                star_transformation_story(Neutron_Star);
                return dynamic_cast(star*, new neutron_star(*this));
            }
            else{
                star_transformation_story(Black_Hole);
                return dynamic_cast(star*, new black_hole(*this));
            }
        }
        else if(relative_helium_mass < 1.6 || get_total_mass() < 1.1) {   
            star_transformation_story(Carbon_Dwarf);	   
            return dynamic_cast(star*, new white_dwarf(*this, Carbon_Dwarf));
        }
        else {         
            star_transformation_story(Oxygen_Dwarf);
            return dynamic_cast(star*, new white_dwarf(*this, Oxygen_Dwarf));
        }
    }
    return this;
}


//		general mass transfer utilities.
// //Increase donor mass and possibly relative_mass of donor.
// Check mass-transfer timescales before use.
real helium_giant::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {
    if (mdot<0) {
        cerr << "helium_giant::add_mass_to_accretor(mdot=" << mdot 
        << ", dt=" << dt << ")"<<endl;
        cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        return 0;
    }
    
    
    if(hydrogen){
        //hydrogen accretion
        // is treated in the same way as helium accretion..
        mdot = accretion_limit(mdot, dt);
        
        // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
        // adjust_accretor_age(mdot);
        envelope_mass += mdot;
        accreted_mass += mdot;
        if (accreted_mass > 0.05 * get_total_mass()){
            //if (is_binary_component()) cout<<get_binary()->get_identity();
	    cerr << "\t WARNING: accreted hydrogen mass more than 5% of helium giant"<<endl;
        }
	
	      
        // For now, rejuvenation of SG, CHeB, AGB or He giant accretor   
	// only if mtot > relative_mass
	if (relative_mass<get_total_mass())  {
	  update_relative_mass(get_total_mass());
	}	  

	
        // only neccessary for AGB & He giant accretor as  
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass      
        adjust_next_update_age();  
        
    }
    else{
        //for the moment assume helium accretion
        // for the moment no adjust_accretor_radius
        mdot = accretion_limit(mdot, dt);
        
        // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
        //adjust_accretor_age(mdot);
        envelope_mass += mdot;
        
      
        // For now, rejuvenation of SG, CHeB, AGB or He giant accretor   
	// only if mtot > relative_mass
	if (relative_mass<get_total_mass())  {
	  update_relative_mass(get_total_mass());
	}	  
	
        // only neccessary for AGB & He giant accretor as  
        // next_update_age is a function of total mass
        // as the maximal core mass can be the total mass
        // when total mass < chandrasekhar mass      
        adjust_next_update_age();  
        
    }
    set_spec_type(Accreting);
    return mdot;
}


real helium_giant::accretion_limit(const real mdot, const real dt) {
  // needed in double_star::zeta(donor, accretor) 
  
	return accretion_limit_eddington(mdot, dt);
}

# if 0
// currently not used
// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void helium_giant::adjust_accretor_age(const real mdot,
				       const bool rejuvenate) {
    
    real m_he_rel_new;
    real m_tot_new = get_total_mass() + mdot;
    if (m_tot_new>relative_helium_mass)
        m_he_rel_new = m_tot_new;
    else m_he_rel_new = relative_helium_mass;
    
    real t_hems_old = helium_main_sequence_time_for_solar_metalicity(relative_helium_mass); 
    real dt_Heg_old = helium_giant_end_time(relative_helium_mass, get_total_mass()) - t_hems_old;

    real t_hems_new = helium_main_sequence_time_for_solar_metalicity(m_he_rel_new);    
    //He star+giant not a function of metalicity for now
    real t_Heg_new = helium_giant_end_time(m_he_rel_new, m_tot_new);
    real dt_Heg_new = t_Heg_new - t_hems_new; 
    
    real dtime = relative_age - t_hems_old;
    
    // (GN+SPZ May  4 1999) update last_update_age
    last_update_age = t_hems_new;
    relative_age = t_hems_new 
    + dtime*(dt_Heg_new/dt_Heg_old);
    
    if (rejuvenate)
        relative_age *= rejuvenation_fraction(mdot/m_tot_new); 
    
    if (relative_age < last_update_age + cnsts.safety(minimum_timestep)){
        cerr<<"In helium_giant::adjust_accretor_age relative age updated on HeGiant, but < last_update_age"<<endl;
    }
    relative_age = max(relative_age, 
                       last_update_age + cnsts.safety(minimum_timestep)); 
    relative_age = min(relative_age, t_Heg_new);

    // next_update_age should not be reset here
    // next_update_age = t_ms_new + t_hg_new;
    
//    real frac = (1-pow(mdot/(get_total_mass()+mdot),
//			     cnsts.parameters(rejuvenation_exponent)));
//      last_update_age *= frac;
//      relative_age *= frac;
}
#endif



//Eq.85
// determine whether the star is a helium giant or helium hertzsprung gap star
 real helium_giant::helium_giant_type(const real lum, const real mass, const real mass_tot, const real z){
    
    real lambda = 500*(2.0+pow(mass_tot,5))/pow(mass_tot, 2.5);
    real Rzhe = helium_star_radius_for_solar_metalicity(mass_tot);
    real L_tHems = terminal_helium_main_sequence_luminosity(mass); 
    real R1 = Rzhe * pow(lum/L_tHems, 0.2) + 0.02*(exp(lum/lambda)-exp(L_tHems/lambda));//He Hg
    real R2 = 0.08*pow(lum, 0.75); //He GB
 
    if (R1 < R2) return 5.;
    else return 6;
}
    


// see helium_star.C
real helium_giant::zeta_adiabatic() {

     real heg_type = helium_giant_type(luminosity, relative_helium_mass, get_total_mass(), metalicity);
     if (heg_type ==  5){ // helium hertzsprung 
     	return 4;
	}
     else { // helium giant
	real z = 0.;
//      Hjellming and Webbink 1987 ApJ, 318, 804
     	real x = core_mass/get_total_mass();
     	real A = -0.220823;
     	real B = -2.84699;
     	real C = 32.0344;
     	real D = -75.6863;
     	real E = 57.8109;
 	z = A + x*(B + x*(C + x*(D + x*E)));
     	return z;
      }
}

real helium_giant::zeta_thermal() {

//	real z = -2;
//	real z = 0; // with this helium giants transfer mass on the nuclear timescale
      real z = -2;
//      if (get_core_mass()<=0.4)  // like a white dwarf
//	z = -cnsts.mathematics(one_third);
	return z;
}

#if 0
real helium_giant::CO_core_mass() {
      // C/O core of helium star grows linearly with time
      // (SPZ+GN:26 Sep 1998)

// (GN+SPZ May  4 1999) core groth totally in helium_star phase: core constant
//  real m_core = get_total_mass()*relative_age/next_update_age;
//  m_core = max(core_mass,m_core);
//    return min(m_core, get_total_mass());

  return min(core_mass, get_total_mass());
}
#endif

#if 0
// end time is a function of the envelope mass, which is a function as the wind
// prescription is dependent on the number of steps
void helium_giant::stellar_wind(const real dt) {

//  PRL(last_update_age);
//  PRL(next_update_age);
//  PRL(relative_age);
//  PRL(previous.relative_age);
//  PRL(dt);
// (GN+SPZ Apr 28 1999) wind for low mass stars per phase
    real end_time = next_update_age - last_update_age;
    real relative_time = min(relative_age - last_update_age, end_time);

//    PRL(end_time);
//    PRL(relative_time);
//    PRL(wind_constant);
    real wind_mass = 0;
    if (relative_time - dt > 0){
    wind_mass = wind_constant 
                   * (pow(relative_time/end_time,
			cnsts.parameters(massive_star_mass_loss_law))
	           -  pow((relative_time-dt)/end_time,
			cnsts.parameters(massive_star_mass_loss_law)));
    }
// (GN+SPZ May  6 1999) try low wind: 
//    wind_mass = 0.;

//    PRL(wind_mass);
//    PRL(envelope_mass);

  if (wind_mass>=envelope_mass) {
    wind_mass = envelope_mass;
    effective_radius = radius = core_radius;
  }

  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else
    reduce_mass(wind_mass);
  return;
}
#endif




real helium_giant::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}


//absidal motion constant
real helium_giant::amc() {

// (SilT 13 Feb 22) 
//    based on Brooke & Olle 1955, for n=3 polytrope
//    return 0.143 
//    based on Claret & Gimenez 1992, 96, 225 the value should be smaller, try:
    return 0.05;
}


// Section 7.2 in Hurley, Pols & Tout 2000
real helium_giant::convective_envelope_mass(){
    return envelope_mass;
}

// Section 2.3.1 in Hurley, Tout & Pols 2002
real helium_giant::convective_envelope_radius(){
    return max(0., radius - co_core_radius()); // function of core_mass
}


void helium_giant::update_wind_constant() {
 
 // (GN+SPZ May  3 1999) helium_giants loose complete envelope in wind-
 //  wind_constant = (1 - cnsts.parameters(helium_star_final_core_fraction))
 //                * get_total_mass(); 
 
 // (GN+SPZ May  7 1999) envelope is about 30% of total mass,
 // we loose 10% of total mass ....
//  wind_constant = 0.3*envelope_mass;
     
// (SilT June 15 2012) based on HPT2000
//    Reimers 1975
//     GB like stars
    real neta = 0.5; 
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();
    // (SilT Jan 2020) metalicity dependence on average with (Z/Z_sun)^0.85
    // Vink & de Koter 2005
    dm_r *= pow(metalicity/cnsts.parameters(solar_metalicity),0.85);

    //HPT2000
    //Reduced WR-like mass loss for small H-envelope mass
    real dm_wr = 1.E-13 * pow(luminosity, 1.5);
    // (SilT Jan 2020) metalicity dependence on average with (Z/Z_sun)^0.85
    // Vink & de Koter 2005
    dm_wr *= pow(metalicity/cnsts.parameters(solar_metalicity),0.85);

    //(AD: 21 Apr 2022) added WR wind-prescription from Sanders & Vink 2020 (2009.01849)
//    real dm_wr = 0;
//    real alpha_sv = 0.32 * log10(metalicity/cnsts.parameters(solar_metalicity)) + 1.40;
//    real loglum10 = -0.87 * log10(metalicity/cnsts.parameters(solar_metalicity)) + 5.06;
//    real logdm10 = -0.75 * log10(metalicity/cnsts.parameters(solar_metalicity)) - 4.06;
//    real lum10 = pow(10, loglum10);
//    real dm10 = pow(10, logdm10);
//    if (luminosity > lum10){
//        dm_wr = dm10 * pow((log10(luminosity/lum10)), alpha_sv) * pow((luminosity/(10*lum10)), 0.75);
//        dm_wr = max(dm_wr,1.E-8);
//    } else {
//        dm_wr = 1.E-8;   
//    }

    
    wind_constant = max(max(dm_wr, dm_r), 0.0);   
	
}


stellar_type helium_giant::get_element_type() {

     stellar_type type = Helium_Giant;
//     if (envelope_mass <= 0)
//	 type = Carbon_Star;

     return type;
     }

real helium_giant::temperature() {

  real T_eff = cnsts.parameters(Tsun)
             * sqrt(sqrt(luminosity)/effective_radius);
  
  return T_eff;

  //return sqrt(33.45*sqrt(luminosity)/effective_radius);
}

//Eq. 84
real helium_giant::helium_giant_luminosity_core_mass_relation(const real time, const real mass, const real z){
    real A_He = AGB_A_He_estimator();
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real l_x = helium_giant_x_luminosity(mass);
    real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
    real l_He;


    if (time  <= t_x){
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        real arg = (p-1)*A_He*D*(t_inf1-time);   
        l_He = D * pow(arg, p/(1-p));      
    }
    else {        
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real t_inf2 = specific_time_limit(A_He, t_x,
                                          B, l_x, q);
        real arg = (q-1)*A_He*B*(t_inf2-time);   
        l_He = B * pow(arg, q/(1-q));
    }  

    return l_He; 
}


real helium_giant::helium_giant_age_core_mass_relation(const real m_core, const real mass){
    real age;  
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real A_He    = AGB_A_He_estimator();
    real t_Hems  = helium_main_sequence_time_for_solar_metalicity(mass);
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    
    real mx = helium_giant_x_mass(mass);
    if (m_core <= mx){        
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        age = t_inf1 - pow(m_core, 1.-p)/A_He/D/(p-1.);
    }
    else{
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real l_x = helium_giant_x_luminosity(mass);
        real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
        real t_inf2 = specific_time_limit(A_He, t_x,
                                           B, l_x, q);
        age = t_inf2 - pow(m_core, 1.-q)/A_He/B/(q-1.);
    }
    return age;
}


void helium_giant::evolve_core_mass(const real time,
                    const real mass, const real mass_tot) {

    real mco_core = helium_giant_core_mass(time, mass);    
    // (GN Oct 26 2010) for mergers
    if (mco_core < core_mass) mco_core = core_mass;

        
    if(!update_COcore_mass(mco_core)) {
        cerr << "Update COcore mass failed in helium_giant()"<<endl;
    }
    if(!update_core_and_envelope_mass(mco_core)) {
        cerr << "Update core mass failed in helium_giant()"<<endl;
    }
}

void helium_giant::evolve_core_mass() {
    evolve_core_mass(relative_age, relative_helium_mass, get_total_mass());
}

//related to Eq. 84
real helium_giant::helium_giant_core_mass(const real time, const real mass){
    
    real A_He = AGB_A_He_estimator();
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real l_x = helium_giant_x_luminosity(mass);
    real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
    real m_core;

    
    if (time  <= t_x){
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        real arg = (p-1)*A_He*D*(t_inf1-time);   
        m_core = pow(arg, 1.0/(1-p));     
    }
    else {
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real t_inf2 = specific_time_limit(A_He, t_x,
                                          B, l_x, q);
        real arg = (q-1)*A_He*B*(t_inf2-time);   
        m_core = pow(arg, 1.0/(1-q));      
    }  

    return m_core;     
}


// Eq.75
real helium_giant::maximum_helium_giant_core_mass(const real mass) {
    real m_Ch = cnsts.parameters(Chandrasekar_mass);
    real mc_max = max(m_Ch, 0.773*mass - 0.35);
    return mc_max;
}


void helium_giant::small_envelope_perturbation(){
    real mu = small_envelope_mu(luminosity, get_total_mass(), core_mass);
    if(mu < 1.){
        real lum_c = small_envelope_core_luminosity();
        luminosity = perturb_luminosity(luminosity, lum_c, get_total_mass(), core_mass, mu);
        real rad_c = small_envelope_core_radius();
        if(rad_c < radius){
            radius = perturb_radius(radius, rad_c, get_total_mass(), core_mass, mu);
        }
    }
}


//Eq.98
real helium_giant::small_envelope_mu(const real lum, const real mass_tot, const real m_core){
    // the function parameter lum is not needed here, 
    // but needed in single_star::small_envelope_mu
    // in this way the correct function is called automatically
    
    real mc_max = min(mass_tot, 1.45*mass_tot-0.31);    
    real mu = 5.*(mc_max-m_core) / mc_max;
    mu = Starlab::max(mu,0.);
    return mu;
}

real helium_giant::co_core_radius(const real m_core){
    // due to small nucleair burning layer 
    // r_c > white_dwarf_radius
    //real r_c = 5.*white_dwarf_radius(m_core, 10000.);
    // (GN+SilT Feb 24 2022) NOTE THIS NEEDS UPDATING FOR NON-DEGENERATE CORES
    real r_c = (2.7 - 1.129*m_core)*white_dwarf_radius(m_core, 10000.);
  
    return r_c;
}

real helium_giant::co_core_radius(){
    return co_core_radius(core_mass);
}


real helium_giant::small_envelope_core_radius(const real m_core){
    real r_c = white_dwarf_radius(m_core, 10000. );
    return r_c;
}

real helium_giant::small_envelope_core_radius(){
    return small_envelope_core_radius(core_mass);
}

real helium_giant::small_envelope_core_luminosity(){
    real l_c = 40.;
    return l_c;
}    





