//
// helium_star.C
//

#include "helium_star.h"
#include "main_sequence.h"
#include "hyper_giant.h"
#include "sub_giant.h"
#include "horizontal_branch.h"
#include "super_giant.h"
#include "hertzsprung_gap.h"

// (AD Oct 4 2022)
helium_star::helium_star(main_sequence & h) : single_star(h) {

    delete &h;
    real m_tot = get_total_mass();
    core_mass = m_tot;
    COcore_mass = 0.;
    envelope_mass = 0;
//    lose_envelope_decent(); //assume all helium has been converted to hydrogen 
    accreted_mass = 0;
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age units
    last_update_age = 0.;
    adjust_next_update_age();
    relative_age = 0.;
    
    instantaneous_element();
    update();
    
    post_constructor();

}

helium_star::helium_star(hertzsprung_gap & h) : single_star(h) {

    delete &h;
    lose_envelope_decent();
    accreted_mass = 0;
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age units
    last_update_age = 0.;
    adjust_next_update_age();
    relative_age = 0.;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = COcore_mass = 0.;
    //final_core_mass = final_CO_core_mass(m_tot);
    //core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();
    
    instantaneous_element();
    update();
    
    post_constructor();

}

helium_star::helium_star(sub_giant & g) : single_star(g) {

    delete &g;
    lose_envelope_decent();
    accreted_mass = 0;

    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age units
    last_update_age = 0.;
    adjust_next_update_age();
    relative_age = 0.;

    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = COcore_mass = 0.;
    //final_core_mass = final_CO_core_mass(m_tot);
    //core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();

    instantaneous_element();
    update();
    
    post_constructor();

}

helium_star::helium_star(horizontal_branch & h) : single_star(h) {

    delete &h;
    lose_envelope_decent();
    accreted_mass = 0;

    relative_age = (relative_age - helium_ignition_time(relative_mass, metalicity)) / core_helium_burning_timescale(relative_mass, metalicity)
     * helium_main_sequence_time_for_solar_metalicity(get_total_mass());
    // if relative_mass is updated, then after previous line
    
    last_update_age = 0.;//relative_age;
    
    // (GN+SPZ May  4 1999) last update age is time of previous type change
    // in NEW relative_age uits
    adjust_next_update_age();
    //relative_age = t_frac* next_update_age;
    
    // core_mass is CO core (SPZ+GN: 27 Jul 2000)
    real m_tot = get_total_mass();
    // (SPZ+GN: 28 Jul 2000) was helium core but is now the CO core.
    core_mass = COcore_mass = 0.;
    //final_core_mass = final_CO_core_mass(m_tot);
    //core_mass = COcore_mass = CO_core_mass();
    envelope_mass =  m_tot - core_mass;

// (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
//    update_wind_constant();

    instantaneous_element();
    update();
    
    post_constructor();

 }

#if 0
void helium_star::adjust_initial_star() {

  update_wind_constant();

  if(relative_age<=0)
    relative_age = max(current_time, 0.0);
}
#endif

void helium_star::adjust_next_update_age() {
//    for a helium star relative_helium_mass = get_total_mass()  
    next_update_age = helium_main_sequence_time_for_solar_metalicity(get_total_mass());
}

// Adjust radius & luminosity at relative_age
void helium_star::instantaneous_element() {
//    for a helium star relative_helium_mass = get_total_mass()  
    luminosity  = helium_main_sequence_luminosity(relative_age, get_total_mass());    
    radius      = helium_main_sequence_radius(relative_age, get_total_mass(), get_total_mass());

    if (accreted_mass > 0.05 * get_total_mass()){
        radius *= 50;    
    }

    //effective_radius = max(effective_radius, radius);
    effective_radius = radius;
}

// Evolve a helium_star upto time argument according to
// the new 2000 models.
void helium_star::evolve_element(const real end_time) {
    real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;

//    for a helium star relative_helium_mass = get_total_mass()  
    if (get_total_mass() < cnsts.parameters(minimum_helium_star)) {   
        // Helium Main_sequence star will not ignite core helium burning.
        cerr<<"Warning: not homogeneous WD"<<endl;
        if(!update_core_and_envelope_mass(get_total_mass())) {
            cerr << "Update core mass failed in helium_star()"<<endl;
        }
        
        star_transformation_story(Helium_Dwarf);
        new white_dwarf(*this, Helium_Dwarf);
        return;
    }
    
    if (relative_age <= next_update_age) {
        instantaneous_element();
    } 
    else {
        star_transformation_story(Helium_Giant);
        new helium_giant(*this);
        return;
    }

    update();
    stellar_wind(dt);
}



// (SilT Dec 10 2020) 
real helium_star::get_evolve_timestep() {

  real timestep = min((next_update_age - last_update_age )/ cnsts.safety(number_of_steps), 
                      next_update_age - relative_age - 0.5 * cnsts.safety(minimum_timestep));   
   
//     (SilT Dec 10 2020) extra safety measure for large L
    real dt_mdot = timestep * pow(10000/luminosity, 0.75) / pow(metalicity/cnsts.parameters(solar_metalicity), 0.85) ;    

   return max(min(timestep, dt_mdot), cnsts.safety(minimum_timestep));
   
}






real helium_star::nucleair_evolution_timescale() {    
//    for a helium star relative_helium_mass = get_total_mass()     
    return helium_main_sequence_time_for_solar_metalicity(get_total_mass());

}



void helium_star::update() {
    
    detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
    effective_radius = radius;

}

#if 0
void helium_star::create_remnant() {

/*
//cout<<"void helium_star::create_remnant() called"<<endl;
       // core_mass = get_total_mass();
        stellar_type type = NAS;
        if (get_total_mass()>=cnsts.parameters(helium2black_hole))
           type = Black_Hole;
//			was 1.3*M_NS
        else if((core_mass>=1.2*cnsts.mass.neutron_star     &&
                 relative_mass>=10) ||
                get_total_mass()>=cnsts.mass.helium2neutron_star)
           type = Neutron_Star;
//        else if(get_total_mass()>=1.3*M_NS && relative_mass>=8) 
//           type = Disintegrated;
        else
           type = White_Dwarf;

//        base_element * s = make(type, this);
       switch (type) {
          case Black_Hole : if (is_binary_component())
	                       get_binary()->dump("sn.dat");
                            else dump("sn.dat");
			    
			    star_transformation_story(this,Black_Hole);
                            new black_hole(*this); 
                            break;
          case Neutron_Star : 
                           if (is_binary_component()) get_binary()->dump("sn.dat");
                           else dump("sn.dat");
			   star_transformation_story(this,Neutron_Star);
			   new neutron_star(*this);
                           break;
          case Disintegrated : if (is_binary_component()) get_binary()->dump("sn.dat");
                               else dump("sn.dat");
			       star_transformation_story(this,Disintegrated);
                               new disintegrated(*this);
                               break;
         case White_Dwarf : star_transformation_story(this,White_Dwarf);
	                    new white_dwarf(*this);
			    break;
          default : ;//error << "helium_star::create_remnant()" <<endl;
                    //error << "star_type not recognized." << endl;
                    //error << "Action: Proceed normally." << endl;
       }

*/
     }
#endif


//used by subtrac_mass_from_donor and double_star::perform_mass_transfer
real helium_star::mdot_limit(const real dt, real mdot){
    //real mdot = relative_mass*dt/get_binary()->get_donor_timescale();
    if (is_binary_component()) {
      mdot = get_total_mass()*dt/get_binary()->get_donor_timescale();
    }
    return mass_ratio_mdot_limit(mdot);
    
}

star* helium_star::subtrac_mass_from_donor(const real dt, real& mdot) {
//    for a helium star relative_helium_mass = get_total_mass()  
    
    mdot = mdot_limit(dt, mdot);
    
    //if (envelope_mass >= mdot)
    //  envelope_mass -= mdot;
    //else {
    //  mdot = envelope_mass;
    //  envelope_mass = 0.;
    //  
    //  //	   star_transformation_story(this,White_Dwarf);
    //  //return dynamic_cast(star*, new white_dwarf(*this));
    //}

    adjust_age_after_mass_loss(mdot, true);
    envelope_mass -= mdot;
    adjust_next_update_age();
    
    if (get_total_mass() < cnsts.parameters(minimum_helium_star)) {
        // Helium Main_sequence star will not continue core helium burning.
        cerr<<"Warning: not homogeneous WD"<<endl;
        if(!update_core_and_envelope_mass(get_total_mass())) {
            cerr << "Update core mass failed in helium_star()"<<endl;
        }
        
        star_transformation_story(Helium_Dwarf);
        return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
    }
    return this;
}

star* helium_star::reduce_mass(const real mdot) {
    adjust_age_after_mass_loss(mdot, true);
    envelope_mass -= mdot;
    adjust_next_update_age();
    
    if (get_total_mass() < cnsts.parameters(minimum_helium_star)) {
        // Helium Main_sequence star will not continue core helium burning.
        cerr<<"Warning: not homogeneous WD"<<endl;
       if(!update_core_and_envelope_mass(get_total_mass())) {
           cerr << "Update core mass failed in helium_star()"<<endl;
       }
        star_transformation_story(Helium_Dwarf);
        return dynamic_cast(star*, new white_dwarf(*this, Helium_Dwarf));
    }
    
    //   On the He-MS there is no defined core yet, 
    //   no CO star/WD can form yet 
    //   if (envelope_mass<=mdot) {
    //         envelope_mass = 0.;
    //         //star_transformation_story(Helium_Star);
    //         //return dynamic_cast(star*, new helium_star(*this));
    //      }
    
    
//    if (envelope_mass >= mdot)
//        envelope_mass -= mdot;
//    else {
//      cerr<<"In helium_star::reduce_mass: what does this part do? only for constructing he star?"<<endl;
//      real mass_reduced = mdot;
//      mass_reduced -= envelope_mass;
//      envelope_mass = 0.;
//      if (core_mass>mass_reduced) {
//        if (core_mass-mass_reduced<=cnsts.parameters(minimum_helium_star)) {
//            core_mass -= mass_reduced;
//    //		 star_transformation_story(this,White_Dwarf);
//    //                 new white_dwarf(*this);
//    //                 return;
//	    }
//        else {
//          core_mass -= mass_reduced;
//          COcore_mass = core_mass;
//        }
//      }
//      else {
//        cerr<<"ERROR!:"<<endl;
//        cerr<<"void helium_star::reduce_mass(mdot="
//            <<mass_reduced<<")"<<endl;
//        cerr<<"mdot exceeds helium core mass ("<<core_mass
//            <<")"<<endl;
//        cerr<<"Decision: Disintegrate helium star!"<<endl;
//        //new disintegrated(*this);
//        //return;
//       }
//    }
    
    return this;
}


// Age adjustment especially for (wind) mass loss
// It is part of the single star evolution, 
// so it can include information from tracks
void helium_star::adjust_age_after_mass_loss(const real mdot,
                                                    const bool rejuvenate=true) {
//    for a helium star relative_helium_mass = get_total_mass()      
    real m_tot_new = get_total_mass() - mdot;
    
    real t_hems_old = helium_main_sequence_time_for_solar_metalicity(get_total_mass()); 
    real t_hems_new = helium_main_sequence_time_for_solar_metalicity(m_tot_new);
    
    relative_age *= (t_hems_new/t_hems_old);
//    if (rejuvenate){
//        real mdot_fr = -1.0 * mdot/m_tot_new;
//        real rejuvenation = (1-pow(mdot_fr,
//                                   cnsts.parameters(rejuvenation_exponent)));
//        relative_age *= rejuvenation; 
//    }
    relative_age = min(relative_age, t_hems_new);

    // next_update_age should not be reset here,
    // is done in add_mass_to_accretor, where also relative_mass
    // is updated
    // next_update_age = t_ms_new; 
}

// add mass to accretor
// is a separate function (see single_star.C) because rejuvenation
real helium_star::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {

    if (mdot<0) {
       cerr << "helium_star::add_mass_to_accretor(mdot=" 
             << mdot << ")"<<endl;
       cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
       return 0;
    }
    
    if(hydrogen){
        //hydrogen accretion
        // is treated in the same way as helium accretion..
        mdot = accretion_limit(mdot, dt);

        adjust_accretor_age(mdot, true);
        // (GN+SPZ May  3 1999) Langer wind: see helium_star::stellar_wind
        //	update_wind_constant();
        envelope_mass += mdot;
        accreted_mass += mdot;
        if (accreted_mass > 0.05 * get_total_mass()){
            //if (is_binary_component()) cout<<get_binary()->get_identity();
	    cerr << "\t WARNING: accreted hydrogen mass more than 5% of helium star"<<endl;
        }
	
        adjust_next_update_age();

//        possible track hydrogen accreting helium star can turn into horizontal branch star
//        not implemented currently
//        core_mass = get_total_mass();
//        envelope_mass = mdot;
//        star_transformation_story(Horizontal_Branch);
//        new horizontal_branch(*this);
//        return mdot;
    }
    else{
        //for the moment assume helium accretion
        // for the moment no adjust_accretor_radius
        mdot = accretion_limit(mdot, dt);
        adjust_accretor_age(mdot, true);
        envelope_mass += mdot;
        adjust_next_update_age();
    }
    set_spec_type(Accreting);
    return mdot;

 }

real helium_star::accretion_limit(const real mdot, const real dt) {
//needed in double_star::zeta(donor, accretor)

	return accretion_limit_eddington(mdot, dt);
}


// Star is rejuvenated by accretion.
// Age adjustment especially for accretion from other stars.
// No information from stellar evolution tracks is included.
void helium_star::adjust_accretor_age(const real mdot,
				      const bool rejuvenate) {
//    for a helium star relative_helium_mass = get_total_mass()  
    
    real m_tot_new = get_total_mass() + mdot;
    real t_hems_old = helium_main_sequence_time_for_solar_metalicity(get_total_mass()); 
    real t_hems_new = helium_main_sequence_time_for_solar_metalicity(m_tot_new);
    
    relative_age *= (t_hems_new/t_hems_old);
    if (rejuvenate)
        relative_age *= rejuvenation_fraction(mdot/m_tot_new); 
    
    relative_age = min(relative_age, t_hems_new);

    
    
    // next_update_age should not be reset here,
    // is done in add_mass_to_accretor, where also relative_mass
    // is updated
    // next_update_age = t_ms_new; 
}



// Helium stars have radiative envelopes.
// this requires more research (SPZ+GN: 3 Oct 1998)
real helium_star::zeta_adiabatic() {

  //real z = 15;
  real z = 4;
  
  return z;

}

real helium_star::zeta_thermal() {

// (GN+SPZ May  3 1999) see evolve_element; mass - radius relation
  if (get_total_mass() < 0.2 ) {
    return -0.19;
  } else {
    return 1;
  }
}
#if 0
real helium_star::final_CO_core_mass(const real initial_mass) {

  // Implementation of Nelemans YPZV 2000 (A&A submitted)
  // Itroduces discontinuity at relative_mass = 2.2
  // bases on Habets 1986 & IT85
  // (SPZ+GN: 27 Jul 2000)
  real final_coremass_fraction;
  if(initial_mass <= 0.8) 
    final_coremass_fraction = 1;
  else if(initial_mass >= cnsts.parameters(helium2neutron_star)) 
    final_coremass_fraction = 0.65;
  else 
    final_coremass_fraction = 1 - 0.32 * (initial_mass - 0.8);

  return final_coremass_fraction*initial_mass;
}
#endif
#if 0
real helium_star::CO_core_mass() {

  real m_core = final_core_mass * relative_age/next_update_age;
  m_core = max(core_mass, m_core);

  return min(m_core, get_total_mass());
}
#endif

#if 0
// end time is a function of the envelope mass, which is a function as the wind
// prescription is dependent on the number of steps
void helium_star::stellar_wind(const real dt) {

//cerr<<"void helium_star::stellar_wind(dt="<<dt<<") = ";
#if 0
	    real kappa = pow(get_total_mass(),2.5);

//	real kappa = wind_constant;
//        real wmlc = wind_constant;
//        wind_constant = kappa;

//		WR Wind
//		Langer, N., 1989, AA 220, 135.
        if(get_total_mass()>core_mass) {
//           kappa = 0.5*(kappa + wmlc);
//           real wind_mass = 0.005*dt*kappa;
// (GN+SPZ May  3 1999) Lager wind (again)
// Langer: M_dot = 5-10 10^-8 M^2.5

	  real wind_mass = 0.05*dt*kappa;
// (GN+SPZ May 12 1999) proper integration
	  real m = get_total_mass();
	  real constant = 0.05;

	  real m_next = pow((pow(m,-1.5) + 1.5*constant*dt),-1/1.5);
	  wind_mass = m - m_next;
#endif

            // (GN Apr 12 2001) Helium stars wind according to
        // Nelemans \& van den Heuvel), PhD thesis, p. 120, 
        real m = get_total_mass();
        real alph = 2.87;
        real constant = 0.0138;
        
        real m_next = pow((pow(m, (1. - alph)) 
                           + (alph-1.) * constant*dt),1./(1. - alph));
        real wind_mass = m - m_next;
        
	// (GN Nov 2012) metalicity dependence on average with sqrt(Z/Z_sun)
	// Vink & de Koter 2005
	wind_mass *= pow(metalicity/cnsts.parameters(solar_metalicity),0.5);

            
            
//	  PRC(m);PRC(m_next);PRL(wind_mass);
//           real wind_mass = 0.0025*dt*kappa;
// (GN+SPZ May  3 1999) low mass helium stars have no wind
// helium_giant looses envlope
	  if (get_total_mass() < 2.5 ) wind_mass = 0.;

#if 0
            if (wind_mass>=envelope_mass) {
              wind_mass = envelope_mass;
              radius = core_radius;
           }
#endif

           if (is_binary_component())
              get_binary()->adjust_binary_after_wind_loss(
                          this, wind_mass, dt);
           else
              reduce_mass(wind_mass);
              return;
        //      }
   } 
#endif

real helium_star::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}



//absidal motion constant
real helium_star::amc() {

// (SilT 13 Feb 22) 
//    based on Brooke & Olle 1955, for n=3 polytrope
    return 0.0144;           

}


void helium_star::update_wind_constant() {
// (GN+SPZ May  3 1999) NOT USED!, but called by post_constructor // (SilT June 15th 2012) Used again!
//  wind_constant = (1 - cnsts.parameters(helium_star_final_core_fraction)) 
//                      * get_total_mass();
 
// (SilT June 15 2012) based on HPT2000
    // Reimers 1975
    // GB like stars
    real neta = 0.5; 
    real dm_r = neta * 4.E-13 * radius * luminosity / get_total_mass();
    dm_r *= pow(metalicity/cnsts.parameters(solar_metalicity),0.85);	
    
    //HPT2000
    //Reduced WR-like mass loss for small H-envelope mass
    real dm_wr = 1.E-13 * pow(luminosity, 1.5);
    // (SilT Jan 2020) metallicity dependence on average with (Z/Z_sun)^0.85
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



stellar_type helium_star::get_element_type() {
//  if (envelope_mass <= 0) 
//    return Carbon_Star;
//  else
    return Helium_Star;
}

real helium_star::temperature() {
  real T_eff = cnsts.parameters(Tsun)
             * sqrt(sqrt(luminosity)/effective_radius);
  return T_eff;
  //return sqrt(33.45*sqrt(luminosity)/effective_radius);
}
