//
// white_dwarf.C
//

#include "white_dwarf.h"
#include "hertzsprung_gap.h"
#include "sub_giant.h"
#include "super_giant.h"
#include "helium_star.h"
#include "helium_giant.h"

white_dwarf::white_dwarf(super_giant & g, stellar_type wd_type) : single_star(g) {
      delete &g;
      white_dwarf_type = wd_type;
    
      real m_tot    = get_total_mass();
//    (GN+SilT Mar 2010) was cnsts.parameters(kanonical_neutron_star_mass)
      core_mass     = min(0.999999*cnsts.parameters(Chandrasekar_mass),
			  core_mass); 

      envelope_mass = m_tot - core_mass;
      accreted_mass = 0;

      lose_envelope_decent();

// (GN+SPZ May  4 1999) last update age is time of previous type change
//      last_update_age = next_update_age;
    last_update_age = 0.;
    relative_age = 0.;//1 + nucleair_evolution_time();

      instantaneous_element();
      update();

      post_constructor();
}

white_dwarf::white_dwarf(sub_giant & s, stellar_type wd_type) : single_star(s) {
      delete &s;
      white_dwarf_type = wd_type;

      real m_tot    = get_total_mass();
      core_mass     = min(0.999999*cnsts.parameters(Chandrasekar_mass),
			  core_mass); 
      envelope_mass = m_tot - core_mass;
      accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
//      last_update_age = next_update_age;
    last_update_age = 0.;

      lose_envelope_decent();

    relative_age = 0.;//1 + nucleair_evolution_time();

      instantaneous_element();
      update();

      post_constructor();

}

white_dwarf::white_dwarf(hertzsprung_gap & s, stellar_type wd_type) : single_star(s) {
      delete &s;
      white_dwarf_type = wd_type;

      real m_tot    = get_total_mass();
      core_mass     = min(0.999999*cnsts.parameters(Chandrasekar_mass),
			  core_mass); 
      envelope_mass = m_tot - core_mass;
      accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
//      last_update_age = next_update_age;
    last_update_age = 0.;

      lose_envelope_decent();

    relative_age = 0.;//1 + nucleair_evolution_time();

      instantaneous_element();
      update();

      post_constructor();

}


white_dwarf::white_dwarf(helium_star & h, stellar_type wd_type) : single_star(h) {
 
        delete &h;
    white_dwarf_type = wd_type;

	real m_tot    = get_total_mass();
	core_mass     = min(0.999999*cnsts.parameters(Chandrasekar_mass),
			    core_mass); 
	envelope_mass = m_tot - core_mass;
	accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
//        last_update_age = next_update_age;
    last_update_age = 0.;

	lose_envelope_decent();
	
    relative_age = 0.;//1 + nucleair_evolution_time();

	instantaneous_element();
	update();

	post_constructor();

}


white_dwarf::white_dwarf(helium_giant & h, stellar_type wd_type) :  single_star(h) {

        delete &h;
    white_dwarf_type = wd_type;

	real m_tot    = get_total_mass();
	core_mass     = min(0.999999*cnsts.parameters(Chandrasekar_mass),
			    core_mass); 
	envelope_mass = m_tot - core_mass;
	accreted_mass = 0;

// (GN+SPZ May  4 1999) last update age is time of previous type change
//        last_update_age = next_update_age;
    last_update_age = 0.;

	lose_envelope_decent();

    relative_age = 0.;//1 + nucleair_evolution_time();

	instantaneous_element();
	update();

	post_constructor();

}
      
void white_dwarf::instantaneous_element() {
//cerr << "white_dwarf::instantaneous_element"<<endl;

        next_update_age = relative_age + cnsts.safety(maximum_timestep);

	luminosity = 40;
	// m_rel may not become larger than M_Ch
//	real m_rel = min(0.999999, core_mass/cnsts.parameters(Chandrasekar_mass));

	// Nauenberg, M, 1972, Apj 175, 417
	// mu=4.039 is the mean molecular weight.
	real mu = 2.;

	effective_radius =
	core_radius      =
	radius           = white_dwarf_radius(core_mass, relative_age);
	  //               min(0.1, (0.0225/mu)*
	  //               sqrt(1./pow(m_rel, cnsts.mathematics(two_third))
	  //	           - pow(m_rel, cnsts.mathematics(two_third))));
}       


void white_dwarf::evolve_element(const real end_time) {
         real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

// done in add_mass_to_accretor
//        accrete_from_envelope(dt);


        if (core_mass > cnsts.parameters(Chandrasekar_mass) ||
	    core_mass <= cnsts.safety(minimum_mass_step)) {
            // || accreted_mass > 0.2) {
	    // (GN Mar 26 1999) test edge lid detanation 
	    // needs more research

	  if (is_binary_component()) 
	    get_binary()->dump("binev.data", false);
	  else
	    dump("binev.data", false);
	  
	  star_transformation_story(Disintegrated);
	  new disintegrated(*this);
	  // new neutron_star(*this);    // AIC
	  return;
        }

//    real t_nuc = nucleair_evolution_time();
        next_update_age = relative_age + cnsts.safety(maximum_timestep);

    if (relative_age <= 0.){
  //      if (t_nuc>=relative_age) {
           luminosity = 40;
        //relative_age = t_nuc;
        relative_age = 0.;
        
        }
        else {
	  // (GN May  4 1999) fit to Driebe et al 1999
	  real fit_mass = min(0.6, max(0.18, get_total_mass()));
	  real l_max = pow(10, (3.83 - 4.77* fit_mass));
	  luminosity = l_max/pow((relative_age), 1.4);
	}

	
//	real m_rel = min(0.999999, core_mass/cnsts.parameters(Chandrasekar_mass));

	   // Nauenberg, M, 1972, Apj 175, 417
	   // mu=4.039 is the mean molecular weight.
           real mu = 2.;
	   core_radius =
	     radius           = white_dwarf_radius(core_mass, relative_age);

//             radius = min(0.1, (0.0225/mu)
//	            * sqrt(1./pow(m_rel, cnsts.mathematics(two_third))
//                    - pow(m_rel, cnsts.mathematics(two_third))));

	// (GN+SPZ May  3 1999) critical mass for nova (Livio 1993; Saas-Fee)
	//   real m_crit = 1.7e-4*pow(get_total_mass(),-0.7)
	//               * pow(69.9*radius,2.8); // radius in 10^9 cm

	//if (envelope_mass >= m_crit) 
	//   thermo_nucleair_flash(dt);


//	   if (envelope_mass >= core_mass
//	                     * cnsts.parameters(thermo_nuclear_flash))
//	     thermo_nucleair_flash(dt);
       
	   update();
}

void white_dwarf::update() {

        detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
        effective_radius = max(effective_radius, radius);

}


star* white_dwarf::subtrac_mass_from_donor(const real dt, real& mdot) {
    
    mdot = mdot_limit(dt, mdot);
    
        if (mdot<=envelope_mass)
           envelope_mass -= mdot;
        else {
	  mdot -= envelope_mass;
	  envelope_mass = 0;
	  if (mdot >= core_mass)
	    mdot = core_mass - cnsts.safety(minimum_mass_step);
	  core_mass -= mdot;
	}

        return this;
}

real white_dwarf::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {

        if (mdot<0) {
           cerr << "white_dwarf::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   mdot = 0;
        }

    //For white dwarfs no difference currently between hydrogen/helium/.. accretion
	
// (GN+SPZ May  3 1999)
	real mu =1;
	if (is_binary_component() && 
	    !get_companion()->hydrogen_envelope_star()) mu =2;
	

	if (mdot >= eddington_limit(radius, dt, mu)) {
	  effective_radius = min(10., get_binary()->roche_radius(this));
	}

        mdot = accretion_limit(mdot, dt, hydrogen);

	// (Madelon+GN May 16 2011)
	// Implement accretion to core via retention efficiencies
	// Everything is dealt with in ::accretion_limit so
	// just add mdot to core

	  core_mass += mdot;
	  accreted_mass += mdot;

        adjust_accretor_age(mdot);
        //envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	set_spec_type(Accreting);

        return mdot;

     }


real  white_dwarf::accretion_limit(const real mdot, const real dt) {
    //needed for double_star::zeta and double_star::perform_mass_transfer
    return accretion_limit(mdot, dt, true);
}


real  white_dwarf::accretion_limit(const real mdot, const real dt, bool hydrogen) {

  if (dt < 0) return mdot;

  // (GN+Madelon May 12 2011) implement retention efficiencies
  //return min(maximum_steady_burning(dt), mdot);  

  real dmdt = 1.e-6*mdot/dt;
  real eta = retention_efficiency(dmdt, get_total_mass(), hydrogen);
  
  return mdot*eta;
}


real white_dwarf::retention_efficiency(real dmdt, real M_WD, bool hydrogen) {

  real eta = 1.;

  if (hydrogen) {

    real eta_H = retention_H(dmdt, M_WD);
    real eta_He = retention_He(eta_H*dmdt, M_WD);
    eta = eta_H*eta_He;
    
  } else {

    eta = retention_He(dmdt, M_WD);
  }

  return eta;
}


// (Madelon+GN May 16 2011)
// Retention efficiencies according to Nomoto
// (SilT Jun 19th 2012) rewritten equations
// (SilT Oct 16th 2013) no wind stripping -> c1 = 0 
real white_dwarf::retention_H(real dmdt, real M_WD) {

    if (dmdt < 1.e-7 || M_WD < 0.6) return 0.;
 
    real logdmdt = log10(dmdt); 
    real eta = 1.;
    
    real M_cr_log = log10(7.5e-7 * (M_WD - 0.4));
    real M_st = 3.1e-7 * (M_WD - 0.54);
    real M_st_log = log10(M_st);
    
    if (logdmdt > M_cr_log) 
        // (SilT Oct 16th 2013) no wind stripping -> c1 = 0 
        eta = pow(10, M_cr_log) / dmdt;   
//        eta = 4. * pow(10,M_cr_log) / (3.*pow(10,M_cr_log) + dmdt);  
    if (logdmdt < M_st_log)
        eta = (logdmdt - 1.e-7) / (M_st_log - 1.e-7);
      
    return eta;
}


real white_dwarf::retention_He(real dmdt, real M_WD) {
  
  real logdmdt = log10(dmdt);
  real eta = 0.;

  if (logdmdt > -7.8 && logdmdt < -5.9) 
    eta = -0.175 * pow(logdmdt+5.35,2) + 1.05;
  if (logdmdt > -5.9 && logdmdt < -5.0) eta =  1.;

  return eta;
}



void white_dwarf::adjust_accretor_age(const real mdot,
				      const bool rejuvenate) {
//        real m_rel_new;
//        real m_tot_new = get_total_mass() + mdot;
//        if (m_tot_new>relative_mass)
//           m_rel_new = m_tot_new;
//        else m_rel_new = relative_mass;

    real t_nuc_old = 0.;//nucleair_evolution_time();
//	real z_new = get_metalicity();
    real t_nuc_new = 0.;//nucleair_evolution_time(m_rel_new, m_tot_new, z_new);

        real dtime = relative_age - t_nuc_old;

        relative_age = t_nuc_new + dtime
	             * (1-pow(mdot/(get_total_mass()+mdot),
		     	      cnsts.parameters(rejuvenation_exponent)));
	
     }

real white_dwarf::zeta_thermal() {

     // Based on white dwarf mass radius relation.
     // zeta = (m/r)dr/dm;
  
     return -cnsts.mathematics(one_third);
}

real white_dwarf::zeta_adiabatic() {

     // Based on white dwarf mass radius relation.
     // zeta = (m/r)dr/dm;

     return -cnsts.mathematics(one_third);
}

star* white_dwarf::merge_elements(star* str) {

     envelope_mass = 0;		//Destroy disc.
	
     real merger_core = str->get_core_mass();

     add_mass_to_accretor(str->get_envelope_mass(), str->hydrogen_envelope_star());

     if (relative_mass<get_total_mass() + merger_core)
       relative_mass=get_total_mass() + merger_core;
     core_mass += merger_core;

     spec_type[Merger]=Merger;
     instantaneous_element();

     return this;
}

star* white_dwarf::reduce_mass(const real mdot) {

      if (envelope_mass<+mdot) {
	real mdot_rest = mdot - envelope_mass;
	envelope_mass = 0;
	if (mdot_rest >= core_mass)
	  mdot_rest = core_mass - cnsts.safety(minimum_mass_step);
	core_mass -= mdot_rest;
      }
      else envelope_mass -= mdot;

      return this;
}


real white_dwarf::gyration_radius_sq() {

//  return cnsts.parameters(homogeneous_sphere_gyration_radius_sq); 
  
// (SilT & AD 13 Feb 22) 
// priv. comm. Antonio Claret, Gabriele Columba, Camilla Danielski (see Claret 2019, 628, 29


    const int size_mbins = 19;
	const int size_tbins = 7;
//	real mass_array [size_mbins] = {0.6, 0.8, 0.85, 0.92, 1., 1.2, 1.4, 1.6,1.8, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10.}; 
    real mass_array [size_mbins] = {0.379888,0.436064,0.449099,0.469484,0.515458,0.534108,0.543948,0.551604,0.558214,0.562217, 0.579593, 0.606072, 0.668300,0.792400,	0.803200,0.854300,0.914900,	0.979500,1.104000};


    int im_l, im_u; 
    if (get_total_mass() <= mass_array[0]) im_l = im_u = 0;
    else if (get_total_mass() >= mass_array[size_mbins-1]) im_l = im_u = size_mbins-1;
    else {
        for (im_u = 0; im_u < size_mbins; im_u++){
                if (get_total_mass() < mass_array[im_u]) break; }
        im_l = im_u-1;
    }
    real m_u = mass_array[im_u];
    real m_l = mass_array[im_l];
//    PRC(im_u);PRL(im_l);
//    PRC(m_u);PRC(m_l);PRL(relative_mass);



    real timestamp_array[size_mbins][size_tbins] =  {{ 
        81955700000.0 , 81956198200.0 , 81956255100.0 , 81956443400.0 , 81992213800.0 , 82975363000.0 , 97780762300.0 }, {
        29018600000.0 , 29018870600.0 , 29018897700.0 , 29018952200.0 , 29050423600.0 , 29835839900.0 , 45607247900.0 }, {
        22928900000.0 , 23022439127.7 , 23031749056.6 , 23047400000.0 , 23080988100.0 , 23856471700.0 , 39915020300.0 }, {
        17023800000.0 , 17028147900.0 , 17029451300.0 , 17029582600.0 , 17049157400.0 , 17405411400.0 , 26097458900.0 }, {
        12462200000.0 , 12462606400.0 , 12462607200.0 , 12462722300.0 , 12477096900.0 , 12787815300.0 , 21345748400.0 }, {
        6513840000.0 , 6513884420.0 , 6513888060.0 , 6513965980.0 , 6526380330.0 , 6816330780.0 , 15391243870.0 }, {
        3867600000.0 , 3867641800.0 , 3867644520.0 , 3867678400.0 , 3878067970.0 , 4145011140.0 , 12808680200.0 }, {
        2536800000.0 , 2536803200.0 , 2536803480.0 , 2536812860.0 , 2546637230.0 , 2811251250.0 , 11533649040.0 }, {
        1798550000.0 , 1798552320.0 , 1798552520.0 , 1798558750.0 , 1807663440.0 , 2071939610.0 , 10814551850.0 }, {
        1495930000.0 , 1495931700.0 , 1495931850.0 , 1495936160.0 , 1505567300.0 , 1771428230.0 , 10523809000.0 }, {
        801259000.0 , 801259993.0 , 801260081.0 , 801262319.0 , 801262786.0 , 801263718.0 , 11612959288.8 }, {
        477720000.0 , 477720605.0 , 477720658.0 , 477721844.0 , 477722072.0 , 477722278.0 , 12910031530.2 }, {
        215327000.0 , 215327100.0 , 215327110.0 , 215337480.0 , 222789470.0 , 535945600.0 , 15499169430.0 }, {
        121464000.0 , 121464060.0 , 121464070.0 , 121523460.0 , 123873980.0 , 216512970.0 , 7255527420.0 }, {
        78723000.0 , 78723069.0 , 78723076.0 , 78828520.0 , 80732909.0 , 150232656.0 , 6371399256.0 }, {
        55922400.0 , 55922571.0 , 55922590.0 , 55992842.0 , 59215221.0 , 200431112.0 , 9377535672.0 }, {
        41928600.0 , 41928905.0 , 41928945.0 , 42015547.0 , 44605609.0 , 157267453.0 , 8909536033.0 }, {
        33221800.0 , 33222117.0 , 33222160.0 , 33319965.0 , 34594968.0 , 96889526.0 , 6445349229.0 }, {
        27271500.0 , 27272600.0 , 27378010.0 , 27399966.0 , 28165618.0 , 74578252.0 , 7962969201.0 }};
        
        
        




    real beta_array[size_mbins][size_tbins] =  {{
        0.0875335 , 0.0330072 , 0.00300989 , 0.0937004 , 0.318307 , 0.379106 , 0.433122 }, {
        0.122052 , 0.0461817 , 0.00199293 , 0.0964927 , 0.336591 , 0.390692 , 0.432942 }, {
        0.132354 , 0.0504088 , 0.00183354 , 0.0991473 , 0.34202 , 0.393371 , 0.432807 }, {
        0.039979 , 0.0396956 , 0.0791797 , 0.118232 , 0.352949 , 0.393973 , 0.429487 }, {
        0.0863324 , 0.0851192 , 0.00203793 , 0.11867 , 0.362515 , 0.400497 , 0.429073 }, {
        0.153656 , 0.0740045 , 0.00205782 , 0.121062 , 0.365874 , 0.402575 , 0.428952 }, {
        0.201359 , 0.085153 , 0.00182936 , 0.0996407 , 0.364753 , 0.402207 , 0.428547 }, {
        0.189108 , 0.0704742 , 0.00124846 , 0.0939188 , 0.365237 , 0.402413 , 0.428113 }, {
        0.203384 , 0.0748965 , 0.00108296 , 0.0957417 , 0.365746 , 0.403024 , 0.428018 }, {
        0.215236 , 0.0817385 , 0.00097948 , 0.101165 , 0.368408 , 0.403696 , 0.428402 }, {
        0.234866 , 0.0891593 , 0.000825986 , 0.106575 , 0.170514 , 0.4278515 , 0.4278515 }, {
        0.252813 , 0.0967703 , 0.000701838 , 0.112472 , 0.171158 , 0.427301 , 0.427301 }, {
        0.2434 , 0.1088 , 0.001026 , 0.1272 , 0.3817 , 0.4092 , 0.4262 }, {
        0.3006 , 0.1699 , 0.001846 , 0.1388 , 0.3972 , 0.414 , 0.4223 }, {
        0.2401 , 0.09153 , 0.005978 , 0.1432 , 0.3949 , 0.4122 , 0.4215 }, {
        0.2333 , 0.09277 , 0.01119 , 0.2047 , 0.4 , 0.4119 , 0.4186 }, {
        0.1875 , 0.07913 , 0.008603 , 0.2002 , 0.3977 , 0.4074 , 0.4138 }, {
        0.1729 , 0.06976 , 0.005729 , 0.2072 , 0.3923 , 0.4013 , 0.4077 }, {
        0.1658 , 0.0141 , 0.0394 , 0.2094 , 0.3839 , 0.3896 , 0.3942 }};


    

    //find time for m_u
    real time_yrs_u = timestamp_array[im_u][0] + relative_age * 1e6; // in yrs
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
    real time_yrs_l = timestamp_array[im_l][0] + relative_age * 1e6; // in yrs
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
        beta = beta_m_l_at_t + (beta_m_u_at_t - beta_m_l_at_t) * (get_total_mass() - m_l) / (m_u - m_l);
        
    
    PRC( beta_array[im_l][it_ll]);PRC(beta_array[im_l][it_lu]);       
    PRC( beta_array[im_u][it_ul]);PRL(beta_array[im_u][it_uu]);       
//    PRC(beta_m_l_at_t);PRC(beta_m_u_at_t);PRL(beta); 
//    PRL(cnsts.parameters(homogeneous_sphere_gyration_radius_sq));
    return beta*beta;           

//  return cnsts.parameters(homogeneous_sphere_gyration_radius_sq); 
  
  
}




//absidal motion constant
real white_dwarf::amc() {

// (SilT & AD 13 Feb 22) 
// priv. comm. Antonio Claret, Gabriele Columba, Camilla Danielski (see Claret 2019, 628, 29

    const int size_mbins = 19;
	const int size_tbins = 7;
//	real mass_array [size_mbins] = {0.6, 0.8, 0.85, 0.92, 1., 1.2, 1.4, 1.6,1.8, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., 10.}; 
    real mass_array [size_mbins] = {0.379888,0.436064,0.449099,0.469484,0.515458,0.534108,0.543948,0.551604,0.558214,0.562217, 0.579593, 0.606072, 0.668300,0.792400,	0.803200,0.854300,0.914900,	0.979500,1.104000};


    int im_l, im_u; 
    if (get_total_mass() <= mass_array[0]) im_l = im_u = 0;
    else if (get_total_mass() >= mass_array[size_mbins-1]) im_l = im_u = size_mbins-1;
    else {
        for (im_u = 0; im_u < size_mbins; im_u++){
                if (get_total_mass() < mass_array[im_u]) break; }
        im_l = im_u-1;
    }
    real m_u = mass_array[im_u];
    real m_l = mass_array[im_l];
//    PRC(m_u);PRC(m_l);PRL(relative_mass);



    real timestamp_array[size_mbins][size_tbins] =  {{ 
        81955700000.0 , 81956198200.0 , 81956255100.0 , 81956443400.0 , 81992213800.0 , 82975363000.0 , 97780762300.0 }, {
        29018600000.0 , 29018870600.0 , 29018897700.0 , 29018952200.0 , 29050423600.0 , 29835839900.0 , 45607247900.0 }, {
        22928900000.0 , 23022439127.7 , 23031749056.6 , 23047400000.0 , 23080988100.0 , 23856471700.0 , 39915020300.0 }, {
        17023800000.0 , 17028147900.0 , 17029451300.0 , 17029582600.0 , 17049157400.0 , 17405411400.0 , 26097458900.0 }, {
        12462200000.0 , 12462606400.0 , 12462607200.0 , 12462722300.0 , 12477096900.0 , 12787815300.0 , 21345748400.0 }, {
        6513840000.0 , 6513884420.0 , 6513888060.0 , 6513965980.0 , 6526380330.0 , 6816330780.0 , 15391243870.0 }, {
        3867600000.0 , 3867641800.0 , 3867644520.0 , 3867678400.0 , 3878067970.0 , 4145011140.0 , 12808680200.0 }, {
        2536800000.0 , 2536803200.0 , 2536803480.0 , 2536812860.0 , 2546637230.0 , 2811251250.0 , 11533649040.0 }, {
        1798550000.0 , 1798552320.0 , 1798552520.0 , 1798558750.0 , 1807663440.0 , 2071939610.0 , 10814551850.0 }, {
        1495930000.0 , 1495931700.0 , 1495931850.0 , 1495936160.0 , 1505567300.0 , 1771428230.0 , 10523809000.0 }, {
        801259000.0 , 801259993.0 , 801260081.0 , 801262319.0 , 801262786.0 , 801263718.0 , 11612959288.8 }, {
        477720000.0 , 477720605.0 , 477720658.0 , 477721844.0 , 477722072.0 , 477722278.0 , 12910031530.2 }, {
        215327000.0 , 215327100.0 , 215327110.0 , 215337480.0 , 222789470.0 , 535945600.0 , 15499169430.0 }, {
        121464000.0 , 121464060.0 , 121464070.0 , 121523460.0 , 123873980.0 , 216512970.0 , 7255527420.0 }, {
        78723000.0 , 78723069.0 , 78723076.0 , 78828520.0 , 80732909.0 , 150232656.0 , 6371399256.0 }, {
        55922400.0 , 55922571.0 , 55922590.0 , 55992842.0 , 59215221.0 , 200431112.0 , 9377535672.0 }, {
        41928600.0 , 41928905.0 , 41928945.0 , 42015547.0 , 44605609.0 , 157267453.0 , 8909536033.0 }, {
        33221800.0 , 33222117.0 , 33222160.0 , 33319965.0 , 34594968.0 , 96889526.0 , 6445349229.0 }, {
        27271500.0 , 27272600.0 , 27378010.0 , 27399966.0 , 28165618.0 , 74578252.0 , 7962969201.0 }};
 


    real log_amc_array[size_mbins][size_tbins] =  {{
        -2.53196 , -3.19008 , -5.98352 , -4.10695 , -1.58973 , -1.2167 , -0.929204 }, {
        -2.14636 , -2.82775 , -6.20791 , -4.07928 , -1.46954 , -1.1511 , -0.929444 }, {
        -2.06209 , -2.74406 , -6.24221 , -4.03151 , -1.43524 , -1.13619 , -0.929984 }, {
        -4.61231 , -4.73476 , -4.26493 , -3.60686 , -1.36834 , -1.13314 , -0.946891 }, {
        -2.48002 , -3.65528 , -6.3724 , -3.60831 , -1.31063 , -1.0971 , -0.94829 }, {
        -1.91073 , -2.45601 , -6.19385 , -3.51276 , -1.29041 , -1.08545 , -0.948388 }, {
        -1.63041 , -2.33308 , -6.27036 , -3.95881 , -1.29676 , -1.08723 , -0.950202 }, {
        -1.69074 , -2.42242 , -6.48243 , -4.1024 , -1.29349 , -1.08566 , -0.951926 }, {
        -1.61766 , -2.36413 , -6.65524 , -4.07945 , -1.29048 , -1.08237 , -0.952366 }, {
        -1.55918 , -2.28595 , -6.73109 , -3.98212 , -1.27498 , -1.07873 , -0.950367 }, {
        -1.46731 , -2.20214 , -6.8227 , -3.89022 , -2.8894 , -0.95265025 , -0.95265025 }, {
        -1.38694 , -2.12327 , -6.9179 , -3.78637 , -2.88675 , -0.9549335 , -0.9549335 }, {
        -1.43 , -2.032 , -6.653 , -3.476 , -1.197 , -1.048 , -0.9595 }, {
        -1.183 , -1.632 , -5.975 , -3.233 , -1.109 , -1.021 , -0.9777 }, {
        -1.406 , -2.161 , -4.962 , -3.177 , -1.12 , -1.028 , -0.9805 }, {
        -1.456 , -2.168 , -4.477 , -2.446 , -1.091 , -1.029 , -0.9935 }, {
        -1.669 , -2.317 , -4.889 , -2.487 , -1.101 , -1.05 , -1.016 }, {
        -1.739 , -2.421 , -5.352 , -2.418 , -1.128 , -1.08 , -1.045 }, {
        -1.795 , -5.657 , -4.448 , -2.389 , -1.168 , -1.137 , -1.111 }};

    

    //find time for m_u
    real time_yrs_u = timestamp_array[im_u][0] + relative_age * 1e6; // in yrs
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
    real time_yrs_l = timestamp_array[im_l][0] + relative_age * 1e6; // in yrs
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
    real log_amc_m_l_at_t, log_amc_m_u_at_t, log_amc;

    if (it_uu == it_ul){
        log_amc_m_u_at_t = log_amc_array[im_u][it_uu];
    }
    else {
        log_amc_m_u_at_t = log_amc_array[im_u][it_ul] + (log_amc_array[im_u][it_uu] - log_amc_array[im_u][it_ul]) * (time_yrs_u - t_ul) / (t_uu - t_ul);
    }
    if (it_lu == it_ll){
        log_amc_m_l_at_t = log_amc_array[im_l][it_lu];
    }
    else {
        log_amc_m_l_at_t = log_amc_array[im_l][it_ll] + (log_amc_array[im_l][it_lu] - log_amc_array[im_l][it_ll]) * (time_yrs_l - t_ll) / (t_lu - t_ll);
    }


    if (im_u==im_l)
        log_amc = log_amc_m_u_at_t;
    else    
        log_amc = log_amc_m_l_at_t + (log_amc_m_u_at_t - log_amc_m_l_at_t) * (get_total_mass() - m_l) / (m_u - m_l);
        
        
//    PRC( log_amc_array[im_l][it_ll]);PRC(log_amc_array[im_l][it_lu]);       
//    PRC( log_amc_array[im_u][it_ul]);PRL(log_amc_array[im_u][it_uu]);       
//    PRC(log_amc_m_l_at_t);PRC(log_amc_m_u_at_t);PRL(log_amc); 
    return pow(10,log_amc);           
  
}
  


//stellar_type white_dwarf::get_element_type(){
//
//  if (core_mass     < cnsts.parameters(helium_dwarf_mass_limit) &&
//      relative_mass < cnsts.parameters(upper_ZAMS_mass_for_degenerate_core))
//    return Helium_Dwarf;
//  
//  else if (relative_mass >= cnsts.parameters(super_giant2neutron_star) &&
//           core_mass > cnsts.parameters(carbon_dwarf_mass_limit))
//    return Oxygen_Dwarf;
//  
//  else 
//    return Carbon_Dwarf;
//}


real white_dwarf::get_evolve_timestep() {
        
    return max(next_update_age - relative_age, 0.0001);
}
