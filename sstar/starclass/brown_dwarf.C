//
// brown_dwarf.C
//
//to be based on Dantona, F., Mazzitelli, I., 1985, ApJ 296, 502
// and  2000astro.ph..5557, Chabrier, G.; Baraffe, I.; Allard, F.; 
// Hauschildt, P.

#include "brown_dwarf.h"
#include "main_sequence.h"
#include "proto_star.h"

brown_dwarf::brown_dwarf(proto_star & p) : single_star(p) {
    
      delete &p;

      real m_tot = get_total_mass();
      core_mass = brown_dwarf_core_mass();
      envelope_mass = m_tot - core_mass;

      last_update_age = 0;
      relative_age = 0;
	  
      instantaneous_element();
      update();

      post_constructor();
}

brown_dwarf::brown_dwarf(main_sequence & m) : single_star(m) {
    
      delete &m;

      real m_tot = get_total_mass();
      core_mass = brown_dwarf_core_mass();
      envelope_mass = m_tot - core_mass;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;
      instantaneous_element();
      update();

      post_constructor();
}

void brown_dwarf::instantaneous_element() {

     luminosity = 1.e-4; 

     core_radius = radius = brown_dwarf_radius();
       
     update();
}

real brown_dwarf::get_evolve_timestep() {
    // (GN+SPZ Apr 28 1999) was a bit too small
    //  return max(next_update_age - relative_age
    //	     -0.5*cnsts.safety(minimum_timestep),
    //	     cnsts.safety(minimum_timestep));
    
    // (GN+SPZ May  5 1999) type change time must be small because of rapid
    // growth of giants at end phase 0.0001 seems to be OK (?)
    // return max(next_update_age - relative_age - (0.5*0.001), 0.001);

    return max(next_update_age - relative_age, 0.0001);
}



void brown_dwarf::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;

        next_update_age = relative_age + cnsts.safety(maximum_timestep);

	//Burrows & Libert 1993, J. Rev. Mod. Phys. 65, 301
	luminosity = 938 * pow(relative_mass, 2.64); 
	if(relative_age>1)
	  luminosity = 938 * pow(relative_mass, 2.64) / pow(relative_age, 1.3);

        core_radius = radius = brown_dwarf_radius();
       
        update();
     }

void brown_dwarf::update() {

     detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
     effective_radius = radius;

     }



// (SilT May 26 2021) 
// fit to Chen & Kipping (2017)
real brown_dwarf::brown_dwarf_radius() {

    real m_earth = cnsts.parameters(Mearth);
    real m_jupiter = cnsts.parameters(Mjupiter);
    real r_jupiter = cnsts.parameters(Rjupiter);           
       
    if (get_total_mass() > 0.414 * m_jupiter) // Jovian planets & brown dwarfs
        return (0.027095 + 1.205219* pow(get_total_mass()/m_jupiter,-0.044))*r_jupiter;
    else if (get_total_mass() >2.04 * m_earth) // Neptunian planets
        return  2.151774*r_jupiter * pow(get_total_mass()/m_jupiter, 0.589); 
    else // Terran planets
        return 0.31832 * pow(get_total_mass(), 0.28); 
    
     }


real brown_dwarf::brown_dwarf_core_mass() {
    
        return 0.01 * get_total_mass();
     }

star* brown_dwarf::subtrac_mass_from_donor(const real dt, real& mdot) {

    mdot = mdot_limit(dt, mdot);
    
        if (mdot<=envelope_mass)
	  envelope_mass -= mdot;
        else if (mdot>envelope_mass) 
	  envelope_mass = 0;

        return this;
     }


real brown_dwarf::add_mass_to_accretor(real mdot, bool, const real dt) {

        if (mdot<0) {
           cerr << "brown_dwarf::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   mdot = 0;
        }

        mdot = accretion_limit(mdot, dt);
 
        envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	set_spec_type(Accreting);
	
        return mdot;
     }

real brown_dwarf::accretion_limit(const real mdot, const real dt) {

  if (dt < 0) return mdot;

        real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

        if(mdot>=eddington)
	  return eddington;

        return mdot;
     }


real brown_dwarf::zeta_thermal() {

     return 0;
}

star* brown_dwarf::merge_elements(star* str) {

     real merger_core = str->get_core_mass();

     add_mass_to_accretor(str->get_envelope_mass(), 
			  cnsts.parameters(spiral_in_time), str->hydrogen_envelope_star());

     if (relative_mass < get_total_mass() + merger_core)
       relative_mass=get_total_mass() + merger_core;
     core_mass += merger_core;

     spec_type[Merger]=Merger;
     instantaneous_element();

     return this;
}

star* brown_dwarf::reduce_mass(const real mdot) {

      if (envelope_mass < mdot)
	envelope_mass = 0;
      else
	envelope_mass -= mdot;

      return this;
}




// (SilT 13 Feb 22) 
// based on 
real brown_dwarf::gyration_radius_sq() {

//    return cnsts.parameters(convective_star_gyration_radius_sq); 

    real time_yrs = pow(10,6) * relative_age;
    const int size_mbins = 14;
	const int size_tbins = 15;
	real mass_array [size_mbins] = {5.00e-04, 1.00e-03, 3.00e-03, 5.00e-03, 0.01,
						0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
						
    real timestamp_array [size_tbins] = {0.000000000000000000e+00, 2.232041677489456255e+06, 6.622961400217147544e+06, 1.526087179673122987e+07, 3.225355231238762662e+07, 6.568191378551624715e+07, 1.314428953886585236e+08, 2.608092929874467552e+08, 5.153016138090125918e+08, 1.015944304332175136e+09, 2.000819222228052855e+09, 3.938286046894550323e+09, 7.749711872118158340e+09, 1.524762958417077827e+10, 2.999769231520001221e+10};

    real beta_array[size_mbins][size_tbins] = {{4.284212710334296403e-01,4.425154488598525004e-01,4.576844442792262213e-01,4.686924379730283974e-01,4.867883325341965106e-01,4.993132026773534071e-01,5.070888182759955010e-01,5.116918614668987120e-01,5.162373251319171130e-01,5.207169049823777707e-01,5.252920168281122182e-01,5.291104710316946180e-01,5.326218968569180756e-01,5.358000000000000540e-01,5.381000000000000227e-01},{4.557940766512104114e-01,4.656934435769591873e-01,4.742127602925355068e-01,4.827284167914714930e-01,4.903676380943810531e-01,4.972195280552020313e-01,5.043582815216275383e-01,5.095282184446862095e-01,5.130889684688720864e-01,5.160970550353507491e-01,5.190182774140561506e-01,5.218466055316945962e-01,5.242999999999999883e-01,5.260000000000000231e-01,5.273999999999999799e-01},{4.845896012524040786e-01,4.901506183820633633e-01,4.966729553649681939e-01,5.010661332383635758e-01,5.037453903134077038e-01,5.055754151186144796e-01,5.073511861043721316e-01,5.091863478675210386e-01,5.103783988158270279e-01,5.115267756184218584e-01,5.126116505070280605e-01,5.134491173316945156e-01,5.141999999999999904e-01,5.150000000000000133e-01,5.159000000000000252e-01},{4.832197862782841025e-01,4.952204341790670261e-01,4.988809887724526271e-01,5.020188588278190478e-01,5.036835193619633166e-01,5.043025129259661732e-01,5.044487552755508863e-01,5.048189693533623457e-01,5.054216308772180311e-01,5.059621709929150724e-01,5.064999999999999503e-01,5.071999999999999842e-01,5.077000000000000401e-01,5.081999999999999851e-01,5.087000000000000410e-01},{4.805186877051476069e-01,4.845582523632057104e-01,4.925049696253366305e-01,4.940246166652708659e-01,4.951450236704253438e-01,4.952930430177609478e-01,4.951204208733632206e-01,4.949784130022125694e-01,4.949760793186089858e-01,4.950999999999999845e-01,4.954000000000000070e-01,4.958000000000000185e-01,4.964000000000000079e-01,4.966999999999999749e-01,4.969000000000000083e-01},{4.645784519465647189e-01,4.687200697163165231e-01,4.754641099230546564e-01,4.787594283405924300e-01,4.833057222882783943e-01,4.838025372181446460e-01,4.839420771201419780e-01,4.838999999999999968e-01,4.839457393886090197e-01,4.839999999999999858e-01,4.842000000000000193e-01,4.843000000000000083e-01,4.845999999999999752e-01,4.843999999999999972e-01,4.838999999999999968e-01},{4.604976552976300375e-01,4.628403761394166982e-01,4.659701866056075725e-01,4.763519833751501120e-01,4.778492431244955041e-01,4.781431887917103607e-01,4.782000000000000139e-01,4.781000000000000250e-01,4.779999999999999805e-01,4.779999999999999805e-01,4.781000000000000250e-01,4.781000000000000250e-01,4.774999999999999800e-01,4.763000000000000012e-01,4.758999999999999897e-01},{4.574172288571127787e-01,4.595454933606252235e-01,4.656578162936672038e-01,4.707894195867353693e-01,4.738999999999999879e-01,4.740969388717103383e-01,4.741000000000000214e-01,4.741000000000000214e-01,4.739999999999999769e-01,4.738999999999999879e-01,4.741000000000000214e-01,4.744999999999999774e-01,4.738999999999999879e-01,4.727999999999999980e-01,4.726000000000000201e-01},{4.529535870413303744e-01,4.566992517750683156e-01,4.642908690318094389e-01,4.664104193116077601e-01,4.710159016224775286e-01,4.712999999999999967e-01,4.712999999999999967e-01,4.712999999999999967e-01,4.712000000000000077e-01,4.717999999999999972e-01,4.728999999999999870e-01,4.732000000000000095e-01,4.717999999999999972e-01,4.702999999999999958e-01,4.701000000000000179e-01},{4.492415680123580946e-01,4.551376674164676062e-01,4.626798727808058564e-01,4.649018270269103525e-01,4.669640310579820097e-01,4.692109288882896911e-01,4.692000000000000060e-01,4.692000000000000060e-01,4.692000000000000060e-01,4.706000000000000183e-01,4.718999999999999861e-01,4.722999999999999976e-01,4.699999999999999734e-01,4.677999999999999936e-01,4.672000000000000042e-01},{4.466379197884707275e-01,4.540503114411651442e-01,4.612549871748234942e-01,4.636473242077385204e-01,4.644516846622477324e-01,4.673607556751309677e-01,4.676000000000000156e-01,4.675000000000000266e-01,4.677000000000000046e-01,4.693254997826434849e-01,4.703999999999999848e-01,4.709999999999999742e-01,4.712999999999999967e-01,4.669407946205844651e-01,4.655000000000000249e-01},{4.451424819141328038e-01,4.531260432917555980e-01,4.595638841633183369e-01,4.623510687035423672e-01,4.633800125644954582e-01,4.648757205351309874e-01,4.662999999999999923e-01,4.662000000000000033e-01,4.663999999999999813e-01,4.676000000000000156e-01,4.681000000000000161e-01,4.682000000000000051e-01,4.682000000000000051e-01,4.682000000000000051e-01,4.682000000000000051e-01},{4.438899601109268511e-01,4.522604197871521459e-01,4.582248436600245833e-01,4.612616713952248926e-01,4.624263792422477204e-01,4.629949464717103558e-01,4.652000000000000024e-01,4.649025058442510239e-01,4.648999999999999799e-01,4.652999999999999914e-01,4.653999999999999804e-01,4.653999999999999804e-01,4.653999999999999804e-01,4.653999999999999804e-01,4.653999999999999804e-01},{4.428240029278892598e-01,4.513416905078856400e-01,4.570269453979082797e-01,4.601854787020318582e-01,4.616303921622477713e-01,4.621000000000000107e-01,4.636000000000000121e-01,4.638984930442510701e-01,4.636000000000000121e-01,4.636000000000000121e-01,4.636000000000000121e-01,4.636000000000000121e-01,4.636000000000000121e-01,4.636000000000000121e-01,4.636000000000000121e-01} };



    
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


    //find time
    int it_u, it_l; 
    if (time_yrs <=  timestamp_array[0]) it_u = it_l = 0 ;
    else if (time_yrs >=  timestamp_array[size_tbins-1]) it_u = it_l = size_tbins-1 ;
    else{
        for (it_u = 0; it_u < size_tbins; it_u++){
                if (time_yrs < timestamp_array[it_u]) break; }
        it_l = it_u-1;        
    } 
    real t_u = timestamp_array[it_u];
    real t_l = timestamp_array[it_l];
//    PRC(t_u);PRC(t_l);PRC(relative_age);PRL(time_yrs);
            

    // for readibility keep this separate from previous two loops            
    //step three, interpolation in mass and time...
    real beta_m_l_at_t, beta_m_u_at_t, beta;
    
    if (it_u == it_l){
        beta_m_l_at_t = beta_array[im_l][it_l];
        beta_m_u_at_t = beta_array[im_u][it_u];
    }
    else {
        beta_m_l_at_t = beta_array[im_l][it_l] + (beta_array[im_l][it_u] - beta_array[im_l][it_l]) * (time_yrs - t_l) / (t_u - t_l);
        beta_m_u_at_t = beta_array[im_u][it_l] + (beta_array[im_u][it_u] - beta_array[im_u][it_l]) * (time_yrs - t_l) / (t_u - t_l);
    }
    
    if (im_u==im_l)
        beta = beta_m_u_at_t;
    else    
        beta = beta_m_l_at_t + (beta_m_u_at_t - beta_m_l_at_t) * (relative_mass - m_l) / (m_u - m_l);
        
        
//    PRC( beta_array[im_l][it_l]);PRC(beta_array[im_l][it_u]);       
//    PRC( beta_array[im_u][it_l]);PRL(beta_array[im_u][it_u]);       
//    PRC(beta_m_l_at_t);PRC(beta_m_u_at_t);PRL(beta); 
//    PRL(cnsts.parameters(convective_star_gyration_radius_sq));
    return beta*beta;           



}


stellar_type brown_dwarf::get_element_type() {

    // (SilT January 18 2021) 
    if (get_total_mass() < cnsts.parameters(brown_dwarf_mass_limit))
//  if (get_total_mass() < 0.1*cnsts.parameters(minimum_main_sequence))
      return Planet;
    else
      return Brown_Dwarf;
}
	
