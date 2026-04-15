
#include "bagel.h"

bagel_type::bagel_type(const prior_function_type& prior_function,
		       const feature_vector_function_type& feature_vector_function,
		       const probability_type& p0,
		       const probability_type& p,
		       const real_type& s,
		       const real_type& lst_nu,
		       const real_type& lst_mu,
		       const int& n)
{
  t = 1;
  max_num_particles = n;
  initial_particle = particle_type(prior_function,feature_vector_function,0,p0,p,s,lst_nu,lst_mu);
}


real_type bagel_type::weight_0_t() const
{
  return particles.front().weight;
}

std::list<std::list<ratio_type> > bagel_type::ratios() const
{
  std::list<std::list<ratio_type> >  particle_ratios;
  for(auto& p : particles)
    {
      particle_ratios.push_back(p.ratios);
    }
  return particle_ratios;

}


real_type weight_0_t(const bagel_type& bagel)
{
  return bagel.weight_0_t();
}


std::list<std::list<ratio_type> > ratios(const bagel_type& bagel)
{
  return bagel.ratios();
}



bagel_type& update(bagel_type& bagel, const real_type& y)
{
  return bagel.update(y);
}


bagel_type& bagel_type::update(const real_type& y)
{

  
  if(t == 1)
    {
      initial_particle = theorem_2(initial_particle,t);
      initial_particle = theorem_1(initial_particle,initial_particle,t);
      initial_particle = theorem_4(initial_particle,t,y);
      initial_particle = theorem_3(initial_particle,t,y);
      initial_particle.weight = initial_particle.p0; 
      particles.push_back(initial_particle);
    }
  else
    {
      // new particle
      initial_particle = particles.front();
      particles.push_back(theorem_2(initial_particle,t));
      real_type sum_of_weights = 0.0;      
      for(auto& p : particles)
	{
	  p = theorem_1(initial_particle,p,t);
	}
      
      for(auto& p : particles)
	{
	  p = theorem_4(p,t,y);
	  p = theorem_3(p,t,y);
	  sum_of_weights += p.weight;
	}
      for(auto& p : particles)
	{
	  p.weight = p.weight/sum_of_weights;
	}

      // prune
      if(particles.size() > max_num_particles && particles.size() > 3)
	{
	  auto it_1 = particles.begin();
	  it_1++;
	  auto it_n_minus_1 = particles.end();
	  it_n_minus_1--;
	  it_n_minus_1--; 
	  auto it_min = std::min_element(it_1,
					 it_n_minus_1,
					 [](auto& x,auto& y){return x.weight < y.weight;});		 
	  auto it_right_of_min = it_min;
	  it_right_of_min++;
	  // update ratios
	  auto combined_weight = it_right_of_min -> weight + it_min -> weight;
	  auto weight = it_min -> weight;
	  std::transform(it_min->ratios.begin(),
			 it_min->ratios.end(),
			 it_min->ratios.begin(),
			 [&combined_weight,&weight](auto& ratio){return ratio*weight/combined_weight;});
	  weight = it_right_of_min -> weight;
	  std::transform(it_right_of_min->ratios.begin(),
			 it_right_of_min->ratios.end(),
			 it_right_of_min->ratios.begin(),
			 [&combined_weight,&weight](auto& ratio){return ratio*weight/combined_weight;});
	  it_right_of_min->ratios.insert(it_right_of_min->ratios.begin(),it_min->ratios.begin(),it_min->ratios.end());
	  
	  // update weights
	  it_right_of_min -> weight += it_min -> weight;

	  // evict the pruned particle
	  particles.erase(it_min);
	}
    }
  
  t = t + 1;


  
  return *this;
}

