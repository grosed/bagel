
#ifndef ___PRUNE_H___
#define ___PRUNE_H___


#include "particle_type.h"
#include "plurality_type.h"
#include <list>


template<typename noise, plurality_type plurality>
std::list<particle_type<noise,plurality> >& prune(std::list<particle_type<noise,plurality> >& particles,const int& max_num_particles)
{

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
  
  return particles;
  
}








#endif
