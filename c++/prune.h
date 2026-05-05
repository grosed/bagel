
#ifndef ___PRUNE_H___
#define ___PRUNE_H___


#include "particle_type.h"
#include "plurality_type.h"
#include <list>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>


#include <iostream>


template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<unknown_variance,plurality_type::univariate> >; }
{
  matrix mu_i = a.post.mu;
  matrix mu_i_plus_1 = b.post.mu;
  matrix sigma_i = a.post.sigma;
  matrix sigma_i_plus_1 = b.post.sigma;
  
  double nu_i = a.noise_structure.nu;
  double nu_i_plus_1 = b.noise_structure.nu;
  double iota_i = a.noise_structure.iota;
  double iota_i_plus_1 = b.noise_structure.iota; 
  matrix I = matrix::Identity(sigma_i.rows(),sigma_i.cols());
  double result = (sigma_i_plus_1.inverse() * sigma_i - I).trace();

  result = result + (nu_i/iota_i)*((mu_i_plus_1 - mu_i).transpose() * sigma_i_plus_1.inverse() * (mu_i_plus_1 - mu_i))(0,0);
  result = result + std::log((sigma_i_plus_1.inverse() * sigma_i).determinant());
  result = result + nu_i_plus_1 * std::log(iota_i/iota_i_plus_1);
  result = result - std::log(boost::math::tgamma(nu_i)/boost::math::tgamma(nu_i_plus_1));
  result = result + (nu_i - nu_i_plus_1)*boost::math::digamma(nu_i);
  result = result - (iota_i - iota_i_plus_1)*(nu_i/iota_i);



  std::cout << std::abs(result) << std::endl;  
  return std::abs(result);
}


template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<known_variance,plurality_type::univariate> >; }
{

  matrix mu_i = a.post.mu;
  matrix mu_i_plus_1 = b.post.mu;
  matrix sigma_i = a.post.sigma;
  matrix sigma_i_plus_1 = b.post.sigma;
  double sigma = a.noise_structure.sigma; 
  matrix I = matrix::Identity(sigma_i.rows(),sigma_i.cols());
  double result = (sigma_i_plus_1.inverse() * sigma_i - I).trace();
  std::cout << sigma_i_plus_1.inverse() << std::endl;

  
  result = result + (1/sigma*sigma) * ((mu_i_plus_1 - mu_i).transpose() * sigma_i_plus_1.inverse() * (mu_i_plus_1 - mu_i))(0,0);
  result = result + std::log((sigma_i_plus_1.inverse() * sigma_i).determinant());
  result = 0.5*result;

  std::cout << std::abs(result) << std::endl;
  return std::abs(result);
  return a.weight;
  
  return result;
}


template<typename noise, plurality_type plurality>
std::list<particle_type<noise,plurality> >& prune(std::list<particle_type<noise,plurality> >& particles,const int& max_num_particles)
{
  // return particles;
  if(particles.size() > max_num_particles && particles.size() > 3)
    {
      auto it_1 = particles.begin();
      it_1++;
      auto it_2 = particles.begin();
      it_2++;it_2++;
      auto it_n_minus_1 = particles.end();
      it_n_minus_1--;
      // generate total variations between adjacent particles
      std::list<double> total_variations;      
      std::transform(it_1,it_n_minus_1,it_2,std::back_inserter(total_variations),[](auto& a,auto& b){return total_variation(a,b);});
      // locate the minimum total variation
      auto it_min_total_variation = std::min_element(total_variations.begin(),total_variations.end());
      // locate the particle to be evicted from the particle population
      auto it_evicted = particles.begin();
      std::advance(it_evicted,std::distance(total_variations.begin(),it_min_total_variation));
      auto it_relocation = it_evicted;
      it_relocation++;
      // update ratios
      auto combined_weight = it_relocation -> weight + it_evicted -> weight;
      auto weight = it_evicted -> weight;
      std::transform(it_evicted->ratios.begin(),
		     it_evicted->ratios.end(),
		     it_evicted->ratios.begin(),
		     [&combined_weight,&weight](auto& ratio){return ratio*weight/combined_weight;});
      weight = it_relocation -> weight;
      std::transform(it_relocation->ratios.begin(),
		     it_relocation->ratios.end(),
		     it_relocation->ratios.begin(),
		     [&combined_weight,&weight](auto& ratio){return ratio*weight/combined_weight;});
      it_relocation->ratios.insert(it_relocation->ratios.begin(),it_evicted->ratios.begin(),it_evicted->ratios.end());
      // update weight
      it_relocation -> weight += it_evicted -> weight;
      // evict the pruned particle
      particles.erase(it_evicted);
      std::cout << particles.size() << std::endl;
    }  
  return particles;
  
}




/*

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
				     [](auto& a,auto& b){return compare(a,b);});
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

*/








#endif
