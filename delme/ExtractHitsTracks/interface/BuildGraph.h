
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "test/ExtractHitsTracks/interface/Data.h"
#include <vector>


template <class T>

std::vector<int> find_all_indices_in_vec(std::vector<T> vec, float num){ 
  std::vector<int> B; 
  typename std::vector<T>::iterator it = vec.begin();
  while ((it = std::find_if(it, vec.end(), [num](int x){return x == num; })) != vec.end())
  {
    B.push_back(std::distance(vec.begin(), it));
    it++;
  }
return B; 
}


template <class S>
std::vector<S> subset_vec(std::vector<S> vec, std::vector<int> indices){ 
  std::vector<S> result(indices.size(), 0); 
  std::transform(indices.begin(), indices.end(), result.begin(), [vec](size_t pos) {return vec[pos];});
  return result; 
}

 
class BuildGraph {

  public:
  BuildGraph() {}  



    
  void select_hits(Ntuple *nt){
  /* Performs pt cut 
   * If there are multiple hits in a layer for a given particle, it selects the one with smallest sim_dxy_sig (should be changed in the long run) 
   */
  
  
    // pt cut 
    
    std::vector<int> accepted_hit_indices;
    std::vector<int> new_accepted_hit_indices;  
    //should be a config variable instead
    int pt_cut = 2;  
  
    auto it = std::find_if(std::begin(nt->sim_pt_), std::end(nt->sim_pt_), [pt_cut](int i){return i > pt_cut;});
    while (it != std::end(nt->sim_pt_)) {
      accepted_hit_indices.emplace_back(std::distance(std::begin(nt->sim_pt_), it));
      it = std::find_if(std::next(it), std::end(nt->sim_pt_), [pt_cut](int i){return i > pt_cut;});
    }
  
  
    // remove duplicate hits in a layer for a given particle 
  
    std::vector<int> pt_cut_particle_id = subset_vec(nt->particle_id_,accepted_hit_indices); 
    // find unqiue particle ids 
    std::unordered_set<int> s(pt_cut_particle_id.begin(), pt_cut_particle_id.end());
    std::vector<int> unique_particle_ids; 
    unique_particle_ids.assign(s.begin(), s.end());  
    for(auto particle: pt_cut_particle_id){
      // find the positions in the list of that particle 
          
      //contains the indices of the partilce in the particle cut list 
      std::vector<int> positions = find_all_indices_in_vec(pt_cut_particle_id, particle); 
      
      //layer ids for a given particle 
      auto particle_layer_ids = subset_vec(nt->layer_id_,positions); 
      std::unordered_set<int> s(particle_layer_ids.begin(), particle_layer_ids.end());
      std::vector<int> unique_particle_layer_ids; 
      unique_particle_layer_ids.assign(s.begin(), s.end());  
      //check if the layer id is repeated 
      std::map<int, int> CountMap;
      for (auto layer = unique_particle_layer_ids.begin(); layer!= unique_particle_layer_ids.end(); layer++)
  	CountMap[*layer]++; 
      
      // for each layer for the particle 
      for (auto count = CountMap.begin(); count!=CountMap.end(); count++){
  	// if the layer id is repeated, find the sim_dxy_sig and return the index of the smallest values for this 
  	if (count->second > 1){
  	  int repeated_layer = count->first; 
  	  std::vector<int> indices_repeated_layers = find_all_indices_in_vec(particle_layer_ids, repeated_layer); 
  	  // the position in the data array will be the index for the particle plus the index for the layer 
            std::for_each(indices_repeated_layers.begin(), indices_repeated_layers.end(), [positions](int &d) { d+=positions.front();});
            
  	  std::vector<float> sim_dxy_sig_multiple_layer_hit = subset_vec(nt->sim_dxy_sig_,indices_repeated_layers); 
  	  //std::vector<int>::iterator smallest_sim_dxy_sig = std::min_element(sim_dxy_sig_for_multiple_layer_hit.begin(), sim_dxy_sig_for_multiple_layer_hit.end());
  	  float smallest_sim_dxy_sig = *std::min_element(sim_dxy_sig_multiple_layer_hit.begin(), sim_dxy_sig_multiple_layer_hit.end());
      	  auto it = std::find(sim_dxy_sig_multiple_layer_hit.begin(), sim_dxy_sig_multiple_layer_hit.end(), smallest_sim_dxy_sig);
  	  int index_min_sim_dxy_sig = std::distance(std::begin(sim_dxy_sig_multiple_layer_hit), it); 
  	  new_accepted_hit_indices.emplace_back(positions.front() + index_min_sim_dxy_sig); 	    	    
   	}     
     
      }//layer id counts for a given particle
   
    }// particle id 
  }
  };




