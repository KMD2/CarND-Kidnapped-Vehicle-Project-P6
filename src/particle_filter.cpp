/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::tuple;
using std::normal_distribution;

// Initialize a random engine
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

  
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  // Generating random Gaussian noise
  std::normal_distribution<double> dist_x (x,std[0]);
  std::normal_distribution<double> dist_y (y,std[1]);
  std::normal_distribution<double> dist_theta (theta,std[2]);
  
  for (int i = 0; i < num_particles; i++ ){
  	// Create a particle
    Particle particle;
    
    // Assign ID to the particle
    particle.id = i;
    
    // Sample x,y and theta from the above normal distributions
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    
    // Initialize the particle weight with 1
    particle.weight = 1.0;
    weights.push_back(1.0);
    
    particles.push_back(particle);
  }
  
  // Change the flag to true since the vector particles has been initialized  
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  for (int i = 0; i < num_particles; i++){
    
  if (yaw_rate == 0){
    // Predicte x and y and update them accordingly. Since the yaw rate = 0, the new yaw = old yaw 
    particles[i].x += velocity * delta_t * cos(particles[i].theta);
    particles[i].y += velocity * delta_t * sin(particles[i].theta); 
}
   else
   {
     // Predicte x,y and theta and update them accordingly
     particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
     particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
     particles[i].theta += yaw_rate * delta_t;    
  }
      
   // Generating random Gaussian noise having the means equal to the updated particles position and the standard deviations as the one given in std_pos  
  std::normal_distribution<double> dist_x (particles[i].x,std_pos[0]);
  std::normal_distribution<double> dist_y (particles[i].y,std_pos[1]);
  std::normal_distribution<double>  dist_theta(particles[i].theta,std_pos[2]);
    
  // Updated particles' position by the noise added 
  particles[i].x = dist_x(gen);
  particles[i].y = dist_y(gen);
  particles[i].theta = dist_theta(gen);
  
  }
  
 
}
    

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
     

    // Loop over all the observed landmarks from the LIDAR sensor
  int id = 0;
  
  for (int i = 0; i < observations.size(); i++){
      
      // Initialize a the min_distance variable with infinity, to be overwritten later on by the minimum distance after comparing between the current observed 			landmark and all the predicted ones.  
      double min_distance = std::numeric_limits<double>::infinity();
    
    // Loop over the predicted landmarks
      for (int j = 0; j < predicted.size(); j++){
        
        double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
        
       	// Overwrites the min_distance if a closer predicted landmark is found, and the id is apdated
        if (min_distance > distance){
          
          min_distance = distance;
          
          id = predicted[j].id;
        }
      }
        
        // Assign the id of the nearest predicted landmark to the current observed landmark (Associate) 
        observations[i].id = id;
      }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // For each particle, get the predicted landmarks within the sensor range
  

  LandmarkObs associated_lm;
  double lm_x, lm_y, obs_x, obs_y;
  int lm_id, obs_id;
  
  
  for (int i = 0; i < particles.size(); i++){
    
    
      weights[i] = 1.0;
      vector <LandmarkObs> predicted;
      vector <LandmarkObs> transformed_obs;
    
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      lm_x =  map_landmarks.landmark_list[j].x_f ;
      lm_y =  map_landmarks.landmark_list[j].y_f ;
      lm_id = map_landmarks.landmark_list[j].id_i ;
      
      if (dist(particles[i].x, particles[i].y, lm_x, lm_y) < sensor_range)
        predicted.push_back(LandmarkObs {lm_id, lm_x, lm_y});
    }
    

  
    
    // Transform the observed landmarks from the vehicle's coordinates to the map's coordinates
    for (int j = 0; j < observations.size(); j++){
        
    obs_x = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
    obs_y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
    obs_id = observations[i].id; 
      
    transformed_obs.push_back(LandmarkObs {obs_id, obs_x, obs_y});     
    }
    

  
  // Apply association 
  dataAssociation(predicted, transformed_obs);

    
    vector<int> associations;
	vector<double> sense_x;
    vector<double> sense_y;
  
  // Calculate the particle final weight
for (int j = 0; j < transformed_obs.size(); j++){
        
  /// Find the associated landmark 
  for (int k = 0; k < predicted.size(); k++){
  
  	if (transformed_obs[j].id == predicted[k].id){
    	associated_lm = predicted[k];
    	break;
    }
    
  }
  
  /// Parameters to pass to SetAssociations for debugging
  associations.push_back(associated_lm.id);
  sense_x.push_back(transformed_obs[j].x);
  sense_y.push_back(transformed_obs[j].y);
  
  
  /// Calculate the  normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);

  /// Calculate the exponent
  double exponent;
  exponent = (pow(transformed_obs[j].x - associated_lm.x, 2) / (2 * pow(std_landmark[0], 2)))
               + (pow(transformed_obs[j].y - associated_lm.y, 2) / (2 * pow(std_landmark[1], 2)));
    
  /// Calculate weight using normalization terms and exponent

  weights[i] *= gauss_norm * exp(-exponent);
  particles[i].weight = weights[i];                                                             
  
 }
    /// Call SetAssociations for debugging
    SetAssociations(particles[i], associations, sense_x, sense_y);
 }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  /// Declare a vector of the size of number of particles to store the resampled particles in 
  vector<Particle> particles_resampled (num_particles);
  
  /// Produces random integers on the interval of the begining and end of the weights vector
  std::discrete_distribution<> d(weights.begin(), weights.end());
  for (int i = 0; i < particles.size(); i++){  
    particles_resampled[i] = particles[d(gen)];
  }
  particles = particles_resampled;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}