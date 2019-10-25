/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 * Modified by Shin-Ying Lu on 12/23/2019 for Udacity class hw
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
using std::normal_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  // Set the number of particles
  num_particles = 100;

  // Initilaize particles with Gaussian distribution
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);  

  for (int i = 0; i < num_particles; i++) {
      
    Particle p;
    p.id = i;
    p.x = dist_x(gen); 
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p); 
  }

  is_initialized = true;    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);



  for (int i = 0; i < num_particles; i++) {
    if (fabs(yaw_rate) < 1E-5) {  // prerventing from devided by zero if yaw rate is zero
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } 
    else {

      particles[i].x += velocity / yaw_rate * ( sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta) );
      particles[i].y += velocity / yaw_rate * ( cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t) );
      particles[i].theta += yaw_rate * delta_t;
    }

    // adding sensor noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (unsigned i = 0; i < observations.size(); i++){
    // initialize min_distance and 
    double min_distance = std::numeric_limits<double>::max(); 
    // initialize nearst neighbor map id
    int map_id = -1;

    for (unsigned j = 0; j < predicted.size(); j++){
      double distance_ij = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if (distance_ij < min_distance){
        min_distance = distance_ij;
        map_id = predicted[j].id;
      }
      observations[i].id = map_id; 

    }
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

  double sensor_range_2 = sensor_range * sensor_range;
  
  for (int i = 0; i < num_particles; i++) { //loop over particles
    
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;

    vector<LandmarkObs> predictions;
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {//loop over landmarks

      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;

      if ( pow(p_x-lm_x,2) + pow(p_y-lm_y,2) <= sensor_range_2 ) {
        predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
      }
    }// loop over landmarks 

    // transform observations from car coordinate to map coordinate
    vector<LandmarkObs> transformed_obs;
    for (unsigned int j = 0; j < observations.size(); j++) {
      double trans_x = cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y + p_x;
      double trans_y = sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y + p_y;
      transformed_obs.push_back(LandmarkObs{ observations[j].id, trans_x, trans_y });
    }

    dataAssociation(predictions, transformed_obs);

    // initialize weight
    particles[i].weight = 1.0;

      for (unsigned int j = 0; j < transformed_obs.size(); j++) {// loop over transformed_obs
        
        double obs_x = transformed_obs[j].x;
        double obs_y = transformed_obs[j].y;
        double lmpre_x, lmpre_y;
        int associated_prediction = transformed_obs[j].id;

        // get the x,y coordinates of the prediction associated with the current observation
        for (unsigned int k = 0; k < predictions.size(); k++) {// loop over predictions
          if (predictions[k].id == associated_prediction) {
            lmpre_x = predictions[k].x;
            lmpre_y = predictions[k].y;
          }
        }// loop over predictions

        double sig_x = std_landmark[0];
        double sig_y = std_landmark[1];
        // calculate weight - normalization term
        double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
        // calculate weight - exponent
        double exponent = (pow(obs_x - lmpre_x, 2) / (2 * pow(sig_x, 2)))+ (pow(obs_y - lmpre_y, 2) / (2 * pow(sig_y, 2)));
        // calculate the weight
        double weight = gauss_norm * exp(-exponent);
        
        particles[i].weight *= weight;

      }// loop over transformed_obs

  } // loop over particles
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // Find max(weight)
  vector<double> weights;
  double maxWeight = std::numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > maxWeight ) {
      maxWeight = particles[i].weight;
    }
  }
  

  // the magic wheel algo ... 
  vector<Particle> resample_particles; 
  std::uniform_int_distribution<int> distInt(0, num_particles - 1);
  int index = distInt(gen);
  std::uniform_real_distribution<double> distDouble(0.0, maxWeight);
  double beta = 0.0;
  for (int i = 0; i < num_particles; i++) {
    beta += distDouble(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resample_particles.push_back(particles[index]);
  }

  particles = resample_particles;

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