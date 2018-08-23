/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	 //TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	 //  x, y, theta and their uncertainties from GPS) and all weights to 1. 
	 // Add random Gaussian noise to each particle.
	 // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    default_random_engine gen;
    num_particles = 100;
    
    std::normal_distribution<double> xd{x, std[0]};
    std::normal_distribution<double> yd{y, std[1]};
    std::normal_distribution<double> thetad{theta, std[2]};

    for (int i = 0; i < num_particles; i++)
    {
        Particle particle;
        particle.id = i;
        particle.x = xd(gen);
        particle.y = yd(gen);
        particle.theta = thetad(gen);
        particle.weight = 1.0;

        particles.push_back(particle); 
        weights.push_back(1);
    } 

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	 //TODO: Add measurements to each particle and add random Gaussian noise.
	 ///NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	 //http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	 //http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;
    num_particles = 100;

    std::normal_distribution<double> xd(0, std_pos[0]);
    std::normal_distribution<double> yd(0, std_pos[1]);
    std::normal_distribution<double> thetad(0, std_pos[2]);

    for (int i = 0; i < num_particles; ++i) 
    {
    
        if (abs(yaw_rate) != 0) {
        // Add measurements to particles
            particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
            particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
            particles[i].theta += yaw_rate * delta_t;
      
        } else {
           // Add measurements to particles
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
           // Theta will stay the same due to no yaw_rate      
    }

    // Add noise to the particles
    particles[i].x += xd(gen);
    particles[i].y += yd(gen);
    particles[i].theta += thetad(gen);
    
  }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	 //TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	 //  observed measurement to this particular landmark.
	 //NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	 //  implement this method and use it as a helper during the updateWeights phase.

     for (unsigned int i = 0; i < observations.size(); i++)
     {
         for (unsigned int j = 0; j < predicted.size(); j++)
	 {
             
         }
     }
 
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  // First term of multi-variate normal Gaussian distribution calculated below
  // It stays the same so can be outside the loop
  const double a = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
  
  // The denominators of the mvGd also stay the same
  const double x_denom = 2 * std_landmark[0] * std_landmark[0];
  const double y_denom = 2 * std_landmark[1] * std_landmark[1];

  // Iterate through each particle
  for (int i = 0; i < num_particles; ++i) {
    
    // For calculating multi-variate Gaussian distribution of each observation, for each particle
    double mvGd = 1.0;
    
    // For each observation
    for (unsigned int j = 0; j < observations.size(); ++j) {
      
      // Transform the observation point (from vehicle coordinates to map coordinates)
      double trans_obs_x, trans_obs_y;
      trans_obs_x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
      trans_obs_y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
      
      // Find nearest landmark
      vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
      vector<double> landmark_obs_dist (landmarks.size());
      for (unsigned int k = 0; k < landmarks.size(); ++k) {
        
        // Down-size possible amount of landmarks to look at by only looking at those in sensor range of the particle
        // If in range, put in the distance vector for calculating nearest neighbor
        double landmark_part_dist = sqrt(pow(particles[i].x - landmarks[k].x_f, 2) + pow(particles[i].y - landmarks[k].y_f, 2));
        if (landmark_part_dist <= sensor_range) {
          landmark_obs_dist[k] = sqrt(pow(trans_obs_x - landmarks[k].x_f, 2) + pow(trans_obs_y - landmarks[k].y_f, 2));

        } else {
          // Need to fill those outside of distance with huge number, or they'll be a zero (and think they are closest)
          landmark_obs_dist[k] = 999999.0;
          
        }
        
      }
      
      // Associate the observation point with its nearest landmark neighbor
      int min_pos = distance(landmark_obs_dist.begin(),min_element(landmark_obs_dist.begin(),landmark_obs_dist.end()));
      float nn_x = landmarks[min_pos].x_f;
      float nn_y = landmarks[min_pos].y_f;
      
      // Calculate multi-variate Gaussian distribution
      double x_diff = trans_obs_x - nn_x;
      double y_diff = trans_obs_y - nn_y;
      double b = ((x_diff * x_diff) / x_denom) + ((y_diff * y_diff) / y_denom);
      mvGd *= a * exp(-b);
      
    }
    
    // Update particle weights with combined multi-variate Gaussian distribution
    particles[i].weight = mvGd;
    weights[i] = particles[i].weight;

  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> new_particles (num_particles);
  
  // Use discrete distribution to return particles by weight
  random_device rd;
  default_random_engine gen(rd());
  for (int i = 0; i < num_particles; ++i) {
    discrete_distribution<int> index(weights.begin(), weights.end());
    new_particles[i] = particles[index(gen)];
    
  }
  
  // Replace old particles with the resampled particles
  particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
