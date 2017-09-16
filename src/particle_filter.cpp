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

#define DEBUGGING (1)

using namespace std;

void DebugPrint(const char* message)
{
	if(DEBUGGING){
		std::cout << "DEBUG: " << message << std::endl;
	}
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	DebugPrint("init: Entering.");
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	default_random_engine gen;

	this->num_particles = 4;
	this->is_initialized = true;
	this->weights.resize( this->num_particles, 1.0 );
	this->particles.resize( this->num_particles );

	// Creates normal (Gaussian) distributions for x, y, theta.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[0]);

	for(int ii = 0; ii < this->num_particles ; ++ii ){

		this->weights[ii] = 1.0;

		this->particles[ii].x = dist_x( gen );
		this->particles[ii].y = dist_y( gen );
		this->particles[ii].theta = dist_theta( gen );

	}

	DebugPrint("init: Leaving.");
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	DebugPrint("prediction: Entering.");
	// Add measurements to each particle and add random Gaussian noise.

	// Generator for the gaussian distributions.
	default_random_engine gen;
	// Creates normal (Gaussian) distributions for x, y, theta.
	normal_distribution<double> dist_x(     0.0, std_pos[0] );
	normal_distribution<double> dist_y(     0.0, std_pos[1] );
	normal_distribution<double> dist_theta( 0.0, std_pos[2] );

	for(int ii = 0; ii < this->num_particles ; ++ii ){

		double x_0 = this->particles[ii].x;

		// Apply the state update to the particle (the formula depends on the yaw_rate being 0 or not)
		if(  yaw_rate < 1e-5 ){
			this->particles[ii].x += velocity * delta_t * cos( this->particles[ii].theta );
			this->particles[ii].y += velocity * delta_t * sin( this->particles[ii].theta );
		}
		else{
			double theta_0 = this->particles[ii].theta;
			double theta_f = theta_0 + delta_t * yaw_rate;

			this->particles[ii].x += velocity/yaw_rate * ( sin( theta_f) - sin( theta_0 ) );
			this->particles[ii].y += velocity/yaw_rate * ( cos( theta_0) - cos( theta_f ) );
			this->particles[ii].theta = theta_f;
		}

		// Add the noise
		this->particles[ii].x += dist_x( gen );
		this->particles[ii].y += dist_y( gen );
		this->particles[ii].theta += dist_theta( gen );
	}

	DebugPrint("prediction: Leaving.");
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predictions, std::vector<LandmarkObs>& observations)
{
	DebugPrint("dataAssociation: Entering.");
	// Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//	predictions <- landmarks_within_sensor_range
	//	observations <- observations_in_map_coordinates

	// TODO
	for (auto & ob: observations) {
		auto closest = min_element(begin(predictions), end(predictions),
								   [& ob] (const LandmarkObs &a, const LandmarkObs &b)
								   { return dist(a.x, a.y, ob.x, ob.y) < dist(b.x, b.y, ob.x, ob.y); } );
		ob.id = distance(begin(predictions), closest);
	}

	DebugPrint("dataAssociation: Leaving.");
}

double MultiVariateGaussianProbability(double x, double y, double mu_x, double mu_y, double std_x, double std_y)
{

	double exponent =   pow(x - mu_x, 2) / (2*pow(std_x, 2))
					  + pow(y - mu_y, 2) / (2*pow(std_y, 2)) ;

	double prob = 1.0 / (2*M_PI*std_x*std_y) * exp( -exponent );

    return prob;
}

LandmarkObs ParticleFilter::convertObservationToMapCoordinates( Particle particle, LandmarkObs car_observation )
{
	// Landmark position in the map (global coordinates)
	LandmarkObs observation_in_map_coordinates;

	observation_in_map_coordinates.id = car_observation.id ; // Landmark ID
	observation_in_map_coordinates.x = particle.x + car_observation.x * cos( particle.theta ) - car_observation.y * sin( particle.theta );
	observation_in_map_coordinates.y = particle.y + car_observation.x * sin( particle.theta ) + car_observation.y * cos( particle.theta );

	return observation_in_map_coordinates;
}

void ParticleFilter::updateWeights(	double sensor_range, double std_landmark[],
									std::vector<LandmarkObs> observations, Map map_landmarks)
{
	DebugPrint("updateWeights: Entering.");
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

    vector<LandmarkObs> landmarks_within_sensor_range;
	// Generate predicted landmark observations (convert each observation to map coordinates)
	std::vector<LandmarkObs> observations_in_map_coordinates;
	observations_in_map_coordinates.resize( observations.size() );

	// For each particle (the "candidate" in each iteration)
	for( int ii = 0; ii < this->num_particles; ++ii ){
		std::cout  << "ii " << ii << std::endl;
		DebugPrint("updateWeights: Start particle cycle.");

		// Clear vector if used in previous cycles.
	    if(landmarks_within_sensor_range.size() > 0){
	    	std::cout << "Clear vector of size " << landmarks_within_sensor_range.size() << std::endl;
	    	landmarks_within_sensor_range.clear();
	    	std::cout << "Cleared vector, now of size " << landmarks_within_sensor_range.size() << std::endl;
	    }

		double xp = particles[ii].x;
		double yp = particles[ii].y;

		// Generate predicted landmark observations (convert each observation to map coordinates)

		for( int il = 0 ; il < observations.size() ; ++il){
			observations_in_map_coordinates[ii] = this->convertObservationToMapCoordinates( particles[ii], observations[il] );
		}

		// Limit the landmarks to be used to those in the range of the sensor from the candidate particle position.
		DebugPrint("updateWeights: Limit the landmarks to be used.");

	    for ( int ilm = 0 ; ilm < map_landmarks.landmark_list.size() ; ++ilm ) {
			double xm = (double)map_landmarks.landmark_list[ilm].x_f;
			double ym = (double)map_landmarks.landmark_list[ilm].y_f;

			if (dist(xm, ym, xp, yp) <= sensor_range) {
				LandmarkObs obs;
				obs.x = xm;
				obs.y = ym;
				obs.id = map_landmarks.landmark_list[ilm].id_i;
				DebugPrint("updateWeights: New landmarks_within_sensor_range.");
				std::cout  << " ilm " << ilm << " x " << obs.x << " y " << obs.y << " id " << obs.id << " size_vector " << landmarks_within_sensor_range.size() << std::endl;
				landmarks_within_sensor_range.push_back(obs);
				std::cout  << "updateWeights: New landmarks_within_sensor_range. Done" << std::endl;
			}
	    }


		DebugPrint("updateWeights: If there are landmarks in range.");
	    if (landmarks_within_sensor_range.size() != 0) {
			// Select landmark associated to the specific observation
			this->dataAssociation( landmarks_within_sensor_range, observations_in_map_coordinates);


			DebugPrint("updateWeights: Calculate final weight.");
			// Calculate final weight using Multivariate-Gaussian probability
			double current_weight = 1;
			for (auto & observation : observations_in_map_coordinates) {
				double xo = observation.x;
				double yo = observation.y;
				double xr = landmarks_within_sensor_range[observation.id].x;
				double yr = landmarks_within_sensor_range[observation.id].y;
				//MultiVariateGaussianProbability(double x, double y, double mu_x, double mu_y, double std_x, double std_y)

				DebugPrint("updateWeights: MultiVariateGaussianProbability.");
				current_weight *= MultiVariateGaussianProbability(xo, yo, xr, yr, std_landmark[0], std_landmark[1]);
			}

			DebugPrint("updateWeights: Save current weight.");
			this->weights[ii] = current_weight;
			DebugPrint("updateWeights: Save current weight. Done");

	    }
	    else{
	    	// No weight update
	    }

		DebugPrint("updateWeights: Loop on particles - End iteration.");
	}


	// Normalize to weights


	DebugPrint("updateWeights: Leaving.");
}

void ParticleFilter::resample()
{
	DebugPrint("resample: Entering.");
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> particle_index_generator(weights.begin(), weights.end());

	vector<Particle> particles_copy = this->particles;
	for (int ii = 0; ii < num_particles; ++ii ){
		this->particles[ii] = particles_copy[particle_index_generator(gen)];
	}

	DebugPrint("resample: Leaving.");
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
