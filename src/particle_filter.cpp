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
static default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//Set No. of particles
	this->num_particles = 100;

	std::default_random_engine gen;

	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);

	//Initilize all particles to the first position
	for (int i = 0; i < this->num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		this->particles.push_back(particle);
		weights.push_back(1);
	}

	//Add random Gaussian noise to each particle

	//Set initilization complete
	this->is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	for (int i = 0; i < num_particles; i++)
	{
		double new_x;
		double new_y;
		double new_theta;

		if (yaw_rate == 0)
		{
			new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			new_theta = particles[i].theta;
		}
		else
		{
			new_x = particles[i].x + ((velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)));
			new_y = particles[i].y + ((velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t)));
			new_theta = particles[i].theta + (yaw_rate*delta_t);
		}

		normal_distribution<double> N_x(new_x, std_pos[0]);
		particles[i].x = N_x(gen);

		normal_distribution<double> N_y(new_y, std_pos[1]);
		particles[i].y = N_y(gen);

		normal_distribution<double> N_theta(new_theta, std_pos[2]);
		particles[i].theta = N_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); i++)
	{
		LandmarkObs obv = observations[i];
		double minDist = numeric_limits<double>::max();
		int mapid = -1;

		for (int j = 0; j < predicted.size(); j++)
		{
			LandmarkObs prd = predicted[j];
			double currentDist = dist(obv.x, obv.y, prd.x, prd.y);
			if(currentDist<minDist)
			{
				minDist = currentDist;
				mapid = prd.id;
			}
		}
		observations[i].id = mapid;
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

	for (int i = 0; i < num_particles; i++) 
	{

		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		//Get landmarks that are within sensor range
		vector<LandmarkObs> predictions;
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) 
		{
			float x_landmark = map_landmarks.landmark_list[j].x_f;
			float y_landmark = map_landmarks.landmark_list[j].y_f;
			int id_landmark = map_landmarks.landmark_list[j].id_i;
			
			if (fabs(x_landmark - x) <= sensor_range && fabs(y_landmark - y) <= sensor_range) {
				predictions.push_back(LandmarkObs{ id_landmark, x_landmark, y_landmark });
			}
		}

		// list of observations transformed from vehicle coordinates to map coordinates
		vector<LandmarkObs> transformation;
		for (unsigned int j = 0; j < observations.size(); j++) {
			double t_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
			double t_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
			transformation.push_back(LandmarkObs{ observations[j].id, t_x, t_y });
		}

		//Perform dataAssociation
		dataAssociation(predictions, transformation);

		//reset weight
		particles[i].weight = 1.0;

		//Calculate weights
		for (unsigned int j = 0; j < transformation.size(); j++) {

			double x_obv, y_obv, x_pred, y_pred;
			x_obv = transformation[j].x;
			y_obv = transformation[j].y;

			int id_lm = transformation[j].id;

			for (unsigned int k = 0; k < predictions.size(); k++) {
				if (predictions[k].id == id_lm) {
					x_pred = predictions[k].x;
					y_pred = predictions[k].y;
				}
			}

			//calculate weight with multivariate Gaussian
			double s_x = std_landmark[0];
			double s_y = std_landmark[1];
			double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(x_pred-x_obv,2)/(2*pow(s_x, 2)) + (pow(y_pred-y_obv,2)/(2*pow(s_y, 2))) ) );

			particles[i].weight *= obs_w;
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//From Chapter: Particle filters, refering python code for resampling
	vector<Particle> p;
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	unsigned int index = uniintdist(gen);

	double beta = 0.0;

	vector<double> w;
	for (int i = 0; i < num_particles; i++) {
		w.push_back(particles[i].weight);
	}

	double mw = *max_element(w.begin(), w.end());

	uniform_real_distribution<double> unirealdist(0.0, mw);
	// spin the resample wheel!
	for (int i = 0; i < num_particles; i++) {
		beta += unirealdist(gen) * 2.0;
		while (beta > w[index]) {
			beta -= w[index];
			index = (index + 1) % num_particles;
		}
		p.push_back(particles[index]);
	}
	particles = p;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
	const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

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
