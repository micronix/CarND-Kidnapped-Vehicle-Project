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

#define EPS 0.0001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;

	particles = std::vector<Particle>(num_particles);

	std::default_random_engine gen;
	std::normal_distribution<double> x_distribution(x, std[0]);
	std::normal_distribution<double> y_distribution(y, std[1]);
	std::normal_distribution<double> theta_distribution(theta, std[2]);

	int id = 0;
	for(auto& particle : particles){
		particle.id = id;
		particle.x = x_distribution(gen);
		particle.y = y_distribution(gen);
		particle.theta = theta_distribution(gen);
		particle.weight = 1.0;
		id++;
	}

	weights = std::vector<double>(num_particles);
	for(int i = 0; i < num_particles; i++){
		weights[i] = 1.0;
	}

	is_initialized = true;

}

/*
void ParticleFilter::PrintState(){
	cout << "particles: " << endl;
	for(auto& particle : particles){
		cout << particle.id << ": " << particle.x << ", " << particle.y << ", " << particle.theta << ", " << particle.weight << endl;
	}

	cout << "weights: " << endl;
	for(int i = 0; i < num_particles; i++){
		cout << weights[i] << endl;
	}
}*/


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//PrintState();

	std::default_random_engine gen;

	double dist = velocity * delta_t;
	double vy = velocity / yaw_rate;
	double delta_theta = yaw_rate * delta_t;

	for(auto& particle : particles){
		if (fabs(yaw_rate) < EPS){
			particle.x += dist * cos(particle.theta);
			particle.y += dist * sin(particle.theta);
		} else {
			double theta = delta_theta  + particle.theta;
			particle.x += vy * (sin(theta) - sin(particle.theta));
			particle.y += vy * (cos(particle.theta) - cos(theta));
			particle.theta += delta_theta;
		}

		// add noise
		particle.x = std::normal_distribution<double>(particle.x, std_pos[0])(gen);
		particle.y = std::normal_distribution<double>(particle.y, std_pos[1])(gen);
		particle.theta = std::normal_distribution<double>(particle.theta, std_pos[2])(gen);
	}

	//PrintState();
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


}

/*
std::vector<Map::single_landmark_s> ParticleFilter::nearbyLandmarks(Map map, double sensor_range, double x, double y){
	std::vector<Map::single_landmark_s> landmarks;
	for(auto& landmark : map.landmark_list){
		if (landmark.x_f > x - sensor_range && landmark.x_f < x + sensor_range && landmark.y_f > y - sensor_range && landmark.y_f < y + sensor_range){
			landmarks.push_back(landmark);
		}
	}
	return landmarks;
}*/

/*std::vector<LandmarkObs> ParticleFilter::transformedObservations(std::vector<LandmarkObs> observations, double x, double y, double theta){
	std::vector<LandmarkObs> tobservations(observations.size());
	//cout << "transforming: " << endl;
	for(int i = 0; i < observations.size(); i++){
		tobservations[i].x = x + observations[i].x * cos(theta) - observations[i].y * sin(theta);
		tobservations[i].y = y + observations[i].x * sin(theta) + observations[i].y * cos(theta);
		//cout << observations[i].x << " -> " << tobservations[i].x << endl;
	}
	return tobservations;
}*/

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
	double std_x2 = std_landmark[0] * std_landmark[0];
	double std_y2 = std_landmark[1] * std_landmark[1];
	double coeff = 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);

	for(auto& particle : particles){
		//cout << "particle: " << particle.x << "," << particle.y << "," << particle.theta << " sensor_range: " << sensor_range << endl;
		double x = particle.x;
		double y = particle.y;
		double theta = particle.theta;

		// select only landmarks within sensor range
		std::vector<Map::single_landmark_s> landmarks;
		for(auto& landmark : map_landmarks.landmark_list){
			if (landmark.x_f > x - sensor_range && landmark.x_f < x + sensor_range && landmark.y_f > y - sensor_range && landmark.y_f < y + sensor_range){
				landmarks.push_back(landmark);
			}
		}

		// transform observations to coordinates of particle
		//std::vector<LandmarkObs> tobservations = transformedObservations(observations, particle.x, particle.y, particle.theta);
		std::vector<LandmarkObs> tobservations(observations.size());
		//cout << "transforming: " << endl;
		for(int i = 0; i < observations.size(); i++){
			tobservations[i].x = x + observations[i].x * cos(theta) - observations[i].y * sin(theta);
			tobservations[i].y = y + observations[i].x * sin(theta) + observations[i].y * cos(theta);
			//cout << observations[i].x << " -> " << tobservations[i].x << endl;
		}

		/*cout << " observations: " << endl;
		for(auto& observation : observations){
			cout << "  " << observation.x << ", " << observation.y << endl;
		}

		cout << " transformed observations: " << endl;
		for(auto& observation : tobservations){
			cout << "  " << observation.x << ", " << observation.y << endl;
		}*/
		// look for the closest landmark for each transformed observation
		vector<int> associations(tobservations.size());
		vector<double> sense_x(tobservations.size());
		vector<double> sense_y(tobservations.size());

		for(int i = 0; i < tobservations.size(); i++){
			double x = tobservations[i].x;
			double y = tobservations[i].y;
			double min_distance = std::numeric_limits<double>::max();
			int best_landmark_id = 0;
			double best_x = 0.0;
			double best_y = 0.0;

			for(auto& landmark : landmarks){
				double distance = dist(x,y,landmark.x_f, landmark.y_f);
				if (distance < min_distance){
					min_distance = distance;
					best_landmark_id = landmark.id_i;
					best_x = landmark.x_f;
					best_y = landmark.y_f;
				}
			}
			associations[i] = best_landmark_id;
			sense_x[i] = best_x;
			sense_y[i] = best_y;
		}

		particle = SetAssociations(particle, associations, sense_x, sense_y);

		double weight = 1.0;
		for(int i = 0; i < associations.size(); i++){
			double x = tobservations[i].x;
			double y = tobservations[i].y;
			double e = pow(x - sense_x[i], 2.0) / (2 * std_x2) + pow(y - sense_y[i], 2.0) / (2 * std_y2);
			double res = exp(-e);
			weight *= res;
		}
		particle.weight = weight;
		weights[particle.id] = weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine gen;
	std::discrete_distribution<> dist(weights.begin(), weights.end());

	std::vector<Particle> sampled_particles(particles.size());
	// cout << "sampled particles: ";
	for(int i = 0; i < particles.size(); i++){
		int index = dist(gen);
		// cout << " " << index;
		sampled_particles[i] = particles[index];
		sampled_particles[i].id = i;
	}
	// cout << endl;
	particles = sampled_particles;
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
