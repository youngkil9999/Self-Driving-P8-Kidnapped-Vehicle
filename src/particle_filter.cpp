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
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    num_particles = 5; // set number of particles

    Particle p;
    particles.resize(num_particles);
    weights.resize(num_particles);

    for (int i = 0; i < num_particles; ++i) {
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta =  dist_theta(gen);
        p.weight = 1.0;
        particles[i] = p;

    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;
    double x, y, theta;
    for (auto& p:particles){
        cout << "yaw_rate: " << yaw_rate << endl;
        if (fabs(yaw_rate) < 1e-4){
            yaw_rate = 1e-4;
            x = p.x + velocity/yaw_rate*(sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
            y = p.y + velocity/yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
            theta = yaw_rate * delta_t + p.theta;

        } else{
            x = p.x + velocity/yaw_rate*(sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
            y = p.y + velocity/yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
            theta = yaw_rate * delta_t + p.theta;
        }
        normal_distribution<double> dist_x(x, std_pos[0]);
        normal_distribution<double> dist_y(y, std_pos[1]);
        normal_distribution<double> dist_yaw(theta, std_pos[2]);

        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_yaw(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

//    find original landmark, which minimizes the distance between transformed obs and original landmark
//    predicted = inrange_landmark, observation = transformed obs

    for (auto& observation:observations){
        double  min_distance = 99999999.0;
        for (const auto& p:predicted) {
            double distance = dist(p.x, p.y, observation.x, observation.y);
            if (distance < min_distance){
                min_distance = distance;
                observation.id = p.id;
            }
        }

    }
}

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
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

//    vector<Map::single_landmark_s> landmark_list;
//    landmark_list.resize(map_landmarks.size());

    vector<LandmarkObs> transformed_obs;
    transformed_obs.resize(observations.size());

    LandmarkObs transformed_ob;

    for (int i = 0; i < num_particles; ++i) {

//      particle address allocated to each particles[i]
        Particle &particle = particles[i];

//      observation measurement(car coordinate) transformed to map coordinate, respect to each particle
        for (int j = 0; j < observations.size(); ++j) {

            transformed_ob.x = observations[j].x * cos(particle.theta) - observations[j].y * sin(particle.theta) + particle.x;
            transformed_ob.y = observations[j].x * sin(particle.theta) + observations[j].y * cos(particle.theta) + particle.y;

            transformed_obs[j] = transformed_ob;
//      map coordinate about landmark
        }

        map<int, Map::single_landmark_s> idxlandmarks;
        vector<LandmarkObs> inRange_landmark;

//      Check whether the landmark is in range (50)
//      for loop Range based
//      var landmark loop in landmark_list

        for (const auto& landmark:map_landmarks.landmark_list){

            double distance = dist(particle.x, particle.y, landmark.x_f, landmark.y_f);

            if (distance <= sensor_range) {
//              landmarkobs class add id, x, y
                inRange_landmark.push_back(LandmarkObs {landmark.id_i, landmark.x_f, landmark.y_f});
//                idxlandmarks.insert(std::make_pair(landmark.id_i, landmark));
            }
        }

        if (inRange_landmark.size() > 0 ){

            dataAssociation(inRange_landmark, transformed_obs);

//            calculate weight
            particle.weight = 1.0;

            for (const auto observation:transformed_obs) {
                int index = observation.id - 1;
                double mu_x = map_landmarks.landmark_list[index].x_f;
                double mu_y = map_landmarks.landmark_list[index].y_f;
                double x = observation.x;
                double y = observation.y;
                double s_x = std_landmark[0];
                double s_y = std_landmark[1];
                double x_diff = (x - mu_x) * (x - mu_x) / (2 * s_x * s_x);
                double y_diff = (y - mu_y) * (y - mu_y) / (2 * s_y * s_y);
                particle.weight *= 1 / (2 * M_PI * s_x * s_y) * exp(-(x_diff + y_diff));

            }
            weights[i] = particle.weight;

        }else{
            weights[i] = 0.0;
        }
    }

    }

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());
    vector<Particle> particles_new;

    particles_new.resize(num_particles);
    for (int n=0; n < num_particles; ++n) {
        particles_new[n] = particles[d(gen)];
    }

    particles = particles_new;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
