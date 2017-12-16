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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	if (!is_initialized) {
		cout << "x: " << x << ", y: " << y << ", theat: " << theta << ", std: " ;
		num_particles = 10;
		// debug
		// num_particles = 1;
		for(int iter = 0; iter < num_particles; iter++) {
			// cout << "particles[iter] : " << particles[iter];
			Particle part = Particle();
			part.id = iter;
			part.x = x;
			part.y = y;
			particles.push_back(part);
		}
		is_initialized = true;
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	for (int i = 0; i < num_particles; ++i) {
		Particle part = particles[i];
		
		// 1. Add measurements to each particle
		// codes modified from lecture: https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/2c318113-724b-4f9f-860c-cb334e6e4ad7/lessons/5c50790c-5370-4c80-aff6-334659d5c0d9/concepts/56d08bf5-8668-42e7-a718-1ef40d444259
		double partx = part.x;
		double party = part.y;
		double parttheta = part.theta;
		partx = partx + (velocity / yaw_rate) * (sin(parttheta + yaw_rate * delta_t) - sin(parttheta));
		party = party + (velocity / yaw_rate) * (-cos(parttheta + yaw_rate * delta_t) + cos(parttheta));
		parttheta = parttheta + yaw_rate;

		// 2. and add random Gaussian noise.
		// codes modified from lecture: https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/2c318113-724b-4f9f-860c-cb334e6e4ad7/lessons/5c50790c-5370-4c80-aff6-334659d5c0d9/concepts/53081ef1-a14c-4ae9-a68f-3ddc258cfd95
		// std_pos[]: 0.3, 0.3, 0.01
		normal_distribution<double> dist_x(partx, std_pos[0]);
		normal_distribution<double> dist_y(party, std_pos[1]);
		normal_distribution<double> dist_theta(parttheta, std_pos[2]);
		double sample_x = dist_x(gen);
		double sample_y = dist_y(gen);
		double sample_theta = dist_theta(gen);
		particles[i].x = sample_x;
		particles[i].y = sample_y;
		particles[i].theta = sample_theta;
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

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

	bool bdebug = true;
	bdebug = false;
	// sensor_range = 50;
	// sigma_landmark [2] = {0.3, 0.3};
	cout << "noisy_observations size: " << observations.size() << endl;
	cout << "len(map_landmarks) : " <<  map_landmarks.landmark_list.size() << endl;
	// 10
	// 42

	// algorithm steps summary:
	// 1. For each particle
	// 2. take a list of particle's observations and transform them into maps coordinates
	// 3. For each transformed observation you calculate a distance to all landmarks:
	// 4. vector dists contains distances to all landmarks for single observation. 
	//    You choose the minimum distance and associate id of landmark for this distance to the observation.
	// 5. when your observation has an associated landmark, you're calculating 
	//    the probability that particular particle's observation saw this particular landmark.
	// 6. particle's total weight is a product of probabilities of each observations

	for(int ni = 0; ni < num_particles; ni++) {
		// 2. take a list of particle's observations and transform them into maps coordinates
		// homogeneous transformation matrix:
		// cos(theta), -sin(theta), Xp 
		// sin(theta),  cos(theta), Yp
		//          0,           0,  1
		double theta = particles[ni].theta;
		double Xp = particles[ni].x;
		double Yp = particles[ni].y;
		const int nobssize = observations.size();
		
		// cout << "Xp : " << Xp << " Yp : " << Yp << endl;
		std::vector<LandmarkObs> observations_transf;
		particles[ni].weight = 1;
		for(int nj = 0; nj < nobssize; nj++) {
			double Xc = observations[nj].x;
			double Yc = observations[nj].y;
			double id = observations[nj].id;

			double tranfx = Xc * cos(theta) + Yc * (-sin(theta)) + Xp;
			double tranfy = Xc * sin(theta) + Yc * ( cos(theta)) + Yp;
			
			LandmarkObs obs_t;
			obs_t.id = observations[nj].id;
			obs_t.x = tranfx;
			obs_t.y = tranfy;
			observations_transf.push_back(obs_t);

			// 3. For each transformed observation you calculate a distance to all landmarks:
			vector <double> dists;
			for (Map::single_landmark_s slm : map_landmarks.landmark_list) {
				double dst = dist(slm.x_f, slm.y_f, tranfx, tranfy);
				dists.push_back(dst);
			}

			// 4. vector dists contains distances to all landmarks for single observation. 
			//    You choose the minimum distance and associate id of landmark for this distance to the observation.
			vector<double>::iterator result = min_element(begin(dists), end(dists));
			Map::single_landmark_s lm = map_landmarks.landmark_list[distance(begin(dists), result)];
			obs_t.id = lm.id_i;

			// 5. when your observation has an associated landmark, you're calculating 
			//    the probability that particular particle's observation saw this particular landmark.

			// 6. particle's total weight is a product of probabilities of each observations
			double sig_x = std_landmark[0];
			double sig_y = std_landmark[1];
			double gauss_norm = (1/(2 * M_PI * sig_x * sig_y) );
			double exponent= pow((obs_t.x - lm.x_f), 2) / (2 * pow(sig_x, 2) )+ 
				pow((obs_t.y - lm.y_f), 2) / (2 * pow(sig_y, 2) );
			double weight = gauss_norm * exp(-exponent);
			// cout << "weight : " << weight << endl;
			// maybe I should only multiple with observations of "obs_t.id" ?
			if (weight > 0.0f) {
				particles[ni].weight = particles[ni].weight * weight;
				// cout << "in loop : particles[ni].weight : " << particles[ni].weight << endl;
			}
		}
		cout << "particles[ " << ni << "].weight : " << particles[ni].weight << endl;

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	bool bdebug = true;
	bdebug = false;
	// codes from Lesson 13: Particle Filters "19. Quiz: New Particle".
	// translate python codes to cpp.
	// i = 1
	// while i <=1000:
	// 	x = random.randint(1,1000) 
	// 	xw = random.random() 
	// 	# if this particle is lucky enough to be >= xw.
	// 	if w[x-1]/wsum >= xw:
	// 		p3.append(p[x-1])
	// 		i = i + 1
	double wsum = 0;
	for (int it = 0; it < num_particles; it++) {
		wsum = wsum + particles[it].weight;
	}
	cout << "wsum : " << wsum << endl;
	
	std::vector<Particle> new_particles;
	int ni = 1;
	while(ni <= num_particles) {
		double idx = rand() % num_particles + 1;
		double thresholdweight = ((double) rand() / (RAND_MAX));
		double prob = particles[idx].weight / wsum;
		if (bdebug) {
			cout << " thresholdweight : " << thresholdweight << endl;
			cout << " idx : " << idx << endl;
			cout << " particles[idx].weight  : " << particles[idx].weight << endl;
			cout << " particles[idx].weight / wsum : " << prob << endl;
		}
		cout << " idx : " << idx << endl;
	 	// if this particle is lucky enough to be >= xw.
		if (prob > thresholdweight) {
			if (bdebug) {
				cout << " particles[idx - 1].weight : " << particles[idx - 1].weight << endl;
			}
			new_particles.push_back(particles[idx - 1] );
			ni = ni + 1;
			cout << "ni :" << ni << endl;
		} else {
			cout << "smaller, prob :"  << prob << " thresholdweight : "<< thresholdweight << endl;
			
		}
	}
	cout << "before particles = new_particles; " << endl;
	particles = new_particles;

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
