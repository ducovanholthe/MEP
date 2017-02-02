#include <iostream>
#include <cmath>
#include <cstdlib>
#include <boost/random.hpp>
#include <ctime>
#include <vector>
#include <string>
#include "print.h"

using namespace std;
const double pi = 4 * atan(1);

//Parameters
const double ki = 0.1; //internal spring constant
const double ko = 0.3; //overlap spring constant
const double dzeta = 1.0; //viscous drag coefficient
const double growthRate = 0.00005; //average growth rate, 1.23 per hour, division rate is about 45400/dt
const double diameter = 0.9; //diameter should equal 0.75 +- 0.0375 micrometer
const double startLength = 2.5;
const double dt = 1.0; //timestep
const int Nmax = 3000; //maximum amount of particles

					   //Random number generator
boost::mt19937 generator(time(0)); //number generator from random list
boost::normal_distribution<> //setup distributions
normalDistGrowth(0.0, 0.05 * growthRate), //growth distribution added to daughter cells, sigma = 0.277
normalDistAngle(0.0, 0.01), //Noise in orientation for daughter cells
normalDistLength(4.54, 0.2); //distribution of maximum length that particles can reach, 4.54, sigma = 0.46
boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
randomMu(generator, normalDistGrowth), //generates growth rate
randomTheta(generator, normalDistAngle), //generates noise in orientation
randomLmax(generator, normalDistLength); //generates division length

struct Particle
	//Structure of particles with its properties
{
	Particle(double, double, double, double, double, double); //constructor

	void grow();
	void forceInternal();
	void move();
	void clear();

	vector<int> neighbours;

	int ID;
	double D, L, mu, Lmax; //size properties: diameter, length, growth rate
	double xcpos, ycpos, theta; //spatial variables
	double x1pos, y1pos, x2pos, y2pos, xvel1, yvel1, xvel2, yvel2, Fx1, Fy1, Fx2, Fy2; //more spatial variables
	double oop; int ram; //sets orientational order parameter as proprty of particle
	int colour;
};

Particle::Particle(double xstart, double ystart, double angle, double length, double diameter, double growth)
{
	ID = 0;
	D = diameter;
	L = length; //rest length of particle (internal spring is at rest)
	mu = growth; //initial growth rate of particle
	Lmax = randomLmax();

	theta = angle;
	x1pos = xstart; y1pos = ystart;
	x2pos = x1pos + L * cos(theta); y2pos = y1pos + L * sin(theta);
	xcpos = (x2pos + x1pos) / 2; ycpos = (y2pos + y1pos) / 2;

	xvel1 = 0.0; yvel1 = 0.0; xvel2 = 0.0; yvel2 = 0.0; //velocities of both ends of particle
	Fx1 = 0.0; Fy1 = 0.0; Fx2 = 0.0; Fy2 = 0.0;

	oop = 0.0; ram = 0;
	colour = 240; //used for visualisation
}

void Particle::grow()
{
	L += (mu * dt);

	x1pos -= ((mu * dt) / 2) * cos(theta); y1pos -= ((mu * dt) / 2) * sin(theta);
	x2pos += ((mu * dt) / 2) * cos(theta); y2pos += ((mu * dt) / 2) * sin(theta);

	xcpos = (x2pos + x1pos) / 2; ycpos = (y2pos + y1pos) / 2;
}

void Particle::forceInternal()
{
	double f_x = ki * ((x2pos - x1pos) - L * cos(theta));
	double f_y = ki * ((y2pos - y1pos) - L * sin(theta));

	Fx1 += f_x; Fy1 += f_y;
	Fx2 -= f_x; Fy2 -= f_y;
}

void Particle::move()
{
	xvel1 = 2 * Fx1 / (dzeta * D); yvel1 = 2 * Fy1 / (dzeta * D); //calculate velocities
	xvel2 = 2 * Fx2 / (dzeta * D); yvel2 = 2 * Fy2 / (dzeta * D);

	x1pos += xvel1 * dt; y1pos += yvel1 * dt; //update positions with velocities
	x2pos += xvel2 * dt; y2pos += yvel2 * dt;

	theta = atan((y2pos - y1pos) / (x2pos - x1pos)); //update orientation
}

void Particle::clear()
{
	xvel1 = 0.0; yvel1 = 0.0; xvel2 = 0.0; yvel2 = 0.0; //velocities of both ends of particle
	Fx1 = 0.0; Fy1 = 0.0; Fx2 = 0.0; Fy2 = 0.0;
}

void divide(Particle &pOld, Particle &pNew)
{
	pNew.L = (pNew.L - pNew.D) / 2;
	pNew.x1pos = pOld.x2pos - pNew.L * cos(pNew.theta); pNew.y1pos = pOld.y2pos - pNew.L * sin(pNew.theta);
	pNew.mu += randomMu();
	pNew.theta += randomTheta();

	pOld.L = (pOld.L - pOld.D) / 2;
	pOld.x2pos = pOld.x1pos + pOld.L * cos(pOld.theta); pOld.y2pos = pOld.y1pos + pOld.L * sin(pOld.theta);
	pOld.mu += randomMu();
}

double inRange(double d0)
{
	double rp = d0;
	if (rp > 1.0) //check if relative point is in range [0,1]
		rp = 1.0; //outer right point of line segment, p.x2pos
	else if (rp < 0.0)
		rp = 0.0; //outer left point of line segment,  p.x1pos
	return rp;
}

void forceOverlap(Particle &p1, Particle &p2, double D1x, double D1y, double D2x, double D2y, double R, double rp1, double rp2)
{
	double f_x = 0.0; double f_y = 0.0;

	f_x = ko * (D2x - D1x) * (1 - R);
	f_y = ko * (D2y - D1y) * (1 - R);

	p1.Fx1 += (1 - rp1) * f_x; p1.Fx2 += rp1 * f_x; //Force is distributed over ends of line segments (xy1pos, xy2pos)
	p1.Fy1 += (1 - rp1) * f_y; p1.Fy2 += rp1 * f_y;
	p2.Fx1 -= (1 - rp2) * f_x; p2.Fx2 -= rp2 * f_x; //Newtons second law f1 = -f2
	p2.Fy1 -= (1 - rp2) * f_y; p2.Fy2 -= rp2 * f_y;
}

void distForce(Particle &p1, Particle &p2, vector<double> &distVector)
{
	double d1 = 0.0, d2; //relative point [0,1] on line segments of particles to determine closest distance, d1 starts as 0 when lines are parallel

						 //calculate dot products, s# corresponds to line segment of particle #, r is vector between first point of s1, s2
	double s1s1 = p1.L*p1.L; double s2s2 = p2.L*p2.L;
	double s1s2 = (p1.x2pos - p1.x1pos) * (p2.x2pos - p2.x1pos) + (p1.y2pos - p1.y1pos) * (p2.y2pos - p2.y1pos);
	double s1r = (p1.x2pos - p1.x1pos) * (p1.x1pos - p2.x1pos) + (p1.y2pos - p1.y1pos) * (p1.y1pos - p2.y1pos);
	double s2r = (p1.x1pos - p2.x1pos) * (p2.x2pos - p2.x1pos) + (p1.y1pos - p2.y1pos) * (p2.y2pos - p2.y1pos);

	double deler = s1s1*s2s2 - s1s2*s1s2; //denominator to determine d1, d2. d1 is set 0 if lines are parallel
	if (deler != 0.0)
	{
		d1 = inRange((s1s2*s2r - s2s2*s1r) / deler); //d1 = (s1s2*s2r - s2s2*s1r) / deler, check if in range
	}

	d2 = (s1s2*d1 + s2r) / s2s2;
	if (d2 > 1.0) //recalculate d1 if d2 does not lie on second line segment (not in range [0,1]), fill in for d2.
	{
		d2 = 1.0; //outer right point of second line segment, p2.x2pos
		d1 = inRange((s1s2 - s1r) / s1s1); //d1 = (s1s2 - s1r) / s1s1, check if in range
	}
	else if (d2 < 0.0)
	{
		d2 = 0.0; //outer left point of second line segment,  p2.x1pos
		d1 = inRange(-s1r / s1s1); //d1 = -s1r / s1s1, check if in range
	}
	//determine points on line segments that give smallest distance between line segments
	double D1x = p1.x1pos + d1*(p1.x2pos - p1.x1pos); double D1y = p1.y1pos + d1*(p1.y2pos - p1.y1pos);
	double D2x = p2.x1pos + d2*(p2.x2pos - p2.x1pos); double D2y = p2.y1pos + d2*(p2.y2pos - p2.y1pos);
	double sqd = (D1x - D2x)*(D1x - D2x) + (D1y - D2y)*(D1y - D2y); //square of distance between segments

	distVector.push_back(sqd);

	double touch = (p1.D + p2.D) / 2; //distance for which particles can touch each other
	if (sqd < (touch * touch)) //check if the particles overlap
	{
		double R = 0.0;
		if (sqd != 0.0)
			R = touch / sqrt(sqd); //gives the relative difference in length
		forceOverlap(p1, p2, D1x, D1y, D2x, D2y, R, d1, d2); //calculate forces between n and nb
	}

	if (sqd < 2 * p1.D*p2.D) //Check if particles are neighbours, if so, append to neighbour vector
	{
		p1.neighbours.push_back(p2.ID); //only this line is needed for ram, otherwise ram is performed 3 times more
		p2.neighbours.push_back(p1.ID); //use only for oop
	}
}

int main()
{
	string dataNumber = "0123456789abcdefghij"; //able to create up to 20 folders during one simulation
	for (int loop = 0;loop < 5;loop++) //multiple loops to gather data, set number after loop < ..., to set number of data folders
	{
		cout << "Type of simulation: loop nr = " << loop << ", nop = " << Nmax << ", growthrate = " << growthRate << ", DF = oop" << endl;
		int nop = 1; //number of particles is one at start
		int ts = 0; //timestep
		int save = 1000; //determines how many times data is saved, for nop = 3000 time runs till 650000.

		string nameData(1, dataNumber[loop]); //setup printing for visualisation and data writing
		Print print;
		print.init(nameData, nop);

		vector<Particle> p(nop, Particle(0.0, 0.0, 0.0, startLength, diameter, growthRate)); //initialise one particle

		while (nop < Nmax + 1) //timestep
		{
			//int nodef = 0; //number of defects

			for (int g = 0;g < nop;g++) //first loop to enable growth
			{
				p[g].grow();
				if (p[g].L > p[g].Lmax) //division
				{
					p[g].L = (p[g].x2pos - p[g].x1pos) / cos(p[g].theta); //Check if legal!! Set new length for dividing
					p.push_back(Particle(p[g].x1pos, p[g].y1pos, p[g].theta, p[g].L, p[g].D, p[g].mu)); //new particle is made with same properties as mother
					divide(p[g], p[p.size() - 1]); //function that sets new properties to daughter particles
					p[p.size() - 1].ID = p.size() - 1; //set ID of new particle
				}
			}

			nop = p.size(); //determine new number of particles after growth
			if (ts % save == 0)
				print.N = nop;

			vector<double> distance; //contains all squared distances between particles

			for (int m = 0;m < nop;m++) //second loop for movement of particles
			{
				if (ts % (save * 4) == 0) //print data at chosen timestep for vid
				{
					//p[m].colour = (int)((p[m].theta * 360 / (2 * pi))*2 + 240)%360; //Use for orientation to determine colour in vid
					print.print_data_js(p[m].x1pos, p[m].y1pos, p[m].x2pos, p[m].y2pos, p[m].D / 2, p[m].colour, p[m].ID);
				}
				p[m].colour = 240; //before finding defects set colour to initial value

				for (int nb = 0;nb < nop;nb++)
				{
					if (nb > m)
					{
						distForce(p[m], p[nb], distance); //calculate distances and overlapping forces between neighbours
					}
				}
				p[m].forceInternal(); //calculate internal spring forces
				p[m].move(); //move the particles according to the net force
				p[m].clear(); //clear values for velocities and forces for next step
			}

			for (int d = 0;d < nop;d++) //third loop to find defects
			{
				int amonb = p[d].neighbours.size(); //amount of neighbours
													//double sumcosThetapn = 0.0; //set sum of angles between particle and its neighbours for sop
													//double sumsinThetapn = 0.0;
				double sumsop = 0.0; //sum of sop over triangles
				int notr = 0; //number of triangles

				if (amonb > 0)
				{
					for (int a = 0;a < amonb;a++) //loop over neighbours
					{
						//use to determine orientational order parameter or sop
						int nb1 = p[d].neighbours[a]; //ID of neighbour 1
													  //double thetapn = sqrt((p[d].theta - p[nb1].theta) * (p[d].theta - p[nb1].theta)); //absolute angle between particle and neighbour
													  //if (thetapn > (pi / 2))
													  //	thetapn = pi - thetapn;
													  //double costhetapn = cos(2 * p[nb1].theta); //use for sop with all neighbours
													  //double sinthetapn = sin(2 * p[nb1].theta);
													  //sumcosThetapn += costhetapn;
													  //sumsinThetapn += sinthetapn;

						for (int b = 0;b < amonb;b++) //use for finding defects with rotational angle method
						{
							if (b > a)
							{
								int nb2 = p[d].neighbours[b]; //ID of neighbour 2
								int posDist = ((nop - 1)*nop / 2) - ((nop - nb1 - 2)*(nop - nb1 - 1) / 2) - (nop - nb2 - 1) - 1; //find position in distance vector of distance between neighbours 1 and 2
								if (distance[posDist] < 2 * p[d].D*p[d].D) //check if neighbours are also neighbours
								{
									double theta12 = (p[d].theta - p[nb1].theta) * (p[d].theta - p[nb1].theta);
									double theta23 = (p[nb1].theta - p[nb2].theta) * (p[nb1].theta - p[nb2].theta);
									double theta31 = (p[nb2].theta - p[d].theta) * (p[nb2].theta - p[d].theta);
									double ts = (pi / 2) * (pi / 2);
									//if only one angle of the three is larger than pi/2 than there is a defect
									if ((theta12 < ts && theta23 < ts && theta31 > ts) || (theta23 < ts && theta31 < ts && theta12 > ts) || (theta12 < ts && theta31 < ts && theta23 > ts))
									{
										p[d].ram = 1; //defect
													  //nodef++;
									}
									notr++;
									double sumcos = cos(2 * p[d].theta) + cos(2 * p[nb1].theta) + cos(2 * p[nb2].theta); //use for sop triangle
									double sumsin = sin(2 * p[d].theta) + sin(2 * p[nb1].theta) + sin(2 * p[nb2].theta);
									sumsop += 1 - 2 * sqrt((sumcos / 6) * (sumcos / 6) + (sumsin / 6) * (sumsin / 6));
								}
							}
						}
					}
					//p[d].oop = sumThetapn / (amonb * (pi / 2)); //gives orientational order parameter between 0 and 1
					//p[d].oop = 1 - 2 * sqrt((sumcosThetapn / (2 * amonb)) * (sumcosThetapn / (2 * amonb)) + (sumsinThetapn / (2 * amonb)) * (sumsinThetapn / (2 * amonb))); //gives sop from 0 to 1 for all nb
					if (notr > 0)
						p[d].oop = sumsop / notr; //gives sop from 0 to 1 for triangle
				}
				if (p[d].ram == 1)
					p[d].colour = 120;
				else
					p[d].colour = (int)(p[d].oop * 160 + 200);

				p[d].neighbours.clear(); //clear neighbours vector
				p[d].ram = 0;

			}
			distance.clear(); //clear distance vector

			if (ts % 500 == 0)
				cout << "At time " << ts << " nop is " << nop  << endl;

			ts++;
		}

		cout << "nop is " << nop << endl;
	}
	system("pause");

	return 0;
}