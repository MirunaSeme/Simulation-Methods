
/* Simulation Course

Assignment 1
Molecular Dynamics Simulation / Brownian Dynamics Simulation

Simulate the motion of the particles using brownian dynamics
Pedestrian crossing problem
No optimization

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define CONSTANT_DENSITY 1.00

struct particle_struct
{
	double x, y;        //x,y coordinate of the particles
	double fx, fy;      //fx,fy forces acting on the particle
	int color;          //this is to distinguish the particles
	int ID;             //ID of a particle
} *particles;

double SX, SY;          //system size x,y direction
double SX2, SY2;        //half of the system size x,y direction
int N;                  //number of particles
int T;					//number of time steps
double dt;              //length of a single time step
int t;                  //time - time steps so far
FILE *moviefile;        //file to store the coordinates of the particles
FILE *statistics_file;

//these are for time keeping purposes to count how many seconds the simulation ran
time_t 	time_start;
time_t 	time_end;

void program_timing_begin()
{
	time(&time_start);
}

void program_timing_end()
{
	double time_difference;
	struct tm *timeinfo;

	time(&time_end);

	time_difference = difftime(time_end, time_start);

	timeinfo = localtime(&time_start);
	printf("Program started at: %s", asctime(timeinfo));

	timeinfo = localtime(&time_end);
	printf("Program ended at: %s", asctime(timeinfo));

	printf("Program running time = %lf seconds\n", time_difference);
	printf("%d %lf\n", N, time_difference);
}

void initialize_particles()
{
	int i, j, ii, overlap;
	double dx, dy, dr, dr2;
	double tempx, tempy;

	SX = sqrt(N * (1 / CONSTANT_DENSITY));
	SY = sqrt(N * (1 / CONSTANT_DENSITY));
	SX2 = SX / 2.0;
	SY2 = SY / 2.0;

	particles = (struct particle_struct *) malloc(N * sizeof(struct particle_struct));

	dt = 0.002;
	ii = 0;

	for (i = 0; i<N; i++)
	{
		particles[i].ID = i;
		do
		{
			overlap = 0;
			tempx = SX * rand() / (RAND_MAX + 1.0);
			tempy = SY * rand() / (RAND_MAX + 1.0);

			for (j = 0; j<i; j++)
			{
				dx = tempx - particles[j].x;
				dy = tempy - particles[j].y;

				//PBC fold back
				if (dx >= SX2) dx -= SX;
				if (dx<-SX2) dx += SX;
				if (dy >= SY2) dy -= SY;
				if (dy<-SY2) dy += SY;

				dr2 = dx*dx + dy*dy;
				dr = sqrt(dr2);
				if (dr<0.2) overlap = 1;
			}
		} while (overlap == 1);

		//a better version would check how many attempts were made
		//after too many attempts quit gracefully (and not get stuck)

		particles[i].x = tempx;
		particles[i].y = tempy;

		//solve the problem: two particles should never be on top of each other!!!
		// 0.2 safe distance, I know the force at 0.2, it's not that big (around 100.0)
		particles[ii].fx = 0.0;
		particles[ii].fy = 0.0;

		if (rand() / (RAND_MAX + 1.0) < 0.5)
			particles[i].color = 0;
		else
			particles[i].color = 1;
		ii++;
	}
}

void calculate_external_forces()
{
	int i;

	for (i = 0; i<N; i++)
	{
		if (particles[i].color == 0)    particles[i].fx += 0.5;
		if (particles[i].color == 1)    particles[i].fx -= 0.5;
	}
}

void calculate_pairwise_forces()
{
	int i, j;
	double dx, dy;
	double dr2, dr;
	double f, fx, fy;

	for (i = 0; i<N - 1; i++)
		for (j = i + 1; j<N; j++)
		{
			dx = particles[i].x - particles[j].x;
			dy = particles[i].y - particles[j].y;

			//PBC check
			//maybe the neighbor cell copy of j is closer
			if (dx>SX2) dx -= SX;
			if (dx<-SX2) dx += SX;
			if (dy>SY2) dy -= SY;
			if (dy<-SY2) dy += SY;

			dr2 = dx*dx + dy*dy;

			//we are calculating the forces directly
			dr = sqrt(dr2);

			if (dr<0.2) { f = 100.0; printf("WARNING!!!\n"); }
			//nice way to do this: give a warning or exit if this happens
			else
				//check if dr>4.0 I can cut off the force
			{
				f = 1 / dr2 * exp(-0.25*dr);
			}

			//project it to the axes get the fx, fy components
			fx = f*dx / dr;
			fy = f*dy / dr;


			particles[i].fx += fx;
			particles[i].fy += fy;

			particles[j].fx -= fx;
			particles[j].fy -= fy;
		}
}

void move_particles()
{
	int i;
	double deltax, deltay;

	for (i = 0; i<N; i++)
	{
		//brownian dynamics
		//the particle is in a highly viscous environment
		deltax = particles[i].fx * dt;
		deltay = particles[i].fy * dt;

		particles[i].x += deltax;
		particles[i].y += deltay;

		//PBC check - check if they left the box
		//Box: 0,0 to SX, SY
		if (particles[i].x > SX) particles[i].x -= SX;
		if (particles[i].y > SY) particles[i].y -= SY;
		if (particles[i].x < 0) particles[i].x += SX;
		if (particles[i].y < 0) particles[i].y += SY;

		particles[i].fx = 0.0;
		particles[i].fy = 0.0;
	}
}

//this is for plot(linux plotter)
void write_cmovie()
{
	int i;
	float floatholder;
	int intholder;

	intholder = N;
	fwrite(&intholder, sizeof(int), 1, moviefile);

	intholder = t;
	fwrite(&intholder, sizeof(int), 1, moviefile);

	for (i = 0; i<N; i++)
	{
		intholder = particles[i].color + 2;
		fwrite(&intholder, sizeof(int), 1, moviefile);
		intholder = i;//ID
		fwrite(&intholder, sizeof(int), 1, moviefile);
		floatholder = (float)particles[i].x;
		fwrite(&floatholder, sizeof(float), 1, moviefile);
		floatholder = (float)particles[i].y;
		fwrite(&floatholder, sizeof(float), 1, moviefile);
		floatholder = 1.0;//cum_disp, cmovie format
		fwrite(&floatholder, sizeof(float), 1, moviefile);
	}

}

void write_statistics()
{
	int i;
	double avg_vx;

	avg_vx = 0.0;
	for (i = 0; i<N; i++)
	{
		avg_vx += particles[i].fx;
	}

	avg_vx = avg_vx / (double)N;

	fprintf(statistics_file, "%d %lf\n", t, avg_vx);

}

/*
 * arg1 - N - number of particles
 * arg2 - T - simulation time
 * arg3 - verbose - whether to print (1) progress or not (0).
 */
int main(int argc, char** argv)
{
	N = atoi(argv[1]);
    T = atoi(argv[2]);
    int verbose = atoi(argv[3]);

	printf("\n Pedestrian crossing simulation. No optimization. N=%d. T=%d.\n", N, T);
	program_timing_begin();

	initialize_particles();

    char movie_filename[100];
    char stats_filename[100];
    sprintf(movie_filename, "nop_%d_results.mvi", N);
    sprintf(stats_filename, "nop_%d_stats.txt", N);

	moviefile = fopen(movie_filename, "w");
	statistics_file = fopen(stats_filename, "wt");

	for (t = 0; t<T; t++)
	{
		calculate_pairwise_forces();
		//calculate_thermal_force();

		calculate_external_forces();

		//right now I have all the information
		//time to calculate some statistics
		write_statistics();

		move_particles();

		if (t % 100 == 0)
			write_cmovie();

		if (verbose && (t % 500 == 0))
		{
			printf("time = %d\n", t);
			fflush(stdout);
		}
	}

	fclose(moviefile);
	fclose(statistics_file);
	program_timing_end();
	return 0;
}