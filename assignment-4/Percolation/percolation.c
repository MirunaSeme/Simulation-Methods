//
//percolation simulation


#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define N_grid 150
#define STEPS 1000
#define MOVIE 0
//#define p_steps 11
#define P_STEPS 30

int grid[N_grid][N_grid];
int cluster_number[N_grid][N_grid];

int actual_cluster;

double p_samples[P_STEPS] = {0.10, 0.20, 0.30, 0.35, 0.40, 0.42, 0.45,
                             0.46, 0.47, 0.48, 0.49, 0.50, 
                             0.51, 0.52, 0.53, 0.54, 0.55, 
                             0.56, 0.57, 0.58, 0.59, 0.60,
                             0.61, 0.62, 0.63, 0.64, 0.65, 
                             0.70, 0.8, 0.9};
//double p_samples[P_STEPS] = {0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9};

double  p;

FILE *moviefile;
FILE *resultsfile;

int t;

void initialize_system()
{
    int i,j;
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            {
            grid[i][j] = 0;
            cluster_number[i][j] = -1;
            }
}

void fill_system_with_probability(double p)
{
    int i,j;
    double r;
    int N_particles;
    
    N_particles = 0;
    
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            {
                cluster_number[i][j] = -1;
                
                r = rand()/(RAND_MAX+1.0);
                //if r<p, this site will be filled
                if (r<p)
                    {
                        grid[i][j] = 1;
                        N_particles++;
                    }
            }
    printf("Filled up the system. Ended up with %d particles \n",N_particles);
}

void recursive_clusternumber(int i,int j)
{
    
    cluster_number[i][j] = actual_cluster;
    
    if ((i+1<N_grid))
        //find a position to the right of the actual position
        if ((cluster_number[i+1][j]==-1)&&(grid[i+1][j]==1))
            recursive_clusternumber(i+1,j);
        //find a position to the left
    if ((i-1>=0))
        if ((cluster_number[i-1][j]==-1)&&(grid[i-1][j]==1))
            recursive_clusternumber(i-1,j);
        //up
    if ((j+1<N_grid))
        if ((cluster_number[i][j+1]==-1)&&(grid[i][j+1]==1))
            recursive_clusternumber(i,j+1);
        //down
    if ((j-1>=0))
        if ((cluster_number[i][j-1]==-1)&&(grid[i][j-1]==1))
            recursive_clusternumber(i,j-1);
}

void clusterize_system()
{
    int i,j;
    
    actual_cluster = 1;
    
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            {
                //filled position, it's not in any cluser
                if ((grid[i][j]==1)&&(cluster_number[i][j]==-1))
                    {
                        //start the recursive algorithm
                        //to find everybody connected to this point
                        recursive_clusternumber(i,j);
                        //the next cluster is going to
                        //be a bigger number
                        actual_cluster++;
                    }
            }
    actual_cluster--;
    
}

int is_spanning_cluster()
{
	
	int spanning_flag =0;

	for (int i=0; i<N_grid; i++)
	{
		for (int j=0; j<N_grid; j++){
			if (cluster_number[0][i] == cluster_number[N_grid-1][j] && cluster_number[0][i]!= -1) spanning_flag = 1;
		}
	}

	return spanning_flag;
}

int size_of_cluster(int cluster)
{
	int sum=0;
	for (int i=0; i < N_grid; i++)
	{
		for (int j=0; j<N_grid; j++)
		{
			if (cluster_number[i][j] == cluster) sum++;
		}
	}

	return sum;
}

int max_size(){
	int max = size_of_cluster(1);
	if (actual_cluster >2)
	{
		for (int i=2; i< actual_cluster; i++)
		{
			int size_i = size_of_cluster(i);
			if (size_i > max) max =  size_i;
		}
	}
	return max;
}

void write_cmovie()
{
    int i,j;
    float floatholder;
    int intholder;
    
    intholder = N_grid*N_grid;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    intholder = t;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
        {
            //color the particles based on which
            //cluster they belong to
            if (grid[i][j]==0) intholder = 2;
            else intholder = 3 + (cluster_number[i][j]%10);
            fwrite(&intholder,sizeof(int),1,moviefile);
            intholder = i*N_grid+j;            //spin ID
            fwrite(&intholder,sizeof(int),1,moviefile);
            floatholder = (float)i;
            fwrite(&floatholder,sizeof(float),1, moviefile);
            floatholder = (float)j;
            fwrite(&floatholder,sizeof(float),1,moviefile);
            floatholder = 1.0;//cum_disp, cmovie format
            fwrite(&floatholder,sizeof(float),1,moviefile);
        }
}


int main(int argc, const char * argv[]) {
   
    printf("Percolation calculation \n");

    printf("SIZE %ld", sizeof(int));


    char results[256];
    snprintf(results, sizeof results, "results_%d_%d.dat", N_grid, P_STEPS);
    resultsfile = fopen(results, "wt");

    for (int i_p=0; i_p< P_STEPS; i_p++)
    {

        //probability of occupying a given site: p
        p = p_samples[i_p];
        t = 0;
        
        char movie[256];

        if (movie)
        {
            snprintf(movie, sizeof movie, "./Data/movie_%d_%.2lf.dat", N_grid, p);
            moviefile = fopen(movie,"wb");
        }

        double total_spannings = 0.0;
        double max_avg = 0.0;
        
        for(t=0;t<STEPS;t++)
            {

                int seed = (int)time(NULL)+t*10;
                srand(seed);
                
                initialize_system();
                fill_system_with_probability(p);
                clusterize_system();

                total_spannings += is_spanning_cluster();
                if (movie)
                {
                    write_cmovie();
                }   

                int max = max_size(); 
                max_avg +=max;

                printf("Step = %d, Max = %d\n", t, max);
            }
        
        total_spannings = (double)total_spannings/STEPS; 
        max_avg = max_avg/(double)STEPS;

        printf("Probability of spanning = %lf \nAverage max cluster = %lf \n", total_spannings, max_avg);
        fprintf(resultsfile, "%lf %lf %lf\n", p, total_spannings, max_avg);

        fclose(moviefile);
    }

    fclose(resultsfile);
    return 0;
}
