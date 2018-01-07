//  dla_fraktal

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_grid 2000
#define PI 3.141592653579892
#define STEPS 100
#define TIME 1000000
#define ARRAY_SIZE 1000

int grid[N_grid][N_grid];

double R_max;

int randomwalker_x, randomwalker_y;
double randomwalker_R;

double p_samples[9] = {0.001, 0.04, 0.01,
                       0.1, 0.3, 0.5,
                       0.7, 0.9, 0.9999};

int t;

FILE *moviefile;
FILE *results_file;


int N_particles;

void write_cmovie(void);

void init_grid()
{
    int i,j;

    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            grid[i][j] = 0;

    grid[N_grid/2][N_grid/2] = 1;

    R_max = 5.0;
    N_particles = 1;
}


void init_randomwalker()
{
    double theta;

    if (( (double)N_grid/2 + 3.0*R_max ) > (double) N_grid )
        {
            printf("Cannot create any more random walkers");
            exit(1);
        }

    randomwalker_R = (R_max+5.0);
    theta = 2*PI*rand()/(RAND_MAX+1.0);

    randomwalker_x = N_grid/2 + (int)floor(randomwalker_R * cos(theta));
    randomwalker_y = N_grid/2 + (int)floor(randomwalker_R * sin(theta));

}

void move_randomwalker()
{
    double r;
    int tempx, tempy;

    tempx = randomwalker_x;
    tempy = randomwalker_y;

    r = rand()/(RAND_MAX+1.0);

    if (r<0.25) tempx++;//right
    else if (r<0.5) tempy++;//up
    else if (r<0.75) tempx--;//left
    else tempy--;//down

    if (grid[tempx][tempy]==0)
        {
            randomwalker_x = tempx;
            randomwalker_y = tempy;
        }



}

void does_it_stick(double probability)
{
    double R_max_candidate2;
    double dx,dy;

    double r;

    r = rand()/(RAND_MAX+1.0);
    //printf("Random number =%f\n",r);


      if (r<probability){
      //if it has a neighbor grid site that is occupied
      //nearest neighbor (right, left,up,down)
      // grid[randomwalker_x+1][randomwalker_y+1]
      //would be next nearest neighbor
      if ((grid[randomwalker_x+1][randomwalker_y]==1)||
          (grid[randomwalker_x-1][randomwalker_y]==1)||
          (grid[randomwalker_x][randomwalker_y+1]==1)||
          (grid[randomwalker_x][randomwalker_y-1]==1))
      {

          //printf("sticking at t=%d\n",t);

          grid[randomwalker_x][randomwalker_y] = 1;
          N_particles++;
          write_cmovie();

          dx = (randomwalker_x-(N_grid/2));
          dy = (randomwalker_y-(N_grid/2));

          R_max_candidate2 = dx*dx + dy*dy;

          if (R_max_candidate2 > R_max*R_max)
              R_max = sqrt(R_max_candidate2);

          init_randomwalker();
      }
    }
}


int the_random_walker_is_lost()
{
    double dx,dy, dr2;

    dx  = randomwalker_x - (N_grid/2.0);
    dy  = randomwalker_y - (N_grid/2.0);

    dr2 = dx*dx + dy*dy;

    if (dr2>9.0*R_max*R_max) return 1;
    else return 0;
}

void write_cmovie()
{
    int i,j,k;
    float floatholder;
    int intholder;

    intholder = N_particles;
    fwrite(&intholder,sizeof(int),1,moviefile);

    intholder = t;
    fwrite(&intholder,sizeof(int),1,moviefile);

    k=0;
    for(i=0;i<N_grid;i++)
        for(j=0;j<N_grid;j++)
            if (grid[i][j])
            {
            intholder = 3; //color of spin
            //else intholder = 3 + (klaszterszam[i][j]%10);
            fwrite(&intholder,sizeof(int),1,moviefile);
            intholder = k++;            //spin ID
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

    printf("Diffusion Limited Aggregation\n");

    char filename[256];
    char results_filename[256];


    for(int i = 0; i < 9; i++){
      double probability = p_samples[i];
      double rmax_array[9][ARRAY_SIZE];

      snprintf(filename, sizeof filename, "dla%.3lf.mvi", probability);
      snprintf(results_filename, sizeof results_filename, "results_%.3lf.txt", probability);

      moviefile = fopen(filename,"wb");
      results_file = fopen(results_filename, "wb");

      init_grid();
      init_randomwalker();

      for(int j=0; j < STEPS; j++){
        for(t=0; t < TIME; t++)
            {
                //srand(1446742268);
                int seed = (int)time(NULL)+t*10;
                // printf"%d seed=%d\n",t,seed);
                srand(seed);
                if (t%100000 == 0) printf("=t%d\n",t);
                move_randomwalker();
                if (the_random_walker_is_lost()) init_randomwalker();
                does_it_stick(probability);
                if (t%1000 == 0){
                  rmax_array[i][t/1000] += R_max;
                }
            }
      }

        for(int k = 0; k < ARRAY_SIZE; k++){
          fprintf(results_file, "%d %lf\n", (k+1)*ARRAY_SIZE, rmax_array[i][k]/STEPS);
        }

      fclose(moviefile);
      fclose(results_file);
    }
    return 0;
}
