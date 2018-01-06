
//  main.c
//  ising_model


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define L_max 100

#define MOVIE 0
#define ECHO 0
#define LOG 0

#define T_samples 19

int L_spin;
int N_spin;
int spin[L_max][L_max];
double T;
long int N_steps_per_T;

long int step;
double M_avg,M2_avg,E_avg,E2_avg;
long int N_rej;

double T_values[T_samples] = {0.2, 0.5, 0.7, 1.0, 1.2, 1.5, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.5, 2.8, 3.0, 3.3, 3.5, 3.8};

FILE *outfile;
FILE *moviefile;
FILE *averages;

void initialize_statistics();
void initialize_spins();
void write_cmovie();
void write_out_results();

void initialize_statistics()
{
    M_avg = 0.0;
    M2_avg = 0.0;
    E_avg = 0.0;
    E2_avg = 0.0;
    N_rej = 0;
}

void initialize_spins()
{
    int i,j;
    int pos,neg;
    
    pos = 0;
    neg = 0;
    
    for(i=0;i<L_spin;i++)
        for(j=0;j<L_spin;j++)
            {
            
            if ( (rand()/(RAND_MAX+1.0)) < 0.5)
                {
                    spin[i][j] = 1;
                    pos++;
                }
            else
                {
                    spin[i][j] = -1;
                    neg++;
                }
            //spin[i][j]=-1;
            }
    printf("Initialization complete\n");
    printf("pos = %d neg = %d\n",pos,neg);
   
}

int Energy_contribution_of_a_spin(int fi,int fj)
{
    int ni,nj; //neighbor, we have 4 neighbors (von Neumann neighborhood)
    int i;
    int delta = 0;
    
    for(i=0;i<4;i++)
        {

            switch(i)
                {
                    case 0: {ni = fi - 1; nj = fj; break;}
                    case 1: {ni = fi + 1; nj = fj; break;}
                    case 2: {ni = fi;     nj = fj + 1; break;}
                    case 3: {ni = fi;     nj = fj - 1; break; }
                    default: {
                                ni=fi;nj=fj;
                                printf("You should never reach this line\n");
                             }
                }
            
            //PBC fold-back
            if (ni<0) ni = ni + L_spin;
            if (ni>=L_spin) ni = ni - L_spin;
            if (nj<0) nj = nj + L_spin;
            if (nj>=L_spin) nj = nj - L_spin;
            
            delta = delta - spin[fi][fj]*spin[ni][nj];
            
        }
    
    return delta;
}

int delta_E(int fi,int fj)
{
    int dE;
    
    dE = - 2 * Energy_contribution_of_a_spin(fi,fj);

    return dE;
}

void flip_spin(flipping_i,flipping_j)
{
    spin[flipping_i][flipping_j] = - spin[flipping_i][flipping_j];
}

void calculate_statistics()
{

    double M,M2,E,E2;
    int i,j;

    M = 0.0;
    E = 0.0;
    
    for(i=0;i<L_spin;i++)
        for(j=0;j<L_spin;j++)
            {
                M = M + (double)spin[i][j];
                E = E + (double)Energy_contribution_of_a_spin(i,j);
            }
    
    M2 = M * M;
    E2 = E * E;
    
    //write out an echo for every 1000 steps
    //so the user would know the program is still running

    if (step%1000==0)
        {
        	if (ECHO)
        	{
        		printf("t= %d M = %lf\n",step,M/N_spin);
        	}
            if (LOG){
            	fprintf(outfile,"%d %lf\n",step,M/(double)N_spin);
        	}
        }

    //M_avg = M_avg * ((double)(step)/(double)(step+1.0)) + M/ (double)(step+1.0);
    M_avg += abs(M)/(double)N_steps_per_T;
    
    //M2 needs fixing because it needs to be divided by N^2
    //so does E2
    //you can just add it all together and then divide by N for M,E and by N^2 for M2,E2
    M2_avg = abs(M2_avg) * ((double)(step)/(double)(step+1.0)) + M2/ (double)(step+1.0);
    E_avg = abs(E_avg) * ((double)(step)/(double)(step+1.0)) + E/ (double)(step+1.0);
    E2_avg = abs(E2_avg) * ((double)(step)/(double)(step+1.0)) + E2/ (double)(step+1.0);
    
}

void running_for_fixed_T()
{
    int flipping_i,flipping_j;
    int dE;
    double p,r;
    
    for(step=0;step<N_steps_per_T;step++)
        {
            //pick the spin at random
            flipping_i = (int)floor((double)L_spin * rand()/(RAND_MAX+1.0));
            flipping_j = (int)floor((double)L_spin * rand()/(RAND_MAX+1.0));
            
            //calculate the delta E that this spin flip will cause
            dE = delta_E(flipping_i,flipping_j);
            
            //Metropolis MC rejection
            if (dE<=0) {
            
                        //if the energy is lowered then ALWAYS make the step
                        flip_spin(flipping_i,flipping_j);
                        //printf("sure\n");
                        }
            
            else        {
                        //the energy is actually raised by this flip
            
                        p = exp(-(double)dE/(double)T);
                        r = rand()/(RAND_MAX+1.0);
                
                        if (r<p) 
							flip_spin(flipping_i,flipping_j);
				 			//this actually does the flip
                        else 
                        	N_rej++;
                    }
            
            calculate_statistics();
            if (step%10000==0 && MOVIE) write_cmovie();
            
        }
    
}

void write_out_results()
{
    M_avg = M_avg / N_spin;
    E_avg = E_avg / N_spin;
    M2_avg = M2_avg / (N_spin*N_spin);
    E2_avg = E2_avg / (N_spin*N_spin);

    printf("T = %lf M = %lf Rejections = %ld\n",T,M_avg,N_rej);
    fprintf(averages,"%lf %lf %lf %lf %lf\n", T, M_avg, M2_avg, E_avg, E2_avg);


}


void write_cmovie()
{
    int i,j;
    float floatholder;
    int intholder;
    
    intholder = N_spin;
    fwrite(&intholder,sizeof(int),1,moviefile);
    
    intholder = (int)step;
    fwrite(&intholder,sizeof(int),1,moviefile);
  
    for(i=0;i<L_spin;i++)
        for(j=0;j<L_spin;j++)
        {
            if (spin[i][j]==-1) intholder = 2; //color of spin
            else intholder = 3;
            fwrite(&intholder,sizeof(int),1,moviefile);
            intholder = i*L_spin+j;            //spin ID
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

    printf("Ising model simulation\n");
    char out[256];
    char movie[256];
    char avr[256];

    // set the size of the system
    L_spin = 50;
    N_spin = L_spin * L_spin;
    
    N_steps_per_T = 1000000;

    snprintf(avr, sizeof avr, "%d%s%s", L_spin, "averages", ".dat");
	averages = fopen(avr,"wt");	


    srand((int)time(NULL));

    for (int i=0; i<T_samples; i++){

    	if (LOG)
	    {
	    	snprintf(out, sizeof out, "./Data/%d%s%d%s", L_spin, "ising_data_", i, ".dat");
	    	outfile = fopen(out,"wt");
    	}
    	if (MOVIE)
    	{
	    	snprintf(movie, sizeof movie, "./Data/%d%s%d%s", L_spin, "ising_movie_", i, ".dat");
	    	moviefile = fopen(movie,"wb");
    	}

    	T=T_values[i];

    	printf("\n RUNNING FOR T = %lf\n",T);

    	initialize_spins();
        initialize_statistics();
        running_for_fixed_T();
        write_out_results();

        if (LOG)
	    {
	    	fclose(outfile);
	    }
	    if(MOVIE)
	    {
	    	fclose(moviefile);
	    }
    }

    fclose(averages);

    return 0;
 
}
