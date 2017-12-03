
/* Simulation Course

 Assignment 1
 Molecular Dynamics Simulation / Brownian Dynamics Simulation

 Simulate the motion of the particles using brownian dynamics
 Pedestrian crossing problem

 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define CONSTANT_DENSITY 1.00

struct particle_struct {
    double drx_so_far, dry_so_far;
    double x, y;                         //x,y coordinate of the particles
    double fx, fy;                       //fx,fy forces acting on the particle
    int color;                          //this is to distinguish the particles
    int ID;                             //ID of a particle
} *particles;

typedef struct cell_t {
    int *ps;            // array of particle indices.
    int n;              // number of particle pointers in array.
} cell;

double SX, SY;          //system size x,y direction
double SX2, SY2;        //half of the system size x,y direction
int N;                  //number of particles
int T;                  //number of time steps
double dt;              //length of a single time step
int t;                  //time - time steps so far
FILE *moviefile;        //file to store the coordinates of the particles
FILE *statistics_file;

//Verlet list variables
int *vlist1 = NULL;
int *vlist2 = NULL;
int N_vlist;

//this flag will tell me whether
//I need to rebuild the Verlet list
int flag_to_rebuild_Verlet;

/*
What the Verlet list stores:

 vlist1[5] vlist2[5] - these are the i,j numbers of 2 particles that interact
 vlist1[17] vlist2[17] - these are the i,j numbers of a different two particles that interact

 this is how I read the i,j from the Verlet list
 i = vlist1[7]
 j = vlist2[7]

 particles[vlist1[7]].x
 particles[vlist2[7]].x
 particles[vlist1[7]].y
 particles[vlist2[7]].y
 */

//variables for tabulating the force
double *tabulated_f_per_r;
int N_tabulated;
double tabulalt_start, tabulalt_lepes;

// Verlet Grid variables
cell **grid;
const int cell_size = 7;
int grid_cells_per_axis;

//these are for time keeping purposes
//to count how many seconds the simulation ran
time_t time_start;
time_t time_end;

void initialize_grid() {
    grid = (cell **) malloc(grid_cells_per_axis * sizeof(cell *));
    for (int i = 0; i < grid_cells_per_axis; i++) {
        grid[i] = (cell *) malloc(grid_cells_per_axis * sizeof(cell));
    }

    for (int i = 0; i < grid_cells_per_axis; i++)
        for (int j = 0; j < grid_cells_per_axis; j++) {
            grid[i][j].ps = (int *) malloc(N * sizeof(int));
            grid[i][j].n = 0;
        }
}

void clear_grid() {
    for (int i = 0; i < grid_cells_per_axis; i++)
        for (int j = 0; j < grid_cells_per_axis; j++) {
            grid[i][j].n = 0;
        }
}

void compute_cell_position(int p, int *x, int *y) {
    *x = (int) (particles[p].x / cell_size);
    *y = (int) (particles[p].y / cell_size);
}

int is_cell_out_of_bound(int x, int y) {
    if ((x >= grid_cells_per_axis) ||
        (x < 0) ||
        (y >= grid_cells_per_axis) ||
        (y < 0)) {
        return 1;
    } else {
        return 0;
    }
}

void get_neighbour_cells(int x, int y, int *ret_cells_x, int *ret_cells_y, int *ret_n) {
    int x_offsets[5] = {0, 0, 1, 1, 1};
    int y_offsets[5] = {0, 1, -1, 0, 1};
    int n = 0;
    for (int i = 0; i < 5; i++) {
        int new_x = x + x_offsets[i];
        int new_y = y + y_offsets[i];
        if (is_cell_out_of_bound(new_x, new_y) == 0) {
            ret_cells_x[n] = new_x;
            ret_cells_y[n] = new_y;
            n += 1;
        }
    }
    *ret_n = n;
}

void add_to_grid(int p) {
    int x, y;
    compute_cell_position(p, &x, &y);

    int n = grid[x][y].n;
    grid[x][y].ps[n] = p;
    grid[x][y].n += 1;
}

void rebuild_grid() {
    clear_grid();

    for (int i = 0; i < N; i++) {
        add_to_grid(i);
    }
}

void program_timing_begin() {
    time(&time_start);
}

void program_timing_end() {
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

void tabulate_forces() {
    int i;
    double x_min, x_max;
    double x2, x;
    double f;

    x_min = 0.1;
    x_max = 6.0;

    N_tabulated = 50000;
    tabulated_f_per_r = (double *) malloc(N_tabulated * sizeof(double));
    for (i = 0; i < N_tabulated; i++) {
        x2 = i * (x_max * x_max - x_min * x_min) / (N_tabulated - 1.0) + x_min * x_min;
        x = sqrt(x2);
        f = 1 / x2 * exp(-0.25 * x);
        tabulated_f_per_r[i] = f / x;
        //printf("%d %lf %lf %lf %lf\n",i,x,x2,f,tabulated_f_per_r[i]);
    }

    tabulalt_start = x_min * x_min;
    tabulalt_lepes = (x_max * x_max - x_min * x_min) / (N_tabulated - 1.0);
    printf("Tabulated start = %lf, steps = %lf\n", tabulalt_start, tabulalt_lepes);
}

void rebuild_verlet_list_grid_optimization() {
    int i, j, ix;
    double dx, dy, dr2;

    N_vlist = 0;
    vlist1 = (int *) realloc(vlist1, N_vlist * sizeof(int));
    vlist2 = (int *) realloc(vlist2, N_vlist * sizeof(int));
    rebuild_grid();

    int neighbour_cells_x[5];
    int neighbour_cells_y[5];
    int neighbour_cell_count;
    for (i = 0; i < N; i++) {
        int cell_x, cell_y;
        compute_cell_position(i, &cell_x, &cell_y);
        get_neighbour_cells(cell_x, cell_y, neighbour_cells_x, neighbour_cells_y, &neighbour_cell_count);
        for (int cell_i = 0; cell_i < neighbour_cell_count; cell_i++) {
            cell current_cell = grid[neighbour_cells_x[cell_i]][neighbour_cells_y[cell_i]];
            for (ix = 0; ix < current_cell.n; ix++) {
                j = current_cell.ps[ix];
                dx = particles[i].x - particles[j].x;
                dy = particles[i].y - particles[j].y;

                //PBC check
                //(maybe the neighbor cell copy is closer)
                if (dx > SX2) dx -= SX;
                if (dx < -SX2) dx += SX;
                if (dy > SY2) dy -= SY;
                if (dy < -SY2) dy += SY;

                dr2 = dx * dx + dy * dy;
                if (dr2 <= 36.0) //instead of 4*4 I will take 6*6
                {
                    N_vlist += 2;
                    vlist1 = (int *) realloc(vlist1, N_vlist * sizeof(int));
                    vlist2 = (int *) realloc(vlist2, N_vlist * sizeof(int));
                    vlist1[N_vlist - 2] = i;
                    vlist1[N_vlist - 1] = j;
                    vlist2[N_vlist - 2] = j;
                    vlist2[N_vlist - 1] = i;
                }
            }
        }
    }

    //once I rebuilt the Verlet list,
    //I can start counting the distances again
    for (i = 0; i < N; i++) {
        particles[i].drx_so_far = 0.0;
        particles[i].dry_so_far = 0.0;
    }
    flag_to_rebuild_Verlet = 0;
}

void rebuild_verlet_list() {
    int i, j;
    double dx, dy, dr2;

    //printf("rebuilding Verlet\n");fflush(stdout);

    N_vlist = 0;
    vlist1 = (int *) realloc(vlist1, N_vlist * sizeof(int));
    vlist2 = (int *) realloc(vlist2, N_vlist * sizeof(int));

    for (i = 0; i < N; i++)
        for (j = i + 1; j < N; j++) {
            dx = particles[i].x - particles[j].x;
            dy = particles[i].y - particles[j].y;

            //PBC check
            //(maybe the neighbor cell copy is closer)
            if (dx > SX2) dx -= SX;
            if (dx < -SX2) dx += SX;
            if (dy > SY2) dy -= SY;
            if (dy < -SY2) dy += SY;

            dr2 = dx * dx + dy * dy;

            //dr=sqrt(dr2)
            //if dr<6.0

            if (dr2 <= 36.0) //instead of 4*4 I will take 6*6
            {
                N_vlist++;
                vlist1 = (int *) realloc(vlist1, N_vlist * sizeof(int));
                vlist2 = (int *) realloc(vlist2, N_vlist * sizeof(int));
                vlist1[N_vlist - 1] = i;
                vlist2[N_vlist - 1] = j;
            }
        }
    /*
     for(i=0;i<N;i++)
     printf("%d %d \n",vlist1[i],vlist2[i]);
     */

    //once I rebuilt the Verlet list,
    //I can start counting the distances again
    for (i = 0; i < N; i++) {
        particles[i].drx_so_far = 0.0;
        particles[i].dry_so_far = 0.0;
    }
    flag_to_rebuild_Verlet = 0;
    //printf("Verlet rebuilt at t=%d\n",t);
}

void initialize_particles() {
    int i, j, ii, overlap;
    double dx, dy, dr, dr2;
    double tempx, tempy;
    double multiplier;

    //multiplier = 512.0;

    SX = sqrt(N * (1 / CONSTANT_DENSITY));
    SY = sqrt(N * (1 / CONSTANT_DENSITY));
    SX2 = SX / 2.0;
    SY2 = SY / 2.0;

    particles = (struct particle_struct *) malloc(N * sizeof(struct particle_struct));

    dt = 0.002;
    ii = 0;

    for (i = 0; i < N; i++) {
        particles[i].ID = i;
        do {
            overlap = 0;
            tempx = SX * rand() / (RAND_MAX + 1.0);
            tempy = SY * rand() / (RAND_MAX + 1.0);

            for (j = 0; j < i; j++) {
                dx = tempx - particles[j].x;
                dy = tempy - particles[j].y;

                //PBC fold back
                if (dx >= SX2) dx -= SX;
                if (dx < -SX2) dx += SX;
                if (dy >= SY2) dy -= SY;
                if (dy < -SY2) dy += SY;

                dr2 = dx * dx + dy * dy;
                dr = sqrt(dr2);
                if (dr < 0.2) overlap = 1;
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

void calculate_external_forces() {
    int i;

    for (i = 0; i < N; i++) {
        if (particles[i].color == 0) particles[i].fx += 0.5;
        if (particles[i].color == 1) particles[i].fx -= 0.5;
    }
}

void calculate_pairwise_forces() {
    int i, j;
    double dx, dy;
    double dr2, dr;
    double f, fx, fy;

    for (i = 0; i < N - 1; i++)
        for (j = i + 1; j < N; j++) {
            dx = particles[i].x - particles[j].x;
            dy = particles[i].y - particles[j].y;

            //PBC check
            //maybe the neighbor cell copy of j is closer
            if (dx > SX2) dx -= SX;
            if (dx < -SX2) dx += SX;
            if (dy > SY2) dy -= SY;
            if (dy < -SY2) dy += SY;

            dr2 = dx * dx + dy * dy;

            //we are calculating the forces directly
            dr = sqrt(dr2);

            if (dr < 0.2) {
                f = 100.0;
                printf("WARNING!!!\n");
            }
                //nice way to do this: give a warning or exit if this happens
            else
                //check if dr>4.0 I can cut off the force
            {
                f = 1 / dr2 * exp(-0.25 * dr);
            }

            //project it to the axes get the fx, fy components
            fx = f * dx / dr;
            fy = f * dy / dr;


            particles[i].fx += fx;
            particles[i].fy += fy;

            particles[j].fx -= fx;
            particles[j].fy -= fy;
        }
}

void calculate_pairwise_forces_with_verlet() {
    int i, j, ii;
    double dx, dy;
    double dr2, dr;
    double f, fx, fy;
    int tab_index;

    for (ii = 0; ii < N_vlist; ii++) {
        i = vlist1[ii];
        j = vlist2[ii];
        //printf("%d %d\n",i,j);
        dx = particles[i].x - particles[j].x;
        dy = particles[i].y - particles[j].y;

        //PBC check
        if (dx > SX2) dx -= SX;
        if (dx < -SX2) dx += SX;
        if (dy > SY2) dy -= SY;
        if (dy < -SY2) dy += SY;

        dr2 = dx * dx + dy * dy;

        //recall the tabulated value of the force


        tab_index = (int) floor((dr2 - tabulalt_start) / (tabulalt_lepes));
        if ((tab_index >= N_tabulated)) {
            //printf("tab_index = %d\n",tab_index);
            //tab_index = N_tabulated-1;
            //exit(1);
            fx = 0.0;
            fy = 0.0;
        } else {
            fx = tabulated_f_per_r[tab_index] * dx;  //f/dr is what I recalled
            fy = tabulated_f_per_r[tab_index] * dy;
        }



        //direct calculation of the force
        /*
        dr = sqrt(dr2); //this is very costly

        if (dr<0.1)
        {
            f = 97.53;
            printf("Warning! Particles %d and %d too close at time %d\n",i,j,t);
        }
        else
            //check if dr>4.0 I can cut off the force
        {
            f = 1/dr2 * exp(-0.25*dr);
        }

        fx = f*dx/dr;
        fy = f*dy/dr;
        //printf("%lf %lf hasonlitva %lf %lf\n",fx,fy,f*dx/dr,f*dy/dr);
        */

        particles[i].fx += fx;
        particles[i].fy += fy;
        particles[j].fx -= fx;
        particles[j].fy -= fy;
    }
}

void move_particles() {
    int i;
    double deltax, deltay;

    for (i = 0; i < N; i++) {
        //brownian dynamics
        //the particle is in a highly viscous environment
        deltax = particles[i].fx * dt;
        deltay = particles[i].fy * dt;

        particles[i].x += deltax;
        particles[i].y += deltay;

        particles[i].drx_so_far += deltax;
        particles[i].dry_so_far += deltay;


        if ((particles[i].drx_so_far * particles[i].drx_so_far +
             particles[i].dry_so_far * particles[i].dry_so_far) >= 4.0)
            flag_to_rebuild_Verlet = 1;

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
void write_cmovie() {
    int i;
    float floatholder;
    int intholder;

    intholder = N;
    fwrite(&intholder, sizeof(int), 1, moviefile);

    intholder = t;
    fwrite(&intholder, sizeof(int), 1, moviefile);

    for (i = 0; i < N; i++) {
        intholder = particles[i].color + 2;
        fwrite(&intholder, sizeof(int), 1, moviefile);
        intholder = i;//ID
        fwrite(&intholder, sizeof(int), 1, moviefile);
        floatholder = (float) particles[i].x;
        fwrite(&floatholder, sizeof(float), 1, moviefile);
        floatholder = (float) particles[i].y;
        fwrite(&floatholder, sizeof(float), 1, moviefile);
        floatholder = 1.0;//cum_disp, cmovie format
        fwrite(&floatholder, sizeof(float), 1, moviefile);
    }

}

void write_statistics() {
    int i;
    double avg_vx;

    avg_vx = 0.0;
    for (i = 0; i < N; i++) {
        avg_vx += particles[i].fx;
    }

    avg_vx = avg_vx / (double) N;

    fprintf(statistics_file, "%d %lf\n", t, avg_vx);

}

/*
 * arg1 - N - number of particles
 * arg2 - T - simulation time
 * arg3 - verbose - whether to print (1) progress or not (0).
 */
int main(int argc, char **argv) {
    N = atoi(argv[1]);
    T = atoi(argv[2]);
    int verbose = atoi(argv[3]);

    grid_cells_per_axis = (int) ceil(sqrt(N) / (double) cell_size);
    initialize_grid();

    printf("\nPedestrian crossing simulation. Verlet List and Grid Lookup optimizations. N=%d. T=%d.\n", N, T);
    program_timing_begin();

    tabulate_forces();

    initialize_particles();

    char movie_filename[100];
    char stats_filename[100];
    sprintf(movie_filename, "vlist_%d_results.mvi", N);
    sprintf(stats_filename, "vlist_%d_stats.txt", N);

    moviefile = fopen(movie_filename, "w");
    statistics_file = fopen(stats_filename, "wt");

    for (t = 0; t < T; t++) {
        calculate_pairwise_forces_with_verlet();

        //calculate_thermal_force();
        //calculate_pairwise_forces();

        calculate_external_forces();

        //right now I have all the information
        //time to calculate some statistics
        write_statistics();

        move_particles();

        if (flag_to_rebuild_Verlet)
            rebuild_verlet_list_grid_optimization();

        if (t % 100 == 0)
            write_cmovie();
        //write_movie_frame();

        if (verbose && (t % 500 == 0)) {
            printf("time = %d\n", t);
            fflush(stdout);
        }
    }

    fclose(moviefile);
    fclose(statistics_file);
    program_timing_end();
    return 0;
}