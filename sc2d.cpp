// -----------------------------------------------------------------
//  Simulation of Single Component Multiphase flow in 2D
//
//  Shan & Chen Model
//
//  Periodic boundary conditions
//
//  Written by: Abhijit Joshi
// -----------------------------------------------------------------

// OpenGL specific headers

#include <GLFW/glfw3.h>

// the usual gang of C++ headers

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC

// calculate pixel colors for the current graphics window, defined by the
// minimum and maximum X and Y coordinates

void showGraphics(const int NX, const int NY, int WIDTH, int HEIGHT, double xmin, double xmax, double ymin, double ymax, const double *rho)
{
    //--------------------------------
    //  OpenGL initialization stuff 
    //--------------------------------

    // select background color to be white
    // R = 1, G = 1, B = 1, alpha = 0
    glClearColor (1.0, 1.0, 1.0, 0.0);

    // initialize viewing values
    glMatrixMode(GL_PROJECTION);

    // replace current matrix with the identity matrix
    glLoadIdentity();

    // set clipping planes in the X-Y-Z coordinate system
    glOrtho(xmin,xmax,ymin,ymax, -1.0, 1.0);

    // clear all pixels
    glClear (GL_COLOR_BUFFER_BIT);

    // calculate pixel size (rectangle to be rendered)
    float dx = (xmax - xmin)/WIDTH;
    float dy = (ymax - ymin)/HEIGHT;

    // buffer used to store what we want to plot
    // 2D array size is identical to the window size in pixels
    float *scalar = new float[WIDTH*HEIGHT];

    // scale factors
    float min_rho = 0.1;   
    float max_rho = 1.7; 

    // loop to fill the buffer that OpenGL will render
    // and assign an appropriate color to that pixel
    for(int i = 0; i < WIDTH-1; i++) 
    {
        for(int j = 0; j < HEIGHT-1; j++) 
        {
            // map pixel coordinate (i,j) to LBM lattice coordinates (x,y)
            int xin = i*NX/WIDTH;
            int yin = j*NY/HEIGHT;

            // get locations of 4 data points inside which this pixel lies
            int idx00 = (xin  )*NY+(yin  );   // point (0,0)
            int idx10 = (xin+1)*NY+(yin  );   // point (1,0)
            int idx01 = (xin  )*NY+(yin+1);   // point (0,1)
            int idx11 = (xin+1)*NY+(yin+1);   // point (1,1)

            // calculate the normalized coordinates of the pixel
            float xfl = (float)i * (float)NX / (float) WIDTH;
            float yfl = (float)j * (float)NY / (float) HEIGHT;
            float x = xfl - (float)xin;
            float y = yfl - (float)yin;

            // bilinear interpolation
            float rho_interp = rho[idx00]*(1.0 - x)*(1.0 - y) + rho[idx10] * x * (1.0 - y) + rho[idx01] * (1.0 - x) * y + rho[idx11] * x * y;

            // this is the value we want to plot at this pixel (should be in the range [0-1])
            scalar[i*WIDTH + j] = (rho_interp - min_rho) / (max_rho - min_rho);                         // normalized density

            float x_actual = xmin + i*dx;   // actual x coordinate
            float y_actual = ymin + j*dy;   // actual y coordinate
            float VAL = scalar[i*WIDTH + j];

    //      printf("VAL = %f  ",VAL);
            float R, G, B;

            if(VAL<=0.5)
            {
                // yellow to blue transition
                R = 2*VAL;
                G = 2*VAL;
                B = 1 - 2*VAL;
            }
            else
            {
                // red to yellow transition
                R = 1;
                G = 2 - 2*VAL;
                B = 0;
            }

            // rendering the pixel with the appropriate color
            glColor3f(R,G,B);
            glRectf (x_actual,y_actual,x_actual+dx,y_actual+dy);
        }
    }

    // free memory
    delete[] scalar;
}

double psi(double x)
{
    const double E = 2.71828;
    const double rho0 = 1.0;
    return rho0 * (1 - pow(E, -x/rho0));
}

void initialize(const int NX, const int NY, const double rhoAvg,
                double* ex, double* ey, double* wt,
                double* rho, double* u, double* v,
                double* f, double* f_new, double* f_eq)
{
    // initialize density and velocity
    // initialize random seed
    srand (time(NULL));

    double rhoVar = 0.01 * rhoAvg;
    for(int i = 0; i < NX; i++)
    {
        for(int j = 0; j < NY; j++)
        {
            int N = i*NY + j;
            rho[N] = rhoAvg - 0.5*rhoVar + rhoVar * rand()/RAND_MAX;
            u[N] = 0.0;
            v[N] = 0.0;
        }
    }

    // initialize distribution functions to their equilibrium value

    for(int i = 0; i < NX; i++)
    {
        for(int j = 0; j < NY; j++)
        {
            int N = i*NY + j;
            double udotu = u[N]*u[N] + v[N]*v[N];

            for(int id = 0; id < 9; id++)
            {
                int index_f = 9*N + id;
                double edotu = ex[id]*u[N] + ey[id]*v[N];
                f_eq[index_f] = wt[id] * rho[N] * (1 + 3*edotu + 4.5*edotu*edotu - 1.5*udotu);
                f[index_f] = f_eq[index_f];
                f_new[index_f] = f_eq[index_f];
            }
        }
    }
}

// streaming 
void streaming(const int NX, int double NY,
               double* ex, double* ey, double tau,
               double* f, double* f_new, double* f_eq)
{
    for(int i = 0; i < NX-1; i++)
    {                                                                                                 
        for(int j = 0; j < NY-1; j++)
        {                                                                                             
            int N = i*NY + j;                                                                         
            for(int id = 0; id < 9; id++)                                                             
            {                                                                                         
                int iflow = i + ex[id];                                                               
                int jflow = j + ey[id];                                                               
       
                // periodic B.C.
                if(iflow == -1) {iflow = NX-2;}                                                       
                if(iflow == NX-1) {iflow = 0;}                                                        
                if(jflow == -1) {jflow = NY-2;}                                                       
                if(jflow == NY-1) {jflow = 0;}                                                        
                                                                                                          
                int Nflow = iflow*NY + jflow;                                                         
           
                int f_index_beg = 9*N + id;                                                           
                int f_index_end = 9*Nflow + id;                                                       
        
                f_new[f_index_end] = f[f_index_beg]                                                   
                                   - (f[f_index_beg] - f_eq[f_index_beg]) / tau;                      
            }
        }
    }
}

void calc_dPdt(const int NX, const int NY, 
               double* ex, double* ey, double* G11,
               double* rho, double* dPdt_x, double* dPdt_y)
{       
    // interparticle forces
    for(int i = 0; i < NX-1; i++)                                                                     
    {
        for(int j = 0; j < NY-1; j++)                                                                 
        {                                                                                             
            int N = i*NY + j;
            double Gsumx = 0;
            double Gsumy = 0;                                                                            
            for(int id = 0; id < 9; id++)                                                             
            {                                                                                         
                int iflow = i + ex[id];                                                               
                int jflow = j + ey[id];                                                               
        
                // periodic B.C.
                if(iflow == -1) {iflow = NX-2;}                                                       
                if(iflow == NX-1) {iflow = 0;}                                                        
                if(jflow == -1) {jflow = NY-2;}                                                       
                if(jflow == NY-1) {jflow = 0;}                                                        
        
                int Nflow = iflow*NY + jflow;                                                         
        
                Gsumx += psi(rho[N]) * psi(rho[Nflow]) * G11[id] * ex[id];                            
                Gsumy += psi(rho[N]) * psi(rho[Nflow]) * G11[id] * ey[id];                            
            }       
            dPdt_x[N] = -Gsumx;                                                                       
            dPdt_y[N] = -Gsumy;                                                                       
        }
    }
}

void updateDensityAndVelocity(const int NX, const int NY,
                              double* ex, double* ey, double* wt, double tau,
                              double* rho, double* u, double* v,
                              double* dPdt_x, double* dPdt_y,
                              double* f)
{
        // update density and velocity
        for(int i = 0; i < NX-1; i++)
        {  
            for(int j = 0; j < NY-1; j++)
            {  
                int N = i*NY + j;
                double f_sum = 0;
                double fex_sum = 0;
                double fey_sum = 0;
                for(int id = 0; id < 9; id++)
                {  
                    f_sum += f[9*N + id];
                    fex_sum += f[9*N + id]*ex[id];
                    fey_sum += f[9*N + id]*ey[id];
                }
                rho[N] = f_sum;
                u[N] = fex_sum / rho[N] + tau * dPdt_x[N] / rho[N];
                v[N] = fey_sum / rho[N] + tau * dPdt_y[N] / rho[N];
            }
        }

        // periodic B.C. for rho
        for(int i = 0; i < NX-1; i++)
        {  
            int N_end = i*NY + (NY-1);
            int N_beg = i*NY + 0;
            rho[N_end] = rho[N_beg];
        }
        for(int j = 0; j < NY-1; j++)
        {  
            int N_end = (NX-1)*NY + j;
            int N_beg = j;
            rho[N_end] = rho[N_beg];
        }
}

// main function
int main(void)
{
    // LBM parameters

    const int NX = 12;         // number of lattice points
    const int NY = 12;         // number of lattice points
    const double GEE11 = -0.45; // Shan & Chen parameter (controls density ratio)
    const double tau = 1.0;     // relaxation time
    const double rhoAvg = 0.693;  // reference density value

    // D2Q9 directions

    double ex[] = {0,1,0,-1,0,1,-1,-1,1}; 
    double ey[] = {0,0,1,0,-1,1,1,-1,-1}; 
    double wt[] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
    double G11[] = {0,GEE11,GEE11,GEE11,GEE11,GEE11/2,GEE11/2,GEE11/2,GEE11/2};

    // define buffers
    double*    rho = new double[NX*NY];   // density
    double*      u = new double[NX*NY];   // velocity x-component
    double*      v = new double[NX*NY];   // velocity y-component
    double* dPdt_x = new double[NX*NY];   // change of momentum along x
    double* dPdt_y = new double[NX*NY];   // change of momentum along x
    double* f      = new double[NX*NY*9]; // PDF
    double* f_eq   = new double[NX*NY*9]; // PDF
    double* f_new  = new double[NX*NY*9]; // PDF

    //--------------------------------
    //   Create a WINDOW using GLFW
    //--------------------------------

    GLFWwindow *window;

    // initialize the library
    if(!glfwInit())
        return -1;

    // window size for displaying graphics
    int WIDTH  = 800;
    int HEIGHT = 800;

    // set the window's display mode
    window = glfwCreateWindow(WIDTH, HEIGHT, "2D SCMP Simulation", NULL, NULL);
    if(!window) 
    {
        glfwTerminate();
	return -1;
    }

    // make the windows context current
    glfwMakeContextCurrent(window);

    // initialize fields
    initialize(NX, NY, rhoAvg, &ex[0], &ey[0], &wt[0], rho, u, v, f, f_new, f_eq);

    // time integration
    int time=0;
    clock_t t0, tN;
    t0 = clock();

    //---------------------------------------
    // Loop until the user closes the window
    //---------------------------------------

    // specify min and max window coordinates
    double xmin = 0, xmax = NX, ymin = 0, ymax = NY;

    while(!glfwWindowShouldClose(window))
    {
        // increment lattice time
        time++;

        streaming(NX, NY, ex, ey, tau, f, f_new, f_eq);

        calc_dPdt(NX, NY, ex, ey, G11, rho, dPdt_x, dPdt_y);

        updateDensityAndVelocity(NX, NY, ex, ey, wt, tau, rho, u, v, dPdt_x, dPdt_y, f);

        // update equilibrium functions
        for(int i = 0; i < NX-1; i++)
        {
            for(int j = 0; j < NY-1; j++)
            {
                int N = i*NY + j;
                double udotu = u[N]*u[N] + v[N]*v[N];

                for(int id = 0; id < 9; id++)
                {
                    int index_f = 9*N + id;
                    double edotu = ex[id]*u[N] + ey[id]*v[N];
                    f_eq[index_f] = wt[id] * rho[N] * (1 + 3*edotu + 4.5*edotu*edotu - 1.5*udotu);
                }
            }
        }

        // transfer fnew back to f

        for(int f_index = 0; f_index < NX*NY*9; f_index++)
        {
            f[f_index] = f_new[f_index];
        }

        // on-the-fly OpenGL graphics
        if(time%10 == 0) 
        {
            showGraphics(NX, NY, WIDTH, HEIGHT, xmin, xmax, ymin, ymax, rho);

            // swap front and back buffers
            glfwSwapBuffers(window);

            // poll for and processs events
            glfwPollEvents();
        }

        // calculate and print the number of lattice time-steps per second
        tN = clock() - t0;
        std::cout << "Lattice time " << time 
                  << " clock ticks " << tN 
                  << " wall clock time " << tN/CLOCKS_PER_SEC 
                  << " lattice time steps per second = " << (float) CLOCKS_PER_SEC * time / (float) tN 
                  << std::endl;
    }

    // clean up

    delete[] rho;
    delete[] u;
    delete[] v;
    delete[] dPdt_x;
    delete[] dPdt_y;
    delete[] f;
    delete[] f_eq;
    delete[] f_new;

    // GLFW clean up
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
