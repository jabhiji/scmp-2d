//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
//    Simulation of Single Component Multiphase flow in 2D
//
//    - Shan & Chen Model for interactive forces 
//    - Kupershtokh & Medvedev Exact Difference Method
//
//    Periodic boundary conditions
//
//    Written by: Abhijit Joshi
//    ------------------------------------------------------------------

//    GLFW library header 

      #include <GLFW/glfw3.h>

//    C++ headers

      #include <iostream>     // cout()
      #include <cmath>        // pow()
      #include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC
      #include <iomanip>      // setw

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
//    calculate pixel colors for the current graphics window

      void showGraphics(const int NX, const int NY, 
                        const int WIDTH, const int HEIGHT, 
                        const double xmin, const double xmax, 
                        const double ymin, const double ymax, 
                        const double *rho)
      {
//      --------------------------------
//        OpenGL initialization stuff 
//      --------------------------------

//      select background color to be white
//      R = 1, G = 1, B = 1, alpha = 0

        glClearColor (1.0, 1.0, 1.0, 0.0);

//      initialize viewing values

        glMatrixMode(GL_PROJECTION);

//      replace current matrix with the identity matrix

        glLoadIdentity();

//      set clipping planes in the X-Y-Z coordinate system

        glOrtho(xmin,xmax,ymin,ymax, -1.0, 1.0);

//      clear all pixels

        glClear (GL_COLOR_BUFFER_BIT);

//      calculate pixel size (rectangle to be rendered)

        float dx = (xmax - xmin)/WIDTH;
        float dy = (ymax - ymin)/HEIGHT;

//      find min and max rho values (for color map)

        float min_rho = 10.0;
        float max_rho = 0.0;
        for(int k = 0; k < NX*NY; k++)
        {
          if(rho[k] > max_rho) max_rho = rho[k];
          if(rho[k] < min_rho) min_rho = rho[k];
        }

//      loop to fill the buffer that OpenGL will render
//      and assign an appropriate color to that pixel

        for(int i = 0; i < WIDTH; i++)
        {
          for(int j = 0; j < HEIGHT; j++)
          {
//          map pixel coordinate (i,j) to LBM lattice coordinates (x,y)

            int xin = i*(NX-1)/(WIDTH-1);
            int yin = j*(NY-1)/(HEIGHT-1);

//          get locations of 4 data points inside which this pixel lies

            int idx00 = (xin  )*NY+(yin  );   // point (0,0)
            int idx10 = (xin+1)*NY+(yin  );   // point (1,0)
            int idx01 = (xin  )*NY+(yin+1);   // point (0,1)
            int idx11 = (xin+1)*NY+(yin+1);   // point (1,1)

//          calculate the normalized coordinates of the pixel

            float xfl = (float)i * (float)(NX-1) / (float) (WIDTH-1);
            float yfl = (float)j * (float)(NY-1) / (float) (HEIGHT-1);
            float x = xfl - (float)xin;
            float y = yfl - (float)yin;

//          bilinear interpolation to get rho value for pixel (i,j)

            float rho_interp = rho[idx00] * (1.0 - x) * (1.0 - y) 
                             + rho[idx10] *     x     * (1.0 - y) 
                             + rho[idx01] * (1.0 - x) * y 
                             + rho[idx11] *     x     * y;

//          normalized rho (should be in the range [0-1])

            float rho_norm = (rho_interp - min_rho) 
                           / (max_rho    - min_rho);

//          get coordinates of bottom-left corner of this pixel

            float x_actual = xmin + i*dx;   // x coordinate
            float y_actual = ymin + j*dy;   // y coordinate

            float R, G, B;

            if(rho_norm<=0.5)
            {
//            yellow to blue transition

              R = 2*rho_norm;
              G = 2*rho_norm;
              B = 1 - 2*rho_norm;
            }
            else
            {
//            red to yellow transition

              R = 1;
              G = 2 - 2*rho_norm;
              B = 0;
            }

            // rendering the pixel with the appropriate color

            glColor3f(R,G,B);
            glRectf (x_actual,y_actual,x_actual+dx,y_actual+dy);
          }
        }
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
//    funtion to calculate effective density in the Shan & Chen model

      double psi(double x)
      {
        const double E = 2.71828;
        const double rho0 = 1.0;
        return rho0 * (1 - pow(E, -x/rho0));
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
      void initialize(const int NX, const int NY, const double rhoAvg,
                      double* ex, double* ey, double* wt,
                      double* rho, double* u, double* v,
                      double* f, double* f_new, double* f_eq)
      {
//      initialize random seed

        srand (time(NULL));

//      initialize density and velocity

        double rhoVar = 0.01 * rhoAvg;
        for(int i = 0; i < NX-1; i++)
        {
          for(int j = 0; j < NY-1; j++)
          {
            int N = i*NY + j;
            rho[N] = rhoAvg - 0.5*rhoVar + rhoVar * rand()/RAND_MAX;
            u[N] = 0.0;
            v[N] = 0.0;
          }
        }

//      initialize distribution functions to their equilibrium value

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
              f_eq[index_f] = wt[id] * rho[N]
                            * (1 + 3*edotu
                                 + 4.5*edotu*edotu - 1.5*udotu);
              f[index_f] = f_eq[index_f];
              f_new[index_f] = f_eq[index_f];
            }
          }
        }
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
      void streaming(const int NX, const int NY,
                     double* ex, double* ey, double tau,
                     double* f, double* f_new, double* f_eq, double* f_eq_tilda)
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
              if(jflow == -1) {jflow = NY-2;}
              if(iflow == NX-1) {iflow = 0;}
              if(jflow == NY-1) {jflow = 0;}

              int Nflow = iflow*NY + jflow;
              int f_index_beg = 9*N + id; 
              int f_index_end = 9*Nflow + id;
        
              f_new[f_index_end] = f[f_index_beg]
                                 - (f[f_index_beg] - f_eq[f_index_beg]) / tau
                                 + f_eq_tilda[f_index_beg] - f_eq[f_index_beg];
            }
          }
        }
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
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
              if(iflow == -1)   iflow = NX-2;
              if(jflow == -1)   jflow = NY-2;
              if(iflow == NX-1) iflow = 0;
              if(jflow == NY-1) jflow = 0;
            
              int Nflow = iflow*NY + jflow;
            
              Gsumx += psi(rho[N]) * psi(rho[Nflow]) * G11[id] * ex[id];
              Gsumy += psi(rho[N]) * psi(rho[Nflow]) * G11[id] * ey[id];
            }
            dPdt_x[N] = -Gsumx;
            dPdt_y[N] = -Gsumy;
          }
        }
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
      void updateDensityAndVelocity(const int NX, const int NY,
                                    double* ex, double* ey, double* wt,
                                    double tau,
                                    double* rho, double* u, double* v,
                                    double* dPdt_x, double* dPdt_y,
                                    double* f, 
                                    float & rho_min,
                                    float & rho_max)
      {
        // update density and velocity
        rho_min = 1000.0;
        rho_max = 0.0;
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
            if (f_sum > rho_max) rho_max = (float) f_sum;
            if (f_sum < rho_min) rho_min = (float) f_sum;
            u[N] = fex_sum / rho[N];
            v[N] = fey_sum / rho[N];
          }
        }

        // periodic B.C. for rho (only used for plotting)
        for(int i = 0; i < NX-1; i++)
        {
          int j = NY-1; // top boundary
          int N_end = i*NY + j;
          int N_beg = i*NY + 0;
          rho[N_end] = rho[N_beg];
        }
        for(int j = 0; j < NY-1; j++)
        {
          int i = NX-1; // right boundary
          int N_end = i*NY + j;
          int N_beg = j;
          rho[N_end] = rho[N_beg];
        }
        rho[NX*NY-1] = rho[0];
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
      void updateEquilibrium(const int NX, const int NY,
                             double* ex, double* ey, double* wt,
                             const double* rho, 
                             const double* u, const double* v,
                             const double* dPdt_x, const double* dPdt_y,
                             double* f_eq, double* f_eq_tilda)
      {
        for(int i = 0; i < NX-1; i++)
        {
          for(int j = 0; j < NY-1; j++)
          {
            int N = i*NY + j;
            double udotu = u[N]*u[N] + v[N]*v[N];
            double u_tilda = u[N] + dPdt_x[N] / rho[N];
            double v_tilda = v[N] + dPdt_y[N] / rho[N];
            double utilda_dot_utilda = u_tilda*u_tilda + v_tilda*v_tilda;
            for(int id = 0; id < 9; id++)
            {
              int index_f = 9*N + id;
              double edotu = ex[id]*u[N] + ey[id]*v[N];
              double edotu_tilda = ex[id]*u_tilda + ey[id]*v_tilda;
              f_eq[index_f] = wt[id] * rho[N] 
                            * (1 + 3*edotu
                                 + 4.5*edotu*edotu - 1.5*udotu);
              f_eq_tilda[index_f] = wt[id] * rho[N] 
                            * (1 + 3*edotu_tilda
                                 + 4.5*edotu_tilda*edotu_tilda
                                 - 1.5*utilda_dot_utilda);
            }
          }
        }
      }

//    789012345678901234567890123456789012345678901234567890123456789012
//    ------------------------------------------------------------------
      int main(void)
      {
//      lattice size

        const int NX = 128;         // number of lattice points along X
        const int NY = 128;         // number of lattice points along Y

        // domain size in lattice units
        // grid spacing is unity along X and Y

        const double xmin = 0;
        const double xmax = NX-1;
        const double ymin = 0;
        const double ymax = NY-1;

//      example where NX = 8 and NY = 8
//
//
//       7 P-----P-----P-----P-----P-----P-----P-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       6 *-----*-----*-----*-----*-----*-----*-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       5 *-----*-----*-----*-----*-----*-----*-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       4 *-----*-----*-----*-----*-----*-----*-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       3 *-----*-----*-----*-----*-----*-----*-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       2 *-----*-----*-----*-----*-----*-----*-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       1 *-----*-----*-----*-----*-----*-----*-----P
//         |     |     |     |     |     |     |     |
//         |     |     |     |     |     |     |     |
//       0 *-----*-----*-----*-----*-----*-----*-----P
//         0     1     2     3     4     5     6     7
//
//
//         * = fields are calculated here
//         P = periodic boundary points

//      LBM parameters

        const double GEE11 = -0.55;   // interaction strength
        const double tau = 0.7;       // relaxation time
        const double rhoAvg = 0.693;  // reference density value

//      D2Q9 directions

        double ex[] = {0,1,0,-1,0,1,-1,-1,1}; 
        double ey[] = {0,0,1,0,-1,1,1,-1,-1}; 
        double wt[] = {4./9,1./9,1./9,1./9,1./9, 
                            1./36,1./36,1./36,1./36};
        double G11[] = {0,GEE11,GEE11,GEE11,GEE11,
                          GEE11/4,GEE11/4,GEE11/4,GEE11/4};

//      define buffers

        double *rho        = new double[NX*NY]; // density
        double *u          = new double[NX*NY]; // velocity x-component
        double *v          = new double[NX*NY]; // velocity y-component
        double *dPdt_x     = new double[NX*NY]; // momentum change along x
        double *dPdt_y     = new double[NX*NY]; // momentum change along y
        double *f          = new double[NX*NY*9]; // PDF
        double *f_eq       = new double[NX*NY*9]; // PDF
        double *f_new      = new double[NX*NY*9]; // PDF
        double *f_eq_tilda = new double[NX*NY*9]; // PDF

//      --------------------------------
//         Create a WINDOW using GLFW
//      --------------------------------

        GLFWwindow *window;

//      initialize the library

        if(!glfwInit()) return -1;

//      window size for displaying graphics

        int WIDTH  = 400;
        int HEIGHT = 400;

//      set the window's display mode

        window = glfwCreateWindow(WIDTH, HEIGHT, 
                                  "2D SCMP Simulation", NULL, NULL);

        if(!window) 
        {
          glfwTerminate();
          return -1;
        }

//      make the context current

        glfwMakeContextCurrent(window);

//      initialize fields

        initialize(NX, NY, rhoAvg, 
                   &ex[0], &ey[0], &wt[0], 
                   rho, u, v, f, f_new, f_eq);

        calc_dPdt(NX, NY, ex, ey, G11, rho, dPdt_x, dPdt_y);

        float rho_min, rho_max;

        updateDensityAndVelocity(NX, NY, ex, ey, wt, tau, 
                                 rho, u, v, dPdt_x, dPdt_y, f_new,
                                 rho_max, rho_max);

        updateEquilibrium(NX, NY, ex, ey, wt, rho, u, v, 
                          dPdt_x, dPdt_y, f_eq, f_eq_tilda);

//      time integration

        int time = 0;
        clock_t t0, tN;
        t0 = clock();


        //---------------------------------------
        // Loop until the user closes the window
        //---------------------------------------

        while(!glfwWindowShouldClose(window))
        {
          time++; // increment lattice time

          streaming(NX, NY, ex, ey, tau, f, f_new, f_eq, f_eq_tilda);

          calc_dPdt(NX, NY, ex, ey, G11, rho, dPdt_x, dPdt_y);

          updateDensityAndVelocity(NX, NY, ex, ey, wt, tau,
                                   rho, u, v, dPdt_x, dPdt_y, f_new,
                                   rho_min, rho_max);

          updateEquilibrium(NX, NY, ex, ey, wt, rho, u, v, 
                            dPdt_x, dPdt_y, f_eq, f_eq_tilda);

//        transfer fnew back to f

          for(int f_index = 0; f_index < NX*NY*9; f_index++)
          {
            f[f_index] = f_new[f_index];
          }

//        on-the-fly OpenGL graphics

          if(time%10 == 0) 
          {
            showGraphics(NX, NY, WIDTH, HEIGHT, 
                         xmin, xmax, ymin, ymax, rho);

            glfwSwapBuffers(window); // swap front and back buffers

            glfwPollEvents(); // poll for and processs events
          }

//        calculate the number of lattice time-steps per second

          tN = clock() - t0;

          std::cout << " lattice time steps per second = " 
                    << std::setw(5)
                    << (int) ((float) CLOCKS_PER_SEC * time / (float) tN)
                    << std::setw(10)
                    << " tau = " << tau
                    << std::setw(10)
                    << " GEE11 = " << GEE11
                    << std::setw(10)
                    << " min density = " << rho_min
                    << std::setw(10)
                    << " max density = " << rho_max
                    << std::setw(10)
                    << " density ratio = " << rho_max / rho_min
                    << std::endl;
        }

//      clean up

        delete[] rho;
        delete[] u;
        delete[] v;
        delete[] dPdt_x;
        delete[] dPdt_y;
        delete[] f;
        delete[] f_eq;
        delete[] f_new;

//      GLFW clean up

        glfwDestroyWindow(window);
        glfwTerminate();

//      main program ends

        return 0;
      }
//    ------------------------------------------------------------------
//    789012345678901234567890123456789012345678901234567890123456789012
