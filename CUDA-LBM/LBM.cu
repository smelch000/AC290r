/*
 * 2D flow around a cylinder
 *
 * compile with nvcc -O2 LBM.cu -o LBMcuda
 *
*/
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#define NBLOCKS_X 4
#define NBLOCKS_Y 4
#define NTHREADS_X 32
#define NTHREADS_Y 16


#define PI2 6.28318530718
//====== Flow parameters definition ==============================================
#define MAXITER 500                  // Total number of time iterations.
#define OUTSTEP 10
#define Re 220.0                 // Reynolds number.
#define NX 520                       // Lattice dimensions and populations.
#define NY 180 
#define LY (NY-1.0)
#define Q 9  
#define Q0 3
#define CX (NX/4)                      // Coordinates and radius of the cylinder.
#define CY (NY/2)
#define R (NY/9) 
#define ULB 0.04                  // Velocity in lattice units.


typedef float real_t;
typedef unsigned int uint;

//----- Lattice Constants -------------------------------

__device__ __constant__ real_t C[Q][2] = {
      { 0., 0.},
      { 0.,-1.},
      { 0., 1.},
      {-1., 0.},
      {-1.,-1.},
      {-1., 1.},
      { 1., 0.},
      { 1.,-1.},
      { 1., 1.}
};
__device__ __constant__ int iC[Q][2] = {
      { 0, 0},
      { 0,-1},
      { 0, 1},
      {-1, 0},
      {-1,-1},
      {-1, 1},
      { 1, 0},
      { 1,-1},
      { 1, 1}
};

real_t C_h[Q][2]; //C on host

//noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)] 
__device__ __constant__ int NOSLIP[Q]={0, 2, 1, 6, 8, 7, 3, 5, 4}; //noslip indexes for C

//i1 = arange(Q)[asarray([ci[0]<0  for ci in c])] # Unknown on right wall.
__device__ __constant__ int I1[Q0] = {3, 4, 5}; // Pops unknown on right wall.
//i2 = arange(Q)[asarray([ci[0]==0 for ci in c])] # Vertical middle.
__device__ __constant__ int I2[Q0] = {0,1,2}; // Vertical middle.
//i3 = arange(Q)[asarray([ci[0]>0  for ci in c])] # Unknown on left wall.
__device__ __constant__ int I3[Q0] = {6, 7, 8}; //Pops Unknown on left wall.


const real_t NULB_h = ULB*(real_t)(R)/(real_t)(Re);
const real_t OMEGA_h = 1.0 / (3.*NULB_h+0.5); // Relaxation parameter.

__device__ __constant__ real_t NULB;
__device__ __constant__ real_t OMEGA; // Relaxation parameter.

//----------------------------------------------------------------------

//========= Functions declaration =================


__host__ __device__ real_t norm2(const real_t * v);

__device__ bool isObstacle(int x, int y);

//convert 2D to 1D array address
// A[i][j]
__device__ __forceinline__ uint i2D(uint i, uint j) {return (NY*i+j);}
//convert 3D to 1D array address
// A[m][i][j]
__device__ __forceinline__ uint i3D(uint m, uint i, uint j) {return (m*NX*NY+NY*i+j);}

__global__ void initialize(real_t * vel, real_t * rho);

__global__ void getEquilibrium(const real_t* rho,
                 const real_t* u,
                 const real_t* t,
                 real_t* feq
                );

__global__ void getHvars(real_t* rho,
            real_t* u,
            const real_t * f
     );

// Right wall: apply outflow condition.
__global__ void outflow(real_t *fin);

// Left wall: compute density from known populations.
__global__ void leftwall(real_t* rho, real_t *u, const real_t* vel, const real_t* fin);

// Left wall: Zou/He boundary condition.
__global__ void zouhe(const real_t* feq, real_t* fin);

// Collision step.
__global__ void collision(const real_t* fin, real_t* fout, const real_t* feq);

//Wall rebound
__global__ void rebound(const real_t * fin, real_t* fout);

// Streaming step.
__global__ void streaming(real_t* fin, const real_t* fout);

//--------------------------------------------------


int main()
{

  // copy NULB_h to NULB on device  
  cudaMemcpyToSymbol(NULB, &NULB_h, sizeof(real_t), 0, cudaMemcpyHostToDevice);  
  // copy OMEGA_h to OMEGA on device
  cudaMemcpyToSymbol(OMEGA, &OMEGA_h, sizeof(real_t), 0, cudaMemcpyHostToDevice);  
  // copy C to C_h on host
  cudaMemcpyFromSymbol(&C_h, C, 2*Q*sizeof(real_t), 0, cudaMemcpyDeviceToHost);

  real_t t[Q];
  real_t* t_d;
  //allocation on device
  cudaMalloc((void **)&t_d, Q * sizeof(real_t));

  real_t rho[NX][NY];
  real_t* rho_d;
  //allocation on device
  cudaMalloc((void **)&rho_d, NX*NY * sizeof(real_t));

  real_t* vel_d;
  //allocation on device
  cudaMalloc((void **)&vel_d, 2*NX*NY * sizeof(real_t));

  real_t u[2][NX][NY];
  real_t* u_d;
  //allocation on device
  cudaMalloc((void **)&u_d, 2*NX*NY * sizeof(real_t));

  real_t* feq_d;
  //allocation on device
  cudaMalloc((void **)&feq_d, Q*NX*NY * sizeof(real_t));

  real_t* fin_d;
  //allocation on device
  cudaMalloc((void **)&fin_d, Q*NX*NY * sizeof(real_t));

  real_t* fout_d;
  //allocation on device
  cudaMalloc((void **)&fout_d, Q*NX*NY * sizeof(real_t));
  
   
  t[0]=4./9.;
  for (int iq=1; iq<Q; iq++)
  {
    if (norm2(&C_h[iq][0])<2.)
    {
      t[iq]=1./9.;
    }
    else
    {
      t[iq]=1./36.;
    }
  }
  cudaMemcpy((real_t *)t_d,t,Q*sizeof(real_t),cudaMemcpyHostToDevice);

  //initial velocity and density setup  
  { 
    dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
    dim3 threads(NTHREADS_X,NTHREADS_Y,1);
    initialize<<<grid,threads>>>(vel_d,rho_d);
  }
  //-------------------------------

  //equilibrium DF setup
  {
    dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
    dim3 threads(NTHREADS_X,NTHREADS_Y,1);
    getEquilibrium<<<grid,threads>>>(rho_d,vel_d,t_d,feq_d);
  }
  //-------------------------------

  //copy feq to fin
  cudaMemcpy( fin_d, feq_d, sizeof(real_t)*Q*NX*NY,cudaMemcpyDeviceToDevice );
  
  //################################################################################
  //###### Main time loop ##########################################################
  
  for (int time=0; time < MAXITER; time++)
  {
    // Right wall: apply outflow condition.
    {
      dim3 grid(1,NBLOCKS_Y,1);
      dim3 threads(1,NTHREADS_Y,1);
      outflow<<<grid,threads>>>(fin_d);
    }
    //---------------------------------

    // Calculate macroscopic density and velocity.
    {
      dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
      dim3 threads(NTHREADS_X,NTHREADS_Y,1);
      getHvars<<<grid,threads>>>(rho_d,u_d,fin_d);
    }
    //---------------------------------

    
    // Left wall: compute density from known populations.
    {
      dim3 grid(1,NBLOCKS_Y,1);
      dim3 threads(1,NTHREADS_Y,1);
      leftwall<<<grid,threads>>>(rho_d, u_d, vel_d, fin_d);
    }
    //---------------------------------
        
    {
      dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
      dim3 threads(NTHREADS_X,NTHREADS_Y,1);
      getEquilibrium<<<grid,threads>>>(rho_d,u_d,t_d,feq_d);
    }
    //---------------------------------
    
    // Left wall: Zou/He boundary condition.
    {
      dim3 grid(1,NBLOCKS_Y,1);
      dim3 threads(1,NTHREADS_Y,1);
      zouhe<<<grid,threads>>>(feq_d,fin_d);
    }
    //---------------------------------
    
    // Collision step.
    {
      dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
      dim3 threads(NTHREADS_X,NTHREADS_Y,1);
      collision<<<grid,threads>>>(fin_d,fout_d,feq_d);
    }
    //---------------------------------

    // Wall "rebound" step.
    {
      dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
      dim3 threads(NTHREADS_X,NTHREADS_Y,1);
      rebound<<<grid,threads>>>(fin_d,fout_d);
    }
    //---------------------------------
    
    // Streaming step.
    {
      dim3 grid(NBLOCKS_X,NBLOCKS_Y,1);
      dim3 threads(NTHREADS_X,NTHREADS_Y,1);
      streaming<<<grid,threads>>>(fin_d, fout_d);
    }    
    //---------------------------------
    
    // Output.
    if (time % OUTSTEP==0)
    {

      //copy back density and macroscopic velocity from device 
      cudaMemcpy(rho,rho_d,NX*NY*sizeof(real_t),cudaMemcpyDeviceToHost);
      cudaMemcpy(u,u_d,2*NX*NY*sizeof(real_t),cudaMemcpyDeviceToHost);

      std::ofstream fout;
      std::stringstream filename("");
      filename <<"out"<<std::setw(5)<<std::setfill('0')<<time<<".csv";
      fout.open(filename.str().c_str(),std::ofstream::out);
      fout << "x,y,rho,ux,uy,|u|" << std::endl;
      for (int x=0; x<NX;x++)
      {
        for (int y=0; y<NY;y++)
        {
          fout << x << "," << y << "," <<rho[x][y] <<"," <<u[0][x][y] << "," << u[1][x][y] << "," 
                << sqrt(norm2(&u[0][x][y])) << std::endl;
        }
      }
      fout.close();
      std::cout << "Step "<<time<<" done.\n";
    }
    
  }
  //=========  free device memory ==========
  cudaFree(t_d);

  cudaFree(rho_d);

  cudaFree(vel_d);

  cudaFree(u_d);

  cudaFree(feq_d);

  cudaFree(fin_d);

  cudaFree(fout_d);
  //----------------------------------------
  
}      
//============================================================
//====== functions definition ================================

__host__ __device__ real_t norm2(const real_t * v)
{
  return (v[0]*v[0]+v[1]*v[1]);
} 

__device__ bool isObstacle(int x, int y)
{
  real_t xx[2];
  xx[0]=x-CX; xx[1]=y-CY;
  return norm2(xx)<R*R;
}

// recurrent declarations of threads indexes inside kernels
#define CUDAHEADER_X \
  const uint num_threads_x = gridDim.x*blockDim.x; \
  const uint tstart_x = blockDim.x*blockIdx.x+threadIdx.x;
  
#define CUDAHEADER_Y \
  const uint num_threads_y = gridDim.y*blockDim.y; \
  const uint tstart_y = blockDim.y*blockIdx.y+threadIdx.y;
  
#define CUDAHEADER \
  CUDAHEADER_X \
  CUDAHEADER_Y
  
//-------------------------------------------------------

__global__ void initialize(real_t * vel, real_t * rho)
{
  CUDAHEADER
  
  for (uint x= tstart_x; x<NX; x += num_threads_x)
  {
    for (uint y= tstart_y; y<NY; y += num_threads_y)
    {
      vel[i3D(0,x,y)]=0.;
      if (x>9 && x<20) vel[i3D(0,x,y)]=ULB;
      vel[i3D(1,x,y)]=0.;
      rho[i2D(x,y)]=1.;
    }
  }
}

// Equilibrium distribution function.
__global__ void getEquilibrium(const real_t* rho,
                 const real_t* u,
                 const real_t* t,
                 real_t* feq
                )
{
  CUDAHEADER
  
  real_t cu;
  real_t uxy[2];
  for (uint x= tstart_x; x<NX; x += num_threads_x)
  {
    for (uint y= tstart_y; y<NY; y += num_threads_y)
    {
      uxy[0]=u[i3D(0,x,y)]; uxy[1]=u[i3D(1,x,y)];
      for (uint iq=0; iq<Q;iq++)
      {
        cu = 3.0*(C[iq][0]*uxy[0]+C[iq][1]*uxy[1]);
        feq[i3D(iq,x,y)] = rho[i2D(x,y)]*t[iq]*(1.+cu+0.5*cu*cu-1.5*norm2(uxy));
      }
    }
  }
}

__global__ void getHvars(real_t* rho,
            real_t* u,
            const real_t * f
     )
{
  
  CUDAHEADER
  
  real_t ff;
  for (uint x= tstart_x; x<NX; x += num_threads_x)
  {
    for (uint y= tstart_y; y<NY; y += num_threads_y)
    {
      real_t& rhxy=rho[i2D(x,y)];
      rhxy=0.;
      u[i3D(0,x,y)]=u[i3D(1,x,y)]=0.;
      
      for (uint iq=0; iq<Q;iq++)
      {
        ff=f[i3D(iq,x,y)];
        rhxy += ff;
        u[i3D(0,x,y)] += C[iq][0]*ff;
        u[i3D(1,x,y)] += C[iq][1]*ff;
      }
      u[i3D(0,x,y)] /= rhxy; u[i3D(1,x,y)] /= rhxy;
      
    }
  }
}

// Right wall: apply outflow condition.
__global__ void outflow(real_t *fin)
{

  CUDAHEADER_Y
  
  for (uint y= tstart_y; y<NY; y += num_threads_y)
  {
    for (uint iq=0; iq<Q0;iq++)
    {
      fin[i3D(I1[iq],NX-1,y)] = fin[i3D(I1[iq],NX-2,y)];
    }
  }
}

// Left wall: compute density from known populations.
__global__ void leftwall(real_t* rho, real_t *u, const real_t* vel, const real_t* fin)
{

  CUDAHEADER_Y

  for (uint y= tstart_y; y<NY; y += num_threads_y)
  {
    u[i3D(0,0,y)] =vel[i3D(0,0,y)]; u[i3D(1,0,y)] =vel[i3D(1,0,y)]; 
    real_t &rh0y = rho[i2D(0,y)];
    rh0y = 0.;
    for (uint iq=0; iq<Q0;iq++)
    {
      rh0y += fin[i3D(I2[iq],0,y)] + 2.*fin[i3D(I1[iq],0,y)];
    }
    rh0y /= (1.-u[i3D(0,0,y)]);
  }
}

// Left wall: Zou/He boundary condition.
__global__ void zouhe(const real_t* feq, real_t* fin)
{
  CUDAHEADER_Y
  
  for (uint y= tstart_y; y<NY; y += num_threads_y)
  {
    for (uint iq=0; iq<Q0;iq++)
    {
      fin[i3D(I3[iq],0,y)] = fin[i3D(I1[iq],0,y)] + feq[i3D(I3[iq],0,y)] - feq[i3D(I1[iq],0,y)];
    }
  }
}

// Collision step.
__global__ void collision(const real_t* fin, real_t* fout, const real_t* feq)
{
  CUDAHEADER
  
  for (uint x= tstart_x; x<NX; x += num_threads_x)
  {
    for (uint y= tstart_y; y<NY; y += num_threads_y)
    {
      for (uint iq=0; iq<Q;iq++)
      {
          uint i=i3D(iq,x,y);
          fout[i] = fin[i] - OMEGA * (fin[i] - feq[i]);
      }
    }
  }
}

//Wall rebound
__global__ void rebound(const real_t * fin, real_t* fout)
{
  CUDAHEADER
  
  for (uint x= tstart_x; x<NX; x += num_threads_x)
  {
    for (uint y= tstart_y; y<NY; y += num_threads_y)
    {
      if (isObstacle(x,y))
      {
        for (uint iq=0; iq<Q;iq++)
        {
          fout[i3D(iq,x,y)] = fin[i3D(NOSLIP[iq],x,y)];
        }
      }
    }
  }
}

// Streaming step.
__global__ void streaming(real_t* fin, const real_t* fout)
{
  CUDAHEADER
  int xout,yout;
  
  for (int x= tstart_x; x<NX; x += num_threads_x)
  {
    for (int y= tstart_y; y<NY; y += num_threads_y)
    {
      fin[i3D(0,x,y)]=fout[i3D(0,x,y)];
      for (int iq=1; iq<Q;iq++)
      {
          //handle periodic conditions
          xout = ((x + iC[iq][0])+NX) % NX; 
          yout = ((y + iC[iq][1])+NY) % NY;
          fin[i3D(iq,xout,yout)]=fout[i3D(iq,x,y)];
      }
    }
  }
}  


