/*************************************************************************\

Ewha womans University

\**************************************************************************/
#include <iostream>//onlyfortest, huangxin


#include <stdlib.h>
#include <math.h>

// PQP
#include <MatVec.h>
#include <RectDist.h>
#include <OBB_Disjoint.h>

#include "C2A/C2A_BV.h"
#include "C2A/C2A_RectDist.h"
#include "C2A/C2A_Tri.h"

C2A_BV::C2A_BV()
: BV()
{
	trilength=0;  
}

C2A_BV::~C2A_BV()
{
}

static
inline 
PQP_REAL 
MaxOfTwo(PQP_REAL a, PQP_REAL b) 
{
  if (a > b) return a;
  return b;
}

void C2A_BV::ComputeCenterOfMassBV(C2A_Tri *tris, int num_tris)
{
    int i;
    PQP_REAL area_x2;
    PQP_REAL total_area;
	
    com[0] = 0.0;
    com[1] = 0.0;
    com[2] = 0.0;
	
    total_area = 0.0;
    for( i = 0; i < num_tris; i++ ) 
	{
		PQP_REAL E0[3], E1[3], S[3];
		VmV(E0, tris[i].p1, tris[i].p2);
		VmV(E1, tris[i].p1, tris[i].p3);
		VcrossV(S, E0, E1);
		area_x2=sqrt(S[0]*S[0] + S[1]*S[1] + S[2]*S[2]);
        total_area += area_x2;
		
        com[0] += (tris[i].p1[0] + tris[i].p2[0] + tris[i].p3[0]) * area_x2;
        com[1] += (tris[i].p1[1] + tris[i].p2[1] + tris[i].p3[1]) * area_x2;
		com[2] += (tris[i].p1[2] + tris[i].p2[2] + tris[i].p3[2]) * area_x2;
    }
	com[0] = com[0] / (3.0 * total_area);
	com[1] = com[1] / (3.0 * total_area);
	com[2] = com[2] / (3.0 * total_area);
	
}

// Calculate the maximum length between the vertex of triangles which belong the C2A_BV
void C2A_BV::ComputeAngularRadius(PQP_REAL O[3][3], C2A_Tri *tris, int num_tris)
{
    int i;
	angularRadius=0.0;
	PQP_REAL d;

	// comRoot is the center of whole model

    comRoot[0] = 0.0;
    comRoot[1] = 0.0;
    comRoot[2] = 0.0;
	
    for( i = 0; i < num_tris; i++ )
	{
		PQP_REAL p[3], pr[3];
		
//		MTxV(P[point],R,tris[i].p1);
		
		p[0] = tris[i].p1[0];//x
		p[1] = tris[i].p1[1];//y
		p[2] = tris[i].p1[2];//z	
		VmV(pr, p, comRoot);
		d = Vlength(pr);		
		if( d > angularRadius ) angularRadius = d;
		
		p[0] = tris[i].p2[0];
		p[1] = tris[i].p2[1];
		p[2] = tris[i].p2[2];		
		VmV(pr, p, comRoot);//pr=p-comRoot
		d = Vlength(pr);//d=|pr|	
		if( d > angularRadius ) angularRadius = d;
		
		p[0] = tris[i].p3[0];
		p[1] = tris[i].p3[1];
		p[2] = tris[i].p3[2];		
		VmV(pr, p, comRoot);
		d = Vlength(pr);	
		if( d > angularRadius ) angularRadius = d;	

    }
}

void
C2A_BV::FitToTris(PQP_REAL O[3][3], C2A_Tri *tris, int num_tris, PQP_REAL comR[3])
{

  VcV(comRoot, comR);//comRoot=comR
  ComputeCenterOfMassBV(tris, num_tris);//calculate the mass center of C2A_BV
  ComputeAngularRadius(O, tris, num_tris);//calculate the maximum distance of points in C2A_BV to the model center
  
  // store orientation  
  McM(R,O);//R=O

  // project points of tris to R coordinates

  int num_points = 3*num_tris;
  PQP_REAL (*P)[3] = new PQP_REAL[num_points][3];
  int point = 0;
  int i;
  for (i = 0; i < num_tris; i++) 
  {
    MTxV(P[point],R,tris[i].p1);
    point++;

    MTxV(P[point],R,tris[i].p2);
    point++;

    MTxV(P[point],R,tris[i].p3);
    point++;
  }

  PQP_REAL minx, maxx, miny, maxy, minz, maxz, c[3];

#if PQP_BV_TYPE & OBB_TYPE // if using OBB as C2A_BV
  minx = maxx = P[0][0];
  miny = maxy = P[0][1];
  minz = maxz = P[0][2];
  for (i = 1; i < num_points; i++)
  {
    if (P[i][0] < minx) minx = P[i][0];
    else if (P[i][0] > maxx) maxx = P[i][0];
    if (P[i][1] < miny) miny = P[i][1];
    else if (P[i][1] > maxy) maxy = P[i][1];
    if (P[i][2] < minz) minz = P[i][2];
    else if (P[i][2] > maxz) maxz = P[i][2];
  }
  c[0] = (PQP_REAL)0.5*(maxx + minx);
  c[1] = (PQP_REAL)0.5*(maxy + miny);
  c[2] = (PQP_REAL)0.5*(maxz + minz); 
  MxV(To,R,c);

  d[0] = (PQP_REAL)0.5*(maxx - minx);
  d[1] = (PQP_REAL)0.5*(maxy - miny);
  d[2] = (PQP_REAL)0.5*(maxz - minz);
#endif
  
#if PQP_BV_TYPE & RSS_TYPE //if use the RSS as C2A_BV

  // compute thickness, which determines radius, and z of rectangle corner
  
  PQP_REAL cz,radsqr;
  minz = maxz = P[0][2];
  for (i = 1; i < num_points; i++) 
  {
    if (P[i][2] < minz) minz = P[i][2];
    else if (P[i][2] > maxz) maxz = P[i][2];
  }
  r = (PQP_REAL)0.5*(maxz - minz);
  radsqr = r*r;
  cz = (PQP_REAL)0.5*(maxz + minz);

  // compute an initial length of rectangle along x direction

  // find minx and maxx as starting points

  int minindex, maxindex;
  minindex = maxindex = 0;
  for (i = 1; i < num_points; i++) 
  {
    if (P[i][0] < P[minindex][0]) minindex = i; 
    else if (P[i][0] > P[maxindex][0]) maxindex = i;
  }
  PQP_REAL x, dz;
  dz = P[minindex][2] - cz;
  minx = P[minindex][0] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
  dz = P[maxindex][2] - cz;
  maxx = P[maxindex][0] - sqrt(MaxOfTwo(radsqr - dz*dz,0));

  // grow minx

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][0] < minx) 
    {
      dz = P[i][2] - cz;
      x = P[i][0] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (x < minx) minx = x;
    }
  }

  // grow maxx

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][0] > maxx) 
    {
      dz = P[i][2] - cz;
      x = P[i][0] - sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (x > maxx) maxx = x;
    }
  }
  
  // compute an initial length of rectangle along y direction

  // find miny and maxy as starting points

  minindex = maxindex = 0;
  for (i = 1; i < num_points; i++) 
  {
    if (P[i][1] < P[minindex][1]) minindex = i;
    else if (P[i][1] > P[maxindex][1]) maxindex = i;
  }
  PQP_REAL y;
  dz = P[minindex][2] - cz;
  miny = P[minindex][1] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
  dz = P[maxindex][2] - cz;
  maxy = P[maxindex][1] - sqrt(MaxOfTwo(radsqr - dz*dz,0));

  // grow miny

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][1] < miny) 
    {
      dz = P[i][2] - cz;
      y = P[i][1] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (y < miny) miny = y;
    }
  }

  // grow maxy

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][1] > maxy) 
    {
      dz = P[i][2] - cz;
      y = P[i][1] - sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (y > maxy) maxy = y;
    }
  }
  
  // corners may have some points which are not covered - grow lengths if
  // necessary
  
  PQP_REAL dx, dy, u, t;
  PQP_REAL a = sqrt((PQP_REAL)0.5);
  for (i = 0; i < num_points; i++) 
  {
    if (P[i][0] > maxx) 
    {
      if (P[i][1] > maxy) 
      {
        dx = P[i][0] - maxx;
        dy = P[i][1] - maxy;
        u = dx*a + dy*a;
        t = (a*u - dx)*(a*u - dx) + 
            (a*u - dy)*(a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          maxx += u*a;
          maxy += u*a;
        }
      }
      else if (P[i][1] < miny) 
      {
        dx = P[i][0] - maxx;
        dy = P[i][1] - miny;
        u = dx*a - dy*a;
        t = (a*u - dx)*(a*u - dx) + 
            (-a*u - dy)*(-a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          maxx += u*a;
          miny -= u*a;
        }
      }
    }
    else if (P[i][0] < minx) 
    {
      if (P[i][1] > maxy) 
      {
        dx = P[i][0] - minx;
        dy = P[i][1] - maxy;
        u = dy*a - dx*a;
        t = (-a*u - dx)*(-a*u - dx) + 
            (a*u - dy)*(a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          minx -= u*a;
          maxy += u*a;
        }     
      }
      else if (P[i][1] < miny) 
      {
        dx = P[i][0] - minx;
        dy = P[i][1] - miny;
        u = -dx*a - dy*a;
        t = (-a*u - dx)*(-a*u - dx) + 
            (-a*u - dy)*(-a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          minx -= u*a; 
          miny -= u*a;
        }
      }
    }
  }

  c[0] = minx;
  c[1] = miny;
  c[2] = cz;
  //huangxinBV
  MxV(Tr,R,c);

  l[0] = maxx - minx;  
  if (l[0] < 0) l[0] = 0;
  l[1] = maxy - miny;
  if (l[1] < 0) l[1] = 0;
#endif

  delete [] P;
}


void
C2A_BV::FitToTris_Corner(PQP_REAL O[3][3], C2A_Tri *tris, int num_tris, PQP_REAL comR[3])
{

  VcV(comRoot, comR);//comRoot=comR
  ComputeCenterOfMassBV(tris, num_tris);//calculate the mass center of C2A_BV
  ComputeAngularRadius(O, tris, num_tris);//calculate the maximum distance of points in C2A_BV to the model center
  trilength=0;
  // store orientation  
  McM(R,O);//R=O

  // project points of tris to R coordinates

  int num_points = 3*num_tris;
  PQP_REAL (*P)[3] = new PQP_REAL[num_points][3];
  int point = 0;
  int i;
  for (i = 0; i < num_tris; i++) 
  {
    MTxV(P[point],R,tris[i].p1);
    point++;

    MTxV(P[point],R,tris[i].p2);
    point++;

    MTxV(P[point],R,tris[i].p3);
    point++;
  }
  if (num_tris==1)
  {
	  PQP_REAL v[3],a,b,c;
	  VmV(v,P[0],P[1]);
      a=Vlength(v);
	  VmV(v,P[0],P[2]);
      b=Vlength(v);
	  VmV(v,P[2],P[1]);
      c=Vlength(v);
	  trilength=a;
	  if (trilength<b) trilength=b;
	  if (trilength<c) trilength=c;
  }

  PQP_REAL minx, maxx, miny, maxy, minz, maxz, c[3];

#if PQP_BV_TYPE & OBB_TYPE // if using OBB as C2A_BV
  minx = maxx = P[0][0];
  miny = maxy = P[0][1];
  minz = maxz = P[0][2];
  for (i = 1; i < num_points; i++)
  {
    if (P[i][0] < minx) minx = P[i][0];
    else if (P[i][0] > maxx) maxx = P[i][0];
    if (P[i][1] < miny) miny = P[i][1];
    else if (P[i][1] > maxy) maxy = P[i][1];
    if (P[i][2] < minz) minz = P[i][2];
    else if (P[i][2] > maxz) maxz = P[i][2];
  }
  c[0] = (PQP_REAL)0.5*(maxx + minx);
  c[1] = (PQP_REAL)0.5*(maxy + miny);
  c[2] = (PQP_REAL)0.5*(maxz + minz); 
  MxV(To,R,c);

  d[0] = (PQP_REAL)0.5*(maxx - minx);
  d[1] = (PQP_REAL)0.5*(maxy - miny);
  d[2] = (PQP_REAL)0.5*(maxz - minz);
#endif
  
#if PQP_BV_TYPE & RSS_TYPE //if use the RSS as C2A_BV

  // compute thickness, which determines radius, and z of rectangle corner
  
  PQP_REAL cz,radsqr;
  minz = maxz = P[0][2];
  for (i = 1; i < num_points; i++) 
  {
    if (P[i][2] < minz) minz = P[i][2];
    else if (P[i][2] > maxz) maxz = P[i][2];
  }
  r = (PQP_REAL)0.5*(maxz - minz);
  radsqr = r*r;
  cz = (PQP_REAL)0.5*(maxz + minz);

  // compute an initial length of rectangle along x direction

  // find minx and maxx as starting points

  int minindex, maxindex;
  minindex = maxindex = 0;
  for (i = 1; i < num_points; i++) 
  {
    if (P[i][0] < P[minindex][0]) minindex = i; 
    else if (P[i][0] > P[maxindex][0]) maxindex = i;
  }
  PQP_REAL x, dz;
  dz = P[minindex][2] - cz;
  minx = P[minindex][0] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
  dz = P[maxindex][2] - cz;
  maxx = P[maxindex][0] - sqrt(MaxOfTwo(radsqr - dz*dz,0));

  // grow minx

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][0] < minx) 
    {
      dz = P[i][2] - cz;
      x = P[i][0] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (x < minx) minx = x;
    }
  }

  // grow maxx

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][0] > maxx) 
    {
      dz = P[i][2] - cz;
      x = P[i][0] - sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (x > maxx) maxx = x;
    }
  }
  
  // compute an initial length of rectangle along y direction

  // find miny and maxy as starting points

  minindex = maxindex = 0;
  for (i = 1; i < num_points; i++) 
  {
    if (P[i][1] < P[minindex][1]) minindex = i;
    else if (P[i][1] > P[maxindex][1]) maxindex = i;
  }
  PQP_REAL y;
  dz = P[minindex][2] - cz;
  miny = P[minindex][1] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
  dz = P[maxindex][2] - cz;
  maxy = P[maxindex][1] - sqrt(MaxOfTwo(radsqr - dz*dz,0));

  // grow miny

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][1] < miny) 
    {
      dz = P[i][2] - cz;
      y = P[i][1] + sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (y < miny) miny = y;
    }
  }

  // grow maxy

  for (i = 0; i < num_points; i++) 
  {
    if (P[i][1] > maxy) 
    {
      dz = P[i][2] - cz;
      y = P[i][1] - sqrt(MaxOfTwo(radsqr - dz*dz,0));
      if (y > maxy) maxy = y;
    }
  }
  
  // corners may have some points which are not covered - grow lengths if
  // necessary
  
  PQP_REAL dx, dy, u, t;
  PQP_REAL a = sqrt((PQP_REAL)0.5);
  for (i = 0; i < num_points; i++) 
  {
    if (P[i][0] > maxx) 
    {
      if (P[i][1] > maxy) 
      {
        dx = P[i][0] - maxx;
        dy = P[i][1] - maxy;
        u = dx*a + dy*a;
        t = (a*u - dx)*(a*u - dx) + 
            (a*u - dy)*(a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          maxx += u*a;
          maxy += u*a;
        }
      }
      else if (P[i][1] < miny) 
      {
        dx = P[i][0] - maxx;
        dy = P[i][1] - miny;
        u = dx*a - dy*a;
        t = (a*u - dx)*(a*u - dx) + 
            (-a*u - dy)*(-a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          maxx += u*a;
          miny -= u*a;
        }
      }
    }
    else if (P[i][0] < minx) 
    {
      if (P[i][1] > maxy) 
      {
        dx = P[i][0] - minx;
        dy = P[i][1] - maxy;
        u = dy*a - dx*a;
        t = (-a*u - dx)*(-a*u - dx) + 
            (a*u - dy)*(a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          minx -= u*a;
          maxy += u*a;
        }     
      }
      else if (P[i][1] < miny) 
      {
        dx = P[i][0] - minx;
        dy = P[i][1] - miny;
        u = -dx*a - dy*a;
        t = (-a*u - dx)*(-a*u - dx) + 
            (-a*u - dy)*(-a*u - dy) +
            (cz - P[i][2])*(cz - P[i][2]);
        u = u - sqrt(MaxOfTwo(radsqr - t,0));
        if (u > 0) 
        {
          minx -= u*a; 
          miny -= u*a;
        }
      }
    }
  }

  PQP_REAL Corner_temp[4][3];

  Corner_temp[0][0] = minx;
  Corner_temp[0][1] = miny;
  Corner_temp[0][2] = cz;
  //
  Corner_temp[1][0] = maxx;
  Corner_temp[1][1] = miny;
  Corner_temp[1][2] = cz;
  //
  Corner_temp[2][0] = minx;
  Corner_temp[2][1] = maxy;
  Corner_temp[2][2] = cz;
  
  Corner_temp[3][0] = maxx;
  Corner_temp[3][1] = maxy;
  Corner_temp[3][2] = cz;
  


  MxV(Corner[0],R,Corner_temp[0]);
  MxV(Corner[2],R,Corner_temp[2]);
  Length_C1 = Vlength(Corner[0]) < Vlength(Corner[2]) ? Vlength(Corner[0]) : Vlength(Corner[2]);

  MxV(Corner[1],R,Corner_temp[1]);
  MxV(Corner[3],R,Corner_temp[3]);
  Length_C2 = Vlength(Corner[1]) < Vlength(Corner[3])? Vlength(Corner[1]) : Vlength(Corner[3]);

  MxV(Corner[2],R,Corner_temp[2]);

  VcV(Tr,Corner[0]);


  max_Length=Length_C1>Length_C2 ? Length_C1:Length_C2;
  C2A_BV_Size=GetSize();

  l[0] = maxx - minx;  
  if (l[0] < 0) l[0] = 0;
  l[1] = maxy - miny;
  if (l[1] < 0) l[1] = 0;
#endif

  delete [] P;
}




int 
C2A_BV_Overlap(PQP_REAL R[3][3], PQP_REAL T[3], C2A_BV *b1, C2A_BV *b2)
{
#if PQP_BV_TYPE & OBB_TYPE
  return (obb_disjoint(R,T,b1->d,b2->d) == 0);
#else
  PQP_REAL dist = RectDist(R,T,b1->l,b2->l);
  if (dist <= (b1->r + b2->r)) return 1;
  return 0;
#endif
}

#if PQP_BV_TYPE & RSS_TYPE

PQP_REAL
C2A_BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], C2A_BV *b1, C2A_BV *b2,PQP_REAL S[3])
{
	PQP_REAL P[3];
	PQP_REAL Q[3];
 	bool bValidPQ;

  PQP_REAL dist = C2ARectDist(R,T,b1->l,b2->l, P, Q, bValidPQ, S);
  dist -= (b1->r + b2->r);
  return (dist < (PQP_REAL)0.0)? (PQP_REAL)0.0 : dist;
}

PQP_REAL
C2A_BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], C2A_BV *b1, C2A_BV *b2,
			PQP_REAL P[3], PQP_REAL Q[3],
			bool &bValidPQ,PQP_REAL S[3])
{
  // The closest points between the two rectangles not the two RSS!!!
  PQP_REAL dist = C2ARectDist(R,T,b1->l,b2->l, P, Q, bValidPQ, S);
  dist -= (b1->r + b2->r);
  return (dist < (PQP_REAL)0.0)? (PQP_REAL)0.0 : dist;
}
#endif



