#include "PQP.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;



#include "C2A/C2A.h"

#include "BVTQ.h"
#include "Build.h"
#include "MatVec.h"
#include "GetTime.h"
#include "TriDist.h"
#include "C2A/ListTri.h"
#include "C2A/InterpMotion.h"
#include "C2A/C2A_Internal.h"



#ifdef _WIN32
#include <float.h>
#define isnan _isnan
#endif

#define Max_Value 1e+30
#define TanslationCCD    

bool bPQPContact_Represenative = true;
double PQPContact_radius = 5;
bool b_TanslationCCD = false;

// Defined in PQP.cpp
int
TriContact(PQP_REAL *P1, PQP_REAL *P2, PQP_REAL *P3,
		   PQP_REAL *Q1, PQP_REAL *Q2, PQP_REAL *Q3); 


PQP_REAL
TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], Tri *t1, Tri *t2,
			PQP_REAL p[3], PQP_REAL q[3]);


int C2A_TimeOfContactStep(CInterpMotion *objmotion1, CInterpMotion *objmotion2, C2A_TimeOfContactResult *res, 
			PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
			PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
			PQP_REAL tolerance_t,PQP_REAL tolerance_d);


enum BUILD_STATE
{ 
  PQP_BUILD_STATE_EMPTY,     // empty state, immediately after constructor
  PQP_BUILD_STATE_BEGUN,     // after BeginModel(), state for adding triangles
  PQP_BUILD_STATE_PROCESSED  // after tree has been built, ready to use
};

static void
ClosestPoints(PQP_REAL VEC[3], 
	  PQP_REAL X[3], PQP_REAL Y[3],             // closest points
          const PQP_REAL P[3], const PQP_REAL A[3], // seg 1 origin, vector
          const PQP_REAL Q[3], const PQP_REAL B[3]) // seg 2 origin, vector
{
  PQP_REAL T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
  PQP_REAL TMP[3];

  VmV(T,Q,P);
  A_dot_A = VdotV(A,A);
  B_dot_B = VdotV(B,B);
  A_dot_B = VdotV(A,B);
  A_dot_T = VdotV(A,T);
  B_dot_T = VdotV(B,T);

  // t parameterizes ray P,A 
  // u parameterizes ray Q,B 

  PQP_REAL t,u;

  // compute t for the closest point on ray P,A to
  // ray Q,B

  PQP_REAL denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;

  t = (A_dot_T*B_dot_B - B_dot_T*A_dot_B) / denom;

  // clamp result so t is on the segment P,A

  if ((t < 0) || isnan(t)) t = 0; else if (t > 1) t = 1;

  // find u for point on ray Q,B closest to point at t

  u = (t*A_dot_B - B_dot_T) / B_dot_B;

  // if u is on segment Q,B, t and u correspond to 
  // closest points, otherwise, clamp u, recompute and
  // clamp t 

  if ((u <= 0) || isnan(u)) {

    VcV(Y, Q);

    t = A_dot_T / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Q, P);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Q, X);
    }
    else {
      VpVxS(X, P, A, t);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  }
  else if (u >= 1) {

    VpV(Y, Q, B);

    t = (A_dot_B + A_dot_T) / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Y, P);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Y, X);
    }
    else {
      VpVxS(X, P, A, t);
      VmV(T, Y, P);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  }
  else {

    VpVxS(Y, Q, B, u);

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(T, Q, X);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    }
    else {
      VpVxS(X, P, A, t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, T) < 0) {
        VxS(VEC, VEC, -1);
      }
    }
  }
}

PQP_REAL 
TriDist(PQP_REAL P[3], PQP_REAL Q[3],
        const PQP_REAL S[3][3], const PQP_REAL T[3][3], C2A_ContactFeature &f1, C2A_ContactFeature &f2, bool &bCollided)  
{
  // Compute vectors along the 6 sides

  PQP_REAL Sv[3][3], Tv[3][3];
  PQP_REAL VEC[3];

  VmV(Sv[0],S[1],S[0]);
  VmV(Sv[1],S[2],S[1]);
  VmV(Sv[2],S[0],S[2]);

  VmV(Tv[0],T[1],T[0]);
  VmV(Tv[1],T[2],T[1]);
  VmV(Tv[2],T[0],T[2]);

  // For each edge pair, the vector connecting the closest points 
  // of the edges defines a slab (parallel planes at head and tail
  // enclose the slab). If we can show that the off-edge vertex of 
  // each triangle is outside of the slab, then the closest points
  // of the edges are the closest points for the triangles.
  // Even if these tests fail, it may be helpful to know the closest
  // points found, and whether the triangles were shown disjoint

  PQP_REAL V[3];
  PQP_REAL Z[3];
  PQP_REAL minP[3], minQ[3], mindd;
  int shown_disjoint = 0;

  mindd = VdistV2(S[0],T[0]) + 1;  // Set first minimum safely high

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      // Find closest points on edges i & j, plus the 
      // vector (and distance squared) between these points

      ClosestPoints(VEC,P,Q,S[i],Sv[i],T[j],Tv[j]);
      
      VmV(V,Q,P);
      PQP_REAL dd = VdotV(V,V);

      // Verify this closest point pair only if the distance 
      // squared is less than the minimum found thus far.

      if (dd <= mindd)
      {
        VcV(minP,P);
        VcV(minQ,Q);
        mindd = dd;

		f1.type = 1;
		f1.fid[0] = i;
		f2.type = 1;
		f2.fid[0] = j;

        VmV(Z,S[(i+2)%3],P);
        PQP_REAL a = VdotV(Z,VEC);
        VmV(Z,T[(j+2)%3],Q);
        PQP_REAL b = VdotV(Z,VEC);

        if ((a <= 0) && (b >= 0)) return sqrt(dd);

        PQP_REAL p = VdotV(V, VEC);

        if (a < 0) a = 0;
        if (b > 0) b = 0;
        if ((p - a + b) > 0) shown_disjoint = 1;	
      }
    }
  }

  // No edge pairs contained the closest points.  
  // either:
  // 1. one of the closest points is a vertex, and the
  //    other point is interior to a face.
  // 2. the triangles are overlapping.
  // 3. an edge of one triangle is parallel to the other's face. If
  //    cases 1 and 2 are not true, then the closest points from the 9
  //    edge pairs checks above can be taken as closest points for the
  //    triangles.
  // 4. possibly, the triangles were degenerate.  When the 
  //    triangle points are nearly colinear or coincident, one 
  //    of above tests might fail even though the edges tested
  //    contain the closest points.

  // First check for case 1

  PQP_REAL Sn[3], Snl;       
  VcrossV(Sn,Sv[0],Sv[1]); // Compute normal to S triangle
  Snl = VdotV(Sn,Sn);      // Compute square of length of normal
  
  // If cross product is long enough,

  if (Snl > 1e-15)  
  {
    // Get projection lengths of T points

    PQP_REAL Tp[3]; 

    VmV(V,S[0],T[0]);
    Tp[0] = VdotV(V,Sn);

    VmV(V,S[0],T[1]);
    Tp[1] = VdotV(V,Sn);

    VmV(V,S[0],T[2]);
    Tp[2] = VdotV(V,Sn);

    // If Sn is a separating direction,
    // find point with smallest projection

    int point = -1;
    if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0))
    {
      if (Tp[0] < Tp[1]) point = 0; else point = 1;
      if (Tp[2] < Tp[point]) point = 2;
    }
    else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0))
    {
      if (Tp[0] > Tp[1]) point = 0; else point = 1;
      if (Tp[2] > Tp[point]) point = 2;
    }

    // If Sn is a separating direction, 

    if (point >= 0) 
    {
      shown_disjoint = 1;

      // Test whether the point found, when projected onto the 
      // other triangle, lies within the face.
    
      VmV(V,T[point],S[0]);
      VcrossV(Z,Sn,Sv[0]);
      if (VdotV(V,Z) > 0)
      {
        VmV(V,T[point],S[1]);
        VcrossV(Z,Sn,Sv[1]);
        if (VdotV(V,Z) > 0)
        {
          VmV(V,T[point],S[2]);
          VcrossV(Z,Sn,Sv[2]);
          if (VdotV(V,Z) > 0)
          {
            // T[point] passed the test - it's a closest point for 
            // the T triangle; the other point is on the face of S

			  f1.type = 2;
			  f1.fid[0] = -1;
			  f2.type = 0;
			  f2.fid[0] = point;
			  
            VpVxS(P,T[point],Sn,Tp[point]/Snl);
            VcV(Q,T[point]);
            return sqrt(VdistV2(P,Q));
          }
        }
      }
    }
  }

  PQP_REAL Tn[3], Tnl;       
  VcrossV(Tn,Tv[0],Tv[1]); 
  Tnl = VdotV(Tn,Tn);      
  
  if (Tnl > 1e-15)  
  {
    PQP_REAL Sp[3]; 

    VmV(V,T[0],S[0]);
    Sp[0] = VdotV(V,Tn);

    VmV(V,T[0],S[1]);
    Sp[1] = VdotV(V,Tn);

    VmV(V,T[0],S[2]);
    Sp[2] = VdotV(V,Tn);

    int point = -1;
    if ((Sp[0] > 0) && (Sp[1] > 0) && (Sp[2] > 0))
    {
      if (Sp[0] < Sp[1]) point = 0; else point = 1;
      if (Sp[2] < Sp[point]) point = 2;
    }
    else if ((Sp[0] < 0) && (Sp[1] < 0) && (Sp[2] < 0))
    {
      if (Sp[0] > Sp[1]) point = 0; else point = 1;
      if (Sp[2] > Sp[point]) point = 2;
    }

    if (point >= 0) 
    { 
      shown_disjoint = 1;

      VmV(V,S[point],T[0]);
      VcrossV(Z,Tn,Tv[0]);
      if (VdotV(V,Z) > 0)
      {
        VmV(V,S[point],T[1]);
        VcrossV(Z,Tn,Tv[1]);
        if (VdotV(V,Z) > 0)
        {
          VmV(V,S[point],T[2]);
          VcrossV(Z,Tn,Tv[2]);
          if (VdotV(V,Z) > 0)
          {
			  f1.type = 0;
			  f1.fid[0] = point;
			  f2.type = 2;
			  f2.fid[0] = -1;
			  
            VcV(P,S[point]);
            VpVxS(Q,S[point],Tn,Sp[point]/Tnl);
            return sqrt(VdistV2(P,Q));
          }
        }
      }
    }
  }

  // Case 1 can't be shown.
  // If one of these tests showed the triangles disjoint,
  // we assume case 3 or 4, otherwise we conclude case 2, 
  // that the triangles overlap.
  
  if (shown_disjoint)
  {
    VcV(P,minP);
    VcV(Q,minQ);
    return sqrt(mindd);
  }
  else
  {
	  bCollided = true;
	   return 0;
  }
  
}


PQP_REAL
C2A_TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], Tri *t1, Tri *t2,
            PQP_REAL p[3], PQP_REAL q[3], C2A_ContactFeature &f1, C2A_ContactFeature &f2, bool &bCollided)
{
	// transform tri 2 into same space as tri 1
	
	PQP_REAL tri1[3][3], tri2[3][3];
	
	VcV(tri1[0], t1->p1);
	VcV(tri1[1], t1->p2);
	VcV(tri1[2], t1->p3);
	MxVpV(tri2[0], R, t2->p1, T);
	MxVpV(tri2[1], R, t2->p2, T);
	MxVpV(tri2[2], R, t2->p3, T);
	
	return TriDist(p,q,tri1,tri2, f1, f2, bCollided);
}
// CCD algorithm based on conservative advancement for polygonal soup models 


void
ContactQueueRecurse(C2A_ContactResult *res, 
                     PQP_REAL R[3][3], PQP_REAL T[3],
                     C2A_Model *o1, int b1,
                     C2A_Model *o2, int b2)
{
	if (res->modelcollided) 
	{
		return;
	}
  BVTQ bvtq(res->qsize);

  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R,R);
  VcV(min_test.T,T);

  while(1) 
  {  
    int l1 = o1->child(min_test.b1)->Leaf();
    int l2 = o2->child(min_test.b2)->Leaf();
    
    if (l1 && l2) 
    {  
      // both leaves.  Test the triangles beneath them.

      res->num_tri_tests++;

      PQP_REAL p[3], q[3];

      Tri *t1 = o1->GetTriangle(-o1->child(min_test.b1)->first_child - 1);
      Tri *t2 = o2->GetTriangle(-o2->child(min_test.b2)->first_child - 1);
	  
	  C2A_ContactFeature f1, f2;
      PQP_REAL d = C2A_TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,res->modelcollided);
  
	  if (d <= res->ctolerance) 
	  {
		  //add the contact feature if < ctolerance
		  f1.tri = t1->id;
		  f2.tri = t2->id;
		  res->Add(f1, f2, d, p, q);
		  //res->closer_than_ctolerance = true;
		  
		  VcV(res->p1, p);         // p already in c.s. 1
		  VcV(res->p2, q);         // q must be transformed 
		  
		  o1->last_tri = t1;
		  o2->last_tri = t2;		  
	  }
    }		 
    else if (bvtq.GetNumTests() == bvtq.GetSize() - 1) 
    {  
      // queue can't get two more tests, recur
      
      ContactQueueRecurse(res,min_test.R,min_test.T,
                           o1,min_test.b1,o2,min_test.b2);
	  if (res->modelcollided) 
	  {
		  return;
	  }
    }
    else 
    {  
      // decide how to descend to children
      
      PQP_REAL sz1 = o1->child(min_test.b1)->GetSize();
      PQP_REAL sz2 = o2->child(min_test.b2)->GetSize();

      res->num_bv_tests += 2;
 
      BVT bvt1,bvt2;
      PQP_REAL Ttemp[3];

      if (l2 || (!l1 && (sz1 > sz2)))	
      {  
        // put new tests on queue consisting of min_test.b2 
        // with children of min_test.b1 
      
        int c1 = o1->child(min_test.b1)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R,o1->child(c1)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c1)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c1)->To);
#endif
        MTxV(bvt1.T,o1->child(c1)->R,Ttemp);
		PQP_REAL S[3];
		
        bvt1.d = C2A_BV_Distance(bvt1.R,bvt1.T,
                            (C2A_BV *)o1->child(bvt1.b1),(C2A_BV *)o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = c2;
        bvt2.b2 = min_test.b2;
        MTxM(bvt2.R,o1->child(c2)->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,o1->child(c2)->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c2)->To);
#endif
        MTxV(bvt2.T,o1->child(c2)->R,Ttemp);
		
		
        bvt2.d = C2A_BV_Distance(bvt2.R,bvt2.T,
                            (C2A_BV *)o1->child(bvt2.b1),(C2A_BV *)o2->child(bvt2.b2),S);
      }
      else 
      {
        // put new tests on queue consisting of min_test.b1 
        // with children of min_test.b2
      
        int c1 = o2->child(min_test.b2)->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R,min_test.R,o2->child(c1)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->Tr,min_test.T);
#else
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->To,min_test.T);
#endif
		PQP_REAL S[3];
		
        bvt1.d = C2A_BV_Distance(bvt1.R,bvt1.T,
                            (C2A_BV *)o1->child(bvt1.b1), (C2A_BV *)o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R,min_test.R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->Tr,min_test.T);
#else
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->To,min_test.T);
#endif
		
		
        bvt2.d = C2A_BV_Distance(bvt2.R,bvt2.T,
                            (C2A_BV *)o1->child(bvt2.b1),(C2A_BV *)o2->child(bvt2.b2),S);
      }

      bvtq.AddTest(bvt1);	
      bvtq.AddTest(bvt2);
    }

    if (bvtq.Empty())
    {
      break;
    }
    else
    {
      min_test = bvtq.ExtractMinTest();

      if ((min_test.d + res->abs_err >= res->ctolerance) && 
         ((min_test.d * (1 + res->rel_err)) >= res->ctolerance)) 
      {
        break;
      }
    }
  }  
}

//
static void DeleteDuplicatePairs(C2A_ContactResult *res, 
							PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
							PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2)
{
	int iDistNum = 0;
	int iPairsNum = res->num_cpairs;
	int *pPairIndex = new int[iPairsNum];
	
	int i,j;
	bool bDistinct;
	int maxPair = 100;

	if(iPairsNum > maxPair)
	{
		// When the parameter for the contact query is very large, the number iPairsNum could be huge: for instance 342000.
		// At this time, just truncate it. 
		static int iWarning = 0;
		if(iWarning < 10)
		{
			iWarning++;
			printf("Warning %d: Wrong with the paramters of contact query. The number of paris is %d\n", iWarning, iPairsNum);
		}

		iPairsNum = maxPair;
	}

	for(i = iPairsNum-1; i >=0; i--)
	{
		//the start of first iterative
		bDistinct = true;
		for(j = i - 1; j >= 0; j--)
		{
			//the start of second iterative
			//go on if the distance and type of the two pairs are the same
			if (fabs(res->cpairs[i].distance - res->cpairs[j].distance) < 1e-10) 
			{
				if ((res->cpairs[i].f1.type == res->cpairs[j].f1.type) &&
					(res->cpairs[i].f2.type == res->cpairs[j].f2.type)) 
				{
					// compare the case of V/F, F/V, E/E
					//1 the case of V/F
					if (res->cpairs[i].f1.type == 0) 
					{
						// compare whether the vertex and face are the same
						int iVert, jVert, iVertTri, jVertTri, iTri, jTri;
						iVert = res->cpairs[i].f1.fid[0];
						jVert = res->cpairs[j].f1.fid[0];
						iVertTri = res->cpairs[i].f1.tri;
						jVertTri = res->cpairs[j].f1.tri;

						PQP_REAL tri1[3][3], tri2[3][3];						
						VcV(tri1[0], o1->trisConst[iVertTri].p1);
						VcV(tri1[1], o1->trisConst[iVertTri].p2);
						VcV(tri1[2], o1->trisConst[iVertTri].p3);

						VcV(tri2[0], o1->trisConst[jVertTri].p1);
						VcV(tri2[1], o1->trisConst[jVertTri].p2);
						VcV(tri2[2], o1->trisConst[jVertTri].p3);
						
						iTri = res->cpairs[i].f2.tri;
						jTri = res->cpairs[j].f2.tri;
						if ((iTri == jTri) && (tri1[iVert][0] == tri2[jVert][0]) &&
							(tri1[iVert][1] == tri2[jVert][1] && (tri1[iVert][2] == tri2[jVert][2]))) 
						{
							bDistinct = false;
							break;
						}
					}

					//2 the case of F/V
					if (res->cpairs[i].f1.type == 2) 
					{
					    // compare whether the face and vertex are the same
						int iVert, jVert, iVertTri, jVertTri, iTri, jTri;
						iVert = res->cpairs[i].f2.fid[0];
						jVert = res->cpairs[j].f2.fid[0];
						iVertTri = res->cpairs[i].f2.tri;
						jVertTri = res->cpairs[j].f2.tri;
						
						PQP_REAL tri1[3][3], tri2[3][3];						
						VcV(tri1[0], o2->trisConst[iVertTri].p1);
						VcV(tri1[1], o2->trisConst[iVertTri].p2);
						VcV(tri1[2], o2->trisConst[iVertTri].p3);
						
						VcV(tri2[0], o2->trisConst[jVertTri].p1);
						VcV(tri2[1], o2->trisConst[jVertTri].p2);
						VcV(tri2[2], o2->trisConst[jVertTri].p3);
						
						iTri = res->cpairs[i].f1.tri;
						jTri = res->cpairs[j].f1.tri;
						if ((iTri == jTri) && (tri1[iVert][0] == tri2[jVert][0]) &&
							(tri1[iVert][1] == tri2[jVert][1] && (tri1[iVert][2] == tri2[jVert][2]))) 
						{
							bDistinct = false;
							break;
						}
						
					}

					//3 the case of E/E
					if (res->cpairs[i].f1.type == 1) 
					{
						// compare whether the edge and edge are the same
						int iEdge1, iEdge2, iTri1, iTri2, jEdge1, jEdge2, jTri1, jTri2;
						iEdge1 = res->cpairs[i].f1.fid[0];
						iEdge2 = res->cpairs[i].f2.fid[0];
						iTri1 = res->cpairs[i].f1.tri;
						iTri2 = res->cpairs[i].f2.tri;

						jEdge1 = res->cpairs[j].f1.fid[0];
						jEdge2 = res->cpairs[j].f2.fid[0];												
						jTri1 = res->cpairs[j].f1.tri;
						jTri2 = res->cpairs[j].f2.tri;
						
						PQP_REAL iVtri1[3][3], iVtri2[3][3];						
						PQP_REAL jVtri1[3][3], jVtri2[3][3];		
						
						VcV(iVtri1[0], o1->trisConst[iTri1].p1);
						VcV(iVtri1[1], o1->trisConst[iTri1].p2);
						VcV(iVtri1[2], o1->trisConst[iTri1].p3);
						
						VcV(iVtri2[0], o2->trisConst[iTri2].p1);
						VcV(iVtri2[1], o2->trisConst[iTri2].p2);
						VcV(iVtri2[2], o2->trisConst[iTri2].p3);

						VcV(jVtri1[0], o1->trisConst[jTri1].p1);
						VcV(jVtri1[1], o1->trisConst[jTri1].p2);
						VcV(jVtri1[2], o1->trisConst[jTri1].p3);
						
						VcV(jVtri2[0], o2->trisConst[jTri2].p1);
						VcV(jVtri2[1], o2->trisConst[jTri2].p2);
						VcV(jVtri2[2], o2->trisConst[jTri2].p3);

						
						//compare the edge on model 1
						if ((((iVtri1[iEdge1][0] == jVtri1[jEdge1][0]) && (iVtri1[iEdge1][1] == jVtri1[jEdge1][1]) && (iVtri1[iEdge1][2] == jVtri1[jEdge1][2])) || 
											((iVtri1[(iEdge1+1)%3][0] == jVtri1[(jEdge1+1)%3][0]) && (iVtri1[(iEdge1+1)%3][1] == jVtri1[(jEdge1+1)%3][1]) && (iVtri1[(iEdge1+1)%3][2] == jVtri1[(jEdge1+1)%3][2]))) ||
											(((iVtri1[iEdge1][0] == jVtri1[(jEdge1+1)%3][0]) && (iVtri1[iEdge1][1] == jVtri1[(jEdge1+1)%3][1]) && (iVtri1[iEdge1][2] == jVtri1[(jEdge1+1)%3][2])) || 
											((iVtri1[(iEdge1+1)%3][0] == jVtri1[jEdge1][0]) && (iVtri1[(iEdge1+1)%3][1] == jVtri1[jEdge1][1]) && (iVtri1[(iEdge1+1)%3][2] == jVtri1[jEdge1][2]))))
								
						{
							//compare the edge on model 2
							if ((((iVtri2[iEdge2][0] == jVtri2[jEdge2][0]) && (iVtri2[iEdge2][1] == jVtri2[jEdge2][1]) && (iVtri2[iEdge2][2] == jVtri2[jEdge2][2])) || 
								((iVtri2[(iEdge2+1)%3][0] == jVtri2[(jEdge2+1)%3][0]) && (iVtri2[(iEdge2+1)%3][1] == jVtri2[(jEdge2+1)%3][1]) && (iVtri2[(iEdge2+1)%3][2] == jVtri2[(jEdge2+1)%3][2]))) ||
								(((iVtri2[iEdge2][0] == jVtri2[(jEdge2+1)%3][0]) && (iVtri2[iEdge2][1] == jVtri2[(jEdge2+1)%3][1]) && (iVtri2[iEdge2][2] == jVtri2[(jEdge2+1)%3][2])) || 
								((iVtri2[(iEdge2+1)%3][0] == jVtri2[jEdge2][0]) && (iVtri2[(iEdge2+1)%3][1] == jVtri2[jEdge2][1]) && (iVtri2[(iEdge2+1)%3][2] == jVtri2[jEdge2][2])))) 
							{
								bDistinct = false;
								break;
							}
							
						}

					}
				}
			}
						
		}
		//the end of second iterative

		if (bDistinct) 
		{
			pPairIndex[iDistNum++] = i;
		}
	}
	//the end of first iterative
	if (iDistNum) 
	{
		C2A_ContactPair * newpairs = new C2A_ContactPair[iDistNum];
		for(int k = 0; k < iDistNum; k++)
		{
			newpairs[k].f1 = res->cpairs[pPairIndex[k]].f1;
			newpairs[k].f2 = res->cpairs[pPairIndex[k]].f2;
			newpairs[k].distance = res->cpairs[pPairIndex[k]].distance;
			VcV(newpairs[k].p1, res->cpairs[pPairIndex[k]].p1);
			VcV(newpairs[k].p2, res->cpairs[pPairIndex[k]].p2);
			
		}
		delete [] res->cpairs;
		res->cpairs = newpairs;
		res->num_cpairs = iDistNum;
		res->num_cpairs_alloced = iDistNum;
		
	}
	delete []pPairIndex;

	if(bPQPContact_Represenative)
	{
		//report the representative pairs by distance clustering
		int repNum = 0;
		int *pRepIndex = new int[iDistNum];
		pRepIndex[0] = 0;
		repNum++;
		PQP_REAL p1[3], p2[3], p3[3], p4[3];
		PQP_REAL square_epislon = res->ctolerance * res->ctolerance * PQPContact_radius * PQPContact_radius;

		//cluster the pairs
		for( i = 1; i < iDistNum; i++)
		{

			VcV(p1, res->cpairs[i].p1);
			VcV(p2, res->cpairs[i].p2);

			bool bSame = false;
			for(int j = 0; j < repNum; j++)
			{
				VcV(p3, res->cpairs[pRepIndex[j]].p1);
				VcV(p4, res->cpairs[pRepIndex[j]].p2);

				// if ((VdistV2(p1, p3) <= epislon ) || (VdistV2(p2, p4) <= epislon )) 
				if (  (VdistV2(p1, p3) <= square_epislon ) && (VdistV2(p2, p4) <= square_epislon ) )
				{
	//				pRepIndex[repNum++] = i;
	//				break;
					bSame = true;
					double d1 = res->cpairs[i].distance;
					double d2 = res->cpairs[pRepIndex[j]].distance;

					if (d1 < d2)
					{
						pRepIndex[j] = i;
					}

					break;
				}
			}

			if (!bSame) 
			{
				pRepIndex[repNum++] = i;			
			}

		}

		//copy the representative pairs
		C2A_ContactPair * reppairs = new C2A_ContactPair[repNum];	
		for(i = 0; i < repNum; i++)
		{
			reppairs[i].f1 = res->cpairs[pRepIndex[i]].f1;
			reppairs[i].f2 = res->cpairs[pRepIndex[i]].f2;
			reppairs[i].distance = res->cpairs[pRepIndex[i]].distance;
			VcV(reppairs[i].p1, res->cpairs[pRepIndex[i]].p1);
			VcV(reppairs[i].p2, res->cpairs[pRepIndex[i]].p2);
		}

		delete [] res->cpairs;
		res->cpairs = reppairs;
		res->num_cpairs = repNum;
		res->num_cpairs_alloced = repNum;
		delete []pRepIndex;
	}
}

void
ContactRecurse(C2A_ContactResult *res,
                PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
                C2A_Model *o1, int b1,
                C2A_Model *o2, int b2)
{
	//return if the models are colliding
  if (res->modelcollided) 
  {
	  return;
  }
  PQP_REAL sz1 = o1->child(b1)->GetSize();
  PQP_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.

    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    Tri *t1 = o1->GetTriangle(-o1->child(b1)->first_child - 1);
    Tri *t2 = o2->GetTriangle(-o2->child(b2)->first_child - 1);

	C2A_ContactFeature f1, f2;	
    PQP_REAL d = C2A_TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,res->modelcollided);
  
	if (res->modelcollided) 
	{
		return;
	}
	
	if (d <= res->ctolerance) 
	{
		//add the contact feature if < ctolerance
		f1.tri = t1->id;
		f2.tri = t2->id;
		//huangxin 071126
		f1.bid = b1;
		f2.bid = b2;
		
		res->Add(f1, f2, d, p, q);
		//res->closer_than_ctolerance = true;
		
		VcV(res->p1, p);         // p already in c.s. 1
		VcV(res->p2, q);         // q must be transformed 

		o1->last_tri = t1;
		o2->last_tri = t2;		
	}


    return;
  }

  // First, perform distance tests on the children. Then traverse 
  // them recursively, but test the closer pair first, the further 
  // pair second.

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child+1;
    c2 = b2;
    
    MTxM(R1,o1->child(a1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(c1)->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,o1->child(c1)->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child+1;

    MxM(R1,R,o2->child(a2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,o2->child(a2)->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif

    MxM(R2,R,o2->child(c2)->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,o2->child(c2)->Tr,T);
#else
    MxVpV(T2,R,o2->child(c2)->To,T);
#endif
  }

  res->num_bv_tests += 2;
  PQP_REAL S[3];
  
  PQP_REAL d1 = C2A_BV_Distance(R1, T1, (C2A_BV *)o1->child(a1), (C2A_BV *)o2->child(a2),S);
  PQP_REAL d2 = C2A_BV_Distance(R2, T2, (C2A_BV *)o1->child(c1), (C2A_BV *)o2->child(c2),S);

  if (d2 < d1)
  {
    if ((d2 < (res->ctolerance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
	
	if (res->modelcollided) 
	{
		return;
	}

    if ((d1 < (res->ctolerance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else 
  {
    if ((d1 < (res->ctolerance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R1, T1, o1, a1, o2, a2);
    }

	if (res->modelcollided) 
	{
		return;
	}

    if ((d2 < (res->ctolerance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->ctolerance)) 
    {      
      ContactRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
  }
}



C2A_ContactResult::C2A_ContactResult()
{
	cpairs = 0;
	num_cpairs = num_cpairs_alloced = 0;
}

C2A_ContactResult::~C2A_ContactResult()
{
	delete [] cpairs;
}

void
C2A_ContactResult::FreePairsList()
{
	num_cpairs = num_cpairs_alloced = 0;
	delete [] cpairs;
	cpairs = 0;
}
// may increase OR reduce mem usage
void
C2A_ContactResult::SizeTo(int n)
{
	C2A_ContactPair *temp;
	
	if (n < num_cpairs) 
	{
		fprintf(stderr, "PQP Error: Internal error in "
			"'C2A_ContactResult::SizeTo(int n)'\n");
		fprintf(stderr, "       n = %d, but num_cpairs = %d\n", n, num_cpairs);
		return;
	}
	
	temp = new C2A_ContactPair[n];
	memcpy(temp, cpairs, num_cpairs*sizeof(C2A_ContactPair));
	delete [] cpairs;
	cpairs = temp;
	num_cpairs_alloced = n;
	return;
}

void C2A_ContactResult::Add(const C2A_ContactFeature& f1, const C2A_ContactFeature& f2, PQP_REAL distance, PQP_REAL p[], PQP_REAL q[])
{
	// Check whether a feature pair is duplicated or not
	bool bDuplicated = false;
	if(bPQPContact_Represenative)
	{
		// The tolerance $\kappa_r$
		PQP_REAL square_epislon = ctolerance * ctolerance * PQPContact_radius * PQPContact_radius;

		//report the representative pairs by distance clustering
		for( int i = 0; i < num_cpairs; i++)
		{
			// It is ok to use the coornidates relative to R1, T1. 
			// MxVpV(p1, R1, res->cpairs[i].p1, T1);
			// MxVpV(p2, R1, res->cpairs[i].p2, T1);

			// if ((VdistV2(p1, p3) <= epislon ) || (VdistV2(p2, p4) <= epislon )) 
			if (  (VdistV2(p, cpairs[i].p1) <= square_epislon ) 
			   && (VdistV2(q, cpairs[i].p2) <= square_epislon ) )
			{
				bDuplicated = true;

				// replace the item by the new item
				if (cpairs[i].distance > distance)
				{
					cpairs[i].f1 = f1;
					cpairs[i].f2 = f2;
					cpairs[i].distance = distance;
					VcV(cpairs[i].p1, p);
					VcV(cpairs[i].p2, q);
				}
			}
		}
	}

	// now proceed as usual
	if(!bDuplicated)
	{
		if (num_cpairs >= num_cpairs_alloced) 
		{
			// allocate more
			
			SizeTo(num_cpairs_alloced*2+8);
		}
		cpairs[num_cpairs].f1 = f1;
		cpairs[num_cpairs].f2 = f2;
		cpairs[num_cpairs].distance = distance;
		VcV(cpairs[num_cpairs].p1, p);
		VcV(cpairs[num_cpairs].p2, q);
	
		num_cpairs++;
	}
}

// toc
double TOCStepRecurse_Dis(PQP_REAL tolerance_t,
						  CInterpMotion *objmotion1, 
						  CInterpMotion *objmotion2,
						  C2A_TimeOfContactResult *res, 
						  double R[3][3], 
						  double T[3], 
						  C2A_Model *o1, 
						  int b1,
						  C2A_Model *o2, 
						  int b2,
						  double delta)
{
	//delta =0.0001;
	int l1 = o1->child(b1)->Leaf();
	int l2 = o2->child(b2)->Leaf();


	PQP_REAL r1[3][3], r2[3][3], tt1[3], tt2[3], result1[3], result2[3], temp1[3], temp2[3],cwc[3];

	(objmotion1->transform).Rotation().Get_Value(r1);
	(objmotion2->transform).Rotation().Get_Value(r2);

	(objmotion1->transform).Translation().Get_Value(tt1);
	(objmotion2->transform).Translation().Get_Value(tt2);

	objmotion1->m_axis.Get_Value(cwc);

	if (l1 && l2)
	{ 
		double p[3], q[3];

		Tri *t1 = o1->GetTriangle(-o1->child(b1)->first_child - 1);
		Tri *t2 = o2->GetTriangle(-o2->child(b2)->first_child - 1);

		PQP_REAL dTri = TriDistance(res->R,res->T,t1,t2,p,q);

		if (dTri <= res->distance)
		{
			res->distance=dTri;
			double S1[3],S2[3];
			MxV(temp1, r1, p);
			VpV(result1, temp1, tt1);
			MxV(temp2, r1, q);
			VpV(result2, temp2, tt1);
			VmV(S1,result2,result1);
			VxS(S2,S1,-1);
			VcV(res->p1,p);
			VcV(res->p2,q); 

			PQP_REAL motionbound1 = objmotion1->computeTOC(dTri, ((C2A_BV*) (o1->child(b1)))->angularRadius, S1); //  
			PQP_REAL motionbound2 = objmotion2->computeTOC(dTri, ((C2A_BV*) (o2->child(b2)))->angularRadius, S2);

			PQP_REAL mint =  (dTri)/(motionbound1+motionbound2);
			if (mint<0.0)
			{
				mint=0.0;
			}
			if (mint <= res->mint) 
				res->mint = mint;

			res->distance=dTri;	
			o1->last_tri = t1;
			o2->last_tri = t2;


		}
		res->num_tri_tests++;

		return dTri;  
	}

	// First, perform distance tests on the children. Then traverse
	// them recursively, but test the closer pair first, the further
	// pair second.

	int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
	double R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

	double sz1 = o1->child(b1)->GetSize();
	double sz2 = o2->child(b2)->GetSize();

	if (l2 || (!l1 && (sz1 > sz2)))
	{
		// visit the children of b1

		a1 = o1->child(b1)->first_child;
		a2 = b2;
		c1 = o1->child(b1)->first_child+1;
		c2 = b2;

		MTxM(R1,o1->child(a1)->R,R);
		VmV(Ttemp,T,o1->child(a1)->Tr);
		MTxV(T1,o1->child(a1)->R,Ttemp);

		MTxM(R2,o1->child(c1)->R,R);
		VmV(Ttemp,T,o1->child(c1)->Tr);
		MTxV(T2,o1->child(c1)->R,Ttemp);
	}
	else
	{
		// visit the children of b2

		a1 = b1;
		a2 = o2->child(b2)->first_child;
		c1 = b1;
		c2 = o2->child(b2)->first_child+1;

		MxM(R1,R,o2->child(a2)->R);
		MxVpV(T1,R,o2->child(a2)->Tr,T);

		MxM(R2,R,o2->child(c2)->R);
		MxVpV(T2,R,o2->child(c2)->Tr,T);
	}

	PQP_REAL a1_clo[3], a2_clo[3];
	PQP_REAL c1_clo[3], c2_clo[3];
	bool bValid_A_clo;
	bool bValid_C_clo;
	PQP_REAL S1[3],S2[3];
	PQP_REAL motion_bound1, motion_bound2;
	PQP_REAL minta,mintb;


	double d1 = C2A_BV_Distance(R1, T1, (C2A_BV *)o1->child(a1), (C2A_BV *)o2->child(a2), a1_clo, a2_clo, bValid_A_clo,S1);
	if (d1!=0.0)
	{
		MxV(temp1, ((C2A_BV*)o1->child(a1))->R_loc, S1);
		MxV(S1, r1, temp1);
		VxS(S2, S1, -1);


		motion_bound1 = objmotion1->computeTOC_MotionBound(T1, d1, ((C2A_BV*)o1->child(a1)), S1); 
		motion_bound2 = objmotion2->computeTOC_MotionBound(T1, d1, ((C2A_BV*)o2->child(a2)), S2); 
		mintb         = (d1 )/(motion_bound1+motion_bound2);
		if (mintb <= 0)
		{
			mintb = 0.0;
		}

	}
	else
		mintb = 0.0;


	double d2 = C2A_BV_Distance(R2, T2, (C2A_BV *)o1->child(c1), (C2A_BV *)o2->child(c2), c1_clo, c2_clo, bValid_C_clo,S1);
	if (d2!=0.0)
	{
		MxV(temp1, ((C2A_BV*)o1->child(c1))->R_loc, S1);
		MxV(S1, r1, temp1);
		VxS(S2, S1, -1);

		motion_bound1 = objmotion1->computeTOC_MotionBound(T2, d2, ((C2A_BV*)o1->child(c1)), S1); 
		motion_bound2 = objmotion2->computeTOC_MotionBound(T2, d2, ((C2A_BV*)o2->child(c2)), S2); 

		minta         = (d2)/(motion_bound1+motion_bound2);

		if (minta <= 0)
		{
			minta = 0.0;
		}
	}
	else
		minta = 0.0;


	res->num_bv_tests+=2;

	if (d2 < d1)
	{	

		if (minta<res->UpboundTOC && ((d2 < (res->distance - res->abs_err)) || 
			(d2*(1 + res->rel_err) < res->distance))) 
		{  
			TOCStepRecurse_Dis(tolerance_t,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);  

		}
		else
		{

			if (minta < res->mint)	
			{
				res->mint=minta;
			}		 
		}

		if (mintb<res->UpboundTOC &&((d1 < (res->distance - res->abs_err)) || 
			(d1*(1 + res->rel_err) < res->distance))) 
		{  
			TOCStepRecurse_Dis(tolerance_t,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
		}
		else
		{

			if (mintb < res->mint)
			{     
				res->mint=mintb;
			}

		}

	}
	else 
	{
		if (mintb<res->UpboundTOC &&((d1 < (res->distance - res->abs_err)) || 
			(d1*(1 + res->rel_err) < res->distance))) 
		{ 

			TOCStepRecurse_Dis(tolerance_t,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
		}
		else
		{	

			if (mintb < res->mint)
			{     
				res->mint=mintb;
			}
		}


		if (minta<res->UpboundTOC &&((d2 < (res->distance - res->abs_err)) || 
			(d2*(1 + res->rel_err) < res->distance))) 
		{	 


			TOCStepRecurse_Dis(tolerance_t,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);
		}
		else
		{

			if (minta < res->mint)
			{     
				res->mint=minta;
			}

		}


	}
	return PQP_OK;

}

//*********************************************************************************************************************//
//  CA_Translation                                                                                                     //
//*********************************************************************************************************************//



double 
TOCStepRecurse_Dis_Translation(PQP_REAL tolerance_t, 
							   CInterpMotion *objmotion1, 
							   CInterpMotion *objmotion2,
							   C2A_TimeOfContactResult *res, 
							   double R[3][3], 
							   double T[3], 
							   C2A_Model *o1, 
							   int b1,
							   C2A_Model *o2, 
							   int b2,
							   double delta)
{

	int l1 = o1->child(b1)->Leaf();
	int l2 = o2->child(b2)->Leaf();


	PQP_REAL r1[3][3], r2[3][3], tt1[3], tt2[3], cwc[3];

	(objmotion1->transform).Rotation().Get_Value(r1);
	(objmotion2->transform).Rotation().Get_Value(r2);

	(objmotion1->transform).Translation().Get_Value(tt1);
	(objmotion2->transform).Translation().Get_Value(tt2);

	objmotion1->m_axis.Get_Value(cwc);

	if (l1 && l2)
	{ 
		PQP_REAL tri1[3][3], tri2[3][3], p[3], q[3];

		Tri *t1 = o1->GetTriangle(-o1->child(b1)->first_child - 1);
		Tri *t2 = o2->GetTriangle(-o2->child(b2)->first_child - 1);

		VcV(tri1[0], t1->p1);
		VcV(tri1[1], t1->p2);
		VcV(tri1[2], t1->p3);
		MxVpV(tri2[0], res->R, t2->p1, res->T);
		MxVpV(tri2[1], res->R, t2->p2, res->T);
		MxVpV(tri2[2], res->R, t2->p3, res->T);

		bool CAtri_flag=objmotion1->CAonNonAdjacentTriangles(r1, tt1, 
			tri1,tri2,&res->mint, &res->distance);

		
		if (CAtri_flag)
		{			
			
			VcV(res->p1,p);
			VcV(res->p2,q); 	
			res->last_triA = t1;
			res->last_triB = t2;

		}
		res->num_tri_tests++;		

		return res->distance;  
	}

	// First, perform distance tests on the children. Then traverse
	// them recursively, but test the closer pair first, the further
	// pair second.

	int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
	double R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

	double sz1 = o1->child(b1)->GetSize();
	double sz2 = o2->child(b2)->GetSize();

	if (l2 || (!l1 && (sz1 > sz2)))
	{
		// visit the children of b1

		a1 = o1->child(b1)->first_child;
		a2 = b2;
		c1 = o1->child(b1)->first_child+1;
		c2 = b2;

		MTxM(R1,o1->child(a1)->R,R);
		VmV(Ttemp,T,o1->child(a1)->Tr);
		MTxV(T1,o1->child(a1)->R,Ttemp);

		MTxM(R2,o1->child(c1)->R,R);
		VmV(Ttemp,T,o1->child(c1)->Tr);
		MTxV(T2,o1->child(c1)->R,Ttemp);
	}
	else
	{
		// visit the children of b2

		a1 = b1;
		a2 = o2->child(b2)->first_child;
		c1 = b1;
		c2 = o2->child(b2)->first_child+1;

		MxM(R1,R,o2->child(a2)->R);
		MxVpV(T1,R,o2->child(a2)->Tr,T);

		MxM(R2,R,o2->child(c2)->R);
		MxVpV(T2,R,o2->child(c2)->Tr,T);
	}


	double d1, d2;
	PQP_REAL minta,mintc;


    d1=d2=1e+30;
	minta=mintc=res->mint;
	bool CARSS_flag_a=objmotion1->CAonRSS(r1, tt1, r2, tt2, R1, T1, (C2A_BV *)o1->child(a1), (C2A_BV *)o2->child(a2),  &minta, &d1);

	bool CARSS_flag_c=objmotion1->CAonRSS(r1, tt1, r2, tt2, R2, T2, (C2A_BV *)o1->child(c1), (C2A_BV *)o2->child(c2),  &mintc, &d2);

	res->num_bv_tests+=2;

	if (d2 < d1)
	{	

		if (mintc<res->mint && ((d2 < (res->distance - res->abs_err)) || 
			(d2*(1 + res->rel_err) < res->distance))) 
		{  
			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);  

		}
		



		if (minta<res->mint &&((d1 < (res->distance - res->abs_err)) || 
			(d1*(1 + res->rel_err) < res->distance))) 
		{  
			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
		}
		

	}
	else 
	{
		if (minta<res->mint &&((d1 < (res->distance - res->abs_err)) || 
			(d1*(1 + res->rel_err) < res->distance))) 
		{ 

			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R1, T1, o1, a1, o2, a2, delta);
		}
		

		if (mintc<res->mint &&((d2 < (res->distance - res->abs_err)) || 
			(d2*(1 + res->rel_err) < res->distance))) 
		{	 


			TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R2, T2, o1, c1, o2, c2, delta);
		}		


	}
	return PQP_OK;

}



void
TOCStepRecurse_Dis_contact(C2A_TimeOfContactResult *res,
						   PQP_REAL R[3][3], 
						   PQP_REAL T[3], // b2 relative to b1
						   C2A_Model *o1, 
						   int b1,
						   C2A_Model *o2, 
						   int b2,
						   int depth)
{

	//if (res->cont_l.size()>30) {
		//printf("+++++++++++++++++++++++depth %d num contact %d, feature size %d\n", depth, res->num_contact, res->cont_l.size());
		//res->num_contact = 0;
		//return;
	//}

  PQP_REAL sz1 = o1->child(b1)->GetSize();
  PQP_REAL sz2 = o2->child(b2)->GetSize();
  int l1 = o1->child(b1)->Leaf();
  int l2 = o2->child(b2)->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.



    PQP_REAL p[3], q[3];

	int trnum1 = -o1->child(b1)->first_child - 1;
	int trnum2 = -o2->child(b2)->first_child - 1;
    C2A_Tri *t1 = o1->GetTriangle(trnum1);//(C2A_Tri*)&((C2A_Tri*)o1->tris)[-o1->child(b1)->first_child - 1];
    C2A_Tri *t2 = o2->GetTriangle(trnum2);//(C2A_Tri*)&((C2A_Tri*)o2->tris)[-o2->child(b2)->first_child - 1];
	C2A_ContactFeature f1,f2;
	bool bCollided=false;
    PQP_REAL d = C2A_TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,bCollided) ;
  
	if (f1.type == -1 || f2.type == -1) return;
    if (d <= res->distance) 
    {
	
			int FeatureType_A = f1.type+1;//1-vertex;2-Edge;3-Face;
			int FeatureType_B = f2.type+1;//1-vertex;2-Edge;3-Face;
			int FeatureID_A[3];
			int FeatureID_B[3];
			
			if (f1.type==0)
			{
				FeatureID_A[0]=t1->Index()[f1.fid[0]];
				
			}
			if(f1.type==1)
			{
				FeatureID_A[0]=t1->Index()[f1.fid[0]];
				FeatureID_A[1]=t1->Index()[(f1.fid[0]+1)%3];		
				
			}
			if(f1.type==2)
			{
				FeatureID_A[0]=t1->Index()[0];
				FeatureID_A[1]=t1->Index()[1];	
				FeatureID_A[2]=t1->Index()[2];		
				
			}
			
			
			if(f2.type==0)
			{
				FeatureID_B[0]=t2->Index()[f2.fid[0]];
				
			}
			if(f2.type==1)
			{
				FeatureID_B[0]=t2->Index()[f2.fid[0]];
				FeatureID_B[1]=t2->Index()[(f2.fid[0]+1)%3];		
				
			}
			if(f2.type==2)
			{
				FeatureID_B[0]=t2->Index()[0];
				FeatureID_B[1]=t2->Index()[1];	
				FeatureID_B[2]=t2->Index()[2];		
				
			}
	
			
			int TriangleID_A = -o1->child(b1)->first_child - 1;
			int TriangleID_B = -o2->child(b2)->first_child - 1;
			
			PQP_REAL P_A[3];
			VcV(P_A,p);
			
			PQP_REAL P_B[3];
			PQP_REAL temp[3];
			VmV(temp,q,res->T);
			MTxV(P_B,res->R,temp);
			
			
			
			
			PQP_REAL Distance=d;
			
			res->cont_l.push_front(ContactF(FeatureType_A, FeatureType_B, FeatureID_A, FeatureID_B,
				TriangleID_A, TriangleID_B,  P_A,  P_B, Distance));
			
			res->num_contact++;

    }

    return;
  }

  // First, perform distance tests on the children. Then traverse 
  // them recursively, but test the closer pair first, the further 
  // pair second.

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = o1->child(b1)->first_child;
    a2 = b2;
    c1 = o1->child(b1)->first_child+1;
    c2 = b2;
    
    MTxM(R1,o1->child(a1)->R,R);

    VmV(Ttemp,T,o1->child(a1)->Tr);


    MTxV(T1,o1->child(a1)->R,Ttemp);

    MTxM(R2,o1->child(c1)->R,R);

    VmV(Ttemp,T,o1->child(c1)->Tr);


    MTxV(T2,o1->child(c1)->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = o2->child(b2)->first_child;
    c1 = b1;
    c2 = o2->child(b2)->first_child+1;

    MxM(R1,R,o2->child(a2)->R);

    MxVpV(T1,R,o2->child(a2)->Tr,T);



    MxM(R2,R,o2->child(c2)->R);

    MxVpV(T2,R,o2->child(c2)->Tr,T);

    

  }

  PQP_REAL S[3];
  

  PQP_REAL d1 = C2A_BV_Distance(R1, T1, (C2A_BV*) o1->child(a1), (C2A_BV*) o2->child(a2), S);
  PQP_REAL d2 = C2A_BV_Distance(R2, T2, (C2A_BV*) o1->child(c1), (C2A_BV*) o2->child(c2), S);

  if (d2 < d1)
  {
    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1.0 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R2, T2, o1, c1, o2, c2, depth+1);      
    }

    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1.0 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R1, T1, o1, a1, o2, a2, depth+1);
    }
  }
  else 
  {
    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1.0 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R1, T1, o1, a1, o2, a2, depth+1);
    }

    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1.0 + res->rel_err) < res->distance)) 
    {      
      TOCStepRecurse_Dis_contact(res, R2, T2, o1, c1, o2, c2, depth+1);      
    }
  }
}

int 
C2A_TimeOfContactStep_Contact(CInterpMotion *objmotion1, 
					CInterpMotion *objmotion2, 
					C2A_TimeOfContactResult *res,
					PQP_REAL R1[3][3], 
					PQP_REAL T1[3], 
					C2A_Model *o1,
					PQP_REAL R2[3][3],
					PQP_REAL T2[3], 
					C2A_Model *o2,
					double threshold)
{


	MTxM(res->R,R1,R2);
	double Ttemp[3];
	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	
	double Rtemp[3][3], R[3][3], T[3];
	
	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);
	
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
	MTxV(T,o1->child(0)->R,Ttemp);

	PQP_REAL dTri;
	dTri=res->distance;
	res->distance= threshold;
		
    res->rel_err=0;
	res->abs_err=0;
	
	   
	res->mint = 1;


    TOCStepRecurse_Dis_contact(res, R,T,o1,0,o2,0,0);

	res->distance=dTri;

	return PQP_OK;


	
}

// every CA iterative
int 
C2A_TimeOfContactStep(CInterpMotion *objmotion1, 
			CInterpMotion *objmotion2, 
			C2A_TimeOfContactResult *res,
			PQP_REAL R1[3][3], 
			PQP_REAL T1[3], 
			C2A_Model *o1,
			PQP_REAL R2[3][3], 
			PQP_REAL T2[3], 
			C2A_Model *o2,
			PQP_REAL tolerance_t,
			PQP_REAL tolerance_d)
{
   double Ttemp[3];
   
	MTxM(res->R,R1,R2);
	
	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	
	double Rtemp[3][3], R[3][3], T[3];
	
	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);
	
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
	MTxV(T,o1->child(0)->R,Ttemp);


	PQP_REAL p[3],q[3];

	Tri *t1, *t2;
	t1=res->last_triA;
	t2=res->last_triB;
	/*res->last_triA = t1;
	res->last_triB = t2;*/


	if (b_TanslationCCD)
	{
		PQP_REAL r1[3][3],tt1[3];
		(objmotion1->transform).Rotation().Get_Value(r1);	

		(objmotion1->transform).Translation().Get_Value(tt1);


		PQP_REAL dTri,mint;

		PQP_REAL tri1[3][3], tri2[3][3];


		VcV(tri1[0], t1->p1);
		VcV(tri1[1], t1->p2);
		VcV(tri1[2], t1->p3);
		MxVpV(tri2[0], res->R, t2->p1, res->T);
		MxVpV(tri2[1], res->R, t2->p2, res->T);
		MxVpV(tri2[2], res->R, t2->p3, res->T);
		mint=1.0;

		bool CAtri_flag=objmotion1->CAonNonAdjacentTriangles(r1, tt1, 
			tri1,tri2,&mint, &dTri);
		if (CAtri_flag)
		{
			res->mint=mint;
			res->distance=dTri;
		}
		else
		{
			res->mint=1.0;
			res->distance=1e+30;

		}
	}
	else
	{
		C2A_ContactFeature f1,f2;
		bool bCollided;

		PQP_REAL dTri =res->distance = C2A_TriDistance(res->R,res->T,t1,t2,p,q,f1,f2,bCollided) ;

		if (res->numCA==0)
		{ 		

			res->mint=1;
		}

		//traverse the C2A_BVH Tree to compute the toc
		//controlled

		if (res->mint<=0.005||res->distance<=0.5||res->numCA>5)
		{ 
			res->abs_err=0;
			res->rel_err=0;


		}
		else
		{

			res->abs_err=Max_Value;

			if(res->numCA<=2)
			{   
				res->rel_err=3;

			}
			else
			{
				res->rel_err=0.5;
			}   	

		}


		res->mint = 1; 
	
	}

	
   
	if (b_TanslationCCD)
	{
		res->abs_err=0;
		res->rel_err=0;

		TOCStepRecurse_Dis_Translation(tolerance_t,objmotion1, objmotion2, res, R,T,o1,0,o2,0,tolerance_d);
		t1=res->last_triA;
		t2=res->last_triB;
		objmotion1->integrate(res->mint, R1, T1);
		objmotion2->integrate(res->mint, R2, T2); 

		MTxM(res->R,R1,R2);

		VmV(Ttemp, T2, T1);  
		MTxV(res->T, R1, Ttemp);

		res->distance=TriDistance(res->R,res->T,t1,t2,p,q) ;
	}
	else
	{
		TOCStepRecurse_Dis(tolerance_t,objmotion1, objmotion2, res, R,T,o1,0,o2,0,tolerance_d);
	}


	



	
	return PQP_OK;
	
}




//get the contact features
PQP_REAL 
C2A_QueryContact(CInterpMotion *objmotion1, 
				 CInterpMotion *objmotion2,
				 C2A_TimeOfContactResult *res,  
				 C2A_Model *o1, 
				 C2A_Model *o2,
				 double threshold)
{
	PQP_REAL toc = 0;

	res->num_contact=0;
	res->UpboundTOC=1;
 
	PQP_REAL R1[3][3];
	PQP_REAL T1[3];
	PQP_REAL R2[3][3];
	PQP_REAL T2[3];

	objmotion1->transform.Rotation().Get_Value(R1);
 	objmotion1->transform.Translation().Get_Value(T1);
	objmotion2->transform.Rotation().Get_Value(R2);
	objmotion2->transform.Translation().Get_Value(T2);	
	
	C2A_TimeOfContactStep_Contact(objmotion1, objmotion2, res, R1,T1,o1,
			R2,T2,o2,threshold);

	return PQP_OK;
		
		
}

PQP_REAL
C2A_QueryContactOnly(C2A_TimeOfContactResult *res,  
					PQP_REAL R1[3][3], 
					PQP_REAL T1[3], 
					C2A_Model *o1,
					PQP_REAL R2[3][3],
					PQP_REAL T2[3], 
					C2A_Model *o2,
					double threshold)
{
	res->cont_l.clear();
	res->num_contact=0;
	res->UpboundTOC=1;
	C2A_TimeOfContactStep_Contact(0, 0, res, R1, T1, o1, R2, T2, o2, threshold);
	return PQP_OK;
}


//calculate the toc 
PQP_REAL 
C2A_QueryTimeOfContact(CInterpMotion *objmotion1, 
					   CInterpMotion *objmotion2,
					   C2A_TimeOfContactResult *res,  
					   C2A_Model *o1, 
					   C2A_Model *o2,
					   PQP_REAL tolerance_d, 
					   PQP_REAL tolerance_t, 
					   int qsize)
{
 

	int tmpi;
	int tmpj;

	PQP_REAL mint;
	PQP_REAL toc = 0;

	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	res->num_contact=0;
	res->UpboundTOC=1;
    res->numCA=0;
	PQP_REAL R1[3][3];
	PQP_REAL T1[3];
	PQP_REAL R2[3][3];
	PQP_REAL T2[3];

	objmotion1->transform.Rotation().Get_Value(R1);
 	objmotion1->transform.Translation().Get_Value(T1);
	objmotion2->transform.Rotation().Get_Value(R2);
	objmotion2->transform.Translation().Get_Value(T2);
    PQP_REAL lamda = 0.0, lastLamda=0, dlamda=0.0;
	PQP_REAL dist=0;
	int nItrs = 0;
	bool T_CA=false;
		
	C2A_TimeOfContactStep(objmotion1, objmotion2, res, R1,T1,o1,
		R2,T2,o2,tolerance_t,tolerance_d);

	nItrs=0;

	dist = res->distance; 
	mint = res->mint;

	if (b_TanslationCCD)
	{
		res->toc=mint;
		if (res->toc>=1.0)
		{
			res->collisionfree = true; 
		}
		else
		{
			res->collisionfree = false;

		}
		toc=mint;



		return toc;
	}
	

		
	res->numCA = 1;
		
		
	if(dist  <= tolerance_d  || mint <= tolerance_t)
	{
		res->collisionfree = false; 			
		res->toc = tolerance_t;
		for(tmpi=0; tmpi<3; tmpi++)
		{	 
			for( tmpj=0; tmpj<3; tmpj++)
			{
				res->R_toc[tmpi][tmpj] = R1[tmpi][tmpj];
			}
			res->T_toc[tmpi] = T1[tmpi];
		}
		if(toc>=1)
			toc=0;	

	}
	lastLamda=mint;


	while (dist > tolerance_d) //not close enough
	{
		nItrs++;

		if(nItrs > 150)
		{
			break;
		}
			// 1. If dist > path_max, no collision during t= [0, 1] 
			// 2. If dist < path_max, look for TOC
		if( mint >= 1.0)
		{
			res->collisionfree = true;
			res->toc = 0;
			toc = 0;

			return 0;
				
		}
		else
		{
			dlamda = mint;
			if (dlamda < tolerance_t)
			{												
				break;
			}	
			lamda += dlamda;   				
			if (lamda >= 1.0)
			{
				res->collisionfree = true;
				res->toc = 0;
				toc = 0;
				return 0;

			}
			lastLamda = lamda;				
			res->numCA++;   
			objmotion1->integrate(lamda, R1, T1);
			objmotion2->integrate(lamda, R2, T2); 
			res->UpboundTOC=1.0-lamda;				

			C2A_TimeOfContactStep(objmotion1, objmotion2, res, R1,T1,o1,
					R2,T2,o2,tolerance_t,tolerance_d);

			dist=res->distance;
			mint = res->mint; 

		}
	}

	if(dist == 0 && mint<0)
	{
		toc = lastLamda- dlamda ;
	}
	else
	{
		toc = lastLamda;
	}

	res->collisionfree = false;
			
	if(toc>=1-tolerance_t)
		{toc=0;}
	res->toc = toc;
	//compute the result R, T at TOC 
	objmotion1->integrate(toc, res->R_toc, res->T_toc);//make res->R_toc equals R1
			
		
	return toc;
		
		
}

//calculate the toc 
PQP_REAL 
C2A_QueryTimeOfContact_TranslateCA(CInterpMotion *objmotion1, 
								   CInterpMotion *objmotion2,
								   C2A_TimeOfContactResult *res,  
								   C2A_Model *o1, 
								   C2A_Model *o2,
								   PQP_REAL tolerance_d, 
								   PQP_REAL tolerance_t, 
								   int qsize)
{


	int tmpi;
	int tmpj;

	PQP_REAL mint;
	PQP_REAL toc = 0;

	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	res->num_contact=0;
	res->UpboundTOC=1;
	res->numCA=0;
	PQP_REAL R1[3][3];
	PQP_REAL T1[3];
	PQP_REAL R2[3][3];
	PQP_REAL T2[3];

	objmotion1->transform.Rotation().Get_Value(R1);
	objmotion1->transform.Translation().Get_Value(T1);
	objmotion2->transform.Rotation().Get_Value(R2);
	objmotion2->transform.Translation().Get_Value(T2);
	PQP_REAL lamda = 0.0, lastLamda=0, dlamda=0.0;
	PQP_REAL dist=0;
	int nItrs = 0;
	bool T_CA=false;



	C2A_TimeOfContactStep(objmotion1, objmotion2, res, R1,T1,o1,
		R2,T2,o2,tolerance_t,tolerance_d);

	nItrs=0;

	dist = res->distance; 
	mint = res->mint;

	res->numCA = 1;


	if(dist  <= tolerance_d  || mint <= tolerance_t)
	{
		res->collisionfree = false; 			
		res->toc = tolerance_t;
		for(tmpi=0; tmpi<3; tmpi++)
		{	 
			for( tmpj=0; tmpj<3; tmpj++)
			{
				res->R_toc[tmpi][tmpj] = R1[tmpi][tmpj];
			}
			res->T_toc[tmpi] = T1[tmpi];
		}
		if(toc>=1)
			toc=0;	

	}


	while (dist > tolerance_d) //not close enough
	{
		nItrs++;

		if(nItrs > 150)
		{
			break;
		}
		// 1. If dist > path_max, no collision during t= [0, 1] 
		// 2. If dist < path_max, look for TOC
		if( mint >= 1.0)
		{
			res->collisionfree = true;
			res->toc = 0;
			toc = 0;

			return 0;

		}
		else
		{
			dlamda = mint;
			if (dlamda < tolerance_t)
			{												
				break;
			}	
			lamda += dlamda;   				
			if (lamda >= 1.0)
			{
				res->collisionfree = true;
				res->toc = 0;
				toc = 0;
				return 0;

			}
			lastLamda = lamda;				
			res->numCA++;   
			objmotion1->integrate(lamda, R1, T1);
			objmotion2->integrate(lamda, R2, T2); 
			res->UpboundTOC=1.0-lamda;				

			C2A_TimeOfContactStep(objmotion1, objmotion2, res, R1,T1,o1,
				R2,T2,o2,tolerance_t,tolerance_d);

			dist=res->distance;
			mint = res->mint; 

		}
	}

	if(dist == 0 && mint<0)
	{
		toc = lastLamda- dlamda ;
	}
	else
	{
		toc = lastLamda;
	}

	res->collisionfree = false;

	if(toc>=1-tolerance_t)
	{toc=0;}
	res->toc = toc;
	//compute the result R, T at TOC 
	objmotion1->integrate(toc, res->R_toc, res->T_toc);//make res->R_toc equals R1


	return toc;


}


inline void Transform2PQP(Transform *t, PQP_REAL R[3][3], PQP_REAL T[3])
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			R[i][j] = t->Rotation()[i][j];
		}
		T[i] = t->Translation()[i];
	}	
}

inline void PQP2Transform(PQP_REAL R[3][3], PQP_REAL T[3], Transform &t)
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			t.Rotation()[i][j]= R[i][j];
		}
		t.Translation()[i]= T[i];
	}	
}

C2A_Result C2A_Solve(Transform *trans00, 
					 Transform *trans01,  
					 C2A_Model * obj1_tested,
					 Transform *trans10, 
					 Transform *trans11,  
					 C2A_Model * obj2_tested,
					 Transform & trans0, 
					 Transform & trans1, 
					 PQP_REAL &time_of_contact, 
					 int& number_of_iteration, 
					 int & number_of_contact,
					 PQP_REAL th_ca,
					 C2A_TimeOfContactResult& dres)
{


	double lamda = 0.0, lastLamda=-1.0, dlamda=0.0;

	number_of_iteration=0;

	
	PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3];

	PQP_REAL dR1[12], dR2[12];
	int k = 0;
	int i,j;

	Transform2PQP(trans00, R1, T1);
	Transform2PQP(trans10, R2, T2);

	k = 0;
	for( i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++)
		{
			dR1[k] = R1[i][j];			
			dR2[k] = R2[i][j];	
			k++;
		}
	}

	

	PQP_REAL R1e[3][3], T1e[3], R2e[3][3], T2e[3];	

	// get the rotation matrix and translational vector about the final state
	Transform2PQP(trans01, R1e, T1e);
	Transform2PQP(trans11, R2e, T2e);
	
	k = 0;
	for( i = 0; i < 3; i++)
	{
		for( j = 0; j < 3; j++)
		{
			dR1[k] = R1e[i][j];			
			dR2[k] = R2e[i][j];	
			k++;
		}
	}
  

   //perform toc query

	CInterpMotion_Linear motion1(R1, T1, R1e, T1e);
	CInterpMotion_Linear motion2(R2, T2, R2e, T2e);


//	C2A_TimeOfContactResult dres;

	PQP_REAL t_delta = 0.0001;
	PQP_REAL d_delta = 0.0001;

	motion1.m_toc_delta = d_delta;
	motion2.m_toc_delta = d_delta;

    //calculate the toc 
	b_TanslationCCD = false;
	if (motion1.m_angVel<1e-8 && motion2.m_angVel<1e-8)
	{
		b_TanslationCCD = true;
	}
	C2A_QueryTimeOfContact(&motion1, 
		&motion2, 
		&dres, 
		obj1_tested, 
		obj2_tested, 
		d_delta, 
		t_delta, 
		0);
	

    dres.cont_l.clear();

	printf("\n dres.distance: %f, collision free %d\n", dres.distance, dres.collisionfree );


	if(!dres.collisionfree)
	{
		PQP_REAL qua[7];
		Quaternion q;
		Coord3D org;		

		motion1.integrate(dres.toc, qua);
		q.W() = qua[0]; q.X() = qua[1]; q.Y() = qua[2]; q.Z() = qua[3]; 
		org.X() = qua[4]; org.Y() = qua[5]; org.Z() = qua[6];
		trans0.Set_Rotation(q);
		trans0.Set_Translation(org);

		motion2.integrate(dres.toc, qua);
		q.W() = qua[0]; q.X() = qua[1]; q.Y() = qua[2]; q.Z() = qua[3]; 
		org.X() = qua[4]; org.Y() = qua[5]; org.Z() = qua[6];
		Coord3D axis(q.X(),q.Y(),q.Z());


		trans1.Set_Rotation(q);
		trans1.Set_Translation(org);

	
		double threshold =2*dres.distance+0.001;		
	    C2A_QueryContact(&motion1,&motion2,&dres,obj1_tested,obj2_tested,threshold);  
	}

	time_of_contact = dres.toc;
    number_of_contact = dres.num_contact;
	number_of_iteration = dres.numCA;	   

	
	return TOCFound;
	
}

void C2A_Model::ComputeCenterOfMass()
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
	//    com /= 3.0 * total_area;
}

void C2A_Model::ComputeRadius()
{
    int i;
	radius=0.0;
	PQP_REAL d;
    for( i = 0; i < num_tris; i++ )
	{
		PQP_REAL p[3], pr[3];

		p[0] = tris[i].p1[0];
		p[1] = tris[i].p1[1];
		p[2] = tris[i].p1[2];		
		VmV(pr, p, com);
		d = Vlength(pr);		
		if( d > radius ) radius = d;

		p[0] = tris[i].p2[0];
		p[1] = tris[i].p2[1];
		p[2] = tris[i].p2[2];		
		VmV(pr, p, com);
		d = Vlength(pr);	
		if( d > radius ) radius = d;
		
		p[0] = tris[i].p3[0];
		p[1] = tris[i].p3[1];
		p[2] = tris[i].p3[2];		
		VmV(pr, p, com);
		d = Vlength(pr);	
		if( d > radius ) radius = d;
    }
}