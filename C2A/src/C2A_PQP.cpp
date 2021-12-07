/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             S. Gottschalk, E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/
#include <stdio.h>
#include <string.h>
#include "PQP.h"
#include "BVTQ.h"
#include "Build.h"
#include "MatVec.h"
#include "GetTime.h"
#include "TriDist.h"
#include "C2A/C2A_Internal.h"
#include "C2A/C2A.h"
#include "C2A/LinearMath.h"



enum BUILD_STATE
{ 
  PQP_BUILD_STATE_EMPTY,     // empty state, immediately after constructor
  PQP_BUILD_STATE_BEGUN,     // after BeginModel(), state for adding triangles
  PQP_BUILD_STATE_PROCESSED  // after tree has been built, ready to use
};

C2A_Model::C2A_Model()
: PQP_Model(), max_clear_conf_()
{

  bMotionCoherence = false;
  motionTri = 0;
  motionBV = 0;

  level = 3;
//  numLocalTris = pow(2.0, level);
//  localTris = new int[numLocalTris];
  levelPID = 0;
}

C2A_Model::~C2A_Model()
{
	if(trisConst != NULL)
		delete [] trisConst;
}


int
C2A_Model::BeginModel(int n)
{
  // reset to initial state if necessary

  if (build_state != PQP_BUILD_STATE_EMPTY) 
  {
    delete [] b;
    delete [] tris;
  
    num_tris = num_bvs = num_tris_alloced = num_bvs_alloced = 0;
  }

  // prepare model for addition of triangles

  if (n <= 0) n = 8;
  num_tris_alloced = n;
  tris = new C2A_Tri[n];
  if (!tris) 
  {
    fprintf(stderr, "PQP Error!  Out of memory for tri array on "
                    "BeginModel() call!\n");
    return PQP_ERR_MODEL_OUT_OF_MEMORY;  
  }

  // give a warning if called out of sequence

  if (build_state != PQP_BUILD_STATE_EMPTY)
  {
    fprintf(stderr,
            "PQP Warning! Called BeginModel() on a PQP_Model that \n"
            "was not empty. This model was cleared and previous\n"
            "triangle additions were lost.\n");
    build_state = PQP_BUILD_STATE_BEGUN;
    return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
  }
 
  build_state = PQP_BUILD_STATE_BEGUN;
  return PQP_OK;
}

int
C2A_Model::AddTri(const PQP_REAL *p1, 
                  const PQP_REAL *p2, 
                  const PQP_REAL *p3, 
                  int id)
{
  if (build_state == PQP_BUILD_STATE_EMPTY)
  {
    BeginModel();
  }
  else if (build_state == PQP_BUILD_STATE_PROCESSED)
  {
    fprintf(stderr,"PQP Warning! Called AddTri() on PQP_Model \n"
                   "object that was already ended. AddTri() was\n"
                   "ignored.  Must do a BeginModel() to clear the\n"
                   "model for addition of new triangles\n");
    return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
  }
        
  // allocate for new triangles

  if (num_tris >= num_tris_alloced)
  {
    C2A_Tri *temp;
    temp = new C2A_Tri[num_tris_alloced*2];
    if (!temp)
    {
      fprintf(stderr, "PQP Error!  Out of memory for tri array on"
	              " AddTri() call!\n");
      return PQP_ERR_MODEL_OUT_OF_MEMORY;  
    }
    memcpy(temp, tris, sizeof(C2A_Tri)*num_tris);
    delete [] tris;
    tris = temp;
    num_tris_alloced = num_tris_alloced*2;
  }
  
  // initialize the new triangle

  C2A_Tri* c2a_tri = GetTriangle(num_tris);

  c2a_tri->p1[0] = p1[0];
  c2a_tri->p1[1] = p1[1];
  c2a_tri->p1[2] = p1[2];

  c2a_tri->p2[0] = p2[0];
  c2a_tri->p2[1] = p2[1];
  c2a_tri->p2[2] = p2[2];

  c2a_tri->p3[0] = p3[0];
  c2a_tri->p3[1] = p3[1];
  c2a_tri->p3[2] = p3[2];



  c2a_tri->id = id;

  
  num_tris += 1;

  return PQP_OK;
}

void
C2A_Model::SetCenterOfMass(PQP_REAL c[3])
{
	VcV(com, c);
}
PQP_REAL
C2A_Model::Max_Cross_Product(PQP_REAL dir[3])
{
	PQP_REAL t[3];
	PQP_REAL cp[3];
	PQP_REAL v=0.0, v_;
	for (int i=0; i<num_tris; i++)
	{
		t[0] = tris[i].p1[0];		
		t[1] = tris[i].p1[1];		
		t[2] = tris[i].p1[2];
		VcrossV(cp, t, dir);
		v_=(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
		
		if( v_ > v ) v=v_;
		
		t[0] = tris[i].p2[0];		
		t[1] = tris[i].p2[1];		
		t[2] = tris[i].p2[2];
		VcrossV(cp, t, dir);
		v_=(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
		
		if( v_ > v ) v=v_;
		
		t[0] = tris[i].p3[0];		
		t[1] = tris[i].p3[1];		
		t[2] = tris[i].p3[2];
		VcrossV(cp, t, dir);
		v_=(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
		
		if( v_ > v ) 
			v = v_;
	}
	
	return sqrt(v);
}

int
C2A_Model::AddTri(const PQP_REAL *p1, 
                  const PQP_REAL *p2, 
                  const PQP_REAL *p3, 
                  int id, int i1,int i2,int i3)
{
	if (build_state == PQP_BUILD_STATE_EMPTY)
	{
		BeginModel();
	}
	else if (build_state == PQP_BUILD_STATE_PROCESSED)
	{
		fprintf(stderr,"PQP Warning! Called AddTri() on C2A_Model \n"
			"object that was already ended. AddTri() was\n"
			"ignored.  Must do a BeginModel() to clear the\n"
			"model for addition of new triangles\n");
		return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
	}
	
	// allocate for new triangles
	
	if (num_tris >= num_tris_alloced)
	{
		C2A_Tri *temp;
		temp = new C2A_Tri[num_tris_alloced*2];
		if (!temp)
		{
			fprintf(stderr, "PQP Error!  Out of memory for tri array on"
				" AddTri() call!\n");
			return PQP_ERR_MODEL_OUT_OF_MEMORY;  
		}
		memcpy(temp, tris, sizeof(C2A_Tri)*num_tris_alloced);
		delete [] tris;
		tris = temp;
		num_tris_alloced = num_tris_alloced*2;
	}
	
	// initialize the new triangle
	C2A_Tri* c2a_tri = GetTriangle(num_tris);

	c2a_tri->p1[0] = p1[0];
	c2a_tri->p1[1] = p1[1];
	c2a_tri->p1[2] = p1[2];
	
	c2a_tri->p2[0] = p2[0];
	c2a_tri->p2[1] = p2[1];
	c2a_tri->p2[2] = p2[2];
	
	c2a_tri->p3[0] = p3[0];
	c2a_tri->p3[1] = p3[1];
	c2a_tri->p3[2] = p3[2];

	
	//C2A_Tri* c2a_tri = (C2A_Tri*)tris;
	c2a_tri->id = id;
	c2a_tri->Index()[0] = i1;
	c2a_tri->Index()[1] = i2;
	c2a_tri->Index()[2] = i3;
	
	num_tris += 1;
	
	return PQP_OK;
}
void ComputeAbsTra(PQP_REAL R[3][3], PQP_REAL T[3], C2A_Model *o, int bv_ind)
{

	C2A_BV *pBV = (C2A_BV *)o->child(bv_ind);

	// R1*(R2*p + T2) + T1 = R1*R2*p + R1*T2 + T1
    MxM(pBV->R_abs,R, pBV->R);
#if PQP_BV_TYPE & OBB_TYPE
    MxVpV(pBV->To_abs, R, pBV->To,T);
#else
    MxVpV(pBV->Tc_abs, R, pBV->Tr,T);
#endif
		
	
	int l = pBV->Leaf();

	// recursive to compute the C2A_BV of the child
	if(!l)
	{
		int c1 = pBV->first_child;
		int c2 = c1 + 1;
#if PQP_BV_TYPE & OBB_TYPE
		ComputeAbsTra(pBV->R_abs, pBV->To_abs, 
					  o, c1);

		ComputeAbsTra(pBV->R_abs, pBV->To_abs, 
					  o, c2);
#else
		ComputeAbsTra(pBV->R_abs, pBV->Tc_abs, 
					  o, c1);

		ComputeAbsTra(pBV->R_abs, pBV->Tc_abs, 
					  o, c2);
#endif
	}
}

int
C2A_BuildModel(C2A_Model *m);

int
C2A_Model::EndModel()
{
	trisConst = new C2A_Tri[num_tris];
	printf("end model 1\n");
	for(int i = 0; i < num_tris; i++)
	{
		trisConst[i] = ((C2A_Tri*)(tris))[i];
	}
	
	if (build_state == PQP_BUILD_STATE_PROCESSED)
  {
    printf("PQP Warning! Called EndModel() on C2A_Model \n"
                   "object that was already ended. EndModel() was\n"
                   "ignored.  Must do a BeginModel() to clear the\n"
                   "model for addition of new triangles\n");
    return PQP_ERR_BUILD_OUT_OF_SEQUENCE;
  }

  // report error is no tris

  if (num_tris == 0)
  {
    printf("PQP Error! EndModel() called on model with"
                   " no triangles\n");
    return PQP_ERR_BUILD_EMPTY_MODEL;
  }

  // shrink fit tris array 

  C2A_Tri *new_tris;
  if (num_tris_alloced > num_tris)
  {
    new_tris = new C2A_Tri[num_tris];
    if (!new_tris) 
    {
      fprintf(stderr, "PQP Error!  Out of memory for tri array "
                      "in EndModel() call!\n");
      return PQP_ERR_MODEL_OUT_OF_MEMORY;  
    }
    memcpy(new_tris, tris, sizeof(C2A_Tri)*num_tris);
    delete [] tris;
    tris = new_tris;
    num_tris_alloced = num_tris;
  }

  printf("end model 2\n");
  //compute the center of mass and the radius of the model
  ComputeCenterOfMass();
  ComputeRadius();
  
  // create an array of C2A_BVs for the model

  b = new C2A_BV[2*num_tris - 1];
  if (!b)
  {
    fprintf(stderr,"PQP Error! out of memory for C2A_BV array "
                   "in EndModel()\n");
    return PQP_ERR_MODEL_OUT_OF_MEMORY;
  }
  num_bvs_alloced = 2*num_tris - 1;
  num_bvs = 0;

  // we should build the model now.

  printf("end model 3\n");
  C2A_BuildModel(this);
  build_state = PQP_BUILD_STATE_PROCESSED;

  printf("end model 4\n");
  last_tri = tris;

  PQP_REAL R[3][3];
  PQP_REAL T[3];

  R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
  R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
  R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;

  T[0] = T[1] = T[2] = 0;

  // starting from the root
  ComputeAbsTra(R, T, this, 0);


  return PQP_OK;
}


int
C2A_Model::MemUsage(int msg)
{
  int mem_bv_list = sizeof(C2A_BV)*num_bvs;
  int mem_tri_list = sizeof(C2A_Tri)*num_tris;

  int total_mem = mem_bv_list + mem_tri_list + sizeof(C2A_Model);

  if (msg) 
  {
    fprintf(stderr,"Total for model %x: %d bytes\n", this, total_mem);
    fprintf(stderr,"BVs: %d alloced, take %d bytes each\n", 
            num_bvs, sizeof(C2A_BV));
    fprintf(stderr,"Tris: %d alloced, take %d bytes each\n", 
            num_tris, sizeof(C2A_Tri));
  }
  
  return total_mem;
}

int
TriContact(PQP_REAL *P1, PQP_REAL *P2, PQP_REAL *P3,
		   PQP_REAL *Q1, PQP_REAL *Q2, PQP_REAL *Q3); 
// Defined in PQP.cpp
PQP_REAL
TriDistance(PQP_REAL R[3][3], PQP_REAL T[3], Tri *t1, Tri *t2,
            PQP_REAL p[3], PQP_REAL q[3]);

// TRIANGLE OVERLAP TEST
   
PQP_REAL
TriDistance_Dirsection(PQP_REAL R[3][3], PQP_REAL T[3], Tri *t1, Tri *t2,
			PQP_REAL p[3], PQP_REAL q[3], PQP_REAL direction[3])
{
	// transform tri 2 into same space as tri 1

	PQP_REAL tri1[3][3], tri2[3][3];

	VcV(tri1[0], t1->p1);
	VcV(tri1[1], t1->p2);
	VcV(tri1[2], t1->p3);
	MxVpV(tri2[0], R, t2->p1, T);
	MxVpV(tri2[1], R, t2->p2, T);
	MxVpV(tri2[2], R, t2->p3, T);



	return TriDist(p,q,tri1,tri2);
}



#if PQP_BV_TYPE & RSS_TYPE // distance/tolerance only available with RSS
                           // unless an OBB distance test is supplied in 
                           // C2A_BV.cpp

// DISTANCE STUFF
//
//--------------------------------------------------------------------------

void
C2ADistanceRecurse(C2A_DistanceResult *res,
                PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
                C2A_Model *o1, int b1,
                C2A_Model *o2, int b2)
{
  PQP_REAL sz1 = (o1->child(b1))->GetSize();
  PQP_REAL sz2 = (o2->child(b2))->GetSize();


  int l1 = (o1->child(b1))->Leaf();
  int l2 = (o2->child(b2))->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.

    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    Tri *t1 = o1->GetTriangle(-(o1->child(b1))->first_child - 1);
    Tri *t2 = o2->GetTriangle(-(o2->child(b2))->first_child - 1);

    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
  
    if (d < res->distance) 
    {
      res->distance = d;

	  res->t1 = t1->id;
	  res->t2 = t2->id;
	  
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed 
                               // into c.s. 2 later
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

    a1 = (o1->child(b1))->first_child;
    a2 = b2;
    c1 = (o1->child(b1))->first_child+1;
    c2 = b2;
    
    MTxM(R1,(o1->child(a1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,(o1->child(a1))->R,Ttemp);

    MTxM(R2,(o1->child(c1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,(o1->child(c1))->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,(o1->child(c1))->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = (o2->child(b2))->first_child;
    c1 = b1;
    c2 = (o2->child(b2))->first_child+1;

    MxM(R1,R,(o2->child(a2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,(o2->child(a2))->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif

    MxM(R2,R,(o2->child(c2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,(o2->child(c2))->Tr,T);
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
    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      C2ADistanceRecurse(res, R2, T2, o1, c1, o2, c2);      
    }

    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      C2ADistanceRecurse(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else 
  {
    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      C2ADistanceRecurse(res, R1, T1, o1, a1, o2, a2);
    }

    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      C2ADistanceRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
  }
}






void
C2ADistanceQueueRecurse(C2A_DistanceResult *res, 
                     PQP_REAL R[3][3], PQP_REAL T[3],
                     C2A_Model *o1, int b1,
                     C2A_Model *o2, int b2)
{
  BVTQ bvtq(res->qsize);

  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R,R);
  VcV(min_test.T,T);

  while(1) 
  {  
    int l1 = (o1->child(min_test.b1))->Leaf();
    int l2 = (o2->child(min_test.b2))->Leaf();
    
    if (l1 && l2) 
    {  
      // both leaves.  Test the triangles beneath them.

      res->num_tri_tests++;

      PQP_REAL p[3], q[3];

      Tri *t1 = o1->GetTriangle(-(o1->child(min_test.b1))->first_child - 1);
      Tri *t2 = o2->GetTriangle(-(o2->child(min_test.b2))->first_child - 1);

      PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
  
      if (d < res->distance)
      {
        res->distance = d;

		res->t1 = t1->id;
		res->t2 = t2->id;
		
        VcV(res->p1, p);         // p already in c.s. 1
        VcV(res->p2, q);         // q must be transformed 
                                 // into c.s. 2 later
        o1->last_tri = t1;
        o2->last_tri = t2;
      }
    }		 
    else if (bvtq.GetNumTests() == bvtq.GetSize() - 1) 
    {  
      // queue can't get two more tests, recur
      
      C2ADistanceQueueRecurse(res,min_test.R,min_test.T,
                           o1,min_test.b1,o2,min_test.b2);
    }
    else 
    {  
      // decide how to descend to children
      
      PQP_REAL sz1 = (o1->child(min_test.b1))->GetSize();
      PQP_REAL sz2 = (o2->child(min_test.b2))->GetSize();

      res->num_bv_tests += 2;
 
      BVT bvt1,bvt2;
      PQP_REAL Ttemp[3];

      if (l2 || (!l1 && (sz1 > sz2)))	
      {  
        // put new tests on queue consisting of min_test.b2 
        // with children of min_test.b1 
      
        int c1 = (o1->child(min_test.b1))->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R,(o1->child(c1))->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,(o1->child(c1))->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c1)->To);
#endif
        MTxV(bvt1.T,(o1->child(c1))->R,Ttemp);
		PQP_REAL S[3];
		
        bvt1.d = C2A_BV_Distance(bvt1.R,bvt1.T,
                            (C2A_BV *)o1->child(bvt1.b1),(C2A_BV *)o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = c2;
        bvt2.b2 = min_test.b2;
        MTxM(bvt2.R,(o1->child(c2))->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,(o1->child(c2))->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c2)->To);
#endif
        MTxV(bvt2.T,(o1->child(c2))->R,Ttemp);
		
		
        bvt2.d = C2A_BV_Distance(bvt2.R,bvt2.T,
                            (C2A_BV *)o1->child(bvt2.b1),(C2A_BV *)o2->child(bvt2.b2),S);
      }
      else 
      {
        // put new tests on queue consisting of min_test.b1 
        // with children of min_test.b2
      
        int c1 = (o2->child(min_test.b2))->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R,min_test.R,(o2->child(c1))->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T,min_test.R,(o2->child(c1))->Tr,min_test.T);
#else
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->To,min_test.T);
#endif
		PQP_REAL S[3];
		
        bvt1.d = C2A_BV_Distance(bvt1.R,bvt1.T,
                            (C2A_BV *)o1->child(bvt1.b1),(C2A_BV *)o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R,min_test.R,(o2->child(c2))->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T,min_test.R,(o2->child(c2))->Tr,min_test.T);
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

      if ((min_test.d + res->abs_err >= res->distance) && 
         ((min_test.d * (1 + res->rel_err)) >= res->distance)) 
      {
        break;
      }
    }
  }  
}	


void PrintTri(FILE* output, Tri &ftri)
{
	fprintf(output,"\n p1:%f, %f, %f", ftri.p1[0], ftri.p1[1],  ftri.p1[2]);
	fprintf(output,"\n p2:%f, %f, %f", ftri.p2[0], ftri.p2[1],  ftri.p2[2]);
	fprintf(output,"\n p3:%f, %f, %f", ftri.p3[0], ftri.p3[1],  ftri.p3[2]);
	
}

void
CollideRecurse(PQP_CollideResult *res,
			   PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
			   C2A_Model *o1, int b1, 
			   C2A_Model *o2, int b2, int flag)
{
	// first thing, see if we're overlapping

	res->num_bv_tests++;

	if (!BV_Overlap(R, T, o1->child(b1), o2->child(b2))) return;

	// if we are, see if we test triangles next

	int l1 = o1->child(b1)->Leaf();
	int l2 = o2->child(b2)->Leaf();

	if (l1 && l2) 
	{
		res->num_tri_tests++;

#if 1
		// transform the points in b2 into space of b1, then compare

		Tri *t1 = o1->GetTriangle(-(o1->child(b1))->first_child - 1);
		Tri *t2 = o2->GetTriangle(-(o2->child(b2))->first_child - 1);
	//	printf("%d; ",t1->id);
		PQP_REAL q1[3], q2[3], q3[3];
		PQP_REAL *p1 = t1->p1;
		PQP_REAL *p2 = t1->p2;
		PQP_REAL *p3 = t1->p3;    
		MxVpV(q1, res->R, t2->p1, res->T);
		MxVpV(q2, res->R, t2->p2, res->T);
		MxVpV(q3, res->R, t2->p3, res->T);
		if (TriContact(p1, p2, p3, q1, q2, q3)) 
		{
			// add this to result

			res->Add(t1->id, t2->id);
		}
#else
		PQP_REAL p[3], q[3];

		Tri *t1 = &o1->tris[-o1->child(b1)->first_child - 1];
		Tri *t2 = &o2->tris[-o2->child(b2)->first_child - 1];

		if (TriDistance(res->R,res->T,t1,t2,p,q) == 0.0)
		{
			// add this to result

			res->Add(t1->id, t2->id);
		}
#endif

		return;
	}

	// we dont, so decide whose children to visit next

	PQP_REAL sz1 = o1->child(b1)->GetSize();
	PQP_REAL sz2 = o2->child(b2)->GetSize();

	PQP_REAL Rc[3][3],Tc[3],Ttemp[3];

	if (l2 || (!l1 && (sz1 > sz2)))
	{
		int c1 = o1->child(b1)->first_child;
		int c2 = c1 + 1;

		MTxM(Rc,o1->child(c1)->R,R);
#if PQP_BV_TYPE & OBB_TYPE
		VmV(Ttemp,T,o1->child(c1)->To);
#else
		VmV(Ttemp,T,o1->child(c1)->Tr);
#endif
		MTxV(Tc,o1->child(c1)->R,Ttemp);
		CollideRecurse(res,Rc,Tc,o1,c1,o2,b2,flag);

		if ((flag == PQP_FIRST_CONTACT) && (res->num_pairs > 0)) return;

		MTxM(Rc,o1->child(c2)->R,R);
#if PQP_BV_TYPE & OBB_TYPE
		VmV(Ttemp,T,o1->child(c2)->To);
#else
		VmV(Ttemp,T,o1->child(c2)->Tr);
#endif
		MTxV(Tc,o1->child(c2)->R,Ttemp);
		CollideRecurse(res,Rc,Tc,o1,c2,o2,b2,flag);
	}
	else 
	{
		int c1 = o2->child(b2)->first_child;
		int c2 = c1 + 1;

		MxM(Rc,R,o2->child(c1)->R);
#if PQP_BV_TYPE & OBB_TYPE
		MxVpV(Tc,R,o2->child(c1)->To,T);
#else
		MxVpV(Tc,R,o2->child(c1)->Tr,T);
#endif
		CollideRecurse(res,Rc,Tc,o1,b1,o2,c1,flag);

		if ((flag == PQP_FIRST_CONTACT) && (res->num_pairs > 0)) return;

		MxM(Rc,R,o2->child(c2)->R);
#if PQP_BV_TYPE & OBB_TYPE
		MxVpV(Tc,R,o2->child(c2)->To,T);
#else
		MxVpV(Tc,R,o2->child(c2)->Tr,T);
#endif
		CollideRecurse(res,Rc,Tc,o1,b1,o2,c2,flag);
	}
}
int 
C2A_Collide(PQP_CollideResult *res,
			PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
			PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
			int flag)
{
	double t1 = GetTime();

	// make sure that the models are built

	if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;

	// clear the stats

	res->num_bv_tests = 0;
	res->num_tri_tests = 0;

	// don't release the memory, but reset the num_pairs counter

	res->num_pairs = 0;

	// Okay, compute what transform [R,T] that takes us from cs1 to cs2.
	// [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
	// First compute the rotation part, then translation part

	MTxM(res->R,R1,R2);
	PQP_REAL Ttemp[3];
	VmV(Ttemp, T2, T1);  
	MTxV(res->T, R1, Ttemp);

	// compute the transform from o1->child(0) to o2->child(0)

	PQP_REAL Rtemp[3][3], R[3][3], T[3];

	MxM(Rtemp,res->R,o2->child(0)->R);
	MTxM(R,o1->child(0)->R,Rtemp);

#if PQP_BV_TYPE & OBB_TYPE
	MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->To);
#else
	MxVpV(Ttemp,res->R,o2->child(0)->Tr,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->Tr);
#endif

	MTxV(T,o1->child(0)->R,Ttemp);

	// now start with both top level BVs  

	CollideRecurse(res,R,T,o1,0,o2,0,flag);

	double t2 = GetTime();
	res->query_time_secs = t2 - t1;

	return PQP_OK; 
}

int 
C2A_Distance(C2A_DistanceResult *res,
             PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err,
             int qsize)
{
  
  double time1 = GetTime();
  
  // make sure that the models are built

  if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;
  if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;

  // Okay, compute what transform [R,T] that takes us from cs2 to cs1.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part

  MTxM(res->R,R1,R2);
  PQP_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);  
  MTxV(res->T, R1, Ttemp);
  
  // establish initial upper bound using last triangles which 
  // provided the minimum distance

  PQP_REAL p[3],q[3];
  res->distance = TriDistance(res->R,res->T,o1->last_tri,o2->last_tri,p,q);
  res->t1 = o1->last_tri->id;
  res->t2 = o2->last_tri->id;
		
  VcV(res->p1,p);
  VcV(res->p2,q);

  // initialize error bounds

  res->abs_err = abs_err;
  res->rel_err = rel_err;
  
  // clear the stats

  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // compute the transform from o1->child(0) to o2->child(0)

  PQP_REAL Rtemp[3][3], R[3][3], T[3];

  MxM(Rtemp,res->R,(o2->child(0))->R);
  MTxM(R,(o1->child(0))->R,Rtemp);
  
#if PQP_BV_TYPE & RSS_TYPE
  MxVpV(Ttemp,res->R,(o2->child(0))->Tr,res->T);
  VmV(Ttemp,Ttemp,(o1->child(0))->Tr);
#else
  MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
  MTxV(T,(o1->child(0))->R,Ttemp);

  // choose routine according to queue size
  
  if (qsize <= 2)
  {
    C2ADistanceRecurse(res,R,T,o1,0,o2,0);    
  }
  else 
  { 
    res->qsize = qsize;

    C2ADistanceQueueRecurse(res,R,T,o1,0,o2,0);
  }

  // res->p2 is in cs 1 ; transform it to cs 2

  PQP_REAL u[3];
  VmV(u, res->p2, res->T);
  MTxV(res->p2, res->R, u);

  double time2 = GetTime();
  res->query_time_secs = time2 - time1;  

  return PQP_OK;
}


void
C2ACollideRecurse(C2A_DistanceResult *res,
                PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
                C2A_Model *o1, int b1,
                C2A_Model *o2, int b2)
{
  PQP_REAL sz1 = (o1->child(b1))->GetSize();
  PQP_REAL sz2 = (o2->child(b2))->GetSize();

  if (!C2A_BV_Overlap(R, T, (C2A_BV *)o1->child(b1), (C2A_BV *)o2->child(b2))) {
	  return;
  }

  int l1 = (o1->child(b1))->Leaf();
  int l2 = (o2->child(b2))->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.

    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    Tri *t1 = o1->GetTriangle(-(o1->child(b1))->first_child - 1);
    Tri *t2 = o2->GetTriangle(-(o2->child(b2))->first_child - 1);

    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
  
    if (d < res->distance) 
    {
      res->distance = d;

	  res->t1 = t1->id;
	  res->t2 = t2->id;
	  
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed 
                               // into c.s. 2 later
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

    a1 = (o1->child(b1))->first_child;
    a2 = b2;
    c1 = (o1->child(b1))->first_child+1;
    c2 = b2;
    
    MTxM(R1,(o1->child(a1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,(o1->child(a1))->R,Ttemp);

    MTxM(R2,(o1->child(c1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,(o1->child(c1))->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,(o1->child(c1))->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = (o2->child(b2))->first_child;
    c1 = b1;
    c2 = (o2->child(b2))->first_child+1;

    MxM(R1,R,(o2->child(a2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,(o2->child(a2))->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif

    MxM(R2,R,(o2->child(c2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,(o2->child(c2))->Tr,T);
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
    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      C2ACollideRecurse(res, R2, T2, o1, c1, o2, c2);      
    }

    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      C2ACollideRecurse(res, R1, T1, o1, a1, o2, a2);
    }
  }
  else 
  {
    if ((d1 < (res->distance - res->abs_err)) || 
        (d1*(1 + res->rel_err) < res->distance)) 
    {      
      C2ACollideRecurse(res, R1, T1, o1, a1, o2, a2);
    }

    if ((d2 < (res->distance - res->abs_err)) || 
        (d2*(1 + res->rel_err) < res->distance)) 
    {      
      C2ACollideRecurse(res, R2, T2, o1, c1, o2, c2);      
    }
  }
}


int 
C2A_Collide(C2A_DistanceResult *res,
             PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err,
             int qsize)
{

  res->distance = 10.0;
  
  double time1 = GetTime();
  
  // make sure that the models are built

  if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;
  if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;

  // Okay, compute what transform [R,T] that takes us from cs2 to cs1.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part

  MTxM(res->R,R1,R2);
  PQP_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);  
  MTxV(res->T, R1, Ttemp);
  
  // establish initial upper bound using last triangles which 
  // provided the minimum distance

  PQP_REAL p[3],q[3];
  res->distance = TriDistance(res->R,res->T,o1->last_tri,o2->last_tri,p,q);
  res->t1 = o1->last_tri->id;
  res->t2 = o2->last_tri->id;
		
  VcV(res->p1,p);
  VcV(res->p2,q);

  // initialize error bounds

  res->abs_err = abs_err;
  res->rel_err = rel_err;
  
  // clear the stats

  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // compute the transform from o1->child(0) to o2->child(0)

  PQP_REAL Rtemp[3][3], R[3][3], T[3];

  MxM(Rtemp,res->R,(o2->child(0))->R);
  MTxM(R,(o1->child(0))->R,Rtemp);
  
#if PQP_BV_TYPE & RSS_TYPE
  MxVpV(Ttemp,res->R,(o2->child(0))->Tr,res->T);
  VmV(Ttemp,Ttemp,(o1->child(0))->Tr);
#else
  MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
  MTxV(T,(o1->child(0))->R,Ttemp);

  // choose routine according to queue size

  C2ACollideRecurse(res,R,T,o1,0,o2,0);    


  // res->p2 is in cs 1 ; transform it to cs 2

  PQP_REAL u[3];
  VmV(u, res->p2, res->T);
  MTxV(res->p2, res->R, u);

  double time2 = GetTime();
  res->query_time_secs = time2 - time1;  

  return PQP_OK;
}

void
C2ADistance_At_ContactSpaceRecurse(C2A_DistanceResult *res,
                PQP_REAL R[3][3], PQP_REAL T[3], // b2 relative to b1
                C2A_Model *o1, int b1,
                C2A_Model *o2, int b2,
				PQP_REAL upper_bound)
{

  if (res->distance < upper_bound) return;

  PQP_REAL sz1 = (o1->child(b1))->GetSize();
  PQP_REAL sz2 = (o2->child(b2))->GetSize();


  int l1 = (o1->child(b1))->Leaf();
  int l2 = (o2->child(b2))->Leaf();

  if (l1 && l2)
  {
    // both leaves.  Test the triangles beneath them.

    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    Tri *t1 = o1->GetTriangle(-(o1->child(b1))->first_child - 1);
    Tri *t2 = o2->GetTriangle(-(o2->child(b2))->first_child - 1);

    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
  
    if (d < res->distance) 
    {
      res->distance = d;

	  res->t1 = t1->id;
	  res->t2 = t2->id;
	  
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed 
                               // into c.s. 2 later
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

    a1 = (o1->child(b1))->first_child;
    a2 = b2;
    c1 = (o1->child(b1))->first_child+1;
    c2 = b2;
    
    MTxM(R1,(o1->child(a1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,o1->child(a1)->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,(o1->child(a1))->R,Ttemp);

    MTxM(R2,(o1->child(c1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,(o1->child(c1))->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,(o1->child(c1))->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = (o2->child(b2))->first_child;
    c1 = b1;
    c2 = (o2->child(b2))->first_child+1;

    MxM(R1,R,(o2->child(a2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,(o2->child(a2))->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif

    MxM(R2,R,(o2->child(c2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,(o2->child(c2))->Tr,T);
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
    if (d2 < res->distance) 
    {      
      C2ADistance_At_ContactSpaceRecurse(res, R2, T2, o1, c1, o2, c2, upper_bound);      
    }

    if (d1 < res->distance) 
    {      
      C2ADistance_At_ContactSpaceRecurse(res, R1, T1, o1, a1, o2, a2, upper_bound);
    }
  }
  else 
  {
    if (d1 < res->distance) 
    {      
      C2ADistance_At_ContactSpaceRecurse(res, R1, T1, o1, a1, o2, a2, upper_bound);
    }

    if (d2 < res->distance) 
    {      
      C2ADistance_At_ContactSpaceRecurse(res, R2, T2, o1, c1, o2, c2, upper_bound);      
    }
  }
}

int 
C2A_Distance_At_ContactSpace(C2A_DistanceResult *res, 
             PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
			 PQP_REAL upper_bound)
{
  double time1 = GetTime();
  
  // make sure that the models are built

  if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;
  if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
    return PQP_ERR_UNPROCESSED_MODEL;

  // Okay, compute what transform [R,T] that takes us from cs2 to cs1.
  // [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
  // First compute the rotation part, then translation part

  MTxM(res->R,R1,R2);
  PQP_REAL Ttemp[3];
  VmV(Ttemp, T2, T1);  
  MTxV(res->T, R1, Ttemp);
  
  // establish initial upper bound using last triangles which 
  // provided the minimum distance

  //PQP_REAL p[3],q[3];
  //res->distance = TriDistance(res->R,res->T,o1->last_tri,o2->last_tri,p,q);
  //res->t1 = o1->last_tri->id;
  //res->t2 = o2->last_tri->id;
		
  //VcV(res->p1,p);
  //VcV(res->p2,q);

  // initialize error bounds

  res->distance = upper_bound;
  res->abs_err = 0.0;
  res->rel_err = 0.0;
  
  // clear the stats

  res->num_bv_tests = 0;
  res->num_tri_tests = 0;
  
  // compute the transform from o1->child(0) to o2->child(0)

  PQP_REAL Rtemp[3][3], R[3][3], T[3];

  MxM(Rtemp,res->R,(o2->child(0))->R);
  MTxM(R,(o1->child(0))->R,Rtemp);
  
#if PQP_BV_TYPE & RSS_TYPE
  MxVpV(Ttemp,res->R,(o2->child(0))->Tr,res->T);
  VmV(Ttemp,Ttemp,(o1->child(0))->Tr);
#else
  MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
  VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
  MTxV(T,(o1->child(0))->R,Ttemp);

  C2ADistance_At_ContactSpaceRecurse(res,R,T,o1,0,o2,0, upper_bound);    


  // res->p2 is in cs 1 ; transform it to cs 2

  PQP_REAL u[3];
  VmV(u, res->p2, res->T);
  MTxV(res->p2, res->R, u);

  double time2 = GetTime();
  res->query_time_secs = time2 - time1;  

  return PQP_OK;
}

// Tolerance Stuff
//
//---------------------------------------------------------------------------
void 
C2AToleranceRecurse(PQP_ToleranceResult *res, 
                 PQP_REAL R[3][3], PQP_REAL T[3],
                 C2A_Model *o1, int b1, C2A_Model *o2, int b2)
{
  PQP_REAL sz1 = (o1->child(b1))->GetSize();
  PQP_REAL sz2 = (o2->child(b2))->GetSize();
  int l1 = (o1->child(b1))->Leaf();
  int l2 = (o2->child(b2))->Leaf();

  if (l1 && l2) 
  {
    // both leaves - find if tri pair within tolerance
    
    res->num_tri_tests++;

    PQP_REAL p[3], q[3];

    Tri *t1 = o1->GetTriangle(-(o1->child(b1))->first_child - 1);
    Tri *t2 = o2->GetTriangle(-(o2->child(b2))->first_child - 1);

    PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
    
    if (d <= res->tolerance)  
    {  
      // triangle pair distance less than tolerance

      res->closer_than_tolerance = 1;
      res->distance = d;
      VcV(res->p1, p);         // p already in c.s. 1
      VcV(res->p2, q);         // q must be transformed 
                               // into c.s. 2 later
    }

    return;
  }

  int a1,a2,c1,c2;  // new bv tests 'a' and 'c'
  PQP_REAL R1[3][3], T1[3], R2[3][3], T2[3], Ttemp[3];

  if (l2 || (!l1 && (sz1 > sz2)))
  {
    // visit the children of b1

    a1 = (o1->child(b1))->first_child;
    a2 = b2;
    c1 = (o1->child(b1))->first_child+1;
    c2 = b2;
    
    MTxM(R1,(o1->child(a1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,(o1->child(a1))->Tr);
#else
    VmV(Ttemp,T,o1->child(a1)->To);
#endif
    MTxV(T1,(o1->child(a1))->R,Ttemp);

    MTxM(R2,(o1->child(c1))->R,R);
#if PQP_BV_TYPE & RSS_TYPE
    VmV(Ttemp,T,(o1->child(c1))->Tr);
#else
    VmV(Ttemp,T,o1->child(c1)->To);
#endif
    MTxV(T2,(o1->child(c1))->R,Ttemp);
  }
  else 
  {
    // visit the children of b2

    a1 = b1;
    a2 = (o2->child(b2))->first_child;
    c1 = b1;
    c2 = (o2->child(b2))->first_child+1;

    MxM(R1,R,(o2->child(a2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T1,R,(o2->child(a2))->Tr,T);
#else
    MxVpV(T1,R,o2->child(a2)->To,T);
#endif
    MxM(R2,R,(o2->child(c2))->R);
#if PQP_BV_TYPE & RSS_TYPE
    MxVpV(T2,R,(o2->child(c2))->Tr,T);
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
    if (d2 <= res->tolerance) C2AToleranceRecurse(res, R2, T2, o1, c1, o2, c2);
    if (res->closer_than_tolerance) return;
    if (d1 <= res->tolerance) C2AToleranceRecurse(res, R1, T1, o1, a1, o2, a2);
  }
  else 
  {
    if (d1 <= res->tolerance) C2AToleranceRecurse(res, R1, T1, o1, a1, o2, a2);
    if (res->closer_than_tolerance) return;
    if (d2 <= res->tolerance) C2AToleranceRecurse(res, R2, T2, o1, c1, o2, c2);
  }
}

void
C2AToleranceQueueRecurse(PQP_ToleranceResult *res,
                      PQP_REAL R[3][3], PQP_REAL T[3],
                      C2A_Model *o1, int b1,
                      C2A_Model *o2, int b2)
{
  BVTQ bvtq(res->qsize);
  BVT min_test;
  min_test.b1 = b1;
  min_test.b2 = b2;
  McM(min_test.R,R);
  VcV(min_test.T,T);

  while(1)
  {  
    int l1 = (o1->child(min_test.b1))->Leaf();
    int l2 = (o2->child(min_test.b2))->Leaf();
    
    if (l1 && l2) 
    {  
      // both leaves - find if tri pair within tolerance
    
      res->num_tri_tests++;

      PQP_REAL p[3], q[3];

      Tri *t1 = o1->GetTriangle(-o1->child(min_test.b1)->first_child - 1);
      Tri *t2 = o2->GetTriangle(-o2->child(min_test.b2)->first_child - 1);

      PQP_REAL d = TriDistance(res->R,res->T,t1,t2,p,q);
    
      if (d <= res->tolerance)  
      {  
        // triangle pair distance less than tolerance

        res->closer_than_tolerance = 1;
        res->distance = d;
        VcV(res->p1, p);         // p already in c.s. 1
        VcV(res->p2, q);         // q must be transformed 
                                 // into c.s. 2 later
        return;
      }
    }
    else if (bvtq.GetNumTests() == bvtq.GetSize() - 1)
    {  
      // queue can't get two more tests, recur
      
      C2AToleranceQueueRecurse(res,min_test.R,min_test.T,
                            o1,min_test.b1,o2,min_test.b2);
      if (res->closer_than_tolerance == 1) return;
    }
    else 
    {  
      // decide how to descend to children
      
      PQP_REAL sz1 = (o1->child(min_test.b1))->GetSize();
      PQP_REAL sz2 = (o2->child(min_test.b2))->GetSize();

      res->num_bv_tests += 2;
      
      BVT bvt1,bvt2;
      PQP_REAL Ttemp[3];

      if (l2 || (!l1 && (sz1 > sz2)))	
      {
	      // add two new tests to queue, consisting of min_test.b2
        // with the children of min_test.b1

        int c1 = (o1->child(min_test.b1))->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = c1;
        bvt1.b2 = min_test.b2;
        MTxM(bvt1.R,(o1->child(c1))->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
        VmV(Ttemp,min_test.T,(o1->child(c1))->Tr);
#else
        VmV(Ttemp,min_test.T,o1->child(c1)->To);
#endif
		PQP_REAL S[3];
		
        MTxV(bvt1.T,(o1->child(c1))->R,Ttemp);
        bvt1.d = C2A_BV_Distance(bvt1.R,bvt1.T,
                            (C2A_BV *)o1->child(bvt1.b1),(C2A_BV *)o2->child(bvt1.b2),S);

	      // init bv test 2

	      bvt2.b1 = c2;
	      bvt2.b2 = min_test.b2;
	      MTxM(bvt2.R,(o1->child(c2))->R,min_test.R);
#if PQP_BV_TYPE & RSS_TYPE
	      VmV(Ttemp,min_test.T,(o1->child(c2))->Tr);
#else
	      VmV(Ttemp,min_test.T,o1->child(c2)->To);
#endif
	      MTxV(bvt2.T,(o1->child(c2))->R,Ttemp);
        bvt2.d = C2A_BV_Distance(bvt2.R,bvt2.T,
                            (C2A_BV *)o1->child(bvt2.b1),(C2A_BV *)o2->child(bvt2.b2),S);
      }
      else 
      {
        // add two new tests to queue, consisting of min_test.b1
        // with the children of min_test.b2

        int c1 = (o2->child(min_test.b2))->first_child;
        int c2 = c1 + 1;

        // init bv test 1

        bvt1.b1 = min_test.b1;
        bvt1.b2 = c1;
        MxM(bvt1.R,min_test.R, (o2->child(c1))->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt1.T,min_test.R,(o2->child(c1))->Tr,min_test.T);
#else
        MxVpV(bvt1.T,min_test.R,o2->child(c1)->To,min_test.T);
#endif
		PQP_REAL S[3];
		
        bvt1.d = C2A_BV_Distance(bvt1.R,bvt1.T,
                            (C2A_BV *)o1->child(bvt1.b1),(C2A_BV *)o2->child(bvt1.b2),S);

        // init bv test 2

        bvt2.b1 = min_test.b1;
        bvt2.b2 = c2;
        MxM(bvt2.R,min_test.R,(o2->child(c2))->R);
#if PQP_BV_TYPE & RSS_TYPE
        MxVpV(bvt2.T,min_test.R,(o2->child(c2))->Tr,min_test.T);
#else
        MxVpV(bvt2.T,min_test.R,o2->child(c2)->To,min_test.T);
#endif
        bvt2.d = C2A_BV_Distance(bvt2.R,bvt2.T,
                            (C2A_BV *)o1->child(bvt2.b1), (C2A_BV *)o2->child(bvt2.b2),S);
      }

      // put children tests in queue

      if (bvt1.d <= res->tolerance) bvtq.AddTest(bvt1);
      if (bvt2.d <= res->tolerance) bvtq.AddTest(bvt2);
    }

    if (bvtq.Empty() || (bvtq.MinTest() > res->tolerance)) 
    {
      res->closer_than_tolerance = 0;
      return;
    }
    else 
    {
      min_test = bvtq.ExtractMinTest();
    }
  }  
}	

int
C2A_Tolerance(PQP_ToleranceResult *res,
              PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
              PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
              PQP_REAL tolerance,
              int qsize)
{
	double time1 = GetTime();
	
	// make sure that the models are built
	
	if (o1->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	if (o2->build_state != PQP_BUILD_STATE_PROCESSED) 
		return PQP_ERR_UNPROCESSED_MODEL;
	
	// Compute the transform [R,T] that takes us from cs2 to cs1.
	// [R,T] = [R1,T1]'[R2,T2] = [R1',-R1'T][R2,T2] = [R1'R2, R1'(T2-T1)]
	
	MTxM(res->R,R1,R2);
	PQP_REAL Ttemp[3];
	VmV(Ttemp, T2, T1);
	MTxV(res->T, R1, Ttemp);
	
	// set tolerance, used to prune the search
	
	if (tolerance < 0.0) tolerance = 0.0;
	res->tolerance = tolerance;
	
	// clear the stats
	
	res->num_bv_tests = 0;
	res->num_tri_tests = 0;
	
	// initially assume not closer than tolerance
	
	res->closer_than_tolerance = 0;
	
	// compute the transform from o1->child(0) to o2->child(0)
	
	PQP_REAL Rtemp[3][3], R[3][3], T[3];
	PQP_REAL S[3];
	
	
	MxM(Rtemp,res->R,(o2->child(0))->R);
	MTxM(R,(o1->child(0))->R,Rtemp);
#if PQP_BV_TYPE & RSS_TYPE
	MxVpV(Ttemp,res->R,(o2->child(0))->Tr,res->T);
	VmV(Ttemp,Ttemp,(o1->child(0))->Tr);
#else
	MxVpV(Ttemp,res->R,o2->child(0)->To,res->T);
	VmV(Ttemp,Ttemp,o1->child(0)->To);
#endif
	MTxV(T,(o1->child(0))->R,Ttemp);
	
	// find a distance lower bound for trivial reject
	
	PQP_REAL d = C2A_BV_Distance(R, T, (C2A_BV *)o1->child(0), (C2A_BV *)o2->child(0),S);
	
	if (d <= res->tolerance)
	{
		// more work needed - choose routine according to queue size
		
		if (qsize <= 2) 
		{
			C2AToleranceRecurse(res, R, T, o1, 0, o2, 0);
		}
		else 
		{
			res->qsize = qsize;
			C2AToleranceQueueRecurse(res, R, T, o1, 0, o2, 0);
		}
	}
	
	// res->p2 is in cs 1 ; transform it to cs 2
	
	PQP_REAL u[3];
	VmV(u, res->p2, res->T);
	MTxV(res->p2, res->R, u);
	
	double time2 = GetTime();
	res->query_time_secs = time2 - time1;
	
	return PQP_OK;
}

#endif

