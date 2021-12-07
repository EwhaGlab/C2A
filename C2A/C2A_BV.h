/*************************************************************************\

  Copyright 2009	Computer Graphics Lab, 
									Ewha Womans University,
									Seoul, Korea.
  All Rights Reserved.

  The authors may be contacted via:

  Mailing Address:     Xinyu Zhang
                       Min Tang
                       Department of Computer Science & Engineering
                       Ewha Womans University
                       Seoul, 120-750, Korea.

  EMail:               zhangxy@ewha.ac.kr
                       tangmin@ewha.ac.kr`


\**************************************************************************/

#ifndef C2A_BV_H
#define C2A_BV_H

#include <math.h>


#include <PQP_Compile.h>
#include <BV.h>

struct C2A_Tri;

class C2A_BV: public BV
{
public:

	PQP_REAL R_abs[3][3]; // abs orientation of RSS & OBB: Liangjun
	PQP_REAL R_loc[3][3];
	PQP_REAL trilength;

#if PQP_BV_TYPE & RSS_TYPE

	PQP_REAL Tr_abs[3];   // abs position of rectangle: Liangjun

	PQP_REAL Corner[3][3];// the four corners of rectangle 
	PQP_REAL Length_C1;// the length of C1
	PQP_REAL Length_C2;// the length of C2

	int deep;


	PQP_REAL max_Length;
#endif

#if PQP_BV_TYPE & OBB_TYPE
	PQP_REAL To_abs[3];   // abs position of obb
#endif



	int parentID; //the id of parent C2A_BV
	PQP_REAL com[3]; // the center of Mass, 
	PQP_REAL comRoot[3]; // the center of model that C2A_BV belongs to, 
	PQP_REAL angularRadius; // the angular radius of the C2A_BV to comRoot, 

	PQP_REAL C2A_BV_Size;

	void ComputeCenterOfMassBV(C2A_Tri *tris, int num_tris);
	void ComputeAngularRadius(PQP_REAL O[3][3], C2A_Tri *tris, int num_tris);


	C2A_BV();
	~C2A_BV();

	void     FitToTris(PQP_REAL O[3][3], C2A_Tri *tris, int num_tris, PQP_REAL comR[3]);
	void     FitToTris_Corner(PQP_REAL O[3][3], C2A_Tri *tris, int num_tris, PQP_REAL comR[3]);
};


int
C2A_BV_Overlap(PQP_REAL R[3][3], PQP_REAL T[3], C2A_BV *b1, C2A_BV *b2);

#if PQP_BV_TYPE & RSS_TYPE
PQP_REAL
C2A_BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], C2A_BV *b1, C2A_BV *b2, PQP_REAL S[3]);

// The closest points between the two rectangles not the two RSS!!!
PQP_REAL
C2A_BV_Distance(PQP_REAL R[3][3], PQP_REAL T[3], C2A_BV *b1, C2A_BV *b2,
				PQP_REAL P[3], PQP_REAL Q[3],
				bool &bValidPQ, PQP_REAL S[3]);
#endif

#endif



