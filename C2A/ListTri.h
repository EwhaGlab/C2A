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
                       tangmin@ewha.ac.kr


\**************************************************************************/
#ifndef LIST_H
#define LIST_H

#include <iostream>
#include <list>
#include <string>
#include <algorithm>

void ContactF:: IcI(int a[3], int b[3])
{
	a[0]=b[0];
	a[1]=b[1];
	a[2]=b[2];
}

void ContactF:: VcV_L(PQP_REAL a[3],PQP_REAL b[3])
{
	a[0]=b[0];
	a[1]=b[1];
	a[2]=b[2];
}

ContactF::ContactF(int FeatureType_A, int FeatureType_B, int FeatureID_A1[3], int FeatureID_B1[3],
int TriangleID_A,	int TriangleID_B,	PQP_REAL P_A1[3],	PQP_REAL P_B1[3],
PQP_REAL Distance) : FeatureType_A(FeatureType_A), FeatureType_B(FeatureType_B),
TriangleID_B(TriangleID_B),TriangleID_A(TriangleID_A),Distance(Distance) 
{ 
	IcI(FeatureID_A,FeatureID_A1);
    IcI(FeatureID_B,FeatureID_B1);
	VcV_L(P_A,P_A1);
	VcV_L(P_B,P_B1);
}
bool operator > (const ContactF& T_A,const ContactF& T_B)
{
	return(T_A.Distance > T_B.Distance);
}

bool operator < (const ContactF& T_A,const ContactF& T_B)
{
	return(T_A.Distance < T_B.Distance);
}


#endif









































