/*************************************************************************\

  Ewha


\**************************************************************************/

#ifndef INTERP_MOTION_H
#define INTERP_MOTION_H

#include "C2A/LinearMath.h"
#include <PQP_compile.h>
#include "C2A/c2a_bv.h"
#include <string>

#include <PQP.h>


enum GMP_INTERP_MODE
{
	GMP_IM_EULER   = 0,   // linear interpolation of the euler angles.
	GMP_IM_LINEAR  = 1,   // constant rotation velocity and constant translation velocity of the origin.
	GMP_IM_SCREW   = 2,   // screw motion
	GMP_IM_SLERP   = 3,   // SLERP motion
	GMP_IM_LINEAR_ONE_CONTACT = 4, // Constained motion (v constaint) 
	GMP_IM_LINEAR_TWO_CONTACTS = 5,
	GMP_IM_LINEAR_THR_CONTACTS = 6,
	GMP_IM_ALLMODE = 7,
};

class CScrewMotion
{
public:
	PQP_REAL S[3]; // the rotation axis
	PQP_REAL p[3]; // a point in the axis S
	PQP_REAL b;    // angular velocity
	PQP_REAL d;    // distance of motion
	bool pure_translation;

public:
	CScrewMotion(
		const PQP_REAL R0[3][3], 
		const PQP_REAL T0[3],						 
		const PQP_REAL R1[3][3], 
		const PQP_REAL T1[3]);
	CScrewMotion(
		const PQP_REAL q0[7], 
		const PQP_REAL q1[7]);
	CScrewMotion();

};

class CInterpMotion
{
public:
  
	CInterpMotion(
		GMP_INTERP_MODE itpMode, 
		const PQP_REAL R0[3][3], 
		const PQP_REAL T0[3],						 
		const PQP_REAL R1[3][3], 
		const PQP_REAL T1[3]);

	CInterpMotion(
		GMP_INTERP_MODE itpMode, 
		const PQP_REAL q0[7], 
		const PQP_REAL q1[7]);
	
	CInterpMotion();
	
	
	void SetSrcDst(const PQP_REAL q0[7], const PQP_REAL q1[7]);

	virtual ~CInterpMotion();

	virtual double computeTOC(PQP_REAL d, PQP_REAL r1);
	
	virtual void velocity(void);

	virtual void velocity_screw(const PQP_REAL R0[3][3],const PQP_REAL T0[3],const PQP_REAL R1[3][3],const PQP_REAL T1[3]);
	virtual void velocity_screw(const PQP_REAL q0[7],const PQP_REAL q1[7]);


	// Interpolation: the results are represented as a quaternion
	// (q[0],... q[4]) Rotation, (q[5], q[6], q[7]) 
	virtual bool integrate(const double dt, PQP_REAL qua[7]) = 0;
	bool integrate(const double dt, PQP_REAL R[3][3], PQP_REAL T[3]);



	//compute the TOC using directional motion bound, //use the C2A_BV's angular radius to compute TOC,
	virtual double computeTOC(
		PQP_REAL d, 
		PQP_REAL r1,
		PQP_REAL S[3]);

	//compute the TOC using directional motion bound, //use the C2A_BV's angular radius to compute TOC,
	virtual double computeTOC_MotionBound(
		PQP_REAL T[3],  
		PQP_REAL d, 
		C2A_BV *V,
		PQP_REAL N[3]);



	virtual bool CAonRSS(
		PQP_REAL r1[3][3], 
		PQP_REAL tt1[3],	  
		PQP_REAL r2[3][3],
		PQP_REAL tt2[3],  
		PQP_REAL R[3][3],
		PQP_REAL T[3], 
		C2A_BV *b1, 
		C2A_BV *b2, 
		PQP_REAL *mint,  
		PQP_REAL *distance);

	virtual bool  CAonNonAdjacentTriangles(
		PQP_REAL R[3][3], 
		PQP_REAL T[3],   
		PQP_REAL triA[3][3],
		PQP_REAL triB[3][3],  
		PQP_REAL *mint, 
		PQP_REAL *distance);

	Quaternion DeltaRt(Real t);    // R_t
	Quaternion AbsoluteRt(Real t); // R(t) = R_t * R(0)

	void LinearAngularVelocity(Coord3D &axis, Real &angVel);
  

public:
	Real ConservD(Real d); // ConservD = d - security_distance_ratio * m_toc_delta 
	Real m_toc_delta; // the d_delta for the CCD query

	GMP_INTERP_MODE m_itpMode;
	Transform transform; // huangxin 12/01/2007
	Transform transform_s; //huangxin 12/01/2007
	Transform transform_t; //huangxin 12/01/2007


	Quaternion quaternion_s; 
	Quaternion quaternion_t; 


	Coord3D m_rel_tra; // The amount translation between q0 and q1 

	Transform m_transform_t_inv; // the inverse of transform_t: avoid the comptuation of the inverse during each integration



	Coord3D cv; // linear translation  velocity

	Coord3D cv_screw;
	Coord3D m_axis; // the axis of the rotation, normalized
	Real m_angVel; // the angular velocity around that axis: counter clockwise
					   // it might be negative, which means the clockwise rotation
	Coord3D point_screw;

	CScrewMotion m_screw;
};


class CInterpMotion_Linear : public CInterpMotion
{
public:

	virtual ~CInterpMotion_Linear();

	virtual void velocity(void);

	virtual void velocity_screw(const PQP_REAL R0[3][3],const PQP_REAL T0[3],const PQP_REAL R1[3][3],const PQP_REAL T1[3]);
	virtual void velocity_screw(const PQP_REAL q0[7],const PQP_REAL q1[7]);

	virtual bool integrate(const double dt, PQP_REAL qua[7]);

	CInterpMotion_Linear(
		const PQP_REAL R0[3][3],
		const PQP_REAL T0[3],
		const PQP_REAL R1[3][3], 
		const PQP_REAL T1[3]);

	CInterpMotion_Linear(const PQP_REAL q0[7], const PQP_REAL q1[7]);

	CInterpMotion_Linear();


	virtual double computeTOC(PQP_REAL d, PQP_REAL r1);



	//use the angular radius computed from the nearest points to compute TOC, more efficient than using C2A_BV's angular radius in computeTOC
	virtual double computeTOC(
		PQP_REAL d, 
		PQP_REAL r1,
		PQP_REAL S[3]);
	
	//compute the TOC using directional motion bound, //use the C2A_BV's angular radius to compute TOC, 
	virtual double computeTOC_MotionBound(
		PQP_REAL T[3],
		PQP_REAL d, 
		C2A_BV *V, 
		PQP_REAL N[3]);//using motion bound of RSS 


	virtual bool CAonRSS(
		PQP_REAL r1[3][3], 
		PQP_REAL tt1[3],	  
		PQP_REAL r2[3][3],
		PQP_REAL tt2[3], 	  
		PQP_REAL R[3][3],
		PQP_REAL T[3], 	  
		C2A_BV *b1, 
		C2A_BV *b2,
		PQP_REAL *mint,  
		PQP_REAL *distance);

	virtual bool  CAonNonAdjacentTriangles(
		PQP_REAL R[3][3], 
		PQP_REAL T[3],   
		PQP_REAL triA[3][3],
		PQP_REAL triB[3][3],  
		PQP_REAL *mint, 
		PQP_REAL *distance);

};


#endif;