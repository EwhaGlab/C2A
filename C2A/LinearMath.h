
#ifndef __LinearMath_h
#define __LinearMath_h

#include <math.h>
#ifdef WIN32
#include <windows.h>
#endif


typedef double Real;


// Forward declarations
class Matrix3x3;
class Transform;
class Quaternion;

//////////////////////////////////////////////////////////////////////////////
// Coord3D
//
// Description:
//      A class to represent a v3D ector.
//////////////////////////////////////////////////////////////////////////////
class Coord3D {
  public:
    inline Coord3D( ) {  }
    inline Coord3D( const Coord3D& t ) { *this = t; }
    inline Coord3D( Real x, Real y, Real z )
                                                    { Set_Value( x, y, z ); }
    inline Coord3D( Real t[] ) { Set_Value( t ); }
    inline Coord3D( const Quaternion& q );
    inline ~Coord3D( ) { }

    // Size of the data
    inline int Size( ) const { return 3; }

    // Get/set routines
    inline const Real& X( ) const { return val[0]; }
    inline const Real& Y( ) const { return val[1]; }
    inline const Real& Z( ) const { return val[2]; }

    inline Real& X( ) { return val[0]; }
    inline Real& Y( ) { return val[1]; }
    inline Real& Z( ) { return val[2]; }


	inline const Real& operator[](int i) const { return val[i%3]; }
	inline Real& operator[](int i) { return val[i%3]; }

    inline const Real* Value( ) const { return val; }
    inline void Get_Value( Real v[] ) const;
    inline void Get_Value( Real& x, Real& y, Real& z ) const;

    inline void Set_X( Real x ) { val[0] = x; }
    inline void Set_Y( Real y ) { val[1] = y; }
    inline void Set_Z( Real z ) { val[2] = z; }
    inline void Set_Value( const Real v[] );

    inline void Set_Value( Real x, Real y, Real z );



    // Distance functions
    inline Real Dist_Sq( const Coord3D& t ) const;
    // Distance to the line defined by the two points.
    inline Real Dist( const Coord3D& t ) const;
    inline Real Length_Sq( ) const;
    inline Real Length( ) const;

    // Min max functions
    // Take this triple and merge using min or max into the other triple(s)
    inline void Min( Coord3D& mn );
    inline void Max( Coord3D& mx );
    inline void Min_Max( Coord3D& mn, Coord3D& mx );

    // Operations on this object
    inline void Normalize( );
    inline void Zero( );
    inline void Identity( );
    inline void Negate( );
    inline void Add( Real v1, Real v2, Real v3 );
    inline void Sub( Real v1, Real v2, Real v3 );

    // Operators that change this object
    inline void operator+=( const Coord3D& t );
    inline void operator-=( const Coord3D& t );
    inline void operator*=( Real s );
    inline void operator/=( Real s );

    // Transformation by a row major matrix (right multiplication)
    // Point transformation (vector transformation is the same)
    inline void operator^=( const Matrix3x3& m );
    // Point transformation
    inline void operator^=( const Transform& m );
    // Vector transformation
    inline void operator&=( const Transform& m );

    inline bool operator==( const Coord3D& t ) const;
    inline bool operator!=( const Coord3D& t ) const;
    inline bool operator<=( const Coord3D& t ) const;
    inline bool operator>=( const Coord3D& t ) const;

    friend Coord3D operator-( const Coord3D& t );
    friend Coord3D operator+( const Coord3D&, const Coord3D& );
    friend Coord3D operator-( const Coord3D&, const Coord3D& );
    friend Coord3D operator*( Real s, const Coord3D& t );
    friend Coord3D operator*( const Coord3D& t, Real s );
    friend Coord3D operator/( const Coord3D& t, Real s );
    friend Real operator*( const Coord3D&, const Coord3D& );
    friend Coord3D operator%( const Coord3D&, const Coord3D& );
    friend Coord3D operator*( const Coord3D&, const Matrix3x3& );
    friend Coord3D operator*( const Matrix3x3&, const Coord3D& );
    // m^T * t
    friend Coord3D operator%( const Matrix3x3&, const Coord3D& );

    friend class Transform;
    // Point transformation
    friend Coord3D operator*( const Transform&,
                                   const Coord3D& );
    // Vector transformation
    friend Coord3D operator&( const Transform&,
                                   const Coord3D& );

  private:
    Real val[3];
};

//Zhangxy, Feb 12, 2005
class Quaternion {
  public:
    inline Quaternion( ) {   }
    inline Quaternion( const Quaternion& q ) { *this = q; }
    inline Quaternion( Real x, Real y, Real z, Real w)
                                                    { Set_Value( x, y, z, w ); }
    inline Quaternion( Real t[] ) { Set_Value( t ); }
    // Creates a quaternion [x,y,z,0]
    inline Quaternion( const Coord3D& t ) { Set_Value( t ); }
    // Axis angle constructor
    inline Quaternion( const Coord3D& t, Real theta )
                                                    { Set_Value( t, theta ); }
    // Euler angle constructor (first rotate about x, then y, then z)
    inline Quaternion( Real x, Real y, Real z )
                                                    { Set_Value( x, y, z ); }
    inline ~Quaternion( ) { }

    // Size of the data
    inline int Size( ) const { return 4; }
    
	// Get/set routines
    inline Real& X( )  { return val[0]; }
    inline Real& Y( )  { return val[1]; }
    inline Real& Z( )  { return val[2]; }
    inline Real& W( )  { return val[3]; }

    inline const Real& X( ) const { return val[0]; }
    inline const Real& Y( ) const { return val[1]; }
    inline const Real& Z( ) const { return val[2]; }
    inline const Real& W( ) const { return val[3]; }

    inline void Set_X( Real x ) { val[0] = x; }
    inline void Set_Y( Real y ) { val[1] = y; }
    inline void Set_Z( Real z ) { val[2] = z; }
    inline void Set_W( Real w ) { val[3] = w; }
    inline void Set_Value( const Real v[] );
    
	// Sets the quaternion to be [x,y,z,0]
    inline void Set_Value( const Coord3D& t );
    
	// Axis angle setter
    inline void Set_Value( const Coord3D& t, Real theta );
    
	// Euler angle setter (first rotate about x, then y', then z'')
    inline void Set_Value( Real x, Real y, Real z );

    inline const Real* Value( ) const { return val; }
    inline void Get_Value( Real v[] ) const
	{ v[0]=val[0], v[1]=val[1], v[2]=val[2], v[3]=val[3]; }
    inline void Get_Value( Real& x, Real& y, Real& z, Real& w ) const
	{ x=val[0], y=val[1], z=val[2], w=val[3]; }
    inline void Set_Value( Real x, Real y, Real z, Real w )
		{ val[0] = x, val[1] = y, val[2] = z, val[3] = w; }


    inline Real Length_Sq( ) const{return val[0]*val[0] + val[1]*val[1] + val[2]*val[2] + val[3]*val[3];}
    inline Real Length( ) const { return Real(sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2] + val[3]*val[3]));}

    // Operations on this object
    inline void Normalize( );
    inline void Invert( );
    inline void Zero( );
    inline void Identity( );
    inline void Negate( );

    inline Quaternion Inverse( ){return Quaternion(val[0], val[1], val[2], -val[3]);}
    inline Real Angle( ) {if (val[3] < -1.0 || val[3]> 1.0) return 0.0; return 2.0 * acos(val[3]);}
 
	
	// Operators that change this object
    inline void operator+=( const Quaternion& t );
    inline void operator-=( const Quaternion& t );
    inline void operator*=( const Quaternion& q );
    inline void operator*=( Real s );
    inline void operator/=( Real s );
    inline Real& operator[]( int i ){	
		switch(i)
		{
		case 0:return val[0];
		case 1:return val[1];
		case 2:return val[2];
		case 3:return val[3];
		default:
			return val[0];
		}
	}
    inline Real operator[]( int i ) const {	
		switch(i)
		{
		case 0:return val[0];
		case 1:return val[1];
		case 2:return val[2];
		case 3:return val[3];
		default:
			return val[0];
		}
	}

	Quaternion Farthest( const Quaternion& q) const{ 
		Quaternion diff,sum;
		diff = *this - q;
		sum = *this + q;
		if( diff*diff > sum*sum )
			return q;
		return (-q);
	}


    // Operators that do not change this object
    inline bool operator==( const Quaternion& t ) const;
    inline bool operator!=( const Quaternion& t ) const;
    inline bool operator<=( const Quaternion& t ) const;
    inline bool operator>=( const Quaternion& t ) const;

    friend Quaternion operator-( const Quaternion& q );
    friend Quaternion operator+( const Quaternion&, const Quaternion& );
    friend Quaternion operator-( const Quaternion&, const Quaternion& );
    friend Quaternion operator*( Real s, const Quaternion& t );
    friend Quaternion operator*( const Quaternion& t, Real s );
    friend Quaternion operator/( const Quaternion& t, Real s );
    friend Real		operator*( const Quaternion&, const Quaternion& );
    friend Quaternion operator%( const Quaternion&, const Quaternion& );

  private:
    Real val[4];
};

//////////////////////////////////////////////////////////////////////////////
// Matrix3x3
//
// Description:
//      A class to represent a matrix.
//////////////////////////////////////////////////////////////////////////////
class Matrix3x3 {
  public:
    inline Matrix3x3( ) { }
    inline Matrix3x3( const Matrix3x3& m ) { *this = m; }
    inline Matrix3x3( Real v[3][3] ) { Set_Value( v ); }
    inline Matrix3x3( Real v[] ) { Set_Value( v ); }
    inline Matrix3x3( Real v0, Real v1, Real v2,
                           Real v3, Real v4, Real v5,
                           Real v6, Real v7, Real v8 )
    {
        val[0] = v0, val[1] = v1, val[2] = v2, val[3] = v3, val[4] = v4;
        val[5] = v5, val[6] = v6, val[7] = v7, val[8] = v8;
    }
    inline ~Matrix3x3( ) { }

    // Size of the data
    inline int Size( ) const { return 9; }

    // Get/set routines
    inline void Get_Value( Real v[] ) const;
    inline void Get_Value( Real v[3][3] ) const;
    inline void Set_Value( const Real v[] );
    inline void Set_Value( const Real v[3][3] );

    inline const Real* Value( ) const { return val; }

	//inline Real* operator[](int i) const { return &(val[(i%3)*3]); }
	inline const Real* operator[](int i) const { return &(val[(i%3)*3]); }
	inline Real* operator[](int i)  { return &(val[(i%3)*3]); }

	inline bool IsZero() const { 
		for (int i = 0 ; i < 9 ; i++) 
			if (val[i] != 0.0) return false;
		return true;
	}

    inline void Set_Value( Real rotx, Real roty, Real rotz )
    {
        Real cos_a = cos( rotx );
        Real sin_a = sin( rotx );
        Matrix3x3 rx( 1.0, 0.0, 0.0, 0.0, cos_a, -sin_a,
                           0.0, sin_a, cos_a );
        cos_a = cos( roty );
        sin_a = sin( roty );
        Matrix3x3 ry( cos_a, 0.0, sin_a, 0.0, 1.0, 0.0,
                           -sin_a, 0.0, cos_a );
        cos_a = cos( rotz );
        sin_a = sin( rotz );
        Matrix3x3 rz( cos_a, -sin_a, 0.0, sin_a, cos_a, 0.0,
                     0.0, 0.0, 1.0 );
        *this = rz * (ry * rx);
    }
	
	//zhangxy Feb 11, 2006-------------
    inline Quaternion Quaternion_( );

    inline void Set_Value(const Quaternion& q );
    inline Matrix3x3 Inverse( ) const;
	Real Cofac(int r1, int c1, int r2, int c2) const{
		return val[r1*3+c1] * val[r2*3+c2] - val[r1*3+c2] * val[r2*3+c1];}
	//--------------zhangxy Feb 11, 2006

    // Operations that change this object
    inline void Zero( );
    inline void Identity( );
    //void Negate( );
    inline void Transpose( );
    //bool Invert( );

    // Operators that change this object
    inline void operator+=( const Matrix3x3& m );
    inline void operator-=( const Matrix3x3& m );
    inline void operator*=( Real s );
    inline void operator/=( Real s );

    // Operators that do not change this object
    //Real Determinant( ) const;

    friend void Coord3D::operator^=( const Matrix3x3& m );

    friend Matrix3x3 operator*( Real s, 
									const Matrix3x3& m );
    friend Matrix3x3 operator*( const Matrix3x3& m,
                                     Real s );
    friend Matrix3x3 operator/( const Matrix3x3& m,
                                     Real s );

    friend Coord3D operator*( const Coord3D& t,
                                   const Matrix3x3& m );
    friend Coord3D operator*( const Matrix3x3& m,
                                   const Coord3D& t );
    // m^T * t
    friend Coord3D operator%( const Matrix3x3& m,
                                   const Coord3D& t );

    friend Matrix3x3 operator-( const Matrix3x3& m );
    friend Matrix3x3 operator+( const Matrix3x3& m1,
                                     const Matrix3x3& m2 );
    friend Matrix3x3 operator-( const Matrix3x3& m1,
                                     const Matrix3x3& m2 );
    friend Matrix3x3 operator*( const Matrix3x3& m1,
                                     const Matrix3x3& m2 );
    // m1 * m2^T
    friend Matrix3x3 operator/( const Matrix3x3& m1,
                                     const Matrix3x3& m2 );
    // m1^T * m2
    friend Matrix3x3 operator%( const Matrix3x3& m1,
                                     const Matrix3x3& m2 );


    friend class Transform;
    friend void Coord3D::operator^=( const Transform& m );
    friend void Coord3D::operator&=( const Transform& m );
    friend Coord3D operator*( const Transform& m,
                                   const Coord3D& t );
    friend Coord3D operator&( const Transform& m,
                                   const Coord3D& t );


  private:
    Real val[9];
};

//////////////////////////////////////////////////////////////////////////////
// Transform
//
// Description:
//      A class to represent a rigid transformation.
//////////////////////////////////////////////////////////////////////////////
class Transform {
  public:
    inline Transform( ) { }
    // Set as a row major 3x4 matrix where the fourth column is the translation
    inline Transform( Real v[] ) { Set_Value( v ); }
    inline Transform( const Matrix3x3& matrix ) : R( matrix ), T( 0.0, 0.0, 0.0 ) { }
    inline Transform( const Matrix3x3& matrix,
					  const Coord3D& t ) : R(matrix), T(t) { }
	inline Transform(const Transform& val): R(val.R), T(val.T) {}
    // Set as a row major 3x4 matrix where the fourth column is the translation
    inline Transform( Real v0, Real v1, Real v2,
                      Real v3, Real v4, Real v5,
					  Real v6, Real v7, Real v8,
					  Real v9, Real v10, Real v11
            ) : R( v0, v1, v2, v4, v5, v6, v8, v9, v10 ), T( v3, v7, v11 ) { }

    inline ~Transform( ) { }

    // Size of the data
    inline int Size( ) const { return 12; }

    // Get/set routines
    inline const Coord3D& Translation() const { return T; }
	inline Coord3D& Translation() { return T; }
    inline const Matrix3x3& Rotation() const { return R; }
	inline Matrix3x3& Rotation() { return R; }
    inline const Quaternion Quaternion_()  { return R.Quaternion_(); }

	
    inline void Set_Value(const Transform& m ) { Set_Rotation(m.Rotation());
															Set_Translation(m.Translation());}//zhangxy, Feb 12, 2006
    inline void Set_Value( const Real v[] );
    inline void Set_Value( const Real R[], const Real T[] );
    inline void Set_Value( const Matrix3x3& r, const Coord3D& t ) { R = r; T = t; }
    inline void Set_Rotation( const Matrix3x3& r ) { R = r; }
    inline void Set_Rotation( const Quaternion& q ) { R.Set_Value(q); }//zhangxy, Feb 11, 2006

    inline void Set_Translation(const Coord3D& t ) { T = t; }

  // Operations that change this object

    inline void Invert( const Transform& m );
    inline void Identity( ) { R.Identity(); T.Identity(); }

    friend void Coord3D::operator^=( const Transform& m );
    friend void Coord3D::operator&=( const Transform& m );
    friend Coord3D operator*( const Transform& m,
                                   const Coord3D& t );
    friend Coord3D operator&( const Transform& m,
                                   const Coord3D& t );


  private:
    Matrix3x3 R;
    Coord3D T;
};

//////////////////////////////////////////////////////////////////////////////
// Auxiliary operator functions
//////////////////////////////////////////////////////////////////////////////

inline Coord3D operator-( const Coord3D& t );
inline Coord3D operator+( const Coord3D& t1, const Coord3D& t2 );
inline Coord3D operator-( const Coord3D& t1, const Coord3D& t2 );
inline Coord3D operator*( Real s, const Coord3D& t );
inline Coord3D operator*( const Coord3D& t, Real s );
inline Coord3D operator/( const Coord3D& t, Real s );
// Dot product
inline Real operator*( const Coord3D& t1, const Coord3D& t2 );
// Cross product
inline Coord3D operator%( const Coord3D& t1, const Coord3D& t2 );
// Vector-Matrix multiplication
inline Coord3D operator*( const Coord3D& t, const Matrix3x3& m );
// Matrix-Vector multiplication
inline Coord3D operator*( const Matrix3x3& m, const Coord3D& t );
// Matrix transpose-Vector multiplication
inline Coord3D operator%( const Matrix3x3& m, const Coord3D& t );


inline Matrix3x3 operator-( const Matrix3x3& m );
inline Matrix3x3 operator+( const Matrix3x3& m1,
                                 const Matrix3x3& m2 );
inline Matrix3x3 operator-( const Matrix3x3& m1,
                                 const Matrix3x3& m2 );
inline Matrix3x3 operator*( Real s, const Matrix3x3& m );
inline Matrix3x3 operator*( const Matrix3x3& m, Real s );
inline Matrix3x3 operator/( const Matrix3x3& m, Real s );
// Matrix multiplication = m1 * m2
inline Matrix3x3 operator*( const Matrix3x3& m1,
                                 const Matrix3x3& m2 );
// Matrix multiplication with transpose of second = m1 * m2^T
inline Matrix3x3 operator/( const Matrix3x3& m1,
                                 const Matrix3x3& m2 );
// Matrix multiplication with transpose of first = m1^T * m2
inline Matrix3x3 operator%( const Matrix3x3& m1, 
                                 const Matrix3x3& m2 );

// Point transformation
inline Coord3D operator*( const Transform& m,
                               const Coord3D& t );
// Vector transformation
inline Coord3D operator&( const Transform& m,
                               const Coord3D& t );


//////////////////////////////////////////////////////////////////////////////
// Coord3D member functions inline definitions
//////////////////////////////////////////////////////////////////////////////



inline void Coord3D::Get_Value( Real v[] ) const
{ v[0] = val[0], v[1] = val[1], v[2] = val[2]; }

inline void Coord3D::Get_Value( Real& x, Real& y,
                                     Real& z ) const
{ x = val[0], y = val[1], z = val[2]; }


inline void Coord3D::Set_Value( const Real v[] )
{ val[0] = v[0], val[1] = v[1], val[2] = v[2]; }

inline void Coord3D::Set_Value( Real x, Real y, Real z )
{ val[0] = x, val[1] = y, val[2] = z; }


inline Real Coord3D::Dist_Sq( const Coord3D& t ) const
{
    Real x, y, z;
    t.Get_Value( x, y, z );
    x -= val[0];
    y -= val[1];
    z -= val[2];
    return x * x + y * y + z * z;
}

inline Real Coord3D::Dist( const Coord3D& t ) const
{
    Real x, y, z;
    t.Get_Value( x, y, z );
    x -= val[0];
    y -= val[1];
    z -= val[2];
    return sqrt( x * x + y * y + z * z );
}

inline Real Coord3D::Length_Sq( ) const
{ return val[0] * val[0] + val[1] * val[1] + val[2] * val[2]; }

inline Real Coord3D::Length( ) const
{ return sqrt( val[0] * val[0] + val[1] * val[1] + val[2] * val[2] ); }

#ifndef min
#define min(a,b) a>b?b:a
#define max(a,b) a>b?a:b
#endif

inline void Coord3D::Min( Coord3D& mn )
{
    mn.Set_X( min( mn.X(), val[0] ) );
    mn.Set_Y( min( mn.Y(), val[1] ) );
    mn.Set_Z( min( mn.Z(), val[2] ) );
}

inline void Coord3D::Max( Coord3D& mx )
{
    mx.Set_X( max( mx.X(), val[0] ) );
    mx.Set_Y( max( mx.Y(), val[1] ) );
    mx.Set_Z( max( mx.Z(), val[2] ) );
}

inline void Coord3D::Min_Max( Coord3D& mn, Coord3D& mx )
{ Min( mn ); Max( mx ); }

inline void Coord3D::Normalize( )
{ if (this->Length() > 1.0e-20) *this /= this->Length(); }

inline void Coord3D::Zero( )
{ val[0] = val[1] = val[2] = 0.0; }

inline void Coord3D::Identity( )
{ Zero(); }

inline void Coord3D::Negate( )
{ val[0] = -val[0]; val[1] = -val[1]; val[2] = -val[2]; }

inline void Coord3D::Add( Real v1, Real v2, Real v3 )
{ val[0] += v1; val[1] += v2; val[2] += v3; }

inline void Coord3D::Sub( Real v1, Real v2, Real v3 )
{ val[0] -= v1; val[1] -= v2; val[2] -= v3; }


inline void Coord3D::operator+=( const Coord3D& t )
{ val[0] += t.X(); val[1] += t.Y(); val[2] += t.Z(); }

inline void Coord3D::operator-=( const Coord3D& t )
{ val[0] -= t.X(); val[1] -= t.Y(); val[2] -= t.Z(); }

inline void Coord3D::operator*=( Real s )
{ val[0] *= s; val[1] *= s; val[2] *= s; }

inline void Coord3D::operator/=( Real s )
{ const Real _s = 1.0/s; val[0] *= _s; val[1] *= _s; val[2] *= _s; }

inline void Coord3D::operator^=( const Matrix3x3& m )
{
    Real temp0 = val[0]*m.val[0] + val[1]*m.val[1] + val[2]*m.val[2];
    Real temp1 = val[0]*m.val[3] + val[1]*m.val[4] + val[2]*m.val[5];
    val[2] = val[0]*m.val[6] + val[1]*m.val[7] + val[2]*m.val[8];
    val[0] = temp0;
    val[1] = temp1;
}

inline void Coord3D::operator^=( const Transform& m )
{
    Real temp0 = val[0]*m.R.val[0] + val[1]*m.R.val[1] +
                       val[2]*m.R.val[2] + m.T.val[0];
    Real temp1 = val[0]*m.R.val[3] + val[1]*m.R.val[4] +
                       val[2]*m.R.val[5] + m.T.val[1];
    val[2] = val[0]*m.R.val[6] + val[1]*m.R.val[7] + val[2]*m.R.val[8] +
                                                     m.T.val[2];
    val[0] = temp0;
    val[1] = temp1;
}

inline void Coord3D::operator&=( const Transform& m )
{
    Real temp0 = val[0]*m.R.val[0] + val[1]*m.R.val[1] +
                       val[2]*m.R.val[2];
    Real temp1 = val[0]*m.R.val[3] + val[1]*m.R.val[4] +
                       val[2]*m.R.val[5];
    val[2] = val[0]*m.R.val[6] + val[1]*m.R.val[7] + val[2]*m.R.val[8];
    val[0] = temp0;
    val[1] = temp1;
}

inline bool Coord3D::operator==( const Coord3D& t ) const
{ return val[0] == t.X() && val[1] == t.Y() && val[2] == t.Z(); }

inline bool Coord3D::operator!=( const Coord3D& t ) const
{ return val[0] != t.X() || val[1] != t.Y() || val[2] != t.Z(); }

inline bool Coord3D::operator<=( const Coord3D& t ) const
{ return val[0] <= t.X() && val[1] <= t.Y() && val[2] <= t.Z(); }

inline bool Coord3D::operator>=( const Coord3D& t ) const
{ return val[0] >= t.X() && val[1] >= t.Y() && val[2] >= t.Z(); }


//////////////////////////////////////////////////////////////////////////////
// Quaternion member functions inline definitions
//////////////////////////////////////////////////////////////////////////////

// Sets the quaternion to be [x,y,z,0]
inline void Quaternion::Set_Value( const Coord3D& t )
{ val[0] = t.X(), val[1] = t.Y(), val[2] = t.Z(), val[3] = 0.0; }

// Axis angle setter
inline void Quaternion::Set_Value( const Coord3D& t,
                                         Real theta )
{
    Coord3D unit_axis = t;
    unit_axis.Normalize();
    unit_axis *= sin( 0.5 * theta );
    unit_axis.Get_Value( val[0], val[1], val[2] );
    val[3] = cos( 0.5 * theta );
}

// Euler angle setter (first rotate about x, then y', then z'')
// (The rotations take place in local coordinates)
inline void Quaternion::Set_Value(
                                    Real x, Real y, Real z )
{
    Set_Value( sin( 0.5*x ), 0.0, 0.0, cos( 0.5*x ) );
    *this *= Quaternion( 0.0, sin( 0.5*y ), 0.0, cos( 0.5*y ) );
    *this *= Quaternion( 0.0, 0.0, sin( 0.5*z ), cos( 0.5*z ) );
}



inline void Quaternion::operator+=( const Quaternion& q )
{
    Real v1, v2, v3, v4;
    q.Get_Value( v1, v2, v3, v4 );
    val[0] += v1; val[1] += v2; val[2] += v3; val[3] += v4;
}

inline void Quaternion::operator-=( const Quaternion& q )
{
    Real v1, v2, v3, v4;
    q.Get_Value( v1, v2, v3, v4 );
    val[0] -= v1; val[1] -= v2; val[2] -= v3; val[3] -= v4;
}

inline void Quaternion::operator*=( Real s )
{ val[0] *= s; val[1] *= s; val[2] *= s; val[3] *= s; }

inline void Quaternion::operator*=( const Quaternion& q )
{
    Coord3D t1( *this );
    Coord3D t2( q );
    Coord3D vec_part = (t1 % t2) + val[3] * t2 + q.val[3] * t1;

    vec_part.Get_Value( val[0], val[1], val[2] );

    val[3] = val[3] * q.val[3] - (t1 * t2);
}

inline void Quaternion::operator/=( Real s )
{ const Real _s = 1.0/s;
  val[0] *= _s; val[1] *= _s; val[2] *= _s; val[3] *= _s; }

inline bool Quaternion::operator==( const Quaternion& q ) const
{
    Real v1, v2, v3, v4;
    q.Get_Value( v1, v2, v3, v4 );
    return val[0] == v1 && val[1] == v2 && val[2] == v3 && val[3] == v4;
}

inline Coord3D::Coord3D( const Quaternion& q )
{
    Real dummy;
    q.Get_Value( val[0], val[1], val[2], dummy );
}

inline void Quaternion::Zero( )
{ val[0] = val[1] = val[2] = val[3] = 0.0; }

inline void Quaternion::Invert( )
{ val[0] = -val[0], val[1] = -val[1], val[2] = -val[2]; }

inline void Quaternion::Negate( ) { Invert( ); }



//////////////////////////////////////////////////////////////////////////////
// Matrix3x3 member functions inline definitions
//////////////////////////////////////////////////////////////////////////////
inline void Matrix3x3::Get_Value( Real v[] ) const
{
    v[0] = val[0], v[1] = val[1], v[2] = val[2];
    v[3] = val[3], v[4] = val[4], v[5] = val[5];
    v[6] = val[6], v[7] = val[7], v[8] = val[8];
}

inline void Matrix3x3::Get_Value( Real v[3][3] ) const
{
    v[0][0] = val[0], v[0][1] = val[1], v[0][2] = val[2];
    v[1][0] = val[3], v[1][1] = val[4], v[1][2] = val[5];
    v[2][0] = val[6], v[2][1] = val[7], v[2][2] = val[8];
}

inline void Matrix3x3::Set_Value( const Real v[] )
{
    val[0] = v[0], val[1] = v[1], val[2] = v[2];
    val[3] = v[3], val[4] = v[4], val[5] = v[5];
    val[6] = v[6], val[7] = v[7], val[8] = v[8];
}

inline void Matrix3x3::Set_Value( const Real v[3][3] )
{
    val[0] = v[0][0], val[1] = v[0][1], val[2] = v[0][2];
    val[3] = v[1][0], val[4] = v[1][1], val[5] = v[1][2];
    val[6] = v[2][0], val[7] = v[2][1], val[8] = v[2][2];
}

inline Quaternion Matrix3x3::Quaternion_( )
{
	Quaternion q;
	Real trace = val[0] + val[4] + val[8];
	
	if (trace > 0.0) 
	{
		Real s = sqrt(trace + 1.0);
		q[3] = s * 0.5;
		if (s == 0.0) exit(0);
		s = 0.5 / s;
		
		q[0] = (val[7] - val[5]) * s;
		q[1] = (val[2] - val[6]) * s;
		q[2] = (val[3] - val[1]) * s;
	} 
	else 
	{
		int i = val[0] < val[4] ? 
			(val[4] < val[8] ? 2 : 1) :
			(val[0] < val[8] ? 2 : 0); 
		int j = (i + 1) % 3;  
		int k = (i + 2) % 3;
		
		Real s = sqrt(val[i*3+i] - val[j*3+j] - val[k*3+k] + 1.0);
		q[i] = s * 0.5;
		s = 0.5 / s;
		
		q[3] = (val[k*3+j] - val[j*3+k]) * s;
		q[j] = (val[j*3+i] + val[i*3+j]) * s;
		q[k] = (val[k*3+i] + val[i*3+k]) * s;
	}

	return q;
}


inline Matrix3x3 Matrix3x3::Inverse( ) const
{
	Coord3D co(Cofac(1, 1, 2, 2), Cofac(1, 2, 2, 0), Cofac(1, 0, 2, 1));
	Real det = val[0] * co.X() +  val[1] * co.Y() + val[2] * co.Z() ;
	//assert(det != double(0.0f));
	if (det < 1.0e-15) exit(0);
	Real s = Real(1.0f) / det;
	return Matrix3x3(co.X() * s, Cofac(0, 2, 2, 1) * s, Cofac(0, 1, 1, 2) * s,
						  co.Y() * s, Cofac(0, 0, 2, 2) * s, Cofac(0, 2, 1, 0) * s,
						  co.Z() * s, Cofac(0, 1, 2, 0) * s, Cofac(0, 0, 1, 1) * s);

}

inline void Matrix3x3::Set_Value(const Quaternion& q )
{
	Real d = q.Length_Sq();
	//assert(d != double(0.0));
	Real s = double(2.0) / d;
	Real xs = q[0] * s,   ys = q[1] * s,   zs = q[2] * s;
	Real wx = q[3] * xs,  wy = q[3] * ys,  wz = q[3] * zs;
	Real xx = q[0] * xs,  xy = q[0] * ys,  xz = q[0] * zs;
	Real yy = q[1] * ys,  yz = q[1] * zs,  zz = q[2] * zs;
	Real v[9];
	v[0]=1.0 - (yy + zz); 
	v[1]=xy - wz; 
	v[2]=xz + wy;
	v[3]=xy + wz;
	v[4]=1.0 - (xx + zz);
	v[5]=yz - wx;
	v[6]=xz - wy;
	v[7]=yz + wx;
	v[8]=1.0 - (xx + yy);

	Set_Value(v);

}


inline void Matrix3x3::Zero( )
{
    val[0] = 0.0; val[1] = 0.0; val[2] = 0.0;
    val[3] = 0.0; val[4] = 0.0; val[5] = 0.0;
    val[6] = 0.0; val[7] = 0.0; val[8] = 0.0;
}

inline void Matrix3x3::Identity( )
{
    val[0] = 1.0; val[1] = 0.0; val[2] = 0.0;
    val[3] = 0.0; val[4] = 1.0; val[5] = 0.0;
    val[6] = 0.0; val[7] = 0.0; val[8] = 1.0;
}

void Matrix3x3::Transpose( )
{
    Real temp;

    temp = val[1];
    val[1] = val[3];
    val[3] = temp;

    temp = val[2];
    val[2] = val[6];
    val[6] = temp;

    temp = val[5];
    val[5] = val[7];
    val[7] = temp;
}   

void Matrix3x3::operator+=( const Matrix3x3& m )
{
    Real v[9];
    m.Get_Value( v );
    val[0] += v[0]; val[1] += v[1]; val[2] += v[2];
    val[3] += v[3]; val[4] += v[4]; val[5] += v[5];
    val[6] += v[6]; val[7] += v[7]; val[8] += v[8];
}

void Matrix3x3::operator-=( const Matrix3x3& m )
{
    Real v[9];
    m.Get_Value( v );
    val[0] -= v[0]; val[1] -= v[1]; val[2] -= v[2];
    val[3] -= v[3]; val[4] -= v[4]; val[5] -= v[5];
    val[6] -= v[6]; val[7] -= v[7]; val[8] -= v[8];
}

void Matrix3x3::operator*=( Real s )
{
    val[0] *= s; val[1] *= s; val[2] *= s; val[3] *= s;
    val[4] *= s; val[5] *= s; val[6] *= s; val[7] *= s;
    val[8] *= s;
}

void Matrix3x3::operator/=( Real s )
{
    const Real _s = 1.0/s;
    val[0] *= _s; val[1] *= _s; val[2] *= _s; val[3] *= _s; val[4] *= _s;
    val[5] *= _s; val[6] *= _s; val[7] *= _s; val[8] *= _s;
}


inline void Transform::Invert( const Transform& m )
{
    R.val[0] = m.R.val[0];
    R.val[1] = m.R.val[3];
    R.val[2] = m.R.val[6];
    R.val[3] = m.R.val[1];
    R.val[4] = m.R.val[4];
    R.val[5] = m.R.val[7];
    R.val[6] = m.R.val[2];
    R.val[7] = m.R.val[5];
    R.val[8] = m.R.val[8];
    T.val[0] = -(m.R.val[0]*m.T.val[0] + m.R.val[3]*m.T.val[1] +
                                         m.R.val[6]*m.T.val[2]);
    T.val[1] = -(m.R.val[1]*m.T.val[0] + m.R.val[4]*m.T.val[1] +
                                         m.R.val[7]*m.T.val[2]);
    T.val[2] = -(m.R.val[2]*m.T.val[0] + m.R.val[5]*m.T.val[1] +
                                         m.R.val[8]*m.T.val[2]);
}

inline void Transform::Set_Value( const Real v[] )
{
    R.val[0] = v[0], R.val[1] = v[1], R.val[2] = v[2];
    R.val[3] = v[4], R.val[4] = v[5], R.val[5] = v[6];
    R.val[6] = v[8], R.val[7] = v[9], R.val[8] = v[10];
    T.val[0] = v[3];
    T.val[1] = v[7];
    T.val[2] = v[11];
}

inline void Transform::Set_Value( const Real r[],
                                             const Real t[] )
{
    R.val[0] = r[0], R.val[1] = r[1], R.val[2] = r[2];
    R.val[3] = r[3], R.val[4] = r[4], R.val[5] = r[5];
    R.val[6] = r[6], R.val[7] = r[7], R.val[8] = r[8];
    T.val[0] = t[0];
    T.val[1] = t[1];
    T.val[2] = t[2];
}


//////////////////////////////////////////////////////////////////////////////
// Auxiliary operator functions inline definitions
//////////////////////////////////////////////////////////////////////////////

// Triple negation
inline Coord3D operator-( const Coord3D& t )
{ return Coord3D( -t.val[0], -t.val[1], -t.val[2] ); }

// Triple addition
inline Coord3D operator+( const Coord3D& t1, const Coord3D& t2 )
{ return Coord3D( t1.val[0] + t2.val[0], t1.val[1] + t2.val[1],
                       t1.val[2] + t2.val[2] ); }

// Triple subtraction
inline Coord3D operator-( const Coord3D& t1, const Coord3D& t2 )
{ return Coord3D( t1.val[0] - t2.val[0], t1.val[1] - t2.val[1],
                       t1.val[2] - t2.val[2] ); }

// Scalar multiplication
inline Coord3D operator*( Real s, const Coord3D& t )
{ return Coord3D( s * t.val[0], s * t.val[1], s * t.val[2] ); }

// Scalar multiplication
inline Coord3D operator*( const Coord3D& t, Real s )
{ return Coord3D( s * t.val[0], s * t.val[1], s * t.val[2] ); }

// Scalar division
inline Coord3D operator/( const Coord3D& t, Real s )
{ const Real _s = 1.0/s;
  return Coord3D( t.val[0] * _s, t.val[1] * _s, t.val[2] * _s ); }

// Dot product
inline Real operator*( const Coord3D& t1, const Coord3D& t2 )
{ return Real(t1.val[0]*t2.val[0] + t1.val[1]*t2.val[1] + t1.val[2]*t2.val[2]); }

// Cross product
inline Coord3D operator%( const Coord3D& t1, const Coord3D& t2 )
{
    return Coord3D( t1.val[1] * t2.val[2] - t1.val[2] * t2.val[1],
                         t1.val[2] * t2.val[0] - t1.val[0] * t2.val[2],
                         t1.val[0] * t2.val[1] - t1.val[1] * t2.val[0] );
}

inline Quaternion operator-( const Quaternion& q )
{
	return Quaternion( -q.val[0], -q.val[1], -q.val[2], -q.val[3] );
}
inline Quaternion operator+( const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion( q1.val[0] + q2.val[0], q1.val[1] + q2.val[1],
						   q1.val[2] + q2.val[2], q1.val[3] + q2.val[3] );
}
inline Quaternion operator-( const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion( q1.val[0] - q2.val[0], q1.val[1] - q2.val[1],
						   q1.val[2] - q2.val[2], q1.val[3] - q2.val[3] );
}
inline Quaternion operator*( Real s, const Quaternion& q )
{
	return Quaternion( s * q.val[0], s * q.val[1], s * q.val[2], s * q.val[3] );
}
inline Quaternion operator*( const Quaternion& q, Real s )
{
	return Quaternion( s * q.val[0], s * q.val[1], s * q.val[2], s * q.val[3] );
}
inline Quaternion operator/( const Quaternion& q, Real s )
{ 
	const Real _s = 1.0/s;
	return Quaternion( q.val[0] * _s, q.val[1] * _s, q.val[2] * _s, q.val[3] * _s ); 
}
// Quaternion product
inline Real operator*( const Quaternion& q1, const Quaternion& q2)
{
	return  q1.val[0] * q2.val[0] + q1.val[1] * q2.val[1] +
			q1.val[2] * q2.val[2] + q1.val[3] * q2.val[3] ;
}
// Refer to http://mathworld.wolfram.com/Quaternion.html
// q=(v, s)=(a1, a2, a3, s) where v=(a1, a2, a3)
// q1q2=(s1v2+s2v1+v1xv2, s1s2-v1v2)
inline Quaternion operator%( const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion( 
		q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1],
		q1[3] * q2[1] + q1[1] * q2[3] + q1[2] * q2[0] - q1[0] * q2[2],
		q1[3] * q2[2] + q1[2] * q2[3] + q1[0] * q2[1] - q1[1] * q2[0],
		q1[3] * q2[3] - q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2]);
}


// Matrix negation
inline Matrix3x3 operator-( const Matrix3x3& m )
{
    return Matrix3x3( -m.val[0], -m.val[1], -m.val[2],
                           -m.val[3], -m.val[4], -m.val[5],
                           -m.val[6], -m.val[7], -m.val[8] );
}

// Matrix addition
inline Matrix3x3 operator+( const Matrix3x3& m1,
                                 const Matrix3x3& m2 )
{
    return Matrix3x3( m1.val[0] + m2.val[0],
                           m1.val[1] + m2.val[1],
                           m1.val[2] + m2.val[2],
                           m1.val[3] + m2.val[3],
                           m1.val[4] + m2.val[4],
                           m1.val[5] + m2.val[5],
                           m1.val[6] + m2.val[6],
                           m1.val[7] + m2.val[7],
                           m1.val[8] + m2.val[8] );
}

// Matrix addition
inline Matrix3x3 operator-( const Matrix3x3& m1,
                                 const Matrix3x3& m2 )
{
    return Matrix3x3( m1.val[0] - m2.val[0],
                           m1.val[1] - m2.val[1],
                           m1.val[2] - m2.val[2],
                           m1.val[3] - m2.val[3],
                           m1.val[4] - m2.val[4],
                           m1.val[5] - m2.val[5],
                           m1.val[6] - m2.val[6],
                           m1.val[7] - m2.val[7],
                           m1.val[8] - m2.val[8] );
}

// Scalar multiplication
inline Matrix3x3 operator*( Real s, const Matrix3x3& m )
{
    return Matrix3x3( s * m.val[0], s * m.val[1], s * m.val[2],
                           s * m.val[3], s * m.val[4], s * m.val[5],
                           s * m.val[6], s * m.val[7], s * m.val[8] );
}

// Scalar multiplication
inline Matrix3x3 operator*( const Matrix3x3& m, Real s )
{
    return Matrix3x3( s * m.val[0], s * m.val[1], s * m.val[2],
                           s * m.val[3], s * m.val[4], s * m.val[5],
                           s * m.val[6], s * m.val[7], s * m.val[8] );
}

// Scalar division
inline Matrix3x3 operator/( const Matrix3x3& m, Real s )
{
    s = 1.0/s;
    return Matrix3x3( s * m.val[0], s * m.val[1], s * m.val[2],
                           s * m.val[3], s * m.val[4], s * m.val[5],
                           s * m.val[6], s * m.val[7], s * m.val[8] );
}

// Matrix vector right multiply
inline Coord3D operator*( const Matrix3x3& m, const Coord3D& t )
{
    return Coord3D(
                   t.val[0]*m.val[0] + t.val[1]*m.val[1] + t.val[2]*m.val[2],
                   t.val[0]*m.val[3] + t.val[1]*m.val[4] + t.val[2]*m.val[5],
                   t.val[0]*m.val[6] + t.val[1]*m.val[7] + t.val[2]*m.val[8] );
}

// Matrix transpose-Vector multiplication
inline Coord3D operator%( const Matrix3x3& m, const Coord3D& t )
{
    return Coord3D(
                   t.val[0]*m.val[0] + t.val[1]*m.val[3] + t.val[2]*m.val[6],
                   t.val[0]*m.val[1] + t.val[1]*m.val[4] + t.val[2]*m.val[7],
                   t.val[0]*m.val[2] + t.val[1]*m.val[5] + t.val[2]*m.val[8] );
}

// Matrix multiplication = m1 * m2
inline Matrix3x3 operator*( const Matrix3x3& m1,
                                 const Matrix3x3& m2 )
{
    return Matrix3x3(
            m1.val[0]*m2.val[0] + m1.val[1]*m2.val[3] + m1.val[2]*m2.val[6],
            m1.val[0]*m2.val[1] + m1.val[1]*m2.val[4] + m1.val[2]*m2.val[7],
            m1.val[0]*m2.val[2] + m1.val[1]*m2.val[5] + m1.val[2]*m2.val[8],
            m1.val[3]*m2.val[0] + m1.val[4]*m2.val[3] + m1.val[5]*m2.val[6],
            m1.val[3]*m2.val[1] + m1.val[4]*m2.val[4] + m1.val[5]*m2.val[7],
            m1.val[3]*m2.val[2] + m1.val[4]*m2.val[5] + m1.val[5]*m2.val[8],
            m1.val[6]*m2.val[0] + m1.val[7]*m2.val[3] + m1.val[8]*m2.val[6],
            m1.val[6]*m2.val[1] + m1.val[7]*m2.val[4] + m1.val[8]*m2.val[7],
            m1.val[6]*m2.val[2] + m1.val[7]*m2.val[5] + m1.val[8]*m2.val[8] );
}

// Matrix multiplication with transpose of second = m1 * m2^T
inline Matrix3x3 operator/( const Matrix3x3& m1,
                                 const Matrix3x3& m2 )
{
    return Matrix3x3(
            m1.val[0]*m2.val[0] + m1.val[1]*m2.val[1] + m1.val[2]*m2.val[2],
            m1.val[0]*m2.val[3] + m1.val[1]*m2.val[4] + m1.val[2]*m2.val[5],
            m1.val[0]*m2.val[6] + m1.val[1]*m2.val[7] + m1.val[2]*m2.val[8],
            m1.val[3]*m2.val[0] + m1.val[4]*m2.val[1] + m1.val[5]*m2.val[2],
            m1.val[3]*m2.val[3] + m1.val[4]*m2.val[4] + m1.val[5]*m2.val[5],
            m1.val[3]*m2.val[6] + m1.val[4]*m2.val[7] + m1.val[5]*m2.val[8],
            m1.val[6]*m2.val[0] + m1.val[7]*m2.val[1] + m1.val[8]*m2.val[2],
            m1.val[6]*m2.val[3] + m1.val[7]*m2.val[4] + m1.val[8]*m2.val[5],
            m1.val[6]*m2.val[6] + m1.val[7]*m2.val[7] + m1.val[8]*m2.val[8] );
}

// Matrix multiplication with transpose of first = m1^T * m2
inline Matrix3x3 operator%( const Matrix3x3& m1,
                                 const Matrix3x3& m2 )
{
    return Matrix3x3(
            m1.val[0]*m2.val[0] + m1.val[3]*m2.val[3] + m1.val[6]*m2.val[6],
            m1.val[0]*m2.val[1] + m1.val[3]*m2.val[4] + m1.val[6]*m2.val[7],
            m1.val[0]*m2.val[2] + m1.val[3]*m2.val[5] + m1.val[6]*m2.val[8],
            m1.val[1]*m2.val[0] + m1.val[4]*m2.val[3] + m1.val[7]*m2.val[6],
            m1.val[1]*m2.val[1] + m1.val[4]*m2.val[4] + m1.val[7]*m2.val[7],
            m1.val[1]*m2.val[2] + m1.val[4]*m2.val[5] + m1.val[7]*m2.val[8],
            m1.val[2]*m2.val[0] + m1.val[5]*m2.val[3] + m1.val[8]*m2.val[6],
            m1.val[2]*m2.val[1] + m1.val[5]*m2.val[4] + m1.val[8]*m2.val[7],
            m1.val[2]*m2.val[2] + m1.val[5]*m2.val[5] + m1.val[8]*m2.val[8] );
}


// Matrix point right multiply
inline Coord3D operator*( const Transform& m,
                               const Coord3D& t )
{
    return Coord3D(
        t.val[0]*m.R.val[0] + t.val[1]*m.R.val[1] + t.val[2]*m.R.val[2] +
                                                    m.T.val[0],
        t.val[0]*m.R.val[3] + t.val[1]*m.R.val[4] + t.val[2]*m.R.val[5] +
                                                    m.T.val[1],
        t.val[0]*m.R.val[6] + t.val[1]*m.R.val[7] + t.val[2]*m.R.val[8] +
                                                    m.T.val[2] );
}

// Matrix vector right multiply
inline Coord3D operator&( const Transform& m,
                               const Coord3D& t )
{
    return Coord3D(
        t.val[0]*m.R.val[0] + t.val[1]*m.R.val[1] + t.val[2]*m.R.val[2],
        t.val[0]*m.R.val[3] + t.val[1]*m.R.val[4] + t.val[2]*m.R.val[5],
        t.val[0]*m.R.val[6] + t.val[1]*m.R.val[7] + t.val[2]*m.R.val[8] );
}

#endif


