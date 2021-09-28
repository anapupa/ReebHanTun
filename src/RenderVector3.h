/*
(c) 2012 Fengtao Fan
*/
#ifndef GFXMATH_VEC3_INCLUDED // -*- C++ -*-
#define GFXMATH_VEC3_INCLUDED

//
// Define the RealTypeForVector3 (ie. default floating point) type
//
//#ifdef GFX_REAL_FLOAT
//typedef float RealTypeForVector3;
//#else
//#define GFX_REAL_DOUBLE
//typedef double RealTypeForVector3;
//#endif
typedef double RealTypeForVector3;
//
#ifndef GFX_NO_AXIS_NAMES
enum Axis {
    X = 0, Y = 1, Z = 2, W = 3
};
#endif
//
#ifndef FEQ_EPS
#define FEQ_EPS 1e-6
#define FEQ_EPS2 1e-12
#endif

// Handle platforms which don't define M_PI in <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/************************************************************************

  3D Vector class.

  $Id: Vector3.h,v 1.6 1997/03/17 22:52:26 garland Exp $

 ************************************************************************/
//#include "std.h"
#include <math.h>
#include <iostream>

using namespace std;

class Vector3 {
private:
    RealTypeForVector3 elt[3];

protected:
    inline void copy(const Vector3 &v);

public:
    //
    // Standard constructors
    //
    Vector3(RealTypeForVector3 x = 0, RealTypeForVector3 y = 0, RealTypeForVector3 z = 0) {
        elt[0] = x;
        elt[1] = y;
        elt[2] = z;
    }

#ifdef GFXMATH_VEC2_INCLUDED
    Vector3(const Vec2& v, RealTypeForVector3 z) { elt[0]=v[0]; elt[1]=v[1]; elt[2]=z; }
#endif

    Vector3(const Vector3 &v) { copy(v); }

    Vector3(const RealTypeForVector3 *v) {
        elt[0] = v[0];
        elt[1] = v[1];
        elt[2] = v[2];
    }

    //
    // Access methods
    //
#ifdef SAFETY
    RealTypeForVector3& operator()(int i)       { assert(i>=0 && i<3); return elt[i]; }
    RealTypeForVector3  operator()(int i) const { assert(i>=0 && i<3); return elt[i]; }
#else

    RealTypeForVector3 &operator()(int i) { return elt[i]; }

    RealTypeForVector3 operator()(int i) const { return elt[i]; }

#endif

    RealTypeForVector3 &operator[](int i) { return elt[i]; }

    RealTypeForVector3 operator[](int i) const { return elt[i]; }

    RealTypeForVector3 *raw() { return elt; }

    const RealTypeForVector3 *raw() const { return elt; }

    //
    // Comparison operators
    //
    inline bool operator==(const Vector3 &v) const;

    inline bool operator!=(const Vector3 &v) const;

    //
    // Assignment and in-place arithmetic methods
    //
    inline void set(RealTypeForVector3 x, RealTypeForVector3 y, RealTypeForVector3 z) {
        elt[0] = x;
        elt[1] = y;
        elt[2] = z;
    }

    inline Vector3 &operator=(const Vector3 &v);

    inline Vector3 &operator+=(const Vector3 &v);

    inline Vector3 &operator-=(const Vector3 &v);

    inline Vector3 &operator*=(RealTypeForVector3 s);

    inline Vector3 &operator/=(RealTypeForVector3 s);

    //
    // Binary arithmetic methods
    //
    inline Vector3 operator+(const Vector3 &v) const;

    inline Vector3 operator-(const Vector3 &v) const;

    inline Vector3 operator-() const;

    inline Vector3 operator*(RealTypeForVector3 s) const;

    inline Vector3 operator/(RealTypeForVector3 s) const;

    inline RealTypeForVector3 operator*(const Vector3 &v) const;

    inline Vector3 operator^(const Vector3 &v) const;

    //
    // set individual component
    //
    inline void setX(RealTypeForVector3 s) {
        elt[0] = s;
    }

    inline void setY(RealTypeForVector3 s) {
        elt[1] = s;
    }

    inline void setZ(RealTypeForVector3 s) {
        elt[2] = s;
    }

    inline void setI(unsigned int i, RealTypeForVector3 s) {
        elt[i] = s;
    }
};



////////////////////////////////////////////////////////////////////////
//
// Method definitions
//

inline void Vector3::copy(const Vector3 &v) {
    elt[0] = v.elt[0];
    elt[1] = v.elt[1];
    elt[2] = v.elt[2];
}

inline bool Vector3::operator==(const Vector3 &v) const {
    RealTypeForVector3 dx = elt[X] - v[X], dy = elt[Y] - v[Y], dz = elt[Z] - v[Z];
    return (dx * dx + dy * dy + dz * dz) < FEQ_EPS2;
}

inline bool Vector3::operator!=(const Vector3 &v) const {
    RealTypeForVector3 dx = elt[X] - v[X], dy = elt[Y] - v[Y], dz = elt[Z] - v[Z];
    return (dx * dx + dy * dy + dz * dz) > FEQ_EPS2;
}

inline Vector3 &Vector3::operator=(const Vector3 &v) {
    copy(v);
    return *this;
}

inline Vector3 &Vector3::operator+=(const Vector3 &v) {
    elt[0] += v[0];
    elt[1] += v[1];
    elt[2] += v[2];
    return *this;
}

inline Vector3 &Vector3::operator-=(const Vector3 &v) {
    elt[0] -= v[0];
    elt[1] -= v[1];
    elt[2] -= v[2];
    return *this;
}

inline Vector3 &Vector3::operator*=(RealTypeForVector3 s) {
    elt[0] *= s;
    elt[1] *= s;
    elt[2] *= s;
    return *this;
}

inline Vector3 &Vector3::operator/=(RealTypeForVector3 s) {
    elt[0] /= s;
    elt[1] /= s;
    elt[2] /= s;
    return *this;
}


inline Vector3 Vector3::operator+(const Vector3 &v) const {
    return Vector3(elt[0] + v[0], elt[1] + v[1], elt[2] + v[2]);
}

inline Vector3 Vector3::operator-(const Vector3 &v) const {
    return Vector3(elt[0] - v[0], elt[1] - v[1], elt[2] - v[2]);
}

inline Vector3 Vector3::operator-() const {
    return Vector3(-elt[0], -elt[1], -elt[2]);
}

inline Vector3 Vector3::operator*(RealTypeForVector3 s) const {
    return Vector3(elt[0] * s, elt[1] * s, elt[2] * s);
}

inline Vector3 Vector3::operator/(RealTypeForVector3 s) const {
    return Vector3(elt[0] / s, elt[1] / s, elt[2] / s);
}

inline RealTypeForVector3 Vector3::operator*(const Vector3 &v) const {
    return elt[0] * v[0] + elt[1] * v[1] + elt[2] * v[2];
}

inline Vector3 Vector3::operator^(const Vector3 &v) const {
    Vector3 w(elt[1] * v[2] - v[1] * elt[2],
              -elt[0] * v[2] + v[0] * elt[2],
              elt[0] * v[1] - v[0] * elt[1]);
    return w;
}

// Make scalar multiplication commutative
inline Vector3 operator*(RealTypeForVector3 s, const Vector3 &v) { return v * s; }



////////////////////////////////////////////////////////////////////////
//
// Primitive function definitions
//

inline RealTypeForVector3 norm(const Vector3 &v) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline RealTypeForVector3 norm2(const Vector3 &v) {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

inline RealTypeForVector3 length(const Vector3 &v) { return norm(v); }


inline RealTypeForVector3 unitize(Vector3 &v) {
    RealTypeForVector3 l = norm2(v);
    if (l != 1.0 && l != 0.0) {
        l = sqrt(l);
        v /= l;
    }
    return l;
}



////////////////////////////////////////////////////////////////////////
//
// Misc. function definitions
//

inline ostream &operator<<(ostream &out, const Vector3 &v) {
    return out << "[" << v[0] << " " << v[1] << " " << v[2] << "]";
}

inline istream &operator>>(istream &in, Vector3 &v) {
    //return in >> "[" >> v[0] >> v[1] >> v[2] >> "]";
    return in;
}

#ifdef GFXGL_INCLUDED
inline void glV(const Vector3& v) { glVertex(v[X], v[Y], v[Z]); }
inline void glN(const Vector3& v) { glNormal(v[X], v[Y], v[Z]); }
inline void glC(const Vector3& v) { glColor(v[X], v[Y], v[Z]); }
#endif


// GFXMATH_VEC3_INCLUDED
#endif
