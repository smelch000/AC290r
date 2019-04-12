/*

New type double3 useful to treat coordinates, velocities as vectors

It implements sum, difference, scalar product between vectors,
as well as product with a scalar (to be put ALWAYS on the righthand side).
  
I/O operators is implemented (as friend methods).

June, 2014, P. Miocchi

*/
#ifndef DOUBLE3_H
#define DOUBLE3_H
#include <sstream>
#include <fstream>
#include <cmath>

class double3 {
  

  public:
  double x, y, z;
  
  inline double3(){}
  inline double3(double var){x=y=z=var;}
  inline double3(double v[]){ x=v[0]; y=v[1]; z=v[2];}
  inline double3(double _x, double _y, double _z){ x=_x; y=_y; z=_z;}

  inline double3 & operator=(double3 const & val) {
    if (this != &val) {
      x = val.x;
      y = val.y;
      z = val.z;
    }
    return *this;
  }
  
  /** just put all components equal to the same scalar
  useful, e.g., to set a vector to 0 */
  
  inline double3 & operator=(double const & val) {x=y=z=val; return *this;}

 inline double3 & operator+=(double3 const & val) {
   x +=val.x;
   y +=val.y;
   z +=val.z;
   return *this;
 }


 inline double3 & operator-=(double3 const & val) {
   x -= val.x;
   y -= val.y;
   z -= val.z;
   return *this;
 }

inline double3 operator+(double3 const & val) const {
    return double3(x+val.x,y+val.y,z+val.z);
  }
inline double3 operator-(double3 const &  val) const {
    return double3(x-val.x,y-val.y,z-val.z);
}

  /** external product
   WARNING: double3 must be the "left" term, i.e.
   double3 vector; double scalar;
   double3 result = vector * scalar */
  inline double3 operator*(double const & val) const {
    return double3(x*val,y*val,z*val);
  }


  /** external division
   WARNING: double3 must be the "left" term, i.e.
   double3 vector; double scalar;
   double3 result = vector / scalar  */
   inline double3 operator/(double const & val) const {
     return double3(x/val,y/val,z/val);
   }

  /** Curl (internal) product
   e.g.: double3 vector1; double vector2;
   double3 result = vector1.Curl(vector2) //==> result=vector1 x vector2 */
  inline double3 Curl(double3 const & val) const {
    return double3(y*val.z-z*val.y,
                   z*val.x-x*val.z, x*val.y-y*val.x);  
  }

  /** Return the versor parallel and equiverse to 'this' */
  inline double3 Versor() const{
      double m=Modulus();
      return double3(x/m,y/m,z/m);
  }
  
  /** scalar product */
  inline double operator*(double3 const & val) const {
    return (x*val.x+y*val.y+z*val.z);
  }

  inline double Modulus() const {
    return sqrt(x*x+y*y+z*z);
  }
  
  inline double Modulus2() const {
    return (x*x+y*y+z*z);
  }
  
  inline void Normalize() {
      double m=Modulus();
      x/=m;y/=m;z/=m;
  }

  inline bool operator==(double3 const & val) const {
      return (x==val.x && y==val.y && z==val.z);
  }

  inline bool operator==(double const & val) const {
      return (x==val && y==val && z==val);
  }
  inline void GetComp(double & xx, double & yy, double & zz) {xx=x;yy=y;zz=z;}
  
  void from_string(std::string const & s, int skipwords=0, bool replaceDecimal=false);
  std::string to_string();

  friend std::istream & operator>> (std::istream & in, double3 & val);

  friend std::ostream & operator<< (std::ostream & out, double3 const & val);

};
#endif
