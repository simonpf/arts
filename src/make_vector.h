/* The make_vector functions can be used to initialize
   VECTORS in the source code. Usage is as simple as:

   VECTOR x        = make_vector(33.2, 17.3, 3.8);

   which will generate a VECTOR with the desired arguments. At the
   moment, up to 20 arguments are possible. If more are needed this
   can be easily extended.

   History:
   SAB 02.08.1999 Created.  
   SAB 28.03.2000 Moved make_array to separate file. */

#ifndef make_vector_h
#define make_vector_h

#include "vecmat.h"

// For VECTORs:
// ------------

VECTOR make_vector();

VECTOR make_vector(const Numeric& a1);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17,
		   const Numeric& a18);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17,
		   const Numeric& a18,
		   const Numeric& a19);

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17,
		   const Numeric& a18,
		   const Numeric& a19,
		   const Numeric& a20);



// For ARRAYs:
// -----------

// template<class T>
// ARRAY<T> make_array();

// template<class T>
// ARRAY<T> make_array(const T& a1);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17,
// 		    const T& a18);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17,
// 		    const T& a18,
// 		    const T& a19);

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17,
// 		    const T& a18,
// 		    const T& a19,
// 		    const T& a20);


#endif  // make_vector_h
