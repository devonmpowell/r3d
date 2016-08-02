/*
 *
 *		rNd_unit_tests.c 
 *		
 *		Definitions for rNd unit tests. 
 *		
 *		Devon Powell
 *		31 August 2015
 *
 *		This program was prepared by Los Alamos National Security, LLC at Los Alamos National
 *		Laboratory (LANL) under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). 
 *		All rights in the program are reserved by the DOE and Los Alamos National Security, LLC.  
 *		Permission is granted to the public to copy and use this software without charge, provided that 
 *		this Notice and any statement of authorship are reproduced on all copies.  Neither the U.S. 
 *		Government nor LANS makes any warranty, express or implied, or assumes any liability 
 *		or responsibility for the use of this software.
 *
 */

#include "test_helpers.h"
#include "utest.h"
#include <math.h>

// numerical tolerances for pass/warn/fail tests
#define TOL_WARN 1.0e-8
#define TOL_FAIL 1.0e-4

// minimum volume allowed for test polyhedra 
#define MIN_VOL 1.0e-8

// number of trials per test
#define NUM_TRIALS 10

// for recursive tests, the stack capacity and 
// maximum recursion depth allowed
#define STACK_SIZE 64
#define MAX_DEPTH 6

// order of polynomial integration for all tests 
#define POLY_ORDER 0


// --                 unit tests                    -- //

void test_split_simplex_thru_centroid() {

	// a very basic sanity check. Splits a simplex through its centroid
	// and checks to see whether the two resulting volumes add up to equal the original

	// variables: the polyhedra and their moments
	rNd_int i;
	rNd_real verts[RND_DIM+1][RND_DIM];
	rNd_poly opoly, poly1, poly2;
	rNd_real om[1], m1[1], m2[1];

	// generate a random tet and clip plane
	rand_simplex_Nd(verts, MIN_VOL);
	rNd_init_simplex(&opoly, verts);
	rNd_plane splane = thru_cent_Nd(&opoly);

	// split the poly by making two copies of the original poly
	// and them clipping them against the same plane, with one
	// oriented oppositely
	poly1 = opoly;
	poly2 = opoly;
	rNd_clip(&poly1, &splane, 1);
	for(i = 0; i < RND_DIM; ++i) splane.n[i] *= -1; 
	splane.d *= -1;
	rNd_clip(&poly2, &splane, 1);

	// reduce the original and its two parts
	rNd_reduce(&opoly, om, POLY_ORDER);
	rNd_reduce(&poly1, m1, POLY_ORDER);
	rNd_reduce(&poly2, m2, POLY_ORDER);

	// check to make sure the two halves are good
	ASSERT_TRUE(rNd_is_good(&poly1));
	ASSERT_TRUE(rNd_is_good(&poly2));

	// make sure the sum of moments equals the original 
	ASSERT_EQ(om[0], m1[0] + m2[0], TOL_FAIL);
	EXPECT_EQ(om[0], m1[0] + m2[0], TOL_WARN);
	
	// make sure neither of the two resulting volumes is larger than the original
	// (within some tolerance)
	ASSERT_LT(m1[0], om[0]*(1.0 + TOL_FAIL));
	EXPECT_LT(m1[0], om[0]*(1.0 + TOL_WARN));
	ASSERT_LT(m2[0], om[0]*(1.0 + TOL_FAIL));
	EXPECT_LT(m2[0], om[0]*(1.0 + TOL_WARN));
}

void test_nonconvex() {


	// Initializes an inside-out simplex to test the reduction operation 
	// under parity inversion 

	// variables: the polyhedra and their moments
	rNd_int i;
	rNd_real tmp, tetvol;
	rNd_real verts[RND_DIM+1][RND_DIM];
	rNd_poly opoly;
	rNd_real om[1];

	// generate a random tet and clip plane
	rand_simplex_Nd(verts, MIN_VOL);
	rNd_init_simplex(&opoly, verts);

	// reduce the original
	rNd_reduce(&opoly, om, POLY_ORDER);
	tetvol = rNd_orient(verts);
	printf("Positive volume = %.5e, should be %.5e\n", om[0], tetvol);

	// make sure the sum of moments equals the original 
	ASSERT_EQ(om[0], tetvol, TOL_FAIL);
	EXPECT_EQ(om[0], tetvol, TOL_WARN);

	// flip the simplex and try to reduce again
	for(i = 0; i < RND_DIM; ++i) {
		tmp = verts[0][i];
		verts[0][i] = verts[1][i];
		verts[1][i] = tmp;
	}
	rNd_init_simplex(&opoly, verts);
	rNd_reduce(&opoly, om, POLY_ORDER);
	tetvol = rNd_orient(verts);
	printf("Inverted volume = %.5e, should be %.5e\n", om[0], tetvol);

	// make sure the sum of moments equals the original 
	ASSERT_EQ(om[0], tetvol, TOL_FAIL);
	EXPECT_EQ(om[0], tetvol, TOL_WARN);


}

void test_recursive_splitting_nondegenerate() {

	// recursively splits a polyhedron (starting from a simplex) through the centroid,
	// checking to see that the resulting volumes add up properly.

	// explicit stack-based implementation
	rNd_int nstack, depth, t, i;
	rNd_poly polystack[STACK_SIZE];
	rNd_int depthstack[STACK_SIZE];

	// variables: the polyhedra and their moments
	rNd_real verts[RND_DIM+1][RND_DIM];
	rNd_poly opoly, poly1, poly2;
	rNd_plane splane;
	rNd_real om[1], m1[1], m2[1];

	// do many trials
	printf("Recursively splitting %d simplices, maximum splits per simplex is %d.\n", NUM_TRIALS, MAX_DEPTH);
	for(t = 0; t < NUM_TRIALS; ++t) {

		// generate a random simplex
		rand_simplex_Nd(verts, MIN_VOL);
		rNd_init_simplex(&opoly, verts);
	
		// push the starting simplex to the stack
		nstack = 0;
		polystack[nstack] = opoly;
		depthstack[nstack] = 0;
		++nstack;	
	
		// recursively split the poly
		while(nstack > 0) {
	
			// pop the stack
			--nstack;
			opoly = polystack[nstack];
			depth = depthstack[nstack];
	
			// generate a randomly oriented plane
			// through the centroid of the poly
			splane = thru_cent_Nd(&opoly);

			// split the poly by making two copies of the original poly
			// and them clipping them against the same plane, with one
			// oriented oppositely
			poly1 = opoly;
			poly2 = opoly;
			rNd_clip(&poly1, &splane, 1);
			for(i = 0; i < RND_DIM; ++i) splane.n[i] *= -1; 
			splane.d *= -1;
			rNd_clip(&poly2, &splane, 1);
		
			// reduce the original and its two parts
			rNd_reduce(&opoly, om, POLY_ORDER);
			rNd_reduce(&poly1, m1, POLY_ORDER);
			rNd_reduce(&poly2, m2, POLY_ORDER);
			
			// check to make sure the two halves are good
			// THIS IS SLOW
			//ASSERT_TRUE(rNd_is_good(&poly1));
			//ASSERT_TRUE(rNd_is_good(&poly2));
		
			// make sure the sum of moments equals the original 
			ASSERT_EQ(om[0], m1[0] + m2[0], TOL_FAIL);
			EXPECT_EQ(om[0], m1[0] + m2[0], TOL_WARN);
			
			// make sure neither of the two resulting volumes is larger than the original
			// (within some tolerance)
			ASSERT_LT(m1[0], om[0]*(1.0 + TOL_FAIL));
			EXPECT_LT(m1[0], om[0]*(1.0 + TOL_WARN));
			ASSERT_LT(m2[0], om[0]*(1.0 + TOL_FAIL));
			EXPECT_LT(m2[0], om[0]*(1.0 + TOL_WARN));

			//printf("nstack = %d, depth = %d, opoly = %.10e, p1 = %.10e, p2 = %.10e, err = %.10e\n", 
					//nstack, depth, om[0], m1[0], m2[0], fabs(1.0 - om[0]/(m1[0] + m2[0])));
	
			// push the children to the stack if they have
			// an acceptably large volume
			if(depth < MAX_DEPTH) {
				if(m1[0] > MIN_VOL) {
					polystack[nstack] = poly1;
					depthstack[nstack] = depth + 1;
					++nstack;	
				}
				if(m2[0] > MIN_VOL) {
					polystack[nstack] = poly2;
					depthstack[nstack] = depth + 1;
					++nstack;	
				}
			}
		}
	}
}


void test_moments() {


	// Initializes an inside-out simplex to test the reduction operation 
	// under parity inversion 
	
	// TODO: just volume for now

	// variables: the polyhedra and their moments
	rNd_int i;
	rNd_real tmp, tetvol;
	rNd_real verts[RND_DIM+1][RND_DIM];
	rNd_poly opoly;
	rNd_real om[1];

	// generate a random tet and clip plane
	rand_simplex_Nd(verts, MIN_VOL);
	rNd_init_simplex(&opoly, verts);
	tetvol = rNd_orient(verts);

#undef NUM_TRIALS
#define NUM_TRIALS 1000
	r3d_int trial;
	printf("Computing moments of order %d, %d trials.\n", POLY_ORDER, NUM_TRIALS);
	for(trial = 0; trial < NUM_TRIALS; ++trial) {

		// reduce the original
		rNd_reduce(&opoly, om, POLY_ORDER);
		//printf("Positive volume = %.5e, should be %.5e\n", om[0], tetvol);
	
		// make sure the sum of moments equals the original 
		ASSERT_EQ(om[0], tetvol, TOL_FAIL);
		EXPECT_EQ(om[0], tetvol, TOL_WARN);

	}
}







// -- user-implemented functions declared in utest.h -- //

void register_all_tests() {

	register_test(test_split_simplex_thru_centroid, "split_simplex_thru_centroid");
	//register_test(test_nonconvex, "nonconvex");
	register_test(test_recursive_splitting_nondegenerate, "recursive_splitting_nondegenerate");
	register_test(test_moments, "moments");

}

void setup() {

	// no print buffer
	setbuf(stdout, NULL);

	// random number seed
	srand(time(NULL));
	//srand(10291986);
}

/////////////////////////////////////////////////////////

