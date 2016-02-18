/*
 *  
 *		r3d.c
 *		
 *		See r3d.h for usage.
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

#include "r3d.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

// useful macros
#define ONE_THIRD 0.333333333333333333333333333333333333333333333333333333
#define ONE_SIXTH 0.16666666666666666666666666666666666666666666666666666667
#define dot(va, vb) (va.x*vb.x + va.y*vb.y + va.z*vb.z)
#define wav(va, wa, vb, wb, vr) {			\
	vr.x = (wa*va.x + wb*vb.x)/(wa + wb);	\
	vr.y = (wa*va.y + wb*vb.y)/(wa + wb);	\
	vr.z = (wa*va.z + wb*vb.z)/(wa + wb);	\
}
#define norm(v) {					\
	r3d_real tmplen = sqrt(dot(v, v));	\
	v.x /= (tmplen + 1.0e-299);		\
	v.y /= (tmplen + 1.0e-299);		\
	v.z /= (tmplen + 1.0e-299);		\
}

void r3d_clip(r3d_poly* poly, r3d_plane* planes, r3d_int nplanes) {

	// direct access to vertex buffer
	r3d_vertex* vertbuffer = poly->verts; 
	r3d_int* nverts = &poly->nverts; 
	if(*nverts <= 0) return;

	// variable declarations
	r3d_int v, p, np, onv, vcur, vnext, vstart, 
			pnext, numunclipped;

	// signed distances to the clipping plane
	r3d_real sdists[R3D_MAX_VERTS];
	r3d_real smin, smax;

	// for marking clipped vertices
	r3d_int clipped[R3D_MAX_VERTS];

	// loop over each clip plane
	for(p = 0; p < nplanes; ++p) {

		// calculate signed distances to the clip plane
		onv = *nverts;
		smin = 1.0e30;
		smax = -1.0e30;
		memset(&clipped, 0, sizeof(clipped));
		for(v = 0; v < onv; ++v) {
			sdists[v] = planes[p].d + dot(vertbuffer[v].pos, planes[p].n);
			if(sdists[v] < smin) smin = sdists[v];
			if(sdists[v] > smax) smax = sdists[v];
			if(sdists[v] < 0.0) clipped[v] = 1;
		}

		// skip this face if the poly lies entirely on one side of it 
		if(smin >= 0.0) continue;
		if(smax <= 0.0) {
			*nverts = 0;
			return;
		}

		// check all edges and insert new vertices on the bisected edges 
		for(vcur = 0; vcur < onv; ++vcur) {
			if(clipped[vcur]) continue;
			for(np = 0; np < 3; ++np) {
				vnext = vertbuffer[vcur].pnbrs[np];
				if(!clipped[vnext]) continue;
				vertbuffer[*nverts].pnbrs[0] = vcur;
				vertbuffer[vcur].pnbrs[np] = *nverts;
				wav(vertbuffer[vcur].pos, -sdists[vnext],
					vertbuffer[vnext].pos, sdists[vcur],
					vertbuffer[*nverts].pos);
				(*nverts)++;
			}
		}

		// for each new vert, search around the faces for its new neighbors
		// and doubly-link everything
		for(vstart = onv; vstart < *nverts; ++vstart) {
			vcur = vstart;
			vnext = vertbuffer[vcur].pnbrs[0];
			do {
				for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
				vcur = vnext;
				pnext = (np+1)%3;
				vnext = vertbuffer[vcur].pnbrs[pnext];
			} while(vcur < onv);
			vertbuffer[vstart].pnbrs[2] = vcur;
			vertbuffer[vcur].pnbrs[1] = vstart;
		}

		// go through and compress the vertex list, removing clipped verts
		// and re-indexing accordingly (reusing `clipped` to re-index everything)
		numunclipped = 0;
		for(v = 0; v < *nverts; ++v) {
			if(!clipped[v]) {
				vertbuffer[numunclipped] = vertbuffer[v];
				clipped[v] = numunclipped++;
			}
		}
		*nverts = numunclipped;
		for(v = 0; v < *nverts; ++v) 
			for(np = 0; np < 3; ++np)
				vertbuffer[v].pnbrs[np] = clipped[vertbuffer[v].pnbrs[np]];
	}
}

	
void r3d_reduce(r3d_poly* poly, r3d_real* moments, r3d_int polyorder) {

	// var declarations
	r3d_real sixv;
	r3d_int np, m, i, j, k, corder;
	r3d_int vstart, pstart, vcur, vnext, pnext;
	r3d_rvec3 v0, v1, v2; 

	// direct access to vertex buffer
	r3d_vertex* vertbuffer = poly->verts; 
	r3d_int* nverts = &poly->nverts; 

	// zero the moments
	for(m = 0; m < R3D_NUM_MOMENTS(polyorder); ++m)
		moments[m] = 0.0;

	if(*nverts <= 0) return;
	
	// for keeping track of which edges have been visited
	r3d_int emarks[*nverts][3];
	memset(&emarks, 0, sizeof(emarks));

	// Storage for coefficients
	// keep two layers of the pyramid of coefficients
	// Note: Uses twice as much space as needed, but indexing is faster this way
	r3d_int prevlayer = 0;
	r3d_int curlayer = 1;
	r3d_real S[polyorder+1][polyorder+1][2];
	r3d_real D[polyorder+1][polyorder+1][2];
	r3d_real C[polyorder+1][polyorder+1][2];

	// loop over all vertices to find the starting point for each face
	for(vstart = 0; vstart < *nverts; ++vstart)
	for(pstart = 0; pstart < 3; ++pstart) {
	
		// skip this face if we have marked it
		if(emarks[vstart][pstart]) continue;

		// initialize face looping
		pnext = pstart; 
		vcur = vstart; 
		emarks[vcur][pnext] = 1;
		vnext = vertbuffer[vcur].pnbrs[pnext];
		v0 = vertbuffer[vcur].pos;
	
		// move to the second edge
		for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
		vcur = vnext;
		pnext = (np+1)%3;
		emarks[vcur][pnext] = 1;
		vnext = vertbuffer[vcur].pnbrs[pnext];
	
		// make a triangle fan using edges
		// and first vertex
		while(vnext != vstart) {
	
			v2 = vertbuffer[vcur].pos;
			v1 = vertbuffer[vnext].pos;
	
			sixv = (-v2.x*v1.y*v0.z + v1.x*v2.y*v0.z + v2.x*v0.y*v1.z
				   	- v0.x*v2.y*v1.z - v1.x*v0.y*v2.z + v0.x*v1.y*v2.z); 

			// calculate the moments
			// using the fast recursive method of Koehl (2012)
			// essentially building a set of trinomial pyramids, one layer at a time

			// base case
			S[0][0][prevlayer] = 1.0;
			D[0][0][prevlayer] = 1.0;
			C[0][0][prevlayer] = 1.0;
			moments[0] += ONE_SIXTH*sixv;

			// build up successive polynomial orders
			for(corder = 1, m = 1; corder <= polyorder; ++corder) {
				for(i = corder; i >= 0; --i)
				for(j = corder - i; j >= 0; --j, ++m) {
					k = corder - i - j;
					C[i][j][curlayer] = 0; 
					D[i][j][curlayer] = 0;  
					S[i][j][curlayer] = 0;  
					if(i > 0) {
						C[i][j][curlayer] += v2.x*C[i-1][j][prevlayer];
						D[i][j][curlayer] += v1.x*D[i-1][j][prevlayer]; 
						S[i][j][curlayer] += v0.x*S[i-1][j][prevlayer]; 
					}
					if(j > 0) {
						C[i][j][curlayer] += v2.y*C[i][j-1][prevlayer];
						D[i][j][curlayer] += v1.y*D[i][j-1][prevlayer]; 
						S[i][j][curlayer] += v0.y*S[i][j-1][prevlayer]; 
					}
					if(k > 0) {
						C[i][j][curlayer] += v2.z*C[i][j][prevlayer]; 
						D[i][j][curlayer] += v1.z*D[i][j][prevlayer]; 
						S[i][j][curlayer] += v0.z*S[i][j][prevlayer]; 
					}
					D[i][j][curlayer] += C[i][j][curlayer]; 
					S[i][j][curlayer] += D[i][j][curlayer]; 
					moments[m] += sixv*S[i][j][curlayer];
				}
				curlayer = 1 - curlayer;
				prevlayer = 1 - prevlayer;
			}

			// move to the next edge
			for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
			vcur = vnext;
			pnext = (np+1)%3;
			emarks[vcur][pnext] = 1;
			vnext = vertbuffer[vcur].pnbrs[pnext];
		}
	}

	// reuse C to recursively compute the leading multinomial coefficients
	C[0][0][prevlayer] = 1.0;
	for(corder = 1, m = 1; corder <= polyorder; ++corder) {
		for(i = corder; i >= 0; --i)
		for(j = corder - i; j >= 0; --j, ++m) {
			k = corder - i - j;
			C[i][j][curlayer] = 0.0; 
			if(i > 0) C[i][j][curlayer] += C[i-1][j][prevlayer];
			if(j > 0) C[i][j][curlayer] += C[i][j-1][prevlayer];
			if(k > 0) C[i][j][curlayer] += C[i][j][prevlayer]; 
			moments[m] /= C[i][j][curlayer]*(corder+1)*(corder+2)*(corder+3);
		}
		curlayer = 1 - curlayer;
		prevlayer = 1 - prevlayer;
	}
}


r3d_int r3d_is_good(r3d_poly* poly) {

	r3d_int v, np, rcur;
	r3d_int nvstack;
	r3d_int va, vb, vc;
	r3d_int vct[R3D_MAX_VERTS];
	r3d_int stack[R3D_MAX_VERTS];
	r3d_int regions[R3D_MAX_VERTS];

	// direct access to vertex buffer
	r3d_vertex* vertbuffer = poly->verts; 
	r3d_int* nverts = &poly->nverts; 

	// consistency check
	memset(&vct, 0, sizeof(vct));
	for(v = 0; v < *nverts; ++v) {

		// return false if two vertices are connected by more than one edge
		// or if any edges are obviously invalid
		for(np = 0; np < 3; ++np) {
			if(vertbuffer[v].pnbrs[np] == vertbuffer[v].pnbrs[(np+1)%3]) {
				printf("Double edge.\n");
				return 0;
			}
			if(vertbuffer[v].pnbrs[np] >= *nverts) {
				printf("Bad pointer.\n");
				return 0;
			}
		}

		vct[vertbuffer[v].pnbrs[0]]++;
		vct[vertbuffer[v].pnbrs[1]]++;
		vct[vertbuffer[v].pnbrs[2]]++;
	}
	
	// return false if any vertices are pointed to 
	// by more or fewer than three other vertices
	for(v = 0; v < *nverts; ++v) if(vct[v] != 3) {
		printf("Bad edge count: count[%d] = %d.\n", v, vct[v]);
		return 0;
	}

	// check for 3-vertex-connectedness
	// this is O(nverts^2)

	// handle multiply-connected polyhedra by testing each 
	// component separately. Flood-fill starting from each vertex
	// to give each connected region a unique ID.
	rcur = 1;
	memset(&regions, 0, sizeof(regions));
	for(v = 0; v < *nverts; ++v) {
		if(regions[v]) continue;
		nvstack = 0;
		stack[nvstack++] = v;
		while(nvstack > 0) {
			vc = stack[--nvstack];
			if(regions[vc]) continue;
			regions[vc] = rcur;
			stack[nvstack++] = vertbuffer[vc].pnbrs[0];
			stack[nvstack++] = vertbuffer[vc].pnbrs[1];
			stack[nvstack++] = vertbuffer[vc].pnbrs[2];
		}
		++rcur;
	}

	// loop over unique pairs of verts
	for(va = 0; va < *nverts; ++va) {
		rcur = regions[va];
		for(vb = va+1; vb < *nverts; ++vb) {
	
			// make sure va and vb are in the same connected component
			if(regions[vb] != rcur) continue;
	
			// pick vc != va && vc != vb 
			// and in the same connected component as va and vb
			for(vc = 0; vc < *nverts; ++vc)
			   if(regions[vc] == rcur && vc != va && vc != vb) break;
	
			// use vct to mark visited verts
			// mask out va and vb
			memset(&vct, 0, sizeof(vct));
			vct[va] = 1;
			vct[vb] = 1;
			
			// flood-fill from vc to make sure the graph is 
			// still connected when va and vb are masked
			nvstack = 0;
			stack[nvstack++] = vc;
			while(nvstack > 0) {
				vc = stack[--nvstack];
				if(vct[vc]) continue;
				vct[vc] = 1;
				stack[nvstack++] = vertbuffer[vc].pnbrs[0];
				stack[nvstack++] = vertbuffer[vc].pnbrs[1];
				stack[nvstack++] = vertbuffer[vc].pnbrs[2];
			}
	
			// if any verts in the region rcur were untouched, 
			// the graph is only 2-vertex-connected and hence an invalid polyhedron
			for(v = 0; v < *nverts; ++v) if(regions[v] == rcur && !vct[v]) {
				printf("Not 3-vertex-connected.\n");
				return 0;
			}
		}
	}	

	return 1;
}

void r3d_rotate(r3d_poly* poly, r3d_real theta, r3d_int axis) {
	r3d_int v;
	r3d_rvec3 tmp;
	r3d_real sine = sin(theta);
	r3d_real cosine = cos(theta);
	for(v = 0; v < poly->nverts; ++v) {
		tmp = poly->verts[v].pos;
		poly->verts[v].pos.xyz[(axis+1)%3] = cosine*tmp.xyz[(axis+1)%3] - sine*tmp.xyz[(axis+2)%3]; 
		poly->verts[v].pos.xyz[(axis+2)%3] = sine*tmp.xyz[(axis+1)%3] + cosine*tmp.xyz[(axis+2)%3]; 
	}
}

void r3d_translate(r3d_poly* poly, r3d_rvec3 shift) {
	r3d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		poly->verts[v].pos.x += shift.x;
		poly->verts[v].pos.y += shift.y;
		poly->verts[v].pos.z += shift.z;
	}
}

void r3d_scale(r3d_poly* poly, r3d_real scale) {
	r3d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		poly->verts[v].pos.x *= scale;
		poly->verts[v].pos.y *= scale;
		poly->verts[v].pos.z *= scale;
	}
}

void r3d_shear(r3d_poly* poly, r3d_real shear, r3d_int axb, r3d_int axs) {
	r3d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		poly->verts[v].pos.xyz[axb] += shear*poly->verts[v].pos.xyz[axs];
	}
}

void r3d_affine(r3d_poly* poly, r3d_real mat[4][4]) {
	r3d_int v;
	r3d_rvec3 tmp;
	r3d_real w;
	for(v = 0; v < poly->nverts; ++v) {
		tmp = poly->verts[v].pos;

		// affine transformation
		poly->verts[v].pos.x = tmp.x*mat[0][0] + tmp.y*mat[0][1] + tmp.z*mat[0][2] + mat[0][3];
		poly->verts[v].pos.y = tmp.x*mat[1][0] + tmp.y*mat[1][1] + tmp.z*mat[1][2] + mat[1][3];
		poly->verts[v].pos.z = tmp.x*mat[2][0] + tmp.y*mat[2][1] + tmp.z*mat[2][2] + mat[2][3];
		w = tmp.x*mat[3][0] + tmp.y*mat[3][1] + tmp.z*mat[3][2] + mat[3][3];
	
		// homogeneous divide if w != 1, i.e. in a perspective projection
		poly->verts[v].pos.x /= w;
		poly->verts[v].pos.y /= w;
		poly->verts[v].pos.z /= w;
	}
}


void r3d_init_tet(r3d_poly* poly, r3d_rvec3 verts[4]) {

	// direct access to vertex buffer
	r3d_vertex* vertbuffer = poly->verts; 
	r3d_int* nverts = &poly->nverts; 
	
	// initialize graph connectivity
	*nverts = 4;
	vertbuffer[0].pnbrs[0] = 1;	
	vertbuffer[0].pnbrs[1] = 3;	
	vertbuffer[0].pnbrs[2] = 2;	
	vertbuffer[1].pnbrs[0] = 2;	
	vertbuffer[1].pnbrs[1] = 3;	
	vertbuffer[1].pnbrs[2] = 0;	
	vertbuffer[2].pnbrs[0] = 0;	
	vertbuffer[2].pnbrs[1] = 3;	
	vertbuffer[2].pnbrs[2] = 1;	
	vertbuffer[3].pnbrs[0] = 1;	
	vertbuffer[3].pnbrs[1] = 2;	
	vertbuffer[3].pnbrs[2] = 0;	

	// copy vertex coordinates
	r3d_int v;
	for(v = 0; v < 4; ++v) vertbuffer[v].pos = verts[v];

}


void r3d_init_box(r3d_poly* poly, r3d_rvec3 rbounds[2]) {

	// direct access to vertex buffer
	r3d_vertex* vertbuffer = poly->verts; 
	r3d_int* nverts = &poly->nverts; 
	
	*nverts = 8;
	vertbuffer[0].pnbrs[0] = 1;	
	vertbuffer[0].pnbrs[1] = 4;	
	vertbuffer[0].pnbrs[2] = 3;	
	vertbuffer[1].pnbrs[0] = 2;	
	vertbuffer[1].pnbrs[1] = 5;	
	vertbuffer[1].pnbrs[2] = 0;	
	vertbuffer[2].pnbrs[0] = 3;	
	vertbuffer[2].pnbrs[1] = 6;	
	vertbuffer[2].pnbrs[2] = 1;	
	vertbuffer[3].pnbrs[0] = 0;	
	vertbuffer[3].pnbrs[1] = 7;	
	vertbuffer[3].pnbrs[2] = 2;	
	vertbuffer[4].pnbrs[0] = 7;	
	vertbuffer[4].pnbrs[1] = 0;	
	vertbuffer[4].pnbrs[2] = 5;	
	vertbuffer[5].pnbrs[0] = 4;	
	vertbuffer[5].pnbrs[1] = 1;	
	vertbuffer[5].pnbrs[2] = 6;	
	vertbuffer[6].pnbrs[0] = 5;	
	vertbuffer[6].pnbrs[1] = 2;	
	vertbuffer[6].pnbrs[2] = 7;	
	vertbuffer[7].pnbrs[0] = 6;	
	vertbuffer[7].pnbrs[1] = 3;	
	vertbuffer[7].pnbrs[2] = 4;	
	vertbuffer[0].pos.x = rbounds[0].x; 
	vertbuffer[0].pos.y = rbounds[0].y; 
	vertbuffer[0].pos.z = rbounds[0].z; 
	vertbuffer[1].pos.x = rbounds[1].x; 
	vertbuffer[1].pos.y = rbounds[0].y; 
	vertbuffer[1].pos.z = rbounds[0].z; 
	vertbuffer[2].pos.x = rbounds[1].x; 
	vertbuffer[2].pos.y = rbounds[1].y; 
	vertbuffer[2].pos.z = rbounds[0].z; 
	vertbuffer[3].pos.x = rbounds[0].x; 
	vertbuffer[3].pos.y = rbounds[1].y; 
	vertbuffer[3].pos.z = rbounds[0].z; 
	vertbuffer[4].pos.x = rbounds[0].x; 
	vertbuffer[4].pos.y = rbounds[0].y; 
	vertbuffer[4].pos.z = rbounds[1].z; 
	vertbuffer[5].pos.x = rbounds[1].x; 
	vertbuffer[5].pos.y = rbounds[0].y; 
	vertbuffer[5].pos.z = rbounds[1].z; 
	vertbuffer[6].pos.x = rbounds[1].x; 
	vertbuffer[6].pos.y = rbounds[1].y; 
	vertbuffer[6].pos.z = rbounds[1].z; 
	vertbuffer[7].pos.x = rbounds[0].x; 
	vertbuffer[7].pos.y = rbounds[1].y; 
	vertbuffer[7].pos.z = rbounds[1].z; 

}

void r3d_init_poly(r3d_poly* poly, r3d_rvec3* vertices, r3d_int numverts, 
					r3d_int** faceinds, r3d_int* numvertsperface, r3d_int numfaces) {

	// dummy vars
	r3d_int v, vprev, vcur, vnext, f, np;

	// direct access to vertex buffer
	r3d_vertex* vertbuffer = poly->verts; 
	r3d_int* nverts = &poly->nverts; 

	// count up the number of faces per vertex
	// and act accordingly
	r3d_int eperv[R3D_MAX_VERTS];
	r3d_int minvperf = R3D_MAX_VERTS;
	r3d_int maxvperf = 0;
	memset(&eperv, 0, sizeof(eperv));
	for(f = 0; f < numfaces; ++f)
		for(v = 0; v < numvertsperface[f]; ++v)
			++eperv[faceinds[f][v]];
	for(v = 0; v < numverts; ++v) {
		if(eperv[v] < minvperf) minvperf = eperv[v];
		if(eperv[v] > maxvperf) maxvperf = eperv[v];
	}

	// clear the poly
	*nverts = 0;

	// return if we were given an invalid poly
	if(minvperf < 3) return;

	if(maxvperf == 3) {

		// simple case with no need for duplicate vertices

		// read in vertex locations
		*nverts = numverts;
		for(v = 0; v < *nverts; ++v) {
			vertbuffer[v].pos = vertices[v];
			for(np = 0; np < 3; ++np) vertbuffer[v].pnbrs[np] = R3D_MAX_VERTS;
		}
	
		// build graph connectivity by correctly orienting half-edges for each vertex 
		for(f = 0; f < numfaces; ++f) {
			for(v = 0; v < numvertsperface[f]; ++v) {
				vprev = faceinds[f][v];
				vcur = faceinds[f][(v+1)%numvertsperface[f]];
				vnext = faceinds[f][(v+2)%numvertsperface[f]];
				for(np = 0; np < 3; ++np) {
					if(vertbuffer[vcur].pnbrs[np] == vprev) {
						vertbuffer[vcur].pnbrs[(np+2)%3] = vnext;
						break;
					}
					else if(vertbuffer[vcur].pnbrs[np] == vnext) {
						vertbuffer[vcur].pnbrs[(np+1)%3] = vprev;
						break;
					}
				}
				if(np == 3) {
					vertbuffer[vcur].pnbrs[1] = vprev;
					vertbuffer[vcur].pnbrs[0] = vnext;
				}
			}
		}
	}
	else {

		// we need to create duplicate, degenerate vertices to account for more than
		// three edges per vertex. This is complicated.

		r3d_int tface = 0;
		for(v = 0; v < numverts; ++v) tface += eperv[v];

		// need more variables
		r3d_int v0, v1, v00, v11, numunclipped;

		// we need a few extra buffers to handle the necessary operations
		r3d_vertex vbtmp[3*R3D_MAX_VERTS];
		r3d_int util[3*R3D_MAX_VERTS];
		r3d_int vstart[R3D_MAX_VERTS];

		// build vertex mappings to degenerate duplicates
		// and read in vertex locations
		*nverts = 0;
		for(v = 0; v < numverts; ++v) {
			vstart[v] = *nverts;
			for(vcur = 0; vcur < eperv[v]; ++vcur) {
				vbtmp[*nverts].pos = vertices[v];
				for(np = 0; np < 3; ++np) vbtmp[*nverts].pnbrs[np] = R3D_MAX_VERTS;
				++(*nverts);
			}	
		}

		// fill in connectivity for all duplicates
		memset(&util, 0, sizeof(util));
		for(f = 0; f < numfaces; ++f) {
			for(v = 0; v < numvertsperface[f]; ++v) {
				vprev = faceinds[f][v];
				vcur = faceinds[f][(v+1)%numvertsperface[f]];
				vnext = faceinds[f][(v+2)%numvertsperface[f]];
				vcur = vstart[vcur] + util[vcur]++;
				vbtmp[vcur].pnbrs[1] = vnext;
				vbtmp[vcur].pnbrs[2] = vprev;
			}
		}

		// link degenerate duplicates, putting them in the correct order
		// use util to mark and avoid double-processing verts
		memset(&util, 0, sizeof(util));
		for(v = 0; v < numverts; ++v) {
			for(v0 = vstart[v]; v0 < vstart[v] + eperv[v]; ++v0) {
				for(v1 = vstart[v]; v1 < vstart[v] + eperv[v]; ++v1) {
					if(vbtmp[v0].pnbrs[2] == vbtmp[v1].pnbrs[1] && !util[v0]) {
						vbtmp[v0].pnbrs[2] = v1;
						vbtmp[v1].pnbrs[0] = v0;
						util[v0] = 1;
					}
				}
			}
		}

		// complete vertex pairs
		memset(&util, 0, sizeof(util));
		for(v0 = 0; v0 < numverts; ++v0)
		for(v1 = v0 + 1; v1 < numverts; ++v1) {
			for(v00 = vstart[v0]; v00 < vstart[v0] + eperv[v0]; ++v00)
			for(v11 = vstart[v1]; v11 < vstart[v1] + eperv[v1]; ++v11) {
				if(vbtmp[v00].pnbrs[1] == v1 && vbtmp[v11].pnbrs[1] == v0 
						&& !util[v00] && !util[v11]) {
					vbtmp[v00].pnbrs[1] = v11;
					vbtmp[v11].pnbrs[1] = v00;
					util[v00] = 1;
					util[v11] = 1;
				}
			}
		}

		// remove unnecessary dummy vertices
		memset(&util, 0, sizeof(util));
		for(v = 0; v < numverts; ++v) {
			v0 = vstart[v];
			v1 = vbtmp[v0].pnbrs[0];
			v00 = vbtmp[v0].pnbrs[2];
			v11 = vbtmp[v1].pnbrs[0];
			vbtmp[v00].pnbrs[0] = vbtmp[v0].pnbrs[1];
			vbtmp[v11].pnbrs[2] = vbtmp[v1].pnbrs[1];
			for(np = 0; np < 3; ++np) if(vbtmp[vbtmp[v0].pnbrs[1]].pnbrs[np] == v0) break;
			vbtmp[vbtmp[v0].pnbrs[1]].pnbrs[np] = v00;
			for(np = 0; np < 3; ++np) if(vbtmp[vbtmp[v1].pnbrs[1]].pnbrs[np] == v1) break;
			vbtmp[vbtmp[v1].pnbrs[1]].pnbrs[np] = v11;
			util[v0] = 1;
			util[v1] = 1;
		}

		// copy to the real vertbuffer and compress
		numunclipped = 0;
		for(v = 0; v < *nverts; ++v) {
			if(!util[v]) {
				vertbuffer[numunclipped] = vbtmp[v];
				util[v] = numunclipped++;
			}
		}
		*nverts = numunclipped;
		for(v = 0; v < *nverts; ++v) 
			for(np = 0; np < 3; ++np)
				vertbuffer[v].pnbrs[np] = util[vertbuffer[v].pnbrs[np]];
	}
}

void r3d_tet_faces_from_verts(r3d_plane* faces, r3d_rvec3* verts) {
	r3d_rvec3 tmpcent;
	faces[0].n.x = ((verts[3].y - verts[1].y)*(verts[2].z - verts[1].z) 
			- (verts[2].y - verts[1].y)*(verts[3].z - verts[1].z));
	faces[0].n.y = ((verts[2].x - verts[1].x)*(verts[3].z - verts[1].z) 
			- (verts[3].x - verts[1].x)*(verts[2].z - verts[1].z));
	faces[0].n.z = ((verts[3].x - verts[1].x)*(verts[2].y - verts[1].y) 
			- (verts[2].x - verts[1].x)*(verts[3].y - verts[1].y));
	norm(faces[0].n);
	tmpcent.x = ONE_THIRD*(verts[1].x + verts[2].x + verts[3].x);
	tmpcent.y = ONE_THIRD*(verts[1].y + verts[2].y + verts[3].y);
	tmpcent.z = ONE_THIRD*(verts[1].z + verts[2].z + verts[3].z);
	faces[0].d = -dot(faces[0].n, tmpcent);

	faces[1].n.x = ((verts[2].y - verts[0].y)*(verts[3].z - verts[2].z) 
			- (verts[2].y - verts[3].y)*(verts[0].z - verts[2].z));
	faces[1].n.y = ((verts[3].x - verts[2].x)*(verts[2].z - verts[0].z) 
			- (verts[0].x - verts[2].x)*(verts[2].z - verts[3].z));
	faces[1].n.z = ((verts[2].x - verts[0].x)*(verts[3].y - verts[2].y) 
			- (verts[2].x - verts[3].x)*(verts[0].y - verts[2].y));
	norm(faces[1].n);
	tmpcent.x = ONE_THIRD*(verts[2].x + verts[3].x + verts[0].x);
	tmpcent.y = ONE_THIRD*(verts[2].y + verts[3].y + verts[0].y);
	tmpcent.z = ONE_THIRD*(verts[2].z + verts[3].z + verts[0].z);
	faces[1].d = -dot(faces[1].n, tmpcent);

	faces[2].n.x = ((verts[1].y - verts[3].y)*(verts[0].z - verts[3].z) 
			- (verts[0].y - verts[3].y)*(verts[1].z - verts[3].z));
	faces[2].n.y = ((verts[0].x - verts[3].x)*(verts[1].z - verts[3].z) 
			- (verts[1].x - verts[3].x)*(verts[0].z - verts[3].z));
	faces[2].n.z = ((verts[1].x - verts[3].x)*(verts[0].y - verts[3].y) 
			- (verts[0].x - verts[3].x)*(verts[1].y - verts[3].y));
	norm(faces[2].n);
	tmpcent.x = ONE_THIRD*(verts[3].x + verts[0].x + verts[1].x);
	tmpcent.y = ONE_THIRD*(verts[3].y + verts[0].y + verts[1].y);
	tmpcent.z = ONE_THIRD*(verts[3].z + verts[0].z + verts[1].z);
	faces[2].d = -dot(faces[2].n, tmpcent);

	faces[3].n.x = ((verts[0].y - verts[2].y)*(verts[1].z - verts[0].z) 
			- (verts[0].y - verts[1].y)*(verts[2].z - verts[0].z));
	faces[3].n.y = ((verts[1].x - verts[0].x)*(verts[0].z - verts[2].z) 
			- (verts[2].x - verts[0].x)*(verts[0].z - verts[1].z));
	faces[3].n.z = ((verts[0].x - verts[2].x)*(verts[1].y - verts[0].y) 
			- (verts[0].x - verts[1].x)*(verts[2].y - verts[0].y));
	norm(faces[3].n);
	tmpcent.x = ONE_THIRD*(verts[0].x + verts[1].x + verts[2].x);
	tmpcent.y = ONE_THIRD*(verts[0].y + verts[1].y + verts[2].y);
	tmpcent.z = ONE_THIRD*(verts[0].z + verts[1].z + verts[2].z);
	faces[3].d = -dot(faces[3].n, tmpcent);
}

void r3d_box_faces_from_verts(r3d_plane* faces, r3d_rvec3* rbounds) {
	faces[0].n.x = 0.0; faces[0].n.y = 0.0; faces[0].n.z = 1.0; faces[0].d = -rbounds[0].z; 
	faces[2].n.x = 0.0; faces[2].n.y = 1.0; faces[2].n.z = 0.0; faces[2].d = -rbounds[0].y; 
	faces[4].n.x = 1.0; faces[4].n.y = 0.0; faces[4].n.z = 0.0; faces[4].d = -rbounds[0].x; 
	faces[1].n.x = 0.0; faces[1].n.y = 0.0; faces[1].n.z = -1.0; faces[1].d = rbounds[1].z; 
	faces[3].n.x = 0.0; faces[3].n.y = -1.0; faces[3].n.z = 0.0; faces[3].d = rbounds[1].y; 
	faces[5].n.x = -1.0; faces[5].n.y = 0.0; faces[5].n.z = 0.0; faces[5].d = rbounds[1].x; 
}

void r3d_poly_faces_from_verts(r3d_plane* faces, r3d_rvec3* vertices, r3d_int numverts, 
						r3d_int** faceinds, r3d_int* numvertsperface, r3d_int numfaces) {

	// dummy vars
	r3d_int v, f;
	r3d_rvec3 p0, p1, p2, centroid;

	// calculate a centroid and a unit normal for each face 
	for(f = 0; f < numfaces; ++f) {

		centroid.x = 0.0;
		centroid.y = 0.0;
		centroid.z = 0.0;
		faces[f].n.x = 0.0;
		faces[f].n.y = 0.0;
		faces[f].n.z = 0.0;
		
		for(v = 0; v < numvertsperface[f]; ++v) {

			// add cross product of edges to the total normal
			p0 = vertices[faceinds[f][v]];
			p1 = vertices[faceinds[f][(v+1)%numvertsperface[v]]];
			p2 = vertices[faceinds[f][(v+2)%numvertsperface[v]]];
			faces[f].n.x += (p1.y - p0.y)*(p2.z - p0.z) - (p1.z - p0.z)*(p2.y - p0.y);
			faces[f].n.y += (p1.z - p0.z)*(p2.x - p0.x) - (p1.x - p0.x)*(p2.z - p0.z);
			faces[f].n.z += (p1.x - p0.x)*(p2.y - p0.y) - (p1.y - p0.y)*(p2.x - p0.x);

			// add the vertex position to the centroid
			centroid.x += p0.x;
			centroid.y += p0.y;
			centroid.z += p0.z;
		}

		// normalize the normals and set the signed distance to origin
		centroid.x /= numvertsperface[f];
		centroid.y /= numvertsperface[f];
		centroid.z /= numvertsperface[f];
		norm(faces[f].n);
		faces[f].d = -dot(faces[f].n, centroid);
	}
}

r3d_real r3d_orient(r3d_rvec3* verts) {
	r3d_real adx, bdx, cdx;
	r3d_real ady, bdy, cdy;
	r3d_real adz, bdz, cdz;
	adx = verts[0].x - verts[3].x;
	bdx = verts[1].x - verts[3].x;
	cdx = verts[2].x - verts[3].x;
	ady = verts[0].y - verts[3].y;
	bdy = verts[1].y - verts[3].y;
	cdy = verts[2].y - verts[3].y;
	adz = verts[0].z - verts[3].z;
	bdz = verts[1].z - verts[3].z;
	cdz = verts[2].z - verts[3].z;
	return -ONE_SIXTH*(adx * (bdy * cdz - bdz * cdy)
			+ bdx * (cdy * adz - cdz * ady)
			+ cdx * (ady * bdz - adz * bdy));
}

void r3d_print(r3d_poly* poly) {
	r3d_int v;
	for(v = 0; v < poly->nverts; ++v) {
		printf("  vertex %d: pos = ( %.10e , %.10e , %.10e ), nbrs = %d %d %d\n", 
				v, poly->verts[v].pos.x, poly->verts[v].pos.y, poly->verts[v].pos.z, poly->verts[v].pnbrs[0], poly->verts[v].pnbrs[1], poly->verts[v].pnbrs[2]);
	}
}
