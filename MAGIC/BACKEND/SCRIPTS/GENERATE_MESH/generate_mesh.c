/* 
 * MUPHY mesh generator from STL file.
 * Written by Mauro Bisson (mauro.bis@gmail.com)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <float.h> // for [FLT|DBL]_M[IN|AX]

#define DIM		3
#define MAX_FNAME	256
#define	MAX_LINE	256

#define MAX_INOUTLET	100

#define MIN(a,b)        (((a)<(b))?(a):(b))
#define MAX(a,b)        (((a)>(b))?(a):(b))

#define TOS(x)		#x

#if defined(READ_DOUBLE)
#define REAL		double
#define FPSUFF(a)       a
#define REAL_SPEC	TOS(%lf)
#else
#define REAL		float
#define FPSUFF(a)       a##f
#define REAL_SPEC 	TOS(%f)
#endif

#define CEIL_CPU(x)     	FPSUFF(ceil)((x))
#define FLOOR_CPU(x)    	FPSUFF(floor)((x))
#define SQRT_CPU(x)     	FPSUFF(sqrt)((x))
#define FABS_CPU(x)		FPSUFF(fabs)((x))
#define ROUND_CPU(x)    	FPSUFF(round)((x))
#define ACOS_CPU(x)		FPSUFF(acos)((x))
#define COPYSIGN_CPU(x,y)	FPSUFF(copysign)((x),(y))

#define EPSILON		FPSUFF(1.0E-6)
#define HALF		FPSUFF(0.5E+0)
#define ZERO		FPSUFF(0.0E+0)
#define ONE             FPSUFF(1.0E+0)
#define TWO             FPSUFF(2.0E+0)
#define THREE		FPSUFF(3.0E+0)

#if defined(READ_DOUBLE)
#define	REAL_MAX	DBL_MAX
#define	REAL_MIN	(-DBL_MAX)
#else
#define	REAL_MAX	FLT_MAX
#define	REAL_MIN	(-FLT_MAX)
#endif

#ifndef M_PI
#define M_PI 		3.14159265358979323846
#endif

#define ALLOC_SSIZE     10000

#define	INTERNAL	1
#define	EXTERNAL	2
#define	INLET		3
#define	OUTLET		4
#define DEAD		5

#define INTERNAL_BITOFF	0
#define EXTERNAL_BITOFF	4
#define INLET_BITOFF	8
#define OUTLET_BITOFF	12
#define DEAD_BITOFF	16

/* in futuro usare questi al posto di INTERNAL, EXTERNAL, ecc... 
 * togliendo la "A" iniziale ovv.te. */
#define	AINTERNAL(subtype)	((subtype+1)<<INTERNAL_BITOFF)
#define	AEXTERNAL(subtype)	((subtype+1)<<EXTERNAL_BITOFF)
#define	AINLET(subtype)		((subtype+1)<<INLET_BITOFF)
#define	AOUTLET(subtype)	((subtype+1)<<OUTLET_BITOFF)
#define	ADEAD(subtype)		((subtype+1)<<DEAD_BITOFF)

#define IS_INTERNAL(type)	(((type)&(0xF<<INTERNAL_BITOFF))!=0)
#define IS_EXTERNAL(type)	(((type)&(0xF<<EXTERNAL_BITOFF))!=0)
#define IS_INLET(type)		(((type)&(0xF<<INLET_BITOFF))!=0)
#define IS_OUTLET(type)		(((type)&(0xF<<OUTLET_BITOFF))!=0)
#define IS_DEAD(type)		(((type)&(0xF<<DEAD_BITOFF))!=0)

#define	ORIGIN		10

#define REALLOC(ptr,type,to_size,ret)	do {\
						type *__tmp_ptr;\
						__tmp_ptr = (type *)realloc((ptr), (to_size));\
						if (__tmp_ptr == NULL) {\
							ret = 0;\
						}\
						else {\
							ptr = __tmp_ptr;\
							ret = 1;\
						}\
					} while(0)

#define LOG(lvl,format,...)	{\
				if (lvl <= verblvl) {\
					fprintf(stderr, format, ## __VA_ARGS__);\
					fflush(stderr);\
				}\
				}

static int verblvl = 0;

typedef struct {
	int i, j, k;
	long long i4;
	char type;
} point_t;

typedef struct {
	point_t		*points;
	long long	pnum;
	long long	pnum_alloc;

	int		imin, imax;
	int		jmin, jmax;
	int		kmin, kmax;

	int		idim, jdim, kdim;
	int		i_tr, j_tr, k_tr;
} mesh_t;

typedef struct {
	REAL v1[DIM];
	REAL v2[DIM];
	REAL v3[DIM];
	REAL n[DIM];
} triangle_t;

typedef struct {
	triangle_t	*triangles;
	long long	tnum;
	long long       tnum_alloc;

	REAL		xmin, xmax;
	REAL		ymin, ymax;
	REAL		zmin, zmax;
} surface_t;

void usage(const char *pname) {

	fprintf(stderr, "Usage: %s -f <geometry_stlfile> [-s <scale_factor>] [-n <normals_file>] "
//			"[-x <geometry_unit>] [-y <physycal_viscosity>] [-z <lattice_viscosity>] "
			"[-i <inlet_stlfile> [-i ...]] [-o <outlet_stlfile> [-o ...]] [-v]\n", pname);
	return;
}

static long long binsearch(point_t *c, long long i4, long long ffree) {

	long long	min = 0;
	long long	max = ffree-1;
	long long	mid = (min + max)>>1;

	while(min <= max) {

		if (c[mid].i4 == i4)	return mid;
		if (c[mid].i4  < i4)	min = mid+1;
		else			max = mid-1;
		mid = (max + min) >> 1;
	}
	return -1;
}
/*
static REAL proddot(REAL *v1, REAL *v2) {

	return SQRT_CPU(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

static void normalizev(REAL *v) {

	REAL length;

	length = proddot(v,v);
	v[0] /= length;
	v[1] /= length;
	v[2] /= length;

	return;
}
*/
static void prodvect(REAL *a, REAL *b, REAL *c) {

	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];

	return;
}

static REAL dvol(REAL *v1, REAL *v2, REAL *v3, REAL *n) {

	REAL	bz;
	REAL	v1proj[3], v2proj[3];
	REAL	tri_proj_area;
	REAL	vprod[3];

	/* z coo of baricentrum */
	bz = (v1[2] + v2[2] + v3[2])/THREE;

	/* v2-v1 projection on plane z=0 */
	v1proj[0] = v2[0] - v1[0];
	v1proj[1] = v2[1] - v1[1];
	v1proj[2] = ZERO;
	
	/* v3-v1 projection on plane z=0 */
	v2proj[0] = v3[0] - v1[0];
	v2proj[1] = v3[1] - v1[1];
	v2proj[2] = ZERO;

	prodvect(v1proj, v2proj, vprod);

	tri_proj_area = FABS_CPU(HALF*(vprod[2])) * bz;

	return COPYSIGN_CPU(tri_proj_area, n[2]);
}

static REAL dsur(REAL *v1, REAL *v2, REAL *v3) {

	REAL    vprod[3];
	REAL	va[3], vb[3];

	va[0] = v2[0] - v1[0];
	va[1] = v2[1] - v1[1];
	va[2] = v2[2] - v1[2];

	vb[0] = v3[0] - v1[0];
	vb[1] = v3[1] - v1[1];
	vb[2] = v3[2] - v1[2];

	prodvect(va, vb, vprod);

	return SQRT_CPU(vprod[0]*vprod[0]+
			vprod[1]*vprod[1]+
			vprod[2]*vprod[2]) / TWO;
}

static char *skip_spaces(char *str) {
	
	for(; (*str == ' ') || (*str == '\t'); str++);
	return str;
}

static char *get_neline(char *buf, int n, FILE *f) {

	char	*rv;

	while(1) {
		rv = fgets(buf, n, f);
		if (rv == NULL) break;
		rv = skip_spaces(buf);
		if ((rv[0] != '\n') && (rv[0] != '\0')) break;
	}
	return rv;
}

/* 
 * Loads the triangle set contained into file fname, scaled by the
 * coefficient scale. 
 * Returns the pointer to the surface_t structure created on success,
 * NULL otherwise.
 */
static surface_t *load_surface(const char *fname, REAL scale) {

	FILE 		*f;
	surface_t	*s;
	char 		buffer[MAX_LINE], *line;
	long long	c;
	REAL		volume = ZERO, area = ZERO;

	if ((f = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open file %s for reading, quitting...\n", __func__, fname);
		return NULL;
	}

	s = (surface_t *)malloc(sizeof(surface_t));
	if (s == NULL) {

		fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n", 
				__func__, sizeof(surface_t));
		fclose(f);
		return NULL;
	}
	s->triangles = (triangle_t *)malloc(sizeof(triangle_t)*ALLOC_SSIZE);
	if (s == NULL) {

		fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n", 
				__func__, sizeof(triangle_t)*ALLOC_SSIZE);
		free(s);
		fclose(f);
		return NULL;
	}
	s->tnum_alloc = ALLOC_SSIZE;

	s->xmin = s->ymin = s->zmin = REAL_MAX;
	s->xmax = s->ymax = s->zmax = REAL_MIN;

//	LOG(0, "\nReading %s...\n", fname);
	c = 0;
	get_neline(buffer, MAX_LINE, f); /* discards "solid ..." line */
	while( (line = get_neline(buffer, MAX_LINE, f)) ) {

		if (c >= s->tnum_alloc) {

			triangle_t	*tmp;
			tmp = (triangle_t *)realloc(s->triangles, (s->tnum_alloc + ALLOC_SSIZE)*sizeof(triangle_t));
			if (tmp == NULL) {
				fprintf(stderr, "%s: cannot realloc to %zu bytes, quitting...\n", 
						__func__, (size_t)(s->tnum_alloc + ALLOC_SSIZE)*sizeof(triangle_t));
				free(s->triangles);
				free(s);
				fclose(f);
				return NULL;
			}
			s->triangles = tmp;
			s->tnum_alloc += ALLOC_SSIZE;
		}

		if (strncmp(line, "endsolid", 8) == 0) break;

		sscanf(line, "facet normal "REAL_SPEC" "REAL_SPEC" "REAL_SPEC"\n",
			&(s->triangles[c].n[0]), &(s->triangles[c].n[1]), &(s->triangles[c].n[2]));

		get_neline(buffer, MAX_LINE, f); /* reads "outer loop ..." */

		line = get_neline(buffer, MAX_LINE, f); /* reads "vertex ..." */
		sscanf(line, "vertex "REAL_SPEC" "REAL_SPEC" "REAL_SPEC"\n", 
			&(s->triangles[c].v1[0]), &(s->triangles[c].v1[1]), &(s->triangles[c].v1[2]));
		s->triangles[c].v1[0] *= scale;
		s->triangles[c].v1[1] *= scale;
		s->triangles[c].v1[2] *= scale;

		line = get_neline(buffer, MAX_LINE, f); /* reads "vertex ..." */
		sscanf(line, "vertex "REAL_SPEC" "REAL_SPEC" "REAL_SPEC"\n",
			&(s->triangles[c].v2[0]), &(s->triangles[c].v2[1]), &(s->triangles[c].v2[2]));
		s->triangles[c].v2[0] *= scale;
		s->triangles[c].v2[1] *= scale;
		s->triangles[c].v2[2] *= scale;

		line = get_neline(buffer, MAX_LINE, f); /* reads "vertex ..." */
		sscanf(line, "vertex "REAL_SPEC" "REAL_SPEC" "REAL_SPEC"\n",
			&(s->triangles[c].v3[0]), &(s->triangles[c].v3[1]), &(s->triangles[c].v3[2]));
		s->triangles[c].v3[0] *= scale;
		s->triangles[c].v3[1] *= scale;
		s->triangles[c].v3[2] *= scale;

		s->xmin = MIN(s->xmin, MIN(s->triangles[c].v1[0], MIN(s->triangles[c].v2[0], s->triangles[c].v3[0])));
		s->ymin = MIN(s->ymin, MIN(s->triangles[c].v1[1], MIN(s->triangles[c].v2[1], s->triangles[c].v3[1])));
		s->zmin = MIN(s->zmin, MIN(s->triangles[c].v1[2], MIN(s->triangles[c].v2[2], s->triangles[c].v3[2])));
		s->xmax = MAX(s->xmax, MAX(s->triangles[c].v1[0], MAX(s->triangles[c].v2[0], s->triangles[c].v3[0])));
		s->ymax = MAX(s->ymax, MAX(s->triangles[c].v1[1], MAX(s->triangles[c].v2[1], s->triangles[c].v3[1])));
		s->zmax = MAX(s->zmax, MAX(s->triangles[c].v1[2], MAX(s->triangles[c].v2[2], s->triangles[c].v3[2])));

		get_neline(buffer, MAX_LINE, f);
		get_neline(buffer, MAX_LINE, f);
		
		volume += dvol(s->triangles[c].v1, s->triangles[c].v2, s->triangles[c].v3, s->triangles[c].n);
		area   += dsur(s->triangles[c].v1, s->triangles[c].v2, s->triangles[c].v3);
		
		c++;
	}
	s->tnum = c;
	LOG(0, "\tread %lld triangles\n", c);
	LOG(0, "\tarea: %f square units (non-scaled %f)\n", area, area/(scale*scale));
	LOG(0, "\tvolume: %f cubic units (non-scaled %f)\n", volume, volume/(scale*scale*scale));
	LOG(1, "\tbounding box: [%E, %E] [%E, %E] [%E, %E]\n",
                s->xmin, s->xmax, s->ymin, s->ymax, s->zmin, s->zmax);
	fclose(f);
/*
for(s->tnum=0; s->tnum<c; s->tnum++) {

	fprintf(stdout, "vertex %E %E %E\n", s->triangles[s->tnum].v1[0], s->triangles[s->tnum].v1[1], s->triangles[s->tnum].v1[2]);
	fprintf(stdout, "vertex %E %E %E\n", s->triangles[s->tnum].v2[0], s->triangles[s->tnum].v2[1], s->triangles[s->tnum].v2[2]);
	fprintf(stdout, "vertex %E %E %E\n", s->triangles[s->tnum].v3[0], s->triangles[s->tnum].v3[1], s->triangles[s->tnum].v3[2]);
}
*/
	return s;
}

/* point_t compare function used for qsort */
static int points_cmp(const void *a, const void *b) {

	point_t *p1 = (point_t *)a;
	point_t *p2 = (point_t *)b;

	if (p1->i4 < p2->i4) return -1;
	if (p1->i4 > p2->i4) return 1;
	return 0;
}

/* 
 * Return in norm[] the normal to the triangle identifed by
 * vectors v1, v2, v3.
 * The directon of the normail is such that v1, v2 and v3 are
 * encountered counter-clockwise.
 */
/*
static void get_normal(REAL *v1, REAL *v2, REAL *v3, REAL *norm) {

	REAL a[3];
	REAL b[3];
	REAL mod;

	a[0] = v1[0] - v2[0];
	a[1] = v1[1] - v2[1];
	a[2] = v1[2] - v2[2];

	b[0] = v2[0] - v3[0];
	b[1] = v2[1] - v3[1];
	b[2] = v2[2] - v3[2];

	norm[0] = a[1] * b[2] - a[2] * b[1];
	norm[1] = a[2] * b[0] - a[0] * b[2];
	norm[2] = a[0] * b[1] - a[1] * b[0];

	mod = SQRT_CPU(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0] /= mod;
	norm[1] /= mod;
	norm[2] /= mod;

	return;
}
*/

/* 
 * Returns 1 in point p is belongs to triangle v1, v2 and v3.
 * 0 otherwise.
 */
static int is_in_triangle(REAL px, REAL py, REAL pz, 
			  REAL *v1, REAL *v2, REAL *v3) {

	REAL p1[3], p2[3], p3[3];
	REAL p1m, p2m, p3m;
	REAL alpha, beta, gamma;

	p1[0] = v1[0] - px;
	p1[1] = v1[1] - py;
	p1[2] = v1[2] - pz;
	p1m = SQRT_CPU(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);

	p2[0] = v2[0] - px;
	p2[1] = v2[1] - py;
	p2[2] = v2[2] - pz;
	p2m = SQRT_CPU(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]);

	p3[0] = v3[0] - px;
	p3[1] = v3[1] - py;
	p3[2] = v3[2] - pz;
	p3m = SQRT_CPU(p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2]);

	alpha = ACOS_CPU((p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]) / (p1m*p2m));
	beta  = ACOS_CPU((p1[0]*p3[0] + p1[1]*p3[1] + p1[2]*p3[2]) / (p1m*p3m));
	gamma = ACOS_CPU((p3[0]*p2[0] + p3[1]*p2[1] + p3[2]*p2[2]) / (p3m*p2m));

//	if ((alpha+beta+gamma) >= (2*M_PI*FPSUFF(0.9E+0))) return 1;
	if ((alpha+beta+gamma) >= (2*M_PI*(ONE-EPSILON))) return 1;
	return 0;
}

static REAL distpp(REAL *p, REAL *q) {

	return SQRT_CPU((p[0]+q[0])*(p[0]+q[0]) + 
			(p[1]+q[1])*(p[1]+q[1]) + 
			(p[2]+q[2])*(p[2]+q[2]));
}

static REAL distpl(REAL *p, REAL *v1, REAL *v2) {

	REAL p1[3], p2[3], n[3];
	REAL n_mod;

	p1[0] = p[0] - v1[0];
	p1[1] = p[1] - v1[1];
	p1[2] = p[2] - v1[2];

	p2[0] = v2[0] - v1[0];
	p2[1] = v2[1] - v1[1];
	p2[2] = v2[2] - v1[2];
/*
	n[0] = p1[1] * p2[2] - p1[2] * p2[1];
        n[1] = p1[2] * p2[0] - p1[0] * p2[2];
        n[2] = p1[0] * p2[1] - p1[1] * p2[0];
*/
	prodvect(p1, p2, n);

	n_mod = SQRT_CPU(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

	return n_mod/distpp(v1, v2);
}

/* 
 * Returns the signed distance between point p andtriangle v1, v2 and v3.
 * If p lies in the semi-space pointed by n the distance is positive,
 * otherwise is negative.
 */
static REAL dist_point_triangle(REAL *p, REAL *v1, REAL *v2, REAL *v3, REAL *n) {

	REAL	pt[3];
	REAL	pv1[3];
	REAL	dotp;
	REAL	ds1, ds2, ds3, dp1, dp2, dp3;
	REAL	ret;

	pv1[0] = p[0] - v1[0];
	pv1[1] = p[1] - v1[1];
	pv1[2] = p[2] - v1[2];

	dotp = pv1[0]*n[0] + pv1[1]*n[1] + pv1[2]*n[2];
	pt[0] = p[0] - dotp*n[0];
	pt[1] = p[1] - dotp*n[1];
	pt[2] = p[2] - dotp*n[2];

	if (is_in_triangle(pt[0], pt[1], pt[2], v1, v2, v3)) 
		return dotp;

	dp1 = distpp(p, v1);
	dp2 = distpp(p, v2);
	dp3 = distpp(p, v3);
	ret = MIN(dp1,MIN(dp2, dp3));

	ds1 = distpl(p, v1, v2);
	ds2 = distpl(p, v2, v3);
	ds3 = distpl(p, v3, v1);
	ret = MIN(ret,MIN(ds1,MIN(ds2, ds3)));

	return COPYSIGN_CPU(ret, dotp);
}

/* 
 * Returns a mest_t pointer containing the mesh approximating
 * the triangle surface s (containing snum elements).
 * The mesh contains all the points at distance less than spc
 * form any of the triangles in s and each point. Moreover
 * the points are tagged as either INTERNAL or EXTERNAL.
 * Please note that the WHOLE internal volume is enclosed
 * into INTERNAL point and EXTERNAL ones are only in contact with
 * the region external to the volume.
 * The mesh returned is translated so that no point has 
 * coordinates less than (1,1,1).
 * The points stored in the returned mesh are sorted for increasing
 * i4.
 */
static mesh_t *generate_mesh_surf(surface_t *s, REAL spc) {

	mesh_t		*m;
	long long	n, c, l;
	int		imin, imax, jmin, jmax, kmin, kmax;
	REAL		xmin, xmax, ymin, ymax, zmin, zmax;
	int		i, j, k;

	//REAL	norm[3];
	
	m = (mesh_t *)malloc(sizeof(mesh_t));
	if (m == NULL) {
		fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n", 
				__func__, sizeof(mesh_t));
		return NULL;
	}

	m->points = (point_t *)malloc(sizeof(point_t)*ALLOC_SSIZE);
	if (m->points == NULL) {
		fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n", 
				__func__, sizeof(point_t)*ALLOC_SSIZE);
		free(m);
		return NULL;
	}
	m->pnum_alloc = ALLOC_SSIZE;

//	LOG(0, "Generating surface mesh...\n");
	c = 0;
	for(n = 0; n < s->tnum; n++) {

		xmin = MIN(s->triangles[n].v1[0], MIN(s->triangles[n].v2[0], s->triangles[n].v3[0]));
		ymin = MIN(s->triangles[n].v1[1], MIN(s->triangles[n].v2[1], s->triangles[n].v3[1]));
		zmin = MIN(s->triangles[n].v1[2], MIN(s->triangles[n].v2[2], s->triangles[n].v3[2]));

		xmax = MAX(s->triangles[n].v1[0], MAX(s->triangles[n].v2[0], s->triangles[n].v3[0]));
		ymax = MAX(s->triangles[n].v1[1], MAX(s->triangles[n].v2[1], s->triangles[n].v3[1]));
		zmax = MAX(s->triangles[n].v1[2], MAX(s->triangles[n].v2[2], s->triangles[n].v3[2]));

		imin = FLOOR_CPU(xmin-EPSILON);
		jmin = FLOOR_CPU(ymin-EPSILON);
		kmin = FLOOR_CPU(zmin-EPSILON);

		imax = CEIL_CPU(xmax+EPSILON);
		jmax = CEIL_CPU(ymax+EPSILON);
		kmax = CEIL_CPU(zmax+EPSILON);
/*
		fprintf(stdout, "Trinagolo: (%E, %E, %E), (%E, %E, %E), (%E, %E, %E) imM=[%d, %d] jmM=[%d, %d] kmM=[%d, %d]\n",
				s->triangles[n].v1[0], s->triangles[n].v1[1], s->triangles[n].v1[2], 
				s->triangles[n].v2[0], s->triangles[n].v2[1], s->triangles[n].v2[2], 
				s->triangles[n].v3[0], s->triangles[n].v3[1], s->triangles[n].v3[2],
				imin, imax, jmin, jmax, kmin, kmax);
*/

/*
		get_normal(s->triangles[n].v1, s->triangles[n].v2, s->triangles[n].v3, norm);
		if ( (FABS_CPU(norm[0] - s->triangles[n].n[0]) > FPSUFF(0.1E-5)) ||\
		     (FABS_CPU(norm[1] - s->triangles[n].n[1]) > FPSUFF(0.1E-5)) ||\
		     (FABS_CPU(norm[2] - s->triangles[n].n[2]) > FPSUFF(0.1E-5)) ) {

			fprintf(stderr, "Normal read from file (%E, %E, %E) differs from computed ((%E, %E, %E)... quitting.\n",
				s->triangles[n].n[0], s->triangles[n].n[1], s->triangles[n].n[2], norm[0], norm[1], norm[2]);
		}
*/
		for(i = imin; i <= imax; i++) {
		for(j = jmin; j <= jmax; j++) {
		for(k = kmin; k <= kmax; k++) {

			/* distance point-triangle */
			REAL dpt;
			REAL p[3];

			p[0] = (REAL)i;
			p[1] = (REAL)j;
			p[2] = (REAL)k;
/*
			fprintf(stdout, "T = {(%E, %E, %E), (%E, %E, %E), (%E, %E, %E)}\n", 
					s->triangles[n].v1[0], s->triangles[n].v1[1], s->triangles[n].v1[2], 
					s->triangles[n].v2[0], s->triangles[n].v2[1], s->triangles[n].v2[2], 
					s->triangles[n].v3[0], s->triangles[n].v3[1], s->triangles[n].v3[2]);
			fprintf(stdout, "d({%E, %E, %E}, T) = %E\n\n", px, py, pz, dpt);
*/
			if (c >= m->pnum_alloc) {

				point_t	*tmp;
				tmp = (point_t *)realloc(m->points, (m->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
				if (tmp == NULL) {
					fprintf(stderr, "%s: cannot realloc to %zu bytes, quitting...\n", 
							__func__, (size_t)(m->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
					free(m->points);
					free(m);
					return NULL;
				}

				m->points = tmp;
				m->pnum_alloc += ALLOC_SSIZE;
				//fprintf(stderr, "ALLOCATO: %lld bytes\n", m->pnum_alloc);
			}
			
			dpt = dist_point_triangle(p, s->triangles[n].v1, s->triangles[n].v2, s->triangles[n].v3, s->triangles[n].n);
			if (FABS_CPU(dpt) < spc) {
		
					m->points[c].i = i; 
					m->points[c].j = j; 
					m->points[c].k = k;
					
					/* type internal(0) - external(1) */
					if (dpt > ZERO) m->points[c].type = EXTERNAL;
					else		m->points[c].type = INTERNAL;
					
					c++;
	
				}
		}}}
	}

	/* compute mesh bbox (this can be moved into if (dpt < spc) {} block */
	imin = imax = m->points[0].i;
	jmin = jmax = m->points[0].j;
	kmin = kmax = m->points[0].k;
	for(n = 1; n < c; n++) {
		imin = MIN(imin, m->points[n].i);
		imax = MAX(imax, m->points[n].i);
		jmin = MIN(jmin, m->points[n].j);
		jmax = MAX(jmax, m->points[n].j);
		kmin = MIN(kmin, m->points[n].k);
		kmax = MAX(kmax, m->points[n].k);
	}
	LOG(1, "\toriginal surface mesh bounding box: [%d, %d] [%d, %d] [%d, %d]\n",
		imin, imax, jmin, jmax, kmin, kmax);

	/* translating mesh to (ORIGIN,ORIGIN,ORIGIN) and computing i4s */
/*	m->imin = ORIGIN;
	m->jmin = ORIGIN;
	m->kmin = ORIGIN;
*/
	imax -= (imin-ORIGIN);
	jmax -= (jmin-ORIGIN);
	kmax -= (kmin-ORIGIN);
	imax += ORIGIN;
	jmax += ORIGIN;
	kmax += ORIGIN;

	/* save tringle frame to mesh frame translation */
	m->i_tr = ORIGIN-imin;
	m->j_tr = ORIGIN-jmin;
	m->k_tr = ORIGIN-kmin;

	LOG(0, "\ttranslation vector: %d,%d,%d\n", m->i_tr, m->j_tr, m->k_tr);

	for(n = 0; n < c; n++) {
		/* i = i - imin + 1 */
		m->points[n].i -= (imin-ORIGIN);
		m->points[n].j -= (jmin-ORIGIN);
		m->points[n].k -= (kmin-ORIGIN);
		m->points[n].i4 = (long long)m->points[n].k*(jmax+1)*(imax+1) + 
				  (long long)m->points[n].j*(imax+1) + 
			  	  (long long)m->points[n].i;
	}

	m->imin = ORIGIN;
	m->jmin = ORIGIN;
	m->kmin = ORIGIN;

	m->imax = imax;
        m->jmax = jmax;
        m->kmax = kmax;

	m->idim = m->imax - m->imin +1;
        m->jdim = m->jmax - m->jmin +1;
        m->kdim = m->kmax - m->kmin +1;

	LOG(1, "\tsorting surface mesh points...\n");
	/* sort points with respect to i4s */
	qsort(m->points, c, sizeof(point_t), points_cmp);

	/* ...remove duplicates */
	l = 1;
	for(n = 1; n < c; n++)
		if (m->points[n].i4 != m->points[n-1].i4)
			m->points[l++] = m->points[n]; /* structure copy, np here */
		else {
			/* when collapsing duplicates in signle points
			 * make them EXTERNAL if there is at least an
			 * EXTERNAL in the duplicates 
			 */
			if (m->points[n].type == INTERNAL)
				m->points[l-1].type = INTERNAL;
/*			fprintf(stderr, "m->points[%d]={%d, %d, %d}[%lld] == m->points[%d]={%d, %d, %d}[%lld]\n",
					n, m->points[n].i, m->points[n].j, m->points[n].k, m->points[n].i4, 
					n-1, m->points[n-1].i, m->points[n-1].j, m->points[n-1].k, m->points[n-1].i4);
*/		}
	m->pnum = l;
	LOG(1, "\tfinal surface mesh bounding box: [%d, %d] [%d, %d] [%d, %d]\n",
			m->imin, m->imax, m->jmin, m->jmax, m->kmin, m->kmax);
	LOG(1, "\tsurface mesh number of points: %lld\n", m->pnum);

	return m;
}

/* 
 * Fills mesh m with the points enclosed by the INTERNAL points.
 * Returns 1 on success, 0 otherwise.
 */
static int fill_mesh(mesh_t *m) {

	int 		i;
       	long long	n, c;
	point_t		prev, curr;

//	LOG(1, "Filling mesh...\n");

	c = m->pnum;
	for(n = 1; n < m->pnum; n++) {

		prev = m->points[n-1];
		curr = m->points[n];

		if ( !(prev.i < (curr.i-1)) ) continue;
		if ( !((prev.type == INTERNAL) && (curr.type == INTERNAL)) ) continue;
		if ( !((prev.j == curr.j) && (prev.k == curr.k)) ) {
			fprintf(stderr, "\tFollowing an internal mesh point possibly inside dead region,\n"
					"\tplease verify the mesh produced.\n");
			/*
					fprintf(stderr, "%s: mesh point (%d %d %d){%lld} tagged as INTERNAL exposed"
						        "to the outside, quitting...\n", 
							__func__, prev.i, prev.j, prev.k, prev.i4);
					return 0;
			*/
			continue;
		}

		/* insert internal points */
		for (i = prev.i+1; i < curr.i; i++) {
			if (c >= m->pnum_alloc) {

				point_t	*tmp;
				tmp = (point_t *)realloc(m->points, (m->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
				if (tmp == NULL) {
					fprintf(stderr, "%s: cannot realloc to %zu bytes, quitting...\n", 
							__func__, (size_t)(m->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
					return 0;
				}
				m->points = tmp;
				m->pnum_alloc += ALLOC_SSIZE;
			}
			m->points[c].i = i;
			m->points[c].j = prev.j;
			m->points[c].k = prev.k;
			m->points[c].type = INTERNAL;
			m->points[c].i4 = prev.i4 + (long long)(i - prev.i);
					/*(long long)k*(m->jmax+1)*(m->imax+1) +
					  (long long)j*(m->imax+1) +
					  (long long)(i);
					  */
			c++;
		}
	}

	LOG(1, "\taddedd %lld INTERNAL mesh points\n", (c - m->pnum));

	/* sort all the points */
	m->pnum = c;
	qsort(m->points, m->pnum, sizeof(point_t), points_cmp);

	/* ...remove duplicates (just as a sanity check) */
	c = 1;
	for(n = 1; n < m->pnum; n++)
		if (m->points[n].i4 != m->points[n-1].i4)
			m->points[c++] = m->points[n]; /* structure copy, np here */

	if (c < m->pnum) {
		LOG(1, "\tselected %lld duplicate points\n", m->pnum - c);
		m->pnum = c;
	}
	else {
		LOG(1, "\tselected no duplicate points\n");
	}

	return 1;
}

/* PROVARE A CICLARE SOLO SUI VICINI A DIST 1 */
/* Smooths mesh m and return the numeber of removed points. */
/* Looping through the 2-neighborhood and discriminating above
 * 125/3 neighbors smooths better meshes of increasing size */
static long long smooth_mesh(mesh_t *m) {

	long long	n, c;
	int		i, j, k;
	long long	i4;

	/* mark as "DEAD" pointes with le 9 neighbors */
	for(n = 0; n < m->pnum; n++) {
		int neighs = 0, delete = 1;
		long long ifl;
		
		for(i = m->points[n].i-1; i <= m->points[n].i+1; i++) {
		for(j = m->points[n].j-1; j <= m->points[n].j+1; j++) {
		for(k = m->points[n].k-1; k <= m->points[n].k+1; k++) {

			if ((i==m->points[n].i)&&(j==m->points[n].j)&&(k==m->points[n].k))
				continue;

			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);

			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type != DEAD) {
					neighs++;
				}
			}
			if (neighs > 9) { /* we'll keep this */
				delete = 0;
				goto out;
			}
		}}}
out:
		if (delete) {
			m->points[n].type = DEAD;
		}
	}
	
	/* remove DEAD-d points */
	c = 0;
	for(n = 0; n < m->pnum; n++)
		if (m->points[n].type != DEAD)
			m->points[c++] = m->points[n];

	n = m->pnum - c;
	fprintf(stderr, "\tSmoothing removed %lld INTERNAL points.\n", m->pnum - c);
	m->pnum = c;

	/* no need to re-sort (we are removing elements from a sorted list) */

	return n;
}

static long long smooth_mesh_aligned(mesh_t *m) {

	long long	n, c;
	int		i, j, k;
	long long	i4;

	/* mark as "DEAD" pointes with le 9 neighbors */
	for(n = 0; n < m->pnum; n++) {
		int neighs_26 = 0, neigh_6 = 0, delete = 1;
		long long ifl;
		
		for(i = -1; i <= 1; i++) {
		for(j = -1; j <= 1; j++) {
		for(k = -1; k <= 1; k++) {

			if ((i==0)&&(j==0)&&(k==0)) continue;

			i4 = ((long long)(m->points[n].k+k))*(m->jmax+1)*(m->imax+1) +
                             ((long long)(m->points[n].j+j))*(m->imax+1) +
                             ((long long)(m->points[n].i+i));

			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type != DEAD) {
					neighs_26++;
					if ((  i  && (!j) && (!k)) ||\
					    ((!i) &&   j  && (!k)) ||\
					    ((!i) && (!j) &&   k))
						neigh_6++;
				}
			}
			if ((neighs_26 > 9) && (neigh_6 >2)) { /* we'll keep this */
				delete = 0;
				goto out;
			}
		}}}
out:
		if (delete) {
			m->points[n].type = DEAD;
		}
	}
	
	/* remove DEAD-d points */
	c = 0;
	for(n = 0; n < m->pnum; n++)
		if (m->points[n].type != DEAD)
			m->points[c++] = m->points[n];

	n = m->pnum - c;
	fprintf(stderr, "\tSmoothing removed %lld INTERNAL points.\n", m->pnum - c);
	m->pnum = c;

	/* no need to re-sort (we are removing elements from a sorted list) */

	return n;
}

static long long smooth_mesh_aligned_oneshot(mesh_t *m) {

	long long	n, c, l;
	int		i, j, k;
	long long	i4;

	long long	*nhood; 
	long long	nhood_end, nhood_alloc;		

	nhood = (long long *)malloc(sizeof(long long)*ALLOC_SSIZE);
	if (nhood == NULL) {
		fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n", 
				__func__, (size_t)(sizeof(long long)*ALLOC_SSIZE));
		return -1;
	}
	nhood_alloc = ALLOC_SSIZE;

	/* mark as "DEAD" pointes with le 9 neighbors */
	for(n = 0; n < m->pnum; n++) {

		if (m->points[n].type == DEAD) continue;

		nhood[0] = n;
		nhood_end = 1;
		for(l = 0; l < nhood_end; l++) {
		
			int neighs_26 = 0, neigh_6 = 0, delete = 1;
		       	long long ifl = nhood[l], nifl;

			if (m->points[ifl].type == DEAD) continue;

			for(i = -1; i <= 1; i++) {
			for(j = -1; j <= 1; j++) {
			for(k = -1; k <= 1; k++) {
	
				if ((i==0)&&(j==0)&&(k==0)) continue;
	
				i4 = ((long long)(m->points[ifl].k+k))*(m->jmax+1)*(m->imax+1) +
	                             ((long long)(m->points[ifl].j+j))*(m->imax+1) +
	                             ((long long)(m->points[ifl].i+i));
	
				nifl = binsearch(m->points, i4, m->pnum);
				if (nifl < 0) continue;
				if (m->points[nifl].type == DEAD) continue;

				if ((nhood_end+neighs_26) >= nhood_alloc) {
	
					long long	*tmp;
					tmp = (long long *)realloc(nhood, (nhood_alloc + ALLOC_SSIZE)*sizeof(long long));
					if (tmp == NULL) {
						fprintf(stderr, "%s: cannot realloc to %zu bytes, quitting...\n", 
								__func__, (size_t)(nhood_alloc + ALLOC_SSIZE)*sizeof(long long));
						free(nhood);
						return -1;
					}
					nhood = tmp;
					nhood_alloc += ALLOC_SSIZE;
				}

				nhood[nhood_end+neighs_26] = nifl;
				neighs_26++;

				if ((  i  && (!j) && (!k)) ||\
				    ((!i) &&   j  && (!k)) ||\
				    ((!i) && (!j) &&   k))
					neigh_6++;

				if ((neighs_26 > 9) && (neigh_6 >2)) delete = 0;
			}}}
			
			if (delete) {
				m->points[ifl].type = DEAD;
				nhood_end += neighs_26;
			}
		} /* for l */
	}
	
	/* remove DEAD-d points */
	c = 0;
	for(n = 0; n < m->pnum; n++)
		if (m->points[n].type != DEAD)
			m->points[c++] = m->points[n];

	n = m->pnum - c;
	LOG(1, "\tsmoothing removed %lld INTERNAL points\n", m->pnum - c);
	m->pnum = c;

	/* no need to re-sort (we are removing elements from a sorted list) */

	free(nhood);
	return n;
}

/* 
 * Inserts the wall by performing the following operations:
 * 1) remove EXTERNAL NODES;
 * 2) smooth the volume formed by the INTERNAL nodes;
 * 3) add a new layer of EXTERNAL nodes.
 * Returns 1 on succcess, 0 otherwise.
 */
static int insert_wall(mesh_t *m) {

	long long	n, c;
	int 		i, j, k;
	long long 	i4;

//	LOG(0, "Finalizing mesh...\n");

	/* remove points tagged EXTERNAL */
	c = 0;
	for(n = 0; n < m->pnum; n++)
		if (m->points[n].type == INTERNAL)
			m->points[c++] = m->points[n];
	/* no need to re-sort (we are removing elements from a sorted list) */

	LOG(1, "\tremoved %lld EXTERNAL points\n", m->pnum - c);
	m->pnum = c;

	/* smooths internal points volume */
/*	n = 0;
	while(smooth_mesh(m)) {
		n++;
		if (n > 5) break;
	}
*/

/*
	n = 0;
	while(smooth_mesh_aligned(m)) {
		n++;
//		if (n > 5) break;
	}
*/
	if (smooth_mesh_aligned_oneshot(m) < 0) 
		return 0;

	c = m->pnum;
	/* add external layer of points */
	for(n = 0; n < m->pnum; n++) {
		
		for(i = m->points[n].i-1; i <= m->points[n].i+1; i++) {
		for(j = m->points[n].j-1; j <= m->points[n].j+1; j++) {
		for(k = m->points[n].k-1; k <= m->points[n].k+1; k++) {

			if ((i==m->points[n].i)&&(j==m->points[n].j)&&(k==m->points[n].k))
				continue;

			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);

			if (binsearch(m->points, i4, m->pnum) < 0) {

				if (c >= m->pnum_alloc) {

					point_t	*tmp;
					tmp = (point_t *)realloc(m->points, (m->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
					if (tmp == NULL) {
						fprintf(stderr, "%s: cannot realloc to %zu bytes, quitting...\n", 
								__func__, (size_t)(m->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
						return 0;
					}
					m->points = tmp;
					m->pnum_alloc += ALLOC_SSIZE;
				}

				m->points[c].i = i;
				m->points[c].j = j;
				m->points[c].k = k;
				m->points[c].i4 = i4;
				m->points[c].type = EXTERNAL;
				c++;
//m->points[n].type = 6;
			}
		}}}
	}
	
	/* sort all the points */
	m->pnum = c;
	qsort(m->points, m->pnum, sizeof(point_t), points_cmp);
	/* ...remove duplicates */
	c = 1;
	for(n = 1; n < m->pnum; n++)
		if (m->points[n].i4 != m->points[n-1].i4)
			m->points[c++] = m->points[n];

	m->pnum = c;
	return 1;	
}

/* final touches to the mesh */
static void finalize_mesh(mesh_t *m) {

	long long n;

	fprintf(stderr, "Finalizing mesh...\n");

	/* compute mesh bbox (this can be moved into if (dpt < spc) {} block */
	m->imin = m->imax = m->points[0].i;
	m->jmin = m->jmax = m->points[0].j;
	m->kmin = m->kmax = m->points[0].k;
	for(n = 1; n < m->pnum; n++) {
		m->imin = MIN(m->imin, m->points[n].i);
		m->imax = MAX(m->imax, m->points[n].i);
		m->jmin = MIN(m->jmin, m->points[n].j);
		m->jmax = MAX(m->jmax, m->points[n].j);
		m->kmin = MIN(m->kmin, m->points[n].k);
		m->kmax = MAX(m->kmax, m->points[n].k);
	}

	for(n = 0; n < m->pnum; n++) {
		/* i = i - imin + 1 */
		m->points[n].i -= (m->imin-1);
		m->points[n].j -= (m->jmin-1);
		m->points[n].k -= (m->kmin-1);
		m->points[n].i4 = (long long)m->points[n].k*(m->jmax+1)*(m->imax+1) + 
				  (long long)m->points[n].j*(m->imax+1) + 
			  	  (long long)m->points[n].i;
	}

	m->idim = m->imax - m->imin + 1;
	m->jdim = m->jmax - m->jmin + 1;
	m->kdim = m->kmax - m->kmin + 1;

	fprintf(stderr, "\tFinal mesh bounding box: [%d, %d] [%d, %d] [%d, %d]\n",
			m->imin, m->imax, m->jmin, m->jmax, m->kmin, m->kmax);
	return;
}

static mesh_t *set_mesh_points(mesh_t *m, int color,
	  		       int mini, int maxi,
			       int minj, int maxj,
			       int mink, int maxk) {
	long long	n;
	long long	c;
	mesh_t		*inout;

	inout = (mesh_t *)malloc(sizeof(mesh_t));
        if (inout == NULL) {
                fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n",
                                __func__, sizeof(mesh_t));
                return NULL;
        }
	/* set unused mesh_t fields */
	inout->imin = inout->imax = 0;
	inout->jmin = inout->jmax = 0;
	inout->kmin = inout->kmax = 0;
	inout->i_tr = inout->j_tr = inout->k_tr = 0;
	inout->idim = inout->jdim = inout->kdim = 0;		

        inout->points = (point_t *)malloc(sizeof(point_t)*ALLOC_SSIZE);
        if (inout->points == NULL) {
                fprintf(stderr, "%s: cannot allocate %zu bytes, quitting...\n",
                                __func__, sizeof(point_t)*ALLOC_SSIZE);
                free(inout);
                return NULL;
        }
        inout->pnum_alloc = ALLOC_SSIZE;

	LOG(1, "\tsetting color %d for mesh nodes in [%d, %d] [%d, %d] [%d, %d]...",
		color, mini, maxi, minj, maxj, mink, maxk);
	c = 0;
	for(n = 0; n < m->pnum; n++) {

		if (m->points[n].type != INTERNAL) continue;

		if ((mini <= m->points[n].i) && (m->points[n].i <= maxi) &&\
		    (minj <= m->points[n].j) && (m->points[n].j <= maxj) &&\
		    (mink <= m->points[n].k) && (m->points[n].k <= maxk)) {

			/* set point color in global mesh */
			m->points[n].type = color;

			/* save point in inlet/outlet mesh */
			if (c >= inout->pnum_alloc) {

                                point_t *tmp;
                                tmp = (point_t *)realloc(inout->points, (inout->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
                                if (tmp == NULL) {
                                        fprintf(stderr, "%s: cannot realloc to %zu bytes, quitting...\n",
                                                        __func__, (size_t)(inout->pnum_alloc + ALLOC_SSIZE)*sizeof(point_t));
                                        free(inout->points);
                                        free(inout);
                                        return NULL;
                                }
                                inout->points = tmp;
                                inout->pnum_alloc += ALLOC_SSIZE;
                        }

			/* save mesh point specific structure */
			inout->points[c].i = m->points[n].i;
			inout->points[c].j = m->points[n].j;
			inout->points[c].k = m->points[n].k;
			inout->points[c].i4 = m->points[n].i4;

			c++;
		}
	}

	inout->pnum = c;
	LOG(1, "%lld points colored.\n", c);

	return inout;
}

static int check_closure(mesh_t *m) {

	long long	n;
	int		i, j, k;
	long long 	i4;
	int		close = 1;

	LOG(0, "\nChecking for mesh closure...");

	for(n = 0; n < m->pnum; n++) {

		if (m->points[n].type == EXTERNAL) continue;

		for(i = m->points[n].i-1; i <= m->points[n].i+1; i++) {
		for(j = m->points[n].j-1; j <= m->points[n].j+1; j++) {
		for(k = m->points[n].k-1; k <= m->points[n].k+1; k++) {
			
			if ((i==m->points[n].i)&&(j==m->points[n].j)&&(k==m->points[n].k)) continue;

			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);

			if (binsearch(m->points, i4, m->pnum) < 0) {
				/*
				fprintf(stderr, "%s: internal point (%d, %d, %d) is near to "
						"a dead point (%d, %d, %d).\n", __func__,
						m->points[n].i, m->points[n].j, m->points[n].k,
						i, j, k);
				*/
				close = 0;
				goto check_closure_out;
			}
		}}}
	}
check_closure_out:
	if (close == 1) {
		LOG(0, "ok\n");
	}
	else {
		LOG(0, "point (%d, %d, %d) touches a dead node (possibly more)\n",
			m->points[n].i, m->points[n].j, m->points[n].k);
	}
	return close;
}

static void print_mesh(const char *fname, mesh_t *m) {

	long long	n;
	FILE		*fp;

	if ((fp = fopen(fname, "w")) == NULL) {
                fprintf(stderr, "%s: cannot open file %s for writing, quitting...\n", __func__, fname);
                return;
        }

	fprintf(fp, "%d %d %d\n", m->imax, m->jmax, m->kmax);
	fprintf(fp, "5\n"); /* historical */ 
	for(n = 0; n < m->pnum; n++)
		fprintf(fp, "%d %d %d %d\n", 
			m->points[n].i, m->points[n].j, m->points[n].k, m->points[n].type);

	fclose(fp);
	return;
}

static void print_inlet_outlet(const char *fname, 
				  mesh_t *inlet_mesh[],  REAL inlet_norm[][3], int inlet_cnt,
				  mesh_t *outlet_mesh[], REAL outlet_norm[][3], int outlet_cnt) {

	int	i;
	FILE	*fp;

	long long	n;

	if ((fp = fopen(fname, "w")) == NULL) {
		fprintf(stderr, "%s: cannot open file %s for writing, quitting...\n", __func__, fname);
		return;
	}
	
	/* write header */
	fprintf(fp, "%d\n", inlet_cnt+outlet_cnt);
	for(i = 0; i < inlet_cnt; i++) 
		fprintf(fp, "%-4d  inlet pressure %3d %3d %3d  %12f %12f %12f  %12f\n", i+1, 
			inlet_mesh[i]->i_tr, inlet_mesh[i]->j_tr, inlet_mesh[i]->k_tr,
			-inlet_norm[i][0], -inlet_norm[i][1], -inlet_norm[i][2], ZERO);
			//(-1.5E-2)*(inlet_mesh[i]->i_tr + inlet_mesh[i]->j_tr + inlet_mesh[i]->k_tr));

	for(i = 0; i < outlet_cnt; i++) 
		fprintf(fp, "%-4d outlet     flow %3d %3d %3d  %12f %12f %12f  %12f\n", inlet_cnt+i+1, 
			outlet_mesh[i]->i_tr, outlet_mesh[i]->j_tr, outlet_mesh[i]->k_tr, 
			outlet_norm[i][0], outlet_norm[i][1], outlet_norm[i][2], FPSUFF(5.0E-3));
			//(-1.0E-2)*(outlet_mesh[i]->i_tr + outlet_mesh[i]->j_tr + outlet_mesh[i]->k_tr));

	/* write inlet points */
	n = 0;
	for(i = 0; i < inlet_cnt; i++)
		n += inlet_mesh[i]->pnum;
	fprintf(fp, "%10lld\n", n);
	for(i = 0; i < inlet_cnt; i++) 
		for(n = 0; n < inlet_mesh[i]->pnum; n++)
			fprintf(fp, "%d %d %d  %d\n",
				inlet_mesh[i]->points[n].i,
				inlet_mesh[i]->points[n].j,
				inlet_mesh[i]->points[n].k,
				i+1);

	/* write outlet points */
	n = 0;
	for(i = 0; i < outlet_cnt; i++)
		n += outlet_mesh[i]->pnum;
	fprintf(fp, "%10lld\n", n);
	for(i = 0; i < outlet_cnt; i++) 
		for(n = 0; n < outlet_mesh[i]->pnum; n++)
			fprintf(fp, "%d %d %d  %d\n",
				outlet_mesh[i]->points[n].i,
				outlet_mesh[i]->points[n].j,
				outlet_mesh[i]->points[n].k,
				inlet_cnt+i+1);

	fclose(fp);
	return;
}

static void print_transform(const char *fname, REAL scale, 
			    int i_tr, int j_tr, int k_tr) {

	FILE		*fp;

	if ((fp = fopen(fname, "w")) == NULL) {
                fprintf(stderr, "%s: cannot open file %s for writing, quitting...\n", __func__, fname);
                return;
        }

	fprintf(fp, "scale_factor: "REAL_SPEC"\n", scale);
	fprintf(fp, "translation_vector: %d,%d,%d\n", i_tr, j_tr, k_tr);

	fclose(fp);
	return;
}

static int find_fluid_inoutlet_dir(mesh_t *m, 
	 			   int mini, int maxi,
				   int minj, int maxj,
			   	   int mink, int maxk,
				   int *di, int *dj, int *dk) {
	long long	i4, ifl;
	int		i, j, k;
	int		ret = -1;

	/*           i j k  internal_face external_face     not_crossed crossed
	 * ijk_sides[0 1 2][            0             1] = [          0       1]
	 */
	int		ijk_sides[3][2] = {{0,0}, {0,0}, {0,0}};

//	LOG(0, "Analyzing inlet/outlet boundary...\n");
	LOG(1, "\tbounding box: [%d, %d] [%d, %d] [%d, %d]\n", 
		mini, maxi, minj, maxj, mink, maxk);

	di[0] = dj[0] = dk[0] = 0;
	/* search through faces parallel to xy */
	for(i = mini; i <= maxi; i++) {
		for(j = minj; j < maxj; j++) {
			
			/* try lower face */	
			i4 = (long long)mink*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);
		
			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type == INTERNAL) {
//					fprintf(stderr, "(0, 0, 1)\n");
					ijk_sides[2][0]++;
				}
			}

			/* try upper face */	
			i4 = (long long)maxk*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);
		
			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type == INTERNAL) {
//					fprintf(stderr, "(0, 0, -1)\n");
					ijk_sides[2][1]++;
				}
			}
		}
	}
	
	/* search through faces parallel to xz */
	for(i = mini; i <= maxi; i++) {
		for(k = mink; k < maxk; k++) {
			
			/* try lower face */	
			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)minj*(m->imax+1) +
                             (long long)(i);
			
			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type == INTERNAL) {
//					fprintf(stderr, "(0, 1, 0)\n");
					ijk_sides[1][0]++;
				}
			}

			/* try upper face */	
			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)maxj*(m->imax+1) +
                             (long long)(i);
	
			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type == INTERNAL) {
//					fprintf(stderr, "(0, -1, 0)\n");
					ijk_sides[1][1]++;
				}
			}
		}
	}
	
	/* search through faces parallel to yz */
	for(j = minj; j <= maxj; j++) {
		for(k = mink; k < maxk; k++) {
			
			/* try lower face */	
			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(mini);
		
			ifl = binsearch(m->points, i4, m->pnum);
			if (ifl >= 0) {
				if (m->points[ifl].type == INTERNAL) {
//					fprintf(stderr, "(1, 0, 0)\n");
					ijk_sides[0][0]++;
				}
			}

			/* try upper face */	
			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(maxi);
	
			ifl = binsearch(m->points, i4, m->pnum);	
			if (ifl >= 0) {
				if (m->points[ifl].type == INTERNAL) {
//					fprintf(stderr, "(-1, 0, 0)\n");
					ijk_sides[0][1]++;
				}
			}
		}
	}

	/* check for errors */
	i = 0;
	if (ijk_sides[0][0]) {
		LOG(0, "\tfluid/boundary interface normal (outward): (1, 0, 0)\n");
		LOG(0, "\tfluid/boundary interface points: %d\n", ijk_sides[0][0]);
		di[0] = 1;
		i++;
		ret = ijk_sides[0][0];
	}
	if (ijk_sides[0][1]) {
		LOG(0, "\tfluid/boundary interface normal (outward): (-1, 0, 0)\n");
		LOG(0, "\tfluid/boundary interface points: %d\n", ijk_sides[0][1]);
		di[0] = -1;
		i++;
		ret = ijk_sides[0][1];
	}
	if (ijk_sides[1][0]) {
		LOG(0, "\tfluid/boundary interface normal (outward): (0, 1, 0)\n");
		LOG(0, "\tfluid/boundary interface points: %d\n", ijk_sides[1][0]);
		dj[0] = 1;
		i++;
		ret = ijk_sides[1][0];
	}
	if (ijk_sides[1][1]) {
		LOG(0, "\tfluid/boundary interface normal (outward): (0, -1, 0)\n");
		LOG(0, "\tfluid/boundary interface points: %d\n", ijk_sides[1][1]);
		dj[0] = -1;
		i++;
		ret = ijk_sides[1][1];
	}
	if (ijk_sides[2][0]) {
		LOG(0, "\tfluid/boundary interface normal (outward): (0, 0, 1)\n");
		LOG(0, "\tfluid/boundary interface points: %d\n", ijk_sides[2][0]);
		dk[0] = 1;
		i++;
		ret = ijk_sides[2][0];
	}
	if (ijk_sides[2][1]) {
		LOG(0, "\tfluid/boundary interface normal (outward): (0, 0, -1)\n");
		LOG(0, "\tfluid/boundary interface points: %d\n", ijk_sides[2][1]);
		dk[0] = -1;
		i++;
		ret = ijk_sides[2][1];
	}

	if (i != 1) {
		fprintf(stderr, "\tInlet/outlet zone cuts the domain with zero or more than one plane!\n");
		return -1;
	} else if (ret == 0) {
		fprintf(stderr, "\tInlet/outlet zone cutted by one plane acroos ZERO nodes!\n");
		return -1;
	}
	return ret;
}

int read_normal(char *fname, char *tgtname, REAL *xn, REAL *yn, REAL *zn) {

	FILE	*f;
	char	buffer[MAX_LINE];
	char	*basename, *line;
	int	found = 0;

	if ((f = fopen(fname, "r")) == NULL) return 0;
	
	basename = tgtname + strlen(tgtname);
	for(; (*basename != '/') && (basename >= tgtname); basename--);
	basename++;

	while((line = get_neline(buffer, MAX_LINE, f))) {
		if (strncmp(line, basename, strlen(basename))) continue;
		sscanf(line + strlen(basename), REAL_SPEC" "REAL_SPEC" "REAL_SPEC, xn, yn, zn);
		found = 1;
		break;
	}

	fclose(f);
	return found;
}

static int find_center_line(mesh_t *m) {

	long long	n, c, nifl;
	int 		i, j, k, curr_layer;
	long long 	i4;

	fprintf(stderr, "Layering mesh...\n");

	/* remove points tagged EXTERNAL, INLET and OUTLET */
	c = 0;
	for(n = 0; n < m->pnum; n++)
		if (m->points[n].type == INTERNAL) {
			m->points[c] = m->points[n];
			m->points[c].type = 0;
			c++;
		}
	/* no need to re-sort (we are removing elements from a sorted list) */

	fprintf(stderr, "\tRemoved %lld EXTERNAL, INLET and OUTLET points.\n", m->pnum - c);
	m->pnum = c;

	/* add external layer of points */
	fprintf(stderr, "\tPeeling layer 1...");
	for(n = 0; n < m->pnum; n++) {
		
		for(i = m->points[n].i-1; i <= m->points[n].i+1; i++) {
		for(j = m->points[n].j-1; j <= m->points[n].j+1; j++) {
		for(k = m->points[n].k-1; k <= m->points[n].k+1; k++) {
			
			if ((i==m->points[n].i)&&(j==m->points[n].j)&&(k==m->points[n].k)) continue;

			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);

			if (binsearch(m->points, i4, m->pnum) < 0) {
				m->points[n].type = 1;
				goto find_center_line_out1;
			}
		}}}
find_center_line_out1:
		continue;
	}
	fprintf(stderr, "done.\n");

	curr_layer = 1;
	do {
		fprintf(stderr, "\tPeeling layer %d...", curr_layer+1);
		c = 0;
		for(n = 0; n < m->pnum; n++) {
			if (m->points[n].type != curr_layer /* untested yet */) continue;
			
			for(i = m->points[n].i-1; i <= m->points[n].i+1; i++) {
			for(j = m->points[n].j-1; j <= m->points[n].j+1; j++) {
	                for(k = m->points[n].k-1; k <= m->points[n].k+1; k++) {

				if ((i==m->points[n].i)&&(j==m->points[n].j)&&(k==m->points[n].k)) continue;

				i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
	                             (long long)j*(m->imax+1) +
				     (long long)(i);

				nifl = binsearch(m->points, i4, m->pnum);
				if (nifl >= 0) {
					if (m->points[nifl].type == 0) {
						m->points[nifl].type = curr_layer + 1;
						c++;
//						goto find_center_line_out2;
					}
				}// else {
//					fprintf(stderr, "\tError: mesh not peeled uniformly!\n");
//					return 0;
//				}
			}}}
//find_center_line_out2:
			continue;
		}
		fprintf(stderr, "done.\n");
		curr_layer++;
	} while (c);

	/* mark points with at least a neighbor in a deeper layer */
/*	for(n = 0; n < m->pnum; n++) {
		int on_center_line = 1;

		for(i = m->points[n].i-1; i <= m->points[n].i+1; i++) {
		for(j = m->points[n].j-1; j <= m->points[n].j+1; j++) {
		for(k = m->points[n].k-1; k <= m->points[n].k+1; k++) {

			if ((i==m->points[n].i)&&(j==m->points[n].j)&&(k==m->points[n].k)) continue;

			i4 = (long long)k*(m->jmax+1)*(m->imax+1) +
                             (long long)j*(m->imax+1) +
                             (long long)(i);

			nifl = binsearch(m->points, i4, m->pnum);
			if (nifl >= 0) {
				if ((m->points[nifl].type & 0x7F) > (m->points[n].type & 0x7F)) {
					on_center_line = 0;
					goto find_center_line_out3;
				}
			}
		}}}
find_center_line_out3:
		if (on_center_line)
			m->points[n].type |= 0x80;
	}
	
	c = 0;
	for(n = 0; n < m->pnum; n++)
		if (m->points[n].type & 0x80) {
			m->points[c] = m->points[n];
			m->points[c].type &= 0x7F;
			c++;
		}
	m->pnum = c;
*/	/* no need to re-sort (we are removing elements from a sorted list) */


	return 1;
}

int main(int argc, char **argv) {

	surface_t	*surface = NULL;
	surface_t	*inlet_surf = NULL;
	surface_t	*outlet_surf = NULL;

	REAL		scale_factor = FPSUFF(1.0E+1);

	mesh_t		*mesh = NULL;
	mesh_t		*inlet_mesh[MAX_INOUTLET];  /* to store inlet points */
	mesh_t		*outlet_mesh[MAX_INOUTLET]; /* to store outlet points */

	int		inlet_cnt = 0, outlet_cnt = 0;
	
	REAL		spacing = ONE;//(SQRT_CPU(THREE));///TWO)+EPSILON;

	char		fname[MAX_FNAME];
	char		normfname[MAX_FNAME];

	char		inlet_fnames[MAX_INOUTLET][MAX_FNAME];
	REAL		inlet_normals[MAX_INOUTLET][3];
	int		inlet_interf[MAX_INOUTLET];
	char		outlet_fnames[MAX_INOUTLET][MAX_FNAME];
	REAL		outlet_normals[MAX_INOUTLET][3];
	int		outlet_interf[MAX_INOUTLET];

	char		myname[MAX_LINE];

	int		c, i, j, k;

	strncpy(myname, argv[0], strlen(argv[0])+1);

	for(c = 0; c < MAX_INOUTLET; c++) {
		inlet_mesh[c] = NULL;
		outlet_mesh[c] = NULL;
	}

	normfname[0] = '\0';

	if (argc < 3) {
		usage(myname);
		goto quit;
	}

	while ( (c = getopt(argc, argv, "f:s:i:o:n:v")) != -1 ) {
		switch (c) {
			case 'f':
				strncpy(fname, (optarg), strlen(optarg)+1);
				break;
			case 's':
				scale_factor = (REAL)atof(optarg);
				break;
			case 'i':
				if (inlet_cnt > MAX_INOUTLET) {
					fprintf(stderr, "Too many inlets, please increase MAX_INOUTLET macro value!\n");
					goto quit;
				}
				strncpy(inlet_fnames[inlet_cnt], (optarg), strlen(optarg)+1);
				inlet_cnt++;
				break;
			case 'o':
				if (outlet_cnt > MAX_INOUTLET) {
					fprintf(stderr, "Too many outlets, please increase MAX_INOUTLET macro value!\n");
					goto quit;
				}
				strncpy(outlet_fnames[outlet_cnt], (optarg), strlen(optarg)+1);
				outlet_cnt++;
				break;
			case 'n':
				strncpy(normfname, (optarg), strlen(optarg)+1);
				break;
			case 'v':
				verblvl++;
				break;
			case '?':
				usage(myname);
				goto quit;
		}
	}

	LOG(0, "Volume file specified:\n\t%s\n", fname);
	fprintf(stderr, "Inlet files specified:\n");
	if (inlet_cnt > 0) {
 		for(c = 0; c < inlet_cnt; c++)
			LOG(0, "\t%s\n", inlet_fnames[c]);
	}
	LOG(0, "Outlet files specified:\n");
	if (outlet_cnt > 0) {
 		for(c = 0; c < outlet_cnt; c++)
			LOG(0, "\t%s\n", outlet_fnames[c]);
	}
	LOG(0, "Normals file specified:\n\t%s\n", normfname);
	LOG(0, "Scaling factor: "REAL_SPEC"\n", scale_factor);

	LOG(0, "\nReading geometry file %s:\n", fname);
	if ((surface = load_surface(fname, scale_factor)) == NULL) {
		goto quit;
	}

	LOG(0, "Generating mesh:\n");
	if ((mesh = generate_mesh_surf(surface, spacing)) == NULL) {
		goto quit;
	}

	LOG(1, "Filling mesh:\n");
	if (!fill_mesh(mesh)) {
		goto quit;
	}

	LOG(1, "Finalizing mesh:\n");
	if (!insert_wall(mesh)) {
		goto quit;
	}

	for(c = 0; c < inlet_cnt; c++) {
		LOG(0, "\nReading inlet file %s:\n", inlet_fnames[c]);
		if ((inlet_surf = load_surface(inlet_fnames[c], scale_factor)) == NULL) {
			goto quit;
		}
		LOG(0, "Analyzing inlet:\n");
		inlet_interf[c] = find_fluid_inoutlet_dir(mesh,
					     MAX(0,mesh->i_tr + ((int)FLOOR_CPU(inlet_surf->xmin))-1),
                                             mesh->i_tr + ((int) CEIL_CPU(inlet_surf->xmax))+1,
                                             MAX(0,mesh->j_tr + ((int)FLOOR_CPU(inlet_surf->ymin))-1),
                                             mesh->j_tr + ((int) CEIL_CPU(inlet_surf->ymax))+1,
                                             MAX(0,mesh->k_tr + ((int)FLOOR_CPU(inlet_surf->zmin))-1),
                                             mesh->k_tr + ((int) CEIL_CPU(inlet_surf->zmax))+1,
					     &i, &j, &k);
		if (inlet_interf[c] < 0) goto quit;

		LOG(1, "Updating mesh:\n");
		inlet_mesh[c] = set_mesh_points(mesh, INLET,
						mesh->i_tr + ((int)FLOOR_CPU(inlet_surf->xmin))-1, 
						mesh->i_tr + ((int) CEIL_CPU(inlet_surf->xmax))+1,
						mesh->j_tr + ((int)FLOOR_CPU(inlet_surf->ymin))-1, 
						mesh->j_tr + ((int) CEIL_CPU(inlet_surf->ymax))+1,
						mesh->k_tr + ((int)FLOOR_CPU(inlet_surf->zmin))-1, 
						mesh->k_tr + ((int) CEIL_CPU(inlet_surf->zmax))+1);
		inlet_mesh[c]->i_tr = i;
		inlet_mesh[c]->j_tr = j;
		inlet_mesh[c]->k_tr = k;

		/* defaults normals to cut plane */
		inlet_normals[c][0] = (REAL)i;
		inlet_normals[c][1] = (REAL)j;
		inlet_normals[c][2] = (REAL)k;

//		fprintf(stderr, "Reading flow normal... ");
		if (!read_normal(normfname, inlet_fnames[c], &inlet_normals[c][0], 
		 	         &inlet_normals[c][1], &inlet_normals[c][2]) ) {
//			fprintf(stderr, "not found, defaulting to interface direction.\n");
			LOG(0, "Flow normal not found, defaulting to interface normal\n");
		} else {
//			fprintf(stderr, "%f %f %f\n", inlet_normals[c][0], inlet_normals[c][1], inlet_normals[c][2]);
			LOG(1, "Flow normal: %f %f %f\n", inlet_normals[c][0], inlet_normals[c][1], inlet_normals[c][2]);
		}
		free(inlet_surf->triangles);
		free(inlet_surf); 
		inlet_surf = NULL;
	}
	for(c = 0; c < outlet_cnt; c++) {
		LOG(0, "\nReading outlet file %s:\n", outlet_fnames[c]);
		if ((outlet_surf = load_surface(outlet_fnames[c], scale_factor)) == NULL) {
			goto quit;
		}
		LOG(0, "Analyzing outlet:\n");
		outlet_interf[c] = find_fluid_inoutlet_dir(mesh,
					     MAX(0,mesh->i_tr + ((int)FLOOR_CPU(outlet_surf->xmin))-1),
                                             mesh->i_tr + ((int) CEIL_CPU(outlet_surf->xmax))+1,
                                             MAX(0,mesh->j_tr + ((int)FLOOR_CPU(outlet_surf->ymin))-1),
                                             mesh->j_tr + ((int) CEIL_CPU(outlet_surf->ymax))+1,
                                             MAX(0,mesh->k_tr + ((int)FLOOR_CPU(outlet_surf->zmin))-1),
                                             mesh->k_tr + ((int) CEIL_CPU(outlet_surf->zmax))+1,
					     &i, &j, &k);
		if (outlet_interf[c] < 0) goto quit;
		
		LOG(1, "Updating mesh:\n");
		outlet_mesh[c] = set_mesh_points(mesh, OUTLET,
						 mesh->i_tr + ((int)FLOOR_CPU(outlet_surf->xmin))-1, 
						 mesh->i_tr + ((int) CEIL_CPU(outlet_surf->xmax))+1,
						 mesh->j_tr + ((int)FLOOR_CPU(outlet_surf->ymin))-1, 
						 mesh->j_tr + ((int) CEIL_CPU(outlet_surf->ymax))+1,
						 mesh->k_tr + ((int)FLOOR_CPU(outlet_surf->zmin))-1, 
						 mesh->k_tr + ((int) CEIL_CPU(outlet_surf->zmax))+1);
		outlet_mesh[c]->i_tr = i;
		outlet_mesh[c]->j_tr = j;
		outlet_mesh[c]->k_tr = k;
		
		/* defaults normals to cut plane */
		outlet_normals[c][0] = (REAL)i;
		outlet_normals[c][1] = (REAL)j;
		outlet_normals[c][2] = (REAL)k;

//		fprintf(stderr, "Reading flow normal... ");
		if (!read_normal(normfname, outlet_fnames[c], &outlet_normals[c][0], 
		 	         &outlet_normals[c][1], &outlet_normals[c][2]) ) {
//			fprintf(stderr, "not found, defaulting to interface direction.\n");
			LOG(0, "Flow normal not found, defaulting to interface normal\n");
		} else {
//			fprintf(stderr, "%f %f %f\n", outlet_normals[c][0], outlet_normals[c][1], outlet_normals[c][2]);
			LOG(1, "Flow normal: %f %f %f\n", outlet_normals[c][0], outlet_normals[c][1], outlet_normals[c][2]);
		}
		free(outlet_surf->triangles);
		free(outlet_surf);
		outlet_surf = NULL;
	}

// Not necessary, leave distance between the mesh and the origin...
//	finalize_mesh(mesh);

	if (!check_closure(mesh)) {
		goto quit;	
	}

	LOG(0, "Writing output files...");
	print_mesh("bgkflag.inp", mesh);
	print_inlet_outlet("inletoutlet.inp",
			      inlet_mesh, inlet_normals, inlet_cnt,
			      outlet_mesh, outlet_normals, outlet_cnt);
	print_transform("mesh_transform.inp", scale_factor, 
			mesh->i_tr, mesh->j_tr, mesh->k_tr);
	LOG(0, "done\n\n");
/* test center line */
//	find_center_line(mesh);
//	print_mesh("pippo.inp", mesh);
quit:
	for(c = 0; c < MAX_INOUTLET; c++) {
		if (inlet_mesh[c]) {
			if (inlet_mesh[c]->points) free(inlet_mesh[c]->points);
			free(inlet_mesh[c]);
		}
		if (outlet_mesh[c]) {
			if (outlet_mesh[c]->points) free(outlet_mesh[c]->points);
			free(outlet_mesh[c]);
		}
	}
	if (inlet_surf) {
		if (inlet_surf->triangles) free(inlet_surf->triangles);
		free(inlet_surf);
	}
	if (outlet_surf) {
		if (outlet_surf->triangles) free(outlet_surf->triangles);
		free(outlet_surf);
	}
	if (mesh) {
		if (mesh->points) free(mesh->points);
		free(mesh);
	}
	return 0;
}
