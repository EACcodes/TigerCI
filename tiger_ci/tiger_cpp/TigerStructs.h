/*
 * @file TigerStructs.h
 * @brief Sets up tiger namespace for generalization of basis set and integral data
 *
 * Author: Francis Ricci
 * Date created: 1/26/15
 *
 */

#ifndef TIGER_STRUCTS
#define TIGER_STRUCTS

namespace tiger
{
enum basis_type {full, minimal, aux};
enum return_val {SUCCESS = 0, INPUT_FAILURE = -1, MO_PARSE_FAILURE = -2, BASIS_FAILURE = -3,
	ORBITAL_FAILURE = -4, CHOLESKY_FAILURE = -5, DF_FAILURE = -6, TIGER_FAILURE = -7};

typedef struct {
	double x;
	double y;
	double z;
} coords_t;

typedef struct {
	coords_t coords;
	int Z;
} atom_t;
}

#endif
