#include "utils.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

void set_vector(vector v, double x, double y, double z) {
	v[0] = x;
	v[1] = y;
	v[2] = z;
}


void cross(vector v1, vector v2, vector res) {
	res[0] = v1[1] * v2[2] - v1[2] * v2[1];
	res[1] = v1[2] * v2[0] - v1[0] * v2[2];
	res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


void normalize(vector v) {
	double norm = sqrt(SCALAR(v, v));
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}


double determinant(double (*m)[3]) {
	double det = 0.;

	det += (m[0][0])*(m[1][1])*(m[2][2]);
	det += (m[0][1])*(m[1][2])*(m[2][0]);
	det += (m[0][2])*(m[1][0])*(m[2][1]);
	det -= (m[0][2])*(m[1][1])*(m[2][0]);
	det -= (m[0][0])*(m[1][2])*(m[2][1]);
	det -= (m[0][1])*(m[1][0])*(m[2][2]);

	return det;
}

void gram_schmidt(vector v1, vector v2, vector v3) {
	double v1_norm2 = SCALAR(v1, v1);
	double v2_v1 = SCALAR(v2, v1);
	double buffer1 = v2_v1 / v1_norm2;

	v2[0] -= buffer1 * v1[0];
	v2[1] -= buffer1 * v1[1];
	v2[2] -= buffer1 * v1[2];

	double v3_v1 = SCALAR(v3, v1);
	double v3_v2 = SCALAR(v3, v2);
	double v2_norm2 = SCALAR(v2, v2);
	buffer1 = v3_v1 / v1_norm2;
	double buffer2 = v3_v2 / v2_norm2;

	v3[0] -= buffer1*v1[0] + buffer2*v2[0];
	v3[1] -= buffer1*v1[1] + buffer2*v2[1];
	v3[2] -= buffer1*v1[2] + buffer2*v2[2];

	normalize(v1);
	normalize(v2);
	normalize(v3);
}



void matrix_matrix_multiplication(matrix m, matrix n, matrix res) {
	res[0][0] = m[0][0]*n[0][0] + m[0][1]*n[1][0] + m[0][2]*n[2][0];
	res[0][1] = m[0][0]*n[0][1] + m[0][1]*n[1][1] + m[0][2]*n[2][1];
	res[0][2] = m[0][0]*n[0][2] + m[0][1]*n[1][2] + m[0][2]*n[2][2];

	res[1][0] = m[1][0]*n[0][0] + m[1][1]*n[1][0] + m[1][2]*n[2][0];
	res[1][1] = m[1][0]*n[0][1] + m[1][1]*n[1][1] + m[1][2]*n[2][1];
	res[1][2] = m[1][0]*n[0][2] + m[1][1]*n[1][2] + m[1][2]*n[2][2];

	res[2][0] = m[2][0]*n[0][0] + m[2][1]*n[1][0] + m[2][2]*n[2][0];
	res[2][1] = m[2][0]*n[0][1] + m[2][1]*n[1][1] + m[2][2]*n[2][1];
	res[2][2] = m[2][0]*n[0][2] + m[2][1]*n[1][2] + m[2][2]*n[2][2];
}



void get_perpendicular_versor(vector v, vector res) {
	random_vector_on_sphere(res);

	double v_norm2 = SCALAR(v, v);
	double res_v = SCALAR(res, v);
	double buffer = res_v / v_norm2;

	res[0] -= buffer * v[0];
	res[1] -= buffer * v[1];
	res[2] -= buffer * v[2];

	normalize(res);
}


double biased_angle(double theta_max) {

	double expo = 1.0 / 3.0;
	double u = drand48();
	return theta_max * pow(u, expo);

}

void random_vector_on_sphere(vector res) {
	double ransq;
	double ran1, ran2;
	double ranh;

	do {
		ran1 = 1. - 2.*drand48();
		ran2 = 1. - 2.*drand48();
		ransq = SQR(ran1) + SQR(ran2);
	} while(ransq >= 1.);

	ranh = 2.*sqrt(1. - ransq);

	res[0] = ran1*ranh;
	res[1] = ran2*ranh;
	res[2] = 1. - 2.*ransq;
}


void unitary_quaternion(quaternion v) {

	int i;
	double u, z, w, norm;

	u = drand48();
	z = drand48();
	w = drand48();

  v[0] = sqrt(1.0-u) * sin(2.0*M_PI*z);
  v[1] = sqrt(1.0-u) * cos(2.0*M_PI*z);
  v[2] = sqrt(u) * sin(2.0*M_PI*w);
  v[3] = sqrt(u) * cos(2.0*M_PI*w);

	norm = SQR(v[0]) + SQR(v[1]) + SQR(v[2]) + SQR(v[3]);

	for (i=0; i<4; i++) v[i] /= norm;

}



void random_orientation(System *syst, PatchyParticle *p) {

	matrix R;

	// THE FIRST APPROACH IS WRONG: it does not sample uniformly the surface of the unitary sphere
	// To properly sample the unitary sphere, the angular variable needs a non-uniform distribution scaling as ~sin^2 theta /2 with theta \in [0, pi]
	// See the end of "On random rotations in R^3", Miles, R.E., Biometrika, 52, 1965, https://www.jstor.org/stable/2333716?origin=crossref&seq=1

	// vector axis;
	// double angle;
	// random_vector_on_sphere(axis);
	// angle = drand48() * 2 * M_PI;
	// get_rotation_matrix(axis, angle, R);

	// This methos samples the surface of the unitary sphere uniformly
	quaternion v;
	unitary_quaternion(v);
	get_rotation_matrix_from_quaternion(v, R);


	if (syst->n_patches == 1) {
		vector v1 = {0, 0, syst->a};
		MATRIX_VECTOR_MULTIPLICATION(R, v1, p->patches[0]);
	}
	if (syst->n_patches == 2) {
		vector v1 = {0, 0, syst->a};
		vector v2 = {0, 0, -syst->a};
		MATRIX_VECTOR_MULTIPLICATION(R, v1, p->patches[0]);
		MATRIX_VECTOR_MULTIPLICATION(R, v2, p->patches[1]);
	}
	if (syst->n_patches == 3) {
		double val = M_PI / 3.0;
		vector v1 = {syst->a*cos(val), syst->a*sin(val), 0};
		vector v2 = {syst->a*cos(2*val), syst->a*sin(2*val), 0};
		vector v3 = {syst->a*cos(3*val), syst->a*sin(3*val), 0};
		MATRIX_VECTOR_MULTIPLICATION(R, v1, p->patches[0]);
		MATRIX_VECTOR_MULTIPLICATION(R, v2, p->patches[1]);
		MATRIX_VECTOR_MULTIPLICATION(R, v3, p->patches[2]);
	}
	if (syst->n_patches == 4) {
		double half_isqrt3 = 1.0 / sqrt(3);
		vector v1 = {-half_isqrt3, -half_isqrt3,  half_isqrt3};
		vector v2 = {half_isqrt3, -half_isqrt3, -half_isqrt3};
		vector v3 = {half_isqrt3,  half_isqrt3,  half_isqrt3};
		vector v4 = {-half_isqrt3,  half_isqrt3, -half_isqrt3};
		MATRIX_VECTOR_MULTIPLICATION(R, v1, p->patches[0]);
		MATRIX_VECTOR_MULTIPLICATION(R, v2, p->patches[1]);
		MATRIX_VECTOR_MULTIPLICATION(R, v3, p->patches[2]);
		MATRIX_VECTOR_MULTIPLICATION(R, v4, p->patches[2]);
	}


}



void get_rotation_matrix(vector axis, double t, matrix rotation_matrix) {
	double st = sin(t);
	double ct = cos(t);
	double uct = 1. - ct;

	rotation_matrix[0][0] = SQR(axis[0]) + (1. - SQR(axis[0]))*ct;
	rotation_matrix[0][1] = axis[0]*axis[1]*uct - axis[2]*st;
	rotation_matrix[0][2] = axis[0]*axis[2]*uct + axis[1]*st;

	rotation_matrix[1][0] = axis[0]*axis[1]*uct + axis[2]*st;
	rotation_matrix[1][1] = SQR(axis[1]) + (1. - SQR(axis[1]))*ct;
	rotation_matrix[1][2] = axis[1]*axis[2]*uct - axis[0]*st;

	rotation_matrix[2][0] = axis[0]*axis[2]*uct - axis[1]*st;
	rotation_matrix[2][1] = axis[1]*axis[2]*uct + axis[0]*st;
	rotation_matrix[2][2] = SQR(axis[2]) + (1. - SQR(axis[2]))*ct;
}


void get_rotation_matrix_from_quaternion(quaternion v, matrix rotation_matrix) {

	double w, x, y, z;
	w = v[0];
	x = v[1];
	y = v[2];
 	z = v[3];

	rotation_matrix[0][0] = 1.0 - 2.0*(SQR(y) + SQR(z));
	rotation_matrix[0][1] = 2.0*(x*y - z*w);
	rotation_matrix[0][2] = 2.0*(x*z + y*w);

	rotation_matrix[1][0] = 2.0*(x*y + z*w);
	rotation_matrix[1][1] = 1.0 - 2.0*(SQR(x) + SQR(z));
	rotation_matrix[1][2] = 2.0*(y*z - x*w);

	rotation_matrix[2][0] = 2.0*(x*z - y*w);
	rotation_matrix[2][1] = 2.0*(y*z + x*w);
	rotation_matrix[2][2] = 1.0 - 2.0*(SQR(x) + SQR(y));

}


void utils_rotate_matrix(matrix orient_old, matrix orient_new, vector axis, double t) {
	matrix rotation_matrix;
	get_rotation_matrix(axis, t, rotation_matrix);

	matrix_matrix_multiplication(orient_old, rotation_matrix, orient_new);
}


void get_rotated_vector(vector v, vector axis, double t, vector res) {
	matrix rotation_matrix;
	get_rotation_matrix(axis, t, rotation_matrix);
	MATRIX_VECTOR_MULTIPLICATION(rotation_matrix, v, res);
}


void rotate_vector(vector v, vector axis, double t) {
	vector tmp;
	get_rotated_vector(v, axis, t, tmp);
	v[0] = tmp[0];
	v[1] = tmp[1];
	v[2] = tmp[2];
}


void set_orientation_around_vector(vector v, matrix orient, double t) {
	vector y_tmp, w_tmp, z;

	y_tmp[0] = -v[0];
	y_tmp[1] = -v[1];
	y_tmp[2] = -v[2];

	random_vector_on_sphere(z);
	random_vector_on_sphere(w_tmp);

	// orthonormalize
	gram_schmidt(y_tmp, z, w_tmp);

	matrix rotation_matrix;
	get_rotation_matrix(z, t, rotation_matrix);

	vector y, w;
	MATRIX_VECTOR_MULTIPLICATION(rotation_matrix, y_tmp, y);
	MATRIX_VECTOR_MULTIPLICATION(rotation_matrix, w_tmp, w);

	matrix cambiamento_base;

	cambiamento_base[0][0] = z[0];
	cambiamento_base[1][0] = z[1];
	cambiamento_base[2][0] = z[2];

	cambiamento_base[0][1] = w[0];
	cambiamento_base[1][1] = w[1];
	cambiamento_base[2][1] = w[2];

	cambiamento_base[0][2] = y[0];
	cambiamento_base[1][2] = y[1];
	cambiamento_base[2][2] = y[2];

	// rotations have det(R) == 1
	if(determinant(cambiamento_base) < 0) {
		cambiamento_base[0][0] = w[0];
		cambiamento_base[1][0] = w[1];
		cambiamento_base[2][0] = w[2];

		cambiamento_base[0][1] = z[0];
		cambiamento_base[1][1] = z[1];
		cambiamento_base[2][1] = z[2];
	}

	int i, j;
	matrix orient_old;
	for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) orient_old[i][j] = orient[i][j];
	matrix_matrix_multiplication(cambiamento_base, orient_old, orient);
}

void utils_reset_acceptance_counters(System *syst) {
	int i;
	for(i = 0; i < N_MOVES; i++) {
		syst->tries[i] = 0;
		syst->accepted[i] = 0;
	}
}
