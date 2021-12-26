#pragma once
#include <cmath>
#include <iostream>

namespace algebra3d{
	struct vector_point { float x, y, z; };
	struct matrix3x3 { vector_point cols[3]; };

	matrix3x3 transpose(const matrix3x3&);
	matrix3x3 multiplyMatrix3x3(const matrix3x3&, const matrix3x3&);
	vector_point multiplyMatrixVector(const matrix3x3&, const vector_point&);

	vector_point multiplyScalar(const vector_point&, float);
	vector_point addVectors(const vector_point&, const vector_point&);
	vector_point addVectors(const vector_point&, const vector_point&, float);

	float dotProduct(const vector_point&, const vector_point&);
	vector_point crossProduct(const vector_point&, const vector_point&);

	float norm2(const vector_point&);	

	void printMatrix(const matrix3x3&);


};

