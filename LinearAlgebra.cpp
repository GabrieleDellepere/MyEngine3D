#include "LinearAlgebra.h"

using namespace algebra3d;

matrix3x3 algebra3d::transpose(const matrix3x3 &matrix) {
	auto SWAP = [](float &a, float &b) {float aux = a;  a = b; b = aux; }; 
	matrix3x3 result = matrix;

	SWAP(result.cols[0].y, result.cols[1].x);
	SWAP(result.cols[0].z, result.cols[2].x);
	SWAP(result.cols[1].z, result.cols[2].y);

	return result;
}

//potrei aggiungere un "multiplyTransponseMatrixVector" ma devo vedere se ne vale la pena

vector_point algebra3d::multiplyMatrixVector(const matrix3x3& matrix, const vector_point& vector) {
	vector_point result;
	result.x = matrix.cols[0].x * vector.x + matrix.cols[1].x * vector.y + matrix.cols[2].x * vector.z;
	result.y = matrix.cols[0].y * vector.x + matrix.cols[1].y * vector.y + matrix.cols[2].y * vector.z;
	result.z = matrix.cols[0].z * vector.x + matrix.cols[1].z * vector.y + matrix.cols[2].z * vector.z;
	return result;
}

matrix3x3 algebra3d::multiplyMatrix3x3(const matrix3x3 &matrix1, const matrix3x3 &matrix2) {
	matrix3x3 result;
	for (int i = 0; i < 3; i++) {
		result.cols[i] = algebra3d::multiplyMatrixVector(matrix1, matrix2.cols[i]);
	}
	return result;
}

float algebra3d::dotProduct(const vector_point& vector1, const vector_point& vector2) {
	float result;
	result = vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
	return result;
}

vector_point algebra3d::crossProduct(const vector_point &vector1, const vector_point &vector2) {
	vector_point result;
	result.x = (vector1.y * vector2.z - vector2.y * vector1.z);		//det
	result.y = -(vector1.x * vector2.z - vector2.x * vector1.z);	//-det
	result.z = (vector1.x * vector2.y - vector2.x * vector1.y);		//det
	return result;
}

float algebra3d::norm2(const vector_point& vector) {
	return sqrtf(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

vector_point algebra3d::multiplyScalar(const vector_point& vector, float scalar) {
	vector_point result;
	result.x = vector.x * scalar;
	result.y = vector.y * scalar;
	result.z = vector.z * scalar;
	return result;
}

vector_point algebra3d::addVectors(const vector_point& vector1, const vector_point& vector2, float lambda) {
	vector_point result;
	result.x = vector1.x + vector2.x * lambda;
	result.y = vector1.y + vector2.y * lambda;
	result.z = vector1.z + vector2.z * lambda;
	return result;
}
vector_point algebra3d::addVectors(const vector_point& vector1, const vector_point& vector2) {
	return algebra3d::addVectors(vector1, vector2, 1);
}


void algebra3d::printMatrix(const matrix3x3 &matrix) {

	for (int i = 0; i < 3; i++) std::cout << matrix.cols[i].x << " ";
	std::cout << '\n';
	for (int i = 0; i < 3; i++) std::cout << matrix.cols[i].y << " ";
	std::cout << '\n';
	for (int i = 0; i < 3; i++) std::cout << matrix.cols[i].z << " ";
	std::cout << '\n';

}


struct triangleNode {
	vector_point triangle[3];
	triangleNode* next;
};

class TriangleList {
protected:

	triangleNode* first = nullptr;
	triangleNode* last = nullptr;
public:

	void push_back(vector_point tri[3]) {
		triangleNode* newTriangle = new triangleNode;
		newTriangle->triangle[0] = tri[0];
		newTriangle->triangle[1] = tri[1];
		newTriangle->triangle[2] = tri[2];
		newTriangle->next = nullptr;
		if (last == nullptr) first = newTriangle;
		else last->next = newTriangle;
		last = newTriangle;
	}
	vector_point* front() {
		if (first->triangle == nullptr) throw "list is empty";
		return first->triangle;
	}
	void pop_front() {
		triangleNode* toDelete = first;
		if (toDelete == nullptr) throw "list is empty";
		first = first->next;
		if (first == nullptr) last = first;

		delete toDelete;
	}

	int size() {
		triangleNode* aux = first;
		int count = 0;
		while (aux != nullptr) {
			count++;
			aux = aux->next;
		}
		return count;
	}

	void clear() {
		triangleNode* aux = first;
		while (aux != nullptr) {
			first = first->next;
			delete aux;
			aux = first;
		}
		last = nullptr;
	}

};