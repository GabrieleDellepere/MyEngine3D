#pragma once
#include "LinearAlgebra.h"

using namespace algebra3d;

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



