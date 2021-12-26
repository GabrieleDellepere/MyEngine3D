// Engine3D_2.0.cpp : Questo file contiene la funzione 'main', in cui inizia e termina l'esecuzione del programma.
//

#include <vector>
#include <string>
#include <fstream>
#include <strstream>
#include "Utils.h"

using namespace algebra3d;

const float tollerance = 0.00001f;

class GameEngine : public olcConsoleGameEngine {
private:

	struct triangle {
		int vertex[3];
		short sym = PIXEL_SOLID;
		short color = FG_WHITE;
	};

	struct mesh {
		std::vector<vector_point> points;
		std::vector<triangle> shape;

		bool loadFromObjectFile(std::string filename) {
		
			std::ifstream f(filename);
			if (!f.is_open()) return false;
				
			while (!f.eof()) {
				char line[128];
				std::strstream s;
				char command;

				f.getline(line, 128);
				command = '\0';
				s << line;

				s >> command;
				if (command == 'v') {	//I hate switch statements btw
					vector_point tmp;
					s >> tmp.x >> tmp.y >> tmp.z;	//I have absolutely no idea how the f it auto converts from string to float but whatever
					points.push_back(tmp);
				}
				else if (command == 'f') {
					triangle tmp;
					s >> tmp.vertex[0] >> tmp.vertex[1] >> tmp.vertex[2];
					tmp.vertex[0]--;		//obj files start counting from 1 instead of 0
					tmp.vertex[1]--;
					tmp.vertex[2]--;	
					std::cout << tmp.vertex[0];
					shape.push_back(tmp);
				}
			}
			return true;
		}
	};

	class Player {
		private:
			matrix3x3 camera;	//col2 è la direzione in cui punta, col0 è la perpendicolare sul piano XZ, col1 è l'altra.

		public:
			vector_point position;
			vector_point velocity;
			vector_point camera_pos;
			vector_point camera_vel;

			Player(vector_point pos) {
				this->position = pos;
				this->camera_pos = pos;
				this->velocity = { 0, 0, 0 };
				this->camera_vel = { 0, 0, 0 };

				setCamera({ 0, 0, 1 });

			}

			matrix3x3 getCamera() { return camera; }

			void setCamera(vector_point dir) {
				float norm = norm2(dir);
				vector_point camera_dir, camera_dirX, camera_dirY;
				camera_dir = { dir.x / norm, dir.y / norm, dir.z / norm };

				camera_dirX.x = sqrtf((camera_dir.z * camera_dir.z) / (camera_dir.x * camera_dir.x + camera_dir.z * camera_dir.z));
				camera_dirX.y = 0;
				camera_dirX.z = sqrtf((camera_dir.x * camera_dir.x) / (camera_dir.z * camera_dir.z + camera_dir.x * camera_dir.x));

				
				if (fabs(dotProduct(camera_dirX, camera_dir)) > tollerance) {
					camera_dirX.x = -camera_dirX.x;
				}

				camera_dirY = crossProduct(camera_dir, camera_dirX);
				if (camera_dirY.y < 0) {
					camera_dirY = multiplyScalar(camera_dirY, -1.0f);
					camera_dirX = multiplyScalar(camera_dirX, -1.0f);
				}

				camera = {camera_dirX, camera_dirY, camera_dir};

			}

			void update() {
				position = addVectors(position, velocity);
				camera_pos = addVectors(camera_pos, camera_vel);
			}

	};

	void check_for_inputs(float fElapsedTime) {
		matrix3x3 rotationXZ = { cosf(3 * fElapsedTime), 0, sin(3 * fElapsedTime), 0, 1, 0, -sin(3 * fElapsedTime), 0, cos(3 * fElapsedTime) };
		matrix3x3 rotationZY = { 1, 0, 0, 0, cosf(3 * fElapsedTime), sin(3 * fElapsedTime), 0, -sin(3 * fElapsedTime), cos(3 * fElapsedTime) };
		rotationZY = multiplyMatrix3x3(player->getCamera(), rotationZY);
		rotationZY = multiplyMatrix3x3(rotationZY, transpose(player->getCamera()));
		vector_point tmp;
		float norm;

		if (m_keys[VK_UP].bHeld) {
			tmp = multiplyMatrixVector(rotationZY, player->getCamera().cols[2]);
			if (!(player->getCamera().cols[2].x * tmp.x < 0) && !(player->getCamera().cols[2].z * tmp.z < 0))
				player->setCamera(tmp);
		}
		if (m_keys[VK_DOWN].bHeld) {
			tmp = multiplyMatrixVector(transpose(rotationZY), player->getCamera().cols[2]);
			if (!(player->getCamera().cols[2].x * tmp.x < 0) && !(player->getCamera().cols[2].z * tmp.z < 0))
				player->setCamera(tmp);
		}
		if (m_keys[VK_RIGHT].bHeld) {
			player->setCamera(multiplyMatrixVector(transpose(rotationXZ), player->getCamera().cols[2]));
		}
		if (m_keys[VK_LEFT].bHeld) {
			player->setCamera(multiplyMatrixVector(rotationXZ, player->getCamera().cols[2]));
		}
		if (m_keys[L'W'].bHeld) {
			tmp = player->getCamera().cols[2];
			tmp.y = 0;
			norm = norm2(tmp);
			tmp.x = tmp.x / norm;
			tmp.z = tmp.z / norm;
			player->velocity = addVectors(player->velocity, multiplyScalar(tmp, 10 * fElapsedTime));
		}
		if (m_keys[L'S'].bHeld) {
			tmp = player->getCamera().cols[2];
			tmp.y = 0;
			norm = norm2(tmp);
			tmp.x = tmp.x / norm;
			tmp.z = tmp.z / norm;
			player->velocity = addVectors(player->velocity, multiplyScalar(tmp, -10 * fElapsedTime));
		}
		if (m_keys[L'D'].bHeld) {
			player->velocity = addVectors(player->velocity, multiplyScalar(player->getCamera().cols[0], 10 * fElapsedTime));
		}
		if (m_keys[L'A'].bHeld) {
			player->velocity = addVectors(player->velocity, multiplyScalar(player->getCamera().cols[0], -10 * fElapsedTime));
		}
		if (m_keys[VK_SPACE].bHeld) {
			player->velocity = addVectors(player->velocity, { 0, -10 * fElapsedTime, 0 });
		}
		if (m_keys[VK_SHIFT].bHeld) {
			player->velocity = addVectors(player->velocity, { 0, 10 * fElapsedTime, 0 });
		}
	}

	vector_point intersectWithPlane(const vector_point& lineStart, const vector_point& lineEnd, const vector_point& normal, const vector_point& p) {

		float d = dotProduct(normal, p);
		vector_point lineVector = addVectors(lineEnd, lineStart, -1);	//this makes it a subtractVectors
		float lambda = (d - dotProduct(normal, lineStart)) / (dotProduct(normal, lineVector));
		return {lineStart.x + lambda * lineVector.x, lineStart.y + lambda * lineVector.y, lineStart.z + lambda * lineVector.z};

	}

	int clipTriangleWithPlane(vector_point triangle[6], vector_point normal, vector_point p) {
		/*auto SWAP = [](vector_point& p1, vector_point& p2, float& a, float& b)
		{vector_point aux = p1; p1 = p2; p2 = aux; float aux2 = a; a = b; b = aux2; };*/

		auto ROTATE = [](vector_point& p1, vector_point& p2, vector_point& p3, float& a, float& b, float& c) {
			vector_point aux = p1; p1 = p2; p2 = p3; p3 = aux; float aux2 = a; a = b; b = c; c = aux2; 
		};
		//normal and p are the components of the plane 

		//do the clipping

		//do the intersection between line and plane
		float d = dotProduct(normal, p);
			
		float a = dotProduct(normal, triangle[0]) - d;
		float b = dotProduct(normal, triangle[1]) - d;
		float c = dotProduct(normal, triangle[2]) - d;

		bool intersects = (a * b < 0 || b * c < 0 || c * a < 0);	//one could be zero i must check them all

		if (!intersects) {
			//two cases: triangle is completely inside or is completely outside. I'm using the normal to decide.
			//Because of this, let's enstablish that the normal should always be pointint INSIDE.
			if (a < 0 || b < 0 || c < 0) {
				//the triangle is outside, let's reject it
				return 0;
			}
			else {	
				//the triangle is completely inside
				return 1;
			}
		}
		else {
			
			/*if (a < b) SWAP(triangle[0], triangle[1], a, b);
			if (b < c) SWAP(triangle[1], triangle[2], b, c);
			if (a < b) SWAP(triangle[0], triangle[1], a, b);*/
			if (a < 0 || (b < 0 && c > 0)) ROTATE(triangle[0], triangle[1], triangle[2], a, b, c);
			if (a < 0) ROTATE(triangle[0], triangle[1], triangle[2], a, b, c);	//by doing this we make sure that c (and occasionaly b) are < 0

			if (b < 0) {
				//two points outside and one inside
				//triangle[0] = triangle[0];
				triangle[1] = intersectWithPlane(triangle[0], triangle[1], normal, p);
				triangle[2] = intersectWithPlane(triangle[0], triangle[2], normal, p);
				return 1;
			}
			else {
				//one point outside and two inside

				//triangle[0] = triangle[0];
				//triangle[1] = triangle[1];
				triangle[5] = intersectWithPlane(triangle[0], triangle[2], normal, p);
				triangle[2] = intersectWithPlane(triangle[1], triangle[2], normal, p);

				triangle[3] = triangle[0];
				triangle[4] = triangle[2];
				

				return 2;
			}
						
		}
			
	}

	void clipTriangles(TriangleList& actual_triangles, int& queueSize, vector_point normal, vector_point p) {
		int newSize = 0;
		for (int i = 0; i < queueSize; i++) {
			vector_point temp[6] = { actual_triangles.front()[0], actual_triangles.front()[1], actual_triangles.front()[2], 0, 0 , 0 };
			actual_triangles.pop_front();
			int num = clipTriangleWithPlane(temp, normal, p);
			if (num > 0) {
				actual_triangles.push_back(temp);
				if (num == 2) { 
					actual_triangles.push_back(temp+3); }
			}
			newSize += num;
		}
		queueSize = newSize;
	}

	void fillTriangle(const vector_point &p1, const vector_point &p2, const vector_point &p3, short c, short col) {
		int x1 = (int)(p1.x + 0.5f); int y1 = (int)(p1.y + 0.5f);
		int x2 = (int)(p2.x + 0.5f); int y2 = (int)(p2.y + 0.5f);
		int x3 = (int)(p3.x + 0.5f); int y3 = (int)(p3.y + 0.5f);
		vector_point p12 = { p2.x - p1.x, p2.y - p1.y, p2.z - p1.z };
		vector_point p13 = { p3.x - p1.x, p3.y - p1.y, p3.z - p1.z };

		auto SWAP = [](int& x, int& y) { int t = x; x = y; y = t; };
		auto drawline = [&](int sx, int ex, int ny) { 
			if (ny < 0 || ny >= ScreenHeight()) return;
			for (int i = sx; i <= ex; i++) {
				
				if (i >= 0 && i < ScreenWidth()) {
					float lambda2, lambda3;
					if (fabs(p12.x) > tollerance) {
						lambda3 = ((ny - p1.y) - (i - p1.x) * p12.y / p12.x) / (p13.y - p13.x * p12.y / p12.x);
						lambda2 = ((i - p1.x) - lambda3 * p13.x) / p12.x;
					}
					else {
						if (fabs(p12.y) <= tollerance || fabs(p13.x) <= tollerance) return;
						lambda3 = (i - p1.x) / p13.x;
						lambda2 = ((ny - p1.y) - lambda3 * p13.y) / p12.y;
					}
					if (lambda2 < tollerance) { lambda2 = tollerance; }
					if (lambda3 < tollerance) { lambda3 = tollerance; }
					if (lambda2 > 1 - tollerance) lambda2 = 1 - tollerance;
					if (lambda3 > 1 - tollerance) lambda3 = 1 - tollerance;
						

					float z = p1.z + lambda2 * p12.z + lambda3 * p13.z;
					if (zBuffer[ny * ScreenWidth() + i] == 0 || zBuffer[ny * ScreenWidth() + i] > z) {
						if (col == FG_RED) {
							std::cout << "";
						}
						zBuffer[ny * ScreenWidth() + i] = z;
						Draw(i, ny, c, col);
					}
				}
				//Draw(i, ny, c, col);
			}
			
		};

		//std::cout << x1 << y1 << x2 << y2 << x3 << y3;

		int t1x, t2x, y, minx, maxx, t1xp, t2xp;
		bool changed1 = false;
		bool changed2 = false;
		int signx1, signx2, dx1, dy1, dx2, dy2;
		int e1, e2;
		// Sort vertices
		if (y1 > y2) { SWAP(y1, y2); SWAP(x1, x2); }
		if (y1 > y3) { SWAP(y1, y3); SWAP(x1, x3); }
		if (y2 > y3) { SWAP(y2, y3); SWAP(x2, x3); }

		t1x = t2x = x1; y = y1;   // Starting points
		dx1 = (int)(x2 - x1); if (dx1 < 0) { dx1 = -dx1; signx1 = -1; }
		else signx1 = 1;
		dy1 = (int)(y2 - y1);

		dx2 = (int)(x3 - x1); if (dx2 < 0) { dx2 = -dx2; signx2 = -1; }
		else signx2 = 1;
		dy2 = (int)(y3 - y1);

		if (dy1 > dx1) {   // swap values
			SWAP(dx1, dy1);
			changed1 = true;
		}
		if (dy2 > dx2) {   // swap values
			SWAP(dy2, dx2);
			changed2 = true;
		}

		e2 = (int)(dx2 >> 1);
		// Flat top, just process the second half
		if (y1 == y2) goto next;
		e1 = (int)(dx1 >> 1);

		for (int i = 0; i < dx1;) {
			t1xp = 0; t2xp = 0;
			if (t1x < t2x) { minx = t1x; maxx = t2x; }
			else { minx = t2x; maxx = t1x; }
			// process first line until y value is about to change
			while (i < dx1) {
				i++;
				e1 += dy1;
				while (e1 >= dx1) {
					e1 -= dx1;
					if (changed1) t1xp = signx1;//t1x += signx1;
					else          goto next1;
				}
				if (changed1) break;
				else t1x += signx1;
			}
			// Move line
		next1:
			// process second line until y value is about to change
			while (1) {
				e2 += dy2;
				while (e2 >= dx2) {
					e2 -= dx2;
					if (changed2) t2xp = signx2;//t2x += signx2;
					else          goto next2;
				}
				if (changed2)     break;
				else              t2x += signx2;
			}
		next2:
			if (minx > t1x) minx = t1x; if (minx > t2x) minx = t2x;
			if (maxx < t1x) maxx = t1x; if (maxx < t2x) maxx = t2x;
			drawline(minx, maxx, y);    // Draw line from min to max points found on the y
										 // Now increase y
			if (!changed1) t1x += signx1;
			t1x += t1xp;
			if (!changed2) t2x += signx2;
			t2x += t2xp;
			y += 1;
			if (y == y2) break;

		}
	next:
		// Second half
		dx1 = (int)(x3 - x2); if (dx1 < 0) { dx1 = -dx1; signx1 = -1; }
		else signx1 = 1;
		dy1 = (int)(y3 - y2);
		t1x = x2;

		if (dy1 > dx1) {   // swap values
			SWAP(dy1, dx1);
			changed1 = true;
		}
		else changed1 = false;

		e1 = (int)(dx1 >> 1);

		for (int i = 0; i <= dx1; i++) {
			t1xp = 0; t2xp = 0;
			if (t1x < t2x) { minx = t1x; maxx = t2x; }
			else { minx = t2x; maxx = t1x; }
			// process first line until y value is about to change
			while (i < dx1) {
				e1 += dy1;
				while (e1 >= dx1) {
					e1 -= dx1;
					if (changed1) { t1xp = signx1; break; }//t1x += signx1;
					else          goto next3;
				}
				if (changed1) break;
				else   	   	  t1x += signx1;
				if (i < dx1) i++;
			}
		next3:
			// process second line until y value is about to change
			while (t2x != x3) {
				e2 += dy2;
				while (e2 >= dx2) {
					e2 -= dx2;
					if (changed2) t2xp = signx2;
					else          goto next4;
				}
				if (changed2)     break;
				else              t2x += signx2;
			}
		next4:

			if (minx > t1x) minx = t1x; if (minx > t2x) minx = t2x;
			if (maxx < t1x) maxx = t1x; if (maxx < t2x) maxx = t2x;
			drawline(minx, maxx, y);
			if (!changed1) t1x += signx1;
			t1x += t1xp;
			if (!changed2) t2x += signx2;
			t2x += t2xp;
			y += 1;
			if (y > y3) return;
		}
	}

	
	const float FOV = 3.14159f / 2.0f;
	const float sinfov2 = sinf(FOV / 2.0f);
	const float cosfov2 = cosf(FOV / 2.0f);

	Player* player;
//temp: will eventually be updated to gameObjects
	std::vector<mesh> objects;

	matrix3x3 screenMatrix;
	
	float* zBuffer;

public:

	GameEngine() {
		player = nullptr;
		zBuffer = nullptr;
		screenMatrix = { 0 };
	}

	virtual bool OnUserCreate() {
		player = new Player({ 0, 0, -5 });
		zBuffer = new float[ScreenWidth() * ScreenHeight()];

		float maxZ = (0.5f * ScreenWidth()) / tanf(FOV / 2.0f);
		screenMatrix.cols[0].x = maxZ;
		screenMatrix.cols[1].y = maxZ;
		screenMatrix.cols[2].z = 1.0f;

		/*mesh cube;
		cube.points = {
			
			{3.0f, 0.0f, 0.0f},   {4.0f, 0.0f, 0.0f},
				{3.0f, 0.0f, 1.0f},   {4.0f, 0.0f, 1.0f},

			{3.0f, 1.0f, 0.0f},   {4.0f, 1.0f, 0.0f},
				{3.0f, 1.0f, 1.0f},   {4.0f, 1.0f, 1.0f}
		};

		cube.shape = {
			{0, 3, 1, FG_DARK_RED}, //TOP
			{3, 2, 0, FG_RED}, //TOP

			{5, 4, 6, FG_DARK_YELLOW}, //BOT
			{6, 7, 5, FG_YELLOW}, //BOT

			{0, 4, 5, FG_DARK_GREEN}, //SOUTH
			{5, 1, 0, FG_GREEN}, //SOUTH

			{6, 2, 3, FG_DARK_BLUE}, //NORTH
			{6, 7, 3, FG_BLUE}, //NORTH

			{4, 0, 2, FG_GREY}, //WEST
			{6, 4, 2, FG_WHITE}, //WEST

			{3, 1, 5, FG_DARK_MAGENTA}, //EAST
			{5, 7, 3, FG_MAGENTA}  //EAST

		};*/

		
		mesh cube;
		cube.loadFromObjectFile("res/cube.obj");
		for (unsigned short i = 0; i < cube.shape.size(); i++) cube.shape[i].color = i+1;
		int blocksize = 2;
		for (int x = -4; x < 5; x++) {
			for (int z = -4; z < 5; z++) {
				mesh block = cube;
				for (vector_point& p : block.points) {
					p.x += x * blocksize;
					p.z += z * blocksize;
				}
				objects.push_back(block);
			}
		}
		
		//objects.push_back(cube);
		return true;
	}

	virtual bool OnUserUpdate(float fElapsedTime) {
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);
		memset(zBuffer, 0, sizeof(float) * ScreenWidth() * ScreenHeight());
		//for (int i = 0; i < ScreenWidth() * ScreenHeight(); i++) zBuffer[i] = 0.0f;

		player->velocity = { 0, 0, 0 };

		check_for_inputs(fElapsedTime);
		
		player->camera_vel = player->velocity;
		player->update();

		matrix3x3 cameraT = transpose(player->getCamera());	//cameraT è trasposta per facilitare i prodotti riga colonna
		std::vector<vector_point> temp_points;
		TriangleList actual_triangles;
		for (auto &object : objects) {
			temp_points.clear();
			for (vector_point &realPos : object.points) {
				vector_point relativePos = {realPos.x - player->camera_pos.x, 
					realPos.y - player->camera_pos.y, realPos.z - player->camera_pos.z };

				vector_point viewPos = multiplyMatrixVector(cameraT, relativePos);
				/*vector_point screenPos = multiplyMatrixVector(screenMatrix, viewPos);
				screenPos.x = screenPos.x / screenPos.z + 0.5f * ScreenWidth();
				screenPos.y = screenPos.y / screenPos.z + 0.5f * ScreenWidth();*/

				temp_points.push_back(viewPos);

			}

			for (auto &triangle : object.shape) {
				int i = triangle.vertex[0];
				int j = triangle.vertex[1];
				int k = triangle.vertex[2];
				vector_point relativeTriangle[3] = { temp_points[i], temp_points[j], temp_points[k] };
				//let's check if the triangle is visible in the first place, meaning that it's normal must be pointing towards us (or better, our field of view)
				vector_point normal = crossProduct(addVectors(relativeTriangle[0], relativeTriangle[1], -1), addVectors(relativeTriangle[2], relativeTriangle[1], -1));
				float normalNorm = norm2(normal);
				normal.x /= -normalNorm;
				normal.y /= -normalNorm;
				normal.z /= -normalNorm;

				float invisibility = dotProduct(normal, relativeTriangle[0]);

				if (invisibility >= 0) {
					//the triangle is not visible
					continue;
					
				}


				actual_triangles.clear();
				//let's clip with the screen
				vector_point viewPoints[6] = { relativeTriangle[0], relativeTriangle[1], relativeTriangle[2], 0, 0, 0 };
				int num = clipTriangleWithPlane(viewPoints, { 0, 0, 1 }, { 0.0f, 0.0f, 0.1f });

				/*for (int n = 0; n < 3*num; n++) {
					viewPoints[n] = multiplyMatrixVector(screenMatrix, viewPoints[n]);	//translate into screenPos
					viewPoints[n].x = viewPoints[n].x / viewPoints[n].z + 0.5f * ScreenWidth();
					viewPoints[n].y = viewPoints[n].y / viewPoints[n].z + 0.5f * ScreenWidth();
				}*/
				int queueSize = 0;
				if (num > 0) {
					//vector_point temp[3] = { viewPoints[0], viewPoints[1], viewPoints[2] };
					actual_triangles.push_back(viewPoints);
					if (num == 2) { 
						//temp[0] = viewPoints[3]; temp[1] = viewPoints[4]; temp[2] = viewPoints[5];
						actual_triangles.push_back(viewPoints + 3);
					}
					queueSize = actual_triangles.size();
					/*clipTriangles(actual_triangles, queueSize, {0, 1, 0}, {0, 0, 0});
					clipTriangles(actual_triangles, queueSize, {0, -1, 0}, {0, ScreenWidth() - 1.0f, 0});
					clipTriangles(actual_triangles, queueSize, {1, 0, 0}, {0, 0, 0});
					clipTriangles(actual_triangles, queueSize, {-1, 0, 0}, { ScreenWidth() - 1.0f, 0, 0});*/
					clipTriangles(actual_triangles, queueSize, { cosfov2, 0, sinfov2 }, { 0, 0, 0 });
					clipTriangles(actual_triangles, queueSize, { -cosfov2, 0, sinfov2 }, { 0, 0, 0 });
					clipTriangles(actual_triangles, queueSize, { 0, cosfov2, sinfov2 }, { 0, 0, 0 });
					clipTriangles(actual_triangles, queueSize, { 0, -cosfov2, sinfov2 }, { 0, 0, 0 });
				}

				for (int i = 0; i < queueSize; i++) {
					vector_point* points = actual_triangles.front();
					for (int n = 0; n < 3; n++) {
						points[n] = multiplyMatrixVector(screenMatrix, points[n]);	//translate into screenPos
						points[n].x = points[n].x / points[n].z + 0.5f * ScreenWidth();
						points[n].y = points[n].y / points[n].z + 0.5f * ScreenWidth();
					}
					fillTriangle(points[0], points[1], points[2], triangle.sym, triangle.color);
					/*DrawLine(points[0].x, points[0].y, points[1].x, points[1].y, PIXEL_SOLID, FG_WHITE);
					DrawLine(points[1].x, points[1].y, points[2].x, points[2].y, PIXEL_SOLID, FG_WHITE);
					DrawLine(points[2].x, points[2].y, points[0].x, points[0].y, PIXEL_SOLID, FG_WHITE);*/
					actual_triangles.pop_front();
				}
				
			}
		}
		return true;
	}


};



int main()
{
	GameEngine* game = new GameEngine();
	game->ConstructConsole(240, 240, 2, 2);
	game->Start();

	/*matrix3x3 prova = {2, 0, 0 , 0, 1, 0, 0, 0, 3};
	matrix3x3 prova2 = prova;
	prova2.cols[0].y = 1;
	prova2.cols[2].x = 2;
	matrix3x3 prova3 = multiplyMatrix3x3(prova, prova2);
	for (auto a : prova3.cols) {
		std::cout << a.x << " " << a.y << " " << a.z << '\n';
	}*/
}

