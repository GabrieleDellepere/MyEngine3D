#include "olcConsoleGameEngine.h"
#include <string>

class GameEngine : public olcConsoleGameEngine {

public:

    struct point {
        float x, y, z;
    };

    struct triangle_pointer {
        short p[3];
    };

    struct triangle {
        short p[3];
        COLOUR texture;
    };

    struct triangle_distance {
        float minZ;
        float maxZ;
        float minDis;

    };

    struct polyhedron {
        std::vector<point> vecPoints;
        std::vector<triangle> shape;

        std::vector<point> tempVecPoints;
        
    };

    struct matrix4x4 {
        float m[4][4];
    };

    class Player {
    public:
        struct angle {
            float yAxis, zxAxis; 
            //yAxis = rotation on y Axis and 
            //zxAxis = rotation on zx Plane 
            //(yeah I kinda did an oopsie and was too lazy to correct the variable name)

        };

        Player() {
            position = { 0.0f, 0.0f, 0.0f };
            velocity = { 0.0f, 0.0f, 0.0f };
            aAngle = { 0.0f, 0.0f };
        }

        void tick(float fElapsedtime) {
            position.x += velocity.x * fElapsedtime;
            position.y += velocity.y * fElapsedtime;
            position.z += velocity.z * fElapsedtime;
        }

        GameEngine::point getPosition() { return position; }

        void setPosition(GameEngine::point p) {  position = p; }

        angle getAngle() { return aAngle; }

        void setAngle(angle aA) {
            aAngle = aA;
            if (aAngle.yAxis < -3.14159f) aAngle.yAxis += 2 * 3.14159f;
            if (aAngle.yAxis > +3.14159f) aAngle.yAxis -= 2 * 3.14159f;

            if (aAngle.zxAxis < -3.14159f / 2.0f) aAngle.zxAxis = -3.14159f / 2.0f;
            if (aAngle.zxAxis > +3.14159f / 2.0f) aAngle.zxAxis = +3.14159f / 2.0f;
            
        }

        void setVelocityX(float fVx) { velocity.x = fVx; }
        void setVelocityY(float fVy) { velocity.y = fVy; }
        void setVelocityZ(float fVz) { velocity.z = fVz; }

        float getVelocityX() { return velocity.x; }
        float getVelocityY() { return velocity.y; }
        float getVelocityZ() { return velocity.z; }

    private:
        GameEngine::point position;
        GameEngine::point velocity; //yeah i'm using a point as a vector, i'm disgusting

        angle aAngle;
    };

    GameEngine() {
        m_sAppName = L"3D Engine";
    }

    int ConstructConsole(int width, int height, int fontw, int fonth) {
        int result = olcConsoleGameEngine::ConstructConsole(width, height, fontw, fonth);
        if (result == 1) { 
            zBuffer.resize(ScreenWidth() * ScreenHeight()); 
        }
        return result;
    }

    virtual bool OnUserCreate() override {

        polyhedron cube;

        cube.vecPoints = {
            {30.0f, 10.0f, 0.0f},   {40.0f, 10.0f, 0.0f},
                {30.0f, 10.0f, 10.0f},   {40.0f, 10.0f, 10.0f},

            {30.0f, 0.0f, 0.0f},   {40.0f, 0.0f, 0.0f},
                {30.0f, 0.0f, 10.0f},   {40.0f, 0.0f, 10.0f}
        };

        cube.shape = {
            {1, 3, 0, FG_DARK_RED}, //TOP
            {0, 3, 2, FG_RED}, //TOP

            {4, 5, 6, FG_DARK_YELLOW}, //BOT
            {6, 5, 7, FG_YELLOW}, //BOT

            {0, 5, 4, FG_DARK_GREEN}, //SOUTH
            {0, 1, 5, FG_GREEN}, //SOUTH

            {2, 3, 6, FG_DARK_BLUE}, //NORTH
            {3, 7, 6, FG_BLUE}, //NORTH

            {0, 2, 4, FG_GREY}, //WEST
            {2, 6, 4, FG_WHITE}, //WEST

            {1, 3, 5, FG_DARK_MAGENTA}, //EAST
            {3, 7, 5, FG_MAGENTA}  //EAST
        
        };

        cube.tempVecPoints.resize(cube.vecPoints.size());

        vecGameObjects.push_back(cube);

        player = new Player();
        return true;
    }

    virtual bool OnUserUpdate(float fElapsedTime) override {

        //get Player data
        point pPlayerPos = player->getPosition();
        Player::angle fPlayerAngle = player->getAngle();

        int nTempX, nTempY;
        nTempX = m_mousePosX;
        nTempY = m_mousePosY;

        if (!m_keys[L'Q'].bHeld) {

            Player::angle tempAngle = player->getAngle();
            tempAngle.yAxis += (-nTempX + nMouseX) * 30.0f * fElapsedTime * 3.14159f / (float)ScreenWidth();
            tempAngle.zxAxis += (-nTempY + nMouseY) * 30.0f * fElapsedTime * 3.14159f / (float)ScreenHeight();


            player->setAngle(tempAngle);
            fPlayerAngle = player->getAngle();

            
        }
        nMouseX = nTempX;
        nMouseY = nTempY;
        if (m_mouse[1].bReleased) {
            int nX, nY, nZ;
            nX = (int)(pPlayerPos.x + 20.0f * cosf(fPlayerAngle.yAxis));
            nY = (int)(pPlayerPos.y + 10.0f * sinf(fPlayerAngle.zxAxis));
            nZ = (int)(pPlayerPos.z + 20.0f * sinf(fPlayerAngle.yAxis));
            polyhedron newCube;
            for (int ix = 0; ix < 2; ix++)
                for (int iy = 0; iy < 2; iy++)
                    for (int iz = 0; iz < 2; iz++) {
                        newCube.vecPoints.push_back({ (float)(nX + 10*ix), (float)(nY + 10*iy), (float)(nZ + 10*iz) });
                    }

            newCube.shape = {
            {1, 3, 0, FG_RED}, //TOP
            {0, 3, 2, FG_RED}, //TOP

            {4, 5, 6, FG_DARK_YELLOW}, //BOT
            {6, 5, 7, FG_DARK_YELLOW}, //BOT

            {0, 5, 4, FG_GREEN}, //SOUTH
            {0, 1, 5, FG_GREEN}, //SOUTH

            {2, 3, 6, FG_BLUE}, //NORTH
            {3, 7, 6, FG_BLUE}, //NORTH

            {0, 2, 4, FG_WHITE}, //WEST
            {2, 6, 4, FG_WHITE}, //WEST

            {1, 3, 5, FG_DARK_MAGENTA}, //EAST
            {3, 7, 5, FG_DARK_MAGENTA}  //EAST

            };

            newCube.tempVecPoints.resize(newCube.vecPoints.size());

            vecGameObjects.push_back(newCube);
        }

        //Check for Player input
        if (m_keys[L'W'].bHeld) {
            player->setVelocityX(20.0f * cosf(fPlayerAngle.yAxis));
            player->setVelocityZ(20.0f * sinf(fPlayerAngle.yAxis));

        }
        else if (m_keys[L'S'].bHeld) {
            player->setVelocityX(20.0f * -cosf(fPlayerAngle.yAxis));
            player->setVelocityZ(20.0f * -sinf(fPlayerAngle.yAxis));
        }
        else {
            player->setVelocityX(0.0f);
            player->setVelocityZ(0.0f);
        }

        if (m_keys[L'A'].bHeld) {
            player->setVelocityX(15.0f * -sinf(fPlayerAngle.yAxis));
            player->setVelocityZ(15.0f * cosf(fPlayerAngle.yAxis));
        }
        else if (m_keys[L'D'].bHeld) {
            player->setVelocityX(15.0f * sinf(fPlayerAngle.yAxis));
            player->setVelocityZ(15.0f * -cosf(fPlayerAngle.yAxis));
        }

        if (m_keys[VK_SPACE].bHeld) {
            player->setVelocityY(30.0f);
        }
        else if (m_keys[VK_SHIFT].bHeld) {
            player->setVelocityY(-30.0f);
        }
        else {
            player->setVelocityY(0.0f);
        }

        player->tick(fElapsedTime);

        

        //adjust points position
        for (polyhedron &poly : vecGameObjects) {

            //We get the posizion of game Objects relative to the player

            ConvertToRelativePos(poly, pPlayerPos, fPlayerAngle);

            //NOW WE HAVE TO CONVERT THE RELATIVE POSITIONS TO SCREEN X AND Y

            ConvertToScreenPoint(poly);

            //scrappedIdea
            //OrderTrianglesByDistance(poly);
            
           
        }

        Render();

        return true;
    }

private:

    float fFov = 3.14159f / 2.0f;   //FOV = 90°
    int nMouseX = 0;
    int nMouseY = 0;

    std::vector<float> zBuffer;

    std::list<polyhedron> vecGameObjects;
    matrix4x4 projectionMatrix;

    Player* player;

    float fAnglesDifference(float fA1, float fA2) {
        float difference = fA1 - fA2;
        if (difference < -3.14159f)
            difference += 2.0f * 3.14159f;
        if (difference > +3.14159f)
            difference -= 2.0f * 3.14159f;

        return difference;
    }

    float fDistance(const point &p1, const point &p2) {
        return sqrtf(
            (p1.x - p2.x) * (p1.x - p2.x) + 
            (p1.y - p2.y) * (p1.y - p2.y) +
            (p1.z - p2.z) * (p1.z - p2.z)
        );
    }

    void multiplyVectorMatrix(const point &p1, point &p2, matrix4x4 matrix) {
        p2.x = p1.x * matrix.m[0][0] + p1.y * matrix.m[1][0] + p1.z * matrix.m[2][0] + matrix.m[3][0];
        p2.y = p1.x * matrix.m[0][1] + p1.y * matrix.m[1][1] + p1.z * matrix.m[2][1] + matrix.m[3][1];
        p2.z = p1.x * matrix.m[0][2] + p1.y * matrix.m[1][2] + p1.z * matrix.m[2][2] + matrix.m[3][2];

        float divider = p1.x * matrix.m[0][3] + p1.y * matrix.m[1][3] + p1.z * matrix.m[2][3] + matrix.m[3][3];

        p2.x /= divider;
        p2.y /= divider;
        //p2.z /= divider;  //why divide z and lose info?

    }

    //we convert each point to relative pos
    void ConvertToRelativePos(polyhedron& poly, const point &pPlayerPos, const Player::angle &fPlayerAngle) {
        int i = 0;
        point relativePos;

        for (point& p : poly.vecPoints) {

            relativePos = { 0, 0, 0 };

            float fX = p.x - pPlayerPos.x;
            float fY = p.y - pPlayerPos.y;
            float fZ = p.z - pPlayerPos.z;

            /*if (fabs(fX) < 0.01f && fabs(fY) < 0.01f && fabs(fZ) < 0.01f) {
                screenPoint.z = -1.0f;
                tempVecPoints.push_back(screenPoint);
                continue;
            }
            else {*/
            float fAngularCoefficient = tanf(fPlayerAngle.yAxis);

            relativePos.x = fabs(fAngularCoefficient * fX - fZ) / sqrtf(1.0f + fAngularCoefficient * fAngularCoefficient);
            float fAngleY = atanf(fZ / fX);
            if (fX < 0) fAngleY += 3.14159f;

            //projection on ZX plane's bisector
            float fProjectionOnZX = sqrtf(fX * fX + fZ * fZ - relativePos.x * relativePos.x);

            if (fAnglesDifference(fAngleY, fPlayerAngle.yAxis) > 0.0f)
                relativePos.x = -relativePos.x;

            if (fabs(fAnglesDifference(fAngleY, fPlayerAngle.yAxis)) > 3.14159f / 2.0f)
                fProjectionOnZX = -fProjectionOnZX;

            //now we can just calculate the Y distance between our point
            //and the line defined by the angular coefficient zx
            //we can use the zx projection as independet variable and
            //y as the dependent one.
            fAngularCoefficient = tanf(fPlayerAngle.zxAxis);

            relativePos.y = (fAngularCoefficient * fProjectionOnZX - fY) / sqrtf(1.0f + fAngularCoefficient * fAngularCoefficient);

            //relative z is just the distance from the projection of our
            //point (projZX, fY) onto the line defined by the angular coefficient
            relativePos.z = sqrtf(fY * fY + fProjectionOnZX * fProjectionOnZX - relativePos.y * relativePos.y);

            if (fAngularCoefficient > 0.01f) {
                //I'm sure I won't remember this part:
                //I'm considering the zx vector as a line whose perpendicular must be below/above (depending
                //on the coefficient sign) our point, if so our relZ is positive otherwise it's negative and the point must be rendered *CAREFULLY*

                if (fY - fProjectionOnZX * (-1.0f / fAngularCoefficient) < 0) {
                    relativePos.z = -relativePos.z;
                }
            }
            else if (fAngularCoefficient < -0.01f) {
                if (fY - fProjectionOnZX * (-1.0f / fAngularCoefficient) > 0) {
                    relativePos.z = -relativePos.z;
                }
            }
            else {  //we avoid dividing by 0;
                if (fProjectionOnZX < 0) relativePos.z = -relativePos.z;
            }

            poly.tempVecPoints[i] = relativePos;
            i++;
        }
    }

    /*
    Scraped idea

    bool isBehind(triangle_pointer& t1, triangle_pointer& t2, std::vector<point>& poly) {
        
    }

    void OrderTrianglesByDistance(polyhedron& poly) {
        poly.tempShape.clear();
        //Future update here: just reorder the triangles instead of deleting and reallocating them, the indices don't change
        
        for (triangle_pointer& tri : poly.shape) {

            auto it = poly.tempShape.begin();
            for (; it != poly.tempShape.end(); it++) {
                if (isBehind(tri, *it, poly.tempVecPoints)) 
                    poly.tempShape.insert(it, tri);

            }
            if (it == poly.tempShape.end()) poly.tempShape.insert(it, tri);

        }
    }*/
    
    //we convert each relative point to screen point, except the negative ones: those will be handled separately
    void ConvertToScreenPoint(polyhedron& poly) {
        int i = 0;
        point screenPoint;

        for (point& relativePos : poly.tempVecPoints) {   
            screenPoint = { 0, 0, 0 };

            //MaxZ is the symbolic Max rendering distance, obtained by multiplying cotan of half Fov with the dimension of the screen
            float fMaxZ = ((float)ScreenWidth() / 2.0f) / tanf(fFov / 2.0f);

            if (relativePos.z >= 0.000001f) {

                //RelX : RelZ = ScreenX : MaxZ  and  RelY : RelZ = ScreenY : MaxZ 

                projectionMatrix = { 0 };
                projectionMatrix.m[0][0] = fMaxZ;
                projectionMatrix.m[1][1] = fMaxZ;
                projectionMatrix.m[2][2] = 1.0f;
                projectionMatrix.m[2][3] = 1.0f;

                multiplyVectorMatrix(relativePos, screenPoint, projectionMatrix);

                //the result x will be between -ScreenWith/2 and ScreenWidth/2, so we translate it by half of the screen size
                screenPoint.x += (float)ScreenWidth() / 2.0f;
                //same for the result y
                screenPoint.y += (float)ScreenHeight() / 2.0f;

            }
            else {

                //The point is added as it is. Unfortunately, everytime a triangle
                //has a negative z point all the points must be re-converted to their
                //relative values before rendering. I'm doing it this way, instead of duplicating the points
                //because i bet most of the points on the screen will have positive z.
                screenPoint = relativePos;
                if (screenPoint.z > -0.000001f) screenPoint.z = -0.000001f;
            }

            poly.tempVecPoints[i] = screenPoint;
            i++;

        }
    }

    //used to convert negative points to screen points
    point FixNegativeScreenPoint(point relP, const point& negPoint) {
        point negScreenPoint = { 0, 0, 0 };

        float fMaxZ = ((float)ScreenWidth() / 2.0f) / tanf(fFov / 2.0f);
        std::cout << relP.z << negPoint.z;
        //NOTE: relP is NOT relative, I'm using fMaxZ to turn it back to relativePos
        relP.x = (relP.x - ScreenWidth() / 2.0f) * relP.z / fMaxZ;
        relP.y = (relP.y - ScreenWidth() / 2.0f) * relP.z / fMaxZ;
        std::cout << relP.z << negPoint.z;
        float fAngularCoefficientX = (relP.x - negPoint.x) / (relP.z - negPoint.z);
        float fAngularCoefficientY = (relP.y - negPoint.y) / (relP.z - negPoint.z);

        float zx, zy;

        if (fAngularCoefficientX == tanf(fFov / 2.0f)) 
            zx = (relP.z * fAngularCoefficientX - relP.x) / (2 * fAngularCoefficientX);
        else if (fAngularCoefficientX == -tanf(fFov / 2.0f)) 
            zx = (relP.z * fAngularCoefficientX - relP.x) / (2 * fAngularCoefficientX);

        else {
            zx = (relP.z * fAngularCoefficientX - relP.x) / (fAngularCoefficientX - tanf(fFov / 2.0f));
            if (zx <= 0 || zx > relP.z)
                zx = (relP.z * fAngularCoefficientX - relP.x) / (fAngularCoefficientX + tanf(fFov / 2.0f));

        }

        if (fAngularCoefficientY == tanf(fFov / 2.0f))
            zy = (relP.z * fAngularCoefficientY - relP.y) / (2 * fAngularCoefficientY);
        else if (fAngularCoefficientY == -tanf(fFov / 2.0f))
            zy = (relP.z * fAngularCoefficientY - relP.y) / (2 * fAngularCoefficientY);

        else {
            zy = (relP.z * fAngularCoefficientY - relP.y) / (fAngularCoefficientY - tanf(fFov / 2.0f));
            if (zy <= 0 || zy > relP.z) 
                zy = (relP.z * fAngularCoefficientY - relP.y) / (fAngularCoefficientY + tanf(fFov / 2.0f));
        }

        //if (zx > 0 && zy > 0 && zx < relP.z && zy < relP.z) {
            if (zx > zy) {
                negScreenPoint.x = ((relP.x + fAngularCoefficientX * (zx - relP.z)) / zx) * fMaxZ + ScreenWidth() / 2.0f;
                negScreenPoint.y = ((relP.y + fAngularCoefficientY * (zx - relP.z)) / zx) * fMaxZ + ScreenWidth() / 2.0f;
                negScreenPoint.z = zx;
            }
            else {
                negScreenPoint.x = ((relP.x + fAngularCoefficientX * (zy - relP.z)) / zy) * fMaxZ + ScreenWidth() / 2.0f;
                negScreenPoint.y = ((relP.y + fAngularCoefficientY * (zy - relP.z)) / zy) * fMaxZ + ScreenWidth() / 2.0f;
                negScreenPoint.z = zy;
            }

        //}
        /*else {
            negScreenPoint.z = negPoint.z;

            if (relP.x > 0) negScreenPoint.x = (float)ScreenWidth();
            else negScreenPoint.x = -1.0f;

            if(relP.y > 0) negScreenPoint.y = (float)ScreenWidth();
            else negScreenPoint.y = -1.0f;
        }*/

        return negScreenPoint;

    }

    //we render all the screen points and fill all the fucking triangles
    void Render() {

        float fMaxZ = ((float)ScreenWidth() / 2.0f) / tanf(fFov / 2.0f);
        int nHorizon = (ScreenHeight() / 2) + (player->getPosition().y + fMaxZ * tanf(player->getAngle().zxAxis));

        Fill(0, 0, ScreenWidth(), nHorizon, PIXEL_SOLID, FG_BLUE);
        Fill(0, nHorizon, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

        for (int i = 0; i < ScreenWidth() * ScreenHeight(); i++) zBuffer[i] = 0.0f;   //set all the values in zBuffer to 0

        for (polyhedron& poly : vecGameObjects) {
            for (triangle tri : poly.shape) {
                //triangle_pointer tri = pair.first;
                point p0, p1, p2;
                p0 = poly.tempVecPoints[tri.p[0]];
                p1 = poly.tempVecPoints[tri.p[1]];
                p2 = poly.tempVecPoints[tri.p[2]];


                FillTriangle3D(p0, p1, p2, tri.texture);

                /*DrawLine(p0, p1);
                DrawLine(p1, p2);
                DrawLine(p2, p0);*/
                

            }
        }
    }

    void DrawLine(point p1, point p2) {
        if (p1.z > 0 && p2.z > 0) {
            if(!(p1.y < 0 || p1.y >= ScreenHeight() || p1.x < 0 || p1.x >= ScreenWidth()
                || p2.y < 0 || p2.y >= ScreenHeight() || p2.x < 0 || p2.x >= ScreenWidth()
                ))
            if(p1.z > zBuffer[(int)p1.y * ScreenWidth() + (int)p1.x] && p2.z > zBuffer[(int)p2.y * ScreenWidth() + (int)p2.x])
            olcConsoleGameEngine::DrawLine((int)p1.x, (int)p1.y, (int)p2.x, (int)p2.y, PIXEL_SOLID, FG_WHITE);
        }
        else {
            if (p1.z <= 0 && p2.z <= 0) return;

            point negScreenPoint;

            //swap points because i don't want to re-write the code for the other scenario
            if (p1.z > 0) { 
                negScreenPoint = FixNegativeScreenPoint(p1, p2);
                //if (negScreenPoint.z < 0) return;
                olcConsoleGameEngine::DrawLine((int)negScreenPoint.x, (int)negScreenPoint.y, (int)p1.x, (int)p1.y, PIXEL_SOLID, FG_WHITE);

            }
            else {
                negScreenPoint = FixNegativeScreenPoint(p2, p1);
                //if (negScreenPoint.z < 0) return;
                olcConsoleGameEngine::DrawLine((int)negScreenPoint.x, (int)negScreenPoint.y, (int)p2.x, (int)p2.y, PIXEL_SOLID, FG_WHITE);
            }

            
        }
    }

    void FillTriangle2D( point p1, point p2, point p3, COLOUR texture) {
        auto SWAP = [](auto& p, auto& q) {auto r = p; p = q; q = r; };
        auto MIN = [](float& f1, float& f2) {return (f2 < f1) ? f2 : f1; };
        auto MAX = [](float& f1, float& f2) {return (f2 > f1) ? f2 : f1; };

        //sort the points by y
        if (p1.y > p2.y) SWAP(p1, p2);
        if (p1.y > p3.y) SWAP(p1, p3);
        if (p2.y > p3.y) SWAP(p2, p3);

        float fAngularCoefficient12, fAngularCoefficient23, fAngularCoefficient31;

        //here must check for when angular coeff -> infinity (rememebr to just set it to k*screenW)
        if (p2.y - p1.y < 0.0001f) fAngularCoefficient12 = (float)(p2.x - p1.x) + 1;
        else fAngularCoefficient12 = float(p2.x - p1.x) / float(p2.y - p1.y);
        if (p3.y - p2.y < 0.0001f) fAngularCoefficient23 = (float)(p3.x - p2.x) + 1;
        else fAngularCoefficient23 = float(p3.x - p2.x) / float(p3.y - p2.y);

        //we could actually return in this scenario, I'm leaving it because
        //this way FillTriangle can also draw a line
        if (p3.y - p1.y < 0.0001f) fAngularCoefficient31 = (float)(p1.x - p3.x) + 1;
        else fAngularCoefficient31 = float(p3.x - p1.x) / float(p3.y - p1.y);

        float fCoefficientZero, fCoefficientEnd;

        //we fill the first half of the triangle
        fCoefficientZero = MIN(fAngularCoefficient12, fAngularCoefficient31);
        fCoefficientEnd = MAX(fAngularCoefficient12, fAngularCoefficient31);

        float lambda1, lambda2, tmpZ;
        point relP[2];
        relP[0] = { p2.x - p1.x, p2.y - p1.y, p2.z - p1.z };
        relP[1] = { p3.x - p1.x, p3.y - p1.y, p3.z - p1.z };

        int index;

        for (int dy = 0; dy < (int)(p2.y - p1.y + 1.0f); dy++) {
            for (int dx = (int)(fCoefficientZero * dy + 0.5f); dx < (int)(fCoefficientEnd * dy + 1.0f); dx++) {
                //we then solve a 2x2 linear system to know which constants make (dx,dy) a linear composition
                //of p1(2d) and p2(2d), i'm calling the costants lambda            
                if (((int)(p1.x + 0.5f) + dx) < 0 || ((int)(p1.x + 0.5f) + dx) >= ScreenWidth()) continue;
                if (((int)(p1.y + 0.5f) + dy) < 0 || ((int)(p1.y + 0.5f) + dy) >= ScreenHeight()) continue;

                lambda2 = ((float)(dy) - (float)(dx)*(relP[0].y / relP[0].x)) / (relP[1].y - relP[1].x * (relP[0].y / relP[0].x));
                lambda1 = ((float)(dx) - lambda2 * relP[1].x) / relP[0].x;
                tmpZ = p1.z + lambda1 * relP[0].z + lambda2 * relP[1].z;

                index = ((int)(p1.y + 0.5f) + dy) * ScreenWidth() + (int)(p1.x + 0.5f) + dx;

                if (tmpZ > 0 && (tmpZ <= zBuffer[index] || zBuffer[index] == 0)) { //if zBuffer[index] == 0 then anything goes   
                    
                    zBuffer[index] = tmpZ;
                    //if every check goes well then we can draw the point
                    Draw(((int)(p1.x + 0.5f) + dx), ((int)(p1.y + 0.5f) + dy), PIXEL_SOLID, texture);
                } 
                
            }
        }

        int nStartX, nOffsetX;

        nOffsetX = (int)((fCoefficientEnd - fCoefficientZero) * float(p2.y - p1.y) + 1);
        if (p2.y - p1.y < 0.0001f) nOffsetX = (int)(abs(p2.x - p1.x));

        if (fCoefficientZero == fAngularCoefficient12) {
            fCoefficientZero = fAngularCoefficient23;
            nStartX = (int)(p2.x + 0.5f);

        }
        else {
            fCoefficientEnd = fAngularCoefficient23;
            nStartX = (int)(p2.x - nOffsetX + 0.5f);

        }

        for (int dy = 0; dy < (int)(p3.y - p2.y + 0.5f); dy++) {
            for (int dx = (int)(fCoefficientZero * dy + 0.5f); dx < (int)(fCoefficientEnd * dy + 0.5f + nOffsetX); dx++) {
                if ((nStartX + dx) < 0 || (nStartX + dx) >= ScreenWidth()) continue;
                if (((int)(p2.y + 0.5f) + dy) < 0 || ((int)(p2.y + 0.5f) + dy) >= ScreenHeight()) continue;

                lambda2 = (relP[0].y + (float)dy - ((float)nStartX + (float)dx - p1.x) * (relP[0].y / relP[0].x)) / (relP[1].y - relP[1].x * relP[0].y / relP[0].x);
                lambda1 = (((float)nStartX + (float)dx - p1.x) - (lambda2 * relP[1].x)) / relP[0].x;
                tmpZ = p1.z + lambda1 * relP[0].z + lambda2 * relP[1].z;

                index = (((int)(p2.y + 0.5f) + dy) * ScreenWidth() + nStartX + dx);
                if (tmpZ > 0 && (tmpZ <= zBuffer[index] || zBuffer[index] == 0)) { //if zBuffer[index] == 0 then anything goes   

                    zBuffer[index] = tmpZ;
                    Draw(nStartX + dx, ((int)(p2.y + 0.5f) + dy), PIXEL_SOLID, texture);
                }
                
            }
        }
    }

    void FillTriangle3D(point p1, point p2, point p3, COLOUR texture=FG_BLACK) {
        auto SWAP = [](point &p, point &q) {point r = p; p = q; q = r; };

        if (p1.z > 0 && p2.z > 0 && p3.z > 0) {
        
            FillTriangle2D(p1, p2, p3, texture);

        }
        else {
            //sort the points by z
            if (p1.z > p2.z) SWAP(p1, p2);
            if (p2.z > p3.z) SWAP(p2, p3);
            if (p1.z > p2.z) SWAP(p1, p2);


            //all points have negative z, no reason to draw the triangle then
            if (p3.z < 0) return;   

            if (p2.z < 0) { //it means both p1.z and p2.z are < 0
                p1 = FixNegativeScreenPoint(p3, p1);
                p2 = FixNegativeScreenPoint(p3, p2);

                //FillTriangle2D(p1, p2, p3, texture);

            }
            else {  //it means only p1.z < 0
               
                point p4 = { 0.0f, 0.0f, 0.0f };
                p4 = FixNegativeScreenPoint(p2, p1);
                p1 = FixNegativeScreenPoint(p3, p1);
                //FillTriangle2D(p1, p2, p3, texture);
                //FillTriangle2D(p4, p1, p3, texture);

            }
        }
    }

};

int main()
{
    GameEngine game;
    game.ConstructConsole(200, 200, 2, 2);
    game.Start();
    return 0;
}
