#include "Polytope.h"
#include <iostream>
#include<vector>
#include<sstream>
#include<string>
#include<cmath>
#include <glut.h>

using namespace std;

const int numPoints = 16;//total no. of  points
vector<vector<float>> vertex_coord = {
		{30, 0},
	{27.07, 15.00},
	{20.61, 23.09},
	{10.61, 26.43},
	{0, 30},
	{-10.61, 26.43},
	{-20.61, 23.09},
	{-27.07, 15.00},
	{-30, 0},
	{-27.07, -15.00},
	{-20.61, -23.09},
	{-10.61, -26.43},
	{0, -30},
	{10.61, -26.43},
	{20.61, -23.09},
	{27.07, -15.00}
	//{5, 5},
	//{ 5,25 },
	//{ 30,25 },
	//{ 30,30 },
	//{ -5,30 },
	//{ -5,5 },
	//{ -25,5 },
	//{ -25,30 },
	//{ -30,30 },
	//{ -30,-5 },
	//{ -5,-5 },
	//{ -5,-25 },
	//{ -30,-25 },
	//{ -30,-30 },
	//{ 5,-30 },
	//{ 5,-5 },
	//{ 25,-5 },
	//{ 25,-30 },
	//{ 30,-30 },
	//{ 30,5 }

};//input data
vector<vector<float>> more_DT = vertex_coord;

vector<Vertex> inputs;
Polytope poly1, poly2;

float Area_Predicate(float v[3][2]) //this is giving area as of now but we only need the sign lol
{
	float area = 0;
	for (int i = 0; i < 3 - 1; i++)
	{

		area = area + (v[i][0] * v[i + 1][1]) - (v[i + 1][0] * v[i][1]); //summation of x(i)y(i+1) - x(i+1)y(i) 
	}
	area = area + (v[2][0] * v[0][1]) - (v[0][0] * v[2][1]);

	area = area / 2;

	return area;
}
//intersection returns if the line  segments intersect
int Intersection(float a[2][2], float b[2][2])
{
	float c1l1t1[3][2] = { {a[0][0], a[0][1]}, {a[1][0], a[1][1]}, {b[0][0], b[0][1]} }; // check 1 and line 1 triangle 1
	float c2l1t2[3][2] = { {a[0][0], a[0][1]}, {a[1][0], a[1][1]}, {b[1][0], b[1][1]} }; // check 2 and line 1 triangle 2
	float c1l2t3[3][2] = { {b[0][0], b[0][1]}, {b[1][0], b[1][1]}, {a[0][0], a[0][1]} }; // check 1 and line  2 triangle 3
	float c2l2t4[3][2] = { {b[0][0], b[0][1]}, {b[1][0], b[1][1]}, {a[1][0], a[1][1]} }; // check 2 and line  2 triangle 4
	float t1 = Area_Predicate(c1l1t1) * Area_Predicate(c2l1t2);

	float t2 = Area_Predicate(c1l2t3) * Area_Predicate(c2l2t4);

	float e = 0.0001; // tolerance
	if (fabs(t1) < e || fabs(t2) < e)
	{
		// one point is on the other line, so that means that it interesects
		if (fabs(t1) < e && fabs(t2) < e)
		//{ // both lines coinside
		//	cout << "both lines coincide\n";
		//	return 1;
		//}
		if (/*(fabs(t1) < e && fabs(t2) < e) ||*/ t2 < 0 || t1 < 0)
		{
			printf("lines intersect at a single point\n");
			return 0; // dioagonal is interesecting with the edge at a single point so it should  be fine
		}
		else
		{
			printf("lines dont intersect\n");
			return 0; // if they dont intersect then its good only
		}
	}
	else if (t1 < 0)
	{
		// line 2 is on the opposite  side of line1
		if (t2 < 0)
		{
			// line1 is on the opposite side of line 2
			printf("lines pakka intersect\n");
			return 1; // lines intersect means it shuold not happen
		}
		else
		{
			printf("lines dont intersect\n");
			return 0; // good only as lines  dont intersect
		}
	}

	// return false;
}

bool EdgeIntersectWithPolygon(Edge* e) {
	float a[2][2] = { {e->v1->x,e->v1->y},{e->v2->x,e->v2->y} };
	bool ans=false;
	float b[2][2];
	for (int i = 0; i < numPoints-1; i++) {
		b[0][0] = vertex_coord[i][0];
		b[0][1] = vertex_coord[i][1];
		b[1][0] = vertex_coord[i+1][0];
		b[1][1] = vertex_coord[i+1][1];

		if (Intersection(a, b) == 1) { //means that they intersect
			ans = true;
			return ans;
		}
		else if(Intersection(a,b)==0){ //means that they dont intersect
			ans = false;
		}
	}
	return ans;

}
int ptInOut(float a[2]) {
	int counter = 0;
	float c[2][2] = { {a[0],a[1]},{a[0],-1000} };
	float b[2][2];
	for (int i = 0; i < numPoints-1; i++) {
		b[0][0] = vertex_coord[i][0];
		b[0][1] = vertex_coord[i][1];
		b[1][0] = vertex_coord[i+1][0];
		b[1][1] = vertex_coord[i+1][1];

		if (Intersection(c, b) == 1) {
			counter++;
		}
	}

	if (counter % 2 == 0) {//event intersection
		return 0; //pt outside
	}
	else {//odd intersection
		return 1; //pt inside
	}
}


int func(vector<vector<float>> vertex)
{
	inputs.clear();
	//Read inputs:
	//string fname = "src/inputs/sphere_less_resolution.txt";
	
	//inputting the 3d points
	float z = 0;
	for (int i = 0; i <vertex.size() ; i++) {
		z = sqrt(pow(vertex[i][0], 2) + pow(vertex[i][1],2) );
		inputs.push_back({ vertex[i][0],vertex[i][1], z});
	}
	//inputParse(fname, inputs);

	/////////////////////////////////////////
	//Actual convex hull generation

	



	return 0;
}

void operationsOnInputsPoly1() {

	poly1.createDoubleTriangle(inputs);
	poly1.addFourthPoint(inputs);
	poly1.addfromList(inputs);


	poly1.computeNormal();
	//poly1.printStatus();
	//poly1.parseToSTL();
}
void operationsOnInputsPoly2() {

	poly2.createDoubleTriangle(inputs);
	poly2.addFourthPoint(inputs);
	poly2.addfromList(inputs);


	poly2.computeNormal();
	//poly1.printStatus();
	//poly1.parseToSTL();
}



void DrawCircle(float cx, float cy, float r, int num_segments)
{
	glBegin(GL_LINE_LOOP);
	for (int ii = 0; ii < num_segments; ii++)
	{
		float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

		float x = r * cosf(theta);//calculate the x component
		float y = r * sinf(theta);//calculate the y component

		glVertex2f(x + cx, y + cy);//output vertex

	}
	glEnd();
}

void updateMore_DT() {
	Edge* e = poly1.getEdgeHead();

	for (int i = 0; i < poly1.nEdge(); i++) {
		if ((*e->adjf1).normal.z < 0 && (*e->adjf2).normal.z < 0) {
			
			float centroid1x = (e->adjf1->v1->x + e->adjf1->v2->x + e->adjf1->v3->x) / (3);
			float centroid1y = (e->adjf1->v1->y + e->adjf1->v2->y + e->adjf1->v3->y) / (3);


			float centroid2x = (e->adjf2->v1->x + e->adjf2->v2->x + e->adjf2->v3->x) / (3);
			float centroid2y = (e->adjf2->v1->y + e->adjf2->v2->y + e->adjf2->v3->y) / (3);

			
			more_DT.push_back({ centroid1x,centroid1y });
			more_DT.push_back({ centroid2x,centroid2y });


		}

		e = e->next;
	}


}
void displayPoly1(void)
{
	glClear(GL_COLOR_BUFFER_BIT);


	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINE_LOOP); //printing the polygon
	for (int i = 0; i < numPoints; i++) {
		glVertex3f(vertex_coord[i][0] / 50, vertex_coord[i][1] / 50, 0.0);
	}
	glEnd();

	int counter = 1;
	Edge* e = poly1.getEdgeHead();

	for (int i = 0; i < poly1.nEdge(); i++) {
		if ((*e->adjf1).normal.z < 0 && (*e->adjf2).normal.z < 0) {
			if(!EdgeIntersectWithPolygon(e)){
				float centroid1x = (e->adjf1->v1->x + e->adjf1->v2->x + e->adjf1->v3->x) / (3 );
				float centroid1y = (e->adjf1->v1->y + e->adjf1->v2->y + e->adjf1->v3->y) / (3 );


				float centroid2x = (e->adjf2->v1->x + e->adjf2->v2->x + e->adjf2->v3->x) / (3 );
				float centroid2y = (e->adjf2->v1->y + e->adjf2->v2->y + e->adjf2->v3->y) / (3 );

				float centroid1[2] = { centroid1x,centroid1y };				float centroid2[2] = { centroid2x,centroid2y };

				cout << "centroid1: " << centroid1x << "," << centroid1y <<"  ptInOut:"<<ptInOut(centroid1)<< endl;
				cout << "centroid2: " << centroid2x << "," << centroid2y << "  ptInOut:" << ptInOut(centroid2) << endl;
				

				if (ptInOut(centroid1) == 1 || ptInOut(centroid2) == 1) {
					std::cout << counter << ". ";
					poly1.customPrint(e);
					glColor3f(1.0, 1.0, 0);
					glBegin(GL_LINES);
					glVertex3f(e->v1->x / 50, e->v1->y / 50, 0.0);
					glVertex3f(e->v2->x / 50, e->v2->y / 50, 0.0);
					glEnd();

					glColor3f(0.0, 1.0, 0);


					//DrawCircle(centroid1x, centroid1y, 0.2, 20);

					//more_DT.push_back({ centroid1x,centroid1y });
					//more_DT.push_back({ centroid2x,centroid2y});
					//
					glPointSize(5.0f);
					glBegin(GL_POINTS);
					if(ptInOut(centroid1)==1){
						glVertex2f(centroid1x / 50, centroid1y / 50);
					}
					if (ptInOut(centroid2) == 1) {
						glVertex2f(centroid2x / 50, centroid2y / 50);
					}
					glEnd();



					counter++;
				}
				
			}
			

		}

		e = e->next;
	}



	glFlush();
}

void displayPoly2(void)
{
	glClear(GL_COLOR_BUFFER_BIT);


	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINE_LOOP); //printing the polygon
	for (int i = 0; i < numPoints; i++) {
		glVertex3f(vertex_coord[i][0] / 50, vertex_coord[i][1] / 50, 0.0);
	}
	glEnd();

	int counter = 1;
	Edge* e = poly2.getEdgeHead();
	
	for (int i = 0; i < poly2.nEdge(); i++) {
		if ((*e->adjf1).normal.z < 0 && (*e->adjf2).normal.z < 0) {
			std::cout << counter << ". ";
			poly2.customPrint(e);
			glColor3f(1.0, 1.0, 0);
			glBegin(GL_LINES);
			glVertex3f(e->v1->x / 50, e->v1->y / 50, 0.0);
			glVertex3f(e->v2->x / 50, e->v2->y / 50, 0.0);
			glEnd();

			glColor3f(0.0, 1.0, 0);
			float centroid1x = (e->adjf1->v1->x + e->adjf1->v2->x + e->adjf1->v3->x)/(3*50);
			float centroid1y = (e->adjf1->v1->y + e->adjf1->v2->y + e->adjf1->v3->y)/(3*50);
			

			float centroid2x = (e->adjf2->v1->x + e->adjf2->v2->x + e->adjf2->v3->x)/(3*50);
			float centroid2y = (e->adjf2->v1->y + e->adjf2->v2->y + e->adjf2->v3->y)/(3*50);
			
			//DrawCircle(centroid1x, centroid1y, 0.2, 20);

			//more_DT.push_back({ centroid1x,centroid1y });
			//more_DT.push_back({ centroid2x,centroid2y});
			//
			glPointSize(5.0f);
			glBegin(GL_POINTS);
			glVertex2f(centroid1x, centroid1y);
			glVertex2f(centroid2x, centroid2y);
			glEnd();

			
			
			counter++;
		
		}

		e = e->next;
	}



	glFlush();
}


void init(void)
{
	/* select clearing color 	*/
	glClearColor(0.0, 0.0, 0.0, 0.0);

	/* initialize viewing values  */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
}

int main(int argc, char** argv)
{
	func(vertex_coord);
	operationsOnInputsPoly1();
	updateMore_DT();
	func(more_DT);
	operationsOnInputsPoly2();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(250, 250);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("hello");
	init();
	glutDisplayFunc(displayPoly1);
	//print_points(vertex_coord);
	//printf("%f",Area_Predicate(vertex_coord));
	glutMainLoop();

	
	return 0;   /* ANSI C requires main to return int. */
}



////function to parse the input file containing coordinates of the points in 3D space
//void inputParse(const std::string& filepath, vector<Vertex>& points)
//{
//	std::ifstream point_source(filepath);
//	std::string line;
//
//	string tmp;
//	double x = 0;
//	double y = 0;
//	double z = 0;
//
//	while (std::getline(point_source, line))
//	{
//		std::stringstream ss;
//		ss << line;
//		while (!ss.eof())
//		{
//			ss >> tmp;
//			x = stof(tmp);
//			ss >> tmp;
//			y = stof(tmp);
//			ss >> tmp;
//			z = stof(tmp);
//		}
//		points.push_back({ x,y,z });
//	}
//
//	return;
//}