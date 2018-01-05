#define Nsimp 90 //Must be multiple of 3
#define Nres 50 
#define PI 3.14159265359
#define EA 1.0
#define EI 1.0
#define ZERO 1e-16
#define TOL 1e-8
#define MAX 1e4
#define rho 1.0
#define Nbez 2

struct point {
	double x,y,z;
};
void printpoint(point a) {std::cout << a.x << " " << a.y << " " << a.z << std::endl;}
double normpoint(point a) {return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);}

typedef double (*FunctionHandle) (double t);
FunctionHandle B[] = {Bern0, Bern1, Bern2, Bern3 };
FunctionHandle dB[] = {dBern0, dBern1, dBern2, dBern3 };
FunctionHandle ddB[] = {ddBern0, ddBern1, ddBern2, ddBern3 };
