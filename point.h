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
#define TOL1 1e-6
#define TOL2 1e-2

struct point {
	double x,y,z;
};
void printpoint(point a) {std::cout << a.x << " " << a.y << " " << a.z << std::endl;}
double normpoint(point a) {return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);}

typedef double (*FunctionHandle) (double t);
FunctionHandle B[] = {Bern0, Bern1, Bern2, Bern3 };
FunctionHandle dB[] = {dBern0, dBern1, dBern2, dBern3 };
FunctionHandle ddB[] = {ddBern0, ddBern1, ddBern2, ddBern3 };

point S(point p[], double t) {
	point pt = { .x=0.0, .y=0.0, .z=0.0 };
	for (int i=0;i<4;i++) {
		pt.x += p[i].x*B[i](t);
		pt.y += p[i].y*B[i](t);
		pt.z += p[i].z*B[i](t);
	}
	return pt;
}

point dS(point p[], double t) {
	point pt = { .x=0.0, .y=0.0, .z=0.0 };
	for (int i=0;i<4;i++) {
		pt.x += p[i].x*dB[i](t);
		pt.y += p[i].y*dB[i](t);
		pt.z += p[i].z*dB[i](t);
	}
	return pt;
}

point ddS(point p[], double t) {
	point pt = { .x=0.0, .y=0.0, .z=0.0 };
	for (int i=0;i<4;i++) {
		pt.x += p[i].x*ddB[i](t);
		pt.y += p[i].y*ddB[i](t);
		pt.z += p[i].z*ddB[i](t);
	}
	return pt;
}

int f_Length(unsigned ndim, const double *t, void *fdata, unsigned fdim, double *fval)
{
	point p[4] = {((point *) fdata)[0],((point *) fdata)[1],((point *) fdata)[2],((point *) fdata)[3]};
	fval[0] = normpoint(dS(p,t[0]));
	return 0; // success
}

int f_AxForce(unsigned ndim, const double *t, void *fdata, unsigned fdim, double *fval)
{
	point p[4] = {((point *) fdata)[0],((point *) fdata)[1],((point *) fdata)[2],((point *) fdata)[3]};
	point dR = dS(p,t[0]);
	double ds = normpoint(dR);
	
	for(int j=0;j<4;j++) fval[j] = (1.0/ds)*dR.x*dB[j](t[0]);
	for(int j=0;j<4;j++) fval[j+4] = (1.0/ds)*dR.y*dB[j](t[0]);
	for(int j=0;j<4;j++) fval[j+8] = (1.0/ds)*dR.z*dB[j](t[0]);
	
	return 0; // success
}

int f_BeEnergy(unsigned ndim, const double *t, void *fdata, unsigned fdim, double *fval)
{
	point p[4] = {((point *) fdata)[0],((point *) fdata)[1],((point *) fdata)[2],((point *) fdata)[3]};
	point dR = dS(p,t[0]),ddR = ddS(p,t[0]);
	
	double num = pow((ddR.z*dR.y - ddR.y*dR.z),2) + pow((ddR.x*dR.z - ddR.z*dR.x),2) + pow((ddR.y*dR.x - ddR.x*dR.y),2);
	double den = pow(normpoint(dR),5);
	fval[0] = num/den;
	
	return 0; // success
}

int f_BeForce(unsigned ndim, const double *t, void *fdata, unsigned fdim, double *fval)
{
	point p[4] = {((point *) fdata)[0],((point *) fdata)[1],((point *) fdata)[2],((point *) fdata)[3]};
	point dR = dS(p,t[0]),ddR = ddS(p,t[0]);
	double ds = normpoint(dR);
	
	double num = pow((ddR.z*dR.y - ddR.y*dR.z),2) + pow((ddR.x*dR.z - ddR.z*dR.x),2) + pow((ddR.y*dR.x - ddR.x*dR.y),2);
	double Nm[3] = {(ddR.z*dR.y - ddR.y*dR.z), (ddR.x*dR.z - ddR.z*dR.x), (ddR.y*dR.x - ddR.x*dR.y)};
	
	for(int j=0;j<4;j++) {
		double dNdxj = 2*Nm[1]*(ddB[j](t[0])*dR.z - ddR.z*dB[j](t[0])) + 2*Nm[2]*(ddR.y*dB[j](t[0]) - ddB[j](t[0])*dR.y);
		double dsdxj = (dR.x/ds)*dB[j](t[0]);
		fval[j] = (ds*dNdxj - 5*num*dsdxj)/pow(ds,6);
	}
	
	for(int j=0;j<4;j++) {
		double dNdyj = 2*Nm[0]*(ddR.z*dB[j](t[0]) - ddB[j](t[0])*dR.z) + 2*Nm[2]*(ddB[j](t[0])*dR.x - ddR.x*dB[j](t[0]));
		double dsdyj = (dR.y/ds)*dB[j](t[0]);
		fval[j+4] = (ds*dNdyj - 5*num*dsdyj)/pow(ds,6);
	}
	
	for(int j=0;j<4;j++) {
		double dNdzj = 2*Nm[0]*(ddB[j](t[0])*dR.y - ddR.y*dB[j](t[0])) + 2*Nm[1]*(ddR.x*dB[j](t[0]) - ddB[j](t[0])*dR.x);
		double dsdzj = (dR.z/ds)*dB[j](t[0]);
		fval[j+8] = (ds*dNdzj - 5*num*dsdzj)/pow(ds,6);
	}
	
	return 0; // success
}

int f_Mass(unsigned ndim, const double *t, void *fdata, unsigned fdim, double *fval)
{
	point p[4] = {((point *) fdata)[0],((point *) fdata)[1],((point *) fdata)[2],((point *) fdata)[3]};
	point dR = dS(p,t[0]);
	double ds = normpoint(dR);
	
	for(int j=0;j<4;j++) for(int k=0;k<4;k++) fval[j*4+k] = ds*B[j](t[0])*B[k](t[0]);
	return 0; // success
}

int f_InterEnergy(unsigned ndim, const double *t, void *fdata, unsigned fdim, double *fval)
{
	point p1[4] = {((point *) fdata)[0],((point *) fdata)[1],((point *) fdata)[2],((point *) fdata)[3]};
	point p2[4] = {((point *) fdata)[4],((point *) fdata)[5],((point *) fdata)[6],((point *) fdata)[7]};
	
	//printpoint(p2[3]);
	
	point R1 = S(p1,t[0]), R2 = S(p2,t[1]);
	point dR1 = dS(p1,t[0]), dR2 = dS(p2,t[1]);
	double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
	fval[0] = pair_energy(d)*normpoint(dR1)*normpoint(dR2);
	
	//fval[0] = pair_energy(d);
	//printf("\n value %.15f %.15f %.15f %.15f %.15f",t[0],t[1],d,pair_energy(d),fval[0]);
	return 0; // success
}
