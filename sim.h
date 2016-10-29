
enum eltype {offset,vdc1_11,endcap,rf,cnt_eltype};
enum axis {X,Y,Z,cnt_axis};
enum fit {offsetfit,ampfit,cnt_fit};
#define maxionanz 2

extern EXPORT void Auswertung(int stepEval);
extern EXPORT void oldAuswertung(int stepEval);


class EXPORT Kibblesim{
public:
	double lebenszeit[maxionanz];
	int simnum;
	Kibblesim();
	~Kibblesim();
	void LaserXYZ(double v[maxionanz][3], double deltat);
	void Laser(double v[maxionanz][3], double deltat, int axis);
	void Freeze(double t, double v[maxionanz][3]);
	void FreezeRad(double t, double v[maxionanz][3]);
	void FreezeAx(double t, double v[maxionanz][3]);
	void Squeeze(double t, double v[maxionanz][3], double x[maxionanz][3], double a[maxionanz][3], double deltat);
	void TestSqueezeKick(double v[maxionanz][3],double v12[maxionanz][3],double x[maxionanz][3],double &t, double h);
	void NoiseHeat();
	void Restart(double &t);
	int Signum(double zahl);
	double RandomReal(double min, double max);
	double RandomGauss(double mean, double sigma);
	double RandomExp(double tau);
	void PhaseScramble(double v[maxionanz][3], double x[maxionanz][3], double a[maxionanz][3]);
	void StepSqueeze(double v[maxionanz][3],double v12[maxionanz][3],double x[maxionanz][3],double t,double h);
	void Force(double t,double x[maxionanz][3],  double v[maxionanz][3],double a[maxionanz][3], double v12[maxionanz][3],double phi);
	void ForcePseudo(double t,double x[maxionanz][3],  double v[maxionanz][3],double a[maxionanz][3], double v12[maxionanz][3],double phi);
	
	void propagateForwardVerlet(double &t,int n,double x[maxionanz][3],double v[maxionanz][3],double *rms=0,int *iternum=0);
	void GetStartPositions(double t,double x[maxionanz][3], double v[maxionanz][3]);
	
	void initSim();
	void Sim();
	double trapFreq(int axis,int anzosci);
	double scatterCount();


	double rfvoltage;
	double rfvoltageX;
	double a[maxionanz][3];
	double x[maxionanz][3];
	double v[maxionanz][3];
	double v12[maxionanz][3];
	double xtemp[maxionanz][3];
	double vtemp[maxionanz][3];
	double t;
	double t0;
	double phi;
	double HCphi;
	int simstepsMT;
	double h;
	
};

