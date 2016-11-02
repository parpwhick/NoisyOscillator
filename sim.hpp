
//#include "stdafx.h"
#include <math.h>


#include <time.h>
#include <boost/thread.hpp>  
#include <boost/date_time.hpp>  
#include "mutex.h"

#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "boost/random/uniform_real_distribution.hpp"

#include <stdio.h>
#include <string>
#include "exports.h"
#include "mcp.h"
#include "Data1DPlotArray.h"
#include "levmarq.h" 
#include "sim.h"
#include <boost/thread/mutex.hpp>
boost::mutex io_mutex; // The iostreams are not guaranteed to be thread-safe!

///gsl includes
#include <stdlib.h>
 #include <stdio.h>
#include <FullGLWidgetCode.h>
#include "ScriptDll.h"
extern EXPORT FullGLWidgetCode &IonTrapDisplay;
using namespace std;
glob(int,numIons,maxionanz);
double Time=4000./1e6;	//30000./1e6; Gesamtdauer der Simulation in sekunden		//49993./1e6;   //1000./1e6;//
						//ACHTUNG!!! Time muss ausreichend klein gewählt werden, ab einem bestimmten größe des Time/steps verhältnisses steigen sonst die fallenfrequenzen an
						// wenn das verhältnis klein genug ist, gehen die fallenfrequenzen "in sättigung", sprich kleinere Time werte verringern die fallenfrequenz nicht weiter.
int const steps=2000000;	//Anzahl der Simulationssteps
glob(double, hglob,5e-9);   //zeitschritt in sekunden


double xStart[maxionanz][3]={};
double vStart[maxionanz][3]={};

	
// Electron charge     
#undef electronCharge
#undef MCa
#undef qDivM
const double  electronCharge=1.602176565e-19;
//Mass of Calcium
const double  MCa=(40.078*1.660538921e-27);
const double  qDivM=electronCharge/MCa;
#define pi 3.141592653589793	
#define epsilon0 8.854187818E-12	/* epsilon0 +- .000000071E-12 F/m */
const static double forceconstant=(electronCharge*electronCharge/(4.*pi*epsilon0))/MCa;


			
#define C 299792458.
#define mm 1./1000.
#define MHz 1000000.
#define hbar 1.054571726e-34		//in [Js]
const double k_B=1.3806488e-23; // Boltzmann Konstante [J/K]



						//ACHTUNG!!! Time muss ausreichend klein gewählt werden, ab einem bestimmten größe des Time/steps verhältnisses steigen sonst die fallenfrequenzen an
						// wenn das verhältnis klein genug ist, gehen die fallenfrequenzen "in sättigung", sprich kleinere Time werte verringern die fallenfrequenz nicht weiter.

const double T=0.001;//0.1;		//in Kelvin
glob(bool,initialdamp,false);

glob(double,OmegaRF,22.75e6);



///gsl includes

//#include <unistd.h>
//#include "linearPaulTrap_edit.h"
glob(bool,runningcalc,true);

bool refine=true;
class D3world;
class D3electrode;

//sollte ungerade anzahl sein da <---- zu faul war für gerade
#define Anzahl 11
#define LOOPS 100
bool maincalled=false;



//include editable scripts here
#ifndef __CINT__
//#include "linearPaulTrap_edit.cxx"
#else
void load(){int retval;gROOT->ProcessLine(".L Traptest.cpp", &retval );}
void unload(){int retval;gROOT->ProcessLine(".U Traptest.cpp", &retval );}
void u() {unload();}
void l() {load();}

#endif

#define Anzahl 11
#define PI 3.141592653589793

struct Xv{
	int i;
	double *xpos;
    double *ypos;
    double *zpos;
    double *vxpos;
    double *vypos;
    double *vzpos;
    double *tpos;
    double *rpos;
};

const int anz=100;
//double arrx[anz];
//double arrxpot[anz];
//double arrxpotfit[anz];
//EXPORT Data1DPlotArray potplotx(anz,arrx,arrxpot);
//EXPORT Data1DPlotArray potplotfitx(anz,arrx,arrxpotfit);

//double arry[anz];
//double arrypot[anz];
//double arrypotfit[anz];
//EXPORT Data1DPlotArray potploty(anz,arry,arrypot);
//EXPORT Data1DPlotArray potplotfity(anz,arry,arrypotfit);

//double arrz[anz];
//double arrzpot[anz];
//double arrzpotfit[anz];
//EXPORT Data1DPlotArray potplotz(anz,arrz,arrzpot);
//EXPORT Data1DPlotArray potplotfitz(anz,arrz,arrzpotfit);






glob(double,dissipationfactor,1000);
glob(double,initialdamptime,0.0001);
glob(double,temperature,0);
glob(bool,displayRad,false);


glob(int,hello_Sam,1);


///////////////////////////////////
struct data3 {
       size_t n;
	   double *x;
       double * y;
       double * sigma;
     };
     

//import capacitor ramps to fit a voltage model used then in sim
	//double fit_val[cnt_eltype][cnt_axis][cnt_fit];



//double xslope[1000000];
//double yslope[1000000];

//EXPORT Data1DPlotArray slope(1000000,xslope,yslope);
//double yslopefit[1000000];
//EXPORT Data1DPlotArray slopefit(1000000,xslope,yslopefit);
//glob(double,fittedoffset,0);
//glob(double,fittedamp,0);
//glob(double,fittedtau,0);
//glob(double,fittedt0,0);
glob(string,runname,"");

glob(string,name,"");

int cntGlobal=0;


//const double tau=7.1e-9;		//Lebensdauer
glob(double, tau,7.1e-9); 
double gamma=1/(2*PI*tau);
EXPORT void setTau(){
	gamma=1/(2*PI*tau);
}
EXPORT void getTau(){
	double tautemp=7.2*0.000000001;
	cout<<"tau "<<tautemp<<endl;
	tauSet(tautemp);
	setTau();
}
glob(double,detune,-1.e7);	//Verstimmung kühllaser
glob(double,sat,5);//=I/I0

EXPORT double setGammaHalf(){
	//detuneSet(-2e8);
	detuneSet(-gamma/2);
	satSet(5);
	return gamma/2;
}



glob(double,rfamp,0);
glob(bool,laseron,true);
glob(bool,laserTrap,false);
glob(bool,dosim,true);
glob(double,startiondist,4.8e-5);
glob(double,startionasym,0.00247);
glob(double,startionshift,7e-7);



glob(int,simsteps,10);
glob(bool,displayiontrap,false);
glob(double,endzeit,0.003);
glob(int,loopanz,100);
int filenum=0;


glob(bool,scatterEvent,true);
glob(bool,useSimSingle,false);
glob(double, reduceAmp, 1.);

glob(int,cyclePts,50);
double cycleArray[50];
double *cycleRadPlot=cycleArray;
double *OmegaRadPlot=cycleArray;


glob(double,freezing,0.99999);
glob(bool, stochKickGlob, true);
glob(bool, doThermalKick, true);
glob(bool, gib_a_aus,false);
glob(bool, verlet, false);
glob(bool,doRMS, false);
glob(double, angle, 11.);
glob(double,Uz,0.16);
glob(double,zOffset,0);  //Center of trap along z axis

glob(double,kicktime,0.);

glob(bool,HeatEngine,false);

double static x0=1.0*0.001;	//Distance rods to trap center
double static z0=4.0*0.001; //Distance endcaps to trap center


#define CPUCNT 8

int ionplotcntmax;
int ionplotcnt[CPUCNT]={};
#define ionplotanz 100000		//100000 geht, bei 1000000 stürzt er ab beim allokieren des pseichers fürs zweite positionsarray
//double tttt[CPUCNT][ionplotanz];
//double yyyy[3][CPUCNT][maxionanz][ionplotanz];

template<class T,int dim> T*makearrayPos(){
  try{
  	return new T[dim];
  }
  catch (exception& e)
  {
	  cout << "Standard exception in file: "<<__FILE__<<" at line "<< __LINE__<<" function: "<<__FUNCTION__<<" problem: "<< e.what() << endl;
  }
}

template<class T,int dim> T*makearrayVel(){
  try{
  	return new T[dim];
  }
  catch (exception& e)
  {
	  cout << "Standard exception in file: "<<__FILE__<<" at line "<< __LINE__<<" function: "<<__FUNCTION__<<" problem: "<< e.what() << endl;
  }
}

template<class T,int dim> T*makearrayTime(){
  try{
  	return new T[dim];
  }
  catch (exception& e)
  {
	  cout << "Standard exception in file: "<<__FILE__<<" at line "<< __LINE__<<" function: "<<__FUNCTION__<<" problem: "<< e.what() << endl;
  }
}

typedef double testtyp[CPUCNT][maxionanz][ionplotanz];
typedef double Tpos[CPUCNT][maxionanz][ionplotanz];
typedef double Tvel[CPUCNT][maxionanz][ionplotanz];
typedef double Ttime[ionplotanz];
Tpos *yyyy=0;//=makearrayPos<Tpos,3>();//--> in function
Tvel *vvvv=0;//=makearrayPos<Tpos,3>();
Ttime *tttt=0;//=makearrayTime<Ttime,CPUCNT>();//--> in function

//delete[] test;


//EXPORT void printptr(){
//  try{
//  	test=new double[3][CPUCNT][maxionanz][ionplotanz];
//  }
//  catch (exception& e)
//  {
//	  cout << "Standard exception in file: "<<__FILE__<<" at line "<< __LINE__<<" function: "<<__FUNCTION__<<" problem: "<< e.what() << endl;
//  }
//}

int dummyLength=10;
double dummyX[10];
double dummyY[10];


//double vvvv[3][CPUCNT][maxionanz][ionplotanz];
double ERad[CPUCNT][maxionanz][ionplotanz];
double ETotal[CPUCNT][maxionanz][ionplotanz];
double ERadMean[maxionanz][ionplotanz];
double ETotMean[maxionanz][ionplotanz];
double yyyyMean[3][maxionanz][ionplotanz];
double HCAmpPlot[ionplotanz];
EXPORT Data1DPlotArray plot1(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot2(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot3(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot4(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot5(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot6(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot7(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot8(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot9(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot10(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot11(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot12(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot13(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot14(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot15(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot16(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot17(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plot18(dummyLength,dummyX, dummyY);

EXPORT Data1DPlotArray plotVoltage(dummyLength,dummyX, dummyY);

EXPORT Data1DPlotArray plotHC(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotMeanX(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotMeanY(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotMeanZ(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotERad(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotETotal(dummyLength,dummyX, dummyY);

EXPORT Data1DPlotArray cyclePlot(cyclePts,OmegaRadPlot,cycleRadPlot);

EXPORT Data1DPlotArray plotX1(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotX2(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotX3(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotX4(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotX5(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotX6(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotX7(dummyLength,dummyX, dummyY);

EXPORT Data1DPlotArray plotY1(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotY2(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotY3(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotY4(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotY5(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotY6(dummyLength,dummyX, dummyY);
EXPORT Data1DPlotArray plotY7(dummyLength,dummyX, dummyY);



double scatterprob[ionplotanz]={};

boost::mutex incmutex;

int incfilenum(){
	boost::mutex::scoped_lock scoped_lock(incmutex);
	int n;
	filenum++;
	n=filenum;
	return n;
}


class semaphore
{
    //The current semaphore count.
    unsigned int count_;

    //mutex_ protects count_.
    //Any code that reads or writes the count_ data must hold a lock on
    //the mutex.

    //Code that increments count_ must notify the condition variable.
    boost::condition_variable condition_;

public:
    boost::mutex mutex_;
    explicit semaphore(unsigned int initial_count) 
       : count_(initial_count),
         mutex_(), 
         condition_()
    {
    }
	void init(int i){count_=i;}

    unsigned int get_count() //for debugging/testing only
    {
        //The "lock" object locks the mutex when it's constructed,
        //and unlocks it when it's destroyed.
        boost::unique_lock<boost::mutex> lock(mutex_);
        return count_;
    }

	void enter() //called "release" in Java
    {
        boost::unique_lock<boost::mutex> lock(mutex_);

        ++count_;

        //Wake up any waiting threads. 
        //Always do this, even if count_ wasn't 0 on entry. 
        //Otherwise, we might not wake up enough waiting threads if we 
        //get a number of signal() calls in a row.
        condition_.notify_all(); 
    }


    void signal() //called "release" in Java
    {
        boost::unique_lock<boost::mutex> lock(mutex_);

        --count_;

        //Wake up any waiting threads. 
        //Always do this, even if count_ wasn't 0 on entry. 
        //Otherwise, we might not wake up enough waiting threads if we 
        //get a number of signal() calls in a row.
        condition_.notify_all(); 
    }

    void wait() //called "acquire" in Java
    {
        boost::unique_lock<boost::mutex> lock(mutex_);
        while (count_ != 0)
        {
             condition_.wait(lock);
        }
    }

};


// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;

base_generator_type generator(std::time(0));




int displaysim=0;


EXPORT void displaysimset(int i){
	if(i<=0) return;
	if(i>CPUCNT) return;

	displaysim=i;
	i--;
	switch(numIons){
			/*case 18:
				plot18.setData(tttt[i],yyyy[2][i][17]);
				plot18.setNum(ionplotcnt[i]);
				plot18.update();
			case 17:
				plot17.setData(tttt[i],yyyy[2][i][16]);
				plot17.setNum(ionplotcnt[i]);
				plot17.update();
			case 16:
				plot16.setData(tttt[i],yyyy[2][i][15]);
				plot16.setNum(ionplotcnt[i]);
				plot16.update();
			case 15:
				plot15.setData(tttt[i],yyyy[2][i][14]);
				plot15.setNum(ionplotcnt[i]);
				plot15.update();
			case 14:
				plot14.setData(tttt[i],yyyy[2][i][13]);
				plot14.setNum(ionplotcnt[i]);
				plot14.update();
			case 13:
				plot13.setData(tttt[i],yyyy[2][i][12]);
				plot13.setNum(ionplotcnt[i]);
				plot13.update();
			case 12:
				plot12.setData(tttt[i],yyyy[2][i][11]);
				plot12.setNum(ionplotcnt[i]);
				plot12.update();
			case 11:
				plot11.setData(tttt[i],yyyy[2][i][10]);
				plot11.setNum(ionplotcnt[i]);
				plot11.update();
			case 10:
				plot10.setData(tttt[i],yyyy[2][i][9]);
				plot10.setNum(ionplotcnt[i]);
				plot10.update();
			case 9:
				plot9.setData(tttt[i],yyyy[2][i][8]);
				plot9.setNum(ionplotcnt[i]);
				plot9.update();
			case 8:
				plot8.setData(tttt[i],yyyy[2][i][7]);
				plot8.setNum(ionplotcnt[i]);
				plot8.update();
			case 7:
				plot7.setData(tttt[i],yyyy[2][i][6]);
				plot7.setNum(ionplotcnt[i]);
				plot7.update();
			case 6:
				plot6.setData(tttt[i],yyyy[2][i][5]);
				plot6.setNum(ionplotcnt[i]);
				plot6.update();
			case 5:
				plot5.setData(tttt[i],yyyy[2][i][4]);
				plot5.setNum(ionplotcnt[i]);
				plot5.update();
			case 4:
				plot4.setData(tttt[i],yyyy[2][i][3]);
				plot4.setNum(ionplotcnt[i]);
				plot4.update();*/
			case 3:
				plot3.setData(tttt[i],yyyy[2][i][2]);
				plot3.setNum(ionplotcnt[i]);
				plot3.update();
			case 2:
				plot2.setData(tttt[i],yyyy[2][i][1]);
				plot2.setNum(ionplotcnt[i]);
				plot2.update();
			case 1:
				plot1.setData(tttt[i],yyyy[2][i][0]);
				plot1.setNum(ionplotcnt[i]);
				plot1.update();
			default: break;
			}

	switch(numIons){
			/*case 7:
				plotX7.setData(tttt[i],yyyy[0][i][6]);
				plotX7.setNum(ionplotcnt[0]);
				plotX7.update();
			case 6:
				plotX6.setData(tttt[i],yyyy[0][i][5]);
				plotX6.setNum(ionplotcnt[0]);
				plotX6.update();
			case 5:
				plotX5.setData(tttt[i],yyyy[0][i][4]);
				plotX5.setNum(ionplotcnt[0]);
				plotX5.update();
			case 4:
				plotX4.setData(tttt[i],yyyy[0][i][3]);
				plotX4.setNum(ionplotcnt[0]);
				plotX4.update();*/
			case 3:
				plotX3.setData(tttt[i],yyyy[0][i][2]);
				plotX3.setNum(ionplotcnt[0]);
				plotX3.update();
			case 2:
				plotX2.setData(tttt[i],yyyy[0][i][1]);
				plotX2.setNum(ionplotcnt[0]);
				plotX2.update();
			case 1:
				plotX1.setData(tttt[i],yyyy[0][i][0]);
				plotX1.setNum(ionplotcnt[0]);
				plotX1.update();
			default: break;
			}

	switch(numIons){
			/*case 7:
				plotY7.setData(tttt[i],yyyy[1][i][6]);
				plotY7.setNum(ionplotcnt[0]);
				plotY7.update();
			case 6:
				plotY6.setData(tttt[i],yyyy[1][i][5]);
				plotY6.setNum(ionplotcnt[0]);
				plotY6.update();
			case 5:
				plotY5.setData(tttt[i],yyyy[1][i][4]);
				plotY5.setNum(ionplotcnt[0]);
				plotY5.update();
			case 4:
				plotY4.setData(tttt[i],yyyy[1][i][3]);
				plotY4.setNum(ionplotcnt[0]);
				plotY4.update();*/
			case 3:
				plotY3.setData(tttt[i],yyyy[1][i][2]);
				plotY3.setNum(ionplotcnt[0]);
				plotY3.update();
			case 2:
				plotY2.setData(tttt[i],yyyy[1][i][1]);
				plotY2.setNum(ionplotcnt[0]);
				plotY2.update();
			case 1:
				plotY1.setData(tttt[i],yyyy[1][i][0]);
				plotY1.setNum(ionplotcnt[0]);
				plotY1.update();
			default: break;
			}
		plotHC.setData(tttt[i],HCAmpPlot);
		plotHC.setNum(ionplotcnt[0]);
		plotHC.update();

		
	//if(displayiontrap) IonTrapDisplay.update();
}




glob(bool,dotrapfreq,false);
glob(double,t0glob,0.002);


EXPORT void clearPlot(){

	if(yyyy==0){
		yyyy=makearrayPos<Tpos,3>();cout<<"yyyy"<<endl;
		vvvv=makearrayVel<Tvel,3>();cout<<"vvvv"<<endl;
		tttt=makearrayTime<Ttime,CPUCNT>();cout<<"tttt"<<endl;
	}

	for(int i=0;i<maxionanz;i++){
		for(int n=0;n<CPUCNT;n++){
			ionplotcnt[n]=0;
			for(int j=0;j<ionplotanz;j++){
				for (int dim=0;dim<3;dim++){
					yyyy[dim][n][i][j]=0;
					vvvv[dim][n][i][j]=0;
					yyyyMean[dim][i][j]=0;
					tttt[n][j]=0;
					HCAmpPlot[j]=0;
				}
			}
		}
	}
	
	switch(maxionanz){
			case 18:
				plot18.setData(tttt[0],yyyy[2][0][17]);
				plot18.setNum(ionplotcnt[0]);
				plot18.update();
			case 17:
				plot17.setData(tttt[0],yyyy[2][0][16]);
				plot17.setNum(ionplotcnt[0]);
				plot17.update();
			case 16:
				plot16.setData(tttt[0],yyyy[2][0][15]);
				plot16.setNum(ionplotcnt[0]);
				plot16.update();
			case 15:
				plot15.setData(tttt[0],yyyy[2][0][14]);
				plot15.setNum(ionplotcnt[0]);
				plot15.update();
			case 14:
				plot14.setData(tttt[0],yyyy[2][0][13]);
				plot14.setNum(ionplotcnt[0]);
				plot14.update();
			case 13:
				plot13.setData(tttt[0],yyyy[2][0][12]);
				plot13.setNum(ionplotcnt[0]);
				plot13.update();
			case 12:
				plot12.setData(tttt[0],yyyy[2][0][11]);
				plot12.setNum(ionplotcnt[0]);
				plot12.update();
			case 11:
				plot11.setData(tttt[0],yyyy[2][0][10]);
				plot11.setNum(ionplotcnt[0]);
				plot11.update();
			case 10:
				plot10.setData(tttt[0],yyyy[2][0][9]);
				plot10.setNum(ionplotcnt[0]);
				plot10.update();
			case 9:
				plot9.setData(tttt[0],yyyy[2][0][8]);
				plot9.setNum(ionplotcnt[0]);
				plot9.update();
			case 8:
				plot8.setData(tttt[0],yyyy[2][0][7]);
				plot8.setNum(ionplotcnt[0]);
				plot8.update();
			case 7:
				plot7.setData(tttt[0],yyyy[2][0][6]);
				plot7.setNum(ionplotcnt[0]);
				plot7.update();
			case 6:
				plot6.setData(tttt[0],yyyy[2][0][5]);
				plot6.setNum(ionplotcnt[0]);
				plot6.update();
			case 5:
				plot5.setData(tttt[0],yyyy[2][0][4]);
				plot5.setNum(ionplotcnt[0]);
				plot5.update();
			case 4:
				plot4.setData(tttt[0],yyyy[2][0][3]);
				plot4.setNum(ionplotcnt[0]);
				plot4.update();
			case 3:
				plot3.setData(tttt[0],yyyy[2][0][2]);
				plot3.setNum(ionplotcnt[0]);
				plot3.update();
			case 2:
				plot2.setData(tttt[0],yyyy[2][0][1]);
				plot2.setNum(ionplotcnt[0]);
				plot2.update();
			case 1:
				plot1.setData(tttt[0],yyyy[2][0][0]);
				plot1.setNum(ionplotcnt[0]);
				plot1.update();
			default: break;
			}

	switch(7){
			case 7:
				plotX7.setData(tttt[0],yyyy[0][0][6]);
				plotX7.setNum(ionplotcnt[0]);
				plotX7.update();
			case 6:
				plotX6.setData(tttt[0],yyyy[0][0][5]);
				plotX6.setNum(ionplotcnt[0]);
				plotX6.update();
			case 5:
				plotX5.setData(tttt[0],yyyy[0][0][4]);
				plotX5.setNum(ionplotcnt[0]);
				plotX5.update();
			case 4:
				plotX4.setData(tttt[0],yyyy[0][0][3]);
				plotX4.setNum(ionplotcnt[0]);
				plotX4.update();
			case 3:
				plotX3.setData(tttt[0],yyyy[0][0][2]);
				plotX3.setNum(ionplotcnt[0]);
				plotX3.update();
			case 2:
				plotX2.setData(tttt[0],yyyy[0][0][1]);
				plotX2.setNum(ionplotcnt[0]);
				plotX2.update();
			case 1:
				plotX1.setData(tttt[0],yyyy[0][0][0]);
				plotX1.setNum(ionplotcnt[0]);
				plotX1.update();
			default: break;
			}

	switch(7){
			case 7:
				plotY7.setData(tttt[0],yyyy[1][0][6]);
				plotY7.setNum(ionplotcnt[0]);
				plotY7.update();
			case 6:
				plotY6.setData(tttt[0],yyyy[1][0][5]);
				plotY6.setNum(ionplotcnt[0]);
				plotY6.update();
			case 5:
				plotY5.setData(tttt[0],yyyy[1][0][4]);
				plotY5.setNum(ionplotcnt[0]);
				plotY5.update();
			case 4:
				plotY4.setData(tttt[0],yyyy[1][0][3]);
				plotY4.setNum(ionplotcnt[0]);
				plotY4.update();
			case 3:
				plotY3.setData(tttt[0],yyyy[1][0][2]);
				plotY3.setNum(ionplotcnt[0]);
				plotY3.update();
			case 2:
				plotY2.setData(tttt[0],yyyy[1][0][1]);
				plotY2.setNum(ionplotcnt[0]);
				plotY2.update();
			case 1:
				plotY1.setData(tttt[0],yyyy[1][0][0]);
				plotY1.setNum(ionplotcnt[0]);
				plotY1.update();
			default: break;
			}

		plotHC.setData(tttt[0],HCAmpPlot);
		plotHC.setNum(ionplotcnt[0]);
		plotHC.update();

	//if(displayiontrap) IonTrapDisplay.update();
}


bool doTestSqueeze=false;
bool TestSqueezeDone[CPUCNT];
bool phaseScrambleDone[CPUCNT];
glob(double, squeezing, 0.);
glob(bool, doSqueeze, false);

double trand[CPUCNT];
bool restarted[CPUCNT];
bool HotLaserOn[CPUCNT];
bool ColdLaserOn[CPUCNT];

int noiseFunction;
double noiseParam[3];
double noiseRF[CPUCNT];


Kibblesim::Kibblesim():simnum(simnum){
	static int simcnt=0;
	simcnt++;
	simnum=simcnt;

	//initSim();	//Das darf hier nicht stehen, sonst schmiewrt alles ab. Konstruktor wird vor der GUI ausgeführt. 
					//Daher werden alle Globs mit ihren startwerten und nicht mit den Werten der Gui gesetzt. 
					//Dies kann zu Problemen führen!!!
}
Kibblesim::~Kibblesim(){}

glob(double, OmegaRes, 2914262.4);
glob(double, epsilon, 1.);
glob(bool,doStepSqueeze,false);

double rfMod[CPUCNT];

int Kibblesim::Signum(double zahl){
	if(zahl>0.)return 1;
	if(zahl<0.)return -1;
	if(zahl==0.)return 0;
}

int Signum2(double zahl){
	if(zahl>0.)return 1;
	if(zahl<0.)return -1;
	if(zahl==0.)return 0;
}

double Kibblesim::RandomReal(double min, double max){
  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::uniform_real<> uni_dist(min,max);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

  std::cout.setf(std::ios::fixed); 
  // calling the generator as a zero-argument function.

return uni();
}

double Kibblesim::RandomGauss(double mean, double sigma){
  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::normal_distribution<> norm_dist(mean,sigma);
  boost::variate_generator<base_generator_type&, boost::normal_distribution<> > norma(generator, norm_dist);

  return norma();
}

double Kibblesim::RandomExp(double tau){
  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::exponential_distribution<> exp_dist(tau);
  boost::variate_generator<base_generator_type&, boost::exponential_distribution<> > expo(generator, exp_dist);

  return expo();
}


glob(bool,phaseScramble,false);

void Kibblesim::PhaseScramble(double v[maxionanz][3], double x[maxionanz][3], double a[maxionanz][3]){
	
	if(phaseScrambleDone[simnum-1]==false && ColdLaserOn[simnum-1]==false){
		phaseScrambleDone[simnum-1]=true;
		double EScr[maxionanz][3]={};
		for(int i=0;i<numIons;i++){
		
			EScr[i][0]=(x[i][0]*x[i][0])*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)) + (v[i][0]*v[i][0]);
			EScr[i][1]=(x[i][1]*x[i][1])*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)) + (v[i][1]*v[i][1]);
			cout<<EScr[i][0]<<"\t"<<x[i][0]*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))<<"\t"<<v[i][0]<<"\t/\t";
			v[i][0]=RandomReal(-sqrt(EScr[i][0]),sqrt(EScr[i][0]));
			x[i][0]=Signum(RandomReal(-1,1))*sqrt(EScr[i][0]-(v[i][0]*v[i][0]))/(2*pi*(x[i][2]*(-2241983470.)+2969444.));
			v[i][1]=RandomReal(-sqrt(EScr[i][1]),sqrt(EScr[i][1]));
			x[i][1]=Signum(RandomReal(-1,1))*sqrt(EScr[i][1]-(v[i][1]*v[i][1]))/(2*pi*(x[i][2]*(-2241983470.)+2969444.));

			//randX=RandomReal(0,pi/2);
			//randY=RandomReal(0,pi/2);
			//v[i][0]=Signum(RandomReal(-1,1))*sqrt(EScr[i][0]*sin(randX));
			//v[i][1]=Signum(RandomReal(-1,1))*sqrt(EScr[i][1]*sin(randY));
			//x[i][0]=Signum(RandomReal(-1,1))*sqrt(EScr[i][0]/((2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)))*cos(randX));
			//x[i][1]=Signum(RandomReal(-1,1))*sqrt(EScr[i][1]/((2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)))*cos(randY));
			a[i][0]=0;
			a[i][1]=0;
			a[i][2]=0;
			//cout<<x[i][0]*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))<<"\t"<<v[i][0]<<endl;
			cout<<(x[i][0]*x[i][0])*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)) + (v[i][0]*v[i][0])<<endl;
		}
	}
}

//double Voltage[ionplotanz];

void Kibblesim::StepSqueeze(double v[maxionanz][3],double v12[maxionanz][3],double x[maxionanz][3],double t,double h){
	
	if((t>endzeit/2.) && (t<endzeit*2./3.)){
		if( Signum(sin(2.*2.*pi*OmegaRes*t)) < Signum(sin(2.*2.*pi*OmegaRes*(t-h))) ){
			rfMod[simnum-1]=epsilon;
			//cout<<"bum"<<endl;
			/*if(simnum-1==0)Voltage[*/
		}
		if( Signum(sin(2.*2.*pi*OmegaRes*t)) > Signum(sin(2.*2.*pi*OmegaRes*(t-h))) ){
			rfMod[simnum-1]=1;
		}
		//if( Signum(sin(OmegaRes*t)) != Signum(sin(OmegaRes*(t-h))) ){	//macht 2 Omega!!!!
		//	for(int ion=0;ion<numIons;ion++){
		//		for(int dim=0;dim<3;dim++){
		//			v[ion][dim]=epsilon*v[ion][dim];
		//			v12[ion][dim]=epsilon*v12[ion][dim];
		//			x[ion][dim]=x[ion][dim]/epsilon;
		//		}
		//	}
		//}
	
	}
	

}



void Kibblesim::TestSqueezeKick(double v[maxionanz][3],double v12[maxionanz][3],double x[maxionanz][3],double &t,double h){
	if(TestSqueezeDone[simnum-1]==false && t>endzeit/3.){
		//cout<<"squeeze"<<endl;
		for(int ion=0;ion<numIons;ion++){
			for(int dim=0;dim<3;dim++){
				v[ion][dim]=v[ion][dim]*squeezing;
				v12[ion][dim]=v12[ion][dim]*squeezing;
				x[ion][dim]=x[ion][dim]/squeezing;
			}
		}
		TestSqueezeDone[simnum-1]=true;
		
	}
	if(t>endzeit/2. && t<endzeit*4./5.){
		Laser(v12,h,0);
		Laser(v12,h,1);
		Laser(v12,h,3);
		Laser(v12,h,4);
	}
	/*if(t<endzeit/4){
		Laser(v12,h,0);
		Laser(v12,h,1);
		Laser(v12,h,3);
		Laser(v12,h,4);
	}*/

}




//glob(double,meanScatteringRate,0);

glob(bool, LaserX1, false);
glob(bool, LaserX2, false);
glob(bool, LaserY1, false);
glob(bool, LaserY2, false);
glob(bool, LaserZ1, false);
glob(bool, LaserZ2, false);

glob(int, PhotonCount,0);

EXPORT void FlipLaser(){
	bool steuerBool=false;
	steuerBool=LaserX1;
	if(LaserX2!=steuerBool || LaserY1!=steuerBool || LaserY2!=steuerBool || LaserZ1!=steuerBool || LaserZ2!=steuerBool){
		LaserX1Set(false);
		LaserX2Set(false);
		LaserY1Set(false);
		LaserY2Set(false);
		LaserZ1Set(false);
		LaserZ2Set(false);
	}
	else{
		LaserX1Set(!LaserX1);
		LaserX2Set(!LaserX2);
		LaserY1Set(!LaserY1);
		LaserY2Set(!LaserY2);
		LaserZ1Set(!LaserZ1);
		LaserZ2Set(!LaserZ2);
	}
}

double satMT[CPUCNT];
double detuneMT[CPUCNT];
double HCAmp[CPUCNT];


glob(double, reduceZZ, 0.9);
glob(bool, Casdorff, false);
glob(bool, MOT_mode, false);

void Kibblesim::Laser(double v[maxionanz][3], double deltat,int axis){//Formeln aus Apl. Phys. B 45, 175
	int r;
	double k[3];					
	double TFreq=C/(396.9592e-9);
	if(MOT_mode)TFreq=C/(780.0e-9);
	double LFreq=TFreq+detuneMT[simnum-1];
	double kBetrag=2.*pi*LFreq/C;
	double kResonanzBetrag=2.*pi*TFreq/C;
	double kL=2.*pi*LFreq/C;			//Wellenvektor Laser(verstimmt)
	double kI=2.*pi*TFreq/C;			//Wellenvektor Übergang Ion
	double gammasqr=gamma*gamma;
	//double temp[maxionanz];
	//double P[maxionanz];
	//double Psp[maxionanz];
	
	switch(axis){
		case 1:
			k[0]=-kBetrag;
			k[1]=0.;
			k[2]=0.;
			break;
		case 4:
			k[0]=kBetrag;
			k[1]=0.;
			k[2]=0.;
			break;
		case 0:
			k[0]=0.;
			k[1]=kBetrag;
			k[2]=0.;
			break;
		case 3:
			k[0]=0.;
			k[1]=-kBetrag;
			k[2]=0.;
			break;
		case 2:
			k[0]=0.;
			k[1]=0.;
			k[2]=kBetrag;
			break;
		case 5:
			k[0]=0.;
			k[1]=0.;
			k[2]=-kBetrag;
			break;
		default: break;
	}

	double skalarprod[maxionanz];
	for(int i=0; i< numIons; i++){
		skalarprod[i]=k[0]*v[i][0]+k[1]*v[i][1]+k[2]*v[i][2];
	}
	if(Casdorff){
		double temp=detuneMT[simnum-1]-skalarprod[0];
		temp*=temp;
		double P=satMT[simnum-1]*gammasqr/(temp+gammasqr);
		double Psp=gamma*deltat*P/(1.+2.*P);
		if(RandomReal(0.,1.)<=Psp){//photonabsorbed
			double u[3];
			double uBetrag;
			for(r=0;r<3;r++){
				u[r]=RandomReal(-1.,1.);
			}
			uBetrag=sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

			for(r=0;r<3;r++){
				u[r]/=uBetrag;
			}
			
			for(r=0;r<3;r++){
				v[0][r]+=hbar/MCa*k[r]  +  hbar*kResonanzBetrag/MCa*u[r];	
			}
	
		}

	}
	else{
		for(int ind=0; ind <numIons; ind++){
			double term=2*(detuneMT[simnum-1]-skalarprod[ind])/gamma;
			double scatteringrate=(satMT[simnum-1]*gamma/2.)/(1.+satMT[simnum-1]+term*term);//http://info.phys.unm.edu/~ideutsch/Classes/Phys500S09/Downloads/handpubl.pdf
			
			if(deltat*scatteringrate>1.) cout <<"timestep too big for laserinteraction" <<endl;

			if(RandomReal(0.,1.)<=(deltat*scatteringrate)){//photonabsorbed
				//lebenszeit[ind]=random.Exp(7e-9);
				scatterEvent=true;
				if(useSimSingle)PhotonCountSet(PhotonCount+1);
				double u[3];
				double thetaL=RandomReal(0,pi);
				double phiL=RandomReal(0,2*pi);
				//u[2]=cos(theta);
				u[2]=cos(thetaL)*reduceZZ;
				u[0]=sin(thetaL)*sin(phiL);
				u[1]=sin(thetaL)*cos(phiL);
				
				for(r=0;r<3;r++){
					v[ind][r]+= hbar*k[r]/MCa   +  hbar*kResonanzBetrag*u[r]/MCa;	
				}
			}
			
		}
	}
}



void Kibblesim::LaserXYZ(double v[maxionanz][3], double deltat){//Formeln aus Apl. Phys. B 45, 175
	int r;
	double k[3];				
	//double P,Psp;		
	double TFreq=C/(396.9592e-9);
	double LFreq=TFreq+detuneMT[simnum-1];
	double kBetrag=2.*pi*LFreq/C;
	double kResonanzBetrag=2.*pi*TFreq/C;
	double kL=2.*pi*LFreq/C;			//Wellenvektor Laser(verstimmt)
	double kI=2.*pi*TFreq/C;			//Wellenvektor Übergang Ion
	k[0]=kBetrag/sqrt(3.);
	k[1]=kBetrag/sqrt(3.);
	k[2]=kBetrag/sqrt(3.);

	double skalarprod[maxionanz];
	for(int i=0; i< numIons; i++){
		skalarprod[i]=k[0]*v[i][0]+k[1]*v[i][1]+k[2]*v[i][2];
	}

	
	double gammasqr=gamma*gamma;
	//double temp[maxionanz];
	//double P[maxionanz];
	//double Psp[maxionanz];
	for(int ind=0; ind <numIons; ind++){
	//ind=7;
		double term=2*(detuneMT[simnum-1]-skalarprod[ind])/gamma;
		double scatteringrate=(satMT[simnum-1]*gamma/2.)/(1.+satMT[simnum-1]+term*term);//http://info.phys.unm.edu/~ideutsch/Classes/Phys500S09/Downloads/handpubl.pdf
		
		//temp[ind]=detune-skalarprod[ind];
		//temp[ind]*=temp[ind];
		//P[ind]=sat*gammasqr/(temp[ind]+gammasqr);
		//Psp[ind]=gamma*deltat*P[ind]/(1.+2.*P[ind]);
		//if(lebenszeit[ind]>0) lebenszeit[ind]-=deltat;
		//if((lebenszeit[ind]<=0) && (random.Uniform(0.,1.)<=Psp[ind])){//photonabsorbed
		if(deltat*scatteringrate>1.) cout <<"timestep too big for laserinteraction" <<endl;
		
		
		//meanScatteringRate+=deltat*scatteringrate;

		
		if(RandomReal(0.,1.)<=deltat*scatteringrate){//photonabsorbed
			//lebenszeit[ind]=random.Exp(7e-9);
			scatterEventSet(true);
			double u[3];
			double thetaL=RandomReal(0,pi);
			double phiL=RandomReal(0,2*pi);
			u[2]=cos(thetaL);
			u[0]=sin(thetaL)*cos(phiL);
			u[1]=sin(thetaL)*sin(phiL);
			
			for(r=0;r<3;r++){
				//k[r]=0;
				v[ind][r]+= hbar*k[r]/(1.*MCa)   +  hbar*kResonanzBetrag*u[r]/(1.*MCa);	
				
			}
		}
		
	}
	
}




//EXPORT void getFittedTaborValues(){
//	double taborOffsettemp=((fittedoffset/2.04)+0.5*(fittedamp/2.04));
//	double taborAmptemp=2*0.5*(fittedamp/2.04);
//	taborAmpSet(taborAmptemp);
//	taborOffsetSet(taborOffsettemp);
//}
	
glob(bool,doScattering,false);
glob(double,tScattering,1e-5);



double Kibblesim::scatterCount(){
	for(int i=0;i<ionplotanz;i++){
		scatterprob[i]=0;
	}
	
	int scatterCountNum=0;
	int axis=2;
	int simnum=0;
	
	doScatteringSet(true);
	
	initSim();
	ionplotcnt[0]=0;
	int oldNumIons=numIons;
	numIons=1;
	int oldsimsteps=simsteps;
	

	
	x[0][0]=0;//RandomReal(-1e6,1.e6);
	x[0][1]=0;//random.Uniform(-1e6,1.e6);
	x[0][2]=1e-6;
	v[0][0]=0;//random.Uniform(-10,10);;
	v[0][1]=0;//random.Uniform(-10,10);;
	v[0][2]=0;//random.Uniform(-10,10);;
	a[0][0]=0;
	a[0][1]=0;
	a[0][2]=0;
	t=0;
	//x[0][axis]=1e-7;//1e-7//RandomReal(-1e4,1e4);
	double valuebefore=1;
	double tstart=-1;
	int oscicnt=0;
	while(doScattering){
		propagateForwardVerlet(t,simsteps,x,v);
		simnum++;
		if(scatterEvent){ 
			scatterCountNum++;
			scatterEvent=false;
		}

		
		double k[3];						
		double TFreq=C/(396.9592e-9);
		double LFreq=TFreq+detune;
		double kBetrag=2.*pi*LFreq/C;
		double kResonanzBetrag=2.*pi*TFreq/C;
		double kL=2.*pi*LFreq/C;			//Wellenvektor Laser(verstimmt)
		double kI=2.*pi*TFreq/C;			//Wellenvektor Übergang Ion
		k[0]=kBetrag/sqrt(3.);
		k[1]=kBetrag/sqrt(3.);
		k[2]=kBetrag/sqrt(3.);
		double skalarprod[maxionanz];
		skalarprod[0]=k[0]*v[0][0]+k[1]*v[0][1]+k[2]*v[0][2];
		double gammasqr=gamma*gamma;
		//double temp[maxionanz];
		//double P[maxionanz];
		//double Psp[maxionanz];
		double term=2*(detune-skalarprod[0])/gamma;
		double scatteringrate=(sat*gamma/2.)/(1.+sat+term*term);
		scatterprob[ionplotcnt[0]]=scatteringrate*h;
		
		tttt[0][ionplotcnt[0]]=t;
		yyyy[2][0][0][ionplotcnt[0]]=x[0][axis];
		
		
		if(t>tScattering){
			doScatteringSet(false);
			break;
		}
		
		valuebefore=yyyy[2][0][0][ionplotcnt[0]];
		ionplotcnt[0]++;
		plot1.setNum(ionplotcnt[0]);
		plot1.update();
		if(displayiontrap) IonTrapDisplay.update();
		if(ionplotcnt[0]>ionplotanz) ionplotcnt[0]=0;
	}
	cout<<simnum<<endl;
	char* filename1 = (char*) calloc(300,sizeof(char));
    sprintf(filename1,"C:\\data\\ergebnis\\ScatterProb.tsv");
    ofstream fname1(filename1);
	for(int i=0;i<simnum;i++){
		fname1<<scatterprob[i]<<endl;
	}
    fname1.close();


	numIons=oldNumIons;
	simsteps=oldsimsteps;
	return double(scatterCountNum)/t/1.0e6;

}



bool dampfinished[CPUCNT];




//for(int i=0;i<CPUCNT;i++)dampfinished[i]=false;


glob(double,phaseTest,0);

void Kibblesim::Force(double t,double x[maxionanz][3], double v[maxionanz][3],double a[maxionanz][3], double v12[maxionanz][3],double phi){
	double tempT=temperature;
	double sigma=sqrt(k_B*temperature/MCa);
	int zahl=0;
	double u[3]={};
	double sigma2=2*k_B*tempT/MCa;
	double alpha=(angle/360.)*(2.*PI);
	double tanangle=tan(alpha);
	double accPot[maxionanz]={};
	double accPotX[maxionanz]={};
	
	double accAx=(Uz+0.5*(rfvoltage-rfvoltageX))/(z0*z0);//Edited to include axial micromotion. Vanishes for asymmetric drive

	//double eta=0;
	
	t0=t0glob;

	if(MOT_mode)rfvoltage=rfamp+noiseRF[simnum-1];
	else{ 
		rfvoltage=(rfamp*rfMod[simnum-1])*cos(2.*pi*OmegaRF*t+phi)+noiseRF[simnum-1];		//rfamp ist glob
		rfvoltageX=(rfamp*rfMod[simnum-1])*cos(2.*pi*OmegaRF*t+phi+phaseTest)+noiseRF[simnum-1];
	}
	// rfMod ist zum squeezen

	//trapping potential
	for(int m=0; m < numIons; m++){
		accPot[m]=rfvoltage/((x0+x[m][2]*tanangle)*(x0+x[m][2]*tanangle));		//EINHEITEN CHECKEN!!!!!!!!!!!!!!!!!!
		accPotX[m]=rfvoltageX/((x0+x[m][2]*tanangle)*(x0+x[m][2]*tanangle));
		
		if(MOT_mode){	
			a[m][X] = -2. * qDivM * accPot[m] * x[m][X];
			a[m][Y] = -2. * qDivM * accPot[m] * x[m][Y];
			a[m][Z] =  qDivM * (2. * accPot[m]/(x0+x[m][Z]*tanangle)*tanangle*(x[m][X]*x[m][X]+x[m][Y]*x[m][Y]) - 2. * accAx * (x[m][Z]-zOffset));
		}
		else{
			a[m][X] = -2. * qDivM * accPotX[m] * x[m][X];
			a[m][Y] =  2. * qDivM * accPot[m] * x[m][Y];
			a[m][Z] =  qDivM * (2. * accPot[m]/(x0+x[m][Z]*tanangle)*tanangle*(x[m][X]*x[m][X]-x[m][Y]*x[m][Y]) - 2. * accAx * (x[m][Z]-zOffset));
		}
	
		//F[2]=  QoverM *(2.*Urf*cos(2.*pi*OmegaRF*t+phi)/pow((x0+x[2][j]*tanalpha),3)*tanalpha* (pow(x[0][j],2)-pow(x[1][j],2)) - 2.*b*(x[2][j]-z0));
	}

	if(numIons!=1){
		//Coulomb interaction
		for(int m=1; m < numIons; m++){	
			for(int k=0; k < m; k++){
				double accel[3]={};
				double r=0;
				for(int dim=0; dim < 3; dim++){
					//force on charge m due to charge k
					double dist=(x[m][dim]-x[k][dim]);
					r+=dist*dist;
					accel[dim]=forceconstant*dist;
				}
				r=sqrt(r);
				r=r*r*r;
				for(int dim=0; dim < 3; dim++) {
					accel[dim]/=r;
					a[m][dim]+=accel[dim];
					a[k][dim]-=accel[dim];
				}
			}
		}
	}

	////stochastic kicks mimicking laser interaction
	//if(verlet && (t>freezingTime) && (stochKickGlob)){
	//	for(int ind=0;ind<numIons;ind++){			
	//		u[2]=RandomGauss(0,sigma);
	//		u[0]=RandomGauss(0,sigma);
	//		u[1]=RandomGauss(0,sigma);
	//		//cout << "u[2]= "<<u[2] <<" u[0]= " << u[0] << " u[1]= " << u[1] <<endl;
	//		//cout << "sigma= " << sigma << endl;
	//		for(int dim=0;dim<3;dim++){
	//			//xStart[ind][dim]=x[ind][dim];
	//			v[ind][dim]+=u[dim];				//ich glaub die muessen auf mm/s umgerechnet werden
	//			v12[ind][dim]+=u[dim];
	//			//a[ind][dim]=0;
	//		}
	//	}
	//}

	
	////initial damping for equilibriums position
	//if((initialdamp==true) && (t<initialdamptime) && (dotrapfreq==false) && (doScattering==false)){
	//	for(int m=0; m < numIons; m++){
	//		for(int dim=0; dim < 3; dim++){			
	//			a[m][dim]=a[m][dim]-v[m][dim]*dissipationfactor;
	//			if(t>(initialdamptime/2.))a[m][dim]=a[m][dim]-v[m][dim]*2.*dissipationfactor;
	//		}
	//	}
	//}

	//thermal kick
	if((t>=kicktime)&&(dampfinished[simnum-1]==false) && doThermalKick){   
		cout<<"KICK"<<endl;
		dampfinished[simnum-1]=true;
		for(int ind=0;ind<numIons;ind++){
			//Würfele nun die Startgeschwindigkeiten
			double u[3];
			/*double theta=random.Uniform(0,pi);
			double phi=random.Uniform(0,2*pi);*/

			/*u[2]=cos(theta)*random.Gaus(0,sigma);
			u[0]=sin(theta)*cos(phi)*random.Gaus(0,sigma);
			u[1]=sin(theta)*sin(phi)*random.Gaus(0,sigma);*/

			
			u[2]=RandomGauss(0,sigma);
			u[0]=RandomGauss(0,sigma);
			u[1]=RandomGauss(0,sigma);
			cout<<"sigma: "<<sigma<<"   z-kick: "<<u[2]<<endl;
			for(int dim=0;dim<3;dim++){
				v[ind][dim]=u[dim];				//ich glaub die muessen auf mm/s umgerechnet werden
				v12[ind][dim]=u[dim];
			}
		}
	}

	//damping
	/*if(!verlet && !HeatEngine){
		for(int m=0;m<numIons;++m){
			for(int dim=0;dim<3;++dim){
				a[m][dim]-=v[m][dim]*eta/MCa;
			}
		}
	}*/

}


void Kibblesim::ForcePseudo(double t,double x[maxionanz][3], double v[maxionanz][3],double a[maxionanz][3], double v12[maxionanz][3],double phi){
	double tempT=temperature;
	double sigma=sqrt(k_B*temperature/MCa);
	int zahl=0;
	double u[3]={};
	double sigma2=2*k_B*tempT/MCa;
	double alpha=(angle/360.)*(2.*PI);
	double tanangle=tan(alpha);


	//pseudo potential
	for(int m=0; m < numIons; m++){
		a[m][X] = -x[m][X]*( 2*pi*(x[m][2]*(-2241983470.)+2969444.))*( 2*pi*(x[m][2]*(-2241983470.)+2969444.));
		a[m][Y] = -x[m][Y]*( 2*pi*(x[m][2]*(-2241983470.)+2969444.))*( 2*pi*(x[m][2]*(-2241983470.)+2969444.));
		a[m][Z] = -x[m][Z]*2*pi*36000*2*pi*36000+ ( 2*pi*(x[m][2]*(-2241983470.)+2969444.))*( 2*pi*(x[m][2]*(-2241983470.)+2969444.))*tanangle*x[m][X] + ( 2*pi*(x[m][2]*(-2241983470.)+2969444.))*( 2*pi*(x[m][2]*(-2241983470.)+2969444.))*tanangle*x[m][Y];
	}

	//Coulomb interaction
	for(int m=1; m < numIons; m++){	
		for(int k=0; k < m; k++){
			double accel[3]={};
			double r=0;
			for(int dim=0; dim < 3; dim++){
				//force on charge m due to charge k
				double dist=(x[m][dim]-x[k][dim]);
				r+=dist*dist;
				accel[dim]=forceconstant*dist;
			}
			r=sqrt(r);
			r=r*r*r;
			for(int dim=0; dim < 3; dim++) {
				accel[dim]/=r;
				a[m][dim]+=accel[dim];
				a[k][dim]-=accel[dim];
			}
		}
	}

	//thermal kick
	if((t>=kicktime)&&(dampfinished[simnum-1]==false) && doThermalKick){   
		//cout<<"KICK"<<endl;
		dampfinished[simnum-1]=true;
		for(int ind=0;ind<numIons;ind++){
			//Würfele nun die Startgeschwindigkeiten
			double u[3];
			u[2]=RandomGauss(0,sigma);
			u[0]=RandomGauss(0,sigma);
			u[1]=RandomGauss(0,sigma);
			for(int dim=0;dim<3;dim++){
				v[ind][dim]=u[dim];				//ich glaub die muessen auf mm/s umgerechnet werden
				v12[ind][dim]=u[dim];
			}
		}
	}
}

glob(bool,doFreeze,false);
glob(bool,verletstoch,false);
glob(bool,RK4_HeatEngine,false);
glob(bool,reduceNoiseZ,false);
glob(bool,RungeKutta4laser,false);
glob(int,rmsion,0);
glob(double,oneminusfreezingforeta,1e-5);


glob(double,HCFreq,35000);
glob(double,HCLimitCold,0.75);
glob(double,HCLimitHot,0.75);
glob(double,StartEngine,0.0001);
glob(double,detuneCold,-1.e7);
glob(double,detuneHot,1.e7);
glob(double,detuneZZ,-1.e7);
glob(double,satZZ,0.1);
glob(bool,coolZ,false);
glob(double,satHEHot,1);
glob(double,satHECold,1);
glob(double,ZZLimit,0.9);
glob(double,coolZZTime,endzeitGet());
glob(double,LaserOnT,0.0);
glob(double,LaserOffT,0.0);
glob(double,TColdBath,0.1);
glob(bool,usePseudoPot,false);
glob(bool,doNoise,false);
glob(double,TimeNoiseStart,0);
glob(double,TimeNoiseStop,0);



void Kibblesim::Restart(double &t){
	if(t>trand[simnum-1] && restarted[simnum-1]==false){
		t=0;
		ionplotcnt[simnum-1]=0;
		restarted[simnum-1]=true;
		//cout<<"restart"<<endl;
	}
}


void Kibblesim::propagateForwardVerlet(double &t,int n,double x[maxionanz][3],double v[maxionanz][3],double *rms, int *iter){
	//int i;
	
	double Ekin=0;
	double eta=MCa*(oneminusfreezingforeta)/h;
	double sigma=sqrt(k_B*temperature/MCa);
	
	double noise_magnitude = sqrt(TColdBath*k_B*(1-freezing)/(MCa*MCa))*sqrt(h);
	
	//variables needed for RK4
	static double kx1[maxionanz][3][CPUCNT];
	static double kx2[maxionanz][3][CPUCNT];
	static double kx3[maxionanz][3][CPUCNT];
	static double kx4[maxionanz][3][CPUCNT];
	static double kx12[maxionanz][3][CPUCNT];
	static double kx22[maxionanz][3][CPUCNT];

	static double kv1[maxionanz][3][CPUCNT];
	static double kv2[maxionanz][3][CPUCNT];
	static double kv3[maxionanz][3][CPUCNT];
	static double kv4[maxionanz][3][CPUCNT];
	static double kv12[maxionanz][3][CPUCNT];
	static double kv22[maxionanz][3][CPUCNT];

	
	//end RK4

	double oldDetune=0;
	double oldSat=0;
	
	
	for(int i=0;i<n;i++){			//n=simsteps
		if(iter) (*iter)++;
		//if(laseron)LaserXYZ(v,h);

		if(!verlet && !HeatEngine && !verletstoch && !RK4_HeatEngine)verletSet(true);

		if(HeatEngine){
			if(verlet) verletSet(false);
			if(RK4_HeatEngine) RK4_HeatEngineSet(false);
			if(verletstoch) verletstochSet(false);

			

			for(int m=0;m<numIons;++m){					//alter Propagator: partitioned Runge Kutta 2. Ordnung
				for(int dim=0;dim<3;++dim){
					v12[m][dim]=v[m][dim]+h/2.*a[m][dim];
					x[m][dim]+=h*v12[m][dim];
				}
			}
			t+=h;
			/*if(usePseudoPot==true)ForcePseudo(t,x,v,a,v12,phi);
			else Force(t,x,v,a,v12,phi);*/
			Force(t,x,v,a,v12,phi);
			HCAmp[simnum-1]=sin(2*pi*HCFreq*t+HCphi);
			
			//if(doTestSqueeze==true && TestSqueezeDone[simnum-1]==false)TestSqueezeKick(v,v12,x,t,h);
			//if(doSqueeze)Squeeze(v,x,h);

			if(t>StartEngine){

				ColdLaserOn[simnum-1]=false;

				if(HCAmp[simnum-1]>HCLimitCold){				
					detuneMT[simnum-1]=detuneCold;
					satMT[simnum-1]=satHECold;
					Laser(v12,h,0);
					Laser(v12,h,1);
					Laser(v12,h,3);
					Laser(v12,h,4);
					if(doFreeze)FreezeRad(t,v12);
					ColdLaserOn[simnum-1]=true;
					phaseScrambleDone[simnum-1]=false;
				}
				
				HotLaserOn[simnum-1]=false;
				
				if(HCAmp[simnum-1]<-HCLimitHot){
					detuneMT[simnum-1]=detuneHot;
					satMT[simnum-1]=satHEHot;
					Laser(v12,h,0);
					Laser(v12,h,1);
					Laser(v12,h,3);
					Laser(v12,h,4);
					//if((sin(2*pi*HCFreq*t+HCphi) < sin(2*pi*HCFreq*(t-h)+HCphi)) && (sin(2*pi*HCFreq*t+HCphi) < sin(2*pi*HCFreq*(t+h)+HCphi)))HotLaserOn[simnum-1]=true;
					HotLaserOn[simnum-1]=true;
					TestSqueezeDone[simnum-1]=false;
				}
				


				if(coolZ && t>coolZZTime){
					if(HCAmp[simnum-1]<-ZZLimit || HCAmp[simnum-1]>ZZLimit){
						detuneMT[simnum-1]=detuneZZ;
						satMT[simnum-1]=satZZ;
						Laser(v12,h,2);
						Laser(v12,h,5);	
						if(doFreeze){
							for(int ind=0; ind<numIons;ind++){
								v12[ind][2]*=reduceZZ;
							}	
						}
					}
				}

				
				
			}
			else if(t>LaserOnT && t<LaserOffT){
				detuneMT[simnum-1]=detune;
				satMT[simnum-1]=sat;
				Laser(v12,h,1);
				Laser(v12,h,4);
				Laser(v12,h,0);
				Laser(v12,h,3);
				Laser(v12,h,5);
				Laser(v12,h,2);

				if(doFreeze)FreezeRad(t,v12);
				if(doFreeze)FreezeAx(t,v12);
			}


			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					v[m][dim]=v12[m][dim]+h/2.*a[m][dim];
					if(ColdLaserOn[simnum-1])v[m][dim]+=noise_magnitude*RandomGauss(0,1);
				}
			}

			
			if(doSqueeze)Squeeze(t,v,x,a,h);
			if(phaseScramble)PhaseScramble(v12,x,a);
			Restart(t);

		}


		else if(verlet){
			RK4_HeatEngine=false;
			verletstoch=false;
			for(int m=0;m<numIons;++m){					//alter Propagator: partitioned Runge Kutta 2. Ordnung
				for(int dim=0;dim<3;++dim){
					v12[m][dim]=v[m][dim]+h/2.*a[m][dim];
					x[m][dim]+=h*v12[m][dim];
				}
			}
			
			t+=h;
			Force(t,x,v,a,v12,phi);

			HCAmp[simnum-1]=sin(2*pi*HCFreq*t+HCphi);

			if(t>LaserOnT && t<LaserOffT){
				detuneMT[simnum-1]=detune;
				satMT[simnum-1]=sat;
				if(laserTrap)LaserXYZ(v12,h);
				if(laseron){
					if(LaserX1)Laser(v12,h,0);
					if(LaserX2)Laser(v12,h,3);
					if(LaserY1)Laser(v12,h,1);
					if(LaserY2)Laser(v12,h,4);
					if(LaserZ1)Laser(v12,h,2);
					if(LaserZ2)Laser(v12,h,5);

					if(doFreeze)FreezeRad(t,v12);
					if(doFreeze)FreezeAx(t,v12);
				}
			}
			
			
			
			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					v[m][dim]=v12[m][dim]+h/2.*a[m][dim];
				}
			}

			if(doTestSqueeze==true && doStepSqueeze==true)StepSqueeze(v,v12,x,t,h);


			if(doTestSqueeze==true && doStepSqueeze==false)TestSqueezeKick(v,v12,x,t,h);
			if(doSqueeze)Squeeze(t,v,x,a,h);

			if(doNoise && t>TimeNoiseStart && t<TimeNoiseStop)NoiseHeat();

			

			//if(doTestSqueeze==true)Restart(t);
			Restart(t);


			/*if(gib_a_aus){
				Ekin=0.5*MCa*v[7][0]*v[7][0];
				cout<<"Ekin= "<<Ekin<<endl;	

			}*/
		}
		
		else if(verletstoch){
			if(RK4_HeatEngine)RK4_HeatEngineSet(false);
			for(int m=0;m<numIons;++m){					//alter Propagator: partitioned Runge Kutta 2. Ordnung
				for(int dim=0;dim<3;++dim){
					v12[m][dim]=v[m][dim]+h/2.*a[m][dim];
					x[m][dim]+=h*v12[m][dim];
				}
			}
			t+=h;
			Force(t,x,v,a,v12,phi);

		
			
			
			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					v[m][dim]=v12[m][dim]+h/2.*a[m][dim];
				}
			}
			

			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					//stochastic kick
					v[m][dim]+=noise_magnitude*RandomGauss(0,1);
				}
			}
			if(gib_a_aus){
				Ekin=0.5*MCa*v[7][0]*v[7][0];
				cout<<"Ekin= "<<Ekin<<endl;	
			}
		}



		else if(RK4_HeatEngine){	//runge kutta 4th order with heat engine
		

			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					kx1[m][dim][simnum-1]=v[m][dim]*h;			//wozu dienen die kx1 und kv1, wenn nur kx12 und kv12 benutzt werden???
					kx12[m][dim][simnum-1]=v[m][dim]*h/2.;
					kv1[m][dim][simnum-1]=a[m][dim]*h;
					kv12[m][dim][simnum-1]=a[m][dim]*h/2.;
					xtemp[m][dim]=x[m][dim]+kx12[m][dim][simnum-1];		//xtemp=x0+v*dt/2
					vtemp[m][dim]=v[m][dim]+kv12[m][dim][simnum-1];
				}
			}
			
			Force(t+h/2.,xtemp,vtemp,a,v12,phi);
			

			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					kx2[m][dim][simnum-1]=vtemp[m][dim]*h;
					kx22[m][dim][simnum-1]=vtemp[m][dim]*h/2.; 
					kv2[m][dim][simnum-1]=a[m][dim]*h;
					kv22[m][dim][simnum-1]=a[m][dim]*h/2.;
					xtemp[m][dim]=x[m][dim]+kx22[m][dim][simnum-1];		//xtemp=x0+vtemp(t/2)*dt/2
					vtemp[m][dim]=v[m][dim]+kv22[m][dim][simnum-1];
				}
			}
			Force(t+h/2.,xtemp,vtemp,a,v12,phi);
			

			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					kx3[m][dim][simnum-1]=vtemp[m][dim]*h;
					kv3[m][dim][simnum-1]=a[m][dim]*h;				//?? kv3=kv2=kv1 (ausser, dass force gewirkt hat. Aber force überschreibt sich selbst
					xtemp[m][dim]=x[m][dim]+kx3[m][dim][simnum-1];
					vtemp[m][dim]=v[m][dim]+kv3[m][dim][simnum-1];
				}
			}
			Force(t+h,xtemp,vtemp,a,v12,phi);
			
			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					kx4[m][dim][simnum-1]=vtemp[m][dim]*h;
					kv4[m][dim][simnum-1]=a[m][dim]*h;
				}
			}

			for(int m=0;m<numIons;++m){
				for(int dim=0;dim<3;++dim){
					x[m][dim]+=1./6.*(kx1[m][dim][simnum-1]+2.*kx2[m][dim][simnum-1]+2.*kx3[m][dim][simnum-1]+kx4[m][dim][simnum-1]);
					v[m][dim]+=1./6.*(kv1[m][dim][simnum-1]+2.*kv2[m][dim][simnum-1]+2.*kv3[m][dim][simnum-1]+kv4[m][dim][simnum-1]);
					//stochastic kick
					/*if(!laseron){
						if(reduceNoiseZ && dim==2) v[m][dim]+=noise_magnitude/100*RandomGauss(0,1);
						else v[m][dim]+=noise_magnitude*RandomGauss(0,1);
					}*/
				}
			}

			t+=h;
			
			Force(t,x,v,a,v12,phi);
			
			
			HCAmp[simnum-1]=sin(2*pi*HCFreq*t+HCphi);
			
			if(t>StartEngine){
				if(HCAmp[simnum-1]>HCLimitCold){				
					detuneMT[simnum-1]=detuneCold;
					satMT[simnum-1]=satHECold;
					Laser(v,h,0);
					Laser(v,h,1);
					Laser(v,h,3);
					Laser(v,h,4);
					if(doFreeze==true)Freeze(t,v);				
				}
				
				HotLaserOn[simnum-1]=false;
				

				if(HCAmp[simnum-1]<-HCLimitHot){
					detuneMT[simnum-1]=detuneHot;
					satMT[simnum-1]=satHEHot;
					Laser(v,h,0);
					Laser(v,h,1);
					Laser(v,h,3);
					Laser(v,h,4);
					//if((sin(2*pi*HCFreq*t+HCphi) < sin(2*pi*HCFreq*(t-h)+HCphi)) && (sin(2*pi*HCFreq*t+HCphi) < sin(2*pi*HCFreq*(t+h)+HCphi)))HotLaserOn[simnum-1]=true;
					HotLaserOn[simnum-1]=true;
					TestSqueezeDone[simnum-1]=false;
				}
				


				if(coolZ && t>coolZZTime){
					if(HCAmp[simnum-1]<-ZZLimit || HCAmp[simnum-1]>ZZLimit){
						detuneMT[simnum-1]=detuneZZ;
						satMT[simnum-1]=satZZ;
						Laser(v,h,2);
						Laser(v,h,5);						
					}
				}

				
				
			}


			else if(t>LaserOnT && t<LaserOffT){
				detuneMT[simnum-1]=detune;
				satMT[simnum-1]=sat;
				if(LaserX1)Laser(v,h,0);
				if(LaserX2)Laser(v,h,3);
				if(LaserY1)Laser(v,h,1);
				if(LaserY2)Laser(v,h,4);
				if(LaserZ1)Laser(v,h,2);
				if(LaserZ2)Laser(v,h,5);
			}

			if(doSqueeze)Squeeze(t,v,x,a,h);
			Restart(t);

		}
			

			
		



		/*if(useSimSingle==true && ((int)(t/h)%simsteps)==0){
			int zahl=(t/h)/simsteps;
			for(int i=0;i<numIons;i++){
				for(int dim=0; dim < 3; dim++){
					frictionCoefficient[i][zahl][dim]+=	(v[i][dim]-v12[i][dim])/(h/2);
					frictionCoefficient[i][zahl][dim]=(-1)*frictionCoefficient[i][zahl][dim] * MCa /(v12[i][dim]);
				}
			}
		}*/


		if(rms && doRMS){
			for(int dim=0;dim<3;++dim){	
				rms[dim]=rms[dim]*(*iter-1)/(*iter) + v[rmsion][dim]*v[rmsion][dim]/(*iter);
			}
		}


	}
}


//EXPORT double GetDamping(){
//	double eta=MCa*(1.-freezing)/h;
//	return eta;
//}




double Tsqueeze[CPUCNT];//initSim null setzen!!!!!
bool startSqueeze[CPUCNT];
double timer[CPUCNT]={};
double timer2[CPUCNT]={};
int SmoothSqueezeCtr[CPUCNT];
glob(bool,potKick,false);
glob(bool,stepKick,false);
glob(bool,scrambleMicroMotion,false);
glob(double,DeltaTSqueeze,8.33333e-08);
glob(int,AnzKickSq,1);
glob(int,SmoothingSq,20);
glob(bool,doubleKick,false);



void Kibblesim::Squeeze(double t, double v[maxionanz][3], double x[maxionanz][3], double a[maxionanz][3], double deltat){
	//calc time hot heat bath is in interaction with ion
	//double RadFreq=2969444.;
	//double eta= squeezing / ((1/HCFreq*acos(HCLimit)/pi) * RadFreq);
	
	if(potKick){
		if(doubleKick){
			if(TestSqueezeDone[simnum-1]==false && HotLaserOn[simnum-1]==false){
				rfMod[simnum-1]=1.+squeezing*(SmoothSqueezeCtr[simnum-1]/SmoothingSq);
				timer[simnum-1]=t+DeltaTSqueeze;	//viertel periode von 1/3MHz
				SmoothSqueezeCtr[simnum-1]++;
				if(SmoothSqueezeCtr[simnum-1]==SmoothingSq+1){
					TestSqueezeDone[simnum-1]=true;
				}
			}
			
			
			if(t>timer[simnum-1]){
				if(SmoothSqueezeCtr[simnum-1]>(-SmoothingSq))SmoothSqueezeCtr[simnum-1]--;
				rfMod[simnum-1]=1.+squeezing*(SmoothSqueezeCtr[simnum-1]/SmoothingSq);
			
			}	

			if(t>(timer[simnum-1]+DeltaTSqueeze)){
				if(SmoothSqueezeCtr[simnum-1]<0)SmoothSqueezeCtr[simnum-1]++;
				rfMod[simnum-1]=1.+squeezing*(SmoothSqueezeCtr[simnum-1]/SmoothingSq);
			
			}	
		}

		else if(TestSqueezeDone[simnum-1]==false && HotLaserOn[simnum-1]==false){
			rfMod[simnum-1]=1.+squeezing*(SmoothSqueezeCtr[simnum-1]/SmoothingSq);
			timer[simnum-1]=t+DeltaTSqueeze;	//viertel periode von 1/3MHz
			SmoothSqueezeCtr[simnum-1]++;
			if(SmoothSqueezeCtr[simnum-1]==SmoothingSq+1){
				TestSqueezeDone[simnum-1]=true;
			}
			
			
		}
		if(t>timer[simnum-1]){
			if(SmoothSqueezeCtr[simnum-1]>0)SmoothSqueezeCtr[simnum-1]--;
			rfMod[simnum-1]=1.+squeezing*(SmoothSqueezeCtr[simnum-1]/SmoothingSq);
			
		}
	}

		


	//else if(stepKick){
	//	if(TestSqueezeDone[simnum-1]==false && HotLaserOn[simnum-1]==false){
	//		rfMod[simnum-1]=1.+squeezing;
	//		timer[simnum-1]=t+DeltaTSqueeze;	//viertel periode von 1/3MHz
	//		TestSqueezeDone[simnum-1]=true;
	//		SmoothSqueezeCtr[simnum-1]=0;
	//	}
	//	if(t>timer[simnum-1] && SmoothSqueezeCtr[simnum-1]<AnzKickSq){
	//		rfMod[simnum-1]=1.;
	//		timer2[simnum-1]=t+DeltaTSqueeze;
	//		timer[simnum-1]=t+2*DeltaTSqueeze;
	//		SmoothSqueezeCtr[simnum-1]++;
	//	}
	//	if(t>timer2[simnum-1] && SmoothSqueezeCtr[simnum-1]<AnzKickSq){
	//		rfMod[simnum-1]=1.+squeezing;
	//	}
	//}


	//else{
	//	if(TestSqueezeDone[simnum-1]==false && HotLaserOn[simnum-1]==false){
	//		for(int ion=0;ion<numIons;ion++){
	//			for(int dim=0;dim<2;dim++){
	//				v[ion][dim]=v[ion][dim]*(1+squeezing);
	//				x[ion][dim]=x[ion][dim]/(1+squeezing);
	//				a[ion][dim]=0;
	//			}
	//		}
	//		TestSqueezeDone[simnum-1]=true;
	//		
	//	}
	//}


}


EXPORT string & noise(int noise, double min, double max, double amp){
	static string uni("Uniform(min,max,amp)");
	static string gauss("Gaussian(x,s,amp)");
	static string expo("Exponential(tau,-,amp)");
	static string warning("N/A choose 0..2");
	noiseFunction=noise;
	noiseParam[0]=min;
	noiseParam[1]=max;
	noiseParam[2]=amp;
	int a;
	switch(a){
		case 0: 
			return uni;
		case 1: 
			return gauss;
		case 2:
			return expo;
		default:
			return warning;
	}
}

void Kibblesim::NoiseHeat(){
	
	switch(noiseFunction){
		case 0:
			noiseRF[simnum-1]=noiseParam[2]*RandomReal(noiseParam[0],noiseParam[1]);
			break;
		case 1:
			noiseRF[simnum-1]=noiseParam[2]*RandomGauss(noiseParam[0],noiseParam[1]);
			break;
		case 2:
			noiseRF[simnum-1]=noiseParam[2]*RandomExp(noiseParam[0]);
			break;
		default:
			noiseRF[simnum-1]=0;
	}


}

double NormalizedFreezing=0;



void Kibblesim::initSim(){
	if(numIons>maxionanz || numIons<=0){
		cout<<"falsche Ionenzahl!!!"<<endl;
		numIons=maxionanz;
	}
	rfamp=rfampGet();
	simstepsMT=simsteps;
	rfvoltage=0;
	t=0;
	PhotonCountSet(0);
	for(int i=0;i<maxionanz;i++) {
		lebenszeit[i]=0;
		a[i][0]=0;a[i][1]=0;a[i][2]=0;
		x[i][0]=0;x[i][1]=0;x[i][2]=0;
		v[i][0]=0;v[i][1]=0;v[i][2]=0;
	}
	SmoothSqueezeCtr[simnum-1]=AnzKickSq;
	h=hglob;
	
	if(scrambleMicroMotion)phi=RandomReal(0,2*pi);
	else phi=0;
	//HCphi=RandomReal(0,2*pi);
	HCphi=0;
	
	t0=t0glob;
	for(int i=0;i<CPUCNT;i++){
		rfMod[i]=1.;
		timer[i]=0.;
		SmoothSqueezeCtr[i]=0;
	}

	NormalizedFreezing=pow(freezing,(h/1.0e-7));
	noiseRF[simnum-1]=0.;

	
	for(int i=0;i<numIons;++i){
		x[i][0]=0;//RandomReal(-1e6,1.e6);
		x[i][1]=0;//RandomReal(-1e6,1.e6);
		if(numIons==1) x[i][2]=zOffset;
		else x[i][2]=((i-0.5*(numIons-1))*startiondist/((numIons-1)*0.5)) + startionshift +  ((i-0.5*(numIons-1))*startiondist/((numIons-1)*0.5))*(i-0.5*(numIons-1))*(i-0.5*(numIons-1))*startionasym;//random.Uniform(-1e4,1e4);

		v[i][0]=0;
		v[i][1]=0;
		v[i][2]=0;
		a[i][0]=0;
		a[i][1]=0;
		a[i][2]=0;	
	}
	t=0;
}

int PlotCountPoints[CPUCNT];
int fileNumber[CPUCNT];
glob(bool, showHCAmp, false);
glob(double, smallstepsT,endzeitGet());
glob(int,smallSimsteps, simsteps);

void Kibblesim::Sim(){
	
	initSim();
	dosimSet(true);
	dotrapfreqSet(false);
	ionplotcnt[simnum-1]=0;
	t=0.;
	initialdamp=true;
	dampfinished[simnum-1]=false;
	detuneMT[simnum-1]=0;
	satMT[simnum-1]=0;
	Tsqueeze[simnum-1]=0;
	startSqueeze[simnum-1]=false;
	TestSqueezeDone[simnum-1]=false;
	restarted[simnum-1]=false;
	TestSqueezeDone[simnum-1]=true;

	trand[simnum-1]=RandomReal(0,(0.05*endzeit));
	//cout<<trand[simnum-1]<<endl;
	
	
	int iter=0;
	double rms[]={0,0,0};
	do{
		propagateForwardVerlet(t,simstepsMT,x,v,rms,&iter);
		//simsteps is number of propagations without drawing a point in the plot
		if(t>smallstepsT)simstepsMT=smallSimsteps;
		//smallsteps --> drwing trajectory with higher resolution
		
		tttt[simnum-1][ionplotcnt[simnum-1]]=t;

		for(int i=0;i<numIons;i++){
			for(int dim=0;dim<3;dim++){
				if(MOT_mode){
					//ERad für MOT_HE
					ERad[simnum-1][i][ionplotcnt[simnum-1]]=( x[i][0]*x[i][0] + x[i][1]*x[i][1] )*( 2*pi*(1632*exp(-x[i][2]/0.00425)+1880))*( 2*pi*(1632*exp(-x[i][2]/0.00425)+1880)) + (v[i][0]*v[i][0] + v[i][1]*v[i][1]);
					ETotal[simnum-1][i][ionplotcnt[simnum-1]]=( x[i][0]*x[i][0] + x[i][1]*x[i][1] )*( 2*pi*(1632*exp(-x[i][2]/0.00425)+1880))*( 2*pi*(1632*exp(-x[i][2]/0.00425)+1880)) + (v[i][0]*v[i][0] + v[i][1]*v[i][1]) + v[i][2]*v[i][2] + (x[i][2]*x[i][2] * (2*pi*HCFreq)*(2*pi*HCFreq));
				}
				else{
					//ERad für HE Ion
					ERad[simnum-1][i][ionplotcnt[simnum-1]]=( x[i][0]*x[i][0] + x[i][1]*x[i][1] )*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)) + (v[i][0]*v[i][0] + v[i][1]*v[i][1]);
					ETotal[simnum-1][i][ionplotcnt[simnum-1]]=( x[i][0]*x[i][0] + x[i][1]*x[i][1] )*( 2*pi*(x[i][2]*(-2241983470.)+2969444.))*( 2*pi*(x[i][2]*(-2241983470.)+2969444.)) + (v[i][0]*v[i][0] + v[i][1]*v[i][1]) + v[i][2]*v[i][2] + (x[i][2]*x[i][2] * (2*pi*36003)*(2*pi*36003));
				}
				if(useSimSingle){
					yyyy[dim][0][i][ionplotcnt[0]]=x[i][dim];
					//vvvv[dim][0][i][ionplotcnt[0]]=v[i][dim];
					HCAmpPlot[ionplotcnt[0]]=HCAmp[0]*reduceAmp;
				}
				else{
					yyyy[dim][simnum-1][i][ionplotcnt[simnum-1]]=x[i][dim];
					vvvv[dim][simnum-1][i][ionplotcnt[simnum-1]]=v[i][dim];
					//HCAmpPlot[ionplotcnt[simnum-1]]=HCAmp[simnum-1]*reduceAmp;
				}
			}
		}
		ionplotcnt[simnum-1]++;
			
		if(displaysim==simnum && useSimSingle==false){
			if(displayiontrap){
				switch(numIons){
					case 5:
						plot5.setData(tttt[0],yyyy[2][0][4]);
						plot5.setNum(ionplotcnt[simnum-1]);
						plot5.update();
					case 4:
						plot4.setData(tttt[0],yyyy[2][0][3]);
						plot4.setNum(ionplotcnt[simnum-1]);
						plot4.update();
					case 3:
						plot3.setData(tttt[0],yyyy[2][0][2]);
						plot3.setNum(ionplotcnt[simnum-1]);
						plot3.update();
					case 2:
						plot2.setData(tttt[0],yyyy[2][0][1]);
						plot2.setNum(ionplotcnt[simnum-1]);
						plot2.update();
					case 1:
						plot1.setData(tttt[0],yyyy[2][0][0]);
						plot1.setNum(ionplotcnt[simnum-1]);
						plot1.update();
					default: break;
				}
				if(displayRad){

					switch(numIons){
						case 7:
							plotX7.setData(tttt[0],yyyy[0][0][6]);
							plotX7.setNum(ionplotcnt[simnum-1]);
							plotX7.update();
						case 6:
							plotX6.setData(tttt[0],yyyy[0][0][5]);
							plotX6.setNum(ionplotcnt[simnum-1]);
							plotX6.update();
						case 5:
							plotX5.setData(tttt[0],yyyy[0][0][4]);
							plotX5.setNum(ionplotcnt[simnum-1]);
							plotX5.update();
						case 4:
							plotX4.setData(tttt[0],yyyy[0][0][3]);
							plotX4.setNum(ionplotcnt[simnum-1]);
							plotX4.update();
						case 3:
							plotX3.setData(tttt[0],yyyy[0][0][2]);
							plotX3.setNum(ionplotcnt[simnum-1]);
							plotX3.update();
						case 2:
							plotX2.setData(tttt[0],yyyy[0][0][1]);
							plotX2.setNum(ionplotcnt[simnum-1]);
							plotX2.update();
						case 1:
							plotX1.setData(tttt[0],yyyy[0][0][0]);
							plotX1.setNum(ionplotcnt[simnum-1]);
							plotX1.update();
						default: break;
						}


					switch(numIons){
						case 7:
							plotY7.setData(tttt[0],yyyy[1][0][6]);
							plotY7.setNum(ionplotcnt[simnum-1]);
							plotY7.update();
						case 6:
							plotY6.setData(tttt[0],yyyy[1][0][5]);
							plotY6.setNum(ionplotcnt[simnum-1]);
							plotY6.update();
						case 5:
							plotY5.setData(tttt[0],yyyy[1][0][4]);
							plotY5.setNum(ionplotcnt[simnum-1]);
							plotY5.update();
						case 4:
							plotY4.setData(tttt[0],yyyy[1][0][3]);
							plotY4.setNum(ionplotcnt[simnum-1]);
							plotY4.update();
						case 3:
							plotY3.setData(tttt[0],yyyy[1][0][2]);
							plotY3.setNum(ionplotcnt[simnum-1]);
							plotY3.update();
						case 2:
							plotY2.setData(tttt[0],yyyy[1][0][1]);
							plotY2.setNum(ionplotcnt[simnum-1]);
							plotY2.update();
						case 1:
							plotY1.setData(tttt[0],yyyy[1][0][0]);
							plotY1.setNum(ionplotcnt[simnum-1]);
							plotY1.update();
						default: break;
						}
				}
				if(showHCAmp){
					plotHC.setData(tttt[0],HCAmpPlot);
					plotHC.setNum(ionplotcnt[simnum-1]);
					plotHC.update();
				}

				IonTrapDisplay.update();
			}
			if(ionplotcnt[simnum-1]>ionplotanz){
				cout <<"array too small!"<<endl;
				ionplotcnt[simnum-1]=0;

			
			}
		}
		
	}while(dosim && t<endzeit);

	simstepsMT=simsteps;

	PlotCountPoints[simnum-1]=ionplotcnt[simnum-1];

	if(doRMS){
		cout <<"-------------"<<endl;
		cout <<"rms x "<<sqrt(rms[0])<<endl;
		cout <<"rms y "<<sqrt(rms[1])<<endl;
		cout <<"rms z "<<sqrt(rms[2])<<endl;
	}
	if(useSimSingle){
		cout<<"generating Plots... "<<ionplotcnt[0]<<endl;
		cout<<"X";
		plotX1.setData(tttt[0],yyyy[0][0][0]);
		plotX1.setNum(ionplotcnt[0]);
		plotX1.update();
		cout<<"Y";
		plotY1.setData(tttt[0],yyyy[1][0][0]);
		plotY1.setNum(ionplotcnt[0]);
		plotY1.update();
		cout<<"Z"<<endl;
		plot1.setData(tttt[0],yyyy[2][0][0]);
		plot1.setNum(ionplotcnt[0]);
		if(showHCAmp){
			cout<<"HC...";
			plotHC.setData(tttt[0],HCAmpPlot);
			plotHC.setNum(ionplotcnt[0]);
			plot1.update();
			plotHC.update();
			cout<<"done"<<endl;
		}
		else plot1.update();
	}

}


EXPORT void UpdatePlot(){
	switch(numIons){
			case 18:
				plot18.setNum(ionplotcnt[0]);
				plot18.update();
			case 17:
				plot17.setNum(ionplotcnt[0]);
				plot17.update();
			case 16:
				plot16.setNum(ionplotcnt[0]);
				plot16.update();
			case 15:
				plot15.setNum(ionplotcnt[0]);
				plot15.update();
			case 14:
				plot14.setNum(ionplotcnt[0]);
				plot14.update();
			case 13:
				plot13.setNum(ionplotcnt[0]);
				plot13.update();
			case 12:
				plot12.setNum(ionplotcnt[0]);
				plot12.update();
			case 11:
				plot11.setNum(ionplotcnt[0]);
				plot11.update();
			case 10:
				plot10.setNum(ionplotcnt[0]);
				plot10.update();
			case 9:
				plot9.setNum(ionplotcnt[0]);
				plot9.update();
			case 8:
				plot8.setNum(ionplotcnt[0]);
				plot8.update();
			case 7:
				plot7.setNum(ionplotcnt[0]);
				plot7.update();
			case 6:
				plot6.setNum(ionplotcnt[0]);
				plot6.update();
			case 5:
				plot5.setNum(ionplotcnt[0]);
				plot5.update();
			case 4:
				plot4.setNum(ionplotcnt[0]);
				plot4.update();
			case 3:
				plot3.setNum(ionplotcnt[0]);
				plot3.update();
			case 2:
				plot2.setNum(ionplotcnt[0]);
				plot2.update();
			case 1:
				plot1.setNum(ionplotcnt[0]);
				plot1.update();
			default: break;
		}
		if(displayRad){

			switch(numIons){
				case 7:
					plotX7.setNum(ionplotcnt[0]);
					plotX7.update();
				case 6:
					plotX6.setNum(ionplotcnt[0]);
					plotX6.update();
				case 5:
					plotX5.setNum(ionplotcnt[0]);
					plotX5.update();
				case 4:
					plotX4.setNum(ionplotcnt[0]);
					plotX4.update();
				case 3:
					plotX3.setNum(ionplotcnt[0]);
					plotX3.update();
				case 2:
					plotX2.setNum(ionplotcnt[0]);
					plotX2.update();
				case 1:
					plotX1.setNum(ionplotcnt[0]);
					plotX1.update();
				default: break;
				}


			switch(numIons){
				case 7:
					plotY7.setNum(ionplotcnt[0]);
					plotY7.update();
				case 6:
					plotY6.setNum(ionplotcnt[0]);
					plotY6.update();
				case 5:
					plotY5.setNum(ionplotcnt[0]);
					plotY5.update();
				case 4:
					plotY4.setNum(ionplotcnt[0]);
					plotY4.update();
				case 3:
					plotY3.setNum(ionplotcnt[0]);
					plotY3.update();
				case 2:
					plotY2.setNum(ionplotcnt[0]);
					plotY2.update();
				case 1:
					plotY1.setNum(ionplotcnt[0]);
					plotY1.update();
				default: break;
				}
		}
		if(showHCAmp){
			//plotHC.setData(tttt[0],HCAmp);
			plotHC.setData(tttt[0],HCAmpPlot);
			plotHC.setNum(ionplotcnt[0]);
			plotHC.update();
		}
		else{
			double Nuller[ionplotanz]={};
			for(int j=0;j<ionplotanz;j++)Nuller[j]=0;

			//plotHC.setData(tttt[0],Nuller);
			plotHC.setNum(0);
			plotHC.update();
		}
	}


glob(double,exciteFreq,1e-7)

double Kibblesim::trapFreq(int axis,int anzosci){
	
	bool sauberMachen=false;
	if(yyyy==0){
		sauberMachen=true;
		yyyy=makearrayPos<Tpos,3>();//--> in function
		vvvv=makearrayPos<Tpos,3>();
		tttt=makearrayTime<Ttime,CPUCNT>();//--> in function
	}

	dotrapfreqSet(true);
	if(axis<0) return 0;
	if(axis>2) return 0;
	initSim();
	ionplotcnt[0]=0;
	int oldNumIons=numIons;
	numIons=1;
	int oldsimsteps=simsteps;
	
	bool HEstatusBefore=HeatEngineGet();
	HeatEngine=false;
	bool laserStatusBefore=laseronGet();
	laseron=false;
	bool verletStatusBefore=verletGet();
	verlet=true;
	
	x[0][0]=0;		
	x[0][1]=0;		
	x[0][2]=zOffset;
	v[0][0]=0;		
	v[0][1]=0;		
	v[0][2]=0;		
	a[0][0]=0;
	a[0][1]=0;
	a[0][2]=0;
	t=0;

	if(axis==2)x[0][axis]=zOffset+exciteFreq;//1e-7
	else x[0][axis]=exciteFreq;

	double valuebefore=0;
	double tstart=-1;
	int oscicnt=0;
	while(dotrapfreq){
		propagateForwardVerlet(t,simsteps,x,v);
		tttt[0][ionplotcnt[0]]=t;
		yyyy[0][0][0][ionplotcnt[0]]=x[0][0];
		yyyy[1][0][0][ionplotcnt[0]]=x[0][1];
		yyyy[2][0][0][ionplotcnt[0]]=x[0][2];
		
		if(axis!=2){
			if(valuebefore<0 && yyyy[axis][0][0][ionplotcnt[0]]>0){
				if(oscicnt==0) tstart=t;
				
				oscicnt++;
				
				if(oscicnt==anzosci){
					dotrapfreqSet(false);
					break;
				}
			}
		}
		
		else{
			if(valuebefore<zOffset && yyyy[axis][0][0][ionplotcnt[0]]>zOffset){
				if(oscicnt==0) tstart=t;
				
				oscicnt++;
				
				if(oscicnt==anzosci){
					dotrapfreqSet(false);
					break;
				}
			}
		}

		valuebefore=yyyy[axis][0][0][ionplotcnt[0]];
		ionplotcnt[0]++;
		
		PlotCountPoints[simnum-1]=ionplotcnt[simnum-1];


		if(displayiontrap){
			plot1.setData(tttt[0],yyyy[2][0][0]);
			plot1.setNum(ionplotcnt[0]);
			plot1.update();
			plotX1.setData(tttt[0],yyyy[0][0][0]);
			plotX1.setNum(ionplotcnt[0]);
			plotX1.update();
			plotY1.setData(tttt[0],yyyy[1][0][0]);
			plotY1.setNum(ionplotcnt[0]);
			plotY1.update();

			IonTrapDisplay.update();
		}
		if(ionplotcnt[0]>ionplotanz) ionplotcnt[0]=0;
	}
	plot1.setData(tttt[0],yyyy[2][0][0]);
	plot1.setNum(ionplotcnt[0]);
	plot1.update();
	plotX1.setData(tttt[0],yyyy[0][0][0]);
	plotX1.setNum(ionplotcnt[0]);
	plotX1.update();
	plotY1.setData(tttt[0],yyyy[1][0][0]);
	plotY1.setNum(ionplotcnt[0]);
	plotY1.update();
	if(showHCAmp){
		plotHC.setData(tttt[0],HCAmpPlot);
			plotHC.setNum(ionplotcnt[0]);
			plot1.update();
			plotHC.update();
	}


	numIons=oldNumIons;
	simsteps=oldsimsteps;

	if(sauberMachen==true){
		delete[] yyyy;
		delete[] vvvv;
		delete[] tttt;
		yyyy=0;
	}

	
	verlet=verletStatusBefore;
	laseron=laserStatusBefore;
	HeatEngine=HEstatusBefore;
	return 1./((t-tstart)/double(anzosci))/1000.;
	
}



void Kibblesim::Freeze(double t, double v[maxionanz][3]){
		
	/*double eta=MCa*(1-freezing)/h;
	double sigma=sqrt(k_B*tempT/MCa);
	
	double noise_magnitude = sqrt(tempT*k_B*eta/(MCa*MCa))*sqrt(h);*/


	for(int ind=0; ind<numIons;ind++){
		for(int dim=0; dim<2; dim++){
			v[ind][dim]*=NormalizedFreezing;
			//else if(t>freezeTime/2.){
			//	v[ind][dim]*=0.9999;
			//	//v[ind][dim]=v[ind][dim]*(1.-(freezeFactor*(t-t0Freeze)/(tauFreeze)));
			//	//cout<<1.-(freezeFactor*(t-t0Freeze)/(tauFreeze))<<endl;
			//}
		}
	}

}

void Kibblesim::FreezeRad(double t, double v[maxionanz][3]){
	for(int ind=0; ind<numIons;ind++){
		v[ind][0]*=NormalizedFreezing;
		v[ind][1]*=NormalizedFreezing;
	}
}

void Kibblesim::FreezeAx(double t, double v[maxionanz][3]){
	for(int ind=0; ind<numIons;ind++){
		v[ind][2]*=NormalizedFreezing;
	}
}


Kibblesim simarray[CPUCNT];//this line if commented out removes the problem
	
	

EXPORT void SimSingle(){
	bool sauberMachen=false;
	if(yyyy==0){
		cout<<"create position arrays..."<<endl;
		sauberMachen=true;
		vvvv=makearrayPos<Tvel,3>();
		yyyy=makearrayPos<Tpos,3>();//--> in function
		tttt=makearrayTime<Ttime,CPUCNT>();//--> in function
	}
	
	if(yyyy!=0 && tttt!=0 && vvvv!=0){
		cout<<"alles OK"<<endl;
		useSimSingleSet(true);
		simarray[0].Sim();
	}
	else cout<<"ABBRUCH"<<endl;

	if(sauberMachen==true){
		delete[] yyyy;
		delete[] vvvv;
		delete[] tttt;
		yyyy=0;
		vvvv=0;
		tttt=0;
	}
	
}

EXPORT double trapFreq(int axis,int anzosci){
		return simarray[0].trapFreq(axis,anzosci);
}

EXPORT double scatterCount(){
		return simarray[0].scatterCount();
}

template< class type>
inline std::string to_string( const type & value,int width=0)
{
    std::ostringstream streamOut;
	if (width!=0)   streamOut <<setfill('0')<<setw(width) <<value;
	else streamOut <<value;
    return streamOut.str();
}
glob(bool,dowrite,false);
glob(bool,dowriteAll,false);

EXPORT void displaystep(int step){
	if(step<0) return;
	if(step>ionplotcnt[0]) return;
	if(yyyy!=0){
		for(int m=0;m<numIons;m++)
			for(int dim=0;dim<3;dim++)
				;//simarray[0].x[m][dim]=yyyy[dim][0][m][step];

	}
	IonTrapDisplay.update();
}


EXPORT int jumptolast(){
	cout<<"jumptolast: "<<ionplotcntmax<<endl;
	//displaystep(ionplotcntmax-1);
	return ionplotcntmax-1;
}

glob(bool,doAuswertung,false);


EXPORT int  readresult(int fnum, int psc){
	boost::mutex::scoped_lock scoped_lock(io_mutex);
	
	bool sauberMachen=false;
	if(yyyy==0){
		cout<<"create position arrays..."<<endl;
		sauberMachen=true;
		yyyy=makearrayPos<Tpos,3>();//--> in function
		vvvv=makearrayVel<Tvel,3>();
		tttt=makearrayTime<Ttime,CPUCNT>();//--> in function
	}
	
	
	string fname(runname);
	fname.append("_");
	size_t found;
	found=name.find_last_of("/\\");
	fname.append(name.substr(found+1));
	if(fnum==0)fname.append(to_string("mean"));
	else if(fnum==-1)fname.append(to_string("ERad"));
	else{
		switch(psc){
			case 0:	fname.append(to_string("X"));
					fname.append(to_string(fnum));
					break;
			case 1: fname.append(to_string("V"));
					fname.append(to_string(fnum));
					break;
			default: break;
		}
	}
	fname.append(string(".tsv"));
	ifstream ifs;
	ifs.open(fname.c_str(),std::ios_base::in);
	cout <<fname<<endl;
	int simpoint=0;
	if(ifs.fail()){
		return 1;
	}
	while(!ifs.eof()){
		ifs>>tttt[0][simpoint];
		for(int i=0;i<numIons;i++) {
			if(fnum==-1) ifs >>ERadMean[i][simpoint];
			else{
				for(int dim=0;dim<3;dim++){
					ifs >>yyyy[dim][0][i][simpoint];				
				}
			}
		}
		simpoint++;
	}
	ifs.close();
	ionplotcnt[0]=simpoint-1;
	if(doAuswertung==false){
		ionplotcntmax=ionplotcnt[0];
		cout<<ionplotcntmax<<endl;
		if(fnum>0)displaysimset(1);
		//displaystep(ionplotcnt[0]-1);
		
		else if(fnum==-1){
			plotERad.setData(tttt[0],ERadMean[0]);
			plotERad.setNum(ionplotcnt[0]);
			plotERad.update();

			plotETotal.setData(tttt[0],ETotMean[0]);
			plotETotal.setNum(ionplotcnt[0]);
			plotETotal.update();
		}
		else if(fnum==0){
			plotMeanZ.setData(tttt[0],yyyy[2][0][0]);
			plotMeanZ.setNum(ionplotcnt[0]);
			plotMeanZ.update();
			plotMeanX.setData(tttt[0],yyyy[0][0][0]);
			plotMeanX.setNum(ionplotcnt[0]);
			plotMeanX.update();
		}
	}
	//cout <<simpoint<<endl;
	if(sauberMachen==true){
		delete[] yyyy;
		delete[] vvvv;
		delete[] tttt;
		yyyy=0;
		vvvv=0;
		tttt=0;
	}
	return 0;
}


void writeresult(int simnum,int fnum){
	boost::mutex::scoped_lock scoped_lock(io_mutex);
	
	string fname(runname);
	fname.append("_");
	size_t found;
	found=name.find_last_of("/\\");
	fname.append(name.substr(found+1));
	fname.append(to_string("X"));
	fname.append(to_string(fnum));
	fname.append(string(".tsv"));
	ofstream ofs;
	ofs.precision(10);
	ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
	int anzsimpoints=ionplotcnt[simnum-1];
	for(int simpoint=0;simpoint<anzsimpoints;simpoint++){
		ofs<<tttt[simnum-1][simpoint]<<"\t";
		for(int i=0;i<numIons;i++) {
			for(int dim=0;dim<3;dim++){
				ofs<<yyyy[dim][simnum-1][i][simpoint]<<"\t";
				//yyyy[dim][simnum-1][i][ionplotcnt[simnum-1]]=x[i][dim];	//nur zum vergleich
			}
		}
		ofs<<endl;
	}
	ofs.close();

	string fname2(runname);
	fname2.append("_");
	size_t found2;
	found2=name.find_last_of("/\\");
	fname2.append(name.substr(found2+1));
	fname2.append(to_string("V"));
	fname2.append(to_string(fnum));
	fname2.append(string(".tsv"));
	ofstream ofs2;
	ofs2.precision(10);
	ofs2.open(fname2.c_str(),std::ios_base::out|std::ios_base::trunc);
	for(int simpoint=0;simpoint<anzsimpoints;simpoint++){
		ofs2<<tttt[simnum-1][simpoint]<<"\t";
		for(int i=0;i<numIons;i++) {
			for(int dim=0;dim<3;dim++){
				ofs2<<vvvv[dim][simnum-1][i][simpoint]<<"\t";
			}
		}
		ofs2<<endl;
	}
	ofs2.close();


	displaysimset(0);
}

string sTime;
bool sGiven=false;
glob(int,numberCycles,1);

EXPORT void DrawCycle(int n, int mittel){
	
	bool sauberMachen=false;
	if(yyyy==0){
		sauberMachen=true;
		yyyy=makearrayPos<Tpos,3>();//--> in function
		vvvv=makearrayPos<Tpos,3>();
		tttt=makearrayTime<Ttime,CPUCNT>();//--> in function
	}
	
	readresult(0,0);
	readresult(-1,0);
	int ion=0;
	int p=0;
	int nstart=0;
	int nstop=0;

	

	OmegaRadPlot=(double*) calloc(cyclePtsGet(),sizeof(double));
	cycleRadPlot=(double*) calloc(cyclePtsGet(),sizeof(double));

	//--> yyyy[dim][0][i][simpoint]
	n-=1;
	while(p>-1){		//gehe von hinten durch yyyy array. Setze stop und start werte für einzelne Schwingungsperioden. Über jede Periode wird ein cycle berechnet
		if (p==numberCycles && yyyy[2][0][ion][n]>yyyy[2][0][ion][n-1] && yyyy[2][0][ion][n]>yyyy[2][0][ion][n+1] && yyyy[2][0][ion][n]>0){//finde maximum
			nstart=n;
			p=-1;
			break;
		}
		if (p!=-1 && yyyy[2][0][ion][n]>yyyy[2][0][ion][n-1] && yyyy[2][0][ion][n]>yyyy[2][0][ion][n+1] && yyyy[2][0][ion][n]>0){//finde zweites maximum
			if(p==0)nstop=n;
			p++;
		}
		n--;
	}	
	
	double periodT=tttt[0][nstop]-tttt[0][nstart];
	int ndiff=nstop-nstart;

	cout<<"Zeit: "<<periodT<<" numberPoints: "<<ndiff<<endl;
	cout<<"Start: "<<nstart<<"\t Stopp: "<<nstop<<endl;

	int ntemp=0;
	double ERadMittel=0;

	for(int r=0;r<cyclePts;r++){
		ntemp=nstop-r*ndiff/cyclePts;
		//cycleZPos[r]=yyyy[2][0][ion][ntemp];

		//////////////////////// hier steht die Radial-Axial Beziehung
		//OmegaRadPlot[r]=yyyy[2][0][ion][ntemp]*(-2241983470.)+2969444.;
		//////////////////////////
		
		//nur mit ORT auf x achse
		OmegaRadPlot[r]=yyyy[2][0][ion][ntemp];


		ERadMittel=0.;
		if(mittel==0){
			cycleRadPlot[r]=ERadMean[ion][ntemp];
		}
		else{
			if(mittel<0)mittel=1;
			for(int l=(ntemp-mittel);l<(ntemp+mittel);l++){
				ERadMittel+=ERadMean[ion][l];
			}
			cycleRadPlot[r]=ERadMittel/(2*mittel);
		}
	}
	cout<<mittel<<endl;

	
	//Draw cycle
	cyclePlot.setData(OmegaRadPlot,cycleRadPlot);
	cyclePlot.setNum(cyclePts);
	cyclePlot.update();

	if(sGiven==false){
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		sTime=to_string(timeinfo->tm_year-100,2)+to_string(timeinfo->tm_mon+1,2)+to_string(timeinfo->tm_mday,2)+"_"+to_string(timeinfo->tm_hour,2)+to_string(timeinfo->tm_min,2)+to_string(timeinfo->tm_sec,2);
		//CopyFile((path+"\\session.mcp").toAscii().data(),(path+"\\session"+s.c_str()+".mcp").toAscii().data(),false);
	}

	
	string fname2(runname);
	fname2.append("_");
	size_t found2;
	found2=name.find_last_of("/\\");
	fname2.append(name.substr(found2+1));
	fname2.append(to_string(sTime));
	fname2.append(to_string("_Cycle"));
	fname2.append(string(".tsv"));
	ofstream ofs2;
	ofs2.precision(10);
	ofs2.open(fname2.c_str(),std::ios_base::out|std::ios_base::trunc);
	for(int r=0;r<cyclePts;r++){
		ofs2<<OmegaRadPlot[r]<<"\t";
		ofs2<<cycleRadPlot[r]<<"\t";
		ofs2<<endl;
	}
	ofs2.close();



	free(OmegaRadPlot);
	free(cycleRadPlot);


	if(sauberMachen==true){
		delete[] yyyy;
		delete[] vvvv;
		delete[] tttt;
		yyyy=0;
		vvvv=0;
		tttt=0;
	}
}



EXPORT void Histogram(int axis, int n, int numHist){
	bool sauberMachen=false;
	
	int fnum=0;
	int readreturn=0;
	int deltaSteps=0;
	int stepEnde=0;
	int stepStart=0;
	int stepTemp=0;
	//double Hist[3][maxionanz]	//alloc!
	doAuswertung=true;

	if(yyyy==0){
		cout<<"create position arrays..."<<endl;
		sauberMachen=true;
		yyyy=makearrayPos<Tpos,3>();//--> in function
		vvvv=makearrayVel<Tvel,3>();
		tttt=makearrayTime<Ttime,CPUCNT>();//--> in function
	}


	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	sTime=to_string(timeinfo->tm_year-100,2)+to_string(timeinfo->tm_mon+1,2)+to_string(timeinfo->tm_mday,2)+"_"+to_string(timeinfo->tm_hour,2)+to_string(timeinfo->tm_min,2)+to_string(timeinfo->tm_sec,2);

	//open filestream for histogram
	string fname(runname);
	fname.append("_");
	size_t found;
	found=name.find_last_of("/\\");
	fname.append(name.substr(found+1));
	fname.append(to_string(sTime));
	fname.append(to_string("_Hist"));
	fname.append(to_string(n));
	fname.append(string(".tsv"));
	ofstream ofs;
	ofs.precision(10);
	ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
	cout<<fname<<endl;

	//write coordinates of all realizations for step n into file
	while(readreturn!=1){
		readreturn=readresult(fnum,0);//--> schreibe x in yyyy.
		
		if(readreturn==1){
			cout<<"File does not exist"<<endl;
			break;
		}
		if(n>ionplotcntmax || n<0){
			cout<<"n not in valid range of steps"<<endl;
			break;
		}
		//cout<<readreturn<<endl;
		//cout<<yyyy[2][0][0][n]<<endl;
		if(fnum==0){
			bool stop=false;
			int step=n;
			double mittelWert=0;
			double mittelWertAlt=0;
			//cout<<"find period"<<endl;
			while(stop==false){
				step--;
				mittelWert=0;
				mittelWertAlt=0;
				for(int cnt=0;cnt<11;cnt++){
					mittelWert+=yyyy[axis][0][0][step+4-cnt];
					mittelWertAlt+=yyyy[axis][0][0][step+5-cnt];
				}
				//if(yyyy[axis][0][0][step+1]<=0 && yyyy[axis][0][0][step]>0 && stepEnde!=0){
				if(Signum2(mittelWertAlt)<Signum2(mittelWert) && stepEnde!=0){
					stepStart=step;
					cout<<"Start "<<stepStart<<"\t"<<"Time ms: "<<tttt[0][step]*1000000<<endl;
					stop=true;
				}
				if(Signum2(mittelWertAlt)<Signum2(mittelWert) && stepEnde==0){
					stepEnde=step;
					cout<<"Stopp "<<stepEnde<<"\t"<<"Time ms: "<<tttt[0][step]*1000000<<endl;
					step-=50;
				}
			}
			deltaSteps=(stepEnde-stepStart)/(numHist+1);
			cout<<"deltaSteps "<<deltaSteps<<endl;
			if(deltaSteps<2){
				cout<<"resulution too bad"<<endl;
				break;
			}
			//for(int HistCount=0;HistCount<numHist;HistCount++){
			//	stepTemp=stepStart+deltaSteps*HistCount;
			//	//cout<<stepTemp<<endl;
			//	for(int dim=0;dim<3;dim++){
			//		ofs<<stepTemp<<"\t"<<stepTemp*h<<"\t"<<tttt[0][stepTemp]<<"\t"<<endl;
			//	}
			//}
		}
		else{			
			//cout<<"writing to file...";
			for(int HistCount=0;HistCount<numHist;HistCount++){
				stepTemp=stepStart+deltaSteps*HistCount;
				//cout<<stepTemp<<endl;
				for(int dim=0;dim<3;dim++){
					ofs<<yyyy[dim][0][0][stepTemp]<<"\t";
				}
			}
			//cout<<" done"<<endl;
			
			ofs<<0<<"\t";

			readreturn=readresult(fnum,1);//--> schreibe v in yyyy.
			//cout<<readreturn<<endl;
			if(readreturn==1)break;
			/*if(n>ionplotcntmax || n<0){
				cout<<"n not in valid range of steps"<<endl;
				break;
			}*/
			for(int HistCount=0;HistCount<numHist;HistCount++){
				for(int dim=0;dim<3;dim++){
					ofs<<yyyy[dim][0][0][stepStart+deltaSteps*HistCount]<<"\t";
				}
			}
			ofs<<endl;
		}

		fnum++;
	}
	
	ofs.close();

	doAuswertung=false;

	if(sauberMachen==true){
		delete[] yyyy;
		delete[] vvvv;
		delete[] tttt;
		yyyy=0;
		vvvv=0;
		tttt=0;
	}

	
}

//int numberThreads=0;
int runcnt=0;
semaphore s(0);



void startsim(int num){
	//cout<<"thread started --"<<num<<endl;
	for(int i=0;i<loopanz;i++){

		//s.enter();
		int fnum=incfilenum();
		fileNumber[num]=fnum;
	
		simarray[num].Sim();
		if(!dosim) break;
	
		cout<<"thread: "<<num<<" run:"<<i<<" file:"<<fnum<<endl;
		if(dowriteAll) writeresult(num+1,fnum);

		//s.signal();
		//s.wait();

		{
			boost::mutex::scoped_lock lock(s.mutex_, boost::try_to_lock);
			if (lock){
				runcnt++;
				//cout<<num<<" added"<<" -> "<<runcnt<<endl;
				for(int ion=0;ion<numIons;ion++){
					for(int simpoint=0;simpoint<PlotCountPoints[num];simpoint++){
						ERadMean[ion][simpoint]+=ERad[num][ion][simpoint];
						ETotMean[ion][simpoint]+=ETotal[num][ion][simpoint];
						if(PlotCountPoints[num]>ionplotanz)PlotCountPoints[num]=ionplotanz;
						for(int dim=0;dim<3;dim++){
							yyyyMean[dim][ion][simpoint]+=yyyy[dim][num][ion][simpoint];
						}						
					}
				}				
			}
			else
			{
				boost::unique_lock<boost::mutex> lock(s.mutex_);

			}
		}
		if(num==0){
			cout<<"generating Plots... "<<ionplotcnt[0]<<endl;
			cout<<"X";
			plotX1.setData(tttt[0],yyyy[0][0][0]);
			plotX1.setNum(PlotCountPoints[0]);
			plotX1.update();
			cout<<"Y";
			plotY1.setData(tttt[0],yyyy[1][0][0]);
			plotY1.setNum(PlotCountPoints[0]);
			plotY1.update();
			cout<<"Z"<<endl;
			plot1.setData(tttt[0],yyyy[2][0][0]);
			plot1.setNum(PlotCountPoints[0]);
			if(showHCAmp){
				plotHC.setData(tttt[0],HCAmpPlot);
				plotHC.setNum(ionplotcnt[0]);
				plot1.update();
				plotHC.update();
			}
			else plot1.update();
		}
		
	}
}


void writeAllStuff(){

	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	sTime=to_string(timeinfo->tm_year-100,2)+to_string(timeinfo->tm_mon+1,2)+to_string(timeinfo->tm_mday,2)+"_"+to_string(timeinfo->tm_hour,2)+to_string(timeinfo->tm_min,2)+to_string(timeinfo->tm_sec,2);

	string fname(runname);
	ofstream ofs;
	size_t found;

	if(dowrite){
		
		//write mean radial energy to textfile
		fname=runname;
		fname.append("_");
		found=name.find_last_of("/\\");
		fname.append(name.substr(found+1));
		fname.append(to_string(sTime));
		fname.append(to_string("_ERad"));
		fname.append(string(".tsv"));
		ofs.precision(10);
		ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
		for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
			ofs<<tttt[0][simpoint]<<"\t";
			for(int i=0;i<numIons;i++) {
				ofs<<ERadMean[i][simpoint]<<"\t";
			}
			ofs<<endl;
		}
		ofs.close();


		//write mean total energy to textfile
		fname=runname;
		fname.append("_");
		found=name.find_last_of("/\\");
		fname.append(name.substr(found+1));
		fname.append(to_string(sTime));
		fname.append(to_string("_ETotal"));
		fname.append(string(".tsv"));
		ofs.precision(10);
		ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
		for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
			ofs<<tttt[0][simpoint]<<"\t";
			for(int i=0;i<numIons;i++) {
				ofs<<ETotMean[i][simpoint]<<"\t";
			}
			ofs<<endl;
		}
		ofs.close();


		//write mean axial position to textfile
		fname=runname;
		fname.append("_");
		found=name.find_last_of("/\\");
		fname.append(name.substr(found+1));
		fname.append(to_string(sTime));
		fname.append(to_string("_mean"));
		fname.append(string(".tsv"));
		ofs.precision(10);
		ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
		for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
			ofs<<tttt[0][simpoint]<<"\t";
			for(int i=0;i<numIons;i++) {
				for(int dim=0;dim<3;dim++){
					ofs<<yyyy[dim][0][i][simpoint]<<"\t";
				}
			}
			ofs<<endl;
		}
		ofs.close();


		
		//write paramters to textfile
		fname=runname;
		fname.append("_");
		found=name.find_last_of("/\\");
		fname.append(name.substr(found+1));
		fname.append(to_string(sTime));
		fname.append(to_string("_parameters"));
		fname.append(string(".tsv"));
		ofs.precision(10);
		ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
		
		ofs<<"trap:\n"<<"OmegaRF:\t"<<OmegaRF<<"\n"<<"rfamp:\t"<<rfamp<<"\n"<<"Uz:\t"<<Uz<<"\n"<<"zOffset:\t"<<zOffset<<"\n"<<"angle:\t"<<angle<<"\n"<<"temperature:\t"<<temperature<<"\n"<<"Casdorff:\t"<<Casdorff<<"\n"<<"scrambleMicroMotion:\t"<<scrambleMicroMotion<<"\n"<<"usePseudoPot:\t"<<usePseudoPot<<"\n"<<endl;
		ofs<<"HeatEngine:\n"<<"HeatEngine:\t"<<HeatEngine<<"\n"<<"RK4_HeatEngine:\t"<<RK4_HeatEngine<<"\n"<<"HCFreq:\t"<<HCFreq<<"\n"<<"HCLimitHot:\t"<<HCLimitHot<<"\n"<<"HCLimitCold:\t"<<HCLimitCold<<"\n"<<"detuneHot:\t"<<detuneHot<<"\n"<<"detuneCold:\t"<<detuneCold<<"\n"<<"satHEHot:\t"<<satHEHot<<"\n"<<"satHECold:\t"<<satHECold<<"\n"<<"DeltaT:\t"<<hglob<<"\n"<<endl;
		ofs<<"Cooling Z:\n"<<"CoolZ:\t"<<coolZ<<"\n"<<"detuneZZ:\t"<<detuneZZ<<"\n"<<"satZZ:\t"<<satZZ<<"\n"<<"ZZLimit:\t"<<ZZLimit<<"\n"<<endl;
		ofs<<"Time protocol:\n"<<"Endzeit:\t"<<endzeit<<"\n"<<"LaserOnT:\t"<<LaserOnT<<"\n"<<"LaserOffT:\t"<<LaserOffT<<"\n"<<"StartEngine:\t"<<StartEngine<<"\n"<<"CoolZZTime:\t"<<coolZZTime<<"\n"<<endl;
		ofs<<"loopanz:\t"<<loopanz<<"\n"<<"detune:\t"<<detune<<"\n"<<"tau:\t"<<tau<<"\n"<<"sat:\t"<<sat<<"\n"<<endl;
		ofs<<"squeezing:\t"<<doSqueeze<<"\n"<<"squeezing:\t"<<squeezing<<"\n"<<"OmegaRes:\t"<<OmegaRes<<"\n epsilon:\t"<<epsilon<<"\n doStepSqueeze:\t"<<doStepSqueeze<<"\n doSqueeze:\t"<<doSqueeze<<"\n potKick:\t"<<potKick<<"\nstepKick:\t"<<stepKick<<"\n DeltaTSqueeze:\t"<<DeltaTSqueeze<<"\nSmoothingSq:\t"<<SmoothingSq<<endl;
		ofs<<"freezing:\t"<<doFreeze<<"\n"<<"freezing:\t"<<freezing<<"\n"<<"TColdBath:\t"<<TColdBath<<"\n"<<endl;
		ofs.close();
	}

	
	//wirte mean axial positions to temporary file for readout
	fname=runname;
	fname.append("_");
	found=name.find_last_of("/\\");
	fname.append(name.substr(found+1));
	fname.append(to_string("mean"));
	fname.append(string(".tsv"));
	ofs.precision(10);
	ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
	for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
		ofs<<tttt[0][simpoint]<<"\t";
		for(int i=0;i<numIons;i++) {
			for(int dim=0;dim<3;dim++){
				ofs<<yyyy[dim][0][i][simpoint]<<"\t";
			}
		}
		ofs<<endl;
	}
	ofs.close();

	

	//write mean radial energy to temporary file for readout
	fname=runname;
	fname.append("_");
	found=name.find_last_of("/\\");
	fname.append(name.substr(found+1));
	fname.append(to_string("ERad"));
	fname.append(string(".tsv"));
	ofs.precision(10);
	ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
	for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
		ofs<<tttt[0][simpoint]<<"\t";
		for(int i=0;i<numIons;i++) {
			ofs<<ERadMean[i][simpoint]<<"\t";
		}
		ofs<<endl;
	}
	ofs.close();

	//write mean total energy to temporary file for readout
	fname=runname;
	fname.append("_");
	found=name.find_last_of("/\\");
	fname.append(name.substr(found+1));
	fname.append(to_string("ETotal"));
	fname.append(string(".tsv"));
	ofs.precision(10);
	ofs.open(fname.c_str(),std::ios_base::out|std::ios_base::trunc);
	for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
		ofs<<tttt[0][simpoint]<<"\t";
		for(int i=0;i<numIons;i++) {
			ofs<<ETotMean[i][simpoint]<<"\t";
		}
		ofs<<endl;
	}
	ofs.close();
	

}


EXPORT void SimParallel(int anz){

	useSimSingleSet(false);
	displayiontrap=false;
	//numberThreads=anz;
	runcnt=0;

	cout<<"initializing arrays..."<<endl;
	
	yyyy=makearrayPos<Tpos,3>();cout<<"yyyy"<<endl;
	vvvv=makearrayVel<Tvel,3>();cout<<"vvvv"<<endl;
	tttt=makearrayTime<Ttime,CPUCNT>();cout<<"tttt"<<endl;


	for(int ion=0;ion<numIons;ion++){
		for(int simpoint=0;simpoint<ionplotanz;simpoint++){
			ERadMean[ion][simpoint]=0;
			ETotMean[ion][simpoint]=0;
			for(int tnum=0;tnum<anz;tnum++){
				for(int dim=0;dim<3;dim++){	
					yyyyMean[dim][ion][simpoint]=0;		
					yyyy[dim][tnum][ion][simpoint]=0;
					tttt[tnum][simpoint]=0;					
				}
			}
		}
	}

	cout<<"initializing arrays... done"<<endl;

	//s.init(0); 

	if(anz>CPUCNT) return;
	if(anz<=0) return;
	filenum=0;
	boost::thread *workerThread[CPUCNT];
	for(int i=0;i<anz;i++){
		workerThread[i]=new boost::thread(startsim,i);  
	}
	for(int i=0;i<anz;i++){
		workerThread[i]->join();
	}
	for(int i=0;i<anz;i++){
		delete workerThread[i];
	}
	cout<<"SimParallel fertig"<<endl;
	for(int ion=0;ion<numIons;ion++){
		if(PlotCountPoints[0]>ionplotanz)PlotCountPoints[0]=ionplotanz;
		for(int simpoint=0;simpoint<PlotCountPoints[0];simpoint++){
			ERadMean[ion][simpoint]=0.5*MCa*ERadMean[ion][simpoint]/(anz*loopanz);
			ETotMean[ion][simpoint]=0.5*MCa*ETotMean[ion][simpoint]/(anz*loopanz);
			for(int dim=0;dim<3;dim++){
				//yyyyMean[dim][ion][simpoint]/=runcnt;
				yyyy[dim][0][ion][simpoint]=yyyyMean[dim][ion][simpoint]/(anz*loopanz);
			}
		}
	}
	//hier externes writefile, da writefile ohne Threads nicht läuft
	cout<<"writing...";
	writeAllStuff();
	cout<<"   done"<<endl;
	
	/*cout<<"PlotZ update... ";
	plot1.setData(tttt[0],yyyy[2][0][0]);
	plot1.setNum(PlotCountPoints[0]);
	plot1.update();
	cout<<"   done"<<endl;
	cout<<"PlotX update... ";
	plotX1.setData(tttt[0],yyyy[0][0][0]);
	plotX1.setNum(PlotCountPoints[0]);
	plotX1.update();
	cout<<"   done"<<endl;
	cout<<"PlotY update... ";
	plotY1.setData(tttt[0],yyyy[1][0][0]);
	plotY1.setNum(PlotCountPoints[0]);
	plotY1.update();
	cout<<"   done"<<endl;*/
	cout<<"PlotERad update... ";
	plotERad.setData(tttt[0],ERadMean[0]);
	plotERad.setNum(PlotCountPoints[0]);
	plotERad.update();
	cout<<"   done"<<endl;
	cout<<"PlotZMean update... ";
	plotMeanZ.setData(tttt[0],yyyy[2][0][0]);
	plotMeanZ.setNum(PlotCountPoints[0]);
	plotMeanZ.update();
	cout<<"   done"<<endl;
	cout<<"PlotETotal update... ";
	plotETotal.setData(tttt[0],ETotMean[0]);
	plotETotal.setNum(PlotCountPoints[0]);
	plotETotal.update();
	cout<<"   done"<<endl;
	cout<<"PlotXMean update... ";
	plotMeanX.setData(tttt[0],yyyy[0][0][0]);
	plotMeanX.setNum(PlotCountPoints[0]);
	plotMeanX.update();
	cout<<"   done"<<endl;
	cout<<"PlotYMean update... ";
	plotMeanY.setData(tttt[0],yyyy[1][0][0]);
	plotMeanY.setNum(PlotCountPoints[0]);
	plotMeanY.update();
	cout<<"   done"<<endl;

	cout<<"draw cycle... "<<endl;

	sGiven=true;
	if(dowrite)DrawCycle((PlotCountPoints[0]-100),50);
	sGiven=false;

	cout<<"   done"<<endl;

	delete[] yyyy;
	delete[] vvvv;
	delete[] tttt;
	yyyy=0;
	vvvv=0;
	tttt=0;
}




EXPORT void TestSqueeze(int testLoops){
	/*double OLDangle=angleGet();
	angleSet(0.);
	double OLDendzeit=endzeitGet();
	endzeitSet(0.01);
	int OLDloopanz=loopanzGet();
	loopanzSet(testLoops);*/
	double OLDangle=angle;
	angle=0.;
	int OLDloopanz=loopanz;
	loopanz=testLoops;

	doTestSqueeze=true;
	//cout<<"blooooob"<<endl;
	SimParallel(8);
	doTestSqueeze=false;

	loopanz=OLDloopanz;
	angle=OLDangle;
	for(int i=0;i<CPUCNT;i++)rfMod[i]=1;

	/*angleSet(OLDangle);
	endzeitSet(OLDendzeit);
	loopanzSet(OLDloopanz);*/
}

