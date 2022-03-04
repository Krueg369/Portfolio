#include <stdlib.h>
#include <time.h>
#include <TRandom.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <tgmath.h>
#include <math.h>
#include <cmath>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"
#include "Math/GenVector/RotationZYX.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/BoostY.h"
#include "Math/GenVector/BoostZ.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"
//#include "Vec3.h"
#define INTERVAL 1000000

using namespace std;
bool OctagonAreaCalc_1(Double_t x, Double_t z);
bool OctagonAreaCalc_2(Double_t x, Double_t z);
Double_t distCalc (Double_t x, Double_t y, Double_t z);
Double_t angleCalc(Double_t x, Double_t y, Double_t z, Double_t r);
Double_t angleCalc2(Double_t dx, Double_t dy, Double_t dz, Double_t x, Double_t y, Double_t z);


void Bool(){

// Define Histograms
TH1D *lengthHisto = new TH1D("Legend","Distance of tracks in box; length(m); Frequency",150,0.0,30.0);
TH2D *DistanceHisto = new TH2D("Legend", "Distance from point to yaxis vs. Distance in Box", 100.0,0.0,100.0, 100,0.0,100.0);
TH2D *AngleHisto = new TH2D("Legend", "Angle Vs. Distance", 300.0,0.0,300.0, 300.0,0.0,300.0);
TH2D *angAngHisto = new TH2D("Legend", "Angle Vs. Angle", 32.0,0.0,32.0, 32.0,0.0,32.0);

// Create random seed object
TRandom2 *q = new TRandom2();
q->SetSeed(0);
Double_t x_i,y_i,z_i,r1,r2,dx_i,dy_i,dz_i,boxpnts,sphpnts,detec, dist,spikedist;
//Double_t x4 = 0,x_4 = 0,y15 = 0,y_15 = 0,z4 = 0,z_4 = 0;
//Double_t _xz5 = 0, xz5 = 0, _xz_5 = 0, xz_5 = 0;
Double_t rays = 0, err = 0;
Double_t TestDist, thetaR, thetaV;
Int_t realboxpnts = 0;
bool Test1 = false;
bool Test2 = false;
r1 = 1,0;
r2 = 60.0;

//_______________________Plane definitions for Detector_______________________//

// D's for each plane equation of form Ax+By+Cz = D
Double_t d1 = 4.0,d2 = -4.0,d3 = 15.0,d4 = -15.0,d5 = 4.0,d6 = -4.0;
// Normal vectors to each plane
TVector3 planeNx1(1.0,0.0,0.0);
TVector3 planeNx2(1.0,0.0,0.0);
TVector3 planeNz1(0.0,0.0,1.0);
TVector3 planeNy1(0.0,1.0,0.0);

//octagon plane cuts
Double_t od1 = 5.65685;
Double_t od2 = -5.65685;
TVector3 planeNOtxz1(-1.0,0.0,1.0); //-x + z = 5
TVector3 planeNOtxz2(1.0,0.0,1.0); // x + z = , d = 5.6667

//----------------------------------------------------------------------------//
//____________________Plane definitions for Reference Box_____________________//




//----------------------------------------------------------------------------//

for(Int_t i = 0; i < INTERVAL; i++){

  Double_t x4 = 0,x_4 = 0,y15 = 0,y_15 = 0,z4 = 0,z_4 = 0;
  Double_t _xz5 = 0, xz5 = 0, _xz_5 = 0, xz_5 = 0;
  Int_t n = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0,
  n5 = 0, n6 =0, n7 = 0, n8 = 0, n9 = 0, n10 = 0;
  TVector3 Temp1(0.0,0.0,0.0);
  TVector3 Temp2(0.0,0.0,0.0);
  Int_t RayTraj[10] = {0,0,0,0,0,0,0,0,0,0};
  //cout << "Hello" << endl;
  rays++;
  //generates directional cosines and starting points on sphere

  q->Sphere(dx_i,dy_i,dz_i,r1);
  q->Sphere(x_i,y_i,z_i,r2);


  // put the starting points and directions into vectors for each iteration
  TVector3 d_i;
  TVector3 p_i;
  TVector3 A(0.0, -40.0, 0.0);
  TVector3 B(0.0, 40.0, 0.0);
  d_i.SetXYZ(dx_i,dy_i,dz_i);
  p_i.SetXYZ(x_i,y_i,z_i);
  //cout << "INITIAL POINT: " << x_i << ", " << y_i << ", " << z_i << endl;
  //Double_t Magn = sqrt(x_i * x_i + y_i*y_i + z_i*z_i);
  //cout << "INITIAL Magnitude: " << Magn << endl;
  // Denomenator to calculate parameter t
  Double_t denom1 = planeNx1.Dot(d_i); // x = -4
  Double_t denom2 = planeNx2.Dot(d_i);// x = 4
  Double_t denom5 = planeNz1.Dot(d_i);
  Double_t denom6 = planeNy1.Dot(d_i);
  Double_t denom7 = planeNOtxz1.Dot(d_i);
  Double_t denom8 = planeNOtxz2.Dot(d_i);

  //if statements to avoid dividing by zero
  if(denom1 == 0) break;;
  if(denom2 == 0) break;;
  if(denom5 == 0) break;;
  if(denom6 == 0) break;;
  if(denom7 == 0) break;;
  if(denom8 == 0) break;;
  //calculate the parameter t using normal to plane, starting point, D and denom
  Double_t t1_i = -(planeNx1.Dot(p_i) - d2) / (planeNx1.Dot(d_i));// x = -4
  Double_t t2_i = -(planeNx2.Dot(p_i) - d1) / (planeNx2.Dot(d_i));// x = 4
  Double_t t3_i = -(planeNy1.Dot(p_i) - d3) / (planeNy1.Dot(d_i));// x = 4
  Double_t t4_i = -(planeNy1.Dot(p_i) - d4) / (planeNy1.Dot(d_i));// x = 4
  Double_t t5_i = -(planeNz1.Dot(p_i) - d5) / (planeNz1.Dot(d_i));
  Double_t t6_i = -(planeNz1.Dot(p_i) - d6) / (planeNz1.Dot(d_i));
  Double_t t7_i = -(planeNOtxz1.Dot(p_i) - od1) / (planeNOtxz1.Dot(d_i));
  Double_t t8_i = -(planeNOtxz2.Dot(p_i) - od1) / (planeNOtxz2.Dot(d_i));
  Double_t t9_i = -(planeNOtxz1.Dot(p_i) - od2) / (planeNOtxz1.Dot(d_i));
  Double_t t10_i = -(planeNOtxz2.Dot(p_i) - od2) / (planeNOtxz2.Dot(d_i));

  //cout << "Paramteres: " << t1_i << ", " << t2_i << t3_i << ", " << t4_i << t5_i << ", " << t6_i << ", " <<  t7_i << endl;

  // Dot product between normal vector and line dir vector to test if perp.
  // if perp. then line is paralell to plane, means no hit on that plane
  if(planeNx1.Dot(d_i) == 0) cout << "MISS" << endl;
  if(planeNx2.Dot(d_i) == 0) cout << "MISS" << endl;
  if(planeNz1.Dot(d_i) == 0) cout << "MISS" << endl;
  if(planeNOtxz1.Dot(d_i) == 0) cout << "MISS" << endl;
  if(planeNOtxz2.Dot(d_i) == 0) cout << "MISS" << endl;
  // calculate intersection with starting point vector, parameter t and dir. vector

  TVector3 intersection1_i = p_i + t1_i * d_i;
  TVector3 intersection2_i = p_i + t2_i * d_i;
  TVector3 intersection3_i = p_i + t3_i * d_i;
  TVector3 intersection4_i = p_i + t4_i * d_i;
  TVector3 intersection5_i = p_i + t5_i * d_i;
  TVector3 intersection6_i = p_i + t6_i * d_i;
  TVector3 intersection7_i = p_i + t7_i * d_i;
  TVector3 intersection8_i = p_i + t8_i * d_i;
  TVector3 intersection9_i = p_i + t9_i * d_i;
  TVector3 intersection10_i = p_i + t10_i * d_i;



  //pull each compenent out of point of intersection, test if it is in range
  // of each plane

  Double_t xx1_i = intersection1_i.X();
  Double_t yy1_i = intersection1_i.Y();
  Double_t zz1_i = intersection1_i.Z();

  Double_t xx2_i = intersection2_i.X();
  Double_t yy2_i = intersection2_i.Y();
  Double_t zz2_i = intersection2_i.Z();

  Double_t xx3_i = intersection3_i.X();
  Double_t yy3_i = intersection3_i.Y();
  Double_t zz3_i = intersection3_i.Z();

  Double_t xx4_i = intersection4_i.X();
  Double_t yy4_i = intersection4_i.Y();
  Double_t zz4_i = intersection4_i.Z();

  Double_t xx5_i = intersection5_i.X();
  Double_t yy5_i = intersection5_i.Y();
  Double_t zz5_i = intersection5_i.Z();

  Double_t xx6_i = intersection6_i.X();
  Double_t yy6_i = intersection6_i.Y();
  Double_t zz6_i = intersection6_i.Z();

  Double_t xx7_i = intersection7_i.X();
  Double_t yy7_i = intersection7_i.Y();
  Double_t zz7_i = intersection7_i.Z();

  Double_t xx8_i = intersection8_i.X();
  Double_t yy8_i = intersection8_i.Y();
  Double_t zz8_i = intersection8_i.Z();

  Double_t xx9_i = intersection9_i.X();
  Double_t yy9_i = intersection9_i.Y();
  Double_t zz9_i = intersection9_i.Z();

  Double_t xx10_i = intersection10_i.X();
  Double_t yy10_i = intersection10_i.Y();
  Double_t zz10_i = intersection10_i.Z();
/*
| Determine if intersection is in region of each face, if so
| up the number of elements in RayTraj array and give value of 1 to elements
| This way if the number of elements in the array is anythig but 2 there is
|  a mistake
*/




  //x = - 4 plane
  if((yy1_i >= -15.0)&&(yy1_i <= 15.0)&&(zz1_i >= -4.0)&&(zz1_i <= 1.65685)){
    n1++;
    RayTraj[0] = 1;
    boxpnts++;
  }

  //x = 4 plane
  if((yy2_i >= -15.0)&&(yy2_i <= 15.0)&&(zz2_i >= -4.0)&&(zz2_i <= 1.65685)){
    n2++;
    RayTraj[1] = 1;
    boxpnts++;
  }

  //y = 15
  if(xx3_i <= 4.0 && xx3_i >= -4.0 && zz3_i <= 4.0 && zz3_i >= -4.0 ){
    Test1 = OctagonAreaCalc_1(xx3_i,zz3_i);
  if(Test1 == true){
    n3++;
    RayTraj[2] = 1;
    boxpnts++;
  }
}
  //y = -15
  if(xx4_i <= 4.0 && xx4_i >= -4.0 && zz4_i <= 4.0 && zz4_i >= -4.0 ){

    Test2 = OctagonAreaCalc_2(xx4_i, zz4_i);
    if(Test2 == true){
      n4++;
      RayTraj[3] = 1;
      boxpnts++;
    }
  }

  //z = 4
  if((xx5_i <= 1.65685)&&(xx5_i >= -1.65685)&&(yy5_i >= -15.0)&&(yy5_i <= 15.0)){
    n5++;
    RayTraj[4] = 1;
    boxpnts++;
  }

  // z = -4
  if((xx6_i <= 1.65685)&&(xx6_i >= -4.0)&&(yy6_i >= -15.0)&&(yy6_i <= 15.0)){
    n6++;
    RayTraj[5] = 1;
    boxpnts++;
  }


  // -x + z = 5.65685
  if((yy7_i >= -15.0)&&(yy7_i <= 15.0) && (xx7_i >= -4.0) && (xx7_i <= -1.65685) && (zz7_i >= 1.65685) && (zz7_i <= 4.0)){
    if(zz7_i > 0 && xx7_i < 0){
      n7++;
      RayTraj[6] = 1;
      boxpnts++;
    }
  }

  //x + z = 5
  if((yy8_i >= -15.0)&&(yy8_i <= 15.0)&&(xx8_i >= 1.65685)&&(xx8_i <= 4.0)&&(zz8_i >= 1.65685)&&(zz8_i <= 4.0)){
    n8++;
    RayTraj[7] = 1;
    boxpnts++;
  }

  // -x + z = -5.6667
  if((yy9_i >= -15.0)&&(yy9_i <= 15.0)&&(xx9_i >= 1.65685)&&(xx9_i <= 4.0)&&(zz9_i >= -4.0)&&(zz9_i <= -1.65685)){
    n9++;
    RayTraj[8] = 1;
    boxpnts++;
  }

  //x + z = -6.6667
  if((yy10_i >= -15.0)&&(yy10_i <= 15.0)&&(xx10_i >= -4.0)&&(xx10_i <= -1.65685)&&(zz10_i >= -4.0)&&(zz10_i <= -1.65685)){
    n10++;
    RayTraj[9] = 1;
    boxpnts++;
  }
/*
if(n3 == 1){
  cout << xx3_i << "," << yy3_i << ", " << zz3_i << endl;
}
*/

  n = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9 + n10;

  if(n == 1){
    err++;
  //cout << n1 << ", " << n2 << ", " << n3 << ", " << n4 << ", " << n5 << ", " << n6 << ", " << n7 << ", " << n8 << ", " << n9  << ", " << n10 << ", " << endl;
  }
/*
if(n10 == 1){
  cout << "n10: "<< xx10_i << ", " << yy10_i << ", " << zz10_i << "," << xx10_i + zz10_i <<  endl;
}
*/

  Int_t normal1 = 0;
  Int_t normal2 = 0;
/*
if(-xx7_i + zz7_i == 5.65685){
  normal1++;
}
if(-xx7_i + zz7_i == -5.65685){
  normal2++;
}
if(-xx7_i + zz7_i == -5.65685){cout << "ERROR ERRORR ERRORR" << endl;}



if(normal1 > 0) cout << "-x + z = +5: " << normal1 << endl;
if(normal2 > 0) cout << "-x + z = -5: " << normal2 << endl;
*/
/*
2 for loops, first to fill a vector with first ray that intersects
once a ray takes the first vector the loop is broken
second for loop will fill the second vector with next intersection
*/
  if(n == 2){
    realboxpnts++;
    for(Int_t j = 0; j < 1; j++){
      if(RayTraj[0] == 1){
        Temp1 = intersection1_i;
        x_4++;
        break;;

      }
      if(RayTraj[1] == 1){
        Temp1 = intersection2_i;
        x4++;
        break;;

      }
      if(RayTraj[2] ==1){
        Temp1 = intersection3_i;
        y15++;
        break;;

      }
      if(RayTraj[3] == 1){
        Temp1 = intersection4_i;
        y_15++;
        break;;
      }

      if(RayTraj[4] == 1){
        Temp1 = intersection5_i;
        z4++;
        break;;
      }
      if(RayTraj[5] == 1){
        Temp1 = intersection6_i;
        z_4++;
        break;;
      }
      if(RayTraj[6] == 1){
        Temp1 = intersection7_i;
        _xz5++;
        break;;
      }
      if(RayTraj[7] == 1){
        Temp1 = intersection8_i;
        xz5++;
        break;;
      }
      if(RayTraj[8] == 1){
        Temp1 = intersection9_i;
        _xz_5++;
        break;;
      }
      if(RayTraj[9] == 1){
        Temp1 = intersection10_i;
        xz_5++;
        break;;
      }
//cout << "LOOP 1" << endl;


  }////

    for(Int_t h = 0; h < 1; h++){
      if((RayTraj[0] == 1)&&(x_4 == 0)){
        Temp2 = intersection1_i;
        break;;
      }
      if((RayTraj[1] == 1)&&(x4 == 0)){
        Temp2 = intersection2_i;
        break;;
      }
      if((RayTraj[2] == 1)&&(y15 == 0)){
        Temp2 = intersection3_i;
        break;;
      }
      if((RayTraj[3] == 1)&&(y_15 == 0)){
        Temp2 = intersection4_i;
        break;;
      }
      if((RayTraj[4] == 1)&&(z4 == 0)){
        Temp2 = intersection5_i;
        break;;
      }
      if((RayTraj[5] == 1)&&(z_4 == 0)){
        Temp2 = intersection6_i;
        break;;
      }

      if((RayTraj[6] == 1)&&(_xz5 == 0)){
        Temp2 = intersection7_i;
        break;;
      }
      if((RayTraj[7] == 1)&&(xz5 == 0)){
        Temp2 = intersection8_i;
        break;;
      }
      if((RayTraj[8] == 1)&&(_xz_5 == 0)){
        Temp2 = intersection9_i;
        break;;
      }
      if((RayTraj[9] == 1)&&(xz_5 == 0)){
        Temp2 = intersection10_i;
        break;;
      }

    }

}

  // Calculate the distance, pull out x,y,z from temp vectors and use 3D pythagorean

  Double_t temp1x_i = Temp1.X();
  Double_t temp1y_i = Temp1.Y();
  Double_t temp1z_i = Temp1.Z();

  Double_t temp2x_i = Temp2.X();
  Double_t temp2y_i = Temp2.Y();
  Double_t temp2z_i = Temp2.Z();

  /*
  if(n == 1 && temp1x_i == 0.0 && temp2x_i == 0.0){
    cout << "ERROR" << endl;
  }
*/
/*
if(n == 1){
cout << "First hit: " << temp1x_i << ", " << temp1y_i << ", " << temp1z_i << endl;
cout << "Second hit: " << temp2x_i << ", " << temp2y_i << ", " << temp2z_i << endl;
}
*/
/*
if((x4 != 0) && ())
cout << "X = 4: " << x4 << endl;
cout << "X = -4: " << x_4 << endl;
cout << "Y = 15: " << y15 << endl;
cout << "Y = -15: " << y_15 << endl;
cout << "Z = 4: " << z4 << endl;
cout << "Z = -4: " << z_4 << endl;
cout << "_xZ = 5.65685: " << _xzO << endl;
*/
  Double_t distance = sqrt((pow((temp2x_i - temp1x_i),2)) + (pow((temp2y_i - temp1y_i),2)) + (pow((temp2z_i - temp1z_i),2)));
//TH2D *angAngHisto = new TH2D("Legend", "Angle Vs. Angle", distance,0.0,32.0, distance,0.0,32.0);

  TestDist = distCalc(x_i, y_i, z_i);
  thetaR = angleCalc(x_i, y_i, z_i, r2);
  thetaV = angleCalc2(dx_i, dy_i, dz_i, x_i, y_i, z_i);
  //cout << thetaR << endl;
  if(distance != 0.0){
    dist++;
    lengthHisto-> Fill(distance);
    DistanceHisto->Fill(TestDist,distance);
    AngleHisto->Fill(thetaR, distance);
    angAngHisto->Fill(thetaV, thetaR);
  }




} //END FOR LOOP

//Int_t realboxpnts = boxpnts / 2;
//make sure number of rays and number of distances match
cout <<  realboxpnts << " Rays Hit the Detector" << endl;
cout << "Number of distances: " << dist << endl;
cout << rays << endl;
cout << "Number of errors" << err << endl;
TCanvas* lengCanvas = new TCanvas("distCanvas","distance Histogram",640,480);
lengthHisto->Draw();

TCanvas* lengthCanvas = new TCanvas("testdistCanvas","test distance Histogram",640,480);
DistanceHisto->Draw("SURF3Z");

TCanvas* angleCanvas = new TCanvas("angleCanvas","angle Vs. Distance",640,480);
AngleHisto->Draw("SURF3Z");

TCanvas* angleAngleCanvas = new TCanvas("angleAngleCanvas","angle Vs. angle",640,480);
angAngHisto->Draw("SURF3Z");

}



//----------------------------------------------------------------------------//
//_______________________________Functions____________________________________//

// Area of octagon for y = 15 and y = -15
bool OctagonAreaCalc_1(Double_t x, Double_t z){
  Double_t areaCALC = 0;
  //bool y15FACE;

  areaCALC =
  //APB
  abs((1.65685*z + x*4.0 + -1.65685*4.0 - 4.0*x - z*-1.65685 - 4.0*1.65685)/2)+
  //BPC
  abs((-1.65685*z + x*1.65685 + -4.0*4.0 - (-4.0)*x - z*-4.0 - 1.65685*-1.65686)/2)+
  //CPD
  abs((-4.0*z + x*-1.65685 + -4.0*1.65685 - 1.65685*x - z*-4.0 - (-1.65685)*(-4.0))/2)+
  //DPE
  abs((-4.0*z + x*-4.0 + -1.65685*-1.65685 - (-1.65685)*x - z*-1.65685 - (-4.0)*(-4.0))/2)+
  //EPF
  abs((-1.65685*z + x*-4.0 + 1.65685*-4.0- (-4.0)*x - z*1.65685 - (-4.0)*(-1.65685))/2)+
  //FPG
  abs((1.65685*z + x*-1.65685 + 4.0*(-4.0) - (-4.0)*x - z*4.0 - (-1.65685)*(1.65685))/2)+
  //GPH
  abs((4.0*z + x*1.65685 + 4.0*(-1.65685) - (-1.65685)*x - z*4.0 - 1.65685*(4.0))/2)+
  //HPA
  abs((4.0*z + x*4.0 + 1.65685*1.65685 - 1.65685*x - z*1.65685 - 4.0*4.0)/2);

  if(areaCALC <= 52.901 && areaCALC > 0){
    return true;
  }
  else return false;
}

bool OctagonAreaCalc_2(Double_t x, Double_t z){
  Double_t areaCALC = 0;
  //bool y15FACE;

  areaCALC =
  //APB
  abs((1.65685*z + x*4.0 + -1.65685*4.0 - 4.0*x - z*-1.65685 - 4.0*1.65685)/2)+
  //BPC
  abs((-1.65685*z + x*1.65685 + -4.0*4.0 - (-4.0)*x - z*-4.0 - 1.65685*-1.65686)/2)+
  //CPD
  abs((-4.0*z + x*-1.65685 + -4.0*1.65685 - 1.65685*x - z*-4.0 - (-1.65685)*(-4.0))/2)+
  //DPE
  abs((-4.0*z + x*-4.0 + -1.65685*-1.65685 - (-1.65685)*x - z*-1.65685 - (-4.0)*(-4.0))/2)+
  //EPF
  abs((-1.65685*z + x*-4.0 + 1.65685*-4.0- (-4.0)*x - z*1.65685 - (-4.0)*(-1.65685))/2)+
  //FPG
  abs((1.65685*z + x*-1.65685 + 4.0*(-4.0) - (-4.0)*x - z*4.0 - (-1.65685)*(1.65685))/2)+
  //GPH
  abs((4.0*z + x*1.65685 + 4.0*(-1.65685) - (-1.65685)*x - z*4.0 - 1.65685*(4.0))/2)+
  //HPA
  abs((4.0*z + x*4.0 + 1.65685*1.65685 - 1.65685*x - z*1.65685 - 4.0*4.0)/2);

  if(areaCALC <= 52.901 && areaCALC > 0){
    return true;
  }
  else return false;
}

Double_t distCalc(Double_t x, Double_t y, Double_t z){

  if( y >= 0.0){
    Double_t x1 = 0.0;
    Double_t y1 = 40.0;
    Double_t z1 = 0.0;


  Double_t testDist = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1) + (z - z1)*(z - z1));
  return testDist;
  }

  if( y < 0.0 ){
    Double_t x2 = 0.0;
    Double_t y2 = -40.0;
    Double_t z2 = 0.0;

    Double_t testDist = sqrt((x - x2)*(x - x2)+ (y - y2)*(y - y2) + (z - z2)*(z - z2));
    return testDist;
  }
}

Double_t angleCalc(Double_t x, Double_t y, Double_t z, Double_t r){
  if( y < 0.0){
    Double_t thetaR = acos((y * (-r))/(sqrt(x*x + y*y + z*z)*(-r))) * 100;
    return thetaR;
  }

  if( y >= 0.0){
    Double_t thetaR = acos((y * (r))/(sqrt(x*x + y*y + z*z)*(r))) * 100;
    return thetaR;
  }
}

Double_t angleCalc2(Double_t dx, Double_t dy, Double_t dz, Double_t x, Double_t y, Double_t z){
  Double_t MagP = sqrt(x*x + y*y + z*z);
  Double_t MagD = sqrt(dx*dx + dy*dy + dz*dz);
  Double_t thetaV = acos((dx * x + dy * y + dz * z)/(MagP * MagD)) * 100;
  return thetaV;
}
