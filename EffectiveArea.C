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
#include "TAttMarker.h"
#include <iostream>
//#include "Vec3.h"
#define INTERVAL 10000


using namespace std;
bool OctagonAreaCalc(Double_t x, Double_t z);
Double_t distCalc (Double_t x, Double_t y, Double_t z);
Double_t angleCalc(Double_t x, Double_t y, Double_t z, Double_t r);
Double_t angleCalc2(Double_t dx, Double_t dy, Double_t dz, Double_t x, Double_t y, Double_t z);
Double_t RefFracCalc(Double_t x, Double_t y, Double_t z);


void testTest(){

  // Define Histograms
  TH1D *lengthHisto = new TH1D("Legend","Distance of tracks in box; length(m); Frequency",100,0.0,35.0);
  TH2D *DistanceHisto = new TH2D("Legend", "Distance in box vs. Distance from Y-Axis", 100.0,0.0,100.0, 100,0.0,100.0);
  TH2D *AngleHisto = new TH2D("Legend", "Angle Vs. Distance", 300.0,0.0,300.0, 300.0,0.0,300.0);
  TH2D *angAngHisto = new TH2D("Legend", "Angle Vs. Angle", 32.0,0.0,32.0, 32.0,0.0,32.0);
  TH2D *AreaVSAngleHisto = new TH2D("Legend", "Area as a Function of Theta, Phi", 200.0,-100.0,100.0, 200.0,0.0,200.0);
  TH2D *DivisionHisto= new TH2D("Legend", "Area Vs Theta", 200.0, -100.0,100.0, 200.0,0.0,200.0);
  TH2D *HitsHistogram = new TH2D("Legend", "How many hit planes at each theta,phi", 200.0, -100.0, 100.0, 200.0, 0.0, 200.0);
  TH2D *EffectiveArea = new TH2D();
  TH1D *AreavsTheta = new TH1D("Legend", "Area vs Theta", 200.0, -100.0, 100.0);
  TH1D *AreavsPhi = new TH1D("Legend", "Area vs Phi", 200.0, 0.0, 200.0);
  EffectiveArea->Multiply(AreaVSAngleHisto,DivisionHisto);

  TGraph2D *Points = new TGraph2D();
  TGraph2D *Points2 = new TGraph2D();
  TGraph2D *Points3 = new TGraph2D();
  TGraph2D *Points4 = new TGraph2D();


  // Create random seed object
  TRandom2 *q = new TRandom2();
  q->SetSeed(0);
  Double_t x_i,y_i,z_i,r1,r2,dx_i,dy_i,dz_i,boxpnts = 0, dist;
  Double_t rays = 0, err = 0;
  Double_t TestDist, thetaR, thetaV;
  Int_t realboxpnts = 0;
  Double_t refpnts = 0;

  // Test variables for  octagon face area, not currently being used
  bool Test1 = false;
  bool Test2 = false;
  r1 = 1.0;
  r2 = 40.0;

  Double_t r, phi, theta;
  Int_t iteration = 500;

  //_______________________Plane definitions for Detector_______________________//

  // D's for each plane equation of form Ax+By+Cz = D
  Double_t d1 = 4.0,d2 = -4.0,d3 = 15.0,d4 = -15.0;
  Double_t dB1 = 8.0, dB2 = -8.0, dB3 = 20.0, dB4 = -20.0;

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

  // Reference points for plane calculations
  TVector3 OPoint(0.01, 0.0, 0.01);
  TVector3 OPoint1(-0.01, 0.0, 0.01);
  TVector3 OPoint2(-0.01, 0.0, -0.01);
  TVector3 OPoint3(0.01, 0.0, -0.01);

  // points to calculate parameter t
  TVector3 Point1(0.0, 0.0, -0.01);
  TVector3 Point2(0.0, 0.0, 0.01);
  TVector3 Point3(0.01, 0.0, 0.0);
  TVector3 Point4(-0.01,0.0,0.0);
  TVector3 Point5(0.0, 0.01, 0.0);
  TVector3 Point6(0.0, -0.01, 0.0);

  //----------------------------------------------------------------------------//
  //____________________Plane definitions for Reference Box_____________________//

  TVector3 planeRBx(1.0,0.0,0.0);
  TVector3 planeRBy(0.0,1.0,0.0);
  TVector3 planeRBz(0.0,0.0,1.0);

  // Temporary Area Vectors to be used for later Calculation
  TVector3 TempAreaVec1(0.0,0.0,0.0);
  TVector3 TempAreaVec2(0.0,0.0,0.0);
  TVector3 TempAreaVec3(0.0,0.0,0.0);

  Double_t SumArea = 0;
  Double_t AvgArea = 0;
  Int_t count = 0;
  Int_t count2 = 0;
  Int_t count3 = 0;
  Int_t Visible = 0;
  Double_t Fraction = 0;

  //----------------------------------------------------------------------------//

  for(Int_t h = 0; h < INTERVAL; h++){

    //generates directional cosines and starting points on sphere
    q->Sphere(x_i,y_i,z_i,r2);
    Points3->SetPoint(count2, x_i, y_i, z_i);
    count2++;

    //calculate polar coordinates to find theta and phi of each point
    r = sqrt(x_i*x_i + y_i * y_i + z_i * z_i);
    theta = atan(y_i/x_i) * 180/M_PI;
    phi = acos(z_i/r) * 180/M_PI;

    //____________________________________________________
    Double_t VisibleArea = RefFracCalc(x_i, y_i, z_i);

    if (VisibleArea > 64){
      SumArea = VisibleArea + SumArea;
      AreaVSAngleHisto->Fill(theta, phi, VisibleArea);
      AreavsTheta -> Fill(theta, VisibleArea);
      AreavsPhi -> Fill(phi, VisibleArea);
      Visible++;
    }

    //____________________________________________________

    for(Int_t i = 0; i < iteration; i++){

      Double_t RBx8 = 0, RBx_8 = 0, RBy20 = 0, RBy_20 = 0, RBz8 = 0, RBz_8 = 0;
      //put phi and theta for each itteration into an array

      q->Sphere(dx_i,dy_i,dz_i,r1);
      Points4->SetPoint(count3, dx_i, dy_i,dz_i);
      Double_t x4 = 0,x_4 = 0,y15 = 0,y_15 = 0,z4 = 0,z_4 = 0;
      Double_t _xz5 = 0, xz5 = 0, _xz_5 = 0, xz_5 = 0;
      Int_t n = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0,
      n5 = 0, n6 =0, n7 = 0, n8 = 0, n9 = 0, n10 = 0,
      n11 = 0, n12 = 0, n13 = 0, n14 = 0, n15 = 0, n16 = 0;

      TVector3 Temp1(0.0,0.0,0.0);
      TVector3 Temp2(0.0,0.0,0.0);
      TVector3 RefTemp1(0.0,0.0,0.0);
      TVector3 RefTemp2(0.0,0.0,0.0);

      Int_t RayTraj[10] = {0,0,0,0,0,0,0,0,0,0};
      Int_t RefTraj[6] = {0,0,0,0,0,0};
      rays++;

      // put the starting points and directions into vectors for each iteration
      TVector3 d_i;
      TVector3 p_i;
      TVector3 A(0.0, -40.0, 0.0);
      TVector3 B(0.0, 40.0, 0.0);

      if(dx_i != 0 && dy_i !=0 && dz_i != 0){
        d_i.SetXYZ(dx_i,dy_i,dz_i);
      }
      if(x_i != 0 && y_i != 0 && z_i != 0){
        p_i.SetXYZ(x_i,y_i,z_i);
      }

      // Denomenator to calculate parameter t
      Double_t denom1 = planeNx1.Dot(d_i); // x = -4
      Double_t denom2 = planeNx2.Dot(d_i);// x = 4
      Double_t denom5 = planeNz1.Dot(d_i);
      Double_t denom6 = planeNy1.Dot(d_i);
      Double_t denom7 = planeNOtxz1.Dot(d_i);
      Double_t denom8 = planeNOtxz2.Dot(d_i);
      Double_t denom9 = planeRBx.Dot(d_i);
      Double_t denom10 = planeRBy.Dot(d_i);
      Double_t denom11 = planeRBz.Dot(d_i);

      //if statements to avoid dividing by zero
      if(denom1 == 0) break;;
      if(denom2 == 0) break;;
      if(denom5 == 0) break;;
      if(denom6 == 0) break;;
      if(denom7 == 0) break;;
      if(denom8 == 0) break;;
      if(denom9 == 0) break;;
      if(denom10 == 0) break;;
      if(denom11 == 0) break;;

      // OLD METHODS COMMENTED OUT DUE TO BUGS
      //calculate the parameter t using normal to plane, starting point, D and denom
      /*
      Double_t t1_i = -(planeNx1.Dot(p_i) - d2) / (planeNx1.Dot(d_i));// x = -4
      Double_t t2_i = -(planeNx2.Dot(p_i) - d1) / (planeNx2.Dot(d_i));// x = 4
      Double_t t3_i = -(planeNy1.Dot(p_i) - d3) / (planeNy1.Dot(d_i));//
      Double_t t4_i = -(planeNy1.Dot(p_i) - d4) / (planeNy1.Dot(d_i));//
      Double_t t5_i = -(planeNz1.Dot(p_i) - d1) / (planeNz1.Dot(d_i));
      Double_t t6_i = -(planeNz1.Dot(p_i) - d2) / (planeNz1.Dot(d_i));
      */

      Double_t t1_i = ((Point4 - p_i).Dot(planeNx1) - d1) / denom1;
      Double_t t2_i = ((Point3 - p_i).Dot(planeNx1) + d1) / denom2;
      Double_t t3_i = ((Point5 - p_i).Dot(planeNy1) + d3) / denom6;
      Double_t t4_i = ((Point6 - p_i).Dot(planeNy1) - d3) / denom6;
      Double_t t5_i = ((Point2 - p_i).Dot(planeNz1) + d1) / denom5;
      Double_t t6_i = ((Point1 - p_i).Dot(planeNz1) - d1) / denom5;

      /*
      Double_t t7_i = -(planeNOtxz1.Dot(p_i) - od1) / (planeNOtxz1.Dot(d_i)); //top left
      Double_t t8_i = -(planeNOtxz2.Dot(p_i) - od1) / (planeNOtxz2.Dot(d_i));
      Double_t t9_i = -(planeNOtxz1.Dot(p_i) - od2) / (planeNOtxz1.Dot(d_i)); // top right
      Double_t t10_i = -(planeNOtxz2.Dot(p_i) - od2) / (planeNOtxz2.Dot(d_i));
      */

      Double_t t7_i = ((OPoint1 - p_i).Dot(planeNOtxz1) + od1) / denom7;
      Double_t t8_i = ((OPoint2 - p_i).Dot(planeNOtxz2) - od1) / denom8;
      Double_t t9_i = ((OPoint3 - p_i).Dot(planeNOtxz1) - od1) / denom7;
      Double_t t10_i = ((OPoint - p_i).Dot(planeNOtxz2) + od1) / denom8;

      //Parameteres for reference box
      Double_t t11_i = -(planeRBx.Dot(p_i) - dB1) / (planeRBx.Dot(d_i));
      Double_t t12_i = -(planeRBx.Dot(p_i) - dB2) / (planeRBx.Dot(d_i));
      Double_t t13_i = -(planeRBy.Dot(p_i) - dB3) / (planeRBy.Dot(d_i));
      Double_t t14_i = -(planeRBy.Dot(p_i) - dB4) / (planeRBy.Dot(d_i));
      Double_t t15_i = -(planeRBz.Dot(p_i) - dB1) / (planeRBz.Dot(d_i));
      Double_t t16_i = -(planeRBz.Dot(p_i) - dB2) / (planeRBz.Dot(d_i));

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
      TVector3 intersection11_i = p_i + t11_i * d_i;
      TVector3 intersection12_i = p_i + t12_i * d_i;
      TVector3 intersection13_i = p_i + t13_i * d_i;
      TVector3 intersection14_i = p_i + t14_i * d_i;
      TVector3 intersection15_i = p_i + t15_i * d_i;
      TVector3 intersection16_i = p_i + t16_i * d_i;

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

      Double_t xx11_i = intersection11_i.X();
      Double_t yy11_i = intersection11_i.Y();
      Double_t zz11_i = intersection11_i.Z();

      Double_t xx12_i = intersection12_i.X();
      Double_t yy12_i = intersection12_i.Y();
      Double_t zz12_i = intersection12_i.Z();

      Double_t xx13_i = intersection13_i.X();
      Double_t yy13_i = intersection13_i.Y();
      Double_t zz13_i = intersection13_i.Z();

      Double_t xx14_i = intersection14_i.X();
      Double_t yy14_i = intersection14_i.Y();
      Double_t zz14_i = intersection14_i.Z();

      Double_t xx15_i = intersection15_i.X();
      Double_t yy15_i = intersection15_i.Y();
      Double_t zz15_i = intersection15_i.Z();

      Double_t xx16_i = intersection16_i.X();
      Double_t yy16_i = intersection16_i.Y();
      Double_t zz16_i = intersection16_i.Z();
/*
| Determine if intersection is in region of each face, if so
| up the number of elements in RayTraj array and give value of 1 to elements
| This way if the number of elements in the array is anythig but 2 there is
|  a mistake
*/

      //x = - 4 plane
      if((yy1_i >= -15.0)&&(yy1_i <= 15.0)&&(zz1_i >= -1.65685)&&(zz1_i <= 1.65685)){
        n1 = 1;
        RayTraj[0] = 1;
        boxpnts++;
      }

      //x = 4 plane
      if((yy2_i >= -15.0)&&(yy2_i <= 15.0)&&(zz2_i >= -1.65685)&&(zz2_i <= 1.65685)){
        n2 = 1;
        RayTraj[1] = 1;
        boxpnts++;
      }

      //y = 15
      Test1 = OctagonAreaCalc(xx3_i,zz3_i);
      if(xx3_i <= 4.0 && xx3_i >= -4.0 && zz3_i <= 4.0 && zz3_i >= -4.0 &&
      (xx3_i + zz3_i <= 5.65686) && (-xx3_i + zz3_i <= 5.65686) && (-xx3_i + zz3_i >= -5.65686) && (xx3_i + zz3_i >= -5.65686)){
        n3 = 1;
        RayTraj[2] = 1;
        boxpnts++;
      }
      //}
      //y = -15
      Test2 = OctagonAreaCalc(xx4_i, zz4_i);
      if(xx4_i <= 4.0 && xx4_i >= -4.0 && zz4_i <= 4.0 && zz4_i >= -4.0 &&
      (xx3_i + zz3_i <= 5.65686) && (-xx3_i + zz3_i <= 5.65686) && (-xx3_i + zz3_i >= -5.65686) && (xx3_i + zz3_i >= -5.65686)){
        n4 = 1;
        RayTraj[3] = 1;
        boxpnts++;
      }
      //}

      //z = 4
      if((xx5_i <= 1.65685)&&(xx5_i >= -1.65685)&&(yy5_i >= -15.0)&&(yy5_i <= 15.0)){
        n5 = 1;
        RayTraj[4] = 1;
        boxpnts++;
      }

      // z = -4
      if((xx6_i <= 1.65685)&&(xx6_i >= -1.65685)&&(yy6_i >= -15.0)&&(yy6_i <= 15.0)){
        n6 = 1;
        RayTraj[5] = 1;
        boxpnts++;
      }


      // -x + z = 5.65685
      if((yy7_i >= -15.0)&&(yy7_i <= 15.0) && (xx7_i >= -4.0) && (xx7_i <= -1.65685) && (zz7_i >= 1.65685) && (zz7_i <= 4.0)){

          n7 = 1;
          RayTraj[6] = 1;
          boxpnts++;

      }

      //x + z = -5
      if((yy8_i >= -15.0)&&(yy8_i <= 15.0)&&(xx8_i <= -1.65685)&&(xx8_i >= -4.0)&&(zz8_i <= -1.65685)&&(zz8_i >= -4.0)){
        n8 = 1;
        RayTraj[7] = 1;
        boxpnts++;
      }



      // -x + z = -5.6667
      if((yy9_i >= -15.0)&&(yy9_i <= 15.0)&&(xx9_i >= 1.65685)&&(xx9_i <= 4.0)&&(zz9_i >= -4.0)&&(zz9_i <= -1.65685)){
        n9 = 1;
        RayTraj[8] = 1;
        boxpnts++;
      }

      //x + z = -6.6667
      if((yy10_i >= -15.0)&&(yy10_i <= 15.0)&&(xx10_i <= 4.0)&&(xx10_i >= 1.65685)&&(zz10_i <= 4.0)&&(zz10_i >= 1.65685)){
        n10 = 1;
        RayTraj[9] = 1;
        boxpnts++;
      }

  //BOUNDARIES FOR REFERENCE BOX_______________________________________________________________________________________________________________

      // x = 8
      if((yy11_i >= -20.0)&&(yy11_i <= 20.0)&&(zz11_i >= -8.0)&&(zz11_i <= 8.0)){
        n11++;
        RefTraj[0] = 1;
        refpnts++;
      }
      // x = -8
      if((yy12_i >= -20.0)&&(yy12_i <= 20.0)&&(zz12_i >= -8.0)&&(zz12_i <= 8.0)){
        n12++;
        RefTraj[1] = 1;
        refpnts++;
      }

      // y = 20
      if((xx13_i >= -8.0)&&(xx13_i <= 8.0)&&(zz13_i >= -8.0)&&(zz13_i <= 8.0)){
        n13++;
        RefTraj[2] = 1;
        refpnts++;
      }

      // y = -20
      if((xx14_i >= -8.0)&&(xx14_i <= 8.0)&&(zz14_i >= -8.0)&&(zz14_i <= 8.0)){
        n14++;
        RefTraj[3] = 1;
        refpnts++;
      }

      // z = 8
      if((xx15_i >= -8.0)&&(xx15_i <= 8.0)&&(yy15_i >= -20.0)&&(yy15_i <= 20.0)){
        n15++;
        RefTraj[4] = 1;
        refpnts++;
      }

      // z = -8
      else if((xx16_i >= -8.0)&&(xx16_i <= 8.0)&&(yy16_i >= -20.0)&&(yy16_i <= 20.0)){
        n16++;
        RefTraj[5] = 1;
        refpnts++;
      }

      n = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9 + n10;

      if(n == 1){
        err++;
      }

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
          //break;;
        }

        else if(RayTraj[1] == 1){
          Temp1 = intersection2_i;
          x4++;
          //break;;
        }

        else if(RayTraj[2] ==1){
          Temp1 = intersection3_i;
          y15++;
          //break;;
        }

        else if(RayTraj[3] == 1){
          Temp1 = intersection4_i;
          y_15++;
          //break;;
        }

        else if(RayTraj[4] == 1){
          Temp1 = intersection5_i;
          z4++;
          //break;;
        }

        else if(RayTraj[5] == 1){
          Temp1 = intersection6_i;
          z_4++;
          //break;;
        }

        else if(RayTraj[6] == 1){
          Temp1 = intersection7_i;
          _xz5++;
          //break;;
        }

        else if(RayTraj[7] == 1){
          Temp1 = intersection8_i;
          xz5++;
          //break;;
        }

        else if(RayTraj[8] == 1){
          Temp1 = intersection9_i;
          _xz_5++;
          //break;;
        }

        else if(RayTraj[9] == 1){
          Temp1 = intersection10_i;
          xz_5++;
          //break;;
        }

        else if(RefTraj[0] == 1){
          RefTemp1 = intersection11_i;
          RBx8++;
          //break;;
        }

        else if(RefTraj[1] == 1){
          RefTemp1 = intersection12_i;
          RBx_8++;
          //break;;
        }

        else if(RefTraj[2] == 1){
          RefTemp1 = intersection13_i;
          RBy20++;
          //break;;
        }

        else if(RefTraj[3] == 1){
          RefTemp1 = intersection14_i;
          RBy_20++;
          //break;;
        }

        else if(RefTraj[4] == 1){
          RefTemp1 = intersection15_i;
          RBz8++;
          //break;;
        }

        else if(RefTraj[5] == 1){
          RefTemp1 = intersection16_i;
          RBz_8++;
          //break;;
        }

      }

      for(Int_t h = 0; h < 1; h++){
        if((RayTraj[0] == 1)&&(x_4 == 0)){
          Temp2 = intersection1_i;
          //break;;
        }

        else if((RayTraj[1] == 1)&&(x4 == 0)){
          Temp2 = intersection2_i;
          //break;;
        }

        else if((RayTraj[2] == 1)&&(y15 == 0)){
          Temp2 = intersection3_i;
          //break;;
        }

        else if((RayTraj[3] == 1)&&(y_15 == 0)){
          Temp2 = intersection4_i;
          //break;;
        }

        else if((RayTraj[4] == 1)&&(z4 == 0)){
          Temp2 = intersection5_i;
          //break;;
        }

        else if((RayTraj[5] == 1)&&(z_4 == 0)){
          Temp2 = intersection6_i;
          //break;;
        }

        else if((RayTraj[6] == 1)&&(_xz5 == 0)){
          Temp2 = intersection7_i;
          //break;;
        }

        else if((RayTraj[7] == 1)&&(xz5 == 0)){
          Temp2 = intersection8_i;
          //break;;
        }

        else if((RayTraj[8] == 1)&&(_xz_5 == 0)){
          Temp2 = intersection9_i;
          //break;;
        }
        else if((RayTraj[9] == 1)&&(xz_5 == 0)){
          Temp2 = intersection10_i;
          //break;;
        }

        else if((RefTraj[0] == 1)&&(RBx8 == 0)){
          RefTemp2 = intersection11_i;
          //break;;
        }

        else if((RefTraj[1] == 1)&&(RBx_8 == 0)){
          RefTemp2 = intersection12_i;
          //break;;
        }

        else if((RefTraj[2] == 1)&&(RBy20 == 0)){
          RefTemp2 = intersection13_i;
          //break;;
        }

        else if((RefTraj[3] == 1)&&(RBy_20 == 0)){
          RefTemp2 = intersection14_i;
          //break;;
        }

        else if((RefTraj[4] == 1)&&(RBz8 == 0)){
          RefTemp2 = intersection15_i;
          //break;;
        }

        else if((RefTraj[5] == 1)&&(RBz_8 == 0)){
          RefTemp2 = intersection16_i;
          //break;;
        }

      }

    }



    //Find which reference box planes were visible
    // Calculate the distance, pull out x,y,z from temp vectors and use 3D pythagorean

    Double_t temp1x_i = Temp1.X();
    Double_t temp1y_i = Temp1.Y();
    Double_t temp1z_i = Temp1.Z();


    Double_t temp2x_i = Temp2.X();
    Double_t temp2y_i = Temp2.Y();
    Double_t temp2z_i = Temp2.Z();

    Points->SetPoint(count, temp1x_i, temp1y_i, temp1z_i);
    Points2->SetPoint(count, temp2x_i, temp2y_i, temp2z_i);
    count++;

    Double_t distance = sqrt((pow((temp2x_i - temp1x_i),2)) + (pow((temp2y_i - temp1y_i),2)) + (pow((temp2z_i - temp1z_i),2)));

    TestDist = distCalc(x_i, y_i, z_i);

    Fraction = Visible * pow(realboxpnts,-1);
    DivisionHisto->Fill(theta, phi, Fraction);

    if(distance > 0.0){
      dist++;
      lengthHisto-> Fill(distance);
      DistanceHisto->Fill(TestDist, distance);
      //AngleHisto->Fill(thetaR, distance);
      //angAngHisto->Fill(thetaV, thetaR);
      HitsHistogram->Fill(theta, phi, realboxpnts);
    }




    } //END FOR LOOP
  }

  AvgArea = SumArea / INTERVAL;
  cout << "average area: " << AvgArea << endl;

  //Int_t realboxpnts = boxpnts / 2;
  //make sure number of rays and number of distances match
  cout << Visible << endl;
  cout <<  realboxpnts << " Rays Hit the Detector" << endl;
  cout << Fraction << endl;

  cout << "Number of distances: " << dist << endl;
  cout << "Number of errors" << err << endl;
  TCanvas* lengCanvas = new TCanvas("distCanvas","distance Histogram",640,480);
  lengthHisto->Draw();

  TCanvas* lengthCanvas = new TCanvas("testdistCanvas","test distance Histogram",640,480);
  //DistanceHisto->SetMarkerSize(100);
  DistanceHisto-> SetMarkerStyle(kFullDotMedium);
  DistanceHisto->SetMarkerColorAlpha(kBlue, 1.0);
  DistanceHisto->Draw("");





  //TCanvas* angleCanvas = new TCanvas("angleCanvas","angle Vs. Distance",640,480);
  //AngleHisto->Draw("SURF3Z");

  //TCanvas* angleAngleCanvas = new TCanvas("angleAngleCanvas","angle Vs. angle",640,480);
  //angAngHisto->Draw("SURF3Z");

  TCanvas* AreaVSAngleCanvas = new TCanvas("Area(Theta,Phi)","Area(Theta, Phi)",640,480);
  /*
  AreaVSPhiHisto-> SetMarkerStyle(kFullDotMedium);
  AreaVSPhiHisto->SetMarkerColorAlpha(kBlue, 1.0);
  AreaVSPhiHisto->Draw("");
  */
  AreaVSAngleHisto->Draw("lego2");
  //AreaVSAngleHisto -> GetXaxis() -> SetTitle("Theta (degrees)");
  //AreaVSAngleHisto -> GetYaxis() -> SetTitle("Phi (degrees)");
  //AreaVSAngleHisto -> GetZaxis() -> SetTitle("Detector Area (m^2)");


  TH2D *h = (TH2D*)AreaVSAngleHisto->DrawClone("surf3 same");
  h->SetLineColorAlpha(kRed,0.);
  gPad->Modified();
  gPad->Update();
  TCanvas* HitsCanvas = new TCanvas("Hits at each Angle","Hits at each angle",640,480);
  HitsHistogram->Draw("COLz");

  TCanvas* DivisionCanvas = new TCanvas("Fraction","Fraction",640,480);
  DivisionHisto->Draw("COLZ");

  TCanvas* AreaThetaCanvas = new TCanvas("Area vs theta","Area vs theta",640,480);
  AreavsTheta->Draw("");

  TCanvas* AreaPhiCanvas = new TCanvas("Area vs phi","Area vs phi",640,480);
  AreavsPhi->Draw("");



  //TCanvas* AreaVSThetaCanvas = new TCanvas("Area vs Theta","Area Vs. Theta",640,480);
  /*
  AreaVSThetaHisto-> SetMarkerStyle(kFullDotMedium);
  AreaVSThetaHisto-> SetMarkerColorAlpha(kBlue, 1.0);
  AreaVSThetaHisto-> Draw("");
  */
  //AreaVSThetaHisto->Draw("COLZ");
  /*
  TCanvas* PointsCanvas = new TCanvas("Points", "Points", 640, 480);
  Points->SetMarkerStyle(20);
  //Points3->SetMarkerStyle(20);
  Points->SetMarkerColorAlpha(kBlue, 0.5);
  Points2->SetMarkerStyle(20);
  Points2->SetMarkerColorAlpha(kRed, 0.2);
  Points->GetXaxis()->SetTitle("X");
  Points->GetYaxis()->SetTitle("Y");
  Points->GetZaxis()->SetTitle("Z");
  //Points3->Draw("P");
  Points->Draw("p");
  Points2->Draw("SAME P");


  TCanvas* Points2Canvas = new TCanvas("Points2", "Points2", 640, 480);
  Points3->SetMarkerStyle(20);
  Points4->SetMarkerStyle(20);
  Points4->SetMarkerColorAlpha(kBlue, 0.4);
  Points3->Draw("P");
  Points4->Draw("SAME P");
  */
  return 0;
}



//----------------------------------------------------------------------------//
//_______________________________Functions____________________________________//

// Area of octagon for y = 15 and y = -15
// Calculates area of octagon with respect to vertexes and the intersection point
// If the point is within the bounds of the detector, the area will be equal to  52.901
// This function has a bug as the absolute value causes the points to only count when they are positive
bool OctagonAreaCalc(Double_t x, Double_t z){

  Double_t areaCALC = 0;

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

  if(areaCALC <= 52.901){
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

  else return 0;
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

  else return 0;
}

Double_t angleCalc2(Double_t dx, Double_t dy, Double_t dz, Double_t x, Double_t y, Double_t z){
  Double_t MagP = sqrt(x*x + y*y + z*z);
  Double_t MagD = sqrt(dx*dx + dy*dy + dz*dz);
  Double_t thetaV = acos((dx * x + dy * y + dz * z)/(MagP * MagD)) * 100;
  return thetaV;
}


Double_t RefFracCalc(Double_t x, Double_t y, Double_t z){

  Int_t VisiblePlanes[6] = {0, 0, 0, 0, 0, 0}; //(xP, xN, yP, yN, zP, zN)
  Double_t VisibleArea = 0;

  // Vectors from center of each plane to point on sphere
  TVector3 xPVec(x - 4.0, y, z); // x = 8
  TVector3 xNVec(x + 4.0, y, z); // x = -8
  TVector3 yPVec(x, y - 15.0, z); // y = 20
  TVector3 yNVec(x, y + 15.0, z); // y = -20
  TVector3 zPVec(x, y, z - 4.0); // z = 8
  TVector3 zNVec(x, y, z + 4.0); // z = -8

  // Reference box Normal Vectors
  TVector3 xPNorm(1.0, 0.0, 0.0);
  TVector3 xNNorm(-1.0, 0.0, 0.0);
  TVector3 yPNorm(0.0, 1.0, 0.0);
  TVector3 yNNorm(0.0, -1.0, 0.0);
  TVector3 zPNorm(0.0, 0.0, 1.0);
  TVector3 zNNorm(0.0, 0.0, -1.0);

  // Area vectors for each Side
  TVector3 xPArea = 4 * xPNorm;
  TVector3 xNArea = 4 * xNNorm;
  TVector3 yPArea = 15 * yPNorm;
  TVector3 yNArea = 15 * yNNorm;
  TVector3 zPArea = 4 * zPNorm;
  TVector3 zNArea = 4 * zNNorm;

  // Initiate Flux to zero for each plane
  Double_t xPFlux = 0;
  Double_t xNFlux = 0;
  Double_t yPFlux = 0;
  Double_t yNFlux = 0;
  Double_t zPFlux = 0;
  Double_t zNFlux = 0;

  // Now need to determine which planes are visible at a given point on sphere
  // If dot product is positive, angle between vectors is between 0 and 90, the plane is visible
  Double_t Test1 = xPVec.Dot(xPNorm);
  Double_t Test2 = xNVec.Dot(xNNorm);
  Double_t Test3 = yPVec.Dot(yPNorm);
  Double_t Test4 = yNVec.Dot(yNNorm);
  Double_t Test5 = zPVec.Dot(zPNorm);
  Double_t Test6 = zNVec.Dot(zNNorm);

  if(Test1 > 0){
    VisiblePlanes[0] = 1;
    xPFlux = xPArea.Dot(xPVec);
  }

  if(Test2 > 0){
    VisiblePlanes[1] = 1;
    xNFlux = xNArea.Dot(xNVec);
  }

  if(Test3 > 0){
    VisiblePlanes[2] = 1;
    yPFlux = yPArea.Dot(yPVec);
  }

  if(Test4 > 0){
    VisiblePlanes[3] = 1;
    yNFlux = yNArea.Dot(yNVec);
  }

  if(Test5 >0){
    VisiblePlanes[4] = 1;
    zPArea.Dot(zPVec);
  }

  if(Test6 >0){
    VisiblePlanes[5] = 1;
    zNArea.Dot(zNVec);
  }

  //Flux is 0 if not visible, this gives sum of visible sides only
  VisibleArea = xPFlux + xNFlux + yPFlux + yNFlux + zPFlux + zNFlux;

  return VisibleArea;
}
