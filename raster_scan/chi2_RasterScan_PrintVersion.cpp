#include <iostream>
#include <Eigen/LU> 
#include <fstream>
#include <cmath>
#include <math.h>
#include <vector>


using namespace Eigen;
using namespace std;


class ScanMin{

public: 
  double chi2Min;
  double s2t;
  double dm2;

  // Function to find the minimum chi2 per each dm2 value.
  void min(double tempChi2Min, double tempS2t){
    if(tempChi2Min<chi2Min){
      chi2Min=tempChi2Min;
      s2t = tempS2t;
    }
  }  
 

  // Function that returns the values of s2t that must be used to
  // perform the interpolation, in order to get the point to draw the
  // contour.
  bool contourPosition(double tempChi2First, double tempChi2Second, double level){
    double rasterLevel = chi2Min + level;
    if (tempChi2First>tempChi2Second){
      if (tempChi2First>rasterLevel && tempChi2Second<rasterLevel){
	return true;
      }
    }
    else if (tempChi2First<tempChi2Second){
      if (tempChi2First<rasterLevel && tempChi2Second>rasterLevel){
	return true;
      }      
    } 
  }


  ScanMin (double cchi2Min=1e6, double ss2t=5): chi2Min(cchi2Min), s2t(ss2t){
  }

};



double dm2Tic(int index, int N){
  double dm2min = 1.e-02;
  double dm2max = 1.e02;
    return pow(10, (log10(dm2min) + index*(log10(dm2max)-log10(dm2min))/(N-1)) );
  }

double s2tTic(int index, int N){
    double s2tmin = 3.e-04;
    double s2tmax = 1.;
    return pow(10, (log10(s2tmin) + index*(log10(s2tmax)-log10(s2tmin))/(N-1)));
  }

double linearInterpolate(double x1, double y1, double x2, double y2, double y){
  return x1 + ((y-y1)*(x2-x1)/(y2-y1));
}

int main(){

  static const int N_s2t = 200;
  static const int N_dm2 = 100;
  static const int rows = N_dm2*N_s2t;
  static const int columns = 3;
  static const int NB = 8;
  static const int numberMC = 9934;
  double binlim[NB+1] = {475., 550., 675., 800., 950., 1100., 1300., 1500., 3000.};
  double daEnu[NB] = {83., 90., 63., 58., 49., 45., 35., 67.};
  double bgEnu[NB] = {70.91, 81.51, 61.94, 58.95, 50.01, 44.49, 33.73, 65.27};


  double fracBgMat[NB][NB];
  ifstream fracBg ("../475data/miniboone_nuebgr_fractcovmatrix.txt");
  if(fracBg.is_open()){
    for(int i=0; i<NB; i++){
      for(int j=0; j<NB; j++){
	fracBg >> fracBgMat[i][j];	
      }
    }
    fracBg.close();
  }
  else cout << "Unable to open file"<< endl;


  double fracSgMat[NB][NB];
  ifstream fracSg ("../475data/miniboone_numunuefullosc_fractcovmatrix.txt");
  if(fracSg.is_open()){
    for(int i=0; i<NB; i++){
      for(int j=0; j<NB; j++){
	fracSg >> fracSgMat[i][j];
      }
    }
    fracSg.close();
  }
  else cout << "Unable to open file "<< endl;


  
  // ntuple file MC data//
  double enuqe[numberMC], enu[numberMC], l[numberMC], wgt[numberMC];
  ifstream myfile ("../475data/miniboone_numunuefullosc_ntuple.txt");
  if (myfile.is_open()){
    for(int i=0; i<numberMC; i++)
      myfile >> enuqe[i] >> enu[i] >> l[i] >> wgt[i];
    myfile.close();
  }
  else cout << "Unable to open file" << endl; 

  
  
  // Background covariance matrix, doesn't depend on the parameters. 
  // Build from the background fractional covariance matrix. 
  double bgMat[NB][NB];
  for(int i=0; i<NB; i++){
    for(int j=0; j<NB; j++){
      bgMat[i][j] = bgEnu[i]*bgEnu[j]* fracBgMat[i][j];
    }
  }

  MatrixXd invBgMat(NB,NB);
  //  double detBgMat;
  for(int i=0; i<NB; i++){
    for(int j=0; j<NB; j++){
      invBgMat(i,j) = bgMat[i][j];
    }
  }
  //  detBgMat = fabs(invBgMat.determinant());
  invBgMat = invBgMat.inverse();


  double sgMat[NB][NB];
  MatrixXd invMat(NB,NB);   
  //  double det;
  MatrixXd totMM(NB,NB);
  
  // Initiation of the parameters sin^2(\theta) and \delta m^2
  static const double level90CL = 1.6448;//I have the impression this should be 2.7055
  static const double level2S = 4.0;//Or maybe this should be 2
  static const double level3S = 9.0;//This 3
  static const double level5S = 25.0;//This 5

  double chi2[N_s2t][N_dm2] = {0.0};
  //  double chi2C[N_s2t][N_dm2]= {0.0};
  //  double chi2Sens[N_s2t][N_dm2] = {0.0};
  double sgEnu[NB];
  double sgEnuIterative[NB];
  // double s2t;
  // double dm2;
  vector <ScanMin> rasterScanMin;
  //  vector <ScanMin> rasterScanMinDet;

  
  double chi2_min  = 1e5;
  double min_position[2];
  //  double chi2C_min = 1e5;
  //  double min_positionC[2];
  double chi2MB_min = 1e5;
  double min_positionMB[2];
  
  int flag;
  int maxIteration = 4;
  double tolerance = 1.e-2;
  double chi2MinOld;
  //  double chi2MinDetOld;
  double s2tMinOld;
  //  double s2tMinDetOld;
  
  // double ntupleMB[rows][columns];
  // ifstream myfile2 ("../475data/chi2surface.txt");
  // if (myfile2.is_open() ){
  //   while ( myfile2.good() ){
  //     for(int i=0; i<rows; i++){
  // 	for(int j=0; j<columns; j++){
  // 	  myfile2 >> ntupleMB[i][j];
  // 	}
  //     }
  //   }
  //   myfile2.close();
  // }
  // else cout << "Unable to open file"<< endl; 


  // double chi2MB[N_s2t][N_dm2];
  // for(int i=0; i<N_dm2; i++){
  //   for(int j=0; j<N_s2t; j++){
  //     chi2MB[j][i]=ntupleMB[i*N_s2t+j][2];
  //   }
  // }

  // for(int k = 0; k<rows; k++){
  //     if (ntupleMB[k][2]<chi2MB_min) {
  // 	chi2MB_min = ntupleMB[k][2];
  // 	min_positionMB[0]=  ntupleMB[k][1];
  // 	min_positionMB[1]=  ntupleMB[k][0]; 
  //     }

  // }



  //  Loop trough a variety of values of s2t and dm2
  for(int kk = 0; kk<N_dm2; kk++){
    rasterScanMin.push_back(ScanMin());
    rasterScanMin[kk].dm2 = dm2Tic(kk, N_dm2);
    //    rasterScanMinDet.push_back(ScanMin());
    //    rasterScanMinDet[kk].dm2=dm2;

    flag = 0;    
    chi2MinOld = 1e5;
    //    chi2MinDetOld = 1e5;
    s2tMinOld = 2; 
    //    s2tMinDetOld = 2;

    while(flag<maxIteration && fabs(chi2MinOld-rasterScanMin[kk].chi2Min)>tolerance ){
      // cout << "flag   " << flag << endl;
      // cout << scientific << "chi2MinOld   " << chi2MinOld << endl;
      // cout << scientific << "s2tMinOld   " << s2tMinOld << endl;
      // cout << scientific << "rasterScanMin[kk].chi2Min   " << rasterScanMin[kk].chi2Min << endl;
      // cout << scientific << "rasterScanMin[kk].s2t   " << rasterScanMin[kk].s2t << endl;

      chi2MinOld = rasterScanMin[kk].chi2Min;
      s2tMinOld = rasterScanMin[kk].s2t;
      //      chi2MinDetOld = rasterScanMinDet[kk].chi2Min;
      //      s2tMinDetOld = rasterScanMinDet[kk].s2t;

      for(int k=0; k<N_s2t; k++){
	chi2[k][kk] = 0.;
      }


      for(int k=0; k<N_s2t; k++){
      
	// Signal of the MC reconstruction, it depends on the parameters
	for(int j=0; j<NB; j++){
	  sgEnu[j] = 0.0;
	}
	for(int i=0; i<numberMC; i++){
	  for(int ii = 0; ii<NB; ii++){
	    if(binlim[ii]< enuqe[i] && enuqe[i] <= binlim[ii+1]) // cambie menor igual por menor
	      sgEnu[ii]+= ((wgt[i]/numberMC)*s2tTic(k, N_s2t)*(sin(1.27*dm2Tic(kk, N_dm2)*(l[i]*0.01)/enu[i]))*(sin(1.27*dm2Tic(kk, N_dm2)*(l[i]*0.01)/enu[i])));
	  }
	}


	if (flag==0){
      	  invMat = invBgMat;  
	  //	  det = detBgMat;

	  // for(int i=0; i<NB; i++){
	  //   for(int j=0; j<NB; j++){
	  //     chi2Sens[k][kk] += (sgEnu[i]*invMat(i,j)*sgEnu[j]);
	  //   }
	  // }

	}
	
	else{

	// Signal of the MC reconstruction for the iterative process.
	for(int j=0; j<NB; j++){
	  sgEnuIterative[j] = 0.0;
	}
	for(int i=0; i<numberMC; i++){
	  for(int ii = 0; ii<NB; ii++){
	    if(binlim[ii]< enuqe[i] && enuqe[i] <= binlim[ii+1]) // cambie menor igual por menor
	      sgEnuIterative[ii]+= ((wgt[i]/numberMC)*rasterScanMin[kk].s2t*(sin(1.27*dm2Tic(kk, N_dm2)*(l[i]*0.01)/enu[i]))*(sin(1.27*dm2Tic(kk, N_dm2)*(l[i]*0.01)/enu[i])));
	  }
	}
      

	  // Build from the signal fractional covariance matrix. 
	  for(int i=0; i<NB; i++){
	    for(int j=0; j<NB; j++){
	      sgMat[i][j] = sgEnuIterative[i]*sgEnuIterative[j]* fracSgMat[i][j];
	    }
	  }
      
	  // Total Matrix, the sum of the Signal and Background matrix
	  // plus a term in the diagonal which accounts for the
	  // statistical uncertainty associated with the nue signal
	  // prediction
	  for (int i=0;i<NB;i++){
	    for(int j=0;j<NB;j++){
	      totMM(i,j) = bgMat[i][j] + sgMat[i][j];
	      if(i==j){
		totMM(i,j) += sgEnuIterative[i]; 
	      }
	    }
	  }




	  //	  det = fabs(totMM.determinant());

	  // if(det<0.00001) {
	  //   cout << "Warning!   Determinant= " << det << endl;
	  //   cout << "Index:   "<< k<<"  "<< kk<< endl;
	  //   cout << "It is not invertible." << endl;
	  // }
      


	  invMat = totMM.inverse();

	}//else "flag!=0"

	// Chi^2
	for(int i=0; i<NB; i++){
	  for(int j=0; j<NB; j++){
	    chi2[k][kk]+= (daEnu[i] - (bgEnu[i] + sgEnu[i])) *invMat(i,j)* (daEnu[j] - ( bgEnu[j] + sgEnu[j]));
	    //////	  chi2[k][kk]= chi2[k][kk] + ((daEnu[i] - bgEnu[i] - sgEnu[i]) *invMat(i,j)* (daEnu[j] - bgEnu[j] - sgEnu[j]));
	  }
	}
	//	chi2C[k][kk] = chi2[k][kk] + log(det);
      
	//Finding the minimum chi2 for each dm2. As must be done for the
	//raster scan method.
	rasterScanMin[kk].min(chi2[k][kk], s2tTic(k, N_s2t));
	//	rasterScanMinDet[kk].min(chi2C[k][kk], s2t);
      
	if (isnan(chi2[k][kk])) 
	  cout<< "chi2 :  "<<chi2[k][kk]<<"   CHANGOS NAN!   Index:  "<< k<<"  "<< kk<< endl;

	// if (isnan(chi2C[k][kk])) 
	//   cout<< "chi2C:  "<<chi2C[k][kk]<<"   CHANGOS NAN!   Index:  "<< k<<"  "<< kk<< endl;
      
	if (chi2[k][kk] <chi2_min ){
	  chi2_min  = chi2[k][kk];
	  min_position[0]=  s2tTic(k, N_s2t);
	  min_position[1]=  dm2Tic(kk, N_dm2); 
	}
      
	// if (chi2C[k][kk]<chi2C_min){
	//   chi2C_min = chi2C[k][kk];
	//   min_positionC[0]=  s2t;
	//   min_positionC[1]=  dm2; 
	// }
      
      }//for k<N_s2t
      flag++;
    }//while: flag && tolerance

    ////    cout << kk << endl;//Counter
    // cout << " Exit flag:  " << flag << endl;

    // cout << scientific << "chi2MinOld   " << chi2MinOld << endl;
    // cout << scientific << "rasterScanMin[kk].chi2Min   " << rasterScanMin[kk].chi2Min << endl;
    // cout << scientific << "Minima difference   " << chi2MinOld-rasterScanMin[kk].chi2Min << endl;
    // cout << scientific << "s2t difference   " << s2tMinOld-rasterScanMin[kk].s2t << endl;

    if(!(kk%10)) cout<< kk << endl;//Counter
  }// for kk

  //      Here finishes the loop s2t dm2


  cout <<"#######################" <<endl;
  cout <<"The Minima are:" <<endl;
  cout <<"chi2_min = "<< chi2_min <<endl;
  cout<< "Located in: "<< endl;
  cout<< "s2t = "<< min_position[0] << "     "<< "dm2 = " << min_position[1] <<endl;
  cout <<"#######################" <<endl;
  // cout <<"chi2C_min = "<< chi2C_min <<endl;
  // cout<< "Located in: "<< endl;
  // cout<< "s2t = "<< min_positionC[0] << "     "<< "dm2 = " << min_positionC[1] <<endl;
  // cout <<"#######################" <<endl;
  // cout <<"chi2MB_min = "<< chi2MB_min <<endl;
  // cout<< "Located in: "<< endl;
  // cout<< "s2t = "<< min_positionMB[0] << "     "<< "dm2 = " << min_positionMB[1] <<endl;
  // cout <<"#######################" <<endl;
  


  // Writing ntuple chi2 and chi2C on a text file.
  // Using the following format:
  // (Dm2 (eV2),  sin2(2th), -2ln(L))
  ofstream myfilentuple ("chi2ntuplePrint.txt");
  if (myfilentuple.is_open())
    {
      myfilentuple.precision(8);
      for(int kk = 0; kk<N_dm2; kk++){
  	for(int k = 0; k<N_s2t; k++){
  	    myfilentuple << scientific <<  dm2Tic(kk, N_dm2) << "\t" << s2tTic(k, N_s2t) << "\t" << chi2[k][kk] <<"\n";
  	}
      }
      myfilentuple.close();
    }
  else cout << "Unable to open file"<<endl;
  cout << "Printing chi2ntuplePrint.txt"<<endl; 
  

  // ofstream myfilentupleC ("chi2ntupleC.txt");
  // if (myfilentupleC.is_open())
  //   {
  //     myfilentupleC.precision(8);
  //     for(int kk = 0; kk<N_dm2; kk++){
  // 	dm2 = dm2Tic(kk, N_dm2);
  // 	for(int k = 0; k<N_s2t; k++){
  // 	  s2t = s2tTic(k, N_s2t);
  //     	  myfilentupleC << fixed <<  dm2 << "\t" << s2t << "\t" << chi2C[k][kk] <<"\n";
  // 	}
  //     }
  //     myfilentupleC.close();
  //   }
  // else cout << "Unable to open file"<<endl;
  // cout << "Printing chi2ntupleC.txt"<<endl; 



  // This part of the code finds the intersection of each chi2(s2t)
  // for each dm2 with the different levels (90CL, 3S, 5S). And prints
  // the intersection (using a linear interpolation).  As it is now it
  // is incapable of printing several intersections (actually it stops
  // looking after it finds the first).
  ofstream myfileRasterContoursIterative ("rasterContoursIterative_PrintVersion.txt");
  if(myfileRasterContoursIterative.is_open()){
    for(int kk=0; kk<N_dm2; kk++){
      myfileRasterContoursIterative << fixed <<  dm2Tic(kk, N_dm2) << "\t" ;    

      for(int k=0; k<N_s2t-1; k++){
	if(rasterScanMin[kk].contourPosition(chi2[k][kk], chi2[k+1][kk], level90CL) ){
	  myfileRasterContoursIterative << linearInterpolate(s2tTic(k, N_s2t), chi2[k][kk], s2tTic(k+1, N_s2t), chi2[k+1][kk], rasterScanMin[kk].chi2Min + level90CL);
	  break;
	}
	else if(k==N_s2t-2)
	  myfileRasterContoursIterative << 1.5; // A point with s2t of 1.5 will not be plotted (lame fix)
      }
      myfileRasterContoursIterative << "\t"; 

      for(int k=0; k<N_s2t-1; k++){
	if(rasterScanMin[kk].contourPosition(chi2[k][kk], chi2[k+1][kk], level2S) ){
	  myfileRasterContoursIterative << linearInterpolate(s2tTic(k, N_s2t), chi2[k][kk], s2tTic(k+1, N_s2t), chi2[k+1][kk], rasterScanMin[kk].chi2Min + level2S);
	  break;
	}
	else if(k==N_s2t-2)
	  myfileRasterContoursIterative << 1.5;
      }
      

      myfileRasterContoursIterative << "\t"; 

      for(int k=0; k<N_s2t-1; k++){
	if(rasterScanMin[kk].contourPosition(chi2[k][kk], chi2[k+1][kk], level3S) ){
	  myfileRasterContoursIterative << linearInterpolate(s2tTic(k, N_s2t), chi2[k][kk], s2tTic(k+1, N_s2t), chi2[k+1][kk], rasterScanMin[kk].chi2Min + level3S);
	  break;
	}
	else if(k==N_s2t-2)
	  myfileRasterContoursIterative << 1.5;
      }
      

      myfileRasterContoursIterative << "\t"; 

      for(int k=0; k<N_s2t-1; k++){
      	if(rasterScanMin[kk].contourPosition(chi2[k][kk], chi2[k+1][kk], level5S) ){
      	  myfileRasterContoursIterative << linearInterpolate(s2tTic(k, N_s2t), chi2[k][kk], s2tTic(k+1, N_s2t), chi2[k+1][kk], rasterScanMin[kk].chi2Min + level5S);
      	  break;
      	    }
      	else if(k==N_s2t-2)
      	  myfileRasterContoursIterative << 1.5;
      }
      myfileRasterContoursIterative << "\n";
    }
    myfileRasterContoursIterative.close();
  }
  else cout << "Unable to open rasterContoursIterative_PrintVersion.txt file"<< endl;
  cout << "Printing rasterContoursIterative_PrintVersion.txt" << endl;



  // ofstream myfileRasterContoursSensibility ("rasterContours_sens.txt");
  // if(myfileRasterContoursSensibility.is_open()){
  //   for(int kk=0; kk<N_dm2; kk++){
  //     myfileRasterContoursSensibility << fixed <<  dm2Tic(kk, N_dm2) << "\t" ;    

  //     for(int k=0; k<N_s2t-1; k++){
  // 	if(rasterScanMin[kk].contourPosition(chi2Sens[k][kk], chi2Sens[k+1][kk], level90CL) ){
  // 	  myfileRasterContoursSensibility << linearInterpolate(s2tTic(k, N_s2t), chi2Sens[k][kk], s2tTic(k+1, N_s2t), chi2Sens[k+1][kk], rasterScanMin[kk].chi2Min + level90CL);
  // 	  break;
  // 	}
  // 	else if(k==N_s2t-2)
  // 	  myfileRasterContoursSensibility << 1.5; // A point with s2t of 1.5 will not be plotted (lame fix)
  //     }
  //     myfileRasterContoursSensibility << "\t"; 

  //     for(int k=0; k<N_s2t-1; k++){
  // 	if(rasterScanMin[kk].contourPosition(chi2Sens[k][kk], chi2Sens[k+1][kk], level3S) ){
  // 	  myfileRasterContoursSensibility << linearInterpolate(s2tTic(k, N_s2t), chi2Sens[k][kk], s2tTic(k+1, N_s2t), chi2Sens[k+1][kk], rasterScanMin[kk].chi2Min + level3S);
  // 	  break;
  // 	}
  // 	else if(k==N_s2t-2)
  // 	  myfileRasterContoursSensibility << 1.5;
  //     }
  //     myfileRasterContoursSensibility << "\t"; 

  //     for(int k=0; k<N_s2t-1; k++){
  // 	if(rasterScanMin[kk].contourPosition(chi2Sens[k][kk], chi2Sens[k+1][kk], level5S) ){
  // 	  myfileRasterContoursSensibility << linearInterpolate(s2tTic(k, N_s2t), chi2Sens[k][kk], s2tTic(k+1, N_s2t), chi2Sens[k+1][kk], rasterScanMin[kk].chi2Min + level5S);
  // 	  break;
  // 	    }
  // 	else if(k==N_s2t-2)
  // 	  myfileRasterContoursSensibility << 1.5;
  //     }
  //     myfileRasterContoursSensibility << "\n";
  //   }
  //   myfileRasterContoursSensibility.close();
  // }
  // else cout << "Unable to open rasterContours_sens.txt file"<< endl;
  // cout << "Printing rasterContours_sens.txt" << endl;




  // ofstream myfileRasterContoursDet ("rasterContoursDet.txt");
  // if(myfileRasterContoursDet.is_open()){
  //   for(int kk=0; kk<N_dm2; kk++){
  //     myfileRasterContoursDet << fixed <<  dm2Tic(kk, N_dm2) << "\t" ;    

  //     for(int k=0; k<N_s2t-1; k++){
  // 	if(rasterScanMinDet[kk].contourPosition(chi2C[k][kk], chi2C[k+1][kk], level90CL) ){
  // 	  myfileRasterContoursDet << linearInterpolate(s2tTic(k, N_s2t), chi2C[k][kk], s2tTic(k+1, N_s2t), chi2C[k+1][kk], rasterScanMinDet[kk].chi2Min + level90CL);
  // 	  break;
  // 	}
  // 	else if(k==N_s2t-2)
  // 	  myfileRasterContoursDet << 1.5;
  //     }
  //     myfileRasterContoursDet << "\t"; 

  //     for(int k=0; k<N_s2t-1; k++){
  // 	if(rasterScanMinDet[kk].contourPosition(chi2C[k][kk], chi2C[k+1][kk], level3S) ){
  // 	  myfileRasterContoursDet << linearInterpolate(s2tTic(k, N_s2t), chi2C[k][kk], s2tTic(k+1, N_s2t), chi2C[k+1][kk], rasterScanMinDet[kk].chi2Min + level3S);
  // 	  break;
  // 	}
  // 	else if(k==N_s2t-2)
  // 	  myfileRasterContoursDet << 1.5;
  //     }
  //     myfileRasterContoursDet << "\t"; 

  //     for(int k=0; k<N_s2t-1; k++){
  // 	if(rasterScanMinDet[kk].contourPosition(chi2C[k][kk], chi2C[k+1][kk], level5S) ){
  // 	  myfileRasterContoursDet << linearInterpolate(s2tTic(k, N_s2t), chi2C[k][kk], s2tTic(k+1, N_s2t), chi2C[k+1][kk], rasterScanMinDet[kk].chi2Min + level5S);
  // 	  break;
  // 	    }
  // 	else if(k==N_s2t-2)
  // 	  myfileRasterContoursDet << 1.5;
  //     }
  //     myfileRasterContoursDet << "\n";
  //   }
  //   myfileRasterContoursDet.close();
  // }
  // else cout << "Unable to open rasterContoursDet.txt file"<< endl;
  // cout << "Printing rasterContoursDet.txt" << endl;





  return 0; 

}


