#include<iostream>
#include<vector>
#include <fstream>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<limits>
#include<cfloat>
#include<algorithm>
//#define DEBUG
using namespace std;

vector<vector<double> > pridownscale_smos(vector<vector<double> > datadownscale, vector<vector<double> > dataupscale, vector<double> sigma, double beta, double no_time){
	int n = datadownscale.size();
	int m = dataupscale.size();
	
	vector<vector<double> > data;
	vector<vector<double> > datanew;
	for(int i=0; i<datadownscale.size(); i++){
		data.push_back(datadownscale[0]);		
	}
	
	for(int i=0; i<dataupscale.size(); i++){
		datanew.push_back(dataupscale[0]);		
	}
	
	for(int t=0; t<no_time; t++){
		vector<double> temp;
		double sum =0.0;
		if(beta == 1){
			for(int i=0; i< datanew.size(); i++){
				for(int j=0; j< data.size(); j++){
					double var1 = pow((datanew[i][1] - data[j][1]),2) / (2 * sigma[0] * sigma[0]);
					double var2 = pow((datanew[i][2] - data[j][2]),2) / (2 * sigma[1] * sigma[1]);
					double var3 = pow((datanew[i][3] - data[j][3]),2) / (2 * sigma[2] * sigma[2]);
					
					double result = exp(-(var1+ var2+ var3));
					temp.push_back(result);
				}
				
				vector<double> temp1;
								
				for(int l=0; l< data[0].size(); l++){
					double s=0.0;
					for(int k=0; k< temp.size(); k++){
						s= s + data[k][l] * temp[k];
					}
					temp1.push_back(s/sum);
				}
				
				if(! isnan(temp1[0])){
					datanew[i]=temp1;
				}
			}
		}
	}
	return datanew;
}

void trim1D(vector<int> iSet, vector<double> set , vector<double>  &trimResult){
	for(int i=0; i< iSet.size(); i++){
  		trimResult.push_back(set[iSet[i]]);
  	}
}

// Function to read multidimentional files with double values
void readDoubleFilesNew(char fileName[], vector<vector<double> > &resultVector ){
	string line;
	ifstream myfile (fileName);
  if (myfile.is_open())
  {
	while (getline (myfile,line) )
    {
    	vector<double> temp;
		char *c;
		c=strtok((char*)line.c_str(),",");
		while(c != NULL){
			temp.push_back(atof(c));
			c=strtok(NULL,",");
		}
		resultVector.push_back(temp);
    }
    myfile.close();	
  }
  else 
  	cout << "Unable to open file"; 
}

// Function to read multidimentional files with NAN values
void readDoubleFilesWithNANValues(char fileName[], vector<vector<double> > &resultVector ){
	string line;
	ifstream myfile (fileName);
  if (myfile.is_open())
  {
	while (getline (myfile,line) )
    {
    	vector<double> temp;
		char *c;
		c=strtok((char*)line.c_str(),",");
		while(c != NULL){
			//cout<<c<<" ";
			temp.push_back(atof(c));
			c=strtok(NULL,",");
		}
		resultVector.push_back(temp);
    }
    myfile.close();	
  }
  else 
  	cout << "Unable to open file"; 
}

// Function to read single dimentional files with double values
void read1DDoubleFiles(char fileName[], vector<double> &data){
	string line;
	ifstream myfile (fileName);
  	if (myfile.is_open())
  	{
		while (getline (myfile,line) )
    	{
    		char *c;	
			c=strtok((char*)line.c_str(),",");
			data.push_back(atof(c));
	}
    myfile.close();	
  }
  else 
  	cout << "Unable to open file "<< myfile; 
}

void minDis(vector<vector<double> > longitude, vector<vector<double> > latitude, double stat_lon, double stat_lat, int &row, int &column){
	int n=longitude.size();
	int m= longitude[0].size();
	double d=DBL_MAX ;
	for(int i=0; i< n; i++ ){
		for(int j=0; j< m ; j++){
			double lon = fabs(longitude[i][j]);
			double lat = fabs (latitude[i][j]);
			double currentD= sqrt(pow(lon - stat_lon,2.0) + pow(lat -stat_lat,2.0)); 
			if(currentD < d){
				row=i;
				column=j;
				d=currentD;
				//cout<<row<<" "<<column;
			}
		}
	}
}

int main(){
	vector<vector<double> > cm;
	readDoubleFilesNew("cm.csv", cm);
	//Line Number 10
	cm.erase(cm.begin());
	//Line Number 11
	reverse(cm.begin(), cm.end());
	/*cout<<"cm size "<<cm.size();
	for(int i=0; i< cm.size(); i++){
		for(int j=0; j< cm[i].size(); j++){
			cout<<cm[i][j]<<" ";
		}
		cout<<endl;
	}*/
	
	//Reading values for Line 15
	vector<vector<double> > latgridB_1km;
	readDoubleFilesNew("latgridB_1km.csv", latgridB_1km );
	/*cout<<"latgridB_1km size "<<latgridB_1km.size();
	for(int i=0; i< latgridB_1km.size(); i++){
		for(int j=0; j< latgridB_1km[i].size(); j++){
			cout<<latgridB_1km[i][j]<<" ";
		}
		cout<<endl;
	}*/
	
	vector<vector<double> > longridB_1km;
	readDoubleFilesNew("longridB_1km.csv", longridB_1km );
	/*cout<<"longridB_1km size "<<longridB_1km.size();
	for(int i=0; i< longridB_1km.size(); i++){
		for(int j=0; j< longridB_1km[i].size(); j++){
			cout<<longridB_1km[i][j]<<" ";
		}
		cout<<endl;
	}*/
	
	
	vector<vector<double> > lc12_2d_1km;
	readDoubleFilesWithNANValues("lc12_2d_1km.csv",lc12_2d_1km);
	cout<<"lc12_2d_1km size "<<lc12_2d_1km.size()<<endl;
	cout<<"lc12_2d_1km columns "<<lc12_2d_1km[0].size()<<endl;
	
	int row, col;
	double stat_lat = 29.17;
	double stat_lon= 53.7;
	minDis(longridB_1km, latgridB_1km, stat_lon, stat_lat, row, col);
	cout<<row<<"  "<<col<<endl;
	//Line Number 16
	int index_val_st= row * col;
	//Line Number  17 . Variable is never used and value of index is not known
	//tx_val_SM = sm_est(index,:)';
	
	//Line Number 40
	vector<int> smos_ind;
	for(int i=1; i<=113; i++){
		smos_ind.push_back(i);
	}
	for(int i=24; i<=240; i++){
		smos_ind.push_back(i);
	}
	
	//Reading files for line 41 onwards
	vector<vector<double> > smos_sm_1km;
	readDoubleFilesNew("smos_sm_1km.csv", smos_sm_1km );
	cout<<"smos_sm_1km size "<<smos_sm_1km.size()<<endl;

	vector<vector<double> > sm_est;
	readDoubleFilesNew("sm_est.csv", sm_est );
	
	//Line Number 41
	vector<vector<double> > sm_smos;
	for(int i=0; i< smos_sm_1km.size(); i++){
		vector<double> temp;
		for(int j=0; j<sm_est.size(); j++){
			temp.push_back(0);
		}
		sm_smos.push_back(temp);
	}
	cout<<sm_smos.size()<<endl;
	cout<<sm_smos[0].size()<<endl;
	
	
	/*
	vector<vector<double> > sm_smos(smos_sm_1km.size(), vector<double>(sm_est.size()));
	cout<<sm_smos.size()<<endl;
	*/
	
	//Line Number 42
	vector<vector<double> > smos_sm_25km;
	readDoubleFilesNew("smos_sm_25km.csv", smos_sm_25km);

	vector<vector<double> > sm_smos_25;
	for(int i=0; i< smos_sm_25km.size(); i++){
		vector<double> temp;
		for(int j=0; j<sm_est.size(); j++){
			temp.push_back(0);
		}
		sm_smos_25.push_back(temp);
	}
	
	//Line Number 43 --> value is never used
	/*
	vector<vector<double> > sm_pri(smos_sm_1km.size(), vector<double>(sm_est.size()));
	cout<<sm_pri.size()<<endl;
	*/
	
	//Line Number 44 and 45
	for(int i=0; i<smos_ind.size(); i++){
		for(int j=0; j< sm_smos.size(); j++){
			sm_smos[j][i]=smos_sm_1km[j][i];
		}
		
		for(int j=0; j< sm_smos_25.size(); j++){
			sm_smos_25[j][i]=smos_sm_25km[j][i];
		}
	}
	//Line Number 47 - 52
	double mini=0.03;
	double maxi=0.6;
	double a1[]={0.1300, 0.5838, 0.2108, 0.3412};
	double a2[]={0.5703, 0.5838, 0.2108, 0.3429};
	double c1[]={0.3795, 0.5798, 0.0319, 0.3412};
	double c2[]={0.8187, 0.5798, 0.0319, 0.3429};
	
	//Line Number 69-71
	vector<double> sigma;
	sigma.push_back(0.01);
	sigma.push_back(0.01);
	sigma.push_back(0.05);
	int beta =2 ;
	int no_time=2;
	
	//Line Number 73
	vector<double> data_lat;
	for(int i=0; i<latgridB_1km[0].size(); i++){
		for(int j=0; j<latgridB_1km.size(); j++){
			data_lat.push_back(latgridB_1km[j][i]);
		}
	}
	//Line Number 74
	vector<double> data_lon;
	for(int i=0; i<longridB_1km[0].size(); i++){
		for(int j=0; j<longridB_1km.size(); j++){
			data_lon.push_back(longridB_1km[j][i]);
		}
	}
	
	//Line Number 77
	cout<<"Calculating lc"<<endl;
	vector<double> lc;
	for(int i=0; i<361; i++){
		for(int j=0; j<309; j++){
			lc.push_back(lc12_2d_1km[j][i]);
		}
	}
	cout<<"lc size "<<lc.size();
	
	//Line Number 78
	vector<int> ind_nonagri;
	for(int i=0; i<lc.size(); i++){
		if(lc[i]!=12)
			ind_nonagri.push_back(i);
	}
	cout<<"ind_nonagri size "<<ind_nonagri.size()<<endl;
	
	//Line Number 79
	vector<int> ind_all;
	for(int i=0; i<lc.size(); i++){
		ind_all.push_back(i);
	}
	//Line Number 80
	for(int i=0; i<sm_est[0].size(); i++){
		//Line number 82
		vector<double> sm_esti;
		vector<vector<double> > datanew;
		for(int l=0; l<data_lat.size(); l++){
			vector<double> temp;
			temp.push_back(data_lat[l]);
			temp.push_back(data_lon[l]);
			temp.push_back(sm_est[l][i]);
			sm_esti.push_back(sm_est[l][i]);
			datanew.push_back(temp);
		}
		
		//Line Number 85
		vector<int> ind_smos;
		vector<double> sm_smosi;
		for(int j=0; j< sm_smos.size(); j++){
			sm_smosi.push_back(sm_smos[j][i]);
			if(sm_smos[j][i] < 0.2)
				ind_smos.push_back(i);
		}
		
		//Finding union between ind_nonagri and ind_smos
		//Line Number 86
		int iind_nonagri=ind_nonagri.size();
		int iind_smos= ind_smos.size();
		int itr1=0;
		int itr2=0;
		vector<int> ind_remove;
		while(itr1<iind_nonagri && itr2 < iind_smos){
			if(ind_nonagri[itr1] < ind_smos[itr2]){
				ind_remove.push_back(itr1);
				itr1++;
			}
			else if(ind_nonagri[itr1] > ind_smos[itr2]){
				ind_remove.push_back(itr2);
				itr2++;
			}
			else{
				ind_remove.push_back(itr1);
				itr1++;
				itr2++;
			}
		}
		
		while(itr1 < iind_nonagri){
			ind_remove.push_back(itr1);
			itr1++;
		}
		
		while(itr2 < iind_smos){
			ind_remove.push_back(itr2);
			itr2++;
		}
		
		//Finding Setdiff between ind_all and ind_remove
		//Line Number 87
		int iind_all=ind_all.size();
		int iind_remove= ind_remove.size();
		itr1=0;
		itr2=0;
		vector<int> ind_net;
		while(itr1 < iind_all && itr2 < iind_remove){
			if(ind_all[itr1] < ind_remove[itr2]){
				ind_net.push_back(itr1);
				itr1++;
			}
			else {
				itr2++;
			}
		}
		
		while(itr1 < iind_all){
			ind_net.push_back(itr1);
			itr1++;
		}
		
		//Finding Setdiff between ind_all and ind_nonagri
		//Line Number 88
		itr1=0;
		itr2=0;
		vector<int> ind_net2;
		while(itr1<iind_all && itr2 < iind_nonagri){
			if(ind_all[itr1] < ind_nonagri[itr2]){
				ind_net2.push_back(itr1);
				itr1++;
			}
			else {
				itr2++;
			}
		}
		
		while(itr1 < iind_all){
			ind_net2.push_back(itr1);
			itr1++;
		}
		
		vector<double> parameter1;
		trim1D(ind_net2, data_lat, parameter1);
		
		vector<double> parameter2;
		trim1D(ind_net2, data_lon, parameter2);
		
		vector<double> parameter3;
		trim1D(ind_net2, sm_esti, parameter3);	
		
		//Line Number 90
		vector<vector<double> > datadownscale;
		for(int j=0; j<ind_net2.size(); j++){
			vector<double> temp;
			temp.push_back(parameter1[j]);
			temp.push_back(parameter2[j]);
			temp.push_back(parameter3[j]);
			datadownscale.push_back(temp);
		}
		
		vector<double> parameter4;
		trim1D(ind_net, data_lat, parameter4);
		
		vector<double> parameter5;
		trim1D(ind_net, data_lon, parameter5);
		
		vector<double> parameter6;
		trim1D(ind_net, sm_smosi, parameter6);	
		
		//Line Number 91
		vector<vector<double> > dataupscale;
		for(int j=0; j<ind_net.size(); j++){
			vector<double> temp;
			temp.push_back(parameter4[j]);
			temp.push_back(parameter5[j]);
			temp.push_back(parameter6[j]);
			dataupscale.push_back(temp);
		}
		
		//Line Number 92
		vector<vector<double> > datanew_pre= pridownscale_smos(datadownscale, dataupscale, sigma, beta, no_time);
		for(int l=0; l<ind_net.size(); l++ ){
			datanew[l][2]= datanew_pre[l][2];
		}
		
		//Line Number 96
		for(int j=0; j<lc.size(); j++){
			if( lc[j] != 12){
				datanew[j][2]=NAN;
			}
		}
		
		//Line Number 101
		for(int k=0; k < sm_smos_25.size(); k++){
			if(sm_smos_25[k][i]< 0.03){
				sm_smos_25[k][i]=NAN;
			}
		}
	}
	
	
}
