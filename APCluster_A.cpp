//APCluster
#include<iostream>
#include<fstream>
#include<algorithm> 
#include <cmath>
#include<vector>
#include<cfloat>
#define DATASIZE 100
using namespace std;

static int dataNum;
double sim[DATASIZE][DATASIZE] = {0};
double oldR[DATASIZE][DATASIZE] = {0};
double oldA[DATASIZE][DATASIZE] = {0};
double newR[DATASIZE][DATASIZE] = {0};
double newA[DATASIZE][DATASIZE] = {0};
struct Point{
       int nID;
       int nType;
       float pX, pY;
}point[DATASIZE];

bool calSimilarity()
{
     vector<double> vecSim;
     for(int i=0; i!=dataNum; ++i)
     {
             for(int j=i+1; j!=dataNum; ++j)
             {  //The similarity matrix is axial symmetry, so only to calculate the half 
                     sim[i][j] = (point[i].pX - point[j].pX) * (point[i].pX - point[j].pX)
                                 + (point[i].pY - point[j].pY) * (point[i].pY - point[j].pY);
                     sim[i][j] = -sim[i][j];
                     vecSim.push_back(sim[i][j]); //the value is added in vector, in order to sort 
             }
     }
     sort(vecSim.begin(), vecSim.end()); //Sort, in order to find the median
     int nSize = vecSim.size();
     double dP = vecSim[nSize/2];  //dP as the similarity value
     if(nSize % 2 == 0) //even data objects
	 {
              dP = (vecSim[nSize/2] + vecSim[nSize/2+1]) / 2.0;
	 }
	 for(int i=0; i!= dataNum; ++i)
	 {
             sim[i][i] = dP;  //sim[i][i] is initialized to the similarity value 
     }

     return true;
}

bool updateRA(float fFactor)
{
     double sumAS[DATASIZE][DATASIZE];
     double tmpAS;
     for(int i=0; i!=dataNum; ++i)
     {
             for(int j=0; j!=dataNum; ++j)
             {
                     oldR[i][j] = newR[i][j];
                     oldA[i][j] = newA[i][j];
                     //Calculate all A(i,j)+S(i,j)
                     tmpAS = (i > j) ? sim[j][i] : sim[i][j];
                     tmpAS += newA[i][j];
                     sumAS[i][j] = tmpAS;
             }
     }
     //Update R(i,k) first
     for(int i=0; i!=dataNum; ++i) //i
     {
             for(int k=0; k!=dataNum; ++k) //k
             {
	         		double maxS_AS;
					if(i == k) //i is equal to k,Using the formula R(k,k)=S(k,k)-max{S(i,j)},k is not equal to j
	                {
	                	maxS_AS = -INT_MAX;
					 	for(int j=0; j!=dataNum; ++j)
					 	{
					 		if(k != j)
					 		{
					 			double tmpS = (i > j) ? sim[j][i] : sim[i][j];
								if(tmpS > maxS_AS)
								{
									maxS_AS = tmpS;
								} 
					 		}
					 	}
					 	//cout<<i<<" : "<<maxS_AS<<endl;
					}
					else //i is not equal to k,Using the formula R(i,k)=S(i,k)-max{A(i,j)+S(i,j)},k is not equal to j
					{
						maxS_AS = -INT_MAX;
						for(int j=0; j!=dataNum; ++j)
					 	{
					 		if(k != j)
					 		{
					 			if(sumAS[i][j] > maxS_AS)
								{
									maxS_AS = sumAS[i][j];
								} 
					 		}
					 	}
					} 
					double tmpS = (i > k) ? sim[k][i] : sim[i][k];
					newR[i][k] = tmpS - maxS_AS;
             } //end for k
     } //end for i
     
    //Iterative attraction degree, complete the final update to R[i][j] 
	for(int i=0; i!=dataNum; ++i)
	{  //The final degree of attraction R = (1-damping factor)*newR + damping factor*oldR  //Retain the effect for oldR
		for(int j=0; j!=dataNum; ++j)
		{
			if(fabs(newR[i][j]) < 0.0000001)
			{
				newR[i][j] = 0;
			}
			newR[i][j] = (1.0 - fFactor) * newR[i][j] + fFactor * oldR[i][j];
		}
	}
	
	//update A(i,k)
	for(int i=0; i!=dataNum; ++i) //i
	{
		for(int k=0; k!=dataNum; ++k) //k
		{
			//Firstly calculate sum{max{0,R(j,k)} 
			double sumTmp = 0.0;
			for(int j=0; j!=dataNum; ++j)
			{
				if(k != j && i!= j)
				{
					sumTmp += (0.0 > newR[j][k]) ? 0.0 : newR[j][k]; 
				}
			}
			if(i == k) //i is equal to k，Update formula：A(k,k)=sum{max{0,R(j,k)}. j is not equal to k
			{
				newA[i][k] = sumTmp;
			}
			else //i is not equal to k，Update formula：A(i,k) = min{0,R(k,k) + sum{max{0,R(j,k)}}
			{
				sumTmp += newR[k][k];
				newA[i][k] = (0.0 < sumTmp) ? 0.0 : sumTmp;
			}
		} //end for k
	} // end for i
	
	//Iterative attraction degree，complete the final update to  A[i][j]
	for(int i=0; i!=dataNum; ++i)
	{  //The final degree of attraction A = (1-damping coefficient)*newA + damping coefficient*oldA  //Retain the  effect for oldA 
		for(int j=0; j!=dataNum; ++j)
		{
			if(fabs(newA[i][j]) < 0.0000001)
			{
				newA[i][j] = 0.0;
			}
			newA[i][j] = (1.0 - fFactor) * newA[i][j] + fFactor * oldA[i][j];
		}
	}

    return true;
}

//To determine the cluster center
int centerJudge(int (&center)[DATASIZE])
{
	int cenNum = 0;
	for(int k=0; k!=dataNum; ++k)
	{
		if(newA[k][k] + newR[k][k] > 0)
		{
			center[cenNum] = k;
			cenNum++;
		}
	}
	return cenNum;
}

//According to the K cluster center, data partitioning
bool partitionClustering(int (&center)[DATASIZE], int cenNum)
{
	if(cenNum <= 0)
	{
		return false;
	}
	for(int i=0; i!=dataNum; ++i) //For each data point
	{
		double dMin = FLT_MAX;
		int cen_I;
		for(int j=0; j!=cenNum; ++j) //For each cluster center
		{
			double dis;
			dis = (point[i].pX - point[center[j]].pX) * (point[i].pX - point[center[j]].pX)
				  + (point[i].pY - point[center[j]].pY) * (point[i].pY - point[center[j]].pY);
			if(dis < dMin)
			{
				dMin = dis;
				cen_I = center[j];
			}
		}
		point[i].nType = cen_I;
	}
	return true;
}

bool cluster_AP(int nIt, float fFactor, int (&center)[DATASIZE], int &cenNum)
{
	//To calculate the similarity between data points 
	calSimilarity();
	
	//To update R and A  
	int i = 0;
	while(i < nIt)
	{
		cout << " Time " << ++i << " iteration." << endl;
		//Update R and A
		cout << "Updating R and A..." << endl;
		if(!updateRA(fFactor))
		{
			cout << "Failed to update R and A !" << endl;
			return false;
		}
	}
	
	//Update the cluster center 
	cenNum = centerJudge(center);
	if(cenNum <= 0)
	{
		cout << "Do not found the cluster !" << endl;
		return false;
	}
	
	//Clustering
	if(!partitionClustering(center, cenNum))
	{
		return false;
	}
	cout << "Cluster number：" << cenNum << endl;
	
	return true;
}

int main()
{
    ifstream inFile;
    inFile.open("./AP_example.txt");
    if(!inFile)
    {
		cout << "Read AP example data failed !" << endl;
		return 0;
	}
	dataNum = 0;
	while(!inFile.eof())
	{
		point[dataNum].nID = dataNum; //The data point label and in the same array corresponding subscript
		point[dataNum].nType = 0;
		inFile >> point[dataNum].pX;
		inFile >> point[dataNum].pY;
		dataNum++;
	}
	inFile.close();
	
	//Clustering
	cout << "Clustering..." << endl;
	float fFactor = 0.5;  //Damping factor is 0.5 
	int nIt = 500;        //Number of iterations 
	int center[DATASIZE];
	int cenNum = 0;
	if(!(cluster_AP(nIt, fFactor, center, cenNum))) //AP Clustering
	{
		cout << "Cluster failed !" << endl;
		return 0;
	}
	cout << "Cluster success !" << endl;
	
	//output the result of the clustering
	ofstream outFile("./cluster_AP.txt"); 
	if(!outFile)
	{
		cout << "Save failed !" << endl;
	}
	outFile << "Data points" << "\t" << "Center point" << endl;
	for(int i=0; i!=dataNum; ++i)
	{
		outFile << point[i].nID << "\t" << point[i].nType + 1 << endl;
	}
	outFile.close();

	cout << "Cluster center：" << endl;
	for(int i=0; i!=cenNum; ++i)
	{
		cout << point[center[i]].nID + 1 << "\t" << point[center[i]].pX << "\t" << point[center[i]].pY;
		cout << endl;
	}
	
    return 0;
}
