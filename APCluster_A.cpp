
#include<iostream>
#include<fstream>
#include<algorithm> 
#include<vector>
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
             {  //相似度矩阵轴对称，因此只计算一半 
                     sim[i][j] = (point[i].pX - point[j].pX) * (point[i].pX - point[j].pX)
                                 + (point[i].pY - point[j].pY) * (point[i].pY - point[j].pY);
                     sim[i][j] = -sim[i][j];
                     vecSim.push_back(sim[i][j]); //值加入vector中，为了排序 
             }
     }
     sort(vecSim.begin(), vecSim.end()); //排序，为了找中值
     int nSize = vecSim.size();
     double dP = vecSim[nSize/2];  //dP为相似度中值 
     if(nSize % 2 == 0) //偶数个数据对象
	 {
              dP = (vecSim[nSize/2] + vecSim[nSize/2+1]) / 2.0;
	 }
	 for(int i=0; i!= dataNum; ++i)
	 {
             sim[i][i] = dP;  //sim[i][i]初始化为相似度中值 
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
                     //计算所有 A(i,j)+S(i,j)
                     tmpAS = i>j?sim[j][i]:sim[i][j];
                     tmpAS += newA[i][j];
                     sumAS[i][j] = tmpAS;
             }
     }
     //先更新吸引度R(i,k)
     for(int i=0; i!=dataNum; ++i) //i
     {
             for(int k=0; k!=dataNum; ++k) //k
             {
	         		double maxS_AS;
					if(i == k) //i等于k,用公式R(k,k)=S(k,k)-max{S(i,j)},k不等于j
	                {
	                	maxS_AS = -INT_MAX;
					 	for(int j=0; j!=dataNum; ++j)
					 	{
					 		if(k != j)
					 		{
					 			int tmpS = i>j?sim[j][i]:sim[i][j];
								if(tmpS > maxS_AS)
								{
									maxS_AS = tmpS;
								} 
					 		}
					 	}
					}
					else //i不等于k,用公式R(i,k)=S(i,k)-max{A(i,j)+S(i,j)},k不等于j
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
					int tmpS = i>k?sim[k][i]:sim[i][k];
					newR[i][k] = tmpS - maxS_AS;
             } //end for k
     } //end for i
     return true;
}

int main()
{
    ifstream inFile;
    inFile.open("../AP_example.txt");
    if(!inFile)
    {
		cout << "读取AP数据例子失败" << endl;
		return 0;
	}
	dataNum = 0;
	while(!inFile.eof())
	{
		point[dataNum].nID = dataNum;
		point[dataNum].nType = 0;
		inFile >> point[dataNum].pX;
		inFile >> point[dataNum].pY;
		dataNum++;
	}
	inFile.close();
	
	calSimilarity();
    return 0;
}
