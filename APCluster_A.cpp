
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
                     tmpAS = (i > j) ? sim[j][i] : sim[i][j];
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
					 			double tmpS = (i > j) ? sim[j][i] : sim[i][j];
								if(tmpS > maxS_AS)
								{
									maxS_AS = tmpS;
								} 
					 		}
					 	}
					 	//cout<<i<<" : "<<maxS_AS<<endl;
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
					double tmpS = (i > k) ? sim[k][i] : sim[i][k];
					newR[i][k] = tmpS - maxS_AS;
             } //end for k
     } //end for i
     
    //迭代吸引度，完成R[i][j]的最终更新 
	for(int i=0; i!=dataNum; ++i)
	{  //最终的吸引度 R = (1-阻尼系数)*newR + 阻尼系数*oldR  //保留oldR的作用 
		for(int j=0; j!=dataNum; ++j)
		{
			if(fabs(newR[i][j]) < 0.0000001)
			{
				newR[i][j] = 0;
			}
			newR[i][j] = (1.0 - fFactor) * newR[i][j] + fFactor * oldR[i][j];
		}
	}
	
	//更新A(i,k)
	for(int i=0; i!=dataNum; ++i) //i
	{
		for(int k=0; k!=dataNum; ++k) //k
		{
			//先计算sum{max{0,R(j,k)} 
			double sumTmp = 0.0;
			for(int j=0; j!=dataNum; ++j)
			{
				if(k != j && i!= j)
				{
					sumTmp += (0.0 > newR[j][k]) ? 0.0 : newR[j][k]; 
				}
			}
			if(i == k) //i等于k，更新公式：A(k,k)=sum{max{0,R(j,k)},其中，j不等于k 
			{
				newA[i][k] = sumTmp;
			}
			else //i不等于k，更新公式：A(i,k) = min{0,R(k,k) + sum{max{0,R(j,k)}}
			{
				sumTmp += newR[k][k];
				newA[i][k] = (0.0 < sumTmp) ? 0.0 : sumTmp;
			}
		} //end for k
	} // end for i
	
	//迭代吸引度，完成A[i][j]的最终更新 
	for(int i=0; i!=dataNum; ++i)
	{  //最终的吸引度 A = (1-阻尼系数)*newA + 阻尼系数*oldA  //保留oldA的作用 
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

//确定簇中心
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

//根据k个聚类中心，划分数据集
bool partitionClustering(int (&center)[DATASIZE], int cenNum)
{
	if(cenNum <= 0)
	{
		return false;
	}
	for(int i=0; i!=dataNum; ++i) //每个数据点 
	{
		double dMin = FLT_MAX;
		int cen_I;
		for(int j=0; j!=cenNum; ++j) //每个簇中心 
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
	//先计算数据点间相似度 
	calSimilarity();
	
	//开始迭代更新R和A 
	int i = 0;
	while(i < nIt)
	{
		cout << "第" << ++i << "次迭代" << endl;
		//更新吸引度和归属度
		cout << "正在更新吸引度和归属度" << endl;
		if(!updateRA(fFactor))
		{
			cout << "更新吸引度和归属度失败" << endl;
			return false;
		}
	}
	
	//更新簇中心 
	cenNum = centerJudge(center);
	if(cenNum <= 0)
	{
		cout << "没有找出簇中心" << endl;
		return false;
	}
	
	//聚类
	if(!partitionClustering(center, cenNum))
	{
		return false;
	}
	cout << "簇数目为：" << cenNum << endl;
	
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
		point[dataNum].nID = dataNum; //数据点标号和在数组对应下标相同 
		point[dataNum].nType = 0;
		inFile >> point[dataNum].pX;
		inFile >> point[dataNum].pY;
		dataNum++;
	}
	inFile.close();
	
	//聚类
	cout << "正在聚类" << endl;
	float fFactor = 0.5;  //阻尼系数为0.5 
	int nIt = 500;        //迭代次数 
	int center[DATASIZE];
	int cenNum = 0;
	if(!(cluster_AP(nIt, fFactor, center, cenNum))) //AP聚类
	{
		cout << "聚类失败" << endl;
		return 0;
	}
	cout << "聚类成功" << endl;
	
	//输出聚类结果
	ofstream outFile("./cluster_AP.txt"); 
	if(!outFile)
	{
		cout << "保存失败" << endl;
	}
	outFile << "数据点" << "\t" << "中心点" << endl;
	for(int i=0; i!=dataNum; ++i)
	{
		outFile << point[i].nID << "\t" << point[i].nType << endl;
	}
	outFile.close();

	cout << "簇中心为：" << endl;
	for(int i=0; i!=cenNum; ++i)
	{
		cout << point[center[i]].nID << "\t" << point[center[i]].pX << "\t" << point[center[i]].pY;
		cout << endl;
	}
	
    return 0;
}
