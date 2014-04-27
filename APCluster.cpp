//AP聚类
#include<iostream>
#include<vector>
#include<list>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <cfloat>
using namespace std;

//表示一个顶点
struct StructPoint{
	int nID;
	int nType;
	float pointX, pointY;
	bool operator == (const StructPoint &p) const{
		return (nID == p.nID);
	}
};

//表示有起点终点权值的一条有向边
struct StructEdge{
    int srcID;
    int desID;
    double fValue;
    StructEdge(){
        srcID = 0;
        desID = 0;
        fValue = 0.0d;
    }
	bool operator == (const StructEdge &e) const{
		return (srcID == e.srcID && desID == e.desID);
	}
};

//计算两个数据点间的欧氏距离 
bool calEuclidDistance(StructPoint &srcP, StructPoint &desP, double &dis){
	dis = (srcP.pointX - desP.pointX) * (srcP.pointX - desP.pointX) + (srcP.pointY - desP.pointY) * (srcP.pointY - desP.pointY);
	dis = - dis;
	return true;
}

//计算所有数据点间的欧氏距离
bool calSimilarity(list<StructPoint> &listPoint, list<StructEdge> &listSim){
	list<StructPoint>::iterator itSrc, itDes;
	StructEdge edgeTmp;
	for(itSrc = listPoint.begin(); itSrc != listPoint.end(); ++itSrc){
		edgeTmp.srcID = itSrc->nID;
		itDes = itSrc;
		for(++itDes; itDes != listPoint.end(); ++itDes){
			edgeTmp.desID = itDes->nID;
			if(!calEuclidDistance(*itSrc, *itDes, edgeTmp.fValue)){
				cout << "计算数据点间欧氏距离出错！" << endl;
				return false;
			}
			listSim.push_back(edgeTmp);
		}
	}
	return true;
}

//更新r(i,k)和a(i,k) 
bool updateRA(list<StructPoint> &listPoint, list<StructEdge> &listSim, list<StructEdge> &listNewR, list<StructEdge> &listNewA, list<StructEdge> &listOldR, list<StructEdge> &listOldA, float fFactor)
{
	//先更新吸引度R(i,k)=S(i,k)-max{A(i,j)+S(i,j)}
	list<StructPoint>::iterator itSrc, itDes, itTmp;
	list<StructEdge>::iterator itEdge, itEdgeTmp;
	StructEdge edge, edgeTmp;
	list<StructEdge> listAS;
	double dMax = -INT_MAX;
	double dMaxTmp = 0;
	double dSum = 0;
	
	//先计算A(i,j)+S(i,j)
	double dATmp = 0;
	double dSTmp = 0;
	for(itSrc=listPoint.begin(); itSrc!=listPoint.end(); ++itSrc){ //i
		for(itDes=listPoint.begin(); itDes!=listPoint.end(); ++itDes){ //j
			edge.srcID = itSrc->nID;
			edge.desID = itDes->nID;
			dATmp = 0;
			if((itEdge=find(listNewA.begin(),listNewA.end(),edge)) != listNewA.end()){ //找到吸引度
				dATmp = itEdge->fValue; //i对j的归属度,A(i,j)
			}
			dSTmp = 0;
			if(itSrc->nID > itDes->nID){ //因为相似矩阵是对称矩阵，只保留了一半
				edge.srcID = itDes->nID;
				edge.desID = itSrc->nID;
			}
			if((itEdge=find(listSim.begin(),listSim.end(),edge)) != listSim.end()){ //找到相似度
				dSTmp = itEdge->fValue;
			}
			edgeTmp.srcID = itSrc->nID;
			edgeTmp.desID = itDes->nID;
			edgeTmp.fValue = (float)(dATmp + dSTmp); //A(i,j)+S(i,j)
			listAS.push_back(edgeTmp);
		}
	}
	
	//更新R(i,k)
	//先清空R
	listNewR.clear();
	for(itSrc=listPoint.begin(); itSrc!=listPoint.end(); ++itSrc){ //i
		for(itDes=listPoint.begin(); itDes!=listPoint.end(); ++itDes){ //k
			if(itDes->nID == itSrc->nID){ //i等于k,用公式R(k,k)=S(k,k)-max{S(i,k')},k不等于k'
				//确定最大S(i,k')
				dMax = -INT_MAX;
				for(itTmp=listPoint.begin(); itTmp!=listPoint.end(); ++itTmp){ //k'
					if(itTmp->nID == itDes->nID){ //k'等于k
						continue;
					}
					edgeTmp.srcID = itSrc->nID; //i
					edgeTmp.desID = itTmp->nID; //k'
					if(itSrc->nID > itTmp->nID){ //相似矩阵是对称矩阵，只保留一半
						edgeTmp.srcID = itTmp->nID; // k'
						edgeTmp.desID = itSrc->nID; // i
					}
					if((itEdge=find(listSim.begin(),listSim.end(),edgeTmp)) != listSim.end()){ //有相似度
						if(itEdge->fValue > dMax){
							dMax = itEdge->fValue;
						}
					}
				}// end for k'
				//确定S(k,k)
				edgeTmp.srcID = itDes->nID; //k
				edgeTmp.desID = itDes->nID; //k
				edgeTmp.fValue = 0;
				if((itEdge=find(listSim.begin(),listSim.end(),edgeTmp)) != listSim.end()){ //有相似度
					edgeTmp.fValue = itEdge->fValue; //S(k,k)
				}
				edge.srcID = itSrc->nID; //i,此时i等于k
				edge.desID = itDes->nID; //k
				edge.fValue = (float)(edgeTmp.fValue - dMax); //S(k,k) - max{S(i,k')}
				if(fabs(edge.fValue) < 0.0000001){ //为0
					edge.fValue = 0;
				}
				listNewR.push_back(edge); //保存新的吸引度
				continue;
			}//end if k'(i == k)

			//确定最大值 max{A(i,j)+S(i,j)} (i != k)
			dMax = -INT_MAX;
			dMaxTmp = 0;
			for(itEdge=listAS.begin(); itEdge!=listAS.end(); ++itEdge){ //每个A(i,j)+S(i,j)
				if(itEdge->srcID == itSrc->nID && itEdge->desID != itDes->nID){ //k'不等于k时
					if(itEdge->fValue > dMax){
						dMax = itEdge->fValue;
					}
				}
			}
			//确定S(i,k)
			edgeTmp.srcID = itSrc->nID; //i
			edgeTmp.desID = itDes->nID; //k
			edgeTmp.fValue = 0;
			if(itSrc->nID > itDes->nID){
				edgeTmp.srcID = itDes->nID;
				edgeTmp.desID = itSrc->nID;
			}
			if((itEdge=find(listSim.begin(),listSim.end(),edgeTmp)) != listSim.end()){ //有相似度
				edgeTmp.fValue = itEdge->fValue;
								
			}
			edge.srcID = itSrc->nID; //i
			edge.desID = itDes->nID; //k
			edge.fValue = (float)(edgeTmp.fValue - dMax); //S(i,k) - max{A(i,j)+S(i,j)}
			if(fabs(edge.fValue) < 0.0000001){ //为0
				edge.fValue = 0;
			}
			listNewR.push_back(edge); //保存新的吸引度
		} //end for k
	}//end for i
	listAS.clear(); //释放内存
	
	//迭代吸引度
	for(itEdge=listNewR.begin(); itEdge!=listNewR.end(); ++itEdge)
	{  //最终的吸引度 R = (1-阻尼系数)*newR + 阻尼系数*oldR  //保留oldR的作用 
		itEdge->fValue *= (float)(1.0 - fFactor);
		if((itEdgeTmp=find(listOldR.begin(),listOldR.end(),*itEdge)) == listOldR.end())
		{
			continue;
		}
		itEdge->fValue += (itEdgeTmp->fValue * fFactor);
	}
	
	//最后更新A(i,k) = min{0,R(k,k) + sum{max{0,R(j,k)}}
	//先清空A
	listNewA.clear();
	for(itSrc=listPoint.begin(); itSrc!=listPoint.end(); ++itSrc) //i
	{
		for(itDes=listPoint.begin(); itDes!=listPoint.end(); ++itDes) //k
		{
			if(itDes->nID == itSrc->nID) //i等于k时，更新公式：A(k,k)=sum{max{0,R(i',k)},其中，i'不等于k
			{
				dMax = 0;
				for(itTmp=listPoint.begin(); itTmp!=listPoint.end(); ++itTmp) //i'
				{
					if(itTmp->nID == itDes->nID) //i'等于k
					{
						continue;
					}
					edgeTmp.srcID = itTmp->nID; //R(i',k), i'
					edgeTmp.desID = itDes->nID; // k
					if((itEdge=find(listNewR.begin(),listNewR.end(),edgeTmp)) != listNewR.end()) //有R(i',k)
					{
						if(itEdge->fValue > 0) // && itEdge->fValue < 0
						{
							dMax += itEdge->fValue;
						}
					}
				}
			
				edge.srcID = itDes->nID; //A(k,k), k
				edge.desID = itDes->nID; // k
				edge.fValue = (float)(dMax); //A(k,k)=max{0,R(i',k)}
				if(fabs(edge.fValue) < 0.0000001) //为0
				{
					edge.fValue = 0;
				}
				listNewA.push_back(edge); //保存更新的A(k,k)
				continue;
			}

			//当i不等于k时
			//求sum{max{0,R(j,k)}
			dSum = 0;
			edge.desID = itDes->nID; //R(j,k)中的k
			for(itTmp=listPoint.begin(); itTmp!=listPoint.end(); ++itTmp) //j
			{
				if(itTmp->nID == itSrc->nID || itTmp->nID == itDes->nID) //j不能等于i,也不能等于k
				{
					continue;
				}
				edge.srcID = itTmp->nID; //R(j,k)中的j
				edge.fValue = 0;
				if((itEdge=find(listNewR.begin(),listNewR.end(),edge)) != listNewR.end()) //没有R(j,k)
				{
					edge.fValue = (itEdge->fValue>0 ? itEdge->fValue : 0);
				}
				dSum += edge.fValue; //累加max{0,R(j,k)}
			} //end j
		
			edge.srcID = itDes->nID; //R(k,k)中的前一个k，实际上都相同，但链表find时需要
			edge.desID = itDes->nID;
			edge.fValue = 0;
			if((itEdge=find(listNewR.begin(),listNewR.end(),edge)) != listNewR.end()) //找到R(k,k)
			{
				dSum += itEdge->fValue; //累加R(k,k),即R(k,k) + sum{max{0,R(j,k)}
			}
			if(dSum < 0) //R(k,k) + sum{max{0,R(j,k)}大于等于0不用记录
			{
				edge.srcID = itSrc->nID; //A(i,k)中的i
				edge.desID = itDes->nID; //A(i,k)中的k
				edge.fValue = (float)(dSum);
				if(fabs(edge.fValue) < 0.0000001) //为0
				{
					edge.fValue = 0;
				}
				listNewA.push_back(edge); //保存更新的A(i,k)
			}
		} //end k
	}//end i

	//迭代归属度A
	for(itEdge=listNewA.begin(); itEdge!=listNewA.end(); ++itEdge)
	{ //最终的归属度 A = (1-阻尼系数)*newA + 阻尼系数*oldA  //保留oldA的作用
		itEdge->fValue *= (float)(1.0 - fFactor);
		if((itEdgeTmp=find(listOldA.begin(),listOldA.end(),*itEdge)) == listOldA.end())
		{
			continue;
		}
		itEdge->fValue += (itEdgeTmp->fValue * fFactor);
	}

	return true;
}

//确定簇中心
bool centerJudge(list<StructEdge>& listNewR, list<StructEdge>& listNewA, set<int>& setCenter)
{
	list<StructEdge>::iterator itEdgeR, itEdgeA;
	setCenter.clear(); //先清空簇中心点ID集合
	float fTmp = 0;
	for(itEdgeR=listNewR.begin(); itEdgeR!=listNewR.end(); ++itEdgeR)
	{
		if(itEdgeR->srcID != itEdgeR->desID) //只判断R(k,k)
		{
			continue;
		}
		fTmp = 0;
		if((itEdgeA=find(listNewA.begin(),listNewA.end(),*itEdgeR)) != listNewA.end())
		{
			fTmp = itEdgeA->fValue;
		}
		if(fTmp + itEdgeR->fValue > 0.000001) //可以作为簇中心
		{
			setCenter.insert(itEdgeR->srcID); //将簇中心加入集合
		}
	}

	return true;
}

//根据k个聚类中心，划分数据集
bool partitionClustering(list<StructPoint>& listPoint, list<StructPoint>& listKPoint)
{
	if(listPoint.empty() || listKPoint.empty())
	{
		return false;
	}
	list<StructPoint>::iterator itStrPoint;
	list<StructPoint>::iterator itKPoint;
	int j = 0;
	for(itStrPoint=listPoint.begin(); itStrPoint!=listPoint.end(); ++itStrPoint) //每个数据
	{
		//计算与每个聚类的距离，选择聚类最小的作为该数据的类
		double dMin = FLT_MAX;
		int nType = itStrPoint->nType; //数据的类别
		j = 1;
		for(itKPoint=listKPoint.begin(); itKPoint!=listKPoint.end(); ++itKPoint,++j) //每个聚类
		{
			double dTmpSum = 0;
			dTmpSum = (itStrPoint->pointX - itKPoint->pointX) * (itStrPoint->pointX - itKPoint->pointX) + (itStrPoint->pointY - itKPoint->pointY) * (itStrPoint->pointY - itKPoint->pointY);
			dTmpSum = sqrt(dTmpSum); //欧式距离
			if(dTmpSum < dMin)
			{
				dMin = dTmpSum;
				nType = j;
			}
		} //end for 每个聚类中心
		itStrPoint->nType = nType; //确定该数据对象的类别
	} //end for 每个数据

	return true;
}

//AP(Affinity Propagation)聚类,首先计算每对数据对象之间的距离，适用于小数据集，因为大数据集内存不能容纳
bool cluster_AP_int(list<StructPoint>& listPoint, list<StructPoint>& listK, int nIt, int nK, float fFactor)
{  //nIt迭代次数； nK数据点数 
	if(listPoint.empty())
	{
		return false;                       
	}

	list<StructEdge> listSim;
	//第一步：计算每对点之间的距离,内存无法容纳每对数据对象之间的距离，因此不能事先计算并保存，只能实时计算
	cout << "正在计算每对数据对象之间的距离" << endl; 
	if(!(calSimilarity(listPoint,listSim)))
	{
		cout << "计算数据对象两两之间的距离出错" << endl;
		return false;
	}
	//用所有距离的中值的一半作为参考度
	cout << "确定参考度" << endl;
	list<StructEdge>::iterator itEdge;
	vector<float> vecSim;  //用来保存所有sim值，然后排序得到中值 
	for(itEdge=listSim.begin(); itEdge!=listSim.end(); ++itEdge)
	{
		vecSim.push_back(itEdge->fValue);
	}
	sort(vecSim.begin(),vecSim.end()); //排序
	int nSize = vecSim.size();
	double dP = vecSim[nSize/2];
	if(nSize % 2 == 0) //偶数个数据对象
	{
		dP = (vecSim[nSize/2] + vecSim[nSize/2+1]) / 2.0;
	}
	vecSim.clear(); //释放内存
	
	list<StructPoint>::iterator itStrPoint;
	StructEdge edge;
	for(itStrPoint=listPoint.begin(); itStrPoint!=listPoint.end(); ++itStrPoint)
	{   //初始化s(k,k)都为dP 
		edge.srcID = itStrPoint->nID;
		edge.desID = itStrPoint->nID;
		edge.fValue = (float)dP;
		listSim.push_back(edge);
	}

	//第二步：迭代寻找合适的聚类中心
	list<StructEdge> listNewR, listNewA, listOldR, listOldA; //新旧吸引度和归属度
	
	int i = 0;
	set<int> setNewCenter, setOldCenter, setTmp;
	insert_iterator<set<int,less<int> > >  res_ins(setTmp, setTmp.begin());
	int t = 0;

	while(i < nIt)
	{
		cout << "第" << ++i << "次迭代" << endl;
		//保存新吸引度和归属度
		listOldR.clear();
		listOldA.clear();
		listOldR.assign(listNewR.begin(), listNewR.end());
		listOldA.assign(listNewA.begin(), listNewA.end());

		//更新吸引度和归属度
		cout << "正在更新吸引度和归属度" << endl;
		if(!updateRA(listPoint, listSim, listNewR, listNewA, listOldR, listOldA, fFactor))
		{
			cout << "更新吸引度和归属度失败" << endl;
			return false;
		}
		
		//确定中心点
		cout << "正在确定中心点" << endl;
		if(!centerJudge(listNewR, listNewA, setNewCenter))
		{
			cout << "确定中心点失败" << endl;
		}
		//更新簇中心
		setOldCenter.clear();
		setOldCenter.insert(setNewCenter.begin(),setNewCenter.end());
	}
	listSim.clear(); //释放数据对象之间的距离占用的内存
	setOldCenter.clear();

	if(setNewCenter.empty())
	{
		cout << "没有找出任何簇中心" << endl;
		return false;
	}
	//第三步：聚类
	StructPoint strPoint;
	//提取聚类中心
	int k = 1;
	for(set<int>::iterator itSet=setNewCenter.begin(); itSet!=setNewCenter.end(); ++itSet,++k)
	{
		for(itStrPoint=listPoint.begin(); itStrPoint!=listPoint.end(); ++itStrPoint)
		{
			if(*itSet == itStrPoint->nID) //簇中心
			{
				strPoint.nID = *itSet;
				strPoint.nType = k;
				strPoint.pointX = itStrPoint->pointX;
				strPoint.pointY = itStrPoint->pointY;
				listK.push_back(strPoint);
			}
		}
	}
	setNewCenter.clear();

	//聚类
	if(!partitionClustering(listPoint,listK))
	{
		cout<<"聚类失败"<<endl;
	}
	cout<<"簇数目为："<<(k-1)<<endl;

	return true;
}



int main()
{
	//AP算法例子
	list<StructPoint> listPoint;
	list<StructPoint> listK; //簇中心
	ifstream inFile;
	inFile.open("./AP_example.txt");
	StructPoint strPoint;
	if(!inFile)
	{
		cout << "读取AP数据例子失败" << endl;
		return 0;
	}
	int k = 0;
	while(!inFile.eof())
	{
		strPoint.nID = ++k;
		strPoint.nType = 0;
		inFile >> strPoint.pointX;
		inFile >> strPoint.pointY;
		listPoint.push_back(strPoint);
	}
	inFile.close();

	//聚类
	cout<<"正在聚类"<<endl;
	float fFactor = 0.5;  //阻尼系数为0.5 
	int nIt = 500;        //迭代次数 
	int nK = 10;		  //数据点数 
	if(!(cluster_AP_int(listPoint, listK, nIt, nK, fFactor))) //AP聚类
	{
		cout<<"聚类失败"<<endl;
		//释放内存
		listPoint.clear();
		return 0;
	}

	cout<<"聚类成功"<<endl;

	//输出聚类结果
	ofstream outFile("./cluster_AP.txt"); 
	if(!outFile)
	{
		cout << "保存失败" << endl;
	}
	for(list<StructPoint>::iterator k=listPoint.begin(); k!=listPoint.end(); ++k)
	{
		outFile << k->nID << "\t" << k->nType << endl;
	}
	outFile.close();

	cout<<"簇中心为："<<endl;
	for(list<StructPoint>::iterator k=listK.begin(); k!=listK.end(); ++k)
	{
		cout << k->nID << "\t" << k->nType << "\t" << k->pointX << "\t" << k->pointY;
		cout<<endl;
	}
	//释放内存
	
	listPoint.clear();
	listK.clear();

	return 0;
}
