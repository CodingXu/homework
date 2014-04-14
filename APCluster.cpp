//AP聚类
#include<iostream>
#include<vector>
#include <list>
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

//表示有起点终点权值的一条边
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
	dis = (srcP.x - desP.x) * (srcP.x - desP.x) + (srcP.y - desP.y) * (srcP.y - desP.y);
	dis = - dis;
	return true;
}

//计算所有数据点间的欧氏距离
bool calSimilarity(list<StructPoint> &listPoint, list<StructEdge> &listSim){
	list<StructPoint>::iterator itSrc, itDes;
	StructEdge tempEdge;
	for(itSrc = listPoint.begin(); itSrc != listPoint.end(); ++itSrc){
		tempEdge.srcID = itSrc->nID;
		itDes = itSrc;
		for(++itDes; itDes != listPoint.end(); ++itDes){
			tempEdge.desID = itDes->nID;
			if(!calEuclidDistance(*itSrc, *itDes, tempEdge.fValue)){
				cout << "计算数据点间欧氏距离出错！" << endl;
				return false;
			}
			listSim.push_back(tempEdge);
		}
	}
	return true;
}

//


int main(){
	
	return 0;
}
