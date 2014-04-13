//AP聚类
#include<iostream>
#include<vector>
using namespace std;

//表示一个顶点
struct StructPoint{
	int id;
	int type;
	float pointX, pointY;
	bool operator == (const StructPoint &a) const{
		return (id == a.id);
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
    StructEdge(int nIDs, int nIDd, float fV)
	{
		srcID = nIDs;
		desID = nIDd;
		fValue = fV;
	}
	bool operator == (const StructEdge &e) const{
		return (srcID == e.srcID && desID == e.desID);
	}
};
int main(){
	
	return 0;
}
