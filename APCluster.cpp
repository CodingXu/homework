//AP����
#include<iostream>
#include<vector>
#include<list>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <cfloat>
using namespace std;

//��ʾһ������
struct StructPoint{
	int nID;
	int nType;
	float pointX, pointY;
	bool operator == (const StructPoint &p) const{
		return (nID == p.nID);
	}
};

//��ʾ������յ�Ȩֵ��һ�������
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

//�����������ݵ���ŷ�Ͼ��� 
bool calEuclidDistance(StructPoint &srcP, StructPoint &desP, double &dis){
	dis = (srcP.pointX - desP.pointX) * (srcP.pointX - desP.pointX) + (srcP.pointY - desP.pointY) * (srcP.pointY - desP.pointY);
	dis = - dis;
	return true;
}

//�����������ݵ���ŷ�Ͼ���
bool calSimilarity(list<StructPoint> &listPoint, list<StructEdge> &listSim){
	list<StructPoint>::iterator itSrc, itDes;
	StructEdge edgeTmp;
	for(itSrc = listPoint.begin(); itSrc != listPoint.end(); ++itSrc){
		edgeTmp.srcID = itSrc->nID;
		itDes = itSrc;
		for(++itDes; itDes != listPoint.end(); ++itDes){
			edgeTmp.desID = itDes->nID;
			if(!calEuclidDistance(*itSrc, *itDes, edgeTmp.fValue)){
				cout << "�������ݵ��ŷ�Ͼ������" << endl;
				return false;
			}
			listSim.push_back(edgeTmp);
		}
	}
	return true;
}

//����r(i,k)��a(i,k) 
bool updateRA(list<StructPoint> &listPoint, list<StructEdge> &listSim, list<StructEdge> &listNewR, list<StructEdge> &listNewA, list<StructEdge> &listOldR, list<StructEdge> &listOldA, float fFactor)
{
	//�ȸ���������R(i,k)=S(i,k)-max{A(i,j)+S(i,j)}
	list<StructPoint>::iterator itSrc, itDes, itTmp;
	list<StructEdge>::iterator itEdge, itEdgeTmp;
	StructEdge edge, edgeTmp;
	list<StructEdge> listAS;
	double dMax = -INT_MAX;
	double dMaxTmp = 0;
	double dSum = 0;
	
	//�ȼ���A(i,j)+S(i,j)
	double dATmp = 0;
	double dSTmp = 0;
	for(itSrc=listPoint.begin(); itSrc!=listPoint.end(); ++itSrc){ //i
		for(itDes=listPoint.begin(); itDes!=listPoint.end(); ++itDes){ //j
			edge.srcID = itSrc->nID;
			edge.desID = itDes->nID;
			dATmp = 0;
			if((itEdge=find(listNewA.begin(),listNewA.end(),edge)) != listNewA.end()){ //�ҵ�������
				dATmp = itEdge->fValue; //i��j�Ĺ�����,A(i,j)
			}
			dSTmp = 0;
			if(itSrc->nID > itDes->nID){ //��Ϊ���ƾ����ǶԳƾ���ֻ������һ��
				edge.srcID = itDes->nID;
				edge.desID = itSrc->nID;
			}
			if((itEdge=find(listSim.begin(),listSim.end(),edge)) != listSim.end()){ //�ҵ����ƶ�
				dSTmp = itEdge->fValue;
			}
			edgeTmp.srcID = itSrc->nID;
			edgeTmp.desID = itDes->nID;
			edgeTmp.fValue = (float)(dATmp + dSTmp); //A(i,j)+S(i,j)
			listAS.push_back(edgeTmp);
		}
	}
	
	//����R(i,k)
	//�����R
	listNewR.clear();
	for(itSrc=listPoint.begin(); itSrc!=listPoint.end(); ++itSrc){ //i
		for(itDes=listPoint.begin(); itDes!=listPoint.end(); ++itDes){ //k
			if(itDes->nID == itSrc->nID){ //i����k,�ù�ʽR(k,k)=S(k,k)-max{S(i,k')},k������k'
				//ȷ�����S(i,k')
				dMax = -INT_MAX;
				for(itTmp=listPoint.begin(); itTmp!=listPoint.end(); ++itTmp){ //k'
					if(itTmp->nID == itDes->nID){ //k'����k
						continue;
					}
					edgeTmp.srcID = itSrc->nID; //i
					edgeTmp.desID = itTmp->nID; //k'
					if(itSrc->nID > itTmp->nID){ //���ƾ����ǶԳƾ���ֻ����һ��
						edgeTmp.srcID = itTmp->nID; // k'
						edgeTmp.desID = itSrc->nID; // i
					}
					if((itEdge=find(listSim.begin(),listSim.end(),edgeTmp)) != listSim.end()){ //�����ƶ�
						if(itEdge->fValue > dMax){
							dMax = itEdge->fValue;
						}
					}
				}// end for k'
				//ȷ��S(k,k)
				edgeTmp.srcID = itDes->nID; //k
				edgeTmp.desID = itDes->nID; //k
				edgeTmp.fValue = 0;
				if((itEdge=find(listSim.begin(),listSim.end(),edgeTmp)) != listSim.end()){ //�����ƶ�
					edgeTmp.fValue = itEdge->fValue; //S(k,k)
				}
				edge.srcID = itSrc->nID; //i,��ʱi����k
				edge.desID = itDes->nID; //k
				edge.fValue = (float)(edgeTmp.fValue - dMax); //S(k,k) - max{S(i,k')}
				if(fabs(edge.fValue) < 0.0000001){ //Ϊ0
					edge.fValue = 0;
				}
				listNewR.push_back(edge); //�����µ�������
				continue;
			}//end if k'(i == k)

			//ȷ�����ֵ max{A(i,j)+S(i,j)} (i != k)
			dMax = -INT_MAX;
			dMaxTmp = 0;
			for(itEdge=listAS.begin(); itEdge!=listAS.end(); ++itEdge){ //ÿ��A(i,j)+S(i,j)
				if(itEdge->srcID == itSrc->nID && itEdge->desID != itDes->nID){ //k'������kʱ
					if(itEdge->fValue > dMax){
						dMax = itEdge->fValue;
					}
				}
			}
			//ȷ��S(i,k)
			edgeTmp.srcID = itSrc->nID; //i
			edgeTmp.desID = itDes->nID; //k
			edgeTmp.fValue = 0;
			if(itSrc->nID > itDes->nID){
				edgeTmp.srcID = itDes->nID;
				edgeTmp.desID = itSrc->nID;
			}
			if((itEdge=find(listSim.begin(),listSim.end(),edgeTmp)) != listSim.end()){ //�����ƶ�
				edgeTmp.fValue = itEdge->fValue;
								
			}
			edge.srcID = itSrc->nID; //i
			edge.desID = itDes->nID; //k
			edge.fValue = (float)(edgeTmp.fValue - dMax); //S(i,k) - max{A(i,j)+S(i,j)}
			if(fabs(edge.fValue) < 0.0000001){ //Ϊ0
				edge.fValue = 0;
			}
			listNewR.push_back(edge); //�����µ�������
		} //end for k
	}//end for i
	listAS.clear(); //�ͷ��ڴ�
	
	//����������
	for(itEdge=listNewR.begin(); itEdge!=listNewR.end(); ++itEdge)
	{  //���յ������� R = (1-����ϵ��)*newR + ����ϵ��*oldR  //����oldR������ 
		itEdge->fValue *= (float)(1.0 - fFactor);
		if((itEdgeTmp=find(listOldR.begin(),listOldR.end(),*itEdge)) == listOldR.end())
		{
			continue;
		}
		itEdge->fValue += (itEdgeTmp->fValue * fFactor);
	}
	
	//������A(i,k) = min{0,R(k,k) + sum{max{0,R(j,k)}}
	//�����A
	listNewA.clear();
	for(itSrc=listPoint.begin(); itSrc!=listPoint.end(); ++itSrc) //i
	{
		for(itDes=listPoint.begin(); itDes!=listPoint.end(); ++itDes) //k
		{
			if(itDes->nID == itSrc->nID) //i����kʱ�����¹�ʽ��A(k,k)=sum{max{0,R(i',k)},���У�i'������k
			{
				dMax = 0;
				for(itTmp=listPoint.begin(); itTmp!=listPoint.end(); ++itTmp) //i'
				{
					if(itTmp->nID == itDes->nID) //i'����k
					{
						continue;
					}
					edgeTmp.srcID = itTmp->nID; //R(i',k), i'
					edgeTmp.desID = itDes->nID; // k
					if((itEdge=find(listNewR.begin(),listNewR.end(),edgeTmp)) != listNewR.end()) //��R(i',k)
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
				if(fabs(edge.fValue) < 0.0000001) //Ϊ0
				{
					edge.fValue = 0;
				}
				listNewA.push_back(edge); //������µ�A(k,k)
				continue;
			}

			//��i������kʱ
			//��sum{max{0,R(j,k)}
			dSum = 0;
			edge.desID = itDes->nID; //R(j,k)�е�k
			for(itTmp=listPoint.begin(); itTmp!=listPoint.end(); ++itTmp) //j
			{
				if(itTmp->nID == itSrc->nID || itTmp->nID == itDes->nID) //j���ܵ���i,Ҳ���ܵ���k
				{
					continue;
				}
				edge.srcID = itTmp->nID; //R(j,k)�е�j
				edge.fValue = 0;
				if((itEdge=find(listNewR.begin(),listNewR.end(),edge)) != listNewR.end()) //û��R(j,k)
				{
					edge.fValue = (itEdge->fValue>0 ? itEdge->fValue : 0);
				}
				dSum += edge.fValue; //�ۼ�max{0,R(j,k)}
			} //end j
		
			edge.srcID = itDes->nID; //R(k,k)�е�ǰһ��k��ʵ���϶���ͬ��������findʱ��Ҫ
			edge.desID = itDes->nID;
			edge.fValue = 0;
			if((itEdge=find(listNewR.begin(),listNewR.end(),edge)) != listNewR.end()) //�ҵ�R(k,k)
			{
				dSum += itEdge->fValue; //�ۼ�R(k,k),��R(k,k) + sum{max{0,R(j,k)}
			}
			if(dSum < 0) //R(k,k) + sum{max{0,R(j,k)}���ڵ���0���ü�¼
			{
				edge.srcID = itSrc->nID; //A(i,k)�е�i
				edge.desID = itDes->nID; //A(i,k)�е�k
				edge.fValue = (float)(dSum);
				if(fabs(edge.fValue) < 0.0000001) //Ϊ0
				{
					edge.fValue = 0;
				}
				listNewA.push_back(edge); //������µ�A(i,k)
			}
		} //end k
	}//end i

	//����������A
	for(itEdge=listNewA.begin(); itEdge!=listNewA.end(); ++itEdge)
	{ //���յĹ����� A = (1-����ϵ��)*newA + ����ϵ��*oldA  //����oldA������
		itEdge->fValue *= (float)(1.0 - fFactor);
		if((itEdgeTmp=find(listOldA.begin(),listOldA.end(),*itEdge)) == listOldA.end())
		{
			continue;
		}
		itEdge->fValue += (itEdgeTmp->fValue * fFactor);
	}

	return true;
}

//ȷ��������
bool centerJudge(list<StructEdge>& listNewR, list<StructEdge>& listNewA, set<int>& setCenter)
{
	list<StructEdge>::iterator itEdgeR, itEdgeA;
	setCenter.clear(); //����մ����ĵ�ID����
	float fTmp = 0;
	for(itEdgeR=listNewR.begin(); itEdgeR!=listNewR.end(); ++itEdgeR)
	{
		if(itEdgeR->srcID != itEdgeR->desID) //ֻ�ж�R(k,k)
		{
			continue;
		}
		fTmp = 0;
		if((itEdgeA=find(listNewA.begin(),listNewA.end(),*itEdgeR)) != listNewA.end())
		{
			fTmp = itEdgeA->fValue;
		}
		if(fTmp + itEdgeR->fValue > 0.000001) //������Ϊ������
		{
			setCenter.insert(itEdgeR->srcID); //�������ļ��뼯��
		}
	}

	return true;
}

//����k���������ģ��������ݼ�
bool partitionClustering(list<StructPoint>& listPoint, list<StructPoint>& listKPoint)
{
	if(listPoint.empty() || listKPoint.empty())
	{
		return false;
	}
	list<StructPoint>::iterator itStrPoint;
	list<StructPoint>::iterator itKPoint;
	int j = 0;
	for(itStrPoint=listPoint.begin(); itStrPoint!=listPoint.end(); ++itStrPoint) //ÿ������
	{
		//������ÿ������ľ��룬ѡ�������С����Ϊ�����ݵ���
		double dMin = FLT_MAX;
		int nType = itStrPoint->nType; //���ݵ����
		j = 1;
		for(itKPoint=listKPoint.begin(); itKPoint!=listKPoint.end(); ++itKPoint,++j) //ÿ������
		{
			double dTmpSum = 0;
			dTmpSum = (itStrPoint->pointX - itKPoint->pointX) * (itStrPoint->pointX - itKPoint->pointX) + (itStrPoint->pointY - itKPoint->pointY) * (itStrPoint->pointY - itKPoint->pointY);
			dTmpSum = sqrt(dTmpSum); //ŷʽ����
			if(dTmpSum < dMin)
			{
				dMin = dTmpSum;
				nType = j;
			}
		} //end for ÿ����������
		itStrPoint->nType = nType; //ȷ�������ݶ�������
	} //end for ÿ������

	return true;
}

//AP(Affinity Propagation)����,���ȼ���ÿ�����ݶ���֮��ľ��룬������С���ݼ�����Ϊ�����ݼ��ڴ治������
bool cluster_AP_int(list<StructPoint>& listPoint, list<StructPoint>& listK, int nIt, int nK, float fFactor)
{  //nIt���������� nK���ݵ��� 
	if(listPoint.empty())
	{
		return false;                       
	}

	list<StructEdge> listSim;
	//��һ��������ÿ�Ե�֮��ľ���,�ڴ��޷�����ÿ�����ݶ���֮��ľ��룬��˲������ȼ��㲢���棬ֻ��ʵʱ����
	cout << "���ڼ���ÿ�����ݶ���֮��ľ���" << endl; 
	if(!(calSimilarity(listPoint,listSim)))
	{
		cout << "�������ݶ�������֮��ľ������" << endl;
		return false;
	}
	//�����о������ֵ��һ����Ϊ�ο���
	cout << "ȷ���ο���" << endl;
	list<StructEdge>::iterator itEdge;
	vector<float> vecSim;  //������������simֵ��Ȼ������õ���ֵ 
	for(itEdge=listSim.begin(); itEdge!=listSim.end(); ++itEdge)
	{
		vecSim.push_back(itEdge->fValue);
	}
	sort(vecSim.begin(),vecSim.end()); //����
	int nSize = vecSim.size();
	double dP = vecSim[nSize/2];
	if(nSize % 2 == 0) //ż�������ݶ���
	{
		dP = (vecSim[nSize/2] + vecSim[nSize/2+1]) / 2.0;
	}
	vecSim.clear(); //�ͷ��ڴ�
	
	list<StructPoint>::iterator itStrPoint;
	StructEdge edge;
	for(itStrPoint=listPoint.begin(); itStrPoint!=listPoint.end(); ++itStrPoint)
	{   //��ʼ��s(k,k)��ΪdP 
		edge.srcID = itStrPoint->nID;
		edge.desID = itStrPoint->nID;
		edge.fValue = (float)dP;
		listSim.push_back(edge);
	}

	//�ڶ���������Ѱ�Һ��ʵľ�������
	list<StructEdge> listNewR, listNewA, listOldR, listOldA; //�¾������Ⱥ͹�����
	
	int i = 0;
	set<int> setNewCenter, setOldCenter, setTmp;
	insert_iterator<set<int,less<int> > >  res_ins(setTmp, setTmp.begin());
	int t = 0;

	while(i < nIt)
	{
		cout << "��" << ++i << "�ε���" << endl;
		//�����������Ⱥ͹�����
		listOldR.clear();
		listOldA.clear();
		listOldR.assign(listNewR.begin(), listNewR.end());
		listOldA.assign(listNewA.begin(), listNewA.end());

		//���������Ⱥ͹�����
		cout << "���ڸ��������Ⱥ͹�����" << endl;
		if(!updateRA(listPoint, listSim, listNewR, listNewA, listOldR, listOldA, fFactor))
		{
			cout << "���������Ⱥ͹�����ʧ��" << endl;
			return false;
		}
		
		//ȷ�����ĵ�
		cout << "����ȷ�����ĵ�" << endl;
		if(!centerJudge(listNewR, listNewA, setNewCenter))
		{
			cout << "ȷ�����ĵ�ʧ��" << endl;
		}
		//���´�����
		setOldCenter.clear();
		setOldCenter.insert(setNewCenter.begin(),setNewCenter.end());
	}
	listSim.clear(); //�ͷ����ݶ���֮��ľ���ռ�õ��ڴ�
	setOldCenter.clear();

	if(setNewCenter.empty())
	{
		cout << "û���ҳ��κδ�����" << endl;
		return false;
	}
	//������������
	StructPoint strPoint;
	//��ȡ��������
	int k = 1;
	for(set<int>::iterator itSet=setNewCenter.begin(); itSet!=setNewCenter.end(); ++itSet,++k)
	{
		for(itStrPoint=listPoint.begin(); itStrPoint!=listPoint.end(); ++itStrPoint)
		{
			if(*itSet == itStrPoint->nID) //������
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

	//����
	if(!partitionClustering(listPoint,listK))
	{
		cout<<"����ʧ��"<<endl;
	}
	cout<<"����ĿΪ��"<<(k-1)<<endl;

	return true;
}



int main()
{
	//AP�㷨����
	list<StructPoint> listPoint;
	list<StructPoint> listK; //������
	ifstream inFile;
	inFile.open("./AP_example.txt");
	StructPoint strPoint;
	if(!inFile)
	{
		cout << "��ȡAP��������ʧ��" << endl;
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

	//����
	cout<<"���ھ���"<<endl;
	float fFactor = 0.5;  //����ϵ��Ϊ0.5 
	int nIt = 500;        //�������� 
	int nK = 10;		  //���ݵ��� 
	if(!(cluster_AP_int(listPoint, listK, nIt, nK, fFactor))) //AP����
	{
		cout<<"����ʧ��"<<endl;
		//�ͷ��ڴ�
		listPoint.clear();
		return 0;
	}

	cout<<"����ɹ�"<<endl;

	//���������
	ofstream outFile("./cluster_AP.txt"); 
	if(!outFile)
	{
		cout << "����ʧ��" << endl;
	}
	for(list<StructPoint>::iterator k=listPoint.begin(); k!=listPoint.end(); ++k)
	{
		outFile << k->nID << "\t" << k->nType << endl;
	}
	outFile.close();

	cout<<"������Ϊ��"<<endl;
	for(list<StructPoint>::iterator k=listK.begin(); k!=listK.end(); ++k)
	{
		cout << k->nID << "\t" << k->nType << "\t" << k->pointX << "\t" << k->pointY;
		cout<<endl;
	}
	//�ͷ��ڴ�
	
	listPoint.clear();
	listK.clear();

	return 0;
}
