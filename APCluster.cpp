//APCluster
#include<iostream>
#include<vector>
using namespace std;

// 
struct StructPoint{
	int id;
	int type;
	float pointX, pointY;
	bool operator == (const StructPoint &a) const{
		return (id == a.id);
	}
};

int main(){
	
	return 0;
}
