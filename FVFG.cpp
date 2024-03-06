/* 
 * FVFG.cpp - Efficient Group based Hilbert Encoding and Decoding Algorithms 
 * 
 * Author:      Lianyin Jia
 *              Dept. of Computer Science
 *              Kunming University of Science & Technology
 * Date:        Mar 6 2024
 * Copyright (c)  Kunming University of Science & Technology
 *
 * Acknowledgement:
 * This implementation is based on the work of Lianyin Jia
 * "Efficient Group based Hilbert Encoding and Decoding Algorithms"
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
using namespace std;

#define ones(T,k) ((((T)2) << (k-1)) - 1) 
#define getBit(X,i,k) (X>>(i-k+1)& ones(halfmask_t,k))  
typedef unsigned long long bitmask_t;
typedef unsigned long halfmask_t;
 
//the four state views, here we use a 1D array to simulate the 3D array 
unsigned *CHM;
unsigned char *CSM;
unsigned *HCM;
unsigned char *HSM;

/*
	en_FVFG, the virtual-filling based group encoding algorithm.	
	X: the first coordinate component
	Y: the second coordinate component 
	order: the number of orders
	groupSize: group size
*/
bitmask_t en_FVFG(halfmask_t X, halfmask_t Y, int order, int groupSize)
{ 
	unsigned groupX = 0, groupY = 0; 
	bitmask_t hcode = 0;
	int pos = 0;
	int groupNum = (order-1)/groupSize +1; // the number of groups
	int startPos=groupNum*groupSize-1; // the encoding start position of each group
	int dimSize = 1<<groupSize;   // the dimension size in each group
	unsigned state = (groupSize*groupNum-order)%2 ? 1 :0;   //computing the state for the first group 
	for (short i = groupNum-1; i >= 0; i--)
	{
		// fetch group value of X and Y
		groupX = getBit(X,startPos,groupSize);
		groupY = getBit(Y,startPos,groupSize);
		startPos -= groupSize;	 
		pos = state * dimSize * dimSize + dimSize * groupX + groupY;		  
		hcode = (hcode << 2*groupSize) | CHM[pos];		 
		state = CSM[pos];
	} 	
	return hcode;
}

/*
	De_FVFG_batch, the batched en_FVFG algorithm.	
	X: the first coordinate component
	Y: the second coordinate component 
	order: the number of orders
	groupSize: group size
	groupNum: the number of groups 
	dimSize: the dimension size in each group 
	startPos: the encoding start position of each group
	state: Hilbert state of the first group 
*/
bitmask_t en_FVFG_batch(halfmask_t X, halfmask_t Y, int order, int groupSize,int groupNum,int dimSize,int startPos, unsigned char state)
{
	unsigned bitX = 0, bitY = 0;
	bitmask_t hcode = 0;
	int pos = 0;  
	for (short i = groupNum-1; i >= 0; i--)
	{
		bitX = getBit(X,startPos,groupSize);
		bitY = getBit(Y,startPos,groupSize);
		startPos-=groupSize;
		pos = state * dimSize * dimSize + dimSize * bitX + bitY; 
		hcode = (hcode << 2*groupSize) | CHM[pos];	 
		state = CSM[pos];
	} 
	return hcode;
}

/*
	de_FVFG, the virtual-filling based group decoding algorithm.	
	X: the first coordinate component
	Y: the second coordinate component 
	order: the number of orders
	groupSize: group size
*/
void de_FVFG(bitmask_t hcode,halfmask_t &X, halfmask_t &Y, int order,int groupSize)
{ 	 
	unsigned bitZ = 0;
	halfmask_t posKey=0;
	int groupNum = (order-1)/groupSize +1;
	int dimSize = 1<<groupSize;
	int startPos=2*groupNum*groupSize-1;
	halfmask_t mask = ones(halfmask_t,groupSize);
	unsigned char state = (groupSize*groupNum-order)%2 ? 1 : 0;	 //computing the state for the first group 
	for (short i = groupNum-1; i >= 0; i--)
	{ 
		bitZ = getBit(hcode,startPos,2*groupSize); 
		startPos-=2*groupSize;
		posKey = HCM[state*dimSize*dimSize + bitZ]; 
		Y = Y <<groupSize | posKey & mask;
		X = X <<groupSize | posKey>>groupSize & mask;
		state = HSM[state*dimSize*dimSize + bitZ];
	} 
}


/*
	de_FVFG_batch, the batched de_FVFG algorithm.	
	X: the first coordinate component
	Y: the second coordinate component 
	order: the number of orders
	groupSize: group size
	groupNum: the number of groups 
	dimSize: the dimension size in each group 
	startPos: the encoding start position of each group
	state: Hilbert state of the first group 
*/
void de_FVFG_batch(bitmask_t hcode,halfmask_t &X, halfmask_t &Y, int order,int groupSize,int groupNum,int dimSize,int startPos,unsigned char state)
{ 	
	unsigned bitZ = 0;
	halfmask_t posKey=0; 	
	halfmask_t mask = ones(halfmask_t,groupSize); 
	for (short i = groupNum-1; i >= 0; i--)
	{ 	
		bitZ = getBit(hcode,startPos,2*groupSize); 
		startPos-=2*groupSize;		
		posKey = HCM[state * dimSize * dimSize + bitZ]; 
		Y = Y <<groupSize | posKey & mask;
		X = X <<groupSize | posKey>>groupSize & mask;	
		state = HSM[state * dimSize * dimSize + bitZ];
	}  
}
/*
reading the state-views from files 
*/
template <typename T> 
bool readSVs(const string& filename, T *arr, int groupSize) {
	ifstream file(filename, ios::binary);
	if (!file.is_open()) {
		cerr << "Error opening file " << filename.c_str()<< endl;
		return false;
	} 
	int dimSize = 1<<groupSize;
	for (int i = 0; i < 4 * dimSize * dimSize; ++i) {
		T val; 
		if (!file.read(reinterpret_cast<char *>(&val), sizeof(T))) {
			cerr << "Error reading file " << filename.c_str() << endl;
			return false;
		} 
		arr[i] = val; 
	}

	file.close();
	return true;
} 

int main()
{
	int n = 21; //the number of orders of the input data. This value should be between 1 to 32.
	int g = 8;  // group size. Do not change this value unless you have generated the state-views for the specific g. 
	unsigned endVal = 10; // we simply test the first endVal*endVal input data here.This value can be changed, however, it should be less than 1<<n.
	int dimSize=1<<g; // the dimension size in each group
	bitmask_t hcode;  // Hilbert code
	halfmask_t x, y;  //coordinate
	int groupNum;  // the number of groups
	int reminder;  // the reminder of the first group
	unsigned char state; // the state of current order
	int startPos; // the encoding or decoding start position of each group

	//Reading state views from file, the four state views are CHM, CSM, HCM and HSM.
	cout<<"_________Reading the state views_____________"<<endl;
	CHM = new unsigned[4*dimSize*dimSize];
	CSM = new unsigned char[4*dimSize*dimSize];
	HCM = new unsigned[4*dimSize*dimSize];
	HSM = new unsigned char[4*dimSize*dimSize];
	readSVs <unsigned>("CHM_8.bin",CHM,g);
	readSVs <unsigned char>("CSM_8.bin",CSM,g);
	readSVs <unsigned>("HCM_8.bin",HCM,g);
	readSVs <unsigned char>("HSM_8.bin",HSM,g);
	 
	cout<<"_________Encoding_____________"<<endl;
	cout<<endl<<"Start test en_FVFG:"<<endl;
	for (x = 0; x < endVal; x++)
	{
		for(y=0;y< endVal;y++)
		{
			hcode= en_FVFG(x, y, n,g);  
			cout<<"The Hilbert code of coordinate ("<<x<<","<<y<<") is: "<<hcode<< endl;
		}
	}   

	cout<<endl<<"Start test en_FVFG_batch:"<<endl;	
	/*for batched algorithms, the groupNum, dimSize, state and startPos need to be computed only once for a batch,
	  so we compute them before the loops start.
	*/
	groupNum = (n-1) /g + 1; 
	reminder = n % g;
	state =reminder && (g-reminder)%2 ? 1:0; 
	dimSize=1<<g;	
	startPos=groupNum*g-1; 
	
	for (x = 0; x < endVal; x++)
	{
		for(y=0;y< endVal;y++)
		{
			hcode= en_FVFG_batch(x, y, n, g, groupNum, dimSize,startPos, state);  
			cout<<"The Hilbert code of coordinate ("<< x <<","<< y <<") is: "<<hcode<< endl;
		}
	}  
	
	cout<<endl<<"Starting test de_FVFG:"<<endl;
	for (hcode = 0; hcode <endVal*endVal; hcode++)
	{
		x =0;
		y =0;					
		de_FVFG(hcode, x, y, n, g);
		cout<<"The coordinate of Hilbert code "<< hcode <<" is: ("<< x <<","<< y <<")"<< endl; 
	} 
  
 
	cout<<endl<<"Starting test de_FVFG_batch:"<<endl;	
	/*for batched algorithms, the groupNum, dimSize, state and startPos need to be computed only once for a batch,
	  so we compute them before the loops start.
	*/
	groupNum = (n-1) / g +1; 
	reminder = n % g;
	state = reminder && (g-reminder)%2 ? 1 : 0;	 //computing the state for the first group
	dimSize=1<<g; 
	startPos=2*groupNum*g-1;
	for (hcode = 0; hcode < endVal*endVal; hcode++)
	{
		x =0;
		y =0;					
		de_FVFG_batch(hcode, x , y, n, g, groupNum, dimSize,startPos, state);
		cout<<"The coordinate of Hilbert code "<< hcode <<" is: ("<< x <<","<< y <<")"<< endl; 
	}  
	return 0;
}