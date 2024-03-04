#include <iostream>
using namespace std;

#define ones(T,k) ((((T)2) << (k-1)) - 1) 
#define getBit(X,i,k) (X>>(i-k+1)& ones(halfmask_t,k)) 
#define ODD(x) ((x)%2==1?1:0)
typedef unsigned long long bitmask_t;
typedef unsigned long halfmask_t;

unsigned CHM[4][256][256];
unsigned char CSM[4][256][256];
unsigned HCM[4][256*256];
unsigned char HSM[4][256*256];

bitmask_t en_FVFG(halfmask_t GridX, halfmask_t GridY, int order, int groupSize, unsigned *arKey, unsigned char *arType)
{
	unsigned nType = 0;
	unsigned bitX = 0, bitY = 0;
	//resKey变量用来存放最终计算结果
	bitmask_t resKey = 0;
	int pos = 0;
	int groupNum = (order-1)/groupSize +1; // 组数
	int startPos=groupNum*groupSize-1;
	int dimSize = 1<<groupSize;
	//int numberOfFirstGroup = k % n;
	if(ODD(groupSize*groupNum-order)) nType =1;  
	for (short i = groupNum-1; i >= 0; i--)
	{
		//取位
		bitX = getBit(GridX,startPos,groupSize);
		bitY = getBit(GridY,startPos,groupSize);
		startPos-=groupSize;
		pos = nType*dimSize*dimSize + dimSize*bitX+bitY;
		//获取当前组
		//resKey = (resKey << 2*groupSize) | arKey[nType][bitX][bitY];
		resKey = (resKey << 2*groupSize) | arKey[pos];
		//查矩阵取type值，用于循环的下一次计算
		nType = arType[pos];
	} 
	//计算完毕，返回结果
	return resKey;
}

bitmask_t en_FVFG_batch(halfmask_t GridX, halfmask_t GridY, int order, int groupNum, int groupSize,int dimSize, unsigned char nType, unsigned *arKey, unsigned char *arType)
{
	//unsigned nType = 0;
	unsigned bitX = 0, bitY = 0;

	//resKey变量用来存放最终计算结果
	bitmask_t resKey = 0;
	int pos = 0; 
	//若startPos指向位置超出32位整型范围,为避免后面用if判断，可把第一组也传过来
	int startPos=groupNum*groupSize-1; 
	//int numberOfFirstGroup = k % n; 
	for (short i = groupNum-1; i >= 0; i--)
	{
		//取位
		bitX = getBit(GridX,startPos,groupSize);
		bitY = getBit(GridY,startPos,groupSize);
		startPos-=groupSize;
		pos = nType*dimSize*dimSize + dimSize*bitX+bitY;
		//获取当前组
		//resKey = (resKey << 2*groupSize) | arKey[nType][bitX][bitY];
		resKey = (resKey << 2*groupSize) | arKey[pos];
		//查矩阵取type值，用于循环的下一次计算
		nType = arType[pos];
	} 
	//计算完毕，返回结果
	return resKey;
}


void de_FVFG(bitmask_t index,halfmask_t &GridX, halfmask_t &GridY, int order,int groupSize,unsigned *invKey, unsigned char *invType)
{ 	 
	unsigned nType = 0;
	unsigned bitZ = 0;
	halfmask_t posKey=0;
	int groupNum = (order-1)/groupSize +1;
	int dimSize = 1<<groupSize;
	int startPos=2*groupNum*groupSize-1;
	halfmask_t mask = ones(halfmask_t,groupSize);
	if(ODD(groupSize*groupNum-order)) nType =1;  
	for (short i = groupNum-1; i >= 0; i--)
	{ 
		//取位
		bitZ = getBit(index,startPos,2*groupSize); 
		startPos-=2*groupSize;
		//查矩阵取值并与上一层结果合并(上一层结果左移两位与本次查矩阵结果求或运算)
		posKey = invKey[nType*dimSize*dimSize + bitZ]; 
		GridY = GridY <<groupSize | posKey & mask;
		GridX = GridX <<groupSize | posKey>>groupSize & mask;
		//查矩阵取type值，用于循环的下一次计算
		nType = invType[nType*dimSize*dimSize + bitZ];
	} 
}

void de_FVFG_batch(bitmask_t index,halfmask_t &GridX, halfmask_t &GridY, int order,int groupNum,int groupSize,unsigned char nType,unsigned *invKey, unsigned char *invType)
{ 	
	unsigned bitZ = 0;
	halfmask_t posKey=0; 
	int dimSize = 1<<groupSize;
	int startPos=2*groupNum*groupSize-1;
	halfmask_t mask = ones(halfmask_t,groupSize); 
	for (short i = groupNum-1; i >= 0; i--)
	{ 
		//取位
		bitZ = getBit(index,startPos,2*groupSize); 
		startPos-=2*groupSize;
		//查矩阵取值并与上一层结果合并(上一层结果左移两位与本次查矩阵结果求或运算)
		posKey = invKey[nType*dimSize*dimSize + bitZ]; 
		GridY = GridY <<groupSize | posKey & mask;
		GridX = GridX <<groupSize | posKey>>groupSize & mask;
		//查矩阵取type值，用于循环的下一次计算
		nType = invType[nType*dimSize*dimSize + bitZ];
	}  
}

 

int main()
{
	int n = 32; //total number of orders
	int g = 8;  // group size
	unsigned endVal = 1000; // endVal should be 
	unsigned *p_key;
	unsigned char *p_type;
	bitmask_t res;  
	p_key = **CHM;
	p_type = **CSM;
	for (halfmask_t i = 0; i <= endVal; i++)
	{
		for(halfmask_t j=0;j<=endVal;j++)
		{
			res= en_FVFG(i, j, n,g,p_key,p_type);  
			cout<<"The result of coordinate ("<<i<<","<<j<<") is"<<res<< endl;
		}
	}  
	return 0;
}