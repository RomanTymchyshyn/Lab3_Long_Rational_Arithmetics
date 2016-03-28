#pragma once

#include "iostream"
#include "fstream"
#include "cmath"
#include "cstdlib"
#include "vector"
#include "ctime"
#include "string"
#include "Support Functions.h"
#include <time.h>


using namespace std;

//vectors:

int max(const int &a, const int &b)
{
	if (a>b) return a;
	else return b;
}

void inverseVector(vector<int> &a)
{
	vector<int> supp=a;
	a.clear();
	for(int i=supp.size()-1; i>=0; --i)
	{
		a.push_back(supp[i]);
	}
	return;
}

void makeVectorCanonical(vector<int> &a, int size)
{
	vector<int>::iterator it;
	it=a.begin();
	a.insert(it,size-a.size(),0);
}

void outputVector(const vector<int> &a)
{
	for(int i=0; i<(int)a.size(); i++)
		cout<<a[i];
	cout<<endl;
}

void outputVectorInStream(ofstream &out, const vector<int> &a)
{
	for(int i=0; i<a.size(); i++)
		out<<a[i];
	out<<endl;
}

bool vectorNotLessThan(const vector<int> &a, const vector<int> &b)
{
	if(a.size()<b.size())
		return false;
	if(a.size()>b.size())
		return true;
	int i=0;
	while(i<a.size() && a[i]==b[i])
		i++;
	if(i>=a.size()) return true;
	if(a[i]>b[i]) return true;
	return false;
}

vector<int> sumVector(vector<int> a, vector<int> b, const int &base)
{
	if(a.size()>b.size())
		makeVectorCanonical(b,a.size());
	else
		makeVectorCanonical(a,b.size());
	vector<int> supp;
	int rest=0, sum=0;
	for(int i=a.size()-1; i>=0; i--)
	{
		sum=a[i]+b[i]+rest;
		rest=sum/base; 
		sum%=base;
		supp.push_back(sum);
	}
	if(rest==1) supp.push_back(1);
	inverseVector(supp);
	return supp;
}

vector<int> mutiplyVectorAndInt(const vector<int> &a, const int &n, const int &base)
{
	vector<int> supp;
	int sum=0, rest=0;
	for(int i=a.size()-1; i>=0; i--)
	{
		sum=a[i]*n+rest;
		rest=sum/base; 
		sum%=base;
		supp.push_back(sum);
	}
	while(rest>0)
	{
		supp.push_back(rest%base);
		rest/=base;
	}
	inverseVector(supp);
	return supp;
}

vector<int> mutiplyVectorAndVector(const vector<int> &a, const vector<int> &b, const int &base)
{
	vector<int> supp,sum;
	vector<int>::iterator it;
	for(int i=b.size()-1; i>=0; i--)
	{
		supp=mutiplyVectorAndInt(a,b[i],base);
		it=supp.end();
		supp.insert(it, b.size()-i-1, 0);
		sum=sumVector(sum,supp,base);
	}
	if(sum[0]==0)
	{
		int u=0;
		while(u<sum.size()-1 && sum[u]==0)
			u++;
		sum.erase(sum.begin(), sum.begin() +u);
	}
	return sum;
}

vector<int> indianPowForVector(const vector<int> &a, int n, const int &base)
{
	vector<int> supp=a;
	vector<int> res(1,1);
	int maxBinPow=31;
	while(!(n>>maxBinPow & 0x1))
		maxBinPow--;
	for(int i=0; i<=maxBinPow; i++)
	{
		if(n>>i & 0x1) 
			res=mutiplyVectorAndVector(res,supp,base);
		supp=mutiplyVectorAndVector(supp,supp,base);
	}
	return res;
}

vector<int> IpowVector(const vector<int> &a, vector<int> n, const vector<int> &roll, const int &base )
{
	vector<int> supp(1,1), _n=n, pows, bin(1,2);
	int maxBinPow=(int)((double)(n.size())/(double)log10((double)2))+1;
	vector<vector<int>> binPow; //=indianPowForVector(bin,maxBinPow,base);
	for(int i=0; i<maxBinPow+1; i++)
	{
		binPow.push_back(supp);
		supp=mutiplyVectorAndInt(supp, 2, base);
	}
	vector<int> res(1,1);
	int j=maxBinPow;
	bool check=false;
	while(j>-1)
	{
		if(vectorNotLessThan(_n,binPow[j])) 
		{
			pows.push_back(1);
			_n=subtractVectors(_n, binPow[j], base);
			check=true;
		}
		else
			if(check==true)
				pows.push_back(0);
		j--;
	}
	supp.clear();
	supp=a;
//	if(pows[pows.size()-1]==1)
//		res=moduloVector(mutiplyVectorAndVector(res, a, base), roll, base);
	for(int i=pows.size()-1; i>-1; i--)
	{
		if(pows[i]==1)
			res=moduloVector(mutiplyVectorAndVector(res, supp, base), roll, base);
		supp=moduloVector(mutiplyVectorAndVector(supp,supp,base), roll, base);
	}

	//cout<<"mult= "<<tM<<endl;
	//cout<<"div= "<<tD<<endl;
	return res;
}

vector<int> subtractVectors(vector<int> a,vector<int> b, const int &base)
{
	makeVectorCanonical(b,a.size());
	vector<int> supp;
	int rest=0, sum=0;
	for(int i=a.size()-1; i>=0; i--)
	{
		if(a[i]-rest-b[i]<0)
		{
			sum=a[i]-rest-b[i]+base;
			rest=1;
		}
		else
		{
			sum=a[i]-rest-b[i];
			rest=0;
		}
		supp.push_back(sum);
	}
	inverseVector(supp);
	int i=2;
	while(i>1 && supp[0]==0)
	{
		supp.erase(supp.begin());
		i=supp.size();
	}
	return supp;
}

vector<int> divideVectorByVector(vector<int> a, vector<int> b, const int &base)
{
	vector<int> res;
	if(!vectorNotLessThan(a,b))
	{
		res.push_back(0);
		return res;
	}
	if (b[0]<base/2)
	{
		int d=base/(b[0]);
		a=mutiplyVectorAndInt(a,d,base);
		b=mutiplyVectorAndInt(b,d,base);
	}
	vector<int> _a=a;
	vector<int> _b;
	int k;
	vector<int> supp=b;
	int shift=a.size()-b.size();
	supp.insert(supp.end(),shift,0);
	if(!vectorNotLessThan(a,supp))
		shift--;	
	while(vectorNotLessThan(_a,b) || shift>-1)
	{
		if(_a.size()==0)
		{
			res.insert(res.end(),shift+1,0);
			return res;
		}
		supp=b;
		supp.insert(supp.end(),shift,0);
		if(_a.size()>supp.size() && _a.size()>1)
			k=(base*_a[0]+_a[1])/supp[0]+1;
		else 
			k=_a[0]/supp[0]+1;
		_b=mutiplyVectorAndInt(supp,k,base);
		while(!vectorNotLessThan(_a,_b))
		{
			_b=subtractVectors(_b,supp,base);
			k--;
		}
		_a=subtractVectors(_a,_b,base);
		res.push_back(k);
		shift--;
	}
	return res;
}

int convertVectorIntoInt(vector<int> a, const int &base)
{
	int res=0;
	int supp=1;
	int i=a.size();
	while(i>0)
	{
		i--;
		res+=a.at(i)*supp;
		supp=supp*base;
	}
	return res;
}

vector<int> convertIntIntoVector(int a, const int &base)
{
	vector<int> res;
	int _a=a;
	int supp=1;
	int i=0;
	while(_a>=base)
	{
		_a=_a/base;
		supp*=base;
		i++;
	}
	int j=0;
	while(j<=i)
	{
		res.push_back((a-a%supp)/supp);
		a=a%supp;
		supp/=base;
		j++;
	}
	return res;
}

vector<int> moduloVector(vector<int> a, vector<int> b, const int &base)
{
	vector<int> temp=divideVectorByVector(a,b,base);
	vector<int> supp=subtractVectors(a,mutiplyVectorAndVector(temp,b,base),base);
	return supp;
}

vector<int> convertInBase(vector<int> a, const int &base)
{

	/*
	vector<int> supp(1,1), _n=a, pows, bin(1,2);
	int maxBasePow=(int)((double)(a.size())/(double)log10((double)base))+2;
	vector<vector<int>> basePow; //=indianPowForVector(bin,maxBinPow,base);
	for(int i=0; i<maxBasePow+1; i++)
	{
		basePow.push_back(supp);
		supp=mutiplyVectorAndInt(supp, base, 10);
	}
	vector<int> res(1,1);
	int j=maxBasePow;
	bool check=false;
	while(j>-1)
	{
		pows.push_back(0);
		while(vectorNotLessThan(_n,basePow[j])) 
		{
			pows[maxBasePow-j]++;
			_n=subtractVectors(_n, basePow[j], 10);
			//check=true;
		}
		j--;
	}
	return pows;

*/

	
	vector<int> res;
	vector<int> _base=convertIntIntoVector(base,10);
	vector<int> supp;
	vector<int> temp=a;
	vector<int> nullVec(1,0);
	while(a!=nullVec)
	{
		temp=divideVectorByVector(a,_base,10);
		supp=subtractVectors(a,mutiplyVectorAndInt(temp,base,10),10);
		//supp=moduloVector(a,_base,10);
		res.push_back(convertVectorIntoInt(supp,10));
		a=temp;
	}
	inverseVector(res);
	return res;
}

vector<int> convertFromBase(vector<int> a, const int &base)
{
	vector<int> res;
	vector<int> _base=convertIntIntoVector(base,10);
	vector<int> supp(1,1);
	int i=a.size();
	while(i>0)
	{
		i--;
		res=sumVector(res, mutiplyVectorAndInt(supp, a[i], 10),10);
		supp=mutiplyVectorAndVector(supp,_base,10);
	}
	return res;
}

void getVectorsFromStream(ifstream &in, vector<int> &a, vector<int> &b)
{
	if(!in) 
	{
		cout<<"input stream error"<<endl;
		return;
	}
	a.clear();
	b.clear();
	string _a,_b;
	in>>_a;
	in>>_b;
	for(int i=0; i<_a.size(); i++)
		a.push_back(_a[i]-48);
	for(int i=0; i<_b.size(); i++)
		b.push_back(_b[i]-48);
	return;	
}

bool soloveiShtrassen(vector<int> m, int k, int base)
{
	if(base!=10)
	{
		base=10;
		convertFromBase(m, base);
	}
	srand(time(NULL));
	for(int j=0; j<k; j++)
	{
	vector<int> w,temp,one(1,1),bin(1,2);
	int size=rand()%(m.size()-2)+1;
	if (size<m.size()/10)
		size+=2*m.size()/10;
	for(int i=0; i<size; i++)
		w.push_back(rand()%9+1);
	if(gcd(m,w,base)!=one)
		return false;
	temp=divideVectorByVector(subtractVectors(m,one,base),bin,base);
	vector<int> q1=IpowVector(w,temp,m,base);
	vector<int> q2=yakobi(w,m,base);
	if(q2!=one && q2!=subtractVectors(m,one,base))
		return false;
	if(q1!=q2)
		return false;
	}
	return true;
}

vector<int> yakobi(vector<int> a, vector< int> P, const int &base)
{
	vector<int> bin(1,2), one(1,1);
	vector<int> res=IpowVector(a,divideVectorByVector(subtractVectors(P,one,base), bin, base), P, base);
	res=moduloVector(res,P,base);
	return res;
}

vector<int> gcd(vector<int> a, vector<int> b, const int &base)
{
	vector<int> nullVec(1,0),temp;
	while(b!=nullVec && b.size()>0)
	{
		a=moduloVector(a,b,base);
		temp=a;
		a=b;
		b=temp;
	}
	return a;
}