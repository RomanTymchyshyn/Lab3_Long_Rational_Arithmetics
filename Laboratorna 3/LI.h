#ifndef _CLASS_LI_
#define _CLASS_LI_


#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <list>
#include <complex>
#include <cmath>
#include <stdlib.h>
#include <time.h>


using namespace std;

enum MultType{Karatsuba,ToomCook,Schonhage,SchonhageShtrassen};


class LI
{
		bool is_fractional;//קט ÷ הנמבמגטל
	public:
		vector<int> long_num;
		static const int BASE=32768;
		LI ();//+
		LI(string a);//+
		LI(int a);//+
		bool isFrac()const {return is_fractional;}//+
		LI operator=(LI a);//+
		LI operator=(int a);//+
		LI operator+(LI a);//+
		LI operator+(int a);//+
		LI operator-(LI & a);//+
		LI operator-(int a);//+
		virtual LI operator*(LI a){return Karatsuba(a);}//+
		LI operator*(int a);//+
		LI operator/(LI a);//+
		LI operator/(int a);//+
		bool operator!=(LI a);//+
		bool operator==(LI a);//+
		bool operator>(LI a);//+
		bool operator<(LI a);//+
		bool operator>=(LI a);//+
		bool operator<=(LI a);//+
		int operator%(int a);//+
		double toDouble();
		double DoubleDivide(LI v);
		LI operator%(LI a);//+
		LI operator^(const LI &n);//+

		vector<int> toDec()const;						//+
		vector<int> toBin()const;						//+
		vector<int> toBin1()const;//רגטהרטי חא toBin	//+
		void inBin(vector<int> b);						//+
		void inBin1(vector<int> b);//רגטהרטי חא inBin	//+
		int last()const {return (long_num.size()>0? long_num.at(long_num.size()-1) : -1);}//+
		void pop_back() { long_num.pop_back(); }//+
		
		bool is_frac()const{return is_fractional;}//+
		void convert_to_LI(string a);//+
		string convert_to_string() const;//+
		void LI_delete_zeros();//+

		LI Primitive (LI a);//+
		LI Karatsuba (LI a);//+
		LI ToomCook  (LI a);//+
		LI Schonhage(LI a);//+
		LI CookDivide(LI v);//+
		LI mod(const int &n);//+
		bool Lemer();//+
		bool isMersenne(int &p);//+
		LI SchonhageShtrassen(LI a);//+
		LI FastPowMod(const LI &n, const LI &m);//+
		bool MillerRabin(const int &r);//+
		bool SolovejShtrassen(int k);
};

template <MultType a>
class LongInt:public LI
{
public:
	virtual LI operator*(LI a);
};

template <MultType b>
LI LongInt<b>::operator*(LI a)
{
	switch(b)
	{
	case 0: return Karatsuba(a);
	case 1: return ToomCook(a);
	case 2: return Schonhage(a);
	case 3: return SchonhageShtrassen(a);
	}
}

void str_delete_zeros(string &a);

vector<complex<double>> FFT(vector<complex<double>>);

vector<complex<double>> FFT_REV(vector<complex<double>>);

long long round(double x);

bool operator>>(istream &i, LI & a);

bool operator<<(ostream &o, LI a);

string FastPow(string x,int n);

string div(string a,int b, string & rest);//a/b

#define ץ Gauss(right)

#endif