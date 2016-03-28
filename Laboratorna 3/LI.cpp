#include "LI.h"
#include "Support Functions.h"

const int MAX_BUF = 700;

const double pi=atan(1.0f)*4;

bool operator>>(istream &i, LI &a)
{
	string s;
	if(!(i>>s))return false;
	a.convert_to_LI(s);
	return true;
}

bool operator<<(ostream & o, LI a)
{
	string s;
	s=a.convert_to_string();
	if(a.isFrac())
	{
		if(!(o<<"0."))return false;
		while(a.last()==0)
		{
			a.pop_back();
			if(!(o<<"0"))return false;
		}
	}
	if(!(o<<s))return false;
	return true;
}

LI::LI()
{
	is_fractional=false;
	long_num.clear();
}

LI::LI(string a)
{
	is_fractional=false;
	if(a!="")
		convert_to_LI(a);
	else 
		long_num.clear();
}

LI::LI(int a)
{
	is_fractional=false;
	if (a!=0)
	{
		while(a>BASE)
		{
			long_num.push_back(a%BASE);
			a/=BASE;
		}
		long_num.push_back(a);
	}
	else long_num.clear();
}

bool LI::operator==(LI a)
{
	LI_delete_zeros();
	a.LI_delete_zeros();
	if (long_num.size()!=a.long_num.size()) return false;
	while(a.long_num.size() < long_num.size()) a.long_num.push_back(0);
	while(a.long_num.size() > long_num.size()) long_num.push_back(0);
	for(int i=0;i<(int)long_num.size();++i)
		if((long_num[i])!=(a.long_num[i])) return false;
	return true;
}

LI LI::operator=(LI a)
{
	if ((*this)==a) return (*this);
	is_fractional=a.is_fractional;
	if (a.long_num.size()>0)
		long_num=a.long_num;
	else long_num.clear();
	return (*this);
}

LI LI::operator=(int a)
{
	LI b(a);
	(*this)=b;
	return *this;
}

LI LI::operator+(LI a)
{
	int shift = 0;
	int temp_result;
	LI result;
	result.is_fractional = is_fractional;
	while (a.long_num.size() < long_num.size()) a.long_num.push_back(0);
	while (a.long_num.size() > long_num.size()) long_num.push_back(0);
	for(int i=0; i<(int)long_num.size(); i++)
	{
		temp_result = long_num[i]+a.long_num[i]+shift;
		result.long_num.push_back(temp_result%BASE);
		shift = temp_result/BASE;
	}
	a.LI_delete_zeros();
	if(shift!=0) result.long_num.push_back(shift);
	LI_delete_zeros();
	return result;
}

LI LI::operator+(int a)
{
	LI temp(a);
	return temp+(*this);
}

LI LI::operator-(LI & a)
{
	a.LI_delete_zeros();
	if((*this)<a) return LI();
	int shift=0;
	int temp_result;
	LI result;
	while (a.long_num.size() < long_num.size()) a.long_num.push_back(0);
	while (a.long_num.size() > long_num.size()) long_num.push_back(0);

	for(int i=0; i<(int)long_num.size(); i++)
	{
		temp_result=long_num[i]-a.long_num[i]-shift;
		if (temp_result<0)
		{
			temp_result+=BASE;
			shift=1;
		}
		else shift=0;
		result.long_num.push_back(temp_result);
	}
	a.LI_delete_zeros();
	LI_delete_zeros();
	return result;
}

LI LI::operator-(int a)
{
	LI temp(a);
	if((*this)<a) return LI ();
	return (*this)-temp;
}

LI LI::operator*(int a)
{
	LI result;
	if (a==0 || (*this)==result) return result;
	int shift=0;
	int temp_result;
	for(int i=0; i<(int)long_num.size(); ++i)
	{
		temp_result=a*long_num[i]+shift;	
		shift=temp_result/BASE;
		result.long_num.push_back(temp_result%BASE);
	}
	if (shift!=0)
	{
		while(shift>=BASE)
		{
			result.long_num.push_back(shift%BASE);
			shift/=BASE;
		}
		result.long_num.push_back(shift);
	}
	return result;
}

void reverse(vector <int> & b)
{
	int temp=0;
	for (unsigned int i=0; i<b.size()/2; ++i)
	{
		temp=b[i];
		b[i]=b[b.size()-1-i];
		b[b.size()-1-i]=temp;
	}
}

LI LI::operator/(int a)
{
	if(a==0) return LI();
	int temp;
	int shift=0;
	LI result;
	for(int  i=long_num.size()-1; i>=0; --i)
	{
		temp=long_num[i]+shift*BASE;
		result.long_num.push_back(temp/a);
		shift=temp%a;
	}
	reverse(result.long_num);
	return result;
}

bool LI::operator!=(LI a)
{
	LI_delete_zeros();
	a.LI_delete_zeros();
	if (long_num.size()!=a.long_num.size()) return true;
	while (a.long_num.size() < long_num.size()) a.long_num.push_back(0);
	while (a.long_num.size() > long_num.size()) long_num.push_back(0);
	if (long_num.size()==0) return false;
	for (int i=0;i<(int)long_num.size();++i)
		if ((long_num[i])!=(a.long_num[i])) return true;
	return false;
}

bool LI::operator<(LI a)
{
	LI_delete_zeros();
	a.LI_delete_zeros();
	if (long_num.size() < a.long_num.size()) return true;
	else if (long_num.size() > a.long_num.size())return false;
	int i = long_num.size() - 1;
	while (i>=0 && long_num[i] == a.long_num[i])--i;
	if (i==-1)return false;
	return long_num[i] < a.long_num[i];
}

bool LI::operator>(LI a)
{
	LI_delete_zeros();
	a.LI_delete_zeros();
	if (long_num.size() > a.long_num.size()) return true;
	else if (long_num.size() < a.long_num.size()) return false;
	int i = long_num.size() - 1;
	while (i>=0 && long_num[i] == a.long_num[i])--i;
	if (i==-1)return false;
	return long_num[i] > a.long_num[i];
}

bool LI::operator<=(LI a)
{
	LI_delete_zeros();
	a.LI_delete_zeros();
	if (long_num.size() < a.long_num.size()) return true;
	else if (long_num.size() > a.long_num.size())return false;
	int i = long_num.size() - 1;
	while (i>=0 && long_num[i] == a.long_num[i])--i;
	if (i==-1) return true;
	return long_num[i] < a.long_num[i];
}

bool LI::operator>=(LI a)
{
	LI_delete_zeros();
	a.LI_delete_zeros();
	if (long_num.size() > a.long_num.size()) return true;
	else if (long_num.size() < a.long_num.size())return false;
	int i = long_num.size() - 1;
	while (i>=0 && long_num[i] == a.long_num[i])--i;
	if (i==-1)return false;
	return long_num[i] < a.long_num[i];
}

void LI::convert_to_LI(string a)
{
	string X=a;
	string rest;
	X=div(a,BASE,rest);
	while(atoi(X.c_str())!=0)
	{
		long_num.push_back(atoi(rest.c_str()));	
		X=div(X,BASE,rest);
	}
	long_num.push_back(atoi(rest.c_str()));
}

string operator+(string b,string a)
{
	string res;
	char buf[MAX_BUF];
	vector<int> dec;
	int temp=0;
	int k=(int)a.size()-(int)b.size();
	if(k>0)//a>b
		for(int i=0;i<k;++i)
			b.insert(0,"0");
	else
	{
		k=-k;
		for(int i=0;i<k;++i)
			a.insert(0,"0");
	}
	for(int i=a.size()-1,j=b.size()-1; i>=0 && j>=0; --i,--j)
	{
		temp+=(int)a[i]-48+(int)b[j]-48;
		if(temp<10) {dec.push_back(temp);temp=0;}
		else {temp=temp-10;dec.push_back(temp);temp=1;}
	}
	if(temp!=0)dec.push_back(temp);
	for(int i=(int)dec.size()-1;i>=0;--i){
		_itoa_s(dec[i],buf,10);
		res.append(buf);
	}
	return res;
}

string operator*(string a, string b)//множенн€ у стовпчик
{
	string result="0";
	string temp="";
	char buf[MAX_BUF];
	if(a.size()<b.size()) a.swap(b);
	//a>b
	int t;
	int shift=0;
	bool c=true;
	for(int i=(int)b.size()-1; i>=0; --i)
	{
		shift=0;
		for(int j=(int)a.size()-1; j>=0; --j)
		{
			t=b[i]-48;
			t=t*(a[j]-48) + shift;
			shift=t/10;
			_itoa_s(t%10,buf,10);
			temp.insert(0,buf);
		}
		_itoa_s(shift,buf,10);
		temp.insert(0,buf);
		for(int t=i+1; t<(int)b.size(); ++t)
			temp.append("0");//реал≥зуЇтьс€ здвиг
		result=result+temp;
		temp="";
	}
	str_delete_zeros(result);
	return result;
}

string div(string a, int b, string & rest)//a/b
{
	int base_iter=1;
	for(int x=b;x/10!=0;x=x/10)
		base_iter++;
	string result;
	char buf[MAX_BUF];
	memset(buf,LI::BASE,0);
	int iter=base_iter;
	int size=(int)a.size();
	int i=0;
	while(1)
	{
		iter=base_iter;
		int X=atoi(a.substr(i,iter).c_str());//переводимо п≥др€док у char*
		if((int)a.size() >= i+iter)
			{if(X<b) X=atoi(a.substr(i,(++iter)).c_str());}
		else 
			{
				rest=a;
				break;
			}
		if(X<b && (int)a.size() > i+iter)
		{
			result.append("0");
			X=atoi(a.substr(i,(++iter)).c_str());
		}
		_itoa_s(X/b,buf,10);
		if(strcmp(buf,"")==0) strcpy_s(buf, "0");
		result.append(buf);
		_itoa_s(X%b,buf,10);
		rest=buf;
		if(i+iter >= (int)a.size()) break;
		if(rest=="")rest="0";
		while((int)rest.size() < base_iter)
			rest=rest.insert(0,"0");
		a.insert(iter+i,rest);//вставл€эмо решту п≥сл€ символ≥в, €к≥ ми вже використовували
		i+=iter;
	}
	return result;
}

string LI::convert_to_string()const
{
	if(long_num.size()==0) return "";
	string result;
	char buf[MAX_BUF];
	_itoa_s(long_num[0],buf,10);
	result.append(buf);
	string temp;
	string base;	
	_itoa_s(BASE,buf,10);
	base=buf;
	string temp_base=base;
	for(int i=1; i<(int)long_num.size(); ++i)
	{
		_itoa_s(long_num[i],buf,10);
		temp=string(buf)*temp_base;
		result=result+temp;
		temp_base=temp_base*base;
	}
	str_delete_zeros(result);
	return result;
}

void LI::LI_delete_zeros()
{
	if (long_num.size() == 0) return;
	while (long_num.size() != 0 && long_num[long_num.size()-1]==0) long_num.pop_back();
}

void str_delete_zeros(string &a)
{
	int i;
	for(i=0; i<(int)a.length() && a[i]=='0'; ++i);
	a.erase(0,i);
}

int LI::operator%(int a)
{
	if (a==0) return -1;
	int temp;
	int shift=0;
	for(int  i=long_num.size()-1; i>=0; --i)
	{
		temp=long_num[i]+shift*BASE;
		shift=temp%a;
	}
	return shift;
}

vector<int> operator+(vector<int> a,vector<int> b)
{
	int shift=0;
	int temp_res;
	vector<int> res;
	while(a.size()<b.size())a.push_back(0);
	while(a.size()>b.size())b.push_back(0);
	for(int i=0;i<(int)a.size();i++)
	{
		temp_res=a[i]+b[i]+shift;
		res.push_back(temp_res%10);
		shift=temp_res/10;
	}
	if(shift!=0)res.push_back(shift);
	return res;
}

vector<int> operator*(vector<int> a,int b)
{
	vector<int> res;
	if(b==0||a==res)return res;
	int shift=0;
	int temp_res;
	for(int i=0; i<(int)a.size(); ++i)
	{
		temp_res=b*a[i]+shift;	
		shift=temp_res/10;
		res.push_back(temp_res%10);
	}
	if(shift!=0)res.push_back(shift);
	return res;
}

vector<int> operator*(vector<int> a,vector<int> b)
{
	int shift=0;
	vector<int> res;
	//==0
	if(a==res||b==res)return res;
	vector<int> temp;
	if(a.size()>b.size())a.swap(b);
	for(int i=0;i<(int)a.size();++i)
	{
		temp=b*a[i];
		for(int j=0;j<i;++j)
			temp.insert(temp.begin(),0);
		res=res+temp;
	}
	return res;
}

vector<int> LI::toDec()const
{																				
	if(long_num.size()==0) return vector<int>();
	vector<int> result;
	int temp;
	temp=long_num[0];
	while(temp/10!=0)
	{
		result.push_back(temp%10);
		temp/=10;
	}
	result.push_back(temp%10);
	vector<int> base;
	temp=BASE;
	while(temp/10!=0)
	{
		base.push_back(temp%10);
		temp/=10;
	}
	base.push_back(temp%10);
	vector<int> temp_result;
	vector<int>temp_base=base;
	for(int i=1; i<(int)long_num.size(); ++i)
	{
		temp=long_num[i];
		while(temp/10!=0)
		{
			temp_result.push_back(temp%10);
			temp/=10;
		}
		temp_result.push_back(temp%10);
		temp_result=temp_result*temp_base;																												
		result=result+temp_result;
		temp_result.clear();
		temp_base=temp_base*base;
	}
	for(int i=0; i<(int)(result.size()-1)/2; ++i)
	{
		temp=result[i];
		result[i]=result[result.size()-1-i];
		result[result.size()-1-i]=temp;
	}
	return result;
}

vector<int> LI::toBin()const
{
	LI X=(*this);
	vector<int> result;
	while(X/2!=0)
	{
		result.push_back(X%2);
		X=X/2;
	}
	result.push_back(X%2);
	return result;
}

vector<int> LI::toBin1()const
{
	LI X=(*this);
	vector<int> result;
	for(int i=0,k=0; i<X.long_num.size(); ++i)
	{
		k=1;
		while(X.long_num[i]/2!=0)
		{
			result.push_back(X.long_num[i]%2);
			X.long_num[i]=X.long_num[i]/2;
			++k;
		}
		result.push_back(X.long_num[i]%2);
		while(k<15){result.push_back(0);++k;}
	}
	while(result[result.size()-1]==0) result.pop_back();
	return result;

}

double LI::toDouble()
{
	if(long_num.size()==0) return 0.;
	double res;
	int temp;
	temp=long_num[0];
	res=temp;
	double base=BASE;
	for(int i=1; i<(int)long_num.size(); ++i)
	{
		res+=long_num[i]*base;
		base=BASE*base;
	}
	return res;
}

double LI::DoubleDivide(LI v)
{
	double x;
	LI u=(*this);
	LI l=u/v;
	x=l.toDouble();
	//if (x>UINT_MAX) return x;
	double x0=0;
	u=u%v;
	if(u==0)return x;
	double y=10;
	double xi;
	for(int i=0; i<16 && u!=0; ++i)
	{
		u=u*10;
		l=u/v;
		xi=l.toDouble();
		x0+=xi/y;
		y*=10;
		u=u%v;
	}
	x=x+x0;
	return x;
}

void LI::inBin(vector<int> b)
{
	long_num.clear();
	LI base=1;
	for(int i=0; i<(int)b.size(); ++i)
	{
		(*this)=(*this)+base*b[i];
		base=base*2;
	}
}

void LI::inBin1(vector<int> b)
{
	long_num.clear();
	int i=0;
	int temp=0;
	int base=1;
	for(i=0; i<b.size(); i+=15)
	{
		base=1;
		temp=0;
		for(int j=0; j<15 && i+j<b.size(); ++j)
		{
			temp=temp+base*b[i+j];
			base*=2;
		}
		long_num.push_back(temp);
	}
}

LI LI::operator^(const LI &n)
{
	vector<int> m;
	//n in 2
	m=n.toBin1();
	LI s=(*this);
	for(int i=(int)m.size()-2; i>=0; --i)
		if(m[i]==0)s=s*s;
		else s=s*s*(*this);
	return s;
}

LI LI::operator/(LI a)
{
	LI result;
	vector<int> c,d;
	LI_delete_zeros();
	a.LI_delete_zeros();
	for(int i=long_num.size()-1; i>=0; --i)
		c.push_back(long_num[i]);
	for(int i=a.long_num.size()-1; i>=0; --i)
		d.push_back(a.long_num[i]);
	vector<int> res=divideVectorByVector(c,d,BASE);
	for(int i=res.size()-1; i>=0; --i)
		result.long_num.push_back(res[i]);
	return result;
}

LI LI::operator%(LI a)
{
	LI result;
	vector<int> c,d;
	LI_delete_zeros();
	a.LI_delete_zeros();
	for(int i=long_num.size()-1; i>=0; --i)
		c.push_back(long_num[i]);
	for(int i=a.long_num.size()-1; i>=0; --i)
		d.push_back(a.long_num[i]);
	vector<int> res=moduloVector(c,d,BASE);
	for(int i=res.size()-1;i>=0;--i)
		result.long_num.push_back(res[i]);
	return result;
}

LI LI::Primitive(LI a)
{
	int shift=0;
	LI result;
	result.is_fractional=is_fractional;
	//==0
	if(a==result||(*this)==result) return result;
	LI temp;
	if(long_num.size() > a.long_num.size()) long_num.swap(a.long_num);
	for(int i=0;i<(int)long_num.size();++i)
	{
		temp=a*long_num[i];
		for(int j=0;j<i;++j)
			temp.long_num.insert(temp.long_num.begin(),0);
		result=result+temp;
	}
	return result;
}

LI LI::Karatsuba(LI a)
{
	while(long_num.size()<a.long_num.size())	
		long_num.push_back(0);
	while(long_num.size()>a.long_num.size())
		a.long_num.push_back(0);

	if(long_num.size()%2!=0)
	{
		a.long_num.push_back(0);
		long_num.push_back(0);
	}

	int n=long_num.size();
	int m=n/2;

	LI A1,A0,B0,B1;
	LI result;
	
	int i=0;
	for(i=0;i<m;++i)
		A0.long_num.push_back(long_num[i]);
	for(;i<n;++i)
		A1.long_num.push_back(long_num[i]);
	for(i=0;i<m;++i)
		B0.long_num.push_back(a.long_num[i]);
	for(;i<n;++i)
		B1.long_num.push_back(a.long_num[i]);
	
	result = A0.Primitive(B0);
	LI res1 = A1.Primitive(B1);
	LI res2 = ((A0+A1).Primitive((B0+B1)) - result - res1);
	for (int i=0; i<m; ++i)
		res2.long_num.insert(res2.long_num.begin(),0);
	for (int i=0; i<n; ++i)
		res1.long_num.insert(res1.long_num.begin(),0);
	result=result+res1+res2;
	return result;
}

LI LI::ToomCook (LI a)
{
	while(long_num.size()<a.long_num.size())	
		long_num.push_back(0);
	while(long_num.size()>a.long_num.size())
		a.long_num.push_back(0);

	LI result;
	int n=long_num.size();
	const LI cod1("-1"),cod2("-2"),cod3("-3"),cod4("-4");
	stack<int> A;
	stack<LI> U,V,C,W;
	int k,Q,R;
	list<int> q,r;
	//T1
	k=1; Q=4; R=2;
	q.push_back(16); q.push_back(16);
	r.push_back(4); r.push_back(4);
	list<int>::iterator i=q.begin();
	//list<int>::iterator j=i-1;
	while((*i)+*(++i)<n)
	{
		++k; 
		Q+=R;
		if((R+1)*(R+1)<=Q)++R;
		q.push_back((int)pow(2.,Q));
		r.push_back((int)pow(2.,R));
	}
	//T2
	int p=(*i)+*(--i);
	//i at end
	while((int)long_num.size()<p)
	{
		a.long_num.push_back(0);
		long_num.push_back(0);
	}
	C.push(cod1);
	C.push(a);
	C.push(*this);
	
	while(true)
	{
		while(k>1)
		{
			--k;
			i=q.end(); --i;
			int s=*i;//q=qk
			p=(*i)+*(--i);
			i=r.end();--i;
			int t=*i;//r=rk
			LI top=C.top();
			//розбити чило в вершины на т+1 чисел довжини с
			vector<LI> u;//разбиение
			for(int i=0, k=0; i<=t; ++i)
			{
				u.push_back(LI());
				for(int j=0; j<s && k<(int)top.long_num.size(); ++j)
				{
					u[i].long_num.push_back(top.long_num[k++]);
				}
				while(int(u[i].long_num.size())<s) 
					u[i].long_num.insert(u[i].long_num.begin(),0);
			}
	
			LI K;
			LI i0;
			for(int j=0; j<=t+t; ++j)
			{
				int l=1;
				K.long_num.clear();
				for(int i=u.size()-1; i>=0; --i)
				{
					K=u[i]*l+K;
					l*=j;
				}
				U.push(K);
			}
			C.pop();
			top=C.top();
			u.clear();
			for(int i=0,k=0; i<=t; ++i)
			{
				u.push_back(LI());
				for(int j=0; j<s && k<(int)top.long_num.size(); ++j)
				{
					u[i].long_num.push_back(top.long_num[k++]);
				}
				while(int(u[i].long_num.size())<s) 
					u[i].long_num.insert(u[i].long_num.begin(),0);
			}
			K.long_num.clear();
			for(int j=0; j<=t+t; ++j)
			{
				int l=1;
				K.long_num.clear();
				for(int i=u.size()-1; i>=0; --i)
				{
					K=u[i]*l+K;
					l*=j;
				}
				V.push(K);
			}
			C.pop();
			C.push(cod2);
			C.push(U.top());
			C.push(V.top());
			U.pop();V.pop();
			for(int i=0;i<t+t;++i)
			{
				C.push(cod3);
				C.push(U.top());
				C.push(V.top());
				U.pop();
				V.pop();
			}
		}

	LI w,u,v;
	u=C.top(); C.pop();
	v=C.top(); C.pop();		//
	w=u.Primitive(v);
		while(k<=1)
		{
			while(C.top()==cod3)
			{//T6
				C.pop();
				W.push(w);
				u=C.top();C.pop();
				v=C.top();C.pop();
				w=u.Primitive(v);
			}
			if(C.top()==cod2)
			{
				k++;
				W.push(w);
				//T7
				C.pop();
				i=q.end();--i;
				int s=*i;
				p=(*i)+*(--i);
				i=r.end();--i;
				int t=*i;
				//trasfer stack W in vector
				vector<LI> l;
				while(!W.empty())
				{
					l.push_back(W.top());
					W.pop();
				}
				LI temp;
				for(int i=0; i<(int)l.size()/2; ++i)
				{
					temp=l[i];
					l[i]=l[l.size()-1-i];
					l[l.size()-1-i]=temp;
				}
				//top of stack in the last elem of l (step 13)
				for(int j=1; j<=t+t; ++j)
				{
					for(int k=t+t; k>=j; --k)
					{
						//T7
						l[k]=(l[k]-l[k-1])/j;
					}
				}
				//T8
				for(int j=t+t-1; j>0; --j)
				{
					for(int k=j; k<t+t; ++k)
					{
						l[k]=l[k]-l[k+1]*j;
					}
				}
				//T9 (step 16)
				int k=1;
				w.long_num.clear();
				for(int i=l.size()-1; i>=0; --i)
				{
					w=l[i]*k+w;//??
					k*=s;
				}
			}
			if(C.top()==cod1)
			{//T10
				return w;
			}
			if(C.top()==cod3||C.top()==cod2) k=1;
		}
	}  
}

LI LI::Schonhage(LI a)
{
	while(long_num.size()<a.long_num.size())	
		long_num.push_back(0);
	while(long_num.size()>a.long_num.size())
		a.long_num.push_back(0);

	int q=1,p=26;
	if((int)long_num.size()<=p) return Primitive(a);
	while((int)long_num.size()>p)
	{
		q=q+q+q-1;
		p=18*q+8;
	}
	int m[6];
//	if(q==1) 
	m[0]=1;
	for(int i=0;i<6*q;++i)
		m[0]=m[0]*2;
	//2^6q
	m[0]=m[0]/2-1; m[1]=m[0]*2-1;	m[2]=m[0]*4-1;	
	m[3]=m[0]*8-1; m[4]=m[0]*32-1;  m[5]=m[0]*128-1;

	int u[6],v[6],w[6];
	for(int i=0;i<6;++i)
	{
		u[i]=(*this)%m[i];
		v[i]=a%m[i];
		w[i]=(u[i]*v[i])%m[i];
	}

	//find w
	int c[6][6];
	for(int i=0; i<6; ++i)
		for(int j=0; j<6; ++j)
		{
			int t=0;
			if(i!=j) for(;(t*m[i])%m[j]!=1;t++);
			c[i][j]=t;
		}	
	int W[6];
	W[0]=w[0];
	for(int i=1; i<6; ++i)
	{
		W[i]=((w[i]-W[0])*c[1][i]%m[i]+m[i])%m[i];//так, бо числа можуть бути в≥д'Їмн≥
		for(int j=1;j<i;++j)	
			W[i]=((W[i]-W[j])*c[j][i]%m[i]+m[i])%m[i];
	}

	LI result;
 	result.long_num.push_back(W[5]);

	for(int i=4;i>=0;--i)
		result=result*m[i]+W[i];

	return result;
}

vector<int> summ(vector<int> a,vector<int> b,const int& base)
{
	int shift=0;
	int temp_result;
	vector<int> result;
	while(a.size()<b.size()) a.push_back(0);
	while(a.size()>b.size()) b.push_back(0);
	for(int i=0; i<(int)a.size(); i++)
	{
		temp_result=a[i]+b[i]+shift;
		result.push_back(temp_result%base);
		shift=temp_result/base;
	}
	if(shift!=0)result.push_back(shift);
	return result;
}

LI LI::mod(const int & n)
{
	vector<int> k;
	k=toBin1();
	vector<int> k1,k2;
	while((int)k.size()>n)
	{
		k1.clear();
		k2.clear();
		for(int i=0; i<n; ++i)
			k1.push_back(k[i]);
		for(int i=n; i<(int)k.size(); ++i)
			k2.push_back(k[i]);
		k.clear();
		k=summ(k1,k2,2);
	}
	LI result;
	result.inBin1(k);
	return result;
}

bool LI::Lemer()
{
	int n=-1;
	if(!isMersenne(n))return false;
	LI s(4);
	for(int i=0; i<n-2; ++i)
		s=(s.Karatsuba(s)-2).mod(n);
	if(s==0) return true;
	LI M(2);
	M=(M^n)-1;
	if(s==M)return true;
	return false;
}

bool LI::isMersenne(int & p)
{
	if(long_num.size()==0)return false;
	LI x=(*this)+1;
	p=1;
	while(x!=2)
	{
		if(x.long_num[0]%2==1)return false;
		x=x/2;
		++p;
	}
	return true;
}

bool isPow2(const int &x)
{
	return ((x!=0)&&((x&(~x+1))==x));
}

LI LI::CookDivide(LI v)
{
	int m=long_num.size();
	int n=v.long_num.size();
	LI rеsult;
	int k=1,j=0;
	while(k<(BASE*n) || k<(BASE*m-BASE*n)) {k*=BASE;++j;}//2
	k=k/BASE;//2
	rеsult=(*this)/v;
	if(k!=BASE-1) return rеsult;

	//дописуЇмо к-n нул≥в
	v.long_num.insert(v.long_num.begin(),k-n,0);
	
	LI a=BASE;//2
	LI d;
	LI a1;
	LI temp;

	for(int i=1;i<j;++i)
	{
		temp=v;
		temp.long_num.erase(temp.long_num.begin(),temp.long_num.begin()+k-i-i);
		a.LI_delete_zeros();
		a1=a;
		a1.long_num.insert(a1.long_num.begin(),3*(int)pow((double)BASE,i),0);//2
		d=a1-a*a*temp;
		a=d;
		a.long_num.erase(a.long_num.begin(),a.long_num.begin()+i-1);
	}

	d.long_num.erase(d.long_num.begin(),d.long_num.end()-k-k);
	a.long_num.erase(a.long_num.begin(),a.long_num.end()-k);
	a1=a;
	a1.long_num.insert(a1.long_num.begin(),k+k,0);

	d=a1-a*a*v;
	LI l=1;
	l.long_num.insert(l.long_num.begin(),k+k-1,0);
	a=(d+l);
	a.long_num.erase(a.long_num.begin(),a.long_num.begin()+k+k-1);
	LI result;
	l=1;
	l.long_num.insert(l.long_num.begin(),k+n-1,0);
	result=a*(*this)+l;
	result=result+l;
	result.long_num.erase(result.long_num.begin(),result.long_num.begin()+k+n-1);
	return rеsult;
}

vector<complex<double>> FFT(vector<complex<double>> long_num)
{
	if(long_num.size()==1) return vector<complex<double>>(1,long_num[0]);
	while(!isPow2(long_num.size()))long_num.push_back(0);	
	int n=long_num.size();
	vector<complex<double>> A;
	vector<complex<double>> B;
	for(int i=0;i<n;++i)
	{
		if(i%2==0) A.push_back(long_num[i]);
		else B.push_back(long_num[i]);
	}
	vector<complex<double>> w;
	double temp;
	for(int i=0;i<n;++i)
	{
		temp=2*pi*i/n;
		w.push_back(complex<double>(cos(temp),sin(temp)));
	}

	vector<complex<double>> A1=FFT(A);
	vector<complex<double>> B1=FFT(B);
	vector<complex<double>> result;
	for(int i=0;i<n;++i)
	{
		result.push_back(A1[i%(n/2)]+w[i]*B1[i%(n/2)]);
	}
	return result;
}

vector<complex<double>> FFT_REV(vector<complex<double>> result)
{
	vector<complex<double>> long_num=FFT(result);
	int n=long_num.size();
	for(int i=0; i<(int)long_num.size(); ++i)
		long_num[i]/=n;
	complex<double> temp=0;
	for(int i=1; i<n/2; ++i)
	{
		temp=long_num[i];
		long_num[i]=long_num[n-i];
		long_num[n-i]=temp;
	}
	return long_num;
}

long long round(double x)
{
	long long d=(long long) x;
	double eps=x-d;
	if(eps<0.5) return d;
	return d+1;
}

LI LI::SchonhageShtrassen(LI a){

	while(long_num.size()<a.long_num.size())	
		long_num.push_back(0);
	while(long_num.size()>a.long_num.size())
		a.long_num.push_back(0);

	vector<complex<double>> u,v;

	for(int i=0;i<(int)long_num.size();++i)
		u.push_back(complex<double>(long_num[i],0));
	for(int i=0;i<(int)a.long_num.size();++i)
		v.push_back(complex<double>(a.long_num[i],0));
	

	for(int i=0;i<(int)long_num.size();++i)
		u.push_back(complex<double>(0,0));
	for(int i=0;i<(int)a.long_num.size();++i)
		v.push_back(complex<double>(0,0));
	
	vector<complex<double>> u1,v1,w;
	u1=FFT(u);v1=FFT(v);
	for(int i=0;i<v1.size();++i){
		w.push_back(u1[i]*v1[i]);
	}
	u.clear();
	u=FFT_REV(w);
	LI result;
	for(int i=0;i<u.size();++i)
		result.long_num.push_back(round(u[i].real()));
	int shift=0;
	for(int i=0;i<result.long_num.size();++i)
	{
		result.long_num[i]+=shift;
		if(result.long_num[i]>BASE)
		{
			shift=result.long_num[i]/BASE;
			result.long_num[i]%=BASE;
		}
		else shift=0;
	}
	while(shift!=0)
	{
		if(shift>BASE)
		{
			result.long_num.push_back(shift%BASE);
			shift/=BASE;}
		else 
		{
			result.long_num.push_back(shift);
			shift =0;
		}
	}
	return result;
}

LI LI::FastPowMod(const LI &n, const LI &m)
{
	vector<int> k;
	//n in 2
	k=n.toBin1();
	//
	LI s=(*this);
	for(int i=(int)k.size()-2;i>=0;--i)
		if(k[i]==0)s=(s*s)%m;
		else {s=(s*s)%m; s=(s*(*this))%m;}
	return s;
}

bool LI::MillerRabin(const int &r)
{
	if(long_num.size()==0)return false;
	LI t=(*this)-1;
	if(t.long_num[0]%2==1)return false;
	int s=1;
	while(t%2==0)
	{
		t=t/2;
		++s;
	}
	srand(unsigned(time(NULL)));
	LI x,a;
	LI m1=(*this);
	m1.long_num[0]-=1;
	for(int i=0;i<r;++i)
	{
		int size;
		if(long_num.size()==1)
		{
			size=1;
			a=rand()%(long_num[0]-2)+2;
		}
		else if(long_num.size()==2)
		{
			size = rand()%2+1;
			a.long_num.push_back(rand()%(long_num[0]-2)+2);
		}
		else
		{
			size=rand()%(long_num.size()-2)+2;
			for(int i=0;i<size;++i)
				a.long_num.push_back((rand()%BASE));
		}		
		x=a.FastPowMod(t,(*this));
		if(x==1 || x==m1) continue;
		int j=0;
		for(j=0; j<s-1; ++j)
		{
			x=(x*x)%(*this);
			if(x==1)return false;
			else if(x==m1)break;
		}
		if(j!=s-1)continue;
		return false;
	}
	return true;
}

bool LI::SolovejShtrassen(int k)
{
	int base=BASE;
	vector<int> m=long_num;
	reverse(m);
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
