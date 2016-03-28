#include "Rational.h"
#include "Support Functions.h"

LI GCD(LI a, LI b)
{
	if (a==0) return b;
	if (b==0) return a;
	LI c;
	while (b!=0)
	{
		a=a%b;
		c=a;
		a=b;
		b=c;
	}
	return a;
}

LI LCM(LI a, LI b)
{
	return (a/GCD(a,b))*b;
}

bool operator>> (istream &i, Rational &a)
{
	LI N,M;
	bool negative=false;
	char sign;
	i>>sign;
	if(sign=='-'){negative=true;}
	else i.unget();
	if(!(i>>N))return false;
	sign=' ';
	i>>sign;
	if(sign=='-'){negative=!negative;}
	else if(sign=='/')
		 {
			if(!(i>>M))return false;
			Rational x(N,M,negative);
			a=x;
			return true;
		 }
		else i.unget();
	Rational p(N,negative);
	a=p;
	return true;
}

bool operator<< (ostream &o, Rational a)
{
	if (a.n == 0) return o<<'0';
	if (a.negative)
		if(!(o<<'-'))return false;
	if(!(o<<a.n))return false;
	if(a.m!=1)
	{
		if(!(o<<'/'))return false;
		if(!(o<<a.m))return false;
	}
	return true;
}

Rational::Rational()
{
	negative=false;
	n=1;
	m=1;
}

Rational::Rational(LI a, LI b, const bool &_negative)
{
	negative=_negative;
	if(b!=0)
	{
		n=a;
		m=b;
	}
	else
	{
		n=1;
		m=1;
	}
	formalise();
}

Rational::Rational(const LI &N,const bool &_negative)
{
	n=N;
	m=1;
	negative=_negative;
}

void Rational::formalise()
{
	LI gcd=GCD(n,m);
	if(gcd==1) return;
	n=n/gcd;
	m=m/gcd;
	return;
}

Rational abs(Rational a)
{
	Rational temp=a;
	temp.negative=false;
	return temp;
}

Rational Rational::operator=(const Rational &a)
{
	n=a.n;
	m=a.m;
	negative=a.negative;
	return (*this);
}

bool Rational::operator== (const Rational & a)
{
	return (negative==a.negative && n==a.n && m==a.m);
}

bool Rational::operator!= (Rational a)
{
	return !((*this)==a);
}

bool Rational::operator< (Rational a)
{
	if (negative && !a.negative) return true;
	if (!negative && a.negative) return false;
	LI gcd=GCD(m, a.m);
	LI N;
	N=n*(a.m/gcd);
	a.n=a.n*(m/gcd);
	if(negative==false)
		return N<a.n;
	return a.n<N;
}

bool Rational::operator> (Rational a)
{
	if (negative && !a.negative) return false;
	if (!negative && a.negative) return true;
	LI gcd=GCD(m, a.m);
	LI N;
	N=n*(a.m/gcd);
	a.n=a.n*(m/gcd);
	if(negative==false)
		return N>a.n;
	return a.n>N;
}

bool Rational::operator<=(Rational a)
{
	return (*this)<a||(*this)==a;
}

bool Rational::operator>=(Rational a)
{
	return !((*this)<a);
}

Rational Rational::operator+(Rational a)
{
	if(negative && !a.negative) return a+(*this);
	LI gcd=GCD(m,a.m);
	LI N;
	N=n*(a.m/gcd);
	a.n=a.n*(m/gcd);
	if(a.negative && !negative)
		if(N==a.n) return Rational(0);
		else if(a.n<N)
			return Rational(N-a.n,(m/gcd)*a.m);
		else return Rational(a.n-N,(m/gcd)*a.m,true);
	return Rational(N+a.n, (m/gcd)*a.m, negative);
}

Rational Rational::operator-(Rational a)
{
	Rational res;
	res=(*this);
	a.negative=!a.negative;
	return res+a;
}

Rational Rational::operator* (int a)
{
	Rational res;
	if (negative) res.negative=true;
	if(a<0)
	{
		if (negative) res.negative=false;
		else res.negative=true;
		a=-a;
	}
	res.n=n*a;
	res.m=m;
	res.formalise();
	return res;
}

Rational Rational::operator* (Rational a)
{
	Rational res;
	res.negative=(negative||a.negative)&&(!negative||!a.negative);
	res.n=n*a.n;
	res.m=m*a.m;
	res.formalise();
	return res;
}

Rational Rational::operator/ (Rational a)
{
	Rational res1;
	res1=(*this);
	Rational res2(a.m,a.n,a.negative);
	return res1*res2;
}

Rational Rational::operator/ (int a)
{
	Rational res;
	if(a<0)
	{
		if (negative) res.negative=false;
		else res.negative=true;
		a=-a;
	}
	res.n=n;
	res.m=m*a;
	res.formalise();
	return res;
}

double Rational::toDouble()
{
	if(negative) return -n.DoubleDivide(m);
	return n.DoubleDivide(m);
}

Rational toRational(double x)
{
	Rational R;
	if (x<=2147483647.)
	{
		long long l=x;
		double dif=x-l;
		Rational a(l);
		dif*=1e+5;
		R.n=LI(long long(dif));
		R.m=1e+5;
		R=R+a;
	}
	else//
	{
		vector<int> a;
		int count=0;
		while (x>10)
		{
			x/=10;
			++count;
		}
		while(count!=-1)
		{
			a.push_back(int(x));
			x-=int(x);
			x*=10;
			count--;
		}
		count=0;
		while(count!=15)
		{
			a.push_back(int(x));
			x-=int(x);
			x*=10;
			count++;
		}
		vector<int> b;
		b=convertInBase(a,LI::BASE);
		a.clear();
		for(int i=(int)b.size()-1; i!=-1; --i)
			a.push_back(b[i]);
		R.n.long_num=a;
		a.clear();
		a.push_back(1);
		for (int i=0; i<15; ++i)
			a.push_back(0);
		b=convertInBase(a,LI::BASE);
		a.clear();
		for(int i=(int)b.size()-1; i!=-1; --i)
			a.push_back(b[i]);
		R.m.long_num=a;
	}
//	R.formalise();
	return R;
}

Rational sqrt(Rational x)
{
	Rational result=x;
	if(x<Rational(0)) x=x*(-1);
	for(int i=0;i<4;++i)
		result=(result+(x)/result)/Rational(2);
	return result;
}

Rational sin(Rational x)
{
	return x - x*x*x/Rational(6) + x*x*x*x*x/Rational(120);
}
Rational cos(Rational x)
{
	return Rational(1) - x*x/Rational(2) + x*x*x*x/Rational(4);
}
Rational tan(Rational x)
{
	return x+x*x*x/Rational(3);
}
Rational atan(Rational x)
{
	return x - x*x*x/Rational(3);
}