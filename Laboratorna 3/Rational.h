#ifndef _RATIONAL_CLASS
#define _RATIONAL_CLASS
#include <iostream>
#include "LI.h"

class Rational	//	n/m
{
	bool negative;
public:
	LI n;
	LI m;
	Rational(LI a, LI b,const bool &_negative=false);//+
	Rational(const LI &N,const bool &_negative=false);//+
	Rational ();//+
	friend bool operator>>(istream &i, Rational &a);//+
	friend bool operator<<(ostream &o, Rational a);//+
	void formalise();//+
	friend Rational abs(Rational a);
	bool operator== (const Rational & a);//+
	bool operator<  (Rational a);//+
	bool operator<= (Rational a);//+
	bool operator>  (Rational a);//+
	bool operator>= (Rational a);//+
	bool operator!= (Rational a);//+
	Rational operator=  (const Rational &a);//+
	Rational operator+ (Rational a);//+
	Rational operator- (Rational a);//+
	Rational operator* (int a);//+
	Rational operator* (Rational a);//+
	Rational operator/ (Rational a);//+
	Rational operator/ (int a);//+
	double toDouble();
	friend Rational toRational(double x);
};


Rational sqrt(Rational x);
Rational sin(Rational x);
Rational cos(Rational x);
Rational tan(Rational x);
Rational atan(Rational x);

#endif