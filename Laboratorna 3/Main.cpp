#include <iostream>
#include "LI.h"
#include "Rational.h"
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include "Matrix.h"

using namespace std;

//програма буде обчислювати трохи довго (в районі 25 секунд, чи менше), бо запущені відразу всі алгоритми, крім того примітивне множення,
//метод Гаусса для раціональних і дійсних, перевід у дійсні і назад

int main()
{
	ifstream fi1("in1.txt");
	ifstream fi2("in2.txt");
	ifstream fi3("in3.txt");
	ifstream fi4("in4.txt");
	ifstream fi5("in5.txt");
	ifstream fi6("in6.txt");
	ifstream fi7("in7.txt");
	ifstream fi8("in8.txt");
	ifstream fi9("in9.txt");
	ofstream fo("out.txt");
	Matrix<Rational> a,b,c;
	fi3>>a;
	fi4>>b;
	double time=(float)clock()/CLOCKS_PER_SEC;
	c=a*b;
	double dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Primitive multiplication (time):\n"<<dif<<endl;
	fo<<"Primitive multiplication:\n";
	fo<<c;
	fo<<endl;
	time=(float)clock()/CLOCKS_PER_SEC;
	c=a.Vinograd(b);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Vinograd multiplication (time):\n"<<dif<<endl;
	fo<<"Vinograd multiplication:\n";
	fo<<c;
	fo<<endl;
	time=(float)clock()/CLOCKS_PER_SEC;
	c=a.Shtrassen(b);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Shtrassen multiplication (time):\n"<<dif<<endl;
	fo<<"Shtrassen multiplication:\n";
	fo<<c;
	fo<<endl;
	Matrix<Rational> d;
	vector<Rational> ans1;
	vector<Rational> right1;
	vector<Rational> answer1;
	fi1>>d;
	fi2>>ans1;
	right1=d.GenRight(ans1);
	fo<<right1;
	fo<<" - right for Rational Gauss\n";
	time=(float)clock()/CLOCKS_PER_SEC;
	answer1=d.Gauss(right1);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Gauss (Rational) (time):\n"<<dif<<endl;
	fo<<"Gauss:\n";
	fo<<answer1;
	fo<<endl;
	Matrix<double> e;
	vector<double> ans2;
	vector<double> right2;
	vector<double> answer2;
	fi5>>e;
	fi6>>ans2;
	right2=e.GenRight(ans2);
	fo<<right2;
	fo<<" - right for Double Gauss, and others\n";
	time=(float)clock()/CLOCKS_PER_SEC;
	answer2=e.Gauss(right2);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Gauss (double) (time):\n"<<dif<<endl;
	fo<<"Gauss:\n";
	fo<<answer2;
	fo<<endl;
	time=(float)clock()/CLOCKS_PER_SEC;
	answer2=e.Givens(right2);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Givens (double) (time):\n"<<dif<<endl;
	fo<<"Givens\n";
	fo<<answer2;
	fo<<endl;
	time=(float)clock()/CLOCKS_PER_SEC;
	answer2=e.HouseHolder(right2);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"HouseHolder (double) (time):\n"<<dif<<endl;
	fo<<"HouseHolder:\n";
	fo<<answer2;
	fo<<endl;
	vector<double> f(5,1);
	time=(float)clock()/CLOCKS_PER_SEC;
	answer2=e.Grad(right2,f,0.01);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Gradients (double) (time):\n"<<dif<<endl;
	fo<<"Gradients:\n";
	fo<<answer2;
	fo<<endl;
	time=(float)clock()/CLOCKS_PER_SEC;
	answer2=e.GramShmidt(right2);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"GramShmidt (double) (time):\n"<<dif<<endl;
	fo<<"GramShmidt:\n";
	fo<<answer2;
	fo<<endl;
	time=(float)clock()/CLOCKS_PER_SEC;
	answer2=e.Holeckyj(right2);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Holeckyj (double) (time):\n"<<dif<<endl;
	fo<<"Holeckyj:\n";
	fo<<answer2;
	fo<<endl;
	Matrix<Rational> u;
	fi7>>u;
	Matrix<Rational> g;
	time=(float)clock()/CLOCKS_PER_SEC;
	g=u^(-1);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Jordan (rational) (time):\n"<<dif<<endl;
	fo<<"Jordan:\n";
	fo<<g;
	fo<<endl;
	Matrix <double> h;
	Matrix <double> t;
	fi8>>h;
	time=(float)clock()/CLOCKS_PER_SEC;
	t=h.Jacobi();//оскільки алгоритм не зовсім точний, то він обчислює позадіагональні елементи не як чисті нулі,
	             //а як значення дуже близькі до нуля
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"Jacobi (double) (time):\n"<<dif<<endl;
	fo<<"Jacobi:\n";
	fo<<t;
	fo<<endl;
	fo<<fixed;
	fo.precision(10);
	Rational p;
	fi9>>p;
	double k;
	time=(float)clock()/CLOCKS_PER_SEC;
	k=p.toDouble();
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"To double (time):"<<dif<<endl;
	fo<<"double(p):\n";
	fo<<k;
	fo<<endl;
	Rational m;
	time=(float)clock()/CLOCKS_PER_SEC;
	m=toRational(k);
	dif=(float)clock()/CLOCKS_PER_SEC-time;
	cout<<"To rational (time):"<<dif<<endl;
	fo<<"Back to rational:\n";
	fo<<m;
	fo.close();
	fi1.close();
	fi2.close();
	fi3.close();
	fi4.close();
	fi5.close();
	fi6.close();
	fi7.close();
	fi8.close();
	fi9.close();
	cin.get();
	return 0;
}
