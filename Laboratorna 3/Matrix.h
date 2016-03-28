#ifndef Matrix____
#define Matrix____

#include "Rational.h"
#include <iostream>
#include <fstream>

using namespace std;

template<class T>
class Matrix
{
	vector<vector<T>> matrix;
public:
	Matrix<T>();
	Matrix<T>(int _size, T val=(T)0);
	int size()const{return (int) matrix.size();}
	template<class T>
	friend Matrix<T> unitMatr(int _size);
	template<class T>
	friend bool operator>>(istream &i, Matrix<T> &a);
	template<class T>
	friend bool operator<<(ostream &o, Matrix<T> a);
	Matrix<T> operator=(Matrix<T> a);
	Matrix<T> operator+(Matrix<T> a);
	Matrix<T> operator-(Matrix<T> a);
	Matrix<T> operator*(int a);
	vector<T> operator*(vector<T> a);//множення на вектор-стовпчик
	Matrix<T> operator*(Matrix<T> a);
	vector<T>& operator[](int i);
	Matrix<T> transpose();
	bool operator==(Matrix<T> a);
	bool operator!=(Matrix<T> a);
	Matrix<T> Vinograd(Matrix<T> b);
	Matrix<T> Shtrassen(Matrix<T> b);
	vector <T> GenRight(vector <T> answer);
	vector <T> Gauss(vector<T> a);
	Matrix<T> operator^(int a);
	vector<T> Givens(vector<T> right);
	vector<T> HouseHolder(vector<T> right);
	vector<T> Grad(vector<T> right, vector <T> x, T eps);
	vector<T> GramShmidt(vector<T> right);
	vector<T> Holeckyj(vector<T> right);
	Matrix<T> Jacobi();
};

//ввід-вивід матриць
template <class T>
bool operator>>(istream &i, Matrix<T> &a)
{
	vector<T> temp;
	T t;
	while(i>>t) temp.push_back(t);
	int size=(int)temp.size();
	size=(int)sqrt((double)size);
	vector<T> y(0);
	int k=0;
	for (int i=0; i<size; ++i)
	{
		a.matrix.push_back(y);
		for (int j=0; j<size; ++j)
			a.matrix[i].push_back(temp[k++]);
	}
	return true;
}

template<class T>
bool operator<<(ostream &o, Matrix<T> a)
{
	if (a.size()==0) return false;
	for (int i=0; i<a.size(); ++i)
	{
		for (int j=0; j<a.size(); ++j)
		{
			o.width(15);
			if (!(o<<a.matrix[i][j]) || !(o<<' ')) return false;
		}
		if(!(o<<endl)) return false;
	}
	return true;
}

//ввід-вивід векторів
template<class T>
bool operator>>(istream & i, vector<T> & a)
{
	a.clear();
	T temp;
	while(i>>temp) a.push_back(temp);
	return true;
}

template<class T>
bool operator<<(ostream & o, vector<T> a)
{
	for(int i=0; i<a.size(); ++i)
	{
		o<<a[i];
		o<<' ';
	}
	return true;
}

//деякі допоміжні операції для векторів
template<class T>
vector<T> operator*(vector<T> a, T b)
{
	for (int i=0; i<(int)a.size(); ++i)
		a[i]=a[i]*b;
	return a;
}

template<class T>
vector<T> operator+(vector<T> a, vector<T> b)
{
	vector<T> result;
	if (a.size()!=b.size()) return result;
	for (int i=0; i<(int)a.size(); ++i)
		result.push_back(a[i]+b[i]);
	return result;
}

template<class T>
vector<T> operator-(vector<T> a, vector<T> b)
{
	vector<T> result;
	if (a.size()!=b.size()) return result;
	for (int i=0; i<(int)a.size(); ++i)
		result.push_back(a[i]-b[i]);
	return result;
}

template<class T>
T operator*(vector<T>a, vector<T> b)
{
	T result(0);
	if (a.size()!=b.size() || a.size()==0) return result;
	for (int i=0; i<(int)a.size(); ++i)
		result=result+a[i]*b[i];
	return result;
}

//матриці основні методи
template<class T>
Matrix<T>::Matrix<T>()
{
	matrix.clear();
}

template<class T>
Matrix<T>::Matrix<T>(int _size, T val)
{
	if (_size==0) matrix.clear();
	vector <T> y(0);
	for (int i=0; i<_size; ++i)
	{
		matrix.push_back(y);
		for (int j=0; j<_size; ++j)
			matrix[i].push_back(val);
	}
}

template<class T>
Matrix<T> Matrix<T>::operator=(Matrix<T> b)
{
	if (b.size()==0) return Matrix<T>();
	matrix=b.matrix;
	return (*this);
}

template<class T>
Matrix<T> unitMatr(int _size)
{
	Matrix<T> E;
	if (_size==0) E.matrix.clear();
	vector <T> y(0);
	for (int i=0; i<_size; ++i)
	{
		E.matrix.push_back(y);
		for (int j=0; j<_size; ++j)
			if (i==j) E.matrix[i].push_back((T)1);
			else E.matrix[i].push_back((T)0);
	}
	return E;
}

template<class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> a)
{
	if(size()==0) return a;
	if(a.size()==0) return (*this);
	if(size()!=a.size()) return Matrix<T>();
	Matrix<T> result(size());
	for(int i=0; i<size(); ++i)
		for(int j=0; j<size(); ++j)
			result.matrix[i][j]=matrix[i][j]+a.matrix[i][j];
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> a)
{
	if(size()==0) return a*(-1);
	if(a.size()==0) return (*this);
	Matrix<T> result(size());
	for(int i=0; i<size(); ++i)
		for(int j=0; j<size(); ++j)
			result.matrix[i][j]=matrix[i][j]-a.matrix[i][j];
	return result;
}

template<class T>
Matrix<T> Matrix<T>::operator*(int a)
{
	if (a==0) return Matrix<T>(size());
	Matrix<T> result(size());
	for(int i=0; i<size(); ++i)
		for(int j=0; j<size(); ++j)
			result.matrix[i][j]=matrix[i][j]*a;
	return result;
}

//множення на вектор стовпчик справа
template<class T>
vector<T> Matrix<T>::operator*(vector<T> a)
{
	
	if (a.size()==0) return vector <T> ();
	if (size()==0||size()!=int(a.size())) return vector<T> ();
	vector<T> result (size(),(T)0);
	for (int i=0; i<size(); ++i)
	{
		T temp=0;
		for (int j=0; j<size(); ++j)
		{
			temp=matrix[i][j]*a[j];
			result[i]=result[i]+temp;
			temp=(T)0;
		}
	}
	return result;
}

//звичайне множення
template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> a)
{
	int n=size();
	if (a.size()==0) return Matrix<T> ();
	if (n==0||n!=int(a.size())) return Matrix<T> ();
	Matrix<T> b;
	b=(*this);
	Matrix<T> result(n);
	for (int i=0; i<n; ++i)
	{
		T temp=0;
		for (int j=0, k=0; k<n; ++j)
		{
			if (j==n){j=0;++k;}
			if (k==n) break;
			temp=b.matrix[i][j]*a.matrix[j][k];
			result.matrix[i][k]=result.matrix[i][k]+temp;
			temp=(T)0;
		}
		temp=(T)0;
	}
	return result;
}

template<class T>
vector<T>& Matrix<T>::operator[](int i)
{
	return matrix[i];
}

template<class T>
Matrix<T> Matrix<T>::transpose()
{
	if (size()==0) return Matrix<T> ();
	if (size()==1) return (*this);
	Matrix<T> result(size());
	for (int i=0; i<size(); ++i)
		for (int j=0; j<size(); ++j)
			result.matrix[i][j]=(*this)[j][i];
	return result;
}

template<class T>
bool Matrix<T>::operator==(Matrix<T> a)
{
	if (size()!=a.size()) return false;
	for (int i=0; i<size(); ++i)
		for (int j=0; j<size(); ++j)
			if (matrix[i][j]!=a.matrix[i][j]) return false;
	return true;
}

template<class T>
bool Matrix<T>::operator!=(Matrix<T> a)
{
	if (size()!=a.size()) return true;
	for (int i=0; i<size(); ++i)
		for (int j=0; j<size(); ++j)
			if (matrix[i][j]!=a.matrix[i][j]) return true;
	return false;
}

//алгоритми
template<class T>
Matrix<T> Matrix<T>::Vinograd(Matrix<T> b)
{
	if (size()!=b.size() || size()==0 || b.size()==0) return Matrix<T>();
	int sz=size();
	Matrix<T> a=(*this);
	bool is_zeros=false;
	if (sz%2)
	{
		is_zeros=true;
		++sz;
		vector <T> temp (sz,T(0));
		b.matrix.push_back(temp);
		a.matrix.push_back(temp);
		for (int i=0; i<sz-1; ++i)
		{
			a.matrix[i].push_back(T(0));
			b.matrix[i].push_back(T(0));
		}
	}
	T *A=new T();
	T *B=new T();
	T *C=new T();
	T *temp=new T();
	Matrix<T> result;
	vector<T> y(0);
	for (int i=0; i<sz; ++i)
	{
		(*A)=(T)0;
		result.matrix.push_back(y);
		for (int k=0; k<sz/2; ++k)
		{
			(*temp)=a.matrix[i][k+k]*a.matrix[i][k+k+1];
			(*A)=(*A)+(*temp);
		}
		for (int j=0; j<sz; ++j)
		{	
			(*B)=(T)0;
			(*C)=(T)0;
			for (int k=0; k<sz/2; ++k)
			{
				(*temp)=b.matrix[k+k][j]*b.matrix[k+k+1][j];
				(*B)=(*B)+(*temp);
				(*temp) = (a.matrix[i][k+k]+b[k+k+1][j]) * (a.matrix[i][k+k+1]+b[k+k][j]);
				(*C)=(*C)+(*temp);
			}
			(*C)=(*C)-(*A)-(*B);
			result.matrix[i].push_back(*C);
		}
	}
	if(is_zeros)
	{
		result.matrix.pop_back();
		for (int i=0; i<(int)result.matrix.size(); ++i)
			result.matrix[i].pop_back();
	}
	delete A;
	delete B;
	delete C;
	delete temp;
	return result;
}

template<class T>
Matrix<T> Matrix<T>::Shtrassen(Matrix<T> b)
{
	if (size()!=b.size() || size()==0 || b.size()==0) return Matrix<T>();
	int size_=size();
	Matrix<T> a=(*this);
	bool is_zeros=false;
	if (size_%2)
	{
		is_zeros=true;
		++size_;
		vector<T> temp(size_,T(0));
		b.matrix.push_back(temp);
		a.matrix.push_back(temp);
		for (int i=0; i<size_; ++i)
		{
			a.matrix[i].push_back(T(0));
			b.matrix[i].push_back(T(0));
		}
	}
	Matrix<T>* A11; Matrix<T> *A12; Matrix<T> *A21; Matrix<T> *A22;
	Matrix<T>* B11; Matrix<T> *B12; Matrix<T> *B21; Matrix<T> *B22;
	A11=new Matrix<T>();
	A12=new Matrix<T>();
	A21=new Matrix<T>();
	A22=new Matrix<T>();
	B11=new Matrix<T>();
	B12=new Matrix<T>();
	B21=new Matrix<T>();
	B22=new Matrix<T>();
	vector<T> temp;
	for (int i=0; i<size_/2; ++i)
	{
		temp.clear();
		for(int j=0; j<size_/2; ++j)
			temp.push_back(a.matrix[i][j]);
		(*A11).matrix.push_back(temp);
	}
	for (int i=0; i<size_/2; ++i)
	{
		temp.clear();
		for(int j=size_/2; j<size_; ++j)
			temp.push_back(a.matrix[i][j]);
		(*A12).matrix.push_back(temp);
	}
	for (int i=size_/2; i<size_; ++i)
	{
		temp.clear();
		for(int j=0; j<size_/2; ++j)
			temp.push_back(a.matrix[i][j]);
		(*A21).matrix.push_back(temp);
	}
	for (int i=size_/2; i<size_; ++i)
	{
		temp.clear();
		for(int j=size_/2; j<size_; ++j)
			temp.push_back(a.matrix[i][j]);
		(*A22).matrix.push_back(temp);
	}
	temp.clear();
	for (int i=0; i<size_/2; ++i)
	{
		temp.clear();
		for(int j=0; j<size_/2; ++j)
			temp.push_back(b.matrix[i][j]);
		(*B11).matrix.push_back(temp);
	}
	for (int i=0; i<size_/2; ++i)
	{
		temp.clear();
		for(int j=size_/2; j<size_; ++j)
			temp.push_back(b.matrix[i][j]);
		(*B12).matrix.push_back(temp);
	}
	for(int i=size_/2; i<size_; ++i)
	{
		temp.clear();
		for(int j=0; j<size_/2; ++j)
			temp.push_back(b.matrix[i][j]);
		(*B21).matrix.push_back(temp);
	}
	for (int i=size_/2; i<size_; ++i)
	{
		temp.clear();
		for(int j=size_/2; j<size_; ++j)
			temp.push_back(b.matrix[i][j]);
		(*B22).matrix.push_back(temp);
	}
	Matrix<T> *P1; Matrix<T> *P2; Matrix<T> *P3; Matrix<T> *P4;
	Matrix<T> *P5; Matrix<T> *P6; Matrix<T> *P7;
	P1=new Matrix<T>(size_/2);
	P2=new Matrix<T>(size_/2);
	P3=new Matrix<T>(size_/2);
	P4=new Matrix<T>(size_/2);
	P5=new Matrix<T>(size_/2);
	P6=new Matrix<T>(size_/2);
	P7=new Matrix<T>(size_/2);
	(*P1)=((*A11)+(*A22)).Vinograd((*B11)+(*B22));
	(*P2)=((*A21)+(*A22)).Vinograd(*B11);
	(*P3)=(*A11).Vinograd((*B12)-(*B22));
	(*P4)=(*A22).Vinograd((*B21)-(*B11));
	(*P5)=((*A11)+(*A12)).Vinograd(*B22);
	(*P6)=((*A21)-(*A11)).Vinograd((*B11)+(*B12));
	(*P7)=((*A12)-(*A22)).Vinograd((*B21)+(*B22));
	Matrix<T> *C11; Matrix<T> *C12; Matrix<T> *C21; Matrix<T> *C22;
	C11=new Matrix<T>(size_/2);
	C12=new Matrix<T>(size_/2);
	C21=new Matrix<T>(size_/2);
	C22=new Matrix<T>(size_/2);
	(*C11)=(*P1)+(*P4)-(*P5)+(*P7);
	(*C12)=(*P3)+(*P5);
	(*C21)=(*P2)+(*P4);
	(*C22)=(*P1)+(*P3)-(*P2)+(*P6);
	Matrix<T> result;
	temp.clear();
	for (int i=0; i<size_/2; ++i)
	{
		temp.clear();
		for(int j=0; j<size_/2; ++j)
			temp.push_back((*C11).matrix[i][j]);
		for(int j=0; j<size_/2; ++j)
			temp.push_back((*C12).matrix[i][j]);
		result.matrix.push_back(temp);
	}
	for (int i=0; i<size_/2; ++i)
	{
		temp.clear();
		for(int j=0; j<size_/2; ++j)
			temp.push_back((*C21).matrix[i][j]);
		for(int j=0; j<size_/2; ++j)
			temp.push_back((*C22).matrix[i][j]);
		result.matrix.push_back(temp);
	}
	if(is_zeros)
	{
		result.matrix.pop_back();
		for (int i=0; i<(int)result.matrix.size(); ++i)
			result.matrix[i].pop_back();
	}
	delete A11; delete A12; delete A21; delete A22;
	delete B11; delete B12; delete B21; delete B22;
	delete C11; delete C12; delete C21; delete C22;
	delete P1; delete P2; delete P3; delete P4; delete P5; delete P6; delete P7;
	return result;
}

template<class T>
vector<T> Matrix<T>::Gauss(vector<T> a)
{
	if (size()!=a.size() || a.size()==0) return vector<T> ();
	Matrix<T> matr=(*this);
	int n = size();
	for (int i=0; i<size(); ++i)
		matr[i].push_back(a[i]);
	for (int col=0; col<n; ++col)
	{
		int row=0;
		for (row=col; row<n && matr[row][col]==(T)0; ++row);
		if (row==n) return vector<T> ();
		swap(matr.matrix[row],matr.matrix[col]);

		T x=matr[col][col];
		for (int j=0; j<n+1; ++j)
			matr[col][j]=matr[col][j]/x;

		vector<T> temp;
		for (int i=col; i<n; ++i)
		{
			if (i!=col)
			{
				temp=matr[col];
				x=matr[i][col];
				for (int k=col; k<n+1; ++k)
				{
					temp[k]=temp[k]*x;
					matr[i][k]=matr[i][k]-temp[k];
				}
			}
		}
	}

	vector<T> ans(n,0);
	T temp(0);
	ans[n-1]=matr[n-1][n];
	for (int i=n-2; i>=0; --i)
	{
		temp=matr[i][n];
		for (int j=n-1; j>i; --j)
			temp = temp - matr[i][j]*ans[j];
		ans[i]=temp;
	}
	return ans;
}

//генерує праву частину за відомим розв'язком
template<class T>
vector<T> Matrix<T>::GenRight(vector<T> answer)
{
	if (size()!=answer.size()) return vector<T> ();
	vector <T> rhs(size(),0);
	Matrix<T> a=(*this);
	for (int i=0; i<a.size(); ++i)
		for (int j=0; j<a.size(); ++j)
			rhs[i]=rhs[i]+a.matrix[i][j]*answer[j];
	return rhs;
}

template<class T>
Matrix<T> Matrix<T>::operator^(int a)
{
	Matrix<T> matr=(*this);
	int n=matr.size();
	if (a!=-1)
	{
		Matrix<T> temp;
		if (a==0)
		{
			temp==unitMatr<T>(n);
			return temp;
		}
		temp=matr;
		for (int i=2; i<=a; ++i)
			temp=temp*matr;
		return temp;
	}
	Matrix<T> L;
	vector<Matrix<T>> LL;
	for (int i=0; i<n; ++i)
	{
		L.matrix.clear();
		L=unitMatr<T>(n);
		for (int j=0; j<n; ++j)
			if (j==i) L[j][i]=T(1)/matr[i][i];
			else L[j][i]=(matr[j][i]/matr[i][i])*(-1);
		matr=L*matr;
		LL.push_back(L);
	}
	Matrix<T> A;
	A=unitMatr<T>(n);
	for (int i=LL.size()-1; i>=0; --i)
		A=A*LL[i];
	return A;
}

template<class T>
vector<T> Matrix<T>::Givens(vector<T> right)
{
	int n=right.size();
	if(n!=size() || n==0) return vector<T>();
	Matrix<T> a;
	a=(*this);
	for (int i=0; i<n; ++i)
		a[i].push_back(right[i]);
	T c, s;
	for (int j=0; j<n; ++j)
	{
		for (int i=j+1; j<n; ++j)
		{
			c=T(0);s=T(0);
			c=a[j][j]/sqrt(a[j][j]*a[j][j] + a[i][j] * a[i][j]);
			s=a[i][j]/sqrt(a[j][j]*a[j][j] + a[i][j] * a[i][j]);
			for (int k=0; k<n+1; ++k)
			{
				a[j][k] = a[j][k]*c + a[i][k]*s;
				a[i][k] = a[i][k]*c - a[j][k]*s;
			}
		}
	}
	vector<T> x(n,T(0));
	T temp=0;
	x[n-1]=a[n-1][n]/a[n-1][n-1];
	for(int i=n-2; i>=0; --i)
	{
		temp=a[i][n];
		for(int j=n-1; j>i; --j)
			temp = temp - a[i][j]*x[j];
		x[i]=temp/a[i][i];
	}
	return х;
}

template<class T>
vector<T> Matrix<T>::HouseHolder(vector<T> right)
{
	int n = size();
	if (n!=right.size() || right.size()==0) return vector<T> ();
	Matrix<T> A(*this);

	vector<T> u, a, w;
	T s, y;
	Matrix<T> A1(*this);
	Matrix<T> I;
	I=unitMatr<T>(n);
	vector<Matrix<T>> PV;
	Matrix<T> P(n);
	Matrix<T> R(n);

	for (int i=0; i<n; ++i)
	{
		a.clear();
		for (int k=0; k<n; ++k)
			a.push_back(A1[k][i]);
		u=a;
		s=sqrt(a*a);
		u[i]=u[i]-s;
		y=T(1)/(s*s-A1[i][i]*s);
		y=sqrt(y);
		w=u*y;
		Matrix<T> W(n);
		for (int j=0; j<n; ++j)
			for (int k=0; k<n; ++k)
				W[j][k]=w[j]*w[k];//множення вектор-стовпчика на вектор-рядок
		P=I-W;
		PV.push_back(P);
		A1=P*A1;
	}

	vector<T> R1;
	R=PV[n-1];
	for (int i=n-2; i>=0; --i)
		R=R*PV[i];
	R1=R*right;
	R=R*A;

	for (int i=0; i<n; ++i)
		R[i].push_back(R1[i]);
	
	vector<T> x(n,0);
	T temp(0);
	x[0]=R[0][n]/R[0][0];
	for (int i=1; i<n; ++i)
	{
		temp=R[i][n];
		for (int j=0; j<i; ++j)
			temp=temp-R[i][j]*x[j];
		x[i]=temp/R[i][i];
	}
	return х;
}

template<class T>
vector<T> Matrix<T>::GramShmidt(vector<T> right)
{
	int n=size();
	if (n!=right.size() || right.size()==0) return vector<T> ();
	Matrix<T> A(*this);

	Matrix<T> Q(n, T(0)), R(n, T(0));
	T temp;

	for (int i=0; i<n; ++i)
	{
		Q[i]=A[i];
		for (int j=0; j<i; ++j)
		{
			R[i][j]=Q[j]*Q[i];
			Q[i]=Q[i]-(Q[j]*R[i][j]);
		}
		R[i][i]=T(0);
		for (int k=0; k<n; ++k)
			R[i][i]=R[i][i]+Q[i]*Q[i];
		R[i][i]=sqrt(R[i][i]);
		if (R[i][i]==T(0)) return vector<T> ();
		temp=(T(1))/R[i][i];
		Q[i]=Q[i]*temp;
	}
	vector<T> x;
	x=((R^(-1))*Q.transpose())*right;
	return х;
}

template<class T>
vector<T> Matrix<T>::Holeckyj(vector<T> right)
{
	int n=size();
	if (n!=right.size() || n==0) return vector<T> ();
	Matrix<T> A(*this);
	Matrix<T> L(n,T(0));
	T sum;

	for (int i=0; i<n; ++i)
	{
		for (int j=0; j<i; ++j)
		{
			sum=T(0);
			for (int k=0; k<j; ++k)
				sum = sum + L[i][k]*L[j][k];
			sum = A[i][j] - sum;
			L[i][j]=sum/L[j][j];
		}
		L[i][i]=T(0);
		for (int k=0; k<i; ++k)
			L[i][i]=L[i][i]+(L[i][k]*L[i][k]);
		L[i][i]=A[i][i]-L[i][i];
		L[i][i]=sqrt(L[i][i]);
	}

	vector<T> y(n), x(n);
	T temp(0);

	y[0]=right[0]/L[0][0];
	for (int i=1; i<n; ++i)
	{
		temp=right[i];
		for (int j=0; j<i; ++j)
			temp=temp-L[i][j]*y[j];
		y[i]=temp/L[i][i];
	}

	L=L.transpose();
	temp=T(0);
	x[n-1]=y[n-1]/L[n-1][n-1];
	for (int i=n-2; i>=0; --i)
	{
		temp=y[i];
		for (int j=n-1; j>i; --j)
			temp=temp-L[i][j]*x[j];
		x[i]=temp/L[i][i];
	}
	return x;
}

template<class T>
T MaxOutDiagonalElem(Matrix<T> a, int & p, int  & q)
{
	T MAX(0);
	for (int i=0; i<a.size(); ++i)
		for (int j=0; j<a.size(); ++j)
			if (i!=j && abs(a[i][j])>MAX) 
			{
				MAX=a[i][j];
				p=i;
				q=j;
			}
	return MAX;
}

template<class T>
T quadEq(T b)
{
	T D;
	D=b*b+T(4);
	if (D<T(0)) return T(-999999);
	T t1, t2;
	t1=(b*(-1)-sqrt(D))/2;
	t2=(b*(-1)+sqrt(D))/2;
	return abs(t1)<abs(t2) ? t1 : t2;
}

template<class T>
Matrix<T> Matrix<T>::Jacobi()
{
	int n=size();
	Matrix<T> A(*this);
	Matrix<T> U(n);
	T MAX(0);
	int p=0, q=0;
	T teta=0;
	T phi, u1, u2;

	for (int p=0; p<n; ++p)
	{
		for (int q=0; q<n; ++q)
		{
			//не шукаємо максимальний позадіагональний елемент, бо так зменшується точність
			//просто проробляємо для всіх елементів
			/*MAX=T(0); p=0; q=0;
			MAX=MaxOutDiagonalElem(A,p,q);*/
			if (p!=q && A[p][q]!=T(0))
			{
				teta=(A[q][q]-A[p][p])/(A[p][q]+A[p][q]);
				teta=quadEq(T(2)*teta);
				//teta=T(1)/(teta+teta);//тут втрачається точність
				phi=atan(teta);
				u1=cos(phi);
				u2=sin(phi);
				for (int i=0; i<n; ++i)
					for (int j=0; j<n; ++j)
						if (i==j) U[i][j]=T(1);
						else U[i][j]=T(0);
				U[p][p]=u1;
				U[q][q]=u1;
				U[p][q]=u2;
				U[q][p]=u2*(-1);
				A=U.transpose()*A*U;
			}
		}
	}
	return A;
}

template<class T>
vector<T> Matrix<T>::Grad(vector<T> right, vector <T> u, T eps)
{
	if(right.size()!=size() || right.size()==0) return vector<T>();
	vector<T> d,x,x1(right),l;
	T s=eps;
	T temp;
	u.clear();

	for(int i=0;i<(int)right.size();++i)
	{
		d.push_back(T(0));
		u.push_back(T(0));
	}

	Matrix<T> A(*this);
	for(int i=0; eps>s; ++i)
	{
		x = A*u - right;	

		d= d*((x*x)/(x1*x1)) - x;

		s = (d*x) / (d*(A*d));

		u = u + d*s;

		x1=x;
	}

	return х;
}

Matrix<Rational> Generate(int n)
{
	Matrix<Rational> result(n,Rational(1));
	result[0][0]=Rational(1);
	LI temp=0;
	for (int i=0; i<10; ++i)
	{
		if (i>0) result[i][0]=result[i-1][1];
		for (int j=1; j<10; ++j)
		{
			temp=result[i][j-1].m;
			result[i][j].m=temp+1;
		}
	}
	return result;
}

#endif