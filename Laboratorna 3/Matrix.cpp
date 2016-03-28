#include "Matrix.h"

#define EPS (Rational)1/10
#define INF vector<Rational> (80, (Rational) 1)

Matrix::Matrix()
{
	matrix.clear();
}

Matrix::Matrix(int _size, Rational val)
{
	if (_size==0) matrix.clear();
	vector <Rational> y(0);
	for (int i=0; i<_size; ++i)
	{
		matrix.push_back(y);
		for (int j=0; j<_size; ++j)
			matrix[i].push_back(val);
	}
}

bool operator>>(istream &i, Matrix &a)
{
	vector<Rational> temp;
	Rational t;
	while(i>>t) temp.push_back(t);
	int size=(int)temp.size();
	size=(int)sqrt((double)size);
	vector<Rational> y(0);
	int k=0;
	for (int i=0; i<size; ++i)
	{
		a.matrix.push_back(y);
		for (int j=0; j<size; ++j)
			a.matrix[i].push_back(temp[k++]);
	}
	return true;
}

bool operator<<(ostream &o, Matrix a)
{
	if (a.size()==0) return false;
	for (int i=0; i<a.size(); ++i)
	{
		for (int j=0; j<a.size(); ++j)
		{
			o.width(10);
			if (!(o<<a.matrix[i][j]) || !(o<<' ')) return false;
		}
		if(!(o<<endl)) return false;
	}
	return true;
}

Matrix Matrix::operator=(Matrix b)
{
	if (b.size()==0) return Matrix();
	matrix=b.matrix;
	return (*this);
}

Matrix unitMatr(int _size)
{
	Matrix g;
	if (_size==0) g.matrix.clear();
	vector <Rational> y(0);
	for (int i=0; i<_size; ++i)
	{
		g.matrix.push_back(y);
		for (int j=0; j<_size; ++j)
			if (i==j) g.matrix[i].push_back((Rational)1);
			else g.matrix[i].push_back((Rational)0);
	}
	return g;
}

Matrix Matrix::operator+(Matrix a)
{
	if(size()==0) return a;
	if(a.size()==0) return (*this);
	if(size()!=a.size()) return Matrix();
	Matrix result(size());
	for(int i=0; i<size(); ++i)
		for(int j=0; j<size(); ++j)
			result.matrix[i][j]=matrix[i][j]+a.matrix[i][j];
	return result;
}

Matrix Matrix::operator-(Matrix a)
{
	if(size()==0) return a*(-1);
	if(a.size()==0) return (*this);
	Matrix result(size());
	for(int i=0; i<size(); ++i)
		for(int j=0; j<size(); ++j)
			result.matrix[i][j]=matrix[i][j]-a.matrix[i][j];
	return result;
}

Matrix Matrix::operator*(int a)
{
	if (a==0) return Matrix(size());
	Matrix result(size());
	for(int i=0; i<size(); ++i)
		for(int j=0; j<size(); ++j)
			result.matrix[i][j]=matrix[i][j]*a;
	return result;
}

vector<Rational> Matrix::operator*(vector<Rational> a)
{
	
	if (a.size()==0) return vector <Rational> ();
	if (size()==0||size()!=int(a.size())) return vector<Rational> ();
	vector <Rational> result (size(),(Rational)0);
	for (int i=0; i<size(); ++i)
	{
		Rational temp=0;
		for (int j=0; j<size(); ++j)
		{
			temp=matrix[i][j]*a[j];
			result[i]=result[i]+temp;
			temp=(Rational)0;
		}
	}
	return result;
}

Matrix Matrix::operator*(Matrix a)
{
	int n=size();
	if (a.size()==0) return Matrix ();
	if (n==0||n!=int(a.size())) return Matrix ();
	Matrix b;
	b=(*this);
	Matrix result(n);
	for (int i=0; i<n; ++i)
	{
		Rational temp=0;
		for (int j=0, k=0; k<n; ++j)
		{
			if (j==n){j=0;++k;}
			if (k==n) break;
			temp=b.matrix[i][j]*a.matrix[j][k];
			result.matrix[i][k]=result.matrix[i][k]+temp;
			temp=(Rational)0;
		}
		temp=(Rational)0;
	}
	return result;
}

vector<Rational>& Matrix::operator[](int i)
{
	return matrix[i];
}

Matrix Matrix::transpose()
{
	if (size()==0) return Matrix();
	if (size()==1) return (*this);
	Matrix result(size());
	for (int i=0; i<size(); ++i)
		for (int j=0; j<size(); ++j)
			result.matrix[i][j]=(*this)[j][i];
	return result;
}

bool Matrix::operator==(Matrix a)
{
	if (size()!=a.size()) return false;
	for (int i=0; i<size(); ++i)
		for (int j=0; j<size(); ++j)
			if (matrix[i][j]!=a.matrix[i][j]) return false;
	return true;
}

bool Matrix::operator!=(Matrix a)
{
	if (size()!=a.size()) return true;
	for (int i=0; i<size(); ++i)
		for (int j=0; j<size(); ++j)
			if (matrix[i][j]!=a.matrix[i][j]) return true;
	return false;
}

Matrix Matrix::Vinograd(Matrix b)
{
	if (size()!=b.size() || size()==0 || b.size()==0) return Matrix();
	int sz=size();
	Matrix a=(*this);
	bool is_zeros=false;
	if (sz%2)
	{
		is_zeros=true;
		++sz;
		vector <Rational> temp (sz,Rational(0));
		b.matrix.push_back(temp);
		a.matrix.push_back(temp);
		for (int i=0; i<sz-1; ++i)
		{
			a.matrix[i].push_back(Rational(0));
			b.matrix[i].push_back(Rational(0));
		}
	}
	Rational *A=new Rational();
	Rational *B=new Rational();
	Rational *C=new Rational();
	Rational *temp=new Rational();
	Matrix result;
	vector<Rational> y(0);
	for (int i=0; i<sz; ++i)
	{
		(*A)=(Rational)0;
		result.matrix.push_back(y);
		for (int k=0; k<sz/2; ++k)
		{
			(*temp)=a.matrix[i][k+k]*a.matrix[i][k+k+1];
			(*A)=(*A)+(*temp);
		}
		for (int j=0; j<sz; ++j)
		{	
			(*B)=(Rational)0;
			(*C)=(Rational)0;
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

Matrix Matrix::Shtrassen(Matrix b)
{
	if (size()!=b.size() || size()==0 || b.size()==0) return Matrix();
	int size_=size();
	Matrix a=(*this);
	bool is_zeros=false;
	if (size_%2)
	{
		is_zeros=true;
		++size_;
		vector <Rational> temp(size_,Rational(0));
		b.matrix.push_back(temp);
		a.matrix.push_back(temp);
		for (int i=0; i<size_; ++i)
		{
			a.matrix[i].push_back(Rational(0));
			b.matrix[i].push_back(Rational(0));
		}
	}
	Matrix* A11; Matrix *A12; Matrix *A21; Matrix *A22;
	Matrix* B11; Matrix *B12; Matrix *B21; Matrix *B22;
	A11=new Matrix();
	A12=new Matrix();
	A21=new Matrix();
	A22=new Matrix();
	B11=new Matrix();
	B12=new Matrix();
	B21=new Matrix();
	B22=new Matrix();
	vector<Rational> temp;
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
	Matrix *P1; Matrix *P2; Matrix *P3; Matrix *P4;
	Matrix *P5; Matrix *P6; Matrix *P7;
	P1=new Matrix(size_/2);
	P2=new Matrix(size_/2);
	P3=new Matrix(size_/2);
	P4=new Matrix(size_/2);
	P5=new Matrix(size_/2);
	P6=new Matrix(size_/2);
	P7=new Matrix(size_/2);
	(*P1)=((*A11)+(*A22)).Vinograd((*B11)+(*B22));
	(*P2)=((*A21)+(*A22)).Vinograd(*B11);
	(*P3)=(*A11).Vinograd((*B12)-(*B22));
	(*P4)=(*A22).Vinograd((*B21)-(*B11));
	(*P5)=((*A11)+(*A12)).Vinograd(*B22);//тут хуйня, А11.сайз=2
	(*P6)=((*A21)-(*A11)).Vinograd((*B11)+(*B12));
	(*P7)=((*A12)-(*A22)).Vinograd((*B21)+(*B22));
	Matrix* C11; Matrix *C12; Matrix *C21; Matrix *C22;
	C11=new Matrix(size_/2);
	C12=new Matrix(size_/2);
	C21=new Matrix(size_/2);
	C22=new Matrix(size_/2);
	(*C11)=(*P1)+(*P4)-(*P5)+(*P7);
	(*C12)=(*P3)+(*P5);
	(*C21)=(*P2)+(*P4);
	(*C22)=(*P1)+(*P3)-(*P2)+(*P6);
	Matrix result;
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

vector<Rational> Matrix::Gauss(vector<Rational> a)
{
	if (size()!=a.size() || a.size()==0) return vector<Rational> ();
	Matrix matr=(*this);
	int n = size();
	for (int i=0; i<size(); ++i)
		matr[i].push_back(a[i]);
	for (int col=0; col<n; ++col)
	{
		int row=0;
		for (row=col; row<n && matr[row][col]==(Rational)0; ++row);
		if (row==n) return vector<Rational> ();
		swap(matr.matrix[row],matr.matrix[col]);

		Rational x=matr[col][col];
		for (int j=0; j<n+1; ++j)
			matr[col][j]=matr[col][j]/x;

		vector<Rational> temp;
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

	vector<Rational> ans(n,0);
	Rational temp(0);
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

vector<Rational> Matrix::GenRight(vector<Rational> answer)
{
	if (size()!=answer.size()) return INF;
	vector <Rational> rhs(size(),0);
	Matrix a=(*this);
	for (int i=0; i<a.size(); ++i)
		for (int j=0; j<a.size(); ++j)
			rhs[i]=rhs[i]+a.matrix[i][j]*answer[j];
	return rhs;
}

Matrix Matrix::operator^(int a)
{
	Matrix matr=(*this);
	int n=matr.size();
	if (a!=-1)
	{
		Matrix temp;
		if (a==0)
		{
			temp==unitMatr(n);
			return temp;
		}
		temp=matr;
		for (int i=2; i<=a; ++i)
			temp=temp*matr;
		return temp;
	}
	Matrix L;
	vector<Matrix> LL;
	for (int i=0; i<n; ++i)
	{
		L.matrix.clear();
		L=unitMatr(n);
		for (int j=0; j<n; ++j)
			if (j==i) L[j][i]=Rational(1)/matr[i][i];
			else L[j][i]=(matr[j][i]/matr[i][i])*(-1);
		matr=L*matr;
		LL.push_back(L);		
	}
	Matrix A;
	A=unitMatr(n);
	for (int i=LL.size()-1; i>=0; --i)
		A=A*LL[i];
	return A;
}
/*
vector<Rational> Matrix::Grad(vector<Rational> right)
{
	vector<Rational> u(size());

}*/

vector<Rational> Matrix::Givens(vector<Rational> right)
{
	int n=right.size();
	if(n!=size() || n==0) return vector<Rational>();
	Matrix a;
	a=(*this);
	for (int i=0; i<n; ++i)
		a[i].push_back(right[i]);
	Rational c, s;
	for (int j=0; j<n; ++j)
	{
		for (int i=j+1; j<n; ++j)
		{
			c=Rational(0);s=Rational(0);
			c=a[j][j]/sqrt(a[j][j]*a[j][j] + a[i][j] * a[i][j]);
			s=a[i][j]/sqrt(a[j][j]*a[j][j] + a[i][j] * a[i][j]);
			for (int k=0; k<n+1; ++k)
			{
				a[j][k] = a[j][k]*c + a[i][k]*s;
				a[i][k] = a[i][k]*c - a[j][k]*s;
			}
		}
	}
	vector<Rational> x(n,Rational(0));
	Rational temp=0;
	x[n-1]=a[n-1][n]/a[n-1][n-1];
	for(int i=n-2; i>=0; --i)
	{
		temp=a[i][n];
		for(int j=n-1; j>i; --j)
			temp = temp - a[i][j]*x[j];
		x[i]=temp/a[i][i];
	}
	return x;

}

vector<Rational> Matrix::HousHolder(vector<Rational> right)
{
	int n = size();
	if (n!=right.size() || right.size()==0) return vector<Rational> ();
	Matrix A(*this);

	vector<Rational> u, a, w;
	Rational s, y;
	Matrix A1(*this);
	Matrix I;
	I=unitMatr(n);
	vector<Matrix> PV;
	Matrix P(n);
	Matrix R(n);

	for (int i=0; i<n; ++i)
	{
		a.clear();
		for (int k=0; k<n; ++k)
			a.push_back(A1[k][i]);
		u=a;
		s=sqrt(a*a);
		u[i]=u[i]-s;
		y=Rational(1)/(s*s-A1[i][i]*s);
		y=sqrt(y);
		w=u*y;
		Matrix W(n);
		for (int j=0; j<n; ++j)
			for (int k=0; k<n; ++k)
				W[j][k]=w[j]*w[k];//множення вектор-стовпчика на вектор-рядок
		P=I-W;
		PV.push_back(P);
		A1=P*A1;
	}

	vector<Rational> R1;
	R=PV[n-1];
	for (int i=n-2; i>=0; --i)
		R=R*PV[i];
	R1=R*right;
	R=R*A;

	for (int i=0; i<n; ++i)
		R[i].push_back(R1[i]);
	
	vector<Rational> answer(n,0);
	Rational temp(0);
	answer[n-1]=R[n-1][n];
	for (int i=n-2; i>=0; --i)
	{
		temp=R[i][n];
		for (int j=n-1; j>i; --j)
			temp=temp-R[i][j]*answer[j];
		answer[i]=temp/R[i][i];
	}
	return answer;
}

vector<Rational> operator*(vector<Rational> a, Rational b)
{
	for (int i=0; i<(int)a.size(); ++i)
		a[i]=a[i]*b;
	return a;
}

vector<Rational> operator+(vector<Rational> a, vector<Rational> b)
{
	vector<Rational> result;
	if (a.size()!=b.size()) return result;
	for (int i=0; i<(int)a.size(); ++i)
		result.push_back(a[i]+b[i]);
	return result;
}

vector<Rational> operator-(vector<Rational> a, vector<Rational> b)
{
	vector<Rational> result;
	if (a.size()!=b.size()) return result;
	for (int i=0; i<(int)a.size(); ++i)
		result.push_back(a[i]-b[i]);
	return result;
}

Rational operator*(vector<Rational>a, vector<Rational> b)
{
	Rational result(0);
	if (a.size()!=b.size() || a.size()==0) return result;
	for (int i=0; i<(int)a.size(); ++i)
		result=result+a[i]*b[i];
	return result;
}

vector<Rational> Matrix::GramShmidt(vector<Rational> right)
{
	int n=size();
	if (n!=right.size() || right.size()==0) return vector<Rational> ();
	Matrix A(*this);

	Matrix Q(n, Rational(0)), R(n, Rational (0));
	Rational temp;

	for (int i=0; i<n; ++i)
	{
		Q[i]=A[i];
		for (int j=0; j<i; ++j)
		{
			R[i][j]=Q[j]*Q[i];
			Q[i]=Q[i]-(Q[j]*R[i][j]);
		}
		R[i][i]=Rational(0);
		for (int k=0; k<n; ++k)
			R[i][i]=R[i][i]+Q[i]*Q[i];
		R[i][i]=sqrt(R[i][i]);
		if (R[i][i]==Rational(0)) return vector<Rational> ();
		temp=(Rational(1))/R[i][i];
		Q[i]=Q[i]*temp;
	}
	vector<Rational> answer;
	answer=((R^(-1))*Q.transpose())*right;
	return answer;
}

vector<Rational> Matrix::Holeckyj(vector<Rational> right)
{
	int n=size();
	if (n!=right.size() || n==0) return vector<Rational> ();
	Matrix A(*this);
	Matrix L(n,Rational(0));
	Rational sum;

	for (int i=0; i<n; ++i)
	{
		for (int j=0; j<i; ++j)
		{
			sum=Rational(0);
			for (int k=0; k<j; ++k)
				sum = sum + L[i][k]*L[j][k];
			sum = A[i][j] - sum;
		}
		L[i][i]=Rational(0);
		for (int k=0; k<i; ++k)
			L[i][i]=L[i][i]+(L[i][k]*L[i][k]);
		L[i][i]=A[i][i]-L[i][i];
		L[i][i]=sqrt(L[i][i]);
	}

	vector<Rational> y(n), x(n);
	Rational temp(0);

	y[0]=right[0]/L[0][0];
	for (int i=1; i<n; ++i)
	{
		temp=right[i];
		for (int j=0; j<i; ++j)
			temp=temp-L[i][j]*y[j];
		y[i]=temp/L[i][i];
	}

	L=L.transpose();
	temp=Rational(0);
	x[n-1]=y[n-1]/L[n-1][n-1];
	for (int i=n-2; i>=0; --i)
	{
		temp=y[i];
		for (int j=n-1; j>i; ++j)
			temp=temp-L[i][j]*x[j];
		x[i]=temp/L[i][i];
	}
	return x;
}

Rational MaxOutDiagonalElem(Matrix a, int & p, int  & q)
{
	Rational MAX(0);
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

Matrix Matrix::Jacobi()
{
	int n=size();
	Matrix A(*this);
	Matrix U(n);
	Rational MAX(0);
	int p=0, q=0;
	Rational teta=0;
	Rational phi, u1, u2;

	for (int i=0; i<n; ++i)
	{
		MAX=Rational(0); p=0; q=0;
		MAX=MaxOutDiagonalElem(A,p,q);
		teta=(A[q][q]-A[p][p])/(A[p][q]+A[p][q]);
		teta=Rational(1)/(teta+teta);
		phi=atan(teta);
		u1=cos(phi);
		u2=sin(phi);
		for (int i=0; i<n; ++i)
			for (int j=0; j<n; ++j)
				if (i==j) U[i][j]=Rational(1);
				else U[i][j]=Rational(0);
		U[p][p]=u1;
		U[q][q]=u1;
		U[p][q]=u2;
		U[q][p]=u2*(-1);
		A=U.transpose()*A*U;
	}
	return A;
}

vector<Rational> Matrix::Grad(vector<Rational> right, vector <Rational> x, Rational eps)
{
	if(right.size()!=size() || right.size()==0) return vector<Rational>();
	vector<Rational> d,g,g1(right),l;
	Rational s=eps;
	Rational temp;
	x.clear();

	for(int i=0;i<(int)right.size();++i)
	{
		d.push_back(Rational(0));
		x.push_back(Rational(0));
	}

	Matrix A(*this);
	for(int i=0; eps>s; ++i)
	{
		g = A*x - right;	

		d= d*((g*g)/(g1*g1)) - g;

		s = (d*g) / (d*(A*d));

		x = x + d*s;

		g1=g;
	}

	return g;
}