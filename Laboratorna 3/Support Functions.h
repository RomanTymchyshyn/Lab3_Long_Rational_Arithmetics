#pragma once

#include <vector>
using namespace std;
vector<int> divideVectorByVector(vector<int> a, vector<int> b, const int &base);
vector<int> moduloVector(vector<int> a, vector<int> b, const int &base);
vector<int> yakobi(vector<int> a, vector< int> P, const int &base);
vector<int> gcd(vector<int> a, vector<int> b, const int &base);
vector<int> subtractVectors(vector<int> a,vector<int> b, const int &base);
vector<int> IpowVector(const vector<int> &a, vector<int> n, const vector<int> &roll, const int &base );
vector<int> convertFromBase(vector<int> a, const int &base);
vector<int> convertInBase(vector<int> a, const int &base);
bool soloveiShtrassen(vector<int> m, int k, int base);