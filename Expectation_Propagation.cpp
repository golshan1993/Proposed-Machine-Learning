//#ifdef _MSC_VER
//#define _CRT_SECURE_NO_WARNINGS
//#endif
//#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <conio.h>
#include <fstream>
//#include <exception>
using namespace std;

vector <int> random(int, float);
vector<int> non_erased(vector<int>);
vector <int> forced(int*, vector <int>, int, int);
int* forced_generetor(int*, vector <int>, vector <int>, int);
int max(vector<int>);
vector <int> multiply(vector<int>, int*, int, int);
vector <int> check_degree(int*, int, int);
vector <int> variable_degree(int*, int, int);
int * Generator(int*, int*, int, int, int);
int * Alist_A(int*, int*, int, int, int);
int * Alist_B(int*, int*, int, int, int); 
int check_chosen(int* , vector <int>, vector <int>, int);
int check_chosen1(int*, vector <int>, vector <int>, int);
int bit_chosen(int*, vector<int>, vector<int>, vector<int>, int, int);
vector<int> check_changes(int*, vector<int>, vector<int>, vector<int>, int, int);
int * Decimation(int*, int*, int*, vector <int>, vector <int>, vector <int>, int, int, int, int, int, int);
int * DecodeMatrix(int*, int*, vector<int>, int, int, int);
vector<int> Parentvariables(int*, vector<int>, vector<int>, int, int);
vector<int> Quantize(int*, int*, vector<int>, vector<int>, vector<int>, vector<int>, int, int, int);
int main()
{
	ofstream myfile;
    myfile.open ("resultGTEP(5,10).txt");
    int n = 1024;
    int k = 512;
    int dc = 10;
    int dv = 5;
	int *A = (int *)malloc(k * dc * sizeof(int));
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < dc; j++) {
			A[i*dc + j] = 0;
		}
	}
	int *B = (int *)malloc(n * dv * sizeof(int));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < dv; j++) {
			B[i*dv + j] = 0;
		}
	}
	int *G = (int *)malloc(k * n * sizeof(int));
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < n; j++) {
			G[i*n + j] = 0;
		}
	}
	FILE * fp1;
	FILE * fp2;
	fp1 = fopen("Alist_A_1024(5,10).txt", "r");
	fp2 = fopen("Alist_B_1024(5,10).txt", "r");
	int c1, c2, a1, a2, b1, b2;
	int cnt = 0;
	while (1)
	{

		fscanf(fp1, "%d", &c1);
		fscanf(fp2, "%d", &c2);

		cnt = cnt + 1;
		a1 = (cnt - 1) / dc;
		b1 = (cnt - 1) % dc;
		a2 = (cnt - 1) / dv;
		b2 = (cnt - 1) % dv;
		B[a2*dv + b2] = c2;
		A[a1*dc + b1] = c1;

		//if (feof(fp1)) {
		//	break;
		//}
		if (cnt == k*dc) {
			break;
		}

	}
	
	fclose(fp1);
	fclose(fp2);
	G = Generator(A, G, k, n, dc);
	double eps;
	int K;
	//vector<double> D;
	for (int i0 = 0; i0 <11; i0++) {
		eps = 0.3+0.02*(i0);
	
		myfile << "\n"<<eps;
		int t1 = 0;
		
		//vector <double> d;
		double S = 0;
		int k1 = 0;
		while (t1 != 1 ) {
			//k1 = 0;
			vector <int> y = random(n, eps);
			vector <int> y1 = non_erased(y);
			vector <int> v1 = forced(A, y, k, dc);
			int a = v1.size();
			int b = y1.size();
			int*G1 = forced_generetor(G, y1, v1, n);
			int*G2 = (int *)malloc(a * b * sizeof(int));
			cnt = 0;
            K = 0;
			//printf("%d\n", b);
			int cn = 0;
			int cn1 = 0;
			vector <int> w;
			vector <int> w1;
			vector <int> W;
			for (int i = 0; i < k; i++) {
				W.push_back(2);
			}
			for (int i = 0; i < a; i++) {
				w.push_back(2);
				w1.push_back(0);
			}
			while (cn != b) {
				cn1 = 0;
				//printf("%d\n", cn);
				std::vector<int> u2 = variable_degree(G1, b, a);
				std::vector<int> u3 = check_degree(G1, b, a);
				int m1 = max(u3);
				int m2 = max(u2);
				int*A1 = (int *)malloc(a * m2 * sizeof(int));
				int*B1 = (int *)malloc(b * m1 * sizeof(int));
				A1 = Alist_A(A1, G1, a, b, m2);
				B1 = Alist_B(B1, G1, a, b, m1);
				cnt = 0;

				int M = check_chosen(B1, u3, w1, m1);
				//printf("%d\n", u3[M]);
			    
				if (u3[M] == 1) {
					W[v1[B1[M*m1 + 0] - 1]] = y[y1[M]];
					w[B1[M*m1 + 0] - 1] = y[y1[M]];
					K++;
				}
				else if (u3[M] > 1) {
				    if (u3[M] > 2){
				    	//printf("%d\n", u3[M]);
				    	//break;
				    	k1++;
					}
					//int M = check_chosen1(B1, u3, w1, m1);
					G2 = DecodeMatrix(G2, B1, u3, M, b, m1);
					w1 = Parentvariables(B1, u3, w1, M, m1);
				}

				a1 = bit_chosen(B1, u3, u2, w1, m1, M);
				if (y[y1[M]] == 1) {
					y = check_changes(A1, y, y1, u2, a1, m2);

				}
				G1 = Decimation(G1, A1, B1, w1, u3, u2, M, a, a1, b, m1, m2);

				cn = 0;
				for (int i = 0; i < b; i++) {
					cnt = 0;
					for (int j = 0; j < a; j++) {
						if (G1[j*b + i] == 1) {
							cnt++;
						}
					}
					if (cnt == 0) {
						cn++;
					}
				}
				for (int i = 0; i < a; i++) {
					if (w1[i] == 2) {
						cn1++;
					}
				}
				//printf("%d\n", cn1);
				//_getch();
				delete[]A1;
				delete[]B1;

			}
			if (K != 0) {
			
				t1++;
				//printf("%d\n",t1);
				vector<int> u4 = check_degree(G2, b, a);
				int m3 = max(u4);

				int*B2 = (int *)malloc(b * m3 * sizeof(int));
				B2 = Alist_B(B2, G2, a, b, m3);
				W = Quantize(G2, B2, u4, v1, w, W, a, b, m3);
				vector<int> Y = multiply(W, G, k, n);

				double s = 0;

				for (int i = 0; i < y1.size(); i++) {
					a = (y[y1[i]] + Y[y1[i]]) % 2;


					if (a == 1) {
						s++;
					}
				}
				double d1 = s / n;
				//printf("%E\n", d1);
				S = S + d1;
			}
		}
			//if (k1 != 0){
				//	myfile << " GTEP";
				    //break;
				//}
		double D = S / (t1);
		//printf("%E\n", D);
		myfile << " "<<D;
	  }
	myfile.close();	
	}

//
// generate random vector:
vector <int> random(int n, float eps)
{
	int e = 2;
	vector <int> y;
	for (int i = 0; i < n; i++) {
		double a = ((double)rand() / (RAND_MAX));
		if (a>(1 - eps)) {
			y.push_back(e);

		}
		else if (a>(1 - eps) / 2 && a<(1 - eps)) {

			y.push_back(0);
		}
		else if (a>0 && a<(1 - eps) / 2) {

			y.push_back(1);
		}

	}
	return y;
}
// non_erased arguments:
vector <int> non_erased(vector<int> y)
{
	int e = 2;
	vector<int> y1;
	for (int i = 0; i < y.size(); i++) {
		if (y[i] != e) {
			y1.push_back(i);
		}
	}
	return y1;
}
// forced variables:
vector <int> forced(int*A, vector<int> y, int k, int dc)
{
	int e = 2;
	vector<int> v1;
	for (int i = 0; i < k; i++) {
		int cnt = 0;
		for (int j = 0; j < dc; j++) {
			if (y[A[i*dc + j] - 1] != e) {
				cnt++;
			}
		}
		if (cnt != 0) {
			v1.push_back(i);
		}
	}
	return v1;
}
// forced generator matrix:
int*forced_generetor(int*G, vector <int> y1, vector <int> v1, int n)
{
	int a = v1.size();
	int b = y1.size();
	int* G1 = (int *)malloc(a * b * sizeof(int));
	
	for (int i = 0; i < a; i++) {

		for (int j = 0; j < b; j++) {
			//G1[i*b + j] = 0;
			if (G[v1[i] * n + y1[j]] == 1) {
				G1[i*b + j] = 1;
			}
			else if (G[v1[i] * n + y1[j]] == 0) {
				G1[i*b + j] = 0;
			}
		}
	}


	return(G1);
}
// finding degrees of check nodes and variable nodes:
vector <int> check_degree(int*G1, int b, int a)
{
	vector <int>  u3;
	for (int i = 0; i < b; i++) {
		int cnt = 0;
		for (int j = 0; j < a; j++) {
			if (G1[j*b + i] == 1) {
				cnt++;
			}

		}
		u3.push_back(cnt);
	}
	return(u3);
}
vector <int> variable_degree(int*G1, int b, int a)
{
	vector <int>  u1;
	for (int i = 0; i < a; i++) {
		int cnt = 0;
		for (int j = 0; j < b; j++) {
			if (G1[i*b + j] == 1) {
				cnt++;
			}

		}
		u1.push_back(cnt);
	}
	return(u1);
}

int * Generator(int *A, int *G,int k,int n,int dc)
{
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < dc; j++) {
			if (A[i*dc + j] != 0) {
				G[i*n + A[i*dc + j] - 1] = 1;
			}
		}
	}
	return(G);
}

// finding maximum of a vector:
int max(vector <int> u1)
{
	int max = u1[0];
	for (int i = 0; i < u1.size(); i++) {
		if (u1[i] > max) {
			max = u1[i];
		}
	}
	return(max);
}
vector<int> multiply(vector<int> W, int*G, int k, int n)
{
	vector<int> Y;
	for (int i = 0; i < n; i++) {
		Y.push_back(0);
		for (int j = 0; j < k; j++) {
			Y[i] = Y[i] + W[j] * G[j*n+i];
		}
	}
	return(Y);
}
// Alist Matrix of forced generator:
int *Alist_A(int*A1, int*G1, int a, int b, int m2)
{
	A1 = (int *)malloc(a * m2 * sizeof(int));
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < m2; j++) {
			A1[i*m2 + j] = 0;
		}
	}
	//A1 = { 0 };
	for (int i = 0; i < a; i++) {

		int cnt0 = 0;
		for (int j = 0; j < b; j++) {
			if (G1[i*b + j] == 1) {
				A1[i*m2 + cnt0] = j + 1;
				cnt0++;
			}
		}
	}
	return(A1);
}

int *Alist_B(int*B1, int*G1, int a, int b, int m1)
{
	B1 = (int *)malloc(b * m1 * sizeof(int));
	//B1 = { 0 };
	for (int i = 0; i < b; i++) {
		for (int j = 0; j < m1; j++) {
			B1[i*m1 + j] = 0;
		}
	}
	for (int i = 0; i < b; i++) {

		int cnt0 = 0;
		for (int j = 0; j < a; j++) {
			if (G1[j*b + i] == 1) {
				B1[i*m1 + cnt0] = j + 1;
				cnt0++;
			}
		}
	}
	return(B1);
}

int check_chosen(int *B1, vector<int> u3, vector<int> w1, int m1)
{
	vector <int> u4;
	for (int i = 0; i < u3.size(); i++) {
		if (u3[i] != 0) {
			u4.push_back(i);
		}
	}
	int min = u3[u4[0]];
	int M;
	for (int i = 0; i < u4.size(); i++) {
		if (min > u3[u4[i]] || min == u3[u4[i]]) {
			min = u3[u4[i]];
			M = u4[i];
		}
	}
	return(M);
}

int check_chosen1(int *B1, vector<int> u3, vector<int> w1, int m1)
{
	vector <int> u4;
	vector <int> u5;
	for (int i = 0; i < u3.size(); i++) {
		
		if (u3[i] != 0) {
			u4.push_back(i);
			
		}
	}

	vector<int> u6;
	for (int i = 0; i < u4.size(); i++) {
		int cn = 0;
		for (int j = 0; j < u3[u4[i]]; j++) {
			if (w1[B1[u4[i] * m1 + j] - 1] != 2) {
				cn++;
			}
		}
		if (cn != 0) {
			u5.push_back(cn);
			u6.push_back(u4[i]);
		}
	}

	//if (u5.size() != 0) {
		int min = u5[0];
		int M = u6[0];
		for (int i = 0; i < u6.size(); i++) {
			if (min > u5[i]) {
				min = u5[i];
				int M = u6[i];
			}
		}
	//}
	return(M);
}



int bit_chosen(int *B1, vector<int> u3, vector<int> u2, vector<int> w1, int m1, int M)
{
	int a1 = B1[M*m1 + 0] - 1;
	int max = u2[a1];
	for (int j = 0; j < u3[M]; j++) {

		if (w1[B1[M*m1 + j] - 1] != 2) {

			if (max < u2[B1[M*m1 + j] - 1] || max == u2[B1[M*m1 + j] - 1]) {
				max = u2[B1[M*m1 + j] - 1];
				a1 = B1[M*m1 + j] - 1;
			}

		}

	}
	return(a1);
}

vector<int> check_changes(int *A1, vector<int>y, vector<int>y1, vector<int>u2, int a1, int m2)
{
	for (int j = 0; j < u2[a1]; j++) {
		y[y1[A1[a1*m2+j] - 1]] = (1 + y[y1[A1[a1*m2+j] - 1]]) % 2;
	}
	return(y);
}

// Decimation process:
int*Decimation(int*G1, int*A1, int*B1, vector <int> w1, vector<int> u3, vector<int> u2, int M, int a, int a1, int b, int m1, int m2)
{
	
    for (int j = 0; j < u3[M]; j++) {

		for (int j0 = 0; j0 < u2[a1]; j0++) {
			//if (A1[a1*m2 + j0] != 0 && B1[M*m1 + j] != 0) {
				G1[(B1[M*m1 + j] - 1)*b + A1[a1*m2 + j0] - 1] = (1 + G1[(B1[M*m1 + j] - 1)*b + A1[a1*m2 + j0] - 1]) % 2;
			//}
		}


	}
	return(G1);
}

int * DecodeMatrix(int *G2, int *B1, vector<int>u3, int M, int b, int m1)
{
	for (int i = 0; i < u3[M]; i++) {
		G2[(B1[M*m1+i]-1)*b+M] = 1;
	}
	return(G2);
}

vector<int> Parentvariables(int *B1, vector<int>u3, vector<int>w1, int M, int m1)
{
	for (int i = 0; i < u3[M]; i++) {
		w1[B1[M*m1 + i] - 1] = 2;
   }
	return(w1);
}

vector<int> Quantize(int *G2, int *B2, vector<int>u4, vector<int> v1, vector<int>w, vector<int>W, int a, int b, int m3)
{
	vector<int> u;
	for (int i = 0; i < u4.size(); i++) {
		if (u4[i] != 0) {
			u.push_back(i);
		}
	}
	int cn1 = 1;
	while (cn1 != 0) {
		cn1 = 0;
		for (int i = 0; i < u.size(); i++) {
			int cn = 0;
			int a3;
			vector<int> v2;
			for (int j = 0; j < u4[u[i]]; j++) {
				if (w[B2[u[i] * m3 + j] - 1] == 2) {
					cn++;
					a3 = j;
				}
				else if (w[B2[u[i] * m3 + j] - 1] != 2) {
					v2.push_back(j);
				}
			}
			if (cn == 1 && v2.size() != 0) {
				cn1++;
				for (int j = 0; j < v2.size(); j++) {
					w[B2[u[i] * m3 + a3] - 1] = w[B2[u[i] * m3 + a3] - 1] + w[B2[u[i] * m3 + v2[j]] - 1];
				}
				w[B2[u[i] * m3 + a3] - 1] = w[B2[u[i] * m3 + a3] - 1] % 2;
				W[v1[B2[u[i] * m3 + a3] - 1]] = w[B2[u[i] * m3 + a3] - 1];
			}
		}
	}
	return(W);
}

