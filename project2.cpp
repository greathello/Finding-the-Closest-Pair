// project2.cpp : 定义控制台应用程序的入口点。
//

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <sstream>
#include <time.h>
#include <math.h>
#include <limits>


#define m 100
const int A = 48271;
const long M = 2147483647;
const int Q = M / A;
const int R = M % A;
long int state = 0;
using namespace std;

class S {
public:
	int ci_j;
	double vi_j;
};

class CP {
public:
	int ImageA;
	int ImageB;
	double dist = 0;
	void set(int a, int b) {
		ImageA = a;
		ImageB = b;
	}
	void set_dist(double d) {
		dist = d;
	}
	int* get() {
		int* a = new int[2];
		a[0] = ImageA;
		a[1] = ImageB;
		return a;
	}
};

int main(int argc, char* argv[]) {

	char* filename = NULL;
	int num_image;
	int demen;

	//编译参数获取
	for (int i = 0; i < argc; i++) {
		if (!strcmp(argv[i], "-f")) {
			if (i  < argc - 1) {
				filename = argv[i + 1];
			}
			else {
				cout << "Wrong!" << endl;
				return 0;
			}
		}
		if (!strcmp(argv[i], "-n")) {
			if (i  < argc - 1) {
				num_image = atoi(argv[i + 1]);
			}
			else {
				cout << "Wrong!" << endl;
				return 0;
			}
		}
		if (!strcmp(argv[i], "-d")) {
			if (i  < argc - 1) {
				demen = atoi(argv[i + 1]);
			}
			else {
				cout << "Wrong!" << endl;
				return 0;
			}
		}
	}

	int** object_list;
	object_list = new int*[num_image];
	for (int i = 0; i < num_image; i++) {
		object_list[i] = new int[demen];
	}


	FILE *file = fopen(filename, "rb");
	if (file == NULL) {
		cout << "The filename is null" << endl;
		abort();
	}

	//使用函数声明
	int change_rule(int);
	int* Random_Projection(int **, int, int, int);

	//文件读取
	int magic, numImage, rows, cols;
	fread(&magic, sizeof(int), 1, file);
	fread(&numImage, sizeof(int), 1, file);
	fread(&rows, sizeof(int), 1, file);
	fread(&cols, sizeof(int), 1, file);
	magic = change_rule(magic);
	numImage = change_rule(numImage);
	rows = change_rule(rows);
	cols = change_rule(cols);
	for (int i = 0; i < num_image; i++) {
		for (int j = 0; j < demen; j++) {
			unsigned char temp;
			fread(&temp, sizeof(unsigned char), 1, file);
			object_list[i][j] = change_rule(int(temp));
		}
	}
	fclose(file);

	//Algorithm 2:Randon Projection
	double point_mul(double*, int*, int);
	S** S1 = new S*[m];
	for (int i = 0; i < m; i++){
		S1[i] = new S[num_image];
	}
	void Normal(double*, int);
	double **a = new double*[m];
	for (int i = 0; i < m; i++) {
		a[i] = new double[demen];
		Normal(a[i], demen);
	}
	int c = 0;
	double v = 0;
	for (int i = 0; i < num_image; i++) {
		for (int j = 0; j < m; j++) {
			c = i;
			v = point_mul(a[j], object_list[i], demen);
			S1[j][i].ci_j = c;
			S1[j][i].vi_j = v;
		}
	}
	delete(a);

	CP cp;
	cp.set(0, 0);
	double min = numeric_limits<double>::max();
	CP Close_Pair_Line(S*, int);

	for (int i = 0; i < m; i++) {
		CP temp;
		int *a1 = new int[2];
		a1 = Close_Pair_Line(S1[i], num_image).get();
		temp.ImageA = a1[0];
		temp.ImageB = a1[1];
		double dist = abs(S1[i][temp.ImageA].vi_j - S1[i][temp.ImageB].vi_j);
		if (dist < min) {
			cp.ImageA = temp.ImageA;
			cp.ImageB = temp.ImageB;
			min = dist;
		}
		delete(a1);
	}

	cout << cp.ImageA << endl;
	for (int i = 0; i < 28; i++) {
		for (int j = 0; j < 28; j++) {
			if (object_list[cp.ImageA][i * 28 + j] == 0) {
				cout << "*";
			} else {
				cout << " ";
			}
		}
		cout << endl;
	}
	cout << endl << endl;

	cout << cp.ImageB << endl;
	for (int i = 0; i < 28; i++) {
		for (int j = 0; j < 28; j++) {
			if (object_list[cp.ImageB][i * 28 + j] == 0) {
				cout << "*";
			}
			else {
				cout << " ";
			}
		}
		cout << endl;
	}

	system("pause");
	return 0;
}

int change_rule(int value) {
	return ((value & 0x000000FF) << 24) |
		((value & 0x0000FF00) << 8) |
		((value & 0x00FF0000) >> 8) |
		((value & 0xFF000000) >> 24);
}

double point_mul(double* a, int* o, int d) {
	double res = 0;
	for (int i = 0; i < d; i++) {
		res += a[i] * o[i];
	}
	return res;
}

double Uniform() {
	while (state == 0) {
		srand((unsigned)time(NULL));
		state = rand();
	}

	double u;

	int tmpState = A * (state % Q) - R * (state / Q);
	if (tmpState >= 0)
		state = tmpState;
	else
		state = tmpState + M;

	u = state / (double)M;
	return u;
}

void Normal(double* nn, int len_nn)
{
	double x1, x2, w;
	int t;
	for (t = 0; 2 * t + 1 < len_nn; t++) {
		w = 2.0;
		while (w > 1.0) {
			x1 = 2.0 * Uniform() - 1.0;
			x2 = 2.0 * Uniform() - 1.0;
			w = x1 * x1 + x2 * x2;
		}
		w = sqrt(-2.0 * log(w) / w);
		nn[2 * t] = x1 * w;
		nn[2 * t + 1] = x2 * w;
	}
	if (len_nn % 2 == 1) {
		w = 2.0;
		while (w > 1.0) {
			x1 = 2.0 * Uniform() - 1.0;
			x2 = 2.0 * Uniform() - 1.0;
			w = x1 * x1 + x2 * x2;
		}
		w = sqrt(-2.0 * log(w) / w);
		nn[len_nn - 1] = x1 * w;
	}
	return;
}

CP Close_Pair_Line(S* s, int k) {
	CP cp, cp1, cp2;
	if (k == 0) {
		cp.set(0, 0);
		cp.set_dist(numeric_limits<double>::max());
		return cp;
	}
	if (k == 1) {
		cp.set(0, 0);
		cp.set_dist(numeric_limits<double>::max());
		return cp;
	}
	if (k == 2) {
		cp.set(s[0].ci_j, s[1].ci_j);
		cp.set_dist(abs(s[0].vi_j - s[1].vi_j));
		return cp;
	}
	int a = int(Uniform() * k);
	int temp_cp = 0;
	double temp_dist = numeric_limits<double>::max();
	int num1 = 0, num2 = 0;
	S* s1,* s2;
	for (int i = 0; i < k; i++) {
		if (i == a) {
		}
		else {
			if (s[i].vi_j < s[a].vi_j) {
				num1++;
			}
			else {
				num2++;
			}
		}
	}
	s1 = new S[num1];
	s2 = new S[num2];
	num1 = 0, num2 = 0;
	for (int i = 0; i < k; i++) {
		if (i == a){
		} else {
			if (s[i].vi_j < s[a].vi_j) {
				if (abs(s[i].vi_j - s[a].vi_j) < temp_dist) {
					temp_dist = abs(s[i].vi_j - s[a].vi_j);
					temp_cp = i;
				}
				s1[num1].ci_j = s[i].ci_j;
				s1[num1].vi_j = s[i].vi_j;
				num1++;
			} else {
				if (abs(s[i].vi_j - s[a].vi_j) < temp_dist) {
					temp_dist = abs(s[i].vi_j - s[a].vi_j);
					temp_cp = i;
				}
				s2[num2].ci_j = s[i].ci_j;
				s2[num2].vi_j = s[i].vi_j;
				num2++;
			}
		}
	}

	CP cp12;
	cp12.set(a, temp_cp);
	cp12.set_dist(temp_dist);
	cp1 = Close_Pair_Line(s1, num1);
	cp2 = Close_Pair_Line(s2, num2);

	if (cp1.dist < cp2.dist) {
		cp = cp1;
	} else {
		cp = cp2;
	}
	if (cp.dist < cp12.dist) {
		cp = cp12;
	}
	delete(s1);
	delete(s2);
	return cp;
}







