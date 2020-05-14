#include<cstdio>
#include<cstring>
#include<iostream>
#include <fstream>
#include<ctime>
#include<random>
#include<cstdlib>
#define NMAX 510
#define INF 10000009
using namespace std;

mt19937 gen(time(NULL));
uniform_real_distribution<double> distribution(0, 1);

double fer[NMAX][NMAX];
double dist[NMAX][NMAX];
bool blacklist[NMAX];
double chance[2][NMAX];

int minRoad[NMAX];
int road[NMAX][NMAX];
int roadDist[NMAX];

void RubishDestroy(int n) {
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			road[i][j] = 0;
		}
		roadDist[i] = 0;
	}
}

double sumForP(int i, int n) {
	const double alpha = 3.645, betta = 0.754;
	double sum = 0.0;
	for (int j = 1; j <= n; j++)
		if (!blacklist[j])
			sum += pow((double) 1 / (double) dist[i][j], alpha) * pow(fer[i][j], betta);
	return sum;
}

double P(int i, int j, double sum) {
	const double alpha = 3.645, betta = 0.754;
	return 100 * (pow((double)1 / (double) dist[i][j], alpha) * pow(fer[i][j], betta)) / sum;
}


void AntRoad(int AntNum, int i, int n, bool elite, int NotCompleteCity) {
	double sum = 0.0;
	sum = sumForP(i, n);
	int countChance = 0;
	for (int j = 1; j <= n; j++) {
		if (!blacklist[j]) {
			chance[0][countChance] = P(i, j, sum);
			chance[1][countChance] = j;
			countChance++;
		}
	}
	if (countChance == 1) {
		road[AntNum][n - NotCompleteCity] = chance[1][0];
		int w = (int)chance[1][0];
		roadDist[AntNum] += dist[i][w];
		blacklist[w] = true;
		roadDist[AntNum] += dist[w][AntNum];
		road[AntNum][n] = AntNum;
		return;
	}

	if (!elite) {
		for (int j = 1; j < countChance; j++) {
			chance[0][j] += chance[0][j - 1];
		}
		double choice = 100 * distribution(gen);

		for (int j = 0; j < countChance; j++) {
			if (choice < chance[0][j]) {
				road[AntNum][n - NotCompleteCity] = chance[1][j];
				int w = (int)chance[1][j];
				roadDist[AntNum] += dist[i][w];
				blacklist[w] = true;
				AntRoad(AntNum, w, n, elite, NotCompleteCity - 1);
				return;
			}
		}
	}
	if (elite) {
		double maxChance = 0.0f;
		int w;
		for (int j = 0; j < countChance; j++) {
			if (chance[0][j] > maxChance) {
				maxChance = chance[0][j];
				w = chance[1][j];
			}
		}
		roadDist[AntNum] += dist[i][w];
		blacklist[w] = true;
		AntRoad(AntNum, w, n, elite, NotCompleteCity - 1);
		return;
	}
}

int AntAlgorithm(int n) {
	int iter = n * 50;
	int min_len = INF;
	const double p = 0.5;
	bool elite;
	while (iter) {
		elite = false;
		for (int i = 1; i <= n; i++) {
			blacklist[i] = true;
			road[i][0] = i;
			AntRoad(i, i, n, elite, n - 1);
			if (min_len > roadDist[i]) {
				min_len = roadDist[i];
				for (int j = 0; j <= n; j++)
					minRoad[j] = road[i][j];
				}
			for (int j = 1; j <= n; j++)
				blacklist[j] = 0;
		}
		for (int i = 1; i <= n; i++)
			for (int j = 1; j <= n; j++)
				fer[i][j] *= (1 - p);

		for (int i = 1; i <= n; i++)
			for (int j = 1; j <= n; j++)
				fer[road[i][j - 1]][road[i][j]] += (1.0 / (double) roadDist[i]);
		RubishDestroy(n);
		iter--;
	}
	return min_len;
}



int main(void) {
	ifstream fin;
	fin.open("att48_d.txt");

	int n, u, v, w, len = 0, schet = 0;
	fin >> n;
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++) {
			fin >> dist[i][j];
			fer[i][j] = 0.525;
		}
	
	len = AntAlgorithm(n);
	
	cout << len << "\n";
	for (int i = 0; i <= n; i++)
		cout << minRoad[i] << " ";

	return 0;
}
