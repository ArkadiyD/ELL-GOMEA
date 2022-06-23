#pragma once

#include <vector>
#include <unordered_map>
#include <fstream>
#include <string>
#include <deque>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <cassert>
#include <iostream>
using namespace std;

class Problem
{
public:
	int numberOfVariables;
	
	Problem(){};
	virtual ~Problem(){};
	virtual void initializeProblem(int numberOfVariables)=0;
	virtual double calculateFitness(vector<bool> &solution)=0;
	virtual double calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals)=0;
};


class concatenatedDeceptiveTrap:public Problem
{
	int k, s;
	bool bimodal;
	vector<vector<int> > trapsForVariable; //determines traps which each variables belongs to

public:
	concatenatedDeceptiveTrap(int k_, int s_): k(k_), s(s_)
	{
		cout<<"creating concatenated Deceptive Trap with trap size=" << k << " and shift=" << s << endl;
	}
	void initializeProblem(int numberOfVariables_);
	double calculateFitness(vector<bool> &solution);
	double calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals)
	{
		cout << "Partial evalautions not implemented!";
		exit(0);
	}
};

struct NKSubfunction
{
	vector<int> variablesPositions;
	vector<double> valuesTable;
};

class ADF:public Problem
{
	string problemInstanceFilename;
	vector<NKSubfunction> subfunctions;
	vector<vector<int> > subfunctionsForVariable;

public:
	ADF(char* problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating ADF\n";
	}

	void initializeProblem(int numberOfVariables_);
	double calculateFitness(vector<bool> &solution);
	double calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals)
	{
		cout << "Partial evalautions not implemented!";
		exit(0);
	}

};

class maxCut:public Problem
{
	string problemInstanceFilename;
	vector<pair<pair<int, int>, double > > edges;
	vector<vector<int> > edgesForVariable;

public:
	maxCut(char* problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating maxCut\n";
	}
	void initializeProblem(int numberOfVariables_);
	double calculateFitness(vector<bool> &solution);
	double calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals)
	{
		cout << "Partial evalautions not implemented!";
		exit(0);
	}
};

class MAXSAT:public Problem
{
	string problemInstanceFilename;
	vector<vector<int> > subfunctions;
	vector<vector<int> > signs;
	
	vector<vector<int> > subfunctionsForVariable;

public:
	MAXSAT(char* problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating MAXSAT " << problemInstanceFilename << endl;
	}

	void initializeProblem(int numberOfVariables_);

	double calculateFitness(vector<bool> &solution);
	double calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals)
	{
		cout << "Partial evalautions not implemented!";
		exit(0);
	}
};

class WalshDecomposition:public Problem
{
	string problemInstanceFilename;
	vector<pair<vector<int>, double> > walsh_coefficients;
	vector<vector<int> > walsh_coefficients_by_index;
	vector<bool> walshCoefficientsTouched;
	
	double walshCoefficientContribution(vector<bool> &solution, vector<int> &coefficient_mask, double coefficient);

public:
	WalshDecomposition(char* problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
	{
		cout<<"creating Walsh Decomposition " << problemInstanceFilename << endl;
	}

	void initializeProblem(int numberOfVariables_);

	double calculateFitness(vector<bool> &solution);
	double calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals);

};

double deceptiveTrap(int unitation, int k);

void createProblemInstance(int problemIndex, int numberOfVariables, Problem **problemInstance, char *instancePath, int k = 1, int s = 1);

