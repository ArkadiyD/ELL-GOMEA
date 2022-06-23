#include "problems.hpp"

void createProblemInstance(int problemIndex, int numberOfVariables, Problem **problemInstance, char *instancePath, int k, int s)
{
	switch (problemIndex)
	{
		case 1: *problemInstance = new concatenatedDeceptiveTrap(k, s); break;
		case 2: *problemInstance = new ADF(instancePath); break;
		case 3: *problemInstance = new maxCut(instancePath); break;		
		case 4: *problemInstance = new MAXSAT(instancePath); break;
		case 5: *problemInstance = new WalshDecomposition(instancePath); break;
		
		default:
		{
			cerr << "No problem with index #" << problemIndex << " installed!\n";
			exit(0);
		};
	}
	(*problemInstance)->initializeProblem(numberOfVariables);
}


double deceptiveTrap(int unitation, int k)
{
	double result;
	if (unitation != k)
		result = k - 1 - unitation;
	else
		result = unitation;
	
	return (double)result;
}

void concatenatedDeceptiveTrap::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
};
	
	
double concatenatedDeceptiveTrap::calculateFitness(vector<bool> &solution)
{
	double res = 0.0;
	for (int i = 0; i < numberOfVariables; i += s)
	{
		int unitation = 0;
		for (int j = i; j < i + k; j++)
		{
			unitation += solution[j % numberOfVariables];
		}

		res += deceptiveTrap(unitation, k);
	}

	return res;
}


void ADF::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	
	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	int N, numFunctions;
	inFile >> N >> numFunctions;
	
	for (int i = 0; i < numFunctions; ++i)
	{
		NKSubfunction subfunction;

		int k;
		inFile >> k;

		for (int j = 0; j < k; ++j)
		{
			int var;
			inFile >> var;
			subfunction.variablesPositions.push_back(var);
		}
		int numCombinations = 1 << k;
		subfunction.valuesTable.resize(numCombinations);
		for (int j = 0; j < numCombinations; ++j)
		{
			string combination;
			inFile >> combination;
			int index = 0, pow = 1;
			for (size_t p = combination.size()-1; p > 0; p--)
			{
				if (combination[p] == '0' || combination[p] == '1')
				{
					index += ((int)(combination[p] == '1') * pow);
					pow *= 2;
				}
			}
			inFile >> subfunction.valuesTable[index];			

		}
		subfunctions.push_back(subfunction);
	}
	inFile.close();
	cout << "init finished\n";
};

double ADF::calculateFitness(vector<bool> &solution)
{
	double res = 0.0;
	
	for (size_t i = 0; i < subfunctions.size(); ++i)
	{
		int index = 0, pow = 1;
		for (int j = subfunctions[i].variablesPositions.size()-1; j >= 0; --j)
		{
			int pos = subfunctions[i].variablesPositions[j];
			index += solution[pos] * pow;
			pow *= 2;
		}
		res += subfunctions[i].valuesTable[index];// * multiplier;
	}
	
	return (double)res;
}

void maxCut::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	edgesForVariable.resize(numberOfVariables);

	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	int N, numEdges;
	inFile >> N >> numEdges;

	for (int i = 0; i < numEdges; ++i)
	{
		int v1, v2;
		double w;
		inFile >> v1 >> v2 >> w;
		edges.push_back(make_pair(make_pair(v1-1, v2-1), w));
		//cout << v1 << " " << v2 << " " << w << endl;
		edgesForVariable[v1-1].push_back(i);
		edgesForVariable[v2-1].push_back(i);		
	}

	inFile.close();
}


double maxCut::calculateFitness(vector<bool> &solution)
{
	double res = 0.0;
	for (size_t i = 0; i < edges.size(); ++i)
	{
		int v1 = edges[i].first.first;
		int v2 = edges[i].first.second;
		double w = edges[i].second;

		if (solution[v1] != solution[v2])
			res += w;
	}

	return res;
}

void MAXSAT::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	subfunctionsForVariable.resize(numberOfVariables);

	ifstream inFile(problemInstanceFilename, ifstream::in);
	if (inFile.fail())
	{
		cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
	
	string line;
	for (int i = 0; i < 8; ++i)
		getline(inFile, line);

	while(true)
	{
		vector<int> subfunction;
		vector<int> subfunctionSigns;
		int var;
		inFile >> var;
		if (inFile.fail())
			break;
		if (var < 0)
		{
			subfunctionSigns.push_back(-1);
			var++;
		}
		else if (var > 0)
		{
			var--;
			subfunctionSigns.push_back(1);
		}
		var = abs(var);
		subfunction.push_back(var);

		while (true)
		{
			inFile >> var;
			if (var == 0)
				break;

			if (var < 0)
			{
				subfunctionSigns.push_back(-1);
				var++;
			}
			else if (var > 0)
			{
				var--;
				subfunctionSigns.push_back(1);
			}

			var = abs(var);
			subfunction.push_back(var);
		}
		for (int i = 0; i < subfunction.size(); ++i)
		{
			//cout << subfunction[i] << " | " << subfunctionSigns[i] << " ";
			subfunctionsForVariable[subfunction[i]].push_back(subfunctions.size());
		}
		//cout << endl;
		subfunctions.push_back(subfunction);
		signs.push_back(subfunctionSigns);

	}
	inFile.close();
}


double MAXSAT::calculateFitness(vector<bool> &solution)
{
	long double res = 0.0;
	
	for (size_t i = 0; i < subfunctions.size(); ++i)
	{
		bool b = false;
		for (int j = 0; j < subfunctions[i].size(); ++j)
		{
			int var = subfunctions[i][j];
			int sign = signs[i][j];
			if (sign > 0 && solution[var] == 1)
				b = true;
			else if (sign < 0 && solution[var] == 0)
				b = true;			
		}
		if (b == false)
			res -= 1;		
	}
	
	return (double)res;
}

void WalshDecomposition::initializeProblem(int numberOfVariables_)
{
	numberOfVariables = numberOfVariables_;
	walsh_coefficients_by_index.resize(numberOfVariables);

	char    c, *string, substring[1000], filename[1000];
  	int     i, j, k, q;
  	double *values;
  	FILE   *file;

  	//sprintf( filename, "%s/walsh_coefficients.dat", folder);
  	file = fopen(problemInstanceFilename.c_str(), "r");
  	if( file == NULL)
  	{
	    cout << "Problem Instance File " << problemInstanceFilename << " does not exist!\n";
		exit(0);
	}
  
	fseek(file, 0 , SEEK_END);
	long fileSize = ftell(file);
	fseek(file, 0 , SEEK_SET);
	string = new char[fileSize * 4];

	pair<vector<int>, double> walsh_coefficient;

	c = fgetc( file );
	k = 0;
	while( c != EOF )
	{
		string[k] = (char) c;
		c      = fgetc( file );
		k++;
	}
	string[k] = '\0';
	//printf("%s %d\n", string, k);
	i = 0;
	while (i < k)
	{
	//printf("i %d\n",i);
		int cur_number_of_variables = 0;
		while(cur_number_of_variables != numberOfVariables)
		{
		  //printf("%c %d %d\n", string[i], i, number_of_variables);
			if( string[i] == '1' )
			{
				walsh_coefficient.first.push_back(cur_number_of_variables);
			    //printf("var %d\n", number_of_variables);

				walsh_coefficients_by_index[cur_number_of_variables].push_back(walsh_coefficients.size());

				cur_number_of_variables++;
			}
			else if( string[i] == '0' )
				cur_number_of_variables++;

			i++;
		}

		q = 0;
		while(string[i] != '\n' && i < k)
		{
			if (string[i] == ' ')
			{
		 		i++;
			continue;
		}

		substring[q] = string[i];
		  //printf("%c", substring[q]);
		q++;
		i++;
	}
	i++;
	substring[q] = '\0';

	walsh_coefficient.second = atof( substring );
	//printf("%lf\n",  walsh_coefficient.second);
	//for (int p = 0; p < walsh_coefficient.first.size();++p)
	  //printf("%d ", walsh_coefficient.first[p]);
	walsh_coefficients.push_back(walsh_coefficient);
	walsh_coefficient.first.clear();

	}
	fclose( file );

	delete[] string;

	walshCoefficientsTouched.resize(walsh_coefficients.size());
}

double WalshDecomposition::walshCoefficientContribution(vector<bool> &solution, vector<int> &coefficient_mask, double coefficient)
{
  int dot_product = 0;
  for (int i = 0; i < coefficient_mask.size(); ++i)
  {
    if (solution[coefficient_mask[i]] == 1)
      dot_product++;
  }

  if (dot_product % 2 == 0)
    dot_product = 1;
  else
    dot_product = -1;

  return (double)dot_product * coefficient;
}

double WalshDecomposition::calculateFitness(vector<bool> &solution)
{
  int    i;
  double result = 0.0;

  for( i = 0; i < walsh_coefficients.size(); i++ )
  {    
    result += walshCoefficientContribution(solution, walsh_coefficients[i].first, walsh_coefficients[i].second);
    //printf("%d %lf\n", i, walshCoefficientContribution(solution, walsh_coefficients[i].first, walsh_coefficients[i].second));
  }
  return result;
}

double WalshDecomposition::calculateFitnessPartialEvaluations(vector<bool> &solution, vector<int> &touched_genes_indices, vector<bool> &genes_before, double fitness_before, double &newEvals)
{
	newEvals = 0;
	fill(walshCoefficientsTouched.begin(), walshCoefficientsTouched.end(), 0);
	double res = fitness_before;

	for(int i = 0; i < touched_genes_indices.size(); i++)
	{
    	int ind = touched_genes_indices[i];
    	for (int j = 0; j < walsh_coefficients_by_index[ind].size(); ++j)
    	{
    		int coef_ind = walsh_coefficients_by_index[ind][j];
    		if (walshCoefficientsTouched[coef_ind] == 0)
    		{
        		res -= walshCoefficientContribution(genes_before, walsh_coefficients[coef_ind].first, walsh_coefficients[coef_ind].second);
        		res += walshCoefficientContribution(solution, walsh_coefficients[coef_ind].first, walsh_coefficients[coef_ind].second);
      		}
      		//cout << walshCoefficientContribution(genes_before, walsh_coefficients[coef_ind].first, walsh_coefficients[coef_ind].second) << " " << walshCoefficientContribution(solution, walsh_coefficients[coef_ind].first, walsh_coefficients[coef_ind].second) << endl;
	        walshCoefficientsTouched[coef_ind] = 1;  
	    }
	}
  
  newEvals = touched_genes_indices.size()/((double)numberOfVariables);

  return res;
}
