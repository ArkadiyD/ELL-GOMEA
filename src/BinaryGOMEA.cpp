/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Header -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * BinaryGOMEA.c
 *
 * Copyright (c) Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 * Genepool Optimal Mixing Evolutionary Algorithm
 *
 * In this implementation, maximization is assumed.
 *
 * Recommended command line for compiling under LINUX:
 * g++ -Wall -O3 -o BinaryGOMEA BinaryGOMEA.c
 */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define OS_WIN
#endif
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <stdio.h>
#include <stdlib.h>
#ifdef OS_WIN
#include <stdint.h>
#endif
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "utils.hpp"
#include "problems.hpp"

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
char     ***BinaryGOMEA_populations,                                                        /* The BinaryGOMEA_populations containing the solutions. */
         ***BinaryGOMEA_offsprings,                                                         /* Offspring solutions (one set per population). */
           *BinaryGOMEA_terminated;                                                         /* Whether a specific BinaryGOMEA with the restart scheme has terminated. */
short       BinaryGOMEA_print_verbose_overview,                                             /* Whether to print a overview of settings (0 = no). */
            BinaryGOMEA_write_FOS_contents_to_file=1,                                         /* Whether to write the contents of the FOS structure to file each generation (0 = no). */
            BinaryGOMEA_write_MI_matrix_to_file,                                            /* Whether to write the contents of the MI matrix to file each generation (0 = no). */
            BinaryGOMEA_elitist_solution_written_before,                                    /* Whether the elitist solution was written to file before. */
            BinaryGOMEA_first_evaluation_ever,                                              /* Whether the next evaluation to be performed is actually the first evaluation ever. */
            BinaryGOMEA_use_VTR,                                                            /* Whether to use the VTR. */
            BinaryGOMEA_use_partial_evaluations,                                            /* Whether to use partial evaluations. */
            BinaryGOMEA_report_every_improvement,                                           /* Whether to report every improvement. */
            BinaryGOMEA_use_predetermined_fos;                                              /* Whether to use a predetermined FOS. */
int         BinaryGOMEA_problem_index,                                                      /* The index of the optimization problem. */
            BinaryGOMEA_number_of_genes,                                                    /* The number of genes (or gene variables) to be optimized. */
           *BinaryGOMEA_number_of_generations,                                              /* The number of generations in the i-th multi-start instance. */
            BinaryGOMEA_number_of_generations_IMS,                                          /* The current generation count in the Interleaved Multi-start Scheme (IMS). */
           *BinaryGOMEA_population_sizes,                                                   /* The number of solutions in each population. */
            BinaryGOMEA_base_population_size,                                               /* The minimum population size used in the smallest BinaryGOMEA instance. */
            BinaryGOMEA_IMS_subgeneration_factor,                                           /* The factor by which the number of subgenerations increases with every new population in the IMS. */
          **BinaryGOMEA_no_improvement_stretchs,                                            /* The number of subsequent generations without an improvement for every individual. */
          **BinaryGOMEA_predetermined_FOS,                                                  /* The family of subsets linkage struture, one per population in the IMS. */
           *BinaryGOMEA_predetermined_FOS_number_of_indices,                                /* The number of variables in each linkage subset, one per population in the IMS. */
            BinaryGOMEA_predetermined_FOS_length,                                           /* The number of linkage subsets, one per population in the IMS. */
         ***BinaryGOMEA_FOSs,                                                               /* The family of subsets linkage struture, one per population in the IMS. */
          **BinaryGOMEA_FOSs_number_of_indices,                                             /* The number of variables in each linkage subset, one per population in the IMS. */
           *BinaryGOMEA_FOSs_length,                                                        /* The number of linkage subsets, one per population in the IMS. */
            BinaryGOMEA_minimum_BinaryGOMEA_index,                                          /* The minimum BinaryGOMEA index that corresponds to the BinaryGOMEA that is still allowed to run (lower ones should be stopped because of average fitness being lower than that of a higher one). */
            BinaryGOMEA_number_of_BinaryGOMEAs,                                             /* The number of BinaryGOMEAs currently running in multipop configuration. */
            BinaryGOMEA_maximum_number_of_BinaryGOMEAs,                                     /* The maximum number of BinaryGOMEAs running in multipop configuration. */
            BinaryGOMEA_elitist_solution_BinaryGOMEA_index,                                 /* The BinaryGOMEA index for which the offspring vector currently holds the best solution ever evaluated. */
            BinaryGOMEA_elitist_solution_offspring_index;                                   /* The offspring index that currently points to the best solution ever evaluated. */
long        BinaryGOMEA_timestamp_start,                                                    /* The time stamp in milliseconds for when the program was started. */
            BinaryGOMEA_timestamp_start_after_init,                                         /* The time stamp in milliseconds for when the algorithm was started (after problem initialization). */
            BinaryGOMEA_elitist_solution_hitting_time_in_milliseconds,                      /* The hitting time of the elitist solution in number of milliseconds. */
            BinaryGOMEA_elitist_solution_last_written_running_time,                         /* The last time the elitist solutions was written to a file. */
            BinaryGOMEA_vtr_hitting_time;                                                   /* The VTR hitting time. */
double      BinaryGOMEA_number_of_evaluations,                                              /* The current number of times a function evaluation was performed. */
            BinaryGOMEA_elitist_solution_hitting_time_in_evaluations,                       /* The hitting time of the elitist solution in number of evaluations. */
            BinaryGOMEA_vtr,                                                                /* The value to reach (fitness of best solution that is feasible). */
          **BinaryGOMEA_fitnesses,                                                          /* BinaryGOMEA_fitnesses of population members. */
          **BinaryGOMEA_fitnesses_offsprings,                                               /* BinaryGOMEA_fitnesses of offspring. */
           *BinaryGOMEA_average_fitnesses,                                                  /* The average fitness of each BinaryGOMEA population. */
         ***BinaryGOMEA_MI_matrices;                                                        /* Mutual information between any two variables. */
long long     BinaryGOMEA_random_seed;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

char folder[100000];
Problem *problem_instance;
char instance_path[1000];
int k,s;

long long timelimitMilliseconds;
std::unordered_map< std::vector<bool>, double, hash_vector > evaluated_solutions;
int SKIP_ALREADY_EVALUATED;

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
void *BinaryGOMEA_Malloc( long size );
int *BinaryGOMEA_MergeSortIntegersDecreasing( int *array, int array_size );
void BinaryGOMEA_MergeSortIntegersDecreasingWithinBounds( int *array, int *sorted, int *tosort, int p, int q );
void BinaryGOMEA_MergeSortIntegersDecreasingMerge( int *array, int *sorted, int *tosort, int p, int r, int q );
void BinaryGOMEA_InterpretCommandLine( int argc, char **argv );
void BinaryGOMEA_ParseCommandLine( int argc, char **argv );
void BinaryGOMEA_ParseOptions( int argc, char **argv, int *index );
void BinaryGOMEA_PrintAllInstalledProblems();
void BinaryGOMEA_OptionError( char **argv, int index );
void BinaryGOMEA_ParseParameters( int argc, char **argv, int *index );
void BinaryGOMEA_PrintUsage();
void BinaryGOMEA_CheckOptions();
void BinaryGOMEA_PrintVerboseOverview();
double BinaryGOMEA_RandomRealUniform01();
int BinaryGOMEA_RandomInt( int maximum );
int *BinaryGOMEA_RandomPermutation( int n );

/* INSTALLED PROBLEMS */
char *BinaryGOMEA_InstalledProblemName();
int BinaryGOMEA_NumberOfInstalledProblems();
void BinaryGOMEA_InstalledProblemEvaluationSpecificBinaryGOMEASpecificOffspringIndex( int BinaryGOMEA_index, int offspring_index, int number_of_touched_genes, int *touched_genes_indices, char *genes_before, double fitness_before );

void BinaryGOMEA_InitializeNewBinaryGOMEA();
void BinaryGOMEA_InitializeNewBinaryGOMEAMemory();
void BinaryGOMEA_InitializeNewBinaryGOMEAPopulationAndFitnessValues();
void BinaryGOMEA_InitializePredeterminedFOS();
void BinaryGOMEA_InitializeRandomNumberGenerator();
void BinaryGOMEA_InitializeProblem();
void BinaryGOMEA_CopyOffspringToPopulationSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_computeAverageFitnessSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_WriteElitistSolutionToFile();
char BinaryGOMEA_CheckTermination();
char BinaryGOMEA_CheckTerminationSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_GenerationalStepAllBinaryGOMEAs();
void BinaryGOMEA_GenerationalStepAllBinaryGOMEAsRecursiveFold( int BinaryGOMEA_index_smallest, int BinaryGOMEA_index_biggest );
void BinaryGOMEA_MakeOffspringSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_LearnLTFOSSpecificBinaryGOMEA( int BinaryGOMEA_index );
int BinaryGOMEA_DetermineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length );
void BinaryGOMEA_ComputeMIMatrixSpecificBinaryGOMEA( int BinaryGOMEA_index );
double *BinaryGOMEA_EstimateParametersForSingleBinaryMarginalSpecificBinaryGOMEA( int BinaryGOMEA_index, int *indices, int number_of_indices, int *factor_size );
void BinaryGOMEA_WriteToFileFOSContentsSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_WriteToFileMIMatrixSpecificBinaryGOMEA( int BinaryGOMEA_index );
double BinaryGOMEA_Log2( double x );
void BinaryGOMEA_GenerateAndEvaluateOffspringSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_GenerateAndEvaluateNewSolutionSpecificBinaryGOMEASpecificOffspringIndex( int BinaryGOMEA_index, int offspring_index );
void BinaryGOMEA_GenerateAndEvaluateNewSolutionSpecificBinaryGOMEASpecificOffspringIndexEvaluationInfoInitializeBackup( int BinaryGOMEA_index, int offspring_index );
void BinaryGOMEA_GenerateAndEvaluateNewSolutionSpecificBinaryGOMEASpecificOffspringIndexEvaluationInfoEzilaitiniBackup();
void BinaryGOMEA_ShuffleFOSSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_EzilaitiniAllBinaryGOMEAs();
void BinaryGOMEA_EzilaitiniSpecificBinaryGOMEA( int BinaryGOMEA_index );
void BinaryGOMEA_EzilaitiniSpecificBinaryGOMEAMemoryForPopulationAndOffspring( int BinaryGOMEA_index );
void BinaryGOMEA_EzilaitiniProblem();
void BinaryGOMEA_RunIMS();
void BinaryGOMEA_Run();

int main( int argc, char **argv );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-= Section Elementary Operations -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Allocates memory and exits the program in case of a memory allocation failure.
 */
void *BinaryGOMEA_Malloc( long size )
{
  void *result;

  result = (void *) malloc( size );
  if( !result )
  {
    printf( "\n" );
    printf( "Error while allocating memory in BinaryGOMEA_Malloc( %ld ), aborting program.", size );
    printf( "\n" );

    exit( 0 );
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Merge Sort -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Sorts an array of integers and returns the sort-order (large to small).
 */
int *BinaryGOMEA_MergeSortIntegersDecreasing( int *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) BinaryGOMEA_Malloc( array_size * sizeof( int ) );
  tosort = (int *) BinaryGOMEA_Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    BinaryGOMEA_MergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void BinaryGOMEA_MergeSortIntegersDecreasingWithinBounds( int *array, int *sorted, int *tosort, int p, int q )
{
  int r;
  if( p < q )
  {
    r = (p + q) / 2;
    BinaryGOMEA_MergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, p, r );
    BinaryGOMEA_MergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, r+1, q );
    BinaryGOMEA_MergeSortIntegersDecreasingMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void BinaryGOMEA_MergeSortIntegersDecreasingMerge( int *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] > array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/




/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
void BinaryGOMEA_InterpretCommandLine( int argc, char **argv )
{
  BinaryGOMEA_ParseCommandLine( argc, argv );
  
  BinaryGOMEA_CheckOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void BinaryGOMEA_ParseCommandLine( int argc, char **argv )
{
  int index;

  index = 1;

  BinaryGOMEA_ParseOptions( argc, argv, &index );
  
  BinaryGOMEA_ParseParameters( argc, argv, &index );
}

/**
 * Parses only the options from the command line.
 */
void BinaryGOMEA_ParseOptions( int argc, char **argv, int *index )
{
  double dummy;

  BinaryGOMEA_use_partial_evaluations    = 0;
  BinaryGOMEA_report_every_improvement   = 0;
  BinaryGOMEA_use_VTR                    = 0;
  BinaryGOMEA_print_verbose_overview     = 0;
  BinaryGOMEA_write_FOS_contents_to_file = 0;
  BinaryGOMEA_use_predetermined_fos      = 0;
  BinaryGOMEA_print_verbose_overview     = 1;
  SKIP_ALREADY_EVALUATED                 = 0;


  for( ; (*index) < argc; (*index)++ )
  {
    if( argv[*index][0] == '-' )
    {
      /* If it is a negative number, the option part is over */
      if( sscanf( argv[*index], "%lf", &dummy ) && argv[*index][1] != '\0' )
        break;

      if( argv[*index][1] == '\0' )
        BinaryGOMEA_OptionError( argv, *index );
      else if( argv[*index][2] != '\0' )
        BinaryGOMEA_OptionError( argv, *index );
      else
      {
        switch( argv[*index][1] )
        {
          case '?': BinaryGOMEA_PrintUsage(); break;
          case 'P': BinaryGOMEA_PrintAllInstalledProblems(); break;
          case 'e': BinaryGOMEA_use_partial_evaluations    = 1; break;
          case 'i': BinaryGOMEA_report_every_improvement   = 1; break;
          case 'r': BinaryGOMEA_use_VTR                    = 1; break;
          case 'v': BinaryGOMEA_print_verbose_overview     = 1; break;
          case 'f': BinaryGOMEA_write_FOS_contents_to_file = 1; break;          
          case 'F': BinaryGOMEA_use_predetermined_fos      = 1; break;          
          case 'S': SKIP_ALREADY_EVALUATED     = 1; break;
          
          default : BinaryGOMEA_OptionError( argv, *index );
        }
      }
    }
    else /* Argument is not an option, so option part is over */
     break;
  }
}

/**
 * Writes the names of all installed problems to the standard output.
 */
void BinaryGOMEA_PrintAllInstalledProblems()
{
  int n, backup;

  n      = BinaryGOMEA_NumberOfInstalledProblems();
  backup = BinaryGOMEA_problem_index;

  printf( "Installed optimization problems:\n" );
  for( BinaryGOMEA_problem_index = 0; BinaryGOMEA_problem_index < n; BinaryGOMEA_problem_index++ )
    printf( "%3d: %s\n", BinaryGOMEA_problem_index, BinaryGOMEA_InstalledProblemName() );

  BinaryGOMEA_problem_index = backup;

  exit( 0 );
}

/**
 * Informs the user of an illegal option and exits the program.
 */
void BinaryGOMEA_OptionError( char **argv, int index )
{
  printf( "Illegal option: %s\n\n", argv[index] );

  BinaryGOMEA_PrintUsage();
}

/**
 * Parses only the parameters from the command line.
 */
void BinaryGOMEA_ParseParameters( int argc, char **argv, int *index )
{
  int noError;

  noError = 1;
  int cnt = 0;
  noError = noError && sscanf( argv[*index+cnt], "%d", &BinaryGOMEA_problem_index );
  cnt++;
  noError = noError && sscanf( argv[*index+cnt], "%d", &BinaryGOMEA_number_of_genes );
  cnt++;
  noError = noError && sscanf( argv[*index+cnt], "%d", &BinaryGOMEA_maximum_number_of_BinaryGOMEAs );
  cnt++;
  noError = noError && sscanf( argv[*index+cnt], "%s", &folder );
  cnt++;
  noError = noError && sscanf( argv[*index+cnt], "%lld", &timelimitMilliseconds );
  cnt++;
  noError = noError && sscanf( argv[*index+cnt], "%lld", &BinaryGOMEA_random_seed );
  cnt++;

  if (BinaryGOMEA_problem_index == 1 || BinaryGOMEA_problem_index == 2)
  {
    noError = noError && sscanf(argv[*index+cnt], "%d", &k);
    cnt++;
    noError = noError && sscanf(argv[*index+cnt], "%d", &s);
    cnt++;
  }

  if (BinaryGOMEA_problem_index == 2 || BinaryGOMEA_problem_index == 3 || BinaryGOMEA_problem_index == 4 || BinaryGOMEA_problem_index == 5)
  {
    noError = noError && sscanf(argv[*index+cnt], "%s", &instance_path);
    cnt++;
  }


  printf("%s\n", folder);

  if( !noError )
  {
    printf( "Error parsing parameters.\n\n" );

    BinaryGOMEA_PrintUsage();
  }
}

/**
 * Prints usage information and exits the program.
 */
void BinaryGOMEA_PrintUsage()
{
  printf( "Usage: BinaryGOMEA [-?] [-P] [-e] [-f] [-i] [-p] [-r] [-v] pro dim vtr mng\n" );
  printf( "   -?: Prints out this usage information.\n" );
  printf( "   -P: Prints out a list of all installed optimization problems.\n" );
  printf( "   -e: Enables use of partial evaluations (if fitness function implementation supports it).\n" );
  printf( "   -i: Enables writing to file every solution that improves fitness.\n" );
  printf( "   -r: Enables use of vtr in termination condition (value-to-reach).\n" );
  printf( "   -v: Enables verbose mode. Prints the settings before starting the run.\n" );
  printf( "\n" );
  printf( "  pro: Index of optimization problem to be solved (minimization).\n" );
  printf( "  dim: Number of genes.\n" );
  //printf( "  vtr: The value to reach. If the objective value of the best feasible solution\n" );
  //printf( "       reaches this value, termination is enforced (if -r is specified).\n" );
  printf( "  mng: Maximum number of GOMEA instances in the Interleaved Multi-start Scheme allowed.\n" );
  printf( "  folder: Folder where GOMEA runs.\n" );
  printf( "  timeLimit: timeLimit in milliseconds.\n" );
  
  
  exit( 0 );
}

/**
 * Checks whether the selected options are feasible.
 */
void BinaryGOMEA_CheckOptions()
{
  if( BinaryGOMEA_number_of_genes < 1 )
  {
    printf( "\n" );
    printf( "Error: number of genes < 1 (read: %d). Require number of genes >= 1.", BinaryGOMEA_number_of_genes);
    printf( "\n\n" );

    exit( 0 );
  }

  if( BinaryGOMEA_maximum_number_of_BinaryGOMEAs <= 0 )
  {
    printf( "\n" );
    printf( "Error: BinaryGOMEA_maximum_number_of_BinaryGOMEAs should be at least 1 (read %d)\n", BinaryGOMEA_maximum_number_of_BinaryGOMEAs);
    printf( "\n\n" );
    exit( 0 );
  }

  if( BinaryGOMEA_InstalledProblemName() == NULL )
  {
    printf( "\n" );
    printf( "Error: unknown index for problem (read index %d).", BinaryGOMEA_problem_index );
    printf( "\n\n" );

    exit( 0 );
  }
}

/**
 * Prints the settings as read from the command line.
 */
void BinaryGOMEA_PrintVerboseOverview()
{
  printf( "### Settings ######################################\n" );
  printf( "#\n" );
  printf( "# FOS contents writing to file : %s\n", BinaryGOMEA_write_FOS_contents_to_file ? "enabled" : "disabled" );
  printf( "# MI matrix writing to file    : %s\n", BinaryGOMEA_write_MI_matrix_to_file ? "enabled" : "disabled" );
  printf( "# Use of VTR                   : %s\n", BinaryGOMEA_use_VTR ? "enabled" : "disabled" );
  printf( "# Verbose mode                 : %s\n", BinaryGOMEA_print_verbose_overview ? "enabled" : "disabled" );
  printf( "#\n" );
  printf( "###################################################\n" );
  printf( "#\n" );
  printf( "# Problem                      = %s\n", BinaryGOMEA_InstalledProblemName());
  printf( "# Problem instance             = %s\n", instance_path);
  printf( "# Number of genes              = %d\n", BinaryGOMEA_number_of_genes);
  printf( "# VTR                          = %f\n", BinaryGOMEA_vtr);
  printf( "# Maximum number of IMS-GOMEAs = %d\n",  BinaryGOMEA_maximum_number_of_BinaryGOMEAs);
  printf( "# Random seed                  = %lld\n", BinaryGOMEA_random_seed);
  printf( "#\n" );
  printf( "###################################################\n" );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Random Numbers -=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns a random double, distributed uniformly between 0 and 1.
 */
double BinaryGOMEA_RandomRealUniform01()
{
  return (double)rand() / (double)RAND_MAX;
}
        
/**
 * Returns a random integer, distributed uniformly between 0 and maximum.
 */
int BinaryGOMEA_RandomInt( int maximum )
{
  return rand() % maximum;
}

/**
 * Returns a random compact (using integers 0,1,...,n-1) permutation
 * of length n using the Fisher-Yates shuffle.
 */
int *BinaryGOMEA_RandomPermutation( int n )
{
  int i, j, dummy, *result;

  result = (int *) BinaryGOMEA_Malloc( n*sizeof( int ) );
  for( i = 0; i < n; i++ )
    result[i] = i;

  for( i = n-1; i > 0; i-- )
  {
    j         = BinaryGOMEA_RandomInt( i+1 );
    dummy     = result[j];
    result[j] = result[i];
    result[i] = dummy;
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns the name of an installed problem.
 */
char *BinaryGOMEA_InstalledProblemName()
{
  switch(BinaryGOMEA_problem_index)
  {
    case  1: return((char *) "Concatenated Deceptive Trap");
    case  2: return((char *) "ADF");
    case  3: return((char *) "MaxCut");
    case  4: return((char *) "MaxSat");
    case  5: return((char *) "Walsh decomposition optimization");            
  }

  return(NULL);
}

/**
 * Returns the number of problems installed.
 */
int BinaryGOMEA_NumberOfInstalledProblems()
{
         int backup;
  static int result = -1;
  
  if( result == -1 )
  {
    backup                    = BinaryGOMEA_problem_index;
    BinaryGOMEA_problem_index = 0;

    while( BinaryGOMEA_InstalledProblemName() != NULL )
      BinaryGOMEA_problem_index++;

    result                    = BinaryGOMEA_problem_index;
    BinaryGOMEA_problem_index = backup;
  }
  
  return( result );
}

/**
 * Evaluates the fitness of a solution.
 */
void BinaryGOMEA_InstalledProblemEvaluationSpecificBinaryGOMEASpecificOffspringIndex( int BinaryGOMEA_index, int offspring_index, int number_of_touched_genes, int *touched_genes_indices, char *genes_before, double fitness_before )
{
  if (GetMilliSecondsRunning(BinaryGOMEA_timestamp_start) > timelimitMilliseconds)
    exit(0);

  double elitist_fitness_before_evaluation = BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index]; // Making a copy of the elitist fitness to be safe against the case where (BinaryGOMEA_index == BinaryGOMEA_elitist_solution_BinaryGOMEA_index) && (offspring_index == BinaryGOMEA_elitist_solution_offspring_index)

  /* Do the actual evaluation */
  
  bool is_found = false;
  double found_value;
  std::vector<bool>solution(BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index], BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index] + BinaryGOMEA_number_of_genes);

  if (!SKIP_ALREADY_EVALUATED)
  {
    archiveRecord record = checkAlreadyEvaluated(solution, evaluated_solutions);
    is_found = record.found;
    found_value = record.value;    
  }

  if (is_found)
  {
    BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] = found_value;
  }
  else
  { 
    double objective_value;
    if (BinaryGOMEA_use_partial_evaluations && genes_before != NULL)
    {
      double newEvals;
      std::vector<bool>solution_backup(genes_before, genes_before + BinaryGOMEA_number_of_genes);
      std::vector<int> touched_genes(touched_genes_indices, touched_genes_indices + number_of_touched_genes);
      objective_value = problem_instance->calculateFitnessPartialEvaluations(solution, touched_genes, solution_backup, fitness_before, newEvals);
      BinaryGOMEA_number_of_evaluations += newEvals;
    }
    else
    {
      objective_value = problem_instance->calculateFitness(solution);
      BinaryGOMEA_number_of_evaluations++;
    }
    
    BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] = objective_value;
    
    if (!SKIP_ALREADY_EVALUATED)
      addToEvaluated(solution, BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index], evaluated_solutions);
   }

  /* Update elitist solution */
  if( (BinaryGOMEA_first_evaluation_ever == 1) || (BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] > elitist_fitness_before_evaluation) )
  {
    BinaryGOMEA_elitist_solution_BinaryGOMEA_index                                 = BinaryGOMEA_index;
    BinaryGOMEA_elitist_solution_offspring_index                                   = offspring_index;
    BinaryGOMEA_elitist_solution_hitting_time_in_milliseconds                      = GetMilliSecondsRunning(BinaryGOMEA_timestamp_start);
    BinaryGOMEA_elitist_solution_hitting_time_in_evaluations                       = BinaryGOMEA_number_of_evaluations;
    //printf("%d\n", BinaryGOMEA_elitist_solution_hitting_time_in_evaluations);
    /* Check the VTR */
    if( BinaryGOMEA_use_VTR > 0 )
    {
      if( BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] >= BinaryGOMEA_vtr )
      {
        BinaryGOMEA_WriteElitistSolutionToFile();
        printf("VTR HIT! %lf %lf\n", BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index], BinaryGOMEA_vtr);
        exit( 0 );
      }
    }

    if( BinaryGOMEA_report_every_improvement )
      BinaryGOMEA_WriteElitistSolutionToFile();
  }

  BinaryGOMEA_first_evaluation_ever = 0;
}


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initialization for a single BinaryGOMEA.
 */
void BinaryGOMEA_InitializeNewBinaryGOMEA()
{
  if( BinaryGOMEA_number_of_BinaryGOMEAs == 0 )
  {
    BinaryGOMEA_populations             = (char ***) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( char ** ) );
    BinaryGOMEA_population_sizes        = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( int ) );
    BinaryGOMEA_fitnesses               = (double **) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( double * ) );
    BinaryGOMEA_offsprings              = (char ***) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( char ** ) );
    BinaryGOMEA_fitnesses_offsprings    = (double **) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( double * ) );
    BinaryGOMEA_average_fitnesses       = (double *) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( double ) );
    BinaryGOMEA_no_improvement_stretchs = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( int * ) );
    BinaryGOMEA_number_of_generations   = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( int ) );
    BinaryGOMEA_terminated              = (char *) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( char ) );
    BinaryGOMEA_FOSs                    = (int ***) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( int ** ) );
    BinaryGOMEA_FOSs_number_of_indices  = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( int * ) );
    BinaryGOMEA_FOSs_length             = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( int ) );

    if( !BinaryGOMEA_use_predetermined_fos )
      BinaryGOMEA_MI_matrices = (double ***) BinaryGOMEA_Malloc( BinaryGOMEA_maximum_number_of_BinaryGOMEAs*sizeof( double ** ) );

    BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs] = BinaryGOMEA_base_population_size;
  }
  else
    BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs] = 2*BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs-1];

  BinaryGOMEA_InitializeNewBinaryGOMEAMemory();

  BinaryGOMEA_InitializeNewBinaryGOMEAPopulationAndFitnessValues();

  BinaryGOMEA_number_of_BinaryGOMEAs++;
}

/**
 * Initializes the memory for a single BinaryGOMEA.
 */
void BinaryGOMEA_InitializeNewBinaryGOMEAMemory()
{
  int i;

  BinaryGOMEA_populations[BinaryGOMEA_number_of_BinaryGOMEAs]             = (char **) BinaryGOMEA_Malloc( BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]*sizeof( char * ) );
  BinaryGOMEA_fitnesses[BinaryGOMEA_number_of_BinaryGOMEAs]               = (double *) BinaryGOMEA_Malloc( BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]*sizeof( double ) );
  BinaryGOMEA_offsprings[BinaryGOMEA_number_of_BinaryGOMEAs]              = (char **) BinaryGOMEA_Malloc( BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]*sizeof( char * ) );
  BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_number_of_BinaryGOMEAs]    = (double *) BinaryGOMEA_Malloc( BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]*sizeof( double ) );
  BinaryGOMEA_no_improvement_stretchs[BinaryGOMEA_number_of_BinaryGOMEAs] = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]*sizeof( int ) );

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]; i++ )
    BinaryGOMEA_populations[BinaryGOMEA_number_of_BinaryGOMEAs][i] = (char *) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( char ) );

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]; i++ )
    BinaryGOMEA_no_improvement_stretchs[BinaryGOMEA_number_of_BinaryGOMEAs][i] = 0;

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]; i++ )
    BinaryGOMEA_offsprings[BinaryGOMEA_number_of_BinaryGOMEAs][i] = (char *) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( char ) );

  if( !BinaryGOMEA_use_predetermined_fos )
  {
    BinaryGOMEA_MI_matrices[BinaryGOMEA_number_of_BinaryGOMEAs] = (double **) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( double * ) );
    for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
      (BinaryGOMEA_MI_matrices[BinaryGOMEA_number_of_BinaryGOMEAs])[i] = (double *) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( double ) );
  }
}

/**
 * Initializes the population and the BinaryGOMEA_fitnesses by randomly
 * generation n solutions.
 */
void BinaryGOMEA_InitializeNewBinaryGOMEAPopulationAndFitnessValues()
{
  int    i, j, number_of_0_to_generate, number_of_1_to_generate;
  double probability_of_generating_0;

  BinaryGOMEA_terminated[BinaryGOMEA_number_of_BinaryGOMEAs]            = 0;
  BinaryGOMEA_number_of_generations[BinaryGOMEA_number_of_BinaryGOMEAs] = 0;
  BinaryGOMEA_FOSs[BinaryGOMEA_number_of_BinaryGOMEAs]                  = NULL;
  if( BinaryGOMEA_predetermined_FOS != NULL )
  {
    BinaryGOMEA_FOSs_length[BinaryGOMEA_number_of_BinaryGOMEAs]            = BinaryGOMEA_predetermined_FOS_length;
    BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_number_of_BinaryGOMEAs] = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_predetermined_FOS_length*sizeof( int ) );
    BinaryGOMEA_FOSs[BinaryGOMEA_number_of_BinaryGOMEAs]                   = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_predetermined_FOS_length*sizeof( int * ) );
    for( i = 0; i < BinaryGOMEA_predetermined_FOS_length; i++ )
    {
      BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_number_of_BinaryGOMEAs][i] = BinaryGOMEA_predetermined_FOS_number_of_indices[i];
      BinaryGOMEA_FOSs[BinaryGOMEA_number_of_BinaryGOMEAs][i]                   = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_predetermined_FOS_number_of_indices[i]*sizeof( int ) );
      for( j = 0; j < BinaryGOMEA_predetermined_FOS_number_of_indices[i]; j++ )
        BinaryGOMEA_FOSs[BinaryGOMEA_number_of_BinaryGOMEAs][i][j] = BinaryGOMEA_predetermined_FOS[i][j];
    }
  }

  if( BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs] == 1 )
  {
    for( j = 0; j < BinaryGOMEA_number_of_genes; j++ )
      BinaryGOMEA_offsprings[BinaryGOMEA_number_of_BinaryGOMEAs][0][j] = BinaryGOMEA_RandomInt( 2 );
  }
  else
  {
    for( j = 0; j < BinaryGOMEA_number_of_genes; j++ )
    {
      number_of_0_to_generate = BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]/2;
      number_of_1_to_generate = BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs] - number_of_0_to_generate;

      for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]; i++ )
      {
        probability_of_generating_0 = ((double) number_of_0_to_generate)/((double) (number_of_0_to_generate+number_of_1_to_generate)); 

        BinaryGOMEA_offsprings[BinaryGOMEA_number_of_BinaryGOMEAs][i][j] = (BinaryGOMEA_RandomRealUniform01() < probability_of_generating_0 ? 0 : 1); 
      
        if( BinaryGOMEA_offsprings[BinaryGOMEA_number_of_BinaryGOMEAs][i][j] == 0 )
          number_of_0_to_generate--;
        else
          number_of_1_to_generate--;
      }
    }
  }

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_number_of_BinaryGOMEAs]; i++ )
    BinaryGOMEA_InstalledProblemEvaluationSpecificBinaryGOMEASpecificOffspringIndex( BinaryGOMEA_number_of_BinaryGOMEAs, i, 0, NULL, NULL, 0 );
  
  BinaryGOMEA_CopyOffspringToPopulationSpecificBinaryGOMEA( BinaryGOMEA_number_of_BinaryGOMEAs );

  BinaryGOMEA_computeAverageFitnessSpecificBinaryGOMEA( BinaryGOMEA_number_of_BinaryGOMEAs );
}

/**
 * Checks to see if a file exists with a predetermined FOS in it.
 */
void BinaryGOMEA_InitializePredeterminedFOS()
{
  char  c, string[100000], filename[1000];
  int   i, j, k;
  FILE *file;

  /* Initialize FOS to NULL */
  BinaryGOMEA_predetermined_FOS = NULL;

  /* Do not go any further if the user did not specificy they wanted to use a predetermined FOS. */
  if( !BinaryGOMEA_use_predetermined_fos )
    return;
  
  /* Open file for reading */
  sprintf( filename, "predetermined_fos.txt" );
  file = fopen( filename, "r" );
  if( file == NULL )
  {
    printf( "\n" );
    printf( "Error: option to use predetermined FOS was specified (-F), but no file \"predetermined_fos.txt\" could be found.\n" );
    printf( "\n\n" );
    exit( 0 );
  }

  /* Count number of FOS elements */
  BinaryGOMEA_predetermined_FOS_length = 0;
  c                              = fgetc( file );
  while( c != EOF )
  {
    if( c == ']' )
      BinaryGOMEA_predetermined_FOS_length++;

    c = fgetc( file );
  }
  printf("%d\n",BinaryGOMEA_predetermined_FOS_length);

  /* Count number of indices in each FOS element */
  rewind( file );
  BinaryGOMEA_predetermined_FOS_number_of_indices = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_predetermined_FOS_length*sizeof( int ) );
  for( i = 0; i < BinaryGOMEA_predetermined_FOS_length; i++ )
    BinaryGOMEA_predetermined_FOS_number_of_indices[i] = 0;
  
  c = fgetc( file );
  i = 0;
  while( c != EOF )
  {
    if( c == ']' )
    {
      BinaryGOMEA_predetermined_FOS_number_of_indices[i]++;
      i++;
    }
    else if( c == ' ' )
    {
      BinaryGOMEA_predetermined_FOS_number_of_indices[i]++;
    }

    c = fgetc( file );
  }

  /* Read indices in FOS elements */
  rewind( file );
  BinaryGOMEA_predetermined_FOS = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_predetermined_FOS_length*sizeof( int * ) );
   for( i = 0; i < BinaryGOMEA_predetermined_FOS_length; i++ )
     BinaryGOMEA_predetermined_FOS[i] = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_predetermined_FOS_number_of_indices[i]*sizeof( int ) );

  c = fgetc( file );
  i = 0;
  j = 0;
  k = 0;
  while( c != EOF )
  {
    if( c == '[' )
    {
      j = 0;
    }
    else if( c == ']' )
    {
      string[k]                           = '\0';
      BinaryGOMEA_predetermined_FOS[i][j] = atof( string );
      printf("%d %d | %d\n", i,j,BinaryGOMEA_predetermined_FOS[i][j]);
      i++;
      j = 0;
      k = 0;
    }
    else if( c == ' ' )
    {
      string[k]                           = '\0';
      BinaryGOMEA_predetermined_FOS[i][j] = atof( string );

      j++;
      k = 0;
    }
    else if( (c >= '0') && (c <= '9') )
    {
      string[k] = c;
      k++;
    }      

    c = fgetc( file );
  }

  /* Close file. */
  printf("FOS read\n");
  for (int i = 0; i < BinaryGOMEA_predetermined_FOS_length; ++i)
  {
    for (int j = 0; j <  BinaryGOMEA_predetermined_FOS_number_of_indices[i]; ++j)
    {
      printf("%d ", BinaryGOMEA_predetermined_FOS[i][j]);
    }
    printf("\n");
  }
  fclose( file );
}

/**
 * Initializes the pseudo-random number generator.
 */
void BinaryGOMEA_InitializeRandomNumberGenerator()
{
  srand(BinaryGOMEA_random_seed);
}

void BinaryGOMEA_InitializeProblem()
{
  createProblemInstance(BinaryGOMEA_problem_index, BinaryGOMEA_number_of_genes, &problem_instance, instance_path, k, s);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*=-=-=-=-=-=-=-=-=-=-= Section Survivor Selection =-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Determines the solutions that finally survive the generation (offspring only).
 */
void BinaryGOMEA_CopyOffspringToPopulationSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int i, j;

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
  {
    for( j = 0; j < BinaryGOMEA_number_of_genes; j++ )
      BinaryGOMEA_populations[BinaryGOMEA_index][i][j] = BinaryGOMEA_offsprings[BinaryGOMEA_index][i][j];
    BinaryGOMEA_fitnesses[BinaryGOMEA_index][i] = BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][i];
  }
}


void BinaryGOMEA_computeAverageFitnessSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int i;

  BinaryGOMEA_average_fitnesses[BinaryGOMEA_index]  = 0;
  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
    BinaryGOMEA_average_fitnesses[BinaryGOMEA_index]  += BinaryGOMEA_fitnesses[BinaryGOMEA_index][i];
  BinaryGOMEA_average_fitnesses[BinaryGOMEA_index] /= (double) (BinaryGOMEA_population_sizes[BinaryGOMEA_index]);
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
void BinaryGOMEA_WriteElitistSolutionToFile()
{
  char       string[10000], filename1[1000], filename2[1000], filename3[1000];
  int        i;
  FILE      *file;

  sprintf(filename1, "%s/elitist.dat", folder);
  BinaryGOMEA_elitist_solution_last_written_running_time = GetMilliSecondsRunning(BinaryGOMEA_timestamp_start);

  if( !BinaryGOMEA_elitist_solution_written_before )
  {
    file = fopen( filename1, "w" );
    BinaryGOMEA_elitist_solution_written_before = 1;
  }
  else
    file = fopen( filename1, "a" );
  /* Write actual values. */
  sprintf( string, "%d %d %.6lf\n", (int)BinaryGOMEA_elitist_solution_hitting_time_in_evaluations, BinaryGOMEA_elitist_solution_hitting_time_in_milliseconds, BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index] );
  fputs( string, file );
  fclose( file );

  sprintf(filename2, "%s/elitistonly.dat", folder);
  file = fopen( filename2, "w" );
  sprintf( string, "%.6lf", BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index] );
  fputs( string, file );
  fclose( file );

  sprintf(filename3, "%s/elitistsolution.dat", folder);
  file = fopen( filename3, "w" );
  for (int i = 0; i < BinaryGOMEA_number_of_genes; ++i)
  {
    sprintf( string, "%d", BinaryGOMEA_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index][i] );
    fputs( string, file );
  }
  fclose( file );

}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced, 0 otherwise.
 */
char BinaryGOMEA_CheckTermination()
{
  int i;

  if( BinaryGOMEA_number_of_BinaryGOMEAs == BinaryGOMEA_maximum_number_of_BinaryGOMEAs )
  {
    for( i = 0; i < BinaryGOMEA_maximum_number_of_BinaryGOMEAs; i++ )
    {
      if( !BinaryGOMEA_terminated[i] )
        return( 0 );
    }

    return( 1 );
  }
  
  return( 0 );
}

char BinaryGOMEA_CheckTerminationSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int i, j;

  for( i = BinaryGOMEA_index+1; i < BinaryGOMEA_number_of_BinaryGOMEAs; i++ )
  {
    if( BinaryGOMEA_average_fitnesses[i] > BinaryGOMEA_average_fitnesses[BinaryGOMEA_index] )
    {
      BinaryGOMEA_minimum_BinaryGOMEA_index = BinaryGOMEA_index+1;

      return( 1 );
    }
  }

  for( i = 1; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
  {
    for( j = 0; j < BinaryGOMEA_number_of_genes; j++ )
    {
      if( BinaryGOMEA_populations[BinaryGOMEA_index][i][j] != BinaryGOMEA_populations[BinaryGOMEA_index][0][j] )
        return( 0 );
    }
  }

  return( 1 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Variation -==-=-=-=-=-=-=-=-=-=-=-=-=-=*/
void BinaryGOMEA_GenerationalStepAllBinaryGOMEAsRecursiveFold( int BinaryGOMEA_index_smallest, int BinaryGOMEA_index_biggest );
void BinaryGOMEA_GenerationalStepAllBinaryGOMEAs()
{
  int BinaryGOMEA_index_smallest, BinaryGOMEA_index_biggest;

  BinaryGOMEA_index_biggest  = BinaryGOMEA_number_of_BinaryGOMEAs-1;
  BinaryGOMEA_index_smallest = 0;
  while( BinaryGOMEA_index_smallest <= BinaryGOMEA_index_biggest )
  {
    if( !BinaryGOMEA_terminated[BinaryGOMEA_index_smallest] )
      break;

    BinaryGOMEA_index_smallest++;
  }

  BinaryGOMEA_GenerationalStepAllBinaryGOMEAsRecursiveFold( BinaryGOMEA_index_smallest, BinaryGOMEA_index_biggest );
}

void BinaryGOMEA_GenerationalStepAllBinaryGOMEAsRecursiveFold( int BinaryGOMEA_index_smallest, int BinaryGOMEA_index_biggest )
{
  int i, BinaryGOMEA_index;

  for( i = 0; i < BinaryGOMEA_IMS_subgeneration_factor-1; i++ )
  {
    for( BinaryGOMEA_index = BinaryGOMEA_index_smallest; BinaryGOMEA_index <= BinaryGOMEA_index_biggest; BinaryGOMEA_index++ )
    {
      if( !BinaryGOMEA_terminated[BinaryGOMEA_index] )
      {
        BinaryGOMEA_terminated[BinaryGOMEA_index] = BinaryGOMEA_CheckTerminationSpecificBinaryGOMEA( BinaryGOMEA_index );
      }

      if( (!BinaryGOMEA_terminated[BinaryGOMEA_index]) && (BinaryGOMEA_index >= BinaryGOMEA_minimum_BinaryGOMEA_index) )
      {
        BinaryGOMEA_MakeOffspringSpecificBinaryGOMEA( BinaryGOMEA_index );

        BinaryGOMEA_CopyOffspringToPopulationSpecificBinaryGOMEA( BinaryGOMEA_index );

        BinaryGOMEA_computeAverageFitnessSpecificBinaryGOMEA( BinaryGOMEA_index );

        BinaryGOMEA_number_of_generations[BinaryGOMEA_index]++;
      }
    }

    for( BinaryGOMEA_index = BinaryGOMEA_index_smallest; BinaryGOMEA_index < BinaryGOMEA_index_biggest; BinaryGOMEA_index++ )
      BinaryGOMEA_GenerationalStepAllBinaryGOMEAsRecursiveFold( BinaryGOMEA_index_smallest, BinaryGOMEA_index );
  }
}

void BinaryGOMEA_MakeOffspringSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
 if( BinaryGOMEA_predetermined_FOS == NULL )
    BinaryGOMEA_LearnLTFOSSpecificBinaryGOMEA( BinaryGOMEA_index );

  if( BinaryGOMEA_write_FOS_contents_to_file )
    BinaryGOMEA_WriteToFileFOSContentsSpecificBinaryGOMEA( BinaryGOMEA_index );
  
  if( BinaryGOMEA_write_MI_matrix_to_file )
    BinaryGOMEA_WriteToFileMIMatrixSpecificBinaryGOMEA( BinaryGOMEA_index );

  BinaryGOMEA_GenerateAndEvaluateOffspringSpecificBinaryGOMEA( BinaryGOMEA_index );
}

/**
 * Learns a linkage tree FOS by means of hierarchical clustering.
 * This implementation follows the reciprocal nearest neighbor approach.
 */
void BinaryGOMEA_LearnLTFOSSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  char     done;
  int      i, j, r0, r1, rswap, *indices, *order,
           BinaryGOMEA_FOSs_index, **mpm, *mpm_number_of_indices, mpm_length,
         **mpm_new, *mpm_new_number_of_indices, mpm_new_length,
          *NN_chain, NN_chain_length, *BinaryGOMEA_FOSs_index_of_mpm_element;
  double **S_matrix, mul0, mul1;

  if( BinaryGOMEA_FOSs[BinaryGOMEA_index] != NULL )
  {
    for( i = 0; i < BinaryGOMEA_FOSs_length[BinaryGOMEA_index]; i++ )
      free( BinaryGOMEA_FOSs[BinaryGOMEA_index][i] );
    free( BinaryGOMEA_FOSs[BinaryGOMEA_index] );
    free( BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index] );
  }

  BinaryGOMEA_FOSs_index_of_mpm_element = NULL; /* Only needed to prevent compiler warnings. */

  /* Compute Mutual Information matrix */
  BinaryGOMEA_ComputeMIMatrixSpecificBinaryGOMEA( BinaryGOMEA_index );

  /* Initialize MPM to the univariate factorization */
  order                 = BinaryGOMEA_RandomPermutation( BinaryGOMEA_number_of_genes );
  mpm                   = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( int * ) );
  mpm_number_of_indices = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( int ) );
  mpm_length            = BinaryGOMEA_number_of_genes;
  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
  {
    if (BinaryGOMEA_predetermined_FOS == NULL)
    {
      indices                  = (int *) BinaryGOMEA_Malloc( 1*sizeof( int ) );
      indices[0]               = order[i];
      mpm[i]                   = indices;
      mpm_number_of_indices[i] = 1;
    }    
  }
  free( order );

  /* Initialize LT to the initial MPM */
  BinaryGOMEA_FOSs_length[BinaryGOMEA_index]            = BinaryGOMEA_number_of_genes+BinaryGOMEA_number_of_genes-1;
  BinaryGOMEA_FOSs[BinaryGOMEA_index]                   = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_FOSs_length[BinaryGOMEA_index]*sizeof( int * ) );
  BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index] = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_FOSs_length[BinaryGOMEA_index]*sizeof( int ) );
  BinaryGOMEA_FOSs_index                                = 0;
  for( i = 0; i < mpm_length; i++ )
  {
    BinaryGOMEA_FOSs[BinaryGOMEA_index][BinaryGOMEA_FOSs_index]                   = mpm[i];
    BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][BinaryGOMEA_FOSs_index] = mpm_number_of_indices[i];
    BinaryGOMEA_FOSs_index++;
  }

  /* Initialize similarity matrix */
  S_matrix = (double **) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( double * ) );
  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
    S_matrix[i] = (double *) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( double ) );
  for( i = 0; i < mpm_length; i++ )
    for( j = 0; j < mpm_length; j++ )
      S_matrix[i][j] = BinaryGOMEA_MI_matrices[BinaryGOMEA_index][mpm[i][0]][mpm[j][0]];
  for( i = 0; i < mpm_length; i++ )
    S_matrix[i][i] = 0;

  NN_chain        = (int *) BinaryGOMEA_Malloc( (BinaryGOMEA_number_of_genes+2)*sizeof( int ) );
  NN_chain_length = 0;
  done            = 0;
  while( !done )
  {
    if( NN_chain_length == 0 )
    {
      NN_chain[NN_chain_length] = BinaryGOMEA_RandomInt( mpm_length );
      NN_chain_length++;
    }

    while( NN_chain_length < 3 )
    {
      NN_chain[NN_chain_length] = BinaryGOMEA_DetermineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      NN_chain_length++;
    }

    while( NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1] )
    {
      NN_chain[NN_chain_length] = BinaryGOMEA_DetermineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      if( ((S_matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length]] == S_matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length-2]])) && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
        NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];
      NN_chain_length++;
      if( NN_chain_length > BinaryGOMEA_number_of_genes )
        break;
    }
    r0 = NN_chain[NN_chain_length-2];
    r1 = NN_chain[NN_chain_length-1];
    if( r0 > r1 )
    {
      rswap = r0;
      r0    = r1;
      r1    = rswap;
    }
    NN_chain_length -= 3;

    if( r1 < mpm_length ) /* This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain */
    {
      indices = (int *) BinaryGOMEA_Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );
  
      i = 0;
      for( j = 0; j < mpm_number_of_indices[r0]; j++ )
      {
        indices[i] = mpm[r0][j];
        i++;
      }
      for( j = 0; j < mpm_number_of_indices[r1]; j++ )
      {
        indices[i] = mpm[r1][j];
        i++;
      }
    
      BinaryGOMEA_FOSs[BinaryGOMEA_index][BinaryGOMEA_FOSs_index]                   = indices;
      BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][BinaryGOMEA_FOSs_index] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      BinaryGOMEA_FOSs_index++;
  
      mul0 = ((double) mpm_number_of_indices[r0])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      mul1 = ((double) mpm_number_of_indices[r1])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      for( i = 0; i < mpm_length; i++ )
      {
        if( (i != r0) && (i != r1) )
        {
          S_matrix[i][r0] = mul0*S_matrix[i][r0] + mul1*S_matrix[i][r1];
          S_matrix[r0][i] = S_matrix[i][r0];
        }
      }
  
      mpm_new                   = (int **) BinaryGOMEA_Malloc( (mpm_length-1)*sizeof( int * ) );
      mpm_new_number_of_indices = (int *) BinaryGOMEA_Malloc( (mpm_length-1)*sizeof( int ) );
      mpm_new_length            = mpm_length-1;
      for( i = 0; i < mpm_new_length; i++ )
      {
        mpm_new[i]                   = mpm[i];
        mpm_new_number_of_indices[i] = mpm_number_of_indices[i];
      }
  
      mpm_new[r0]                   = indices;
      mpm_new_number_of_indices[r0] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      if( r1 < mpm_length-1 )
      {
        mpm_new[r1]                   = mpm[mpm_length-1];
        mpm_new_number_of_indices[r1] = mpm_number_of_indices[mpm_length-1];
  
        for( i = 0; i < r1; i++ )
        {
          S_matrix[i][r1] = S_matrix[i][mpm_length-1];
          S_matrix[r1][i] = S_matrix[i][r1];
        }
  
        for( j = r1+1; j < mpm_new_length; j++ )
        {
          S_matrix[r1][j] = S_matrix[j][mpm_length-1];
          S_matrix[j][r1] = S_matrix[r1][j];
        }
      }
  
      for( i = 0; i < NN_chain_length; i++ )
      {
        if( NN_chain[i] == mpm_length-1 )
        {
          NN_chain[i] = r1;
          break;
        }
      }
  
      free( mpm );
      free( mpm_number_of_indices );
      mpm                   = mpm_new;
      mpm_number_of_indices = mpm_new_number_of_indices;
      mpm_length            = mpm_new_length;
  
      if( mpm_length == 1 )
        done = 1;
    }
  }

  free( NN_chain );

  free( mpm_new );
  free( mpm_number_of_indices );

  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
    free( S_matrix[i] );
  free( S_matrix );

  free( BinaryGOMEA_FOSs_index_of_mpm_element );
}

/**
 * Determines nearest neighbour according to similarity values.
 */
int BinaryGOMEA_DetermineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length )
{
  int i, result;

  result = 0;
  if( result == index )
    result++;
  for( i = 1; i < mpm_length; i++ )
  {
    if( ((S_matrix[index][i] > S_matrix[index][result]) || ((S_matrix[index][i] == S_matrix[index][result]) && (mpm_number_of_indices[i] < mpm_number_of_indices[result]))) && (i != index) )
      result = i;
  }

  return( result );
}

void BinaryGOMEA_ComputeMIMatrixSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int    i, j, k, *indices, factor_size;
  double p, *cumulative_probabilities;

  /* Compute joint entropy matrix */
  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
  {
    for( j = i+1; j < BinaryGOMEA_number_of_genes; j++ )
    {
      indices                  = (int *) BinaryGOMEA_Malloc( 2*sizeof( int ) );
      indices[0]               = i;
      indices[1]               = j;
      cumulative_probabilities = BinaryGOMEA_EstimateParametersForSingleBinaryMarginalSpecificBinaryGOMEA( BinaryGOMEA_index, indices, 2, &factor_size );

      BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j] = 0.0;
      for( k = 0; k < factor_size; k++ )
      {
        if( k == 0 )
          p = cumulative_probabilities[k];
        else
          p = cumulative_probabilities[k]-cumulative_probabilities[k-1];
        if( p > 0 )
          BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j] += -p*BinaryGOMEA_Log2(p);
      }

      BinaryGOMEA_MI_matrices[BinaryGOMEA_index][j][i] = BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j];

      free( indices );
      free( cumulative_probabilities );
    }
    indices                  = (int *) BinaryGOMEA_Malloc( 1*sizeof( int ) );
    indices[0]               = i;
    cumulative_probabilities = BinaryGOMEA_EstimateParametersForSingleBinaryMarginalSpecificBinaryGOMEA( BinaryGOMEA_index, indices, 1, &factor_size );

    BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][i] = 0.0;
    for( k = 0; k < factor_size; k++ )
    {
      if( k == 0 )
        p = cumulative_probabilities[k];
      else
        p = cumulative_probabilities[k]-cumulative_probabilities[k-1];
      if( p > 0 )
       BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][i] += -p*BinaryGOMEA_Log2(p);
    }

    free( indices );
    free( cumulative_probabilities );
  }

  /* Then transform into mutual information matrix MI(X,Y)=H(X)+H(Y)-H(X,Y) */
  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
    for( j = i+1; j < BinaryGOMEA_number_of_genes; j++ )
    {
      BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j] = BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][i] + BinaryGOMEA_MI_matrices[BinaryGOMEA_index][j][j] - BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j];
      BinaryGOMEA_MI_matrices[BinaryGOMEA_index][j][i] = BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j];
    }
}

/**
 * Estimates the cumulative probability distribution of a
 * single binary marginal.
 */
double *BinaryGOMEA_EstimateParametersForSingleBinaryMarginalSpecificBinaryGOMEA( int BinaryGOMEA_index, int *indices, int number_of_indices, int *factor_size )
{
  int     i, j, index, power_of_two;
  double *result;

  *factor_size = (int) pow( 2, number_of_indices );
  result       = (double *) BinaryGOMEA_Malloc( (*factor_size)*sizeof( double ) );

  for( i = 0; i < (*factor_size); i++ )
    result[i] = 0.0;

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
  {
    index        = 0;
    power_of_two = 1;
    for( j = number_of_indices-1; j >= 0; j-- )
    {
      index += BinaryGOMEA_populations[BinaryGOMEA_index][i][indices[j]] ? power_of_two : 0;
      power_of_two *= 2;
    }

    result[index] += 1.0;
  }

  for( i = 0; i < (*factor_size); i++ )
    result[i] /= (double) BinaryGOMEA_population_sizes[BinaryGOMEA_index];

  for( i = 1; i < (*factor_size); i++ )
    result[i] += result[i-1];

  result[(*factor_size)-1] = 1.0;

  return( result );
}

void BinaryGOMEA_WriteToFileFOSContentsSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int   i, j;
  char  string[1000];
  FILE *file;

  sprintf( string, "FOS_BinaryGOMEA_%02d_generation_%02d.dat", BinaryGOMEA_index, BinaryGOMEA_number_of_generations[BinaryGOMEA_index]  );
  file = fopen( string, "w" );

  for( i = 0; i < BinaryGOMEA_FOSs_length[BinaryGOMEA_index]; i++ )
  {
    sprintf( string, "[" );
    fputs( string, file );
    for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
    {
      sprintf( string, "%d", BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j] );
      fputs( string, file );
      if( j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]-1 )
      {
        sprintf( string, " " );
        fputs( string, file );
      }
    }
    sprintf( string, "]\n" );
    fputs( string, file );
  }
  
  fclose( file );
}

void BinaryGOMEA_WriteToFileMIMatrixSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int   i, j;
  char  string[1000];
  FILE *file;

  if( BinaryGOMEA_use_predetermined_fos )
    return;

  sprintf( string, "MI_Matrix_BinaryGOMEA_%02d_generation_%02d.dat", BinaryGOMEA_index, BinaryGOMEA_number_of_generations[BinaryGOMEA_index]  );
  file = fopen( string, "w" );

  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
  {
    fputs( string, file );
    for( j = 0; j < BinaryGOMEA_number_of_genes; j++ )
    {
      sprintf( string, "%37.30e", BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i][j] );
      fputs( string, file );
      if( j < BinaryGOMEA_number_of_genes )
      {
        sprintf( string, " " );
        fputs( string, file );
      }
    }
    sprintf( string, "\n" );
    fputs( string, file );
  }
  
  fclose( file );
}

/**
 * Computes the two-log of x.
 */
double math_log_two = log(2.0);
double BinaryGOMEA_Log2( double x )
{
  return( log(x) / math_log_two );
}

/**
 * Generates new solutions.
 */
void BinaryGOMEA_GenerateAndEvaluateOffspringSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int i;

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
  {
    BinaryGOMEA_GenerateAndEvaluateNewSolutionSpecificBinaryGOMEASpecificOffspringIndex( BinaryGOMEA_index, i );

    if( !(BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][i] > BinaryGOMEA_fitnesses[BinaryGOMEA_index][i]) )
      BinaryGOMEA_no_improvement_stretchs[BinaryGOMEA_index][i]++;
    else
      BinaryGOMEA_no_improvement_stretchs[BinaryGOMEA_index][i] = 0;
  }
}

/**
 * Performs Genepool Optimal Mixing (for one solution in the population).
 */
void BinaryGOMEA_GenerateAndEvaluateNewSolutionSpecificBinaryGOMEASpecificOffspringIndex( int BinaryGOMEA_index, int offspring_index )
{
  char   *backup, donor_genes_are_the_same, this_is_the_elitist_solution;
  short   solution_has_changed;
  int     i, j, donor_index;
  double  obj_backup;

  solution_has_changed = 0;
  this_is_the_elitist_solution = BinaryGOMEA_elitist_solution_BinaryGOMEA_index == BinaryGOMEA_index && BinaryGOMEA_elitist_solution_offspring_index == offspring_index;

  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
    BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][i] = BinaryGOMEA_populations[BinaryGOMEA_index][offspring_index][i];
  BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] = BinaryGOMEA_fitnesses[BinaryGOMEA_index][offspring_index];

  backup = (char *) BinaryGOMEA_Malloc( BinaryGOMEA_number_of_genes*sizeof( char ) );
  for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
    backup[i] = BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][i];

  obj_backup = BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index];

  
  /* Phase 1: optimal mixing with random donors */
  BinaryGOMEA_ShuffleFOSSpecificBinaryGOMEA( BinaryGOMEA_index );

  for( i = 0; i < BinaryGOMEA_FOSs_length[BinaryGOMEA_index]; i++ )
  {
    if( BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i] == BinaryGOMEA_number_of_genes )
      continue;

    do
    {
      donor_index = BinaryGOMEA_RandomInt( BinaryGOMEA_population_sizes[BinaryGOMEA_index] );
    } while( donor_index == offspring_index );

    for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
      BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] = BinaryGOMEA_populations[BinaryGOMEA_index][donor_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]];

    /* Test if the change is for the better */
    donor_genes_are_the_same = 1;
    for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
    {
      if( backup[BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] != BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] )
      {
        donor_genes_are_the_same = 0;
        break;
      }
    }
    if( !donor_genes_are_the_same )
    {
      BinaryGOMEA_InstalledProblemEvaluationSpecificBinaryGOMEASpecificOffspringIndex( BinaryGOMEA_index, offspring_index, BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i], BinaryGOMEA_FOSs[BinaryGOMEA_index][i], backup, obj_backup );
     
      // accept the change if this solution is not the elitist and the fitness is at least equally good (allows random walk in neutral fitness landscape)
      // however, if this is the elitist solution, only accept strict improvements, to avoid convergence problems
      if( (!this_is_the_elitist_solution && BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] >= obj_backup) || 
          (this_is_the_elitist_solution && BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] > obj_backup) )   
      {
        for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
          backup[BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] = BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]];

        obj_backup = BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index];

        
        solution_has_changed = 1;
      }
      else
      {
        for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
          BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] = backup[BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]];

        BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] = obj_backup;

      }
    }
  }

  /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
  if( (!solution_has_changed) || (BinaryGOMEA_no_improvement_stretchs[BinaryGOMEA_index][offspring_index] > (1+(log(BinaryGOMEA_population_sizes[BinaryGOMEA_index])/log(10)))) )
  {
    BinaryGOMEA_ShuffleFOSSpecificBinaryGOMEA( BinaryGOMEA_index );

    solution_has_changed = 0;
    for( i = 0; i < BinaryGOMEA_FOSs_length[BinaryGOMEA_index]; i++ )
    {
      /* Convert elite solution to binary representation and set factor variables. */
      for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
        BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] = BinaryGOMEA_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]];

      /* Test if the change is for the better */
      donor_genes_are_the_same = 1;
      for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
      {
        if( backup[BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] != BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] )
        {
          donor_genes_are_the_same = 0;
          break;
        }
      }
      if( !donor_genes_are_the_same )
      {
        BinaryGOMEA_InstalledProblemEvaluationSpecificBinaryGOMEASpecificOffspringIndex( BinaryGOMEA_index, offspring_index, BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i], BinaryGOMEA_FOSs[BinaryGOMEA_index][i], backup, obj_backup );

        if( BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] > obj_backup )
        {
          for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
            backup[BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] = BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]];

          obj_backup = BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index];
          
          solution_has_changed = 1;
        }
        else
        {
          for( j = 0; j < BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][i]; j++ )
            BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]] = backup[BinaryGOMEA_FOSs[BinaryGOMEA_index][i][j]];

          BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] = obj_backup;

        }
      }
      if( solution_has_changed )
        break;
    }

    if( !solution_has_changed )
    {
      if( BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index] > BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] )
        solution_has_changed = 1;

      for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
        BinaryGOMEA_offsprings[BinaryGOMEA_index][offspring_index][i] = BinaryGOMEA_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index][i];

      BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index][offspring_index] = BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_elitist_solution_BinaryGOMEA_index][BinaryGOMEA_elitist_solution_offspring_index];
    }
  }

  free( backup );
  
}


/**
 * Shuffles the FOS (ordering), but not the contents
 * of the linkage subsets themselves.
 */
void BinaryGOMEA_ShuffleFOSSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int i, *order, **BinaryGOMEA_FOSs_new, *BinaryGOMEA_FOSs_number_of_indices_new;

  BinaryGOMEA_FOSs_new                   = (int **) BinaryGOMEA_Malloc( BinaryGOMEA_FOSs_length[BinaryGOMEA_index]*sizeof( int * ) );
  BinaryGOMEA_FOSs_number_of_indices_new = (int *) BinaryGOMEA_Malloc( BinaryGOMEA_FOSs_length[BinaryGOMEA_index]*sizeof( int ) );
  order                                  = BinaryGOMEA_RandomPermutation( BinaryGOMEA_FOSs_length[BinaryGOMEA_index] );
  for( i = 0; i < BinaryGOMEA_FOSs_length[BinaryGOMEA_index]; i++ )
  {
    BinaryGOMEA_FOSs_new[i]                   = BinaryGOMEA_FOSs[BinaryGOMEA_index][order[i]];
    BinaryGOMEA_FOSs_number_of_indices_new[i] = BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index][order[i]];
  }
  free( BinaryGOMEA_FOSs[BinaryGOMEA_index] );
  free( BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index] );
  BinaryGOMEA_FOSs[BinaryGOMEA_index]                   = BinaryGOMEA_FOSs_new;
  BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index] = BinaryGOMEA_FOSs_number_of_indices_new;

  free( order );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section BinaryGOMEA_Ezilaitini -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Undoes BinaryGOMEA initializations.
 */
void BinaryGOMEA_EzilaitiniAllBinaryGOMEAs()
{
  int i;

  for( i = 0; i < BinaryGOMEA_number_of_BinaryGOMEAs; i++ )
    BinaryGOMEA_EzilaitiniSpecificBinaryGOMEA( i );

  free( BinaryGOMEA_FOSs_length );
  free( BinaryGOMEA_FOSs_number_of_indices );
  free( BinaryGOMEA_FOSs );
  free( BinaryGOMEA_populations );
  free( BinaryGOMEA_fitnesses );
  free( BinaryGOMEA_average_fitnesses );
  free( BinaryGOMEA_terminated );
  free( BinaryGOMEA_offsprings );
  free( BinaryGOMEA_fitnesses_offsprings );
  free( BinaryGOMEA_no_improvement_stretchs );
  free( BinaryGOMEA_population_sizes );
  if( BinaryGOMEA_predetermined_FOS != NULL )
  {
    free( BinaryGOMEA_predetermined_FOS_number_of_indices );
    for( i = 0; i < BinaryGOMEA_predetermined_FOS_length; i++ )
      free( BinaryGOMEA_predetermined_FOS[i] );
    free( BinaryGOMEA_predetermined_FOS );
  }
  if( !BinaryGOMEA_use_predetermined_fos )
    free( BinaryGOMEA_MI_matrices );
}


void BinaryGOMEA_EzilaitiniSpecificBinaryGOMEA( int BinaryGOMEA_index )
{
  int i;

  if( BinaryGOMEA_FOSs[BinaryGOMEA_index] != NULL )
  {
    for( i = 0; i < BinaryGOMEA_FOSs_length[BinaryGOMEA_index]; i++ )
      free( BinaryGOMEA_FOSs[BinaryGOMEA_index][i] );
    free( BinaryGOMEA_FOSs[BinaryGOMEA_index] );
    free( BinaryGOMEA_FOSs_number_of_indices[BinaryGOMEA_index] );
  }

  if( !BinaryGOMEA_use_predetermined_fos )
  {
    for( i = 0; i < BinaryGOMEA_number_of_genes; i++ )
      free( BinaryGOMEA_MI_matrices[BinaryGOMEA_index][i] );
    free( BinaryGOMEA_MI_matrices[BinaryGOMEA_index] );
  }

  BinaryGOMEA_EzilaitiniSpecificBinaryGOMEAMemoryForPopulationAndOffspring( BinaryGOMEA_index );
}

/**
 * Initializes the memory for a single BinaryGOMEA instance, for the population only.
 */
void BinaryGOMEA_EzilaitiniSpecificBinaryGOMEAMemoryForPopulationAndOffspring( int BinaryGOMEA_index )
{
  int i;

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
    free( BinaryGOMEA_offsprings[BinaryGOMEA_index][i] );

  for( i = 0; i < BinaryGOMEA_population_sizes[BinaryGOMEA_index]; i++ )
    free( BinaryGOMEA_populations[BinaryGOMEA_index][i] );
  
  free( BinaryGOMEA_populations[BinaryGOMEA_index] );
  free( BinaryGOMEA_fitnesses[BinaryGOMEA_index] );
  free( BinaryGOMEA_offsprings[BinaryGOMEA_index] );
  free( BinaryGOMEA_fitnesses_offsprings[BinaryGOMEA_index] );
}


void BinaryGOMEA_EzilaitiniProblem()
{
  switch( BinaryGOMEA_problem_index )
  {
    case  0: break;
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/



/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Runs BinaryGOMEA in the Interleaved Multistart Scheme.
 */
void BinaryGOMEA_RunIMS()
{
  BinaryGOMEA_IMS_subgeneration_factor                   = 4;
  BinaryGOMEA_base_population_size                       = 1;
  BinaryGOMEA_number_of_BinaryGOMEAs                     = 0;
  BinaryGOMEA_number_of_generations_IMS                  = 0;
  BinaryGOMEA_number_of_evaluations                      = 0;
  BinaryGOMEA_minimum_BinaryGOMEA_index                  = 0;
  BinaryGOMEA_elitist_solution_written_before            = 0;
  BinaryGOMEA_elitist_solution_last_written_running_time = 0;
  BinaryGOMEA_first_evaluation_ever                      = 1;
  BinaryGOMEA_elitist_solution_BinaryGOMEA_index         = 0;
  BinaryGOMEA_elitist_solution_offspring_index           = 0;

  while( !BinaryGOMEA_CheckTermination() )
  {
    if( BinaryGOMEA_number_of_BinaryGOMEAs < BinaryGOMEA_maximum_number_of_BinaryGOMEAs )
      BinaryGOMEA_InitializeNewBinaryGOMEA();

    BinaryGOMEA_GenerationalStepAllBinaryGOMEAs();

    BinaryGOMEA_number_of_generations_IMS++;

    //BinaryGOMEA_WriteElitistSolutionToFile();
  }

  BinaryGOMEA_EzilaitiniAllBinaryGOMEAs();
}

/**
 * Initializes the random number generator and the problem and runs the BinaryGOMEA.
 */
void BinaryGOMEA_Run()
{
  BinaryGOMEA_timestamp_start = GetCurrentTimeStampInMilliSeconds();

  BinaryGOMEA_InitializePredeterminedFOS();
  
  BinaryGOMEA_InitializeRandomNumberGenerator();
  
  BinaryGOMEA_InitializeProblem();

  if (BinaryGOMEA_use_VTR == 1)
    BinaryGOMEA_vtr = load_vtr(folder);

  if( BinaryGOMEA_print_verbose_overview )
    BinaryGOMEA_PrintVerboseOverview();

  BinaryGOMEA_timestamp_start_after_init = GetCurrentTimeStampInMilliSeconds();

  BinaryGOMEA_RunIMS();

  BinaryGOMEA_EzilaitiniProblem();
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Main -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
int main( int argc, char **argv )
{
  BinaryGOMEA_InterpretCommandLine( argc, argv );
  BinaryGOMEA_Run();

  return( 0 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
