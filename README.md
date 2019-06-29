# LinKernighanAlgorithm
My implementation of the Lin-Kernighan-algorithm as described in Combinatorial optimization by Korte and Vygen with 
improvements by Keld Helsgaun as part of a course at the Uni Bonn

## Compiling

    cmake --build cmake-build-debug --target LinKernighanAlgorithm -- -j 6

## Usage:
    LinKernighanAlgorithm --help
    LinKernighanAlgorithm tsplib_problem.tsp [options]

## Options:
    --number-of-trials=integer
        Sets the maximum number of trials (default: 50)
    --candidate-edges=[ALL|NEAREST|ALPHA_NEAREST|OPT_ALPHA_NEAREST]
        Set the choice of candidate edges for each vertex (default: OPT_ALPHA_NEAREST)
            ALL: all incident edges
            NEAREST: shortest k incident edges
            ALPHA_NEAREST: shortest k incident edges by alpha distance
            OPT_ALPHA_NEAREST: shortest k incident edges by alpha distance after subgradient optimization
    --number-of-candidate-edges=integer
        Set the number of candidate edges k for each vertex. Ignored for --candidate-edges=ALL. (default: 5)
    --dont-store-distances
        The distances between all vertices are not computed once and then stored in a matrix, but are computed when
        needed. This leads to fewer memory usage but also to a considerable decrease in performance.
    --optimum-tour-length=integer
        Set the optimum tour length for the problem. The output will include the percentage error of the solution.
    --acceptable-error=double
        Set the acceptable error compared to the optimum length in percent. The program will stop when a found tour
        is within the acceptable length range. Ignored if --optimum-value is not given. (default: 0)
    --output-to-file
        Output the best tour found to "tsplib_problem.lk.tour"
    --verbose
        Output debug output, useful to estimate the remaining running time. If not given there is no output until the
        result is computed.
        
The data structure for the tour cannot be changed with a command line flag. To set it you must change the last line in 
`Tour.h` to one of these options

```objectivec
using Tour = TwoLevelTreeTour;
using Tour = ArrayTour;
```
and recompile the source code.

## Example:
    LinKernighanAlgorithm ExampleProblems/si1032.tsp --optimum-tour-length=92650 --acceptable-error=0.1 --verbose

## Restrictions:
The dimension of the TSPLIB problem may not be smaller than 3 or the number of candidate edges plus 1.
