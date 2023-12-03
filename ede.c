#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POP_SIZE 80         // Maximum population size
#define NUM_DECISION_VARS 6     // Number of decision variables

// Function declarations
double rand_uniform(double min, double max);  // Function to generate a random number in a given range
void initialize_population(double pop[][NUM_DECISION_VARS], double VarMin[], double VarMax[], int nPop);
double calculate_cost(double x[], double target_field[], double meas_point_set[][3], int nVar);
void mutate_and_crossover(double pop[][NUM_DECISION_VARS], double target_field[], double meas_point_set[][3],
                           double VarMin[], double VarMax[], int nPop, int nVar,
                           double BestSol[], double BestCost, double beta_min, double beta_max, double pCR);

int main() {
    // Problem Definition
    double target_field[] = { /* Initialize with your target field values */ };
    double meas_point_set[][3] = { /* Initialize with your measured point set values */ };

    // DE Parameters
    int MaxIt = 50;                // Maximum number of iterations
    int nPop = MAX_POP_SIZE;       // Population size
    double beta_min = 0.5;         // Minimum value for mutation parameter beta
    double beta_max = 1;           // Maximum value for mutation parameter beta
    double pCR = 0.9;              // Crossover probability

    // Variable Bounds
    double VarMin[NUM_DECISION_VARS] = {-30, -30, -30, -500, -500, -500};  // Minimum values for decision variables
    double VarMax[NUM_DECISION_VARS] = {30, 30, 30, 500, 500, 500};        // Maximum values for decision variables

    // Initialization
    double pop[MAX_POP_SIZE][NUM_DECISION_VARS];   // Population matrix
    double BestSol[NUM_DECISION_VARS];             // Best solution vector
    double BestCost = INFINITY;                    // Best cost initialized to infinity

    initialize_population(pop, VarMin, VarMax, nPop);  // Initialize the population with random values

    // Main Loop
    for (int it = 0; it < MaxIt; ++it) {
        mutate_and_crossover(pop, target_field, meas_point_set, VarMin, VarMax, nPop, NUM_DECISION_VARS, BestSol, BestCost, beta_min, beta_max, pCR);

        // Update Best Cost
        BestCost = calculate_cost(BestSol, target_field, meas_point_set, NUM_DECISION_VARS);

        // Show Iteration Information
        if (it % 10 == 0) {
            printf("Iteration %d: Best Cost = %lf\n", it, BestCost);
        }

        // Termination Condition
        if (BestCost <= 1e-30) {
            printf("Termination condition met.\n");
            break;
        }
    }

    // Show Results
    printf("Best Solution:\n");
    for (int i = 0; i < NUM_DECISION_VARS; ++i) {
        printf("%lf ", BestSol[i]);
    }
    printf("\nBest Cost: %lf\n", BestCost);

    return 0;
}

// Function to initialize the population with random values within variable bounds
void initialize_population(double pop[][NUM_DECISION_VARS], double VarMin[], double VarMax[], int nPop) {
    for (int i = 0; i < nPop; ++i) {
        for (int j = 0; j < NUM_DECISION_VARS; ++j) {
            pop[i][j] = rand_uniform(VarMin[j], VarMax[j]);
        }
    }
}

// Function to calculate the cost of a solution based on the provided cost function
double calculate_cost(double x[], double target_field[], double meas_point_set[][3], int nVar) {
    // Placeholder: Implement your cost function here
    // For example, calculate the sum of squared differences between model and target fields
    double cost = 0.0;
    for (int i = 0; i < nVar; ++i) {
        // Modify this line based on your actual cost function
        cost += pow(x[i] - target_field[i], 2);
    }
    return cost;
}

// Function implementing mutation and crossover steps of the DE algorithm
void mutate_and_crossover(double pop[][NUM_DECISION_VARS], double target_field[], double meas_point_set[][3],
                           double VarMin[], double VarMax[], int nPop, int nVar,
                           double BestSol[], double BestCost, double beta_min, double beta_max, double pCR) {
    for (int i = 0; i < nPop; ++i) {
        // Mutation and Crossover logic
        int a, b, c;
        // Ensure a, b, c are distinct and not equal to i
        do {
            a = rand() % nPop;
        } while (a == i);
        do {
            b = rand() % nPop;
        } while (b == i || b == a);
        do {
            c = rand() % nPop;
        } while (c == i || c == a || c == b);

        // Mutation
        double beta1 = rand_uniform(beta_min, beta_max);
        double beta2 = rand_uniform(beta_min, beta_max);
        double rmax = 1.0;
        double rmin = 0.0;
        double r = rmax - (i / MaxIt) * (rmax - rmin);
        double rr = rand_uniform(0, 1);

        double y[NUM_DECISION_VARS];
        for (int j = 0; j < nVar; ++j) {
            // DE/current-to-best/1
            y[j] = pop[i][j] + beta1 * (pop[b][j] - pop[c][j]) + beta2 * (BestSol[j] - pop[i][j]);
            // Ensure the mutated value is within bounds
            y[j] = fmax(VarMin[j], fmin(y[j], VarMax[j]));
        }

        // Crossover
        double z[NUM_DECISION_VARS];
        int j0 = rand() % nVar;
        for (int j = 0; j < nVar; ++j) {
            if (j == j0 || rand_uniform(0, 1) <= pCR) {
                z[j] = y[j];
            } else {
                z[j] = pop[i][j];
            }
        }

        // Evaluate the cost of the new solution
        double current_cost = calculate_cost(z, target_field, meas_point_set, nVar);

        // Update BestSol if a better solution is found
        if (current_cost < BestCost) {
            for (int j = 0; j < nVar; ++j) {
                BestSol[j] = z[j];
            }
            BestCost = current_cost;
        }
    }
}

// Function to generate a random number between min and max (uniform distribution)
double rand_uniform(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min);
}
