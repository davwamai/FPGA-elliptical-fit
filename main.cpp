#include <stdio.h>
#include <stdlib.h>
#include "fitzgibbon.h"

#include <stdio.h>
#include <stdlib.h>

#define MAX_POINTS 300

// static memory allocation for x and y values
static double x_values[MAX_POINTS];
static double y_values[MAX_POINTS];
static double a [6];

int read_csv(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return -1;
    }

    char buffer[256];  // buffer to hold each line
    int index = 0;

    // skip header line
    fgets(buffer, sizeof(buffer), file);

    while (fgets(buffer, sizeof(buffer), file) && index < MAX_POINTS) {
    if (sscanf(buffer, "%lf,%lf", &x_values[index], &y_values[index]) != 2) {
        // error handling for malformed lines
        printf("Error reading line: %s\n", buffer);
        fclose(file);
        return -1;
    }
    index++;
}

    fclose(file);
    return index;  // returns the number of points read
}

int main() {
    int num_points = read_csv("ellipse_points.csv");
    if (num_points < 0) {
        printf("Error reading CSV file.\n");
        return -1;
    }

    // prints first 10 points to verify
    for (int i = 0; i < 10; i++) {
        printf("Point %d: (%lf, %lf)\n", i+1, x_values[i], y_values[i]);
    }
    printf("...\n\n");
    
    fitzgibbon(x_values, y_values, num_points, a);
    return 0;
}