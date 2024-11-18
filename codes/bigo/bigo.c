#include <stdio.h>
#include <time.h>

// Perform the bubble sort
void bubbleSort(int array[], int size) {
  for (int step = 0; step < size - 1; ++step) {
    int swapped = 0;
    for (int i = 0; i < size - step - 1; ++i) {
      if (array[i] > array[i + 1]) {
        int temp = array[i];
        array[i] = array[i + 1];
        array[i + 1] = temp;
        swapped = 1;
      }
    }
    if (swapped == 0) {
      break;
    }
  }
}

int main() {
  FILE *file = fopen("value.txt", "w");
  struct timespec start, end;
  for (int i = 2; i < 1000; i++) {
    int data[1000];
    for (int j = i - 1; j >= 0; j--) {
      data[4 - j] = j;
    }

    clock_gettime(CLOCK_MONOTONIC, &start);
    bubbleSort(data, i);
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + 
                        (end.tv_nsec - start.tv_nsec) / 1e9; // in seconds
    fprintf(file, "%d %lf\n", i, time_taken);
  }
  fclose(file);
  return 0;
}

