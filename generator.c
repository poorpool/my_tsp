#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int nodes[] = {10, 40, 70, 140, 250, 333, 1017, 2579, 5298, 7748, 10000};

int main(int argc, char *argv[]) {
  srand(time(NULL));
  char filename[25];
  int dataset_num = sizeof(nodes) / sizeof(int);
  for (int i = 0; i < dataset_num; i++) {
    int n = nodes[i];
    sprintf(filename, "datasets/my_data%d.in", i);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d\n", n);
    for (int j = 0; j < n; j++) {
      int px = rand() % 10001;
      int py = rand() % 10001;
      fprintf(fp, "%d %d\n", px, py);
    }
    fclose(fp);
  }
  return 0;
}