#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

int n;
int **distances;
int *solution;

int64_t GetUs() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (int64_t)tv.tv_sec * 1000000 + (int64_t)tv.tv_usec;
}

double GetDistance(int i, int j) { return 1.0 * distances[i][j]; }

// 获得一个初始解
void InitSolution() {
  int *ids = malloc(n * sizeof(int)); // 还没进入初始序列的解
  for (int i = 0; i < n; i++) {
    ids[i] = i;
  }
  int lastpos = n - 1;  // ids 最后一个有效位置
  int nxt = rand() % n; // 下个节点在 ids 中的序号
  solution[0] = ids[nxt];
  ids[nxt] = ids[lastpos];
  lastpos--;
  for (int i = 1; i < n; i++) {
    int lst = solution[i - 1]; // 当前节点编号
    double mxdis = INFINITY;
    nxt = lst;
    for (int j = 0; j <= lastpos; j++) {
      double dis = GetDistance(lst, ids[j]);
      if (dis < mxdis) {
        mxdis = dis;
        nxt = j;
      }
    }
    solution[i] = ids[nxt];
    ids[nxt] = ids[lastpos];
    lastpos--;
  }
  free(ids);
}

// 输出最终 solution 及其结果
void PrintMyResults(int64_t start_us) {
  double ans = 0.0;
  for (int i = 0; i < n; i++) {
    // printf("%d ", solution[i]);
    ans += GetDistance(solution[i], solution[(i + 1) % n]);
  }
  int64_t used_us = GetUs() - start_us;
  printf("\nFor %d node atsp, my ans is %.4f. cost %.2f s\n", n, ans,
         used_us / 1000000.0);
}

// 模拟退火解决 ATSP 问题 主函数
void SimulatedAnnealing() {
  int64_t start_us = GetUs();

  InitSolution();
  
  // 模拟退火
  // TODO(cyx):...

  PrintMyResults(start_us);
}

int main(int argc, char *argv[]) {
  srand(time(NULL));
  if (argc != 2) {
    printf("Usage: %s matrix_file\n", argv[0]);
    return 0;
  }
  FILE *fp = fopen(argv[1], "r");
  if (fp == NULL) {
    printf("cannot open %s\n", argv[1]);
    return 0;
  }
  fscanf(fp, "%d", &n);
  distances = malloc(n * sizeof(int *));
  for (int i = 0; i < n; i++) {
    distances[i] = malloc(n * sizeof(int));
    for (int j = 0; j < n; j++) {
      fscanf(fp, "%d", &distances[i][j]);
    }
  }
  fclose(fp);
  solution = malloc(n * sizeof(int));
  // calc
  SimulatedAnnealing();
  // clean
  for (int i = 0; i < n; i++) {
    free(distances[i]);
  }
  free(distances);
  free(solution);
  return 0;
}