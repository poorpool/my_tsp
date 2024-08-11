#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

const int64_t SA_MAX_RUN_US = 18.5 * 1000000; // 最多允许模拟退火多少秒
const double INIT_ACCEPT_P = 0.75; // 一开始能以多大概率接收差解
const double SA_RATE = 0.99;       // 模拟退火降温率
const int FINAL_ACCEPT_RATE_DAO = 100000; // FINAL_ACCEPT_RATE 的倒数
const double FINAL_ACCEPT_RATE =
    1.0 / FINAL_ACCEPT_RATE_DAO; // 退火过程中温度不得低于多少

int n;              // 节点总数
int **distances;    // 距离矩阵（中小规模使用）
int *best_solution; // 历史最好解
double best_dis;    // 历史最好距离

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
  best_solution[0] = ids[nxt];
  ids[nxt] = ids[lastpos];
  lastpos--;
  for (int i = 1; i < n; i++) {
    int lst = best_solution[i - 1]; // 当前节点编号
    double mxdis = INFINITY;
    nxt = lst;
    for (int j = 0; j <= lastpos; j++) {
      double dis = GetDistance(lst, ids[j]);
      if (dis < mxdis) {
        mxdis = dis;
        nxt = j;
      }
    }
    best_solution[i] = ids[nxt];
    ids[nxt] = ids[lastpos];
    lastpos--;
  }
  free(ids);
  best_dis = 0.0;
  for (int i = 0; i < n; i++) {
    best_dis += GetDistance(best_solution[i], best_solution[(i + 1) % n]);
  }
  printf("at first, best_dis %.4f\n", best_dis);
}

int chk_rand(double p) {
  int acceptance = p * FINAL_ACCEPT_RATE_DAO;
  int rnd = (long long)rand() * rand() % FINAL_ACCEPT_RATE_DAO;
  if (rnd < acceptance) {
    return 1;
  }
  return 0;
}

// 输出最终 best_solution 及其结果
void PrintMyResults(int64_t start_us) {
  double ans = 0.0;
  for (int i = 0; i < n; i++) {
    // printf("%d ", best_solution[i]);
    ans += GetDistance(best_solution[i], best_solution[(i + 1) % n]);
  }
  int64_t used_us = GetUs() - start_us;
  printf("\nFor %d node atsp, my ans is %.4f best_dis %.4f. cost %.2f s\n", n,
         ans, best_dis, used_us / 1000000.0);
}

// 模拟退火解决 ATSP 问题 主函数
void SimulatedAnnealing() {
  int64_t start_us = GetUs();
  int64_t should_end_us = start_us + SA_MAX_RUN_US;

  // 贪心初始化一个解
  InitSolution();

  int *curr_solution = malloc(n * sizeof(int));
  int *tmp_solution = malloc(n * sizeof(int));

  // 模拟退火
  while (GetUs() < should_end_us) {
    // 复制历史最好解，按这个退火
    // TODO(cyx): 不按照历史最好解，而是随机一点会怎么样？
    double curr_dis = best_dis;
    memcpy(curr_solution, best_solution, n * sizeof(int));

    double temprature = -1.0; // 当前温度。起初为 -1，后面根据第一次结果生成
    do {
      // printf("curr temp %.3f\n", temprature);
      // 构造一个新的 curr_solution。即：随机 swap i、j
      int i = rand() % n;
      int j = rand() % n;
      if (i == j) {
        continue;
      }
      if (i > j) {
        int tmp = j;
        j = i;
        i = j;
      }
      // NOTE(cyx): 在移植以后，这个可能还会变
      double new_dis = curr_dis;
      int pre_i = (i - 1 + n) % n;
      int nxt_i = (i + 1) % n;
      int pre_j = (j - 1 + n) % n;
      int nxt_j = (j + 1) % n;
      // 缩水版 2-opt
      if (nxt_i != j) {
        new_dis -= GetDistance(curr_solution[pre_i], curr_solution[i]);
        new_dis -= GetDistance(curr_solution[i], curr_solution[nxt_i]);
        new_dis -= GetDistance(curr_solution[pre_j], curr_solution[j]);
        new_dis -= GetDistance(curr_solution[j], curr_solution[nxt_j]);
        new_dis += GetDistance(curr_solution[pre_i], curr_solution[j]);
        new_dis += GetDistance(curr_solution[j], curr_solution[nxt_i]);
        new_dis += GetDistance(curr_solution[pre_j], curr_solution[i]);
        new_dis += GetDistance(curr_solution[i], curr_solution[nxt_j]);
      } else {
        new_dis -= GetDistance(curr_solution[pre_i], curr_solution[i]);
        new_dis -= GetDistance(curr_solution[i], curr_solution[j]);
        new_dis -= GetDistance(curr_solution[j], curr_solution[nxt_j]);
        new_dis += GetDistance(curr_solution[pre_i], curr_solution[j]);
        new_dis += GetDistance(curr_solution[j], curr_solution[i]);
        new_dis += GetDistance(curr_solution[i], curr_solution[nxt_j]);
      }
      // Metropolis准则：若Δt′<0则接受S′作为新的当前解S，否则以概率exp（-Δt′/T）接受S′作为新的当前解S
      double delta_dis = new_dis - curr_dis;
      if (temprature < -0.5) {
        // exp（-Δt′/T)=p -> -Δt′/t = ln(p) -> t = -Δt′/ln(p)
        temprature = -fabs(delta_dis) / log(INIT_ACCEPT_P);
      }
      if (delta_dis < 0.0 || chk_rand(exp(-delta_dis / temprature))) {
        curr_dis = new_dis;
        int tmp = curr_solution[i];
        curr_solution[i] = curr_solution[j];
        curr_solution[j] = tmp;
        if (curr_dis < best_dis) {
          best_dis = curr_dis;
          memcpy(best_solution, curr_solution, n * sizeof(int));
        }
      }
      temprature *= SA_RATE;
    } while (temprature > FINAL_ACCEPT_RATE);
  }

  free(tmp_solution);
  free(curr_solution);
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
  best_solution = malloc(n * sizeof(int));
  // calc
  SimulatedAnnealing();
  // clean
  for (int i = 0; i < n; i++) {
    free(distances[i]);
  }
  free(distances);
  free(best_solution);
  return 0;
}