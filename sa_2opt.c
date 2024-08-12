#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

const int64_t SA_MAX_RUN_US = 18.5 * 1000000; // 最多允许模拟退火多少秒
const double INIT_ACCEPT_P = 0.25; // 一开始能以多大概率接收差解
const double SA_RATE = 0.99;       // 模拟退火降温率
const int FINAL_ACCEPT_RATE_DAO = 100000; // FINAL_ACCEPT_RATE 的倒数
const double FINAL_ACCEPT_RATE =
    1.0 / FINAL_ACCEPT_RATE_DAO; // 退火过程中温度不得低于多少

int n; // 节点总数
int sqrt_n;
int *point_x;       // 坐标 X
int *point_y;       // 坐标 Y
int *best_solution; // 历史最好解
double best_dis;    // 历史最好距离

static inline int64_t GetUs() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (int64_t)tv.tv_sec * 1000000 + (int64_t)tv.tv_usec;
}

static inline double GetDistance(int i, int j) {
  return sqrt((long long)(point_x[i] - point_x[j]) * (point_x[i] - point_x[j]) +
              (long long)(point_y[i] - point_y[j]) * (point_y[i] - point_y[j]));
}

static inline void SwapInt(int *i, int *j) {
  int tmp = *i;
  *i = *j;
  *j = tmp;
}

// 获得一个初始解
void InitSolution() {
  int *ids = (int *)malloc(n * sizeof(int)); // 还没进入初始序列的解
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

static inline int ChkRand(double p) {
  int acceptance = p * FINAL_ACCEPT_RATE_DAO;
  int rnd = (long long)rand() * rand() % FINAL_ACCEPT_RATE_DAO;
  if (rnd < acceptance) {
    return 1;
  }
  return 0;
}

double CalcSolutionDis(int *sol) {
  double ans = 0.0;
  for (int i = 0; i < n; i++) {
    ans += GetDistance(sol[i], sol[(i + 1) % n]);
  }
  return ans;
}

// 输出最终 best_solution 及其结果
void PrintMyResults(int64_t start_us) {
  double ans = CalcSolutionDis(best_solution);
  int64_t used_us = GetUs() - start_us;
  printf("\nFor %d node atsp, my ans is %.4f best_dis %.4f. cost %.2f s\n", n,
         ans, best_dis, used_us / 1000000.0);
}

// 尝试交换两节点。如果更好，则更新到 best_solution 中
void Opt2PointExchange(int *curr_solution, int *tmp_solution, double *curr_dis,
                       double *temprature) {
  // 构造一个新的 curr_solution。即：随机 swap i、j
  int i = rand() % n;
  int j = rand() % n;
  if (i == j) {
    return;
  }
  if (i > j) {
    SwapInt(&i, &j);
  }
  // NOTE(cyx): 在移植以后，这个可能还会变
  double new_dis = *curr_dis;
  int pre_i = (i - 1 + n) % n;
  int nxt_i = (i + 1) % n;
  int pre_j = (j - 1 + n) % n;
  int nxt_j = (j + 1) % n;
  // 缩水版 2-opt
  if (i == 0 && j == n - 1) {
    new_dis -= GetDistance(curr_solution[pre_j], curr_solution[j]);
    new_dis -= GetDistance(curr_solution[j], curr_solution[i]);
    new_dis -= GetDistance(curr_solution[i], curr_solution[nxt_i]);
    new_dis += GetDistance(curr_solution[pre_j], curr_solution[i]);
    new_dis += GetDistance(curr_solution[i], curr_solution[j]);
    new_dis += GetDistance(curr_solution[j], curr_solution[nxt_i]);
  } else if (nxt_i != j) {
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
  double delta_dis = new_dis - *curr_dis;
  if (*temprature < -0.5) {
    // exp（-Δt′/T)=p -> -Δt′/t = ln(p) -> t = -Δt′/ln(p)
    double wanna_dis = best_dis * 0.001;
    *temprature = -fabs(wanna_dis) / log(INIT_ACCEPT_P);
  }
  if (0) {
    printf("After 2point xchange %d %d, delta_dis %.2f new_dis %.2f best_dis "
           "%.2f temp %.2f p %.2f\n",
           i, j, delta_dis, new_dis, best_dis, *temprature,
           exp(-delta_dis / *temprature));
  }
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / *temprature))) {
    *curr_dis = new_dis;
    SwapInt(&curr_solution[i], &curr_solution[j]);
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}

// 尝试翻转一个 range。如果更好，则更新到 best_solution 中
// 具体而言，有头、中、尾实现
void OptRangeReverse(int *curr_solution, int *tmp_solution, double *curr_dis,
                     double *temprature) {
  // int wanna_len = 3 + rand() % sqrt_n;
  // int wanna_len = sqrt_n;
  // int idx_start = rand() % (n - wanna_len + 1);
  // int idx_end = idx_start + wanna_len - 1;

  int idx_start = rand() % n;
  int idx_end = rand() % n;
  if (idx_start > idx_end) {
    SwapInt(&idx_start, &idx_end);
  }
  if (idx_end - idx_start + 1 == n || idx_end == idx_start) { // 没法搞了
    return;
  }
  if (n > 1200) {
    int wanna_len = rand() % 600 + 3;
    idx_start = rand() % (n - wanna_len + 1);
    idx_end = idx_start + wanna_len - 1;
  }
  // if (sqrt_n > 40 && idx_start + 10 * sqrt_n < idx_end) {
  //   idx_end = idx_start + rand() % (10 * sqrt_n);
  // }
  // 概率地决定是反转最前面、最后面，还是中间
  // 这么生成期望长度是 n / 3
  {
    int rnd = rand() % 1000;
    int len = idx_end - idx_start + 1;
    if (rnd < 200) { // 20% 概率反转开头
      idx_start = 0;
      idx_end = len - 1;
    } else if (rnd < 400) { // 20% 概率反转结尾
      idx_start = n - len;
      idx_end = n - 1;
    }
  }
  double new_dis = *curr_dis;
  // 先内部 new_dis 变化
  for (int i = idx_start + 1; i <= idx_end; i++) {
    new_dis -= GetDistance(curr_solution[i - 1], curr_solution[i]);
    new_dis += GetDistance(curr_solution[i], curr_solution[i - 1]);
  }
  // 再外部 new_dis 变化
  int prev = (idx_start - 1 + n) % n;
  new_dis -= GetDistance(curr_solution[prev], curr_solution[idx_start]);
  new_dis += GetDistance(curr_solution[prev], curr_solution[idx_end]);
  int nxt = (idx_end + 1) % n;
  new_dis -= GetDistance(curr_solution[idx_end], curr_solution[nxt]);
  new_dis += GetDistance(curr_solution[idx_start], curr_solution[nxt]);
  // Metropolis准则：若Δt′<0则接受S′作为新的当前解S，否则以概率exp（-Δt′/T）接受S′作为新的当前解S
  double delta_dis = new_dis - *curr_dis;
  if (*temprature < -0.5) {
    // exp（-Δt′/T)=p -> -Δt′/t = ln(p) -> t = -Δt′/ln(p)
    double wanna_dis = best_dis * 0.001;
    *temprature = -fabs(wanna_dis) / log(INIT_ACCEPT_P);
  }
  if (0) {
    printf("After range xchange %d %d, delta_dis %.2f new_dis %.2f best_dis "
           "%.2f temp %.2f p %.2f\n",
           idx_start, idx_end, delta_dis, new_dis, best_dis, *temprature,
           exp(-delta_dis / *temprature));
  }
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / *temprature))) {
    *curr_dis = new_dis;
    int swapi = idx_start;
    int swapj = idx_end;
    while (swapi < swapj) {
      SwapInt(&curr_solution[swapi], &curr_solution[swapj]);
      swapi++;
      swapj--;
    }
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}

// 尝试插入一个
void OptInsert(int *curr_solution, int *tmp_solution, double *curr_dis,
               double *temprature) {
  // 把 idx_i 挪到 idx_j 后头
  int idx_i = rand() % n;
  int idx_j = rand() % n;
  if (idx_i > idx_j) {
    SwapInt(&idx_i, &idx_j);
  }
  if (idx_j - idx_i + 1 == n || idx_j == idx_i) { // 没法搞了
    return;
  }
  // if (n > 1200) {
  //   int wanna_len = rand() % 600 + 3;
  //   idx_i = rand() % (n - wanna_len + 1);
  //   idx_j = idx_i + wanna_len - 1;
  // }
  int is_qian = rand() % 2; // 为 1 则把 idx_i 挪到 idx_j 后头
  double new_dis = *curr_dis;
  if (is_qian) {
    int pre_i = (idx_i - 1 + n) % n;
    int nxt_i = (idx_i + 1) % n;
    int nxt_j = (idx_j + 1) % n;
    new_dis -= GetDistance(curr_solution[pre_i], curr_solution[idx_i]);
    new_dis -= GetDistance(curr_solution[idx_i], curr_solution[nxt_i]);
    new_dis -= GetDistance(curr_solution[idx_j], curr_solution[nxt_j]);
    new_dis += GetDistance(curr_solution[pre_i], curr_solution[nxt_i]);
    new_dis += GetDistance(curr_solution[idx_j], curr_solution[idx_i]);
    new_dis += GetDistance(curr_solution[idx_i], curr_solution[nxt_j]);
  } else {
    int pre_i = (idx_i - 1 + n) % n;
    int pre_j = (idx_j - 1 + n) % n;
    int nxt_j = (idx_j + 1) % n;
    new_dis -= GetDistance(curr_solution[pre_i], curr_solution[idx_i]);
    new_dis -= GetDistance(curr_solution[pre_j], curr_solution[idx_j]);
    new_dis -= GetDistance(curr_solution[idx_j], curr_solution[nxt_j]);
    new_dis += GetDistance(curr_solution[pre_i], curr_solution[idx_j]);
    new_dis += GetDistance(curr_solution[idx_j], curr_solution[idx_i]);
    new_dis += GetDistance(curr_solution[pre_j], curr_solution[nxt_j]);
  }
  // Metropolis准则：若Δt′<0则接受S′作为新的当前解S，否则以概率exp（-Δt′/T）接受S′作为新的当前解S
  double delta_dis = new_dis - *curr_dis;
  if (*temprature < -0.5) {
    // exp（-Δt′/T)=p -> -Δt′/t = ln(p) -> t = -Δt′/ln(p)
    double wanna_dis = best_dis * 0.001;
    *temprature = -fabs(wanna_dis) / log(INIT_ACCEPT_P);
  }
  if (0) {
    printf("After range xchange %d %d, delta_dis %.2f new_dis %.2f best_dis "
           "%.2f temp %.2f p %.2f\n",
           idx_i, idx_j, delta_dis, new_dis, best_dis, *temprature,
           exp(-delta_dis / *temprature));
  }
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / *temprature))) {
    *curr_dis = new_dis;
    if (is_qian) {
      int tmp = curr_solution[idx_i];
      for (int i = idx_i + 1; i <= idx_j; i++) {
        curr_solution[i - 1] = curr_solution[i];
      }
      curr_solution[idx_j] = tmp;
    } else {
      int tmp = curr_solution[idx_j];
      for (int i = idx_j - 1; i >= idx_i; i--) {
        curr_solution[i + 1] = curr_solution[i];
      }
      curr_solution[idx_i] = tmp;
    }
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      // poorpool debug
      // double qwq_dis = CalcSolutionDis(curr_solution);
      // if (fabs(qwq_dis - *curr_dis) > 0.5) {
      //   printf("strange qwq_dis %.4f curr_dis %.4f idx_i %d idx_j %d\n",
      //          qwq_dis, *curr_dis, idx_i, idx_j);
      // }
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}
// 随机生成排列
void RandomGeneratePermutation(int *arr) {
  for (int i = 0; i < n; i++) {
    arr[i] = i;
  }
  for (int i = n - 1; i > 0; i--) {
    int j = rand() % (i + 1);
    SwapInt(&arr[i], &arr[j]);
  }
}

// 模拟退火解决 ATSP 问题 主函数
void SimulatedAnnealing() {
  int64_t start_us = GetUs();
  int64_t should_end_us = start_us + SA_MAX_RUN_US;

  // 贪心初始化一个解
  InitSolution();

  // 两个缓冲区
  int *curr_solution = (int *)malloc(n * sizeof(int));
  int *tmp_solution = (int *)malloc(n * sizeof(int));

  // 模拟退火
  while (GetUs() < should_end_us) {
    // 复制历史最好解，按这个退火
    // 不按照历史最好解，而是随机一点会怎么样？实测效果不咋样。。。
    double curr_dis = best_dis;
    memcpy(curr_solution, best_solution, n * sizeof(int));

    // // TODO(cyx): 调试随机
    // RandomGeneratePermutation(curr_solution);
    // curr_dis = CalcSolutionDis(curr_solution);

    double temprature = -1.0; // 当前温度。起初为 -1，后面根据第一次结果生成
    do {
      int rnd = rand() % 1000; // 决定操作概率
      if (rnd < 500) {         // 50% 概率交换两点
        Opt2PointExchange(curr_solution, tmp_solution, &curr_dis, &temprature);
      } else if (rnd < 750) { // 25%
        // 概率反转区间。不过这个概率好像也不至于这么大。主要是节点越多，越不适合进行这个操作！这个操作是
        // O(n) 的
        // 1000 的时候，这玩意占 80% 还不错。10000
        // 的时候，这玩意要么设置成0，要么长度限制在sqrt(n)
        // 长度不是越长就越好的
        OptRangeReverse(curr_solution, tmp_solution, &curr_dis, &temprature);
      } else {
        OptInsert(curr_solution, tmp_solution, &curr_dis, &temprature);
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
  sqrt_n = sqrt(n);
  point_x = (int *)malloc(n * sizeof(int));
  point_y = (int *)malloc(n * sizeof(int));
  for (int i = 0; i < n; i++) {
    fscanf(fp, "%d %d", &point_x[i], &point_y[i]);
  }
  fclose(fp);
  best_solution = (int *)malloc(n * sizeof(int));
  // calc
  SimulatedAnnealing();
  // clean
  free(point_x);
  free(point_y);
  free(best_solution);
  return 0;
}