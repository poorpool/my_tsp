#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

// 最多允许模拟退火多少微秒
const int64_t SA_MAX_RUN_US = 18.5 * 1000000;
// 一开始能以多大概率接收差解
// 0.25 效果最好，0.5 略次，0.75 不太好
const double INIT_ACCEPT_P = 0.25;
// 模拟退火降温率
// 0.97/98/99 都差不多
// 0.999 不太好
const double SA_COOL_RATE = 0.97;
// 起初能以 INIT_ACCEPT_P 的概率接收比最好值差多少的解
// 0.001 就不错。更高了就不太好
const double SA_TOLERATE_RATE = 0.001;
// 退火结束时，温度为初始温度的比值
// 比如 exp(-Δt/T)=0.25 的话，exp(-Δt/(T*0.1)) 就很小很小了
const double MIN_SA_COOL_RATIO = 0.1;

int n;              // 节点总数
int *point_x;       // X 坐标数组
int *point_y;       // Y 坐标数组
int *best_solution; // 历史最好解数组
double best_dis;    // 历史最好距离

// 读取当前微秒数
static inline int64_t GetUs() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (int64_t)tv.tv_sec * 1000000 + (int64_t)tv.tv_usec;
}

// 获取两点间距离，不论距离到底是如何计算的
static inline double GetDistance(int i, int j) {
  return sqrt((long long)(point_x[i] - point_x[j]) * (point_x[i] - point_x[j]) +
              (long long)(point_y[i] - point_y[j]) * (point_y[i] - point_y[j]));
}

// 交换两个整数
static inline void SwapInt(int *i, int *j) {
  int tmp = *i;
  *i = *j;
  *j = tmp;
}

// 根据概率 p 返回 1 或 0
static inline int ChkRand(double p) {
  int acceptance = p * 1000000;
  int rnd = (long long)rand() * rand() % 1000000;
  if (rnd < acceptance) {
    return 1;
  }
  return 0;
}

// 获得一条路径的总距离
static inline double CalcTourDistance(int *tour) {
  double ans = 0.0;
  for (int i = 0; i < n; i++) {
    ans += GetDistance(tour[i], tour[(i + 1) % n]);
  }
  return ans;
}

// 贪心获得一个初始解
void InitSolution() {
  // 初始化还没进入解序列的节点编号
  int *ids = (int *)malloc(n * sizeof(int));
  for (int i = 0; i < n; i++) {
    ids[i] = i;
  }

  // 随机设置起始节点
  int lastpos = n - 1;  // ids 最后一个有效位置
  int nxt = rand() % n; // 下个节点在 ids 中的序号
  best_solution[0] = ids[nxt];
  ids[nxt] = ids[lastpos];
  lastpos--;

  // 贪心获得后续节点
  for (int i = 1; i < n; i++) {
    int lst = best_solution[i - 1]; // 当前节点编号
    double mxdis = DBL_MAX;
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

  // 计算贪心解法的总距离
  best_dis = CalcTourDistance(best_solution);
  printf("Init greedy solution, best_dis %.4f\n", best_dis);
}

// 输出最终 best_solution 的结果及耗时
void PrintMyResults(int64_t start_us) {
  double ans = CalcTourDistance(best_solution);
  int64_t used_us = GetUs() - start_us;
  printf("For %d node Euclid TSP, calc_ans %.4f, best_dis %.4f, time %.2f s\n",
         n, ans, best_dis, used_us / 1000000.0);
}

// 随机两个下标 i、j 并尝试交换它们
void Opt2PointExchange(int *curr_solution, double *curr_dis,
                       double temperature) {
  int i = rand() % n;
  int j = rand() % n;
  if (i == j) {
    return;
  }
  if (i > j) {
    SwapInt(&i, &j);
  }

  double new_dis = *curr_dis;
  int pre_i = (i - 1 + n) % n;
  int nxt_i = (i + 1) % n;
  int pre_j = (j - 1 + n) % n;
  int nxt_j = (j + 1) % n;
  // 一头一尾
  if (i == 0 && j == n - 1) {
    new_dis -= GetDistance(curr_solution[pre_j], curr_solution[j]);
    new_dis -= GetDistance(curr_solution[j], curr_solution[i]);
    new_dis -= GetDistance(curr_solution[i], curr_solution[nxt_i]);
    new_dis += GetDistance(curr_solution[pre_j], curr_solution[i]);
    new_dis += GetDistance(curr_solution[i], curr_solution[j]);
    new_dis += GetDistance(curr_solution[j], curr_solution[nxt_i]);
  } else if (nxt_i != j) { // 不相连
    new_dis -= GetDistance(curr_solution[pre_i], curr_solution[i]);
    new_dis -= GetDistance(curr_solution[i], curr_solution[nxt_i]);
    new_dis -= GetDistance(curr_solution[pre_j], curr_solution[j]);
    new_dis -= GetDistance(curr_solution[j], curr_solution[nxt_j]);
    new_dis += GetDistance(curr_solution[pre_i], curr_solution[j]);
    new_dis += GetDistance(curr_solution[j], curr_solution[nxt_i]);
    new_dis += GetDistance(curr_solution[pre_j], curr_solution[i]);
    new_dis += GetDistance(curr_solution[i], curr_solution[nxt_j]);
  } else { // 相连
    new_dis -= GetDistance(curr_solution[pre_i], curr_solution[i]);
    new_dis -= GetDistance(curr_solution[i], curr_solution[j]);
    new_dis -= GetDistance(curr_solution[j], curr_solution[nxt_j]);
    new_dis += GetDistance(curr_solution[pre_i], curr_solution[j]);
    new_dis += GetDistance(curr_solution[j], curr_solution[i]);
    new_dis += GetDistance(curr_solution[i], curr_solution[nxt_j]);
  }

  // Metropolis 准则：若 Δt<0 则接受 S′ 作为新的当前解 S
  // 否则以概率 exp(-Δt/T) 接受 S′ 作为新的当前解S
  double delta_dis = new_dis - *curr_dis;
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / temperature))) {
    SwapInt(&curr_solution[i], &curr_solution[j]);
    *curr_dis = new_dis;
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}

// 尝试随机翻转一个子序列
// 具体而言，有头、中、尾翻转
void OptRangeReverse(int *curr_solution, double *curr_dis, double temperature) {
  // 生成期望长度为 n/3 的线段
  int idx_start = rand() % n;
  int idx_end = rand() % n;
  if (idx_start > idx_end) {
    SwapInt(&idx_start, &idx_end);
  }
  // 长度为 1 或 n 均没必要翻转
  if (idx_end - idx_start + 1 == n || idx_start == idx_end) {
    return;
  }
  // 限制操作长度，降低单次操作耗时
  if (n > 1200) { // magic
    int wanna_len = rand() % 600 + 3;
    idx_start = rand() % (n - wanna_len + 1);
    idx_end = idx_start + wanna_len - 1;
  }
  // 概率地决定是翻转最前面、最后面，还是中间
  {
    int rnd = rand() % 1000;
    int len = idx_end - idx_start + 1;
    if (rnd < 200) { // 20% 概率翻转开头
      idx_start = 0;
      idx_end = len - 1;
    } else if (rnd < 400) { // 20% 概率翻转结尾
      idx_start = n - len;
      idx_end = n - 1;
    } // 60% 概率翻转中间
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

  // Metropolis 准则：若 Δt<0 则接受 S′ 作为新的当前解 S
  // 否则以概率 exp(-Δt/T) 接受 S′ 作为新的当前解S
  double delta_dis = new_dis - *curr_dis;
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / temperature))) {
    int swapi = idx_start;
    int swapj = idx_end;
    while (swapi < swapj) {
      SwapInt(&curr_solution[swapi], &curr_solution[swapj]);
      swapi++;
      swapj--;
    }
    *curr_dis = new_dis;
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}

// 尝试将一个元素随机插到另一个元素前/后
// 该操作相比 range 翻转更符合 ATSP 局部性原则
// 此外，该操作计算 new_dis 为 O(1) 的，因此不需要限制操作长度
// 毕竟只有更改被接受了才会产生 O(n) 的调整
void OptInsert(int *curr_solution, double *curr_dis, double temperature) {
  int idx_i = rand() % n;
  int idx_j = rand() % n;
  if (idx_i > idx_j) {
    SwapInt(&idx_i, &idx_j);
  }
  // 长度为 1 或 n 均没必要插入（shift）
  if (idx_j - idx_i + 1 == n || idx_i == idx_j) {
    return;
  }

  double new_dis = *curr_dis;
  int is_qian = rand() % 2;
  if (is_qian) { // is_qian 为 1 则把 idx_i 挪到 idx_j 后头
    int pre_i = (idx_i - 1 + n) % n;
    int nxt_i = (idx_i + 1) % n;
    int nxt_j = (idx_j + 1) % n;
    new_dis -= GetDistance(curr_solution[pre_i], curr_solution[idx_i]);
    new_dis -= GetDistance(curr_solution[idx_i], curr_solution[nxt_i]);
    new_dis -= GetDistance(curr_solution[idx_j], curr_solution[nxt_j]);
    new_dis += GetDistance(curr_solution[pre_i], curr_solution[nxt_i]);
    new_dis += GetDistance(curr_solution[idx_j], curr_solution[idx_i]);
    new_dis += GetDistance(curr_solution[idx_i], curr_solution[nxt_j]);
  } else { // is_qian 为 0 则把 idx_j 挪到 idx_i 前头
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

  // Metropolis 准则：若 Δt<0 则接受 S′ 作为新的当前解 S
  // 否则以概率 exp(-Δt/T) 接受 S′ 作为新的当前解S
  double delta_dis = new_dis - *curr_dis;
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / temperature))) {
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
    *curr_dis = new_dis;
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}

// 有序交换两个区间
void OptSwapRange(int *curr_solution, double *curr_dis, double temperature) {
  // 生成 [idx_a, idx_b] 和 [idx_c, idx_d] 两个区间
  int idx_a;
  int idx_b;
  int idx_c;
  int idx_d;
  {
    int max_range_len = (n - 1) / 2;
    int range_len = rand() % max_range_len + 1;
    // 两个中心端点都在 [len-1, n-len] 之间，即 len-1+rand{0,n-2len+1}
    idx_b = rand() % (n - 2 * range_len + 1) + range_len - 1;
    idx_c = rand() % (n - 2 * range_len + 1) + range_len - 1;
    if (idx_b > idx_c) {
      SwapInt(&idx_b, &idx_c);
    }
    idx_a = idx_b - range_len + 1;
    idx_d = idx_c + range_len - 1;
  }

  double new_dis = *curr_dis;
  if (idx_a == 0 && idx_d == n - 1) { // 顶天立地
    int nxt_b = (idx_b + 1) % n;
    int prev_c = (idx_c - 1 + n) % n;
    new_dis -= GetDistance(curr_solution[idx_b], curr_solution[nxt_b]);
    new_dis -= GetDistance(curr_solution[prev_c], curr_solution[idx_c]);
    new_dis -= GetDistance(curr_solution[idx_d], curr_solution[idx_a]);
    new_dis += GetDistance(curr_solution[idx_d], curr_solution[nxt_b]);
    new_dis += GetDistance(curr_solution[prev_c], curr_solution[idx_a]);
    new_dis += GetDistance(curr_solution[idx_b], curr_solution[idx_c]);
  } else if (idx_b + 1 == idx_c) { // 区间贴贴
    int prev_a = (idx_a - 1 + n) % n;
    int nxt_d = (idx_d + 1) % n;
    new_dis -= GetDistance(curr_solution[prev_a], curr_solution[idx_a]);
    new_dis -= GetDistance(curr_solution[idx_b], curr_solution[idx_c]);
    new_dis -= GetDistance(curr_solution[idx_d], curr_solution[nxt_d]);
    new_dis += GetDistance(curr_solution[prev_a], curr_solution[idx_c]);
    new_dis += GetDistance(curr_solution[idx_d], curr_solution[idx_a]);
    new_dis += GetDistance(curr_solution[idx_b], curr_solution[nxt_d]);
  } else { // 河汉迢迢
    int prev_a = (idx_a - 1 + n) % n;
    int nxt_b = (idx_b + 1) % n;
    int prev_c = (idx_c - 1 + n) % n;
    int nxt_d = (idx_d + 1) % n;
    new_dis -= GetDistance(curr_solution[prev_a], curr_solution[idx_a]);
    new_dis -= GetDistance(curr_solution[idx_b], curr_solution[nxt_b]);
    new_dis -= GetDistance(curr_solution[prev_c], curr_solution[idx_c]);
    new_dis -= GetDistance(curr_solution[idx_d], curr_solution[nxt_d]);
    new_dis += GetDistance(curr_solution[prev_a], curr_solution[idx_c]);
    new_dis += GetDistance(curr_solution[idx_d], curr_solution[nxt_b]);
    new_dis += GetDistance(curr_solution[prev_c], curr_solution[idx_a]);
    new_dis += GetDistance(curr_solution[idx_b], curr_solution[nxt_d]);
  }

  // Metropolis 准则：若 Δt<0 则接受 S′ 作为新的当前解 S
  // 否则以概率 exp(-Δt/T) 接受 S′ 作为新的当前解S
  double delta_dis = new_dis - *curr_dis;
  if (delta_dis < 0.0 || ChkRand(exp(-delta_dis / temperature))) {
    int len = idx_b - idx_a + 1;
    for (int i = 0; i < len; i++) {
      SwapInt(&curr_solution[idx_a + i], &curr_solution[idx_c + i]);
    }
    *curr_dis = new_dis;
    if (*curr_dis < best_dis) {
      best_dis = *curr_dis;
      memcpy(best_solution, curr_solution, n * sizeof(int));
    }
  }
}

// 用 DP 精准解决小规模【开环】TSP
// 时间复杂度：2**n * n**2
// 空间复杂度：2**n * n
void ExactDp() {
  int n__2 = 1 << n;
  double **dp = (double **)malloc(n__2 * sizeof(double *));
  int **parent = (int **)malloc(n__2 * sizeof(int *));
  for (int i = 0; i < n__2; i++) {
    dp[i] = (double *)malloc(n * sizeof(double));
    parent[i] = (int *)malloc(n * sizeof(int));
    for (int j = 0; j < n; j++) {
      dp[i][j] = DBL_MAX; // 初始化为无穷大
      parent[i][j] = -1;  // 初始化为无
    }
  }

  // 从 0 开始
  dp[1 << 0][0] = 0.0;
  for (int S = 1; S < n__2; S++) {
    for (int i = 0; i < n; i++) {
      if (!(S & (1 << i))) {
        continue;
      }
      for (int j = 0; j < n; j++) {
        if (i == j || !(S & (1 << j))) {
          continue;
        }
        int prev = S & ~(1 << i);
        double dis = dp[prev][j] + GetDistance(j, i);
        if (dis < dp[S][i]) {
          dp[S][i] = dis;
          parent[S][i] = j;
        }
      }
    }
  }

  // 恢复答案
  best_dis = DBL_MAX;
  int lst = -1;
  int S = n__2 - 1;
  for (int i = 0; i < n; i++) {
    if (dp[S][i] < best_dis) {
      best_dis = dp[S][i];
      lst = i;
    }
  }
  best_solution[n - 1] = lst;
  for (int i = n - 2; i >= 0; i--) {
    best_solution[i] = parent[S][lst];
    S &= ~(1 << lst);
    lst = best_solution[i];
  }

  for (int i = 0; i < n__2; i++) {
    free(dp[i]);
    free(parent[i]);
  }
  free(dp);
  free(parent);
}

// 模拟退火解决 ATSP 问题 主函数
void SimulatedAnnealing() {
  int64_t start_us = GetUs();

  if (n <= 15) {
    ExactDp(); // 解决的是【开环】TSP
    for (int i = 0; i < n; i++) {
      printf("%d ", best_solution[i]);
    }
    printf("as result. Note that this is for OTSP!\n");
    PrintMyResults(start_us);
    return;
  }

  int64_t should_end_us = start_us + SA_MAX_RUN_US;

  // 贪心初始化一个解
  InitSolution();

  // 操作时候所用的缓冲区
  int *curr_solution = (int *)malloc(n * sizeof(int));

  int sa_cnt = 0;
  double last_best_dis = best_dis;
  // 模拟退火
  while (GetUs() < should_end_us) {
    // 复制历史最好解，按这个退火
    // 随机一个解退火效果很差
    double curr_dis = best_dis;
    memcpy(curr_solution, best_solution, n * sizeof(int));

    // 初始化温度、结束温度
    double tolerant_dis = best_dis * SA_TOLERATE_RATE;
    // exp(-Δt/T)=p -> -Δt/t = ln(p) -> t = -Δt/ln(p)
    double temperature = -fabs(tolerant_dis) / log(INIT_ACCEPT_P);
    double min_temperature = temperature * MIN_SA_COOL_RATIO;
    do {
      int rnd = rand() % 1000; // 决定操作概率
      if (rnd < 300) {
        // 30% 概率交换两点
        Opt2PointExchange(curr_solution, &curr_dis, temperature);
      } else if (rnd < 300) { // 25%或者10%？
        // 0% 概率反转区间。这个操作可能不太适合 ATSP
        // 额，感觉不写他更好。。。
        OptRangeReverse(curr_solution, &curr_dis, temperature);
      } else if (rnd < 1000) {
        // 70% 概率随机插入
        // 这个操作好像有点牛逼。设置成30%/10%/60%的时候效果挺好
        // 好像把它设置成 100% 都行
        OptInsert(curr_solution, &curr_dis, temperature);
      } else {
        // 0% 概率交换两个区间
        // 额，感觉不写他更好。。。
        OptSwapRange(curr_solution, &curr_dis, temperature);
      }
      temperature *= SA_COOL_RATE;
      sa_cnt++;
      if (best_dis < last_best_dis - 0.5) {
        // printf("after #%d sa, best_dis %.4f\n", sa_cnt, best_dis);
        last_best_dis = best_dis;
      }
    } while (temperature > min_temperature);
  }

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