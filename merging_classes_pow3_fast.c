/*
 * merging_classes_pow3_fast.c
 * ============================
 * merging_classes_fast.c を (n, n+1) ペアから (m, 3^p*m+x) ペアに一般化した版。
 *
 * §13で見つかった「m と 3^p*m+x が同じステップ数で同じ値に合流する」現象
 * (diff は常に -p) について、§5b と同じ手法(合流が"証明可能"な剰余類
 * r (mod 2^k) の割合を調べる)を、より大きな k まで押して、
 * §5b で見つかった「収束せず対数的に増加し続ける」という挙動が
 * この power-of-3 ペアでも成り立つかどうかを確認する。
 *
 * ビルド:
 *   gcc -O3 -march=native -fopenmp merging_classes_pow3_fast.c -o merging_classes_pow3_fast
 *
 * 実行:
 *   ./merging_classes_pow3_fast <p> <x> <k_min> <k_max>
 *   例 (§13のp=2,x=1ケースをk=27まで押す):
 *     ./merging_classes_pow3_fast 2 1 3 27
 *
 * 注意: unsigned __int128 を使用。k=30程度までは安全のはず
 * (power_of_3_merging_classes_fast.py と同じロジックだが、C+OpenMPでさらに高速)。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef unsigned __int128 u128;

static u128 ipow(u128 base, int exp) {
    u128 r = 1;
    for (int i = 0; i < exp; i++) r *= base;
    return r;
}

static int collatz_traj(u128 n, u128 *traj, int max_steps) {
    int len = 0;
    traj[len++] = n;
    while (n != 1 && len < max_steps) {
        if (n % 2 == 0) n = n / 2;
        else n = 3 * n + 1;
        traj[len++] = n;
    }
    return len;
}

// mの軌道とn=c*m+xの軌道を辿り、最初に共通する値が現れる
// (m基準ステップ, n基準ステップ) を求める。見つからなければ0を返す。
static int find_merge(u128 m, u128 n, int max_steps, int *merge_a, int *merge_b) {
    static _Thread_local u128 traj_m[400];
    static _Thread_local u128 traj_n[400];
    int len_m = collatz_traj(m, traj_m, max_steps);
    int len_n = collatz_traj(n, traj_n, max_steps);

    for (int i = 0; i < len_m; i++) {
        for (int j = 0; j < len_n; j++) {
            if (traj_m[i] == traj_n[j]) {
                *merge_a = i;
                *merge_b = j;
                return 1;
            }
        }
    }
    return 0;
}

static int is_merging_class(u128 r, u128 modulus, u128 c, u128 x,
                             int num_samples, int max_steps) {
    int first_a = -1, first_b = -1;
    for (int s = 0; s < num_samples; s++) {
        u128 j = (u128)(s + 2);
        u128 m = r + j * modulus;
        if (m < 1) continue;
        u128 n = c * m + x;
        int a, b;
        if (!find_merge(m, n, max_steps, &a, &b)) {
            return 0;
        }
        if (first_a == -1) {
            first_a = a; first_b = b;
        } else if (a != first_a || b != first_b) {
            return 0;
        }
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s <p> <x> <k_min> <k_max>\n", argv[0]);
        return 1;
    }
    int p = atoi(argv[1]);
    long long x_arg = atoll(argv[2]);
    int k_min = atoi(argv[3]);
    int k_max = atoi(argv[4]);
    u128 c = ipow(3, p);
    u128 x = (u128)x_arg;
    int num_samples = 4;
    int max_steps = 300;

    printf("c=3^%d=%lld, x=%lld, threads=%d\n", p, (long long)c, x_arg, omp_get_max_threads());
    fflush(stdout);

    for (int k = k_min; k <= k_max; k++) {
        u128 modulus = (u128)1 << k;
        long long total = (long long)1 << k;
        long long count = 0;

        double t0 = omp_get_wtime();

        #pragma omp parallel for schedule(dynamic, 1024) reduction(+:count)
        for (long long r = 0; r < total; r++) {
            if (is_merging_class((u128)r, modulus, c, x, num_samples, max_steps)) {
                count++;
            }
        }

        double elapsed = omp_get_wtime() - t0;
        double frac = (double)count / (double)total;
        printf("k=%2d modulus=%lld : merging=%lld/%lld = %.4f%%   (%.1fs)\n",
               k, (long long)modulus, count, total, 100.0 * frac, elapsed);
        fflush(stdout);
    }

    return 0;
}
