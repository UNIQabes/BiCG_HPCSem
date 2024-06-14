#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000
#define p 0.1
#define MAX_L 1000

int main(int argc, char *argv[])
{
	printf("\"繰り返し回数\",\"BiCG_C\"\n");

	// CRS形式の行列の初期化--------------------
	int i, j;
	double val[3 * N - 2];
	int col_ind[3 * N - 2];
	int row_ptr[N + 1];
	val[0] = 2;
	val[1] = 1;
	for (i = 1; i < N - 1; i++)
	{
		val[3 * i - 1] = p;
		val[3 * i] = 2;
		val[3 * i + 1] = 1;
	}
	val[3 * N - 4] = p;
	val[3 * N - 3] = 2;

	col_ind[0] = 0;
	col_ind[1] = 1;
	for (i = 1; i < N - 1; i++)
	{
		col_ind[3 * i - 1] = i - 1;
		col_ind[3 * i] = i;
		col_ind[3 * i + 1] = i + 1;
	}
	col_ind[3 * N - 4] = N - 2;
	col_ind[3 * N - 3] = N - 1;

	row_ptr[0] = 0;
	row_ptr[1] = 2;
	for (i = 2; i < N; i++)
	{
		row_ptr[i] = 3 * i - 1;
	}
	// CRS形式を扱うループの処理を簡略化するために、最後の部分にval配列のサイズを入れる
	row_ptr[N] = 3 * N - 2;
	// (終)CRS形式の行列の初期化--------------------

	// b(=Ax)の値を算出
	double b[N];
	for (i = 0; i < N; i++)
	{
		b[i] = 0;
		for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
		{
			b[i] += val[j];
		}
	}
	// (終)b(=Ax)の値を算出

	double r0[N], r1[N], p0[N], p1[N], q0[N], x0[N], x1[N];
	double R0[N], R1[N], P0[N], P1[N], Q0[N];
	double b_norm, r_norm;
	b_norm = 0;
	for (i = 0; i < N; i++)
	{
		x0[i] = 0;
		r0[i] = b[i];
		R0[i] = r0[i];
		p0[i] = r0[i];
		P0[i] = R0[i];
		b_norm += b[i] * b[i];
	}
	int count = 0;
	double alfa, beta, tmp1, tmp2;
	for (int l = 0; l < MAX_L; l++)
	{

		count++;
		r_norm = 0;
		for (i = 0; i < N; i++)
		{
			r_norm += r0[i] * r0[i];
		}
		printf("%d, %.15lf\n", count - 1, sqrt(r_norm) / sqrt(b_norm));

		if (sqrt(r_norm / b_norm) < 1e-12)
		{
			break;
		}

		// q=Apk----------------
		for (i = 0; i < N; i++)
		{
			q0[i] = 0;
			Q0[i] = 0;
		}
		for (i = 0; i < N; i++)
		{
			for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
			{
				q0[i] += val[j] * p0[col_ind[j]];
				Q0[col_ind[j]] += val[j] * P0[i];
			}
		}
		// (終)q=Apk----------------

		/*
		printf("pstar_k:");
		for (int i = 0; i < 5; i++)
		{
			printf("%lf ", P0[i]);
		}
		printf("\n");
		*/

		/*
		printf("qstar_k:");
		for (int i = 0; i < 5; i++)
		{
			printf("%lf ", Q0[i]);
		}
		printf("\n");
		*/

		// alphak=(rstark,rk)/(pstark,qk)-----------------------
		tmp1 = 0;
		tmp2 = 0;
		for (i = 0; i < N; i++)
		{
			tmp1 += R0[i] * r0[i];
			tmp2 += P0[i] * q0[i];
		}
		alfa = tmp1 / tmp2;
		// printf("alpha:%lf\n", alfa);
		//   (終)alphak=(rstark,rk)/(pstark,qk)-----------------------

		// xk1=xk+alphak*pk  rk1=rk-alphak*qk-----------------------
		for (i = 0; i < N; i++)
		{
			x1[i] = x0[i] + alfa * p0[i];
			r1[i] = r0[i] - alfa * q0[i];
			R1[i] = R0[i] - alfa * Q0[i];
		}
		// (終)xk1=xk+alphak*pk  rk1=rk-alphak*qk-----------------------

		/*
		printf("rstar_k1:");
		for (int i = 0; i < 5; i++)
		{
			printf("%lf ", R1[i]);
		}
		printf("\n");
		*/

		// betak=(rstar_k1,r_k1)/(rstar_k,r_k) -------------------
		tmp1 = 0;
		tmp2 = 0;
		for (i = 0; i < N; i++)
		{
			tmp1 += R1[i] * r1[i];
			tmp2 += R0[i] * r0[i];
		}
		beta = tmp1 / tmp2;
		// printf("beta:%lf\n", beta);
		//   (終)betak=(rstar_k1,r_k1)/(rstar_k,r_k) -------------------

		// pk1=rk1+betak*pk -------------------
		for (i = 0; i < N; i++)
		{
			p1[i] = r1[i] + beta * p0[i];
			P1[i] = R1[i] + beta * P0[i];
		}
		// (終)pk1=rk1+betak*pk -------------------

		for (i = 0; i < N; i++)
		{
			r0[i] = r1[i];
			p0[i] = p1[i];
			x0[i] = x1[i];
			R0[i] = R1[i];
			P0[i] = P1[i];
		}
	}
	// printf("count = %d\n", count);
}