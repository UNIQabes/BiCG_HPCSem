#include <stdio.h>
#include <vector>

using namespace std;

#define N 10
#define GAMMA 0.5

int main()
{
	// 行列Aの初期化
	vector<double> val_CRS(3 * N - 2);
	vector<int> col_ind_CRS(3 * N - 2);
	vector<int> row_ptr_CRS(N);
	for (int r = 0; r < N; r++)
	{
		row_ptr_CRS[r] = 3 * r - 1;
		if (r == 0)
		{
			row_ptr_CRS[r] = 0;
		}
		for (int i = 0; i < 3; i++)
		{
			int c = r - 1 + i;
			if (c < 0 || c >= N)
			{
				continue;
			}
			int counter = 3 * r + i - 1;
			switch (i)
			{
			case 0:
				val_CRS[counter] = GAMMA;
				break;
			case 1:
				val_CRS[counter] = 2;
				break;
			case 2:
				val_CRS[counter] = 1;
				break;
			default:
				break;
			}
			col_ind_CRS[counter] = c;
		}
	}

	// 行列Aの内容を確認
	for (int r = 0; r < N; r++)
	{
		int nextval_CRSindex = row_ptr_CRS[r];
		int nextcolumn_CRSindex = val_CRS.size();
		if (r < N - 1)
		{
			nextcolumn_CRSindex = row_ptr_CRS[r + 1];
		}

		for (int c = 0; c < N; c++)
		{
			double dispNum = 0;
			if (nextval_CRSindex < nextcolumn_CRSindex && col_ind_CRS[nextval_CRSindex] == c)
			{
				dispNum = val_CRS[nextval_CRSindex];
				nextval_CRSindex++;
			}
			printf(" %4lf ", dispNum);
		}
		printf("\n");
	}

	vector<double> initialX(N, 0);
	vector<double> b(N, 1);
}
