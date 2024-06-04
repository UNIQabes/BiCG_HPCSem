#include <stdio.h>
#include <vector>

using namespace std;

#define N 10
#define GAMMA 1

vector<double> vec_numtimes(double num, vector<double> vec)
{
	int vecSize = vec.size();
	vector<double> retValue(vecSize);
	for (int i = 0; i < vecSize; i++)
	{
		retValue[i] = num * vec[i];
	}
	return retValue;
}

double vecDot(vector<double> v1, vector<double> v2)
{
	int vecSize = v1.size();
	double retValue = 0;
	for (int i = 0; i < vecSize; i++)
	{
		retValue += v1[i] * v2[i];
	}
	return retValue;
}

double vec_norm(vector<double> vec)
{
	return sqrt(vecDot(vec, vec));
}

vector<double> VecAddition(vector<double> v1, vector<double> v2)
{
	int vecSize = v1.size();
	vector<double> retValue(vecSize, 0);
	for (int i = 0; i < vecSize; i++)
	{
		retValue[i] = v1[i] + v2[i];
	}
	return retValue;
}

vector<double> MatvecProduct(vector<double> val_CRS, vector<int> col_ind_CRS, vector<int> row_ptr_CRS, vector<double> vec)
{
	int valNum = val_CRS.size();
	int rowNum = row_ptr_CRS.size();
	int colNum = vec.size();
	vector<double> retValue(colNum, 0);
	for (int r = 0; r < N; r++)
	{
		int startCRSindex = row_ptr_CRS[r];
		int nextCRSindex = r == N - 1 ? val_CRS.size() : row_ptr_CRS[r + 1];
		for (int crs_ind = startCRSindex; crs_ind < nextCRSindex; crs_ind++)
		{
			retValue[r] += val_CRS[crs_ind] * vec[col_ind_CRS[crs_ind]];
		}
	}
	return retValue;
}

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
	vector<double> vecFilled1(N, 1);
	vector<double> b = MatvecProduct(val_CRS, col_ind_CRS, row_ptr_CRS, vecFilled1);

	vector<double> x_k = initialX;
	vector<double> r_k = VecAddition(b, MatvecProduct(val_CRS, col_ind_CRS, row_ptr_CRS, x_k)); // xの初期値が0ベクトルなので、bでも良い
	vector<double> p_k = r_k;

	/*
	for (int r = 0; r < N; r++)
	{
		printf(" %lf ", r_k[r]);
	}
	printf("\n");
	*/

	int counter = 0;
	while (vec_norm(r_k) / vec_norm(b) >= 1e-12)
	{
		vector<double> q_k = MatvecProduct(val_CRS, col_ind_CRS, row_ptr_CRS, p_k);
		double alpha_k = vecDot(r_k, r_k) / vecDot(p_k, q_k);

		x_k = VecAddition(x_k, vec_numtimes(alpha_k, p_k));

		vector<double> r_k1 = VecAddition(r_k, vec_numtimes(-alpha_k, q_k));
		double beta_k = vecDot(r_k1, r_k1) / vecDot(r_k, r_k);
		p_k = VecAddition(r_k1, vec_numtimes(beta_k, p_k));

		r_k = r_k1;
		counter++;
		if (counter % 1 == 0)
		{
			printf("%d:", counter);
			for (int r = 0; r < N; r++)
			{
				printf(" %4lf ", x_k[r]);
			}
			printf("\n");
		}
	}

	for (int r = 0; r < N; r++)
	{
		printf(" %4lf ", x_k[r]);
	}
	printf("\n");
}
