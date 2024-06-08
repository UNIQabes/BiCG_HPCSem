#include <stdio.h>
#include <vector>
#include <iostream>

using namespace std;

#define N 1000
#define GAMMA 0.9
#define ITERLIMIT 1000

// CRS形式の行列
struct CRSMat
{
	int rowNum;
	int colNum;
	vector<double> val;
	vector<int> column_index;
	vector<int> rowTop_ind;

	CRSMat(int rowNum,
		   int colNum,
		   vector<double> val,
		   vector<int> column_index,
		   vector<int> rowTop_ind) : rowNum(rowNum), colNum(colNum), val(val), column_index(column_index), rowTop_ind(rowTop_ind)
	{
		// もし行数または、データ数に不整合があった場合
		if (rowTop_ind.size() != rowNum || val.size() != column_index.size())
		{
			// エラー落ちする
			exit(1);
		}
	}
};
template <typename T>
void printVector(string title, vector<T> vec)
{
	printf("%s:", title.c_str());
	for (int i = 0; i < vec.size(); i++)
	{
		cout << " ";
		cout << vec[i];
	}
	cout << "\n";
}

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

vector<double> MatvecProduct(CRSMat mat, vector<double> vec)
{
	int valNum = mat.val.size();
	int rowNum = mat.rowNum;
	int colNum = mat.colNum;
	vector<double> retValue(colNum, 0);
	for (int r = 0; r < rowNum; r++)
	{

		int startCRSindex = mat.rowTop_ind[r];
		int nextCRSindex = r == rowNum - 1 ? valNum : mat.rowTop_ind[r + 1];
		for (int crs_ind = startCRSindex; crs_ind < nextCRSindex; crs_ind++)
		{
			retValue[r] += mat.val[crs_ind] * vec[mat.column_index[crs_ind]];
		}
	}
	return retValue;
}

vector<double> TransMatvecProduct(CRSMat TransMat, vector<double> vec)
{
	int valNum = TransMat.val.size();
	int rowNum = TransMat.colNum;
	int colNum = TransMat.rowNum;
	vector<double> retValue(colNum, 0);
	for (int c = 0; c < colNum; c++)
	{
		int startCRSindex = TransMat.rowTop_ind[c];
		int nextCRSindex = c == colNum - 1 ? valNum : TransMat.rowTop_ind[c + 1];
		for (int crs_ind = startCRSindex; crs_ind < nextCRSindex; crs_ind++)
		{
			int r = TransMat.column_index[crs_ind];
			double TransMatVal_rc = TransMat.val[crs_ind];
			retValue[r] += TransMatVal_rc * vec[c];
		}
		// printf("/");
	}
	// printf("\n");
	/*
	for (int i = 0; i < 5; i++)
	{
		printf("%lf", retValue[i]);
	}
	printf("\n");
	*/
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
	/*
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
	*/
	CRSMat A_CRS(N, N, val_CRS, col_ind_CRS, row_ptr_CRS);

	vector<double> initialX(N, 0);
	vector<double> vecFilled1(N, 1);
	vector<double> b = MatvecProduct(A_CRS, vecFilled1);

	vector<double> x_k = initialX;
	vector<double> r_k = VecAddition(b, vec_numtimes(-1, MatvecProduct(A_CRS, x_k))); // xの初期値が0ベクトルなので、bでも良い
	vector<double> rstar_k = r_k;
	vector<double> p_k = r_k;
	vector<double> pstar_k = rstar_k;

	int counter = 0;
	while (counter < ITERLIMIT && vec_norm(r_k) / vec_norm(b) >= 1e-12)
	{
		printf("\"%d\", \"%.15lf\"\n", counter, vec_norm(r_k) / vec_norm(b));

		vector<double> q_k = MatvecProduct(A_CRS, p_k);
		vector<double> qstar_k = TransMatvecProduct(A_CRS, pstar_k);

		double alpha_k = vecDot(rstar_k, r_k) / vecDot(pstar_k, q_k);
		// printf("alpha:%lf\n", alpha_k);

		x_k = VecAddition(x_k, vec_numtimes(alpha_k, p_k));

		vector<double> r_k1 = VecAddition(r_k, vec_numtimes(-alpha_k, q_k));
		vector<double> rstar_k1 = VecAddition(rstar_k, vec_numtimes(-alpha_k, qstar_k));

		double beta_k = vecDot(rstar_k1, r_k1) / vecDot(rstar_k, r_k);
		// printf("beta:%lf\n", beta_k);
		p_k = VecAddition(r_k1, vec_numtimes(beta_k, p_k));
		pstar_k = VecAddition(rstar_k1, vec_numtimes(beta_k, pstar_k));

		// q,alpha,betaは次のループで使われない。pとxはその場で更新できる。そのため、rのみ最後に更新
		r_k = r_k1;
		rstar_k = rstar_k1;

		counter++;
	}
	printf("\"%d\", \"%.15lf\"\n", counter, vec_norm(r_k) / vec_norm(b));
}