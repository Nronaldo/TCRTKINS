#include "Math_utility.h"
#include "earth_utility.h"
double norm_array(const double * a, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i] * a[i];
	}
	return sqrt(sum);
}
int getbitu(const char *buff, int pos, int len) {
	int bits = 0;
	int i;
	for (i = pos; i < pos + len; i++)
		bits = (bits << 1) + ((buff[i / 8] >> (7 - i % 8)) & 1u);
	return bits;
}
void setbitu(char *buff, int pos, int len, int data) {
	int mask = 1u << (len - 1);
	int i;
	if (len <= 0 || 32 < len) return;
	for (i = pos; i < pos + len; i++, mask >>= 1) {
		if (data&mask)
			buff[i / 8] |= 1u << (7 - i % 8);
		else
			buff[i / 8] &= ~(1u << (7 - i % 8));
	}
	return;
}
void matrixtofile(const Eigen::MatrixXd &mat, std::string matname,char mode) {
	std::string filename = "../matrix/" + matname + ".txt";
	std::ofstream out(filename, std::ios::out);
	if(mode == 'N')
		out <<setprecision(16)<< mat << std::endl;
		//out << mat << std::endl;
	else if (mode == 'T') {
		Eigen::MatrixXd matT = mat.transpose();
		out << setprecision(16) << matT << std::endl;
		//out << matT << std::endl;
	}
	out.close();
	return;
}
void arraymattofile(const double *A, int rows, int columns, int sr, int sc, int er, int ec, char mode, const char f[]) {
	if (sr < 0 || er < sr || er >= rows) {
		printf("行数不对\n");
		return;
	}
	else if (sc < 0 || ec < sc || ec >= columns) {
		printf("列数不对\n");
		return;
	}
	else {
		if (mode == 'F') {
			char filename[20];
			strcpy_s(filename, "../matrix/");
			strcat_s(filename, f);
			strcat_s(filename, ".txt");
			//strcat(filename, "/0");
			FILE*fp = fopen(filename, "w");
			if (!fp)
				printf("open failed\n");
			for (int i = sr; i <= er; i++) {
				for (int j = sc; j <= ec; j++) {
					fprintf(fp, "%-10lf ", A[i + j * rows]);
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
		else {
			for (int i = sr; i <= er; i++) {
				for (int j = sc; j <= ec; j++) {
					printf("%-10lf ", A[i + j * rows]);
				}
				printf("\n");
			}
		}

	}

}
int lambda(int nb, int solnum, const Eigen::VectorXd &a, const Eigen::MatrixXd & Q, Eigen::MatrixXd &F, vector<double> &s) {
	int info = 0;
	if (nb <= 0 || solnum <= 0) return -1;
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(nb, nb);
	Eigen::VectorXd D = Eigen::VectorXd::Zero(nb);
	Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(nb,nb);
	Eigen::VectorXd z = Eigen::VectorXd::Zero(nb);
	Eigen::MatrixXd E = Eigen::MatrixXd::Zero(nb,solnum);
	if (numiterator == debug_stop) {
		//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}
	if (!(info = LD(nb, Q, L, D))) {
		//matrixtofile(L, "L", 'N');
		//matrixtofile(D, "D", 'N');
		if (numiterator == debug_stop) {
			//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
			int aaaa = 1;
		}
		reduction(nb, L, D, Z);
		/*matrixtofile(L, "L", 'N');
		matrixtofile(D, "D", 'N');
		matrixtofile(Z, "Z", 'N');*/
		z = Z.transpose()*a;
		
		/*if (numiterator == debug_stop) {
			matrixtofile(z, "z", 'N');
			int aaaa = 1;
		}*/
		if (numiterator == debug_stop) {
			//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
			int aaaa = 1;
		}
		if (!(info = mysearch(nb, solnum, L, D, z, E, s))) {
			//matrixtofile(E, "E", 'N');
			if (numiterator == debug_stop) {
				//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
				int aaaa = 1;
			}
			F = (Z.inverse()).transpose()*E; /* F=inv(Zt)E */ //反变换
			//matrixtofile(F, "F", 'N');
		}
	}
	return info;
}
int LD(int n, const Eigen::MatrixXd & Q, Eigen::MatrixXd & L, Eigen::VectorXd &D) {
	int i, j, k, info = 0;
	double a;
	Eigen:: MatrixXd A = Q;
	if (numiterator == debug_stop) {
		//arraymattofile(RTK_P, nx, nx, 0, 0, 8, 8, 'F', "P-a");
		int aaaa = 1;
	}
	for (i = n - 1; i >= 0; i--) {
		if ((D(i) = A(i, i)) <= 0.0) {
			info = -1;
			break;
		}
		a = sqrt(D(i));
		L.block(i, 0, 1, i+1) = A.block(i, 0, 1, i+1) / a;
		//matrixtofile(L, "L", 'N');
		for (j = 0; j <=i - 1; j++) {
			A.block(j, 0, 1, j+1) = A.block(j, 0, 1, j+1) - L.block(i, 0, 1, j+1)*L(i, j);
		}
		L.block(i, 0, 1, i+1) = L.block(i, 0, 1, i+1) / L(i, i);	
		//matrixtofile(L, "L", 'N');
	}
	//for (i = n - 1; i >= 0; i--) {
	//	if ((D(i) = A(i, i)) <= 0.0) {
	//		info = -1;
	//		break;
	//	}
	//	a = sqrt(D(i));
	//	for (j = 0; j <= i; j++)
	//		L(i,j) = A(i,j) / a;
	//	//matrixtofile(L, "L", 'N');
	//	for (j = 0; j <= i - 1; j++) {
	//		for (k = 0; k <= j; k++) // 0 到j列
	//			A(j,k) -= L(i,k) * L(i,j);
	//	}
	//	for (j = 0; j <= i; j++)
	//		L(i,j) /= L(i,i);//L(i,1:i) = L(i,1:i)/L(i,i);
	//	//matrixtofile(L, "L", 'N');
	//}
	if (info)
		cout << "LD factorization error" << endl;
	return info;
}
void reduction(int n, Eigen::MatrixXd & L, Eigen::VectorXd &D, Eigen::MatrixXd &Z) {
	int i, j, k;
	double del;
	j = n - 2;
	k = n - 2;
	while (j >= 0) {// 如果此条件不满足说明，上一个循环中没有任何两组数据需要重排列
		if (j <= k)// 记录当前已经到哪一列进行高斯变换了 
			for (i = j + 1; i < n; i++) // j列 i 行。从右到左逐列消相关，最后一列用不着消。
				gauss(n, L, Z, i, j);
		/*matrixtofile(L, "L", 'N');
		matrixtofile(D, "D", 'N');
		matrixtofile(Z, "Z", 'N');*/
		del = D(j) + L(j + 1 ,j) * L(j + 1,j) * D(j + 1);//这个实际上就是计算的更新后的dj+1
		if (del + 1E-6 < D(j + 1)) { /* compared considering numerical error */// 满足条件则发生交换，也就是更新数据。
			/*matrixtofile(L, "L", 'N');
			matrixtofile(D, "D", 'N');
			matrixtofile(Z, "Z", 'N');*/
			myperm(n, L, D, j, del, Z);
			/*matrixtofile(L, "L", 'N');
			matrixtofile(D, "D", 'N');
			matrixtofile(Z, "Z", 'N');*/
			k = j;// 记录当前已经高斯变换的列，下一次循环从这列开始进行高斯变换。 
			j = n - 2;//只要发生交换就要再从最后一列开始判断是否要重排列。
		}
		else j--;
	}


}
void gauss(int n, Eigen::MatrixXd & L, Eigen::MatrixXd& Z, int i, int j) {
	int k, mu;
	if ((mu = round(L(i , j ))) != 0) {// mu = round(L(j,i));if mu % if mu not equal to 0
		//for (k = i; k < n; k++)	//L(j:n,i) = L(j:n,i) - mu * L(j:n,j); 
			L.block(i,j,n-i,1) = L.block(i, j, n - i, 1)-(double)mu*L.block(i ,i ,n-i,1);
		//for (k = 0; k < n; k++)
			Z.block(0,j ,n ,1)= Z.block(0, j, n, 1)-(double)mu*Z.block(0 ,i, n,1);// iZt(:,j) = iZt(:,j) + mu * iZt(:,i);  %iZt is inv(Zt) matrix 
	}
}
void myperm(int n, Eigen::MatrixXd & L, Eigen::VectorXd &D, int j, double del, Eigen::MatrixXd&Z) {
	int k;
	double eta, lam, a0, a1;

	eta = D(j) / del;
	lam = D(j + 1) * L(j + 1, j) / del;
	D(j) = eta * D(j + 1); // 这里实际上是为了更新D
	D(j + 1) = del;
	// 具体公式可以参考参考文献1 2 或者matlab代码。都一样
	for (k = 0; k <= j - 1; k++) {// 这里实际上是为了更新L21（2行j-1列 （这里索引j其实相当于处理的是j+1列）） k是L21的每一列
		a0 = L(j , k); //原始L21 每一列的第一行
		a1 = L(j + 1 ,k);//原始L21 每一列的第二行
		L(j , k) = -L(j + 1 ,j) * a0 + a1;//更新L21 每一列的第一行
		L(j + 1,k) = eta * a0 + lam * a1;//更新L21 每一列的第二行
	}
	//L(i:i+1,1:i-1) = [ -L(i+1,i) 1 ; eta lambda ] * L(i:i+1,1:i-1);
	L(j + 1 ,j) = lam;// 这里是为了更新L22 注意 经过高斯变换后，L22的对角本就是1，不需要变化
	for (k = j + 2; k < n; k++) // 这里实际上是为了更新L32 n-(j+1)行 2列 交换两列数字
		swap(L(k,j), L(k,(j + 1)));
	for (k = 0; k < n; k++) // 顺带更新Z。
		swap(Z(k ,j), Z(k ,(j + 1)));

}
int mysearch(int n, int m, const Eigen::MatrixXd &L, const Eigen::VectorXd &D,
	const Eigen::VectorXd &zs, Eigen::MatrixXd &zn, vector<double> &s) {
	int i, j, k, c, nn = 0, imax = 0;
	double newdist, maxdist = 1E99, y;// maxdist 相当于门限值。在选够前m个候选解的时候 maxdist 也会更新，保持是当前的最大累加值，选够之后，maxdist改为当前m个候选解中对应的最大累加值。
	// 此后若能找到累加值小于maxdist的候选解，代替原候选解，并更新当前m个候选解对应的maxdist
	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
	Eigen::VectorXd dist = Eigen::VectorXd::Zero(n,1);
	Eigen::VectorXd zb = Eigen::VectorXd::Zero(n, 1);
	Eigen::VectorXd z = Eigen::VectorXd::Zero(n, 1);
	Eigen::VectorXd step = Eigen::VectorXd::Zero(n, 1);
	k = n - 1;// n元素存在下标n-1 ，为了统一表示我们还称为n维度
	dist[k] = 0.0;// 存放截止到此维度的距离。注意，还没算上此level的残余 ，实际上只反映了截止到上一维度的累加距离
	zb[k] = zs[k];// n level上 杠z 和 帽z 一致。
	z[k] = round(zb[k]);// 取最接近整数
	y = zb[k] - z[k];// 计算n维度残余
	step[k] = sgn(y);  /* step towards closest integer */// 计算n维度寻找下一整数的方向
	for (c = 0; c < 10000; c++) {
		newdist = dist[k] + y * y / D[k];  /* newdist=sum(((z(j)-zb(j))^2/d(j))) */ // 截止到当前维度的距离加上当前的残余，这样计算的就是当前维度的累加距离，但没有更新当前维度dist[k]
		if (newdist < maxdist) {// 截止到当前维度的累加距离仍然小于最大距离，要么继续往下一个维度搜索，如果已经到了维度1，就可以更新当前的p个候选解
			/* Case 1: move down */
			if (k != 0) {//说明没到维度1
				dist[--k] = newdist;// 当前维度的累加距离，是下一维度的dist,
				for (i = 0; i <= k; i++)
					S(k ,i) = S(k + 1,i) + (z[k + 1] - zb[k + 1])*L(k + 1, i);// 此处是为了计算参考文献二的公式17 ，展开会发现一模一样。声明了Sn*n是为了好计算累加
				zb[k] = zs[k] + S(k, k );//计算当前维度的 杠z。
				z[k] = round(zb[k]); /* next valid integer */ // 计算当前维度的最近整数
				y = zb[k] - z[k]; // 计算当前维度的残余
				step[k] = sgn(y);// 计算当前维度的下一个寻找整数方向，和残余有关系
			}
			/* Case 2: store the found candidate and try next valid integer */
			else {//说明到维度1
				if (nn < m) {  /* store the first m initial points */ // 如果还没存够m个候选解
					if (nn == 0 || newdist > s[imax])// 虽然还没存够，但是只要遇到新候选解累加距离更大，就更新一下当前的最大距离索引
						imax = nn;
					for (i = 0; i < n; i++)
						zn(i, nn) = z[i];// 当前候选解添加到zn中
					s[nn++] = newdist;// 当前候选解对应的累加距离保存到s中
				}
				else { // 已经存够m个候选解
					if (newdist < s[imax]) { // 如果新候选解的累加距离更小，小于当前m个候选解中的最大累加距离。
						for (i = 0; i < n; i++)
							zn(i ,imax) = z[i];// 将累加距离对应的索引对应的候选解更新。
						s[imax] = newdist;// 累加距离更新
						for (i = imax = 0; i < m; i++)// 统计当前m个候选解中最大的累加距离，作为新的maxdist 
							if (s[imax] < s[i])
								imax = i;
					}
					maxdist = s[imax];
				}
				z[0] += step[0]; /* next valid integer */ // 到了这一步，要寻找新的候选解，从维度1开始，寻找维度1的下一个最优整数，注意，这个时候对应的z已经跟新
				y = zb[0] - z[0];// 计算当前整数对应的残余
				step[0] = -step[0] - sgn(step[0]);// 如果寻找下下一个最优整数，要在已经更新的z上寻找。
				// 这里这样理解，如果更新的整数仍满足距离小于maxdist，直接就是一个候选解，
				// 如果不满足，会在当前的维度向上寻找。Case 3
			}
		}
		/* Case 3: exit or move up */
		else {
			if (k == n - 1) // 找到最顶层，还找不到小于maxdist的候选解，就停止不再找了。
				break;
			else {// 否则往上一level寻找
				k++;  /* move up */
				z[k] += step[k];  /* next valid integer */ // 类似，寻找当前维度的下一整数可能。更新。
				y = zb[k] - z[k];// 计算当前维度当前整数对应的残余
				step[k] = -step[k] - sgn(step[k]);// 如果寻找下下一个最优整数，要在已经更新的z上寻找。
			}
		}
	}
	for (i = 0; i < m - 1; i++) { /* sort by s */ // m是选出来的最优解数目 
		for (j = i + 1; j < m; j++) { // 冒泡排序把
			if (s[i] < s[j]) // s是累加距离。
				continue;
			swap(s[i], s[j]);// 让累加距离从小到大排序，对应的也就是m个解中最优解到最差解
			for (k = 0; k < n; k++) // 同时也要把最优解存放矩阵zn中的一列和一列交换。
				swap(zn(k ,i), zn(k ,j));// 这样最后zn中的元素从左往右就是最优解到最差解。
		}
	}

	if (c >= 10000) {
		cout<<"search loop count overflow"<<endl;
		return -1;
	}
	return 0;



}
double sgn(double x) { return x <= 0 ? -1.0 : 1.0; }
void swap(double &a, double &b) { double tmp; tmp = a; a = b; b = tmp; }
void filetovector(Eigen::VectorXd &vec, int rows, string name) {
	std::string filename = "../test/" + name + ".txt";
	std::ifstream fin(filename, std::ios::binary);
	//fin >> vec;
	for (int i = 0; i < rows; i++)
		fin >> vec(i);
	fin.close();
}
void filetomatrix(Eigen::MatrixXd &mat, int rows, int columns, string name) {
	std::string filename = "../test/" + name + ".txt";
	std::ifstream fin(filename, std::ios::binary);
	for (int i = 0; i < rows; i++) {
		for (int  j= 0; j < columns; j++) {
			fin >> mat(i, j);
		}
	}
	fin.close();
}
double sqvar(double covar) { return covar < 0.0 ? -sqrt(-covar) : sqrt(covar); }
