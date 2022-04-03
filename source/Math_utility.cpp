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
		printf("��������\n");
		return;
	}
	else if (sc < 0 || ec < sc || ec >= columns) {
		printf("��������\n");
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
			F = (Z.inverse()).transpose()*E; /* F=inv(Zt)E */ //���任
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
	//		for (k = 0; k <= j; k++) // 0 ��j��
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
	while (j >= 0) {// ���������������˵������һ��ѭ����û���κ�����������Ҫ������
		if (j <= k)// ��¼��ǰ�Ѿ�����һ�н��и�˹�任�� 
			for (i = j + 1; i < n; i++) // j�� i �С����ҵ�����������أ����һ���ò�������
				gauss(n, L, Z, i, j);
		/*matrixtofile(L, "L", 'N');
		matrixtofile(D, "D", 'N');
		matrixtofile(Z, "Z", 'N');*/
		del = D(j) + L(j + 1 ,j) * L(j + 1,j) * D(j + 1);//���ʵ���Ͼ��Ǽ���ĸ��º��dj+1
		if (del + 1E-6 < D(j + 1)) { /* compared considering numerical error */// ������������������Ҳ���Ǹ������ݡ�
			/*matrixtofile(L, "L", 'N');
			matrixtofile(D, "D", 'N');
			matrixtofile(Z, "Z", 'N');*/
			myperm(n, L, D, j, del, Z);
			/*matrixtofile(L, "L", 'N');
			matrixtofile(D, "D", 'N');
			matrixtofile(Z, "Z", 'N');*/
			k = j;// ��¼��ǰ�Ѿ���˹�任���У���һ��ѭ�������п�ʼ���и�˹�任�� 
			j = n - 2;//ֻҪ����������Ҫ�ٴ����һ�п�ʼ�ж��Ƿ�Ҫ�����С�
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
	D(j) = eta * D(j + 1); // ����ʵ������Ϊ�˸���D
	D(j + 1) = del;
	// ���幫ʽ���Բο��ο�����1 2 ����matlab���롣��һ��
	for (k = 0; k <= j - 1; k++) {// ����ʵ������Ϊ�˸���L21��2��j-1�� ����������j��ʵ�൱�ڴ������j+1�У��� k��L21��ÿһ��
		a0 = L(j , k); //ԭʼL21 ÿһ�еĵ�һ��
		a1 = L(j + 1 ,k);//ԭʼL21 ÿһ�еĵڶ���
		L(j , k) = -L(j + 1 ,j) * a0 + a1;//����L21 ÿһ�еĵ�һ��
		L(j + 1,k) = eta * a0 + lam * a1;//����L21 ÿһ�еĵڶ���
	}
	//L(i:i+1,1:i-1) = [ -L(i+1,i) 1 ; eta lambda ] * L(i:i+1,1:i-1);
	L(j + 1 ,j) = lam;// ������Ϊ�˸���L22 ע�� ������˹�任��L22�ĶԽǱ�����1������Ҫ�仯
	for (k = j + 2; k < n; k++) // ����ʵ������Ϊ�˸���L32 n-(j+1)�� 2�� ������������
		swap(L(k,j), L(k,(j + 1)));
	for (k = 0; k < n; k++) // ˳������Z��
		swap(Z(k ,j), Z(k ,(j + 1)));

}
int mysearch(int n, int m, const Eigen::MatrixXd &L, const Eigen::VectorXd &D,
	const Eigen::VectorXd &zs, Eigen::MatrixXd &zn, vector<double> &s) {
	int i, j, k, c, nn = 0, imax = 0;
	double newdist, maxdist = 1E99, y;// maxdist �൱������ֵ����ѡ��ǰm����ѡ���ʱ�� maxdist Ҳ����£������ǵ�ǰ������ۼ�ֵ��ѡ��֮��maxdist��Ϊ��ǰm����ѡ���ж�Ӧ������ۼ�ֵ��
	// �˺������ҵ��ۼ�ֵС��maxdist�ĺ�ѡ�⣬����ԭ��ѡ�⣬�����µ�ǰm����ѡ���Ӧ��maxdist
	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
	Eigen::VectorXd dist = Eigen::VectorXd::Zero(n,1);
	Eigen::VectorXd zb = Eigen::VectorXd::Zero(n, 1);
	Eigen::VectorXd z = Eigen::VectorXd::Zero(n, 1);
	Eigen::VectorXd step = Eigen::VectorXd::Zero(n, 1);
	k = n - 1;// nԪ�ش����±�n-1 ��Ϊ��ͳһ��ʾ���ǻ���Ϊnά��
	dist[k] = 0.0;// ��Ž�ֹ����ά�ȵľ��롣ע�⣬��û���ϴ�level�Ĳ��� ��ʵ����ֻ��ӳ�˽�ֹ����һά�ȵ��ۼӾ���
	zb[k] = zs[k];// n level�� ��z �� ñz һ�¡�
	z[k] = round(zb[k]);// ȡ��ӽ�����
	y = zb[k] - z[k];// ����nά�Ȳ���
	step[k] = sgn(y);  /* step towards closest integer */// ����nά��Ѱ����һ�����ķ���
	for (c = 0; c < 10000; c++) {
		newdist = dist[k] + y * y / D[k];  /* newdist=sum(((z(j)-zb(j))^2/d(j))) */ // ��ֹ����ǰά�ȵľ�����ϵ�ǰ�Ĳ��࣬��������ľ��ǵ�ǰά�ȵ��ۼӾ��룬��û�и��µ�ǰά��dist[k]
		if (newdist < maxdist) {// ��ֹ����ǰά�ȵ��ۼӾ�����ȻС�������룬Ҫô��������һ��ά������������Ѿ�����ά��1���Ϳ��Ը��µ�ǰ��p����ѡ��
			/* Case 1: move down */
			if (k != 0) {//˵��û��ά��1
				dist[--k] = newdist;// ��ǰά�ȵ��ۼӾ��룬����һά�ȵ�dist,
				for (i = 0; i <= k; i++)
					S(k ,i) = S(k + 1,i) + (z[k + 1] - zb[k + 1])*L(k + 1, i);// �˴���Ϊ�˼���ο����׶��Ĺ�ʽ17 ��չ���ᷢ��һģһ����������Sn*n��Ϊ�˺ü����ۼ�
				zb[k] = zs[k] + S(k, k );//���㵱ǰά�ȵ� ��z��
				z[k] = round(zb[k]); /* next valid integer */ // ���㵱ǰά�ȵ��������
				y = zb[k] - z[k]; // ���㵱ǰά�ȵĲ���
				step[k] = sgn(y);// ���㵱ǰά�ȵ���һ��Ѱ���������򣬺Ͳ����й�ϵ
			}
			/* Case 2: store the found candidate and try next valid integer */
			else {//˵����ά��1
				if (nn < m) {  /* store the first m initial points */ // �����û�湻m����ѡ��
					if (nn == 0 || newdist > s[imax])// ��Ȼ��û�湻������ֻҪ�����º�ѡ���ۼӾ�����󣬾͸���һ�µ�ǰ������������
						imax = nn;
					for (i = 0; i < n; i++)
						zn(i, nn) = z[i];// ��ǰ��ѡ����ӵ�zn��
					s[nn++] = newdist;// ��ǰ��ѡ���Ӧ���ۼӾ��뱣�浽s��
				}
				else { // �Ѿ��湻m����ѡ��
					if (newdist < s[imax]) { // ����º�ѡ����ۼӾ����С��С�ڵ�ǰm����ѡ���е�����ۼӾ��롣
						for (i = 0; i < n; i++)
							zn(i ,imax) = z[i];// ���ۼӾ����Ӧ��������Ӧ�ĺ�ѡ����¡�
						s[imax] = newdist;// �ۼӾ������
						for (i = imax = 0; i < m; i++)// ͳ�Ƶ�ǰm����ѡ���������ۼӾ��룬��Ϊ�µ�maxdist 
							if (s[imax] < s[i])
								imax = i;
					}
					maxdist = s[imax];
				}
				z[0] += step[0]; /* next valid integer */ // ������һ����ҪѰ���µĺ�ѡ�⣬��ά��1��ʼ��Ѱ��ά��1����һ������������ע�⣬���ʱ���Ӧ��z�Ѿ�����
				y = zb[0] - z[0];// ���㵱ǰ������Ӧ�Ĳ���
				step[0] = -step[0] - sgn(step[0]);// ���Ѱ������һ������������Ҫ���Ѿ����µ�z��Ѱ�ҡ�
				// ����������⣬������µ��������������С��maxdist��ֱ�Ӿ���һ����ѡ�⣬
				// ��������㣬���ڵ�ǰ��ά������Ѱ�ҡ�Case 3
			}
		}
		/* Case 3: exit or move up */
		else {
			if (k == n - 1) // �ҵ���㣬���Ҳ���С��maxdist�ĺ�ѡ�⣬��ֹͣ�������ˡ�
				break;
			else {// ��������һlevelѰ��
				k++;  /* move up */
				z[k] += step[k];  /* next valid integer */ // ���ƣ�Ѱ�ҵ�ǰά�ȵ���һ�������ܡ����¡�
				y = zb[k] - z[k];// ���㵱ǰά�ȵ�ǰ������Ӧ�Ĳ���
				step[k] = -step[k] - sgn(step[k]);// ���Ѱ������һ������������Ҫ���Ѿ����µ�z��Ѱ�ҡ�
			}
		}
	}
	for (i = 0; i < m - 1; i++) { /* sort by s */ // m��ѡ���������Ž���Ŀ 
		for (j = i + 1; j < m; j++) { // ð�������
			if (s[i] < s[j]) // s���ۼӾ��롣
				continue;
			swap(s[i], s[j]);// ���ۼӾ����С�������򣬶�Ӧ��Ҳ����m���������Ž⵽����
			for (k = 0; k < n; k++) // ͬʱҲҪ�����Ž��ž���zn�е�һ�к�һ�н�����
				swap(zn(k ,i), zn(k ,j));// �������zn�е�Ԫ�ش������Ҿ������Ž⵽���⡣
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
