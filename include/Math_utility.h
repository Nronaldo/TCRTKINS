#ifndef _MATH_UTILITY_H
#define _MATH_UTILITY_H
#include <vector>
#include <math.h>
#include <cmath>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <iostream>
#include "Load_Configuration.h"
using namespace std;
template <typename T>	//模板需要在头文件中直接定义

T dot(const std::vector <T> &v1, int indexs1, const std::vector <T> &v2, int indexs2,int length) {
	double c = 0.0;
	int i = indexs1;
	int j = indexs2;
	for (int count = 0; count < length; i++, count++, j++) {
		c += v1[i] * v2[j];
	}
	return c;
}
template <typename T>
T norm_vector(const std::vector <T> &v, int indexs,int length){
	return sqrt(dot(v, indexs, v, indexs, length));
}
template <typename T>
void divide_vector( std::vector <T> &v, int indexs, int length,T domination) {
	for (int i = indexs, count = 0; count < length; i++, count++)
		v[i] /= domination;

	return;
}
double norm_array(const double * a,int n);
int getbitu(const char *buff, int pos, int len);
void setbitu(char *buff, int pos, int len, int value);
void arraymattofile(const double *A, int rows, int columns, int sr, int sc, int er, int ec, char mode,const char filename[]);
void matrixtofile(const Eigen::MatrixXd &mat,std::string matname,char mode);
void filetovector(Eigen::VectorXd &mat, int rows, string name);
void filetomatrix(Eigen::MatrixXd &mat,int rows,int columns,string name);
int lambda(int nb,int solnum,const Eigen::VectorXd &y,const Eigen::MatrixXd & Qb,Eigen::MatrixXd &b,vector<double> &s);
int LD(int n, const Eigen::MatrixXd & Q, Eigen::MatrixXd & L, Eigen::VectorXd &D);
void reduction(int n, Eigen::MatrixXd & L, Eigen::VectorXd &D, Eigen::MatrixXd &Z);
void gauss(int n, Eigen::MatrixXd & L, Eigen::MatrixXd &Z, int i, int j);
void myperm(int n, Eigen::MatrixXd & L, Eigen::VectorXd &D, int j, double del, Eigen::MatrixXd&Z);
void swap(double &a, double &b);
int mysearch(int n, int m, const Eigen::MatrixXd &L, const Eigen::VectorXd &D,
	const Eigen::VectorXd &zs, Eigen::MatrixXd &zn, vector<double> &s);
//int mysolve(const Eigen::MatrixXd & Z, const Eigen::MatrixXd & E, int n, int m, Eigen::MatrixXd & F);
double sgn(double x);
double sqvar(double covar);
#endif
#pragma once
