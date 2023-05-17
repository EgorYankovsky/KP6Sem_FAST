#pragma once
#include<iostream>
#include<fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

namespace solution
{
   typedef double basisF(double t, double t0, double t1);
   
   double Integ(basisF funct, double a, double b);

   struct Matrix {
      vector <int> ig, jg;
      vector <double> ggl, ggu, di; 

      void sum_elems_to_matrix(int i, int j, double elems) {
         if (i == j)
            di[i - 1] += elems;
         else {
            if (i > j) {
               for (int k = ig[i - 1]; k < ig[i]; k++)
                  if (jg[k - 1] == j)
                     ggl[k - 1] += elems;
            }
            else {
               for (int k = ig[j - 1]; k < ig[j]; k++)
                  if (jg[k - 1] == i)
                     ggu[k - 1] += elems;
            }
         }
      }

      void multiply_by_coef(double coef) {
         for (int i = 0; i < di.size(); i++)
            di[i] *= coef;
         for (int i = 0; i < ggl.size(); i++) {
            ggl[i] *= coef;
            ggu[i] *= coef;
         }
      }

      void copy_matrix(Matrix& C) {
         for (int i = 0; i < di.size(); i++)
            di[i] = C.di[i];
         for (int i = 0; i < ggl.size(); i++) {
            ggl[i] = C.ggl[i];
            ggu[i] = C.ggu[i];
         }
      }

      void clear_matrix() {
         for (int i = 0; i < di.size(); i++)
            di[i] = 0;
         for (int i = 0; i < ggl.size(); i++) {
            ggl[i] = 0;
            ggu[i] = 0;
         }
      }

      vector<double> multiply_matrix_by_vector(vector<double> v, int size) {
         vector<double> res;
         res.resize(v.size());
         for (int i = 0; i < size; i++) {
            int i0 = ig[i] - 1;
            int i1 = ig[i + 1] - 1;
            res[i] = di[i] * v[i];
            for (int k = i0; k < i1; k++) {
               int j = jg[k];
               res[i] += ggl[k] * v[j - 1];
               res[j - 1] += ggu[k] * v[i];
            }
         }
         return res;
      }

      void sum_matrix(Matrix& C) {
         for (int i = 0; i < di.size(); i++)
            di[i] += C.di[i];
         for (int i = 0; i < ggl.size(); i++) {
            ggl[i] += C.ggl[i];
            ggu[i] += C.ggu[i];
         }
      }
   };

   class Difur
   {
   private:
      vector<vector<double>> b2, elems; // ������ ������� ������� (����� ��������, ��� ��������� ����� � ��������), �������� ��������� � ������
      vector<pair<double, double>> RZ, mat, b1; // ���������� ����� (r,z), �������� (������, �����), ������ ������� ������� (������ ����� ����, �������� �������)
      vector<double> hr, hz, d, gl, gu, B;
      double r1, r2, z1, z2, kr, kz, eps = 1e-15; // ������� �� R � Z, ������������ ����������� �� R � Z, �������
      size_t Nr, Nz, Ne, Nn, max_iter = 10'000; // ���-�� ��������� �� R � �� Z, ���-�� ���������, ���������� �����, ������ ���-�� ��������
   public:
      Matrix globalMatrix, globalMatrixM, globalMatrixG;
      vector<double> t, P;
      vector<vector<double>> q_h;
      void TimeGenerator(double k, size_t numTimes, double t0, double t1);
      void Grid_Generator(); // ���������� �����
      void Read(); // ��������� � �������������� ������
      void Portrait(Matrix &matrix); // ��� �������� �������� ig
      void GaMbo(); // ���������� �������
      void InitiateConditions(double time);
      double Func(double r, double z, double t) //������� ������ �����
      {
         return 2.0 * t;
      }
      double FuncReal(double r, double z, double t)
      {
         return t * t;
      }
      void Boundary_conditions(double time); // ���� ������� �������
      void LU();
      void LOS();
      void output();
      void localG(Matrix& matrix);
      void localM(Matrix& matrix);
      void localVR(double time);
      void ImplicitScheme(int layer);

      // ���� ��������������� ������� ��� ������� ���� ����� LU+���
      void mult_A(vector<double>& vect, vector<double>& res);
      void Direct(vector<double>& vec, vector<double>& res);
      void Reverse(vector<double>& vec, vector<double>& res);
   };
}