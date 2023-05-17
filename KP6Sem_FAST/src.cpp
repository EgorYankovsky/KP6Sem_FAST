#include "Header.h"

namespace solution
{
   double R1R1(double r, double r0, double r1)
   {
      return r * (r1 - r) / (r1 - r0) * (r1 - r) / (r1 - r0);
   }

   double R1R2(double r, double r0, double r1)
   {
      return r * (r1 - r) / (r1 - r0) * (r - r0) / (r1 - r0);
   }

   double R2R2(double r, double r0, double r1)
   {
      return r * (r - r0) / (r1 - r0) * (r - r0) / (r1 - r0);
   }

   double dR1dR1(double r, double r0, double r1)
   {
      return -1.0 * r / (r1 - r0) * (-1.0 / (r1 - r0));
   }

   double dR1dR2(double r, double r0, double r1)
   {
      return -1.0 / (r1 - r0) * r / (r1 - r0);
   }
   
   double dR2dR2(double r, double r0, double r1)
   {
      return r / (r1 - r0) * 1.0 / (r1 - r0);
   }

   double Z1Z1(double z, double z0, double z1)
   {
      return (z1 - z) / (z1 - z0) * (z1 - z) / (z1 - z0);
   }

   double Z1Z2(double z, double z0, double z1)
   {
      return (z1 - z) / (z1 - z0) * (z - z0) / (z1 - z0);
   }

   double Z2Z2(double z, double z0, double z1)
   {
      return (z - z0) / (z1 - z0) * (z - z0) / (z1 - z0);
   }

   double dZ1dZ1(double z, double z0, double z1)
   {
      return -1.0 / (z1 - z0) * (-1.0 / (z1 - z0));
   }

   double dZ1dZ2(double z, double z0, double z1)
   {
      return -1.0 / (z1 - z0) * 1.0 / (z1 - z0);
   }

   double dZ2dZ2(double z, double z0, double z1)
   {
      return 1.0 / (z1 - z0) * 1.0 / (z1 - z0);
   }

   double Integ(basisF funct, double a, double b)
   {
      double value = 0.0;

      vector<double> w = { 0.0271524594117541,
                           0.0622535239386479,
                           0.0951585116824928,
                           0.1246289712555339,
                           0.1495959888165767,
                           0.1691565193950025,
                           0.1826034150449236,
                           0.1894506104550685,
                           0.1894506104550685,
                           0.1826034150449236,
                           0.1691565193950025,
                           0.1495959888165767,
                           0.1246289712555339,
                           0.0951585116824928,
                           0.0622535239386479,
                           0.0271524594117541
      };
      vector<double> x = { -0.9894009349916499,
                           -0.9445750230732326,
                           -0.8656312023878318,
                           -0.7554044083550030,
                           -0.6178762444026438,
                           -0.4580167776572274,
                           -0.2816035507792589,
                           -0.0950125098376374,
                           0.0950125098376374,
                           0.2816035507792589,
                           0.4580167776572274,
                           0.6178762444026438,
                           0.7554044083550030,
                           0.8656312023878318,
                           0.9445750230732326,
                           0.9894009349916499,

      };

      for (int i = 0; i < 16; i++)
         value += w[i] * funct(0.5 * (b - a) * x[i] + 0.5 * (b + a), a, b);

      return 0.5 * (b - a) * value;
   }

   void Difur::TimeGenerator(double kt, size_t numTimes, double t0, double t1)
   {
      ofstream fout("time.txt");
      fout << numTimes + 1 << endl;

      double sum = 0;
      for (size_t i(0); i < numTimes; i++)
         sum += pow(kt, i);
      double th0 = (t1 - t0) / sum;

      double t = t0;
      for (size_t i(0); i <= numTimes; i++)
      {
         fout << t << endl;
         t += th0 * pow(kt, i);
      }
      fout.close();
   }

   void Difur::Grid_Generator()
   {
      ifstream in("boundaries.txt");
      in >> r1 >> r2 >> kr >> z1 >> z2 >> kz >> Nr >> Nz;
      in.close();

      Ne = Nr * Nz;
      hr.resize(Nr + 1);
      hz.resize(Nz + 1);
      double sum = 0;
      for (int j = 0; j <= Nr - 1; j++)
         sum += pow(kr, j);
      double w = (r2 - r1) / sum;
      hr[0] = r1;
      double nstu;
      for (int i = 0; i < Nr; i++)
      {
         hr[i + 1] = hr[i] + w * pow(kr, i);
      }
      sum = 0;
      for (int j = 0; j <= Nz - 1; j++)
         sum += pow(kz, j);
      w = (z2 - z1) / sum;
      hz[0] = z1;
      for (int i = 0; i < Nz; i++)
      {
         hz[i + 1] = hz[i] + w * pow(kz, i);
      }
      int Nn = (Nr + 1) * (Nz + 1);
      RZ.resize(Nn);
      ofstream out("RZ.txt");
      out << Nn << '\n';
      for (int q = 0; q < Nn; q++)
      {
         RZ[q].first = hr[q % (Nr + 1)];
         RZ[q].second = hz[q / (Nr + 1)];
         out << RZ[q].first << ' ' << RZ[q].second << '\n';
      }
      out.close();
      elems.resize(Ne);
      int o = 1, u = 0;
      for (int k = 0; k < Nz; k++)
      {
         for (int j = 0; j < Nr; j++)
         {
            if (o % (Nr + 1) == 0)
            {
               o++;
            }
            elems[u].resize(5);
            elems[u][0] = o;
            elems[u][1] = o + 1;
            elems[u][2] = o + Nr + 1;
            elems[u][3] = o + Nr + 2;
            elems[u][4] = 1;
            o++;
            u++;
         }
      }
      out.open("elems.txt");
      out << elems.size() << '\n';
      for (int k = 0; k < Ne; k++)
         out << elems[k][0] << ' ' << elems[k][1] <<
         ' ' << elems[k][2] << ' ' << elems[k][3] <<
         ' ' << elems[k][4] << endl;
      out.close();
   }

   void Difur::Read()
   {
      ifstream fin("RZ.txt");
      fin >> Nn;
      P.resize(Nn);
      RZ.resize(Nn);
      for (int i = 0; i < Nn; i++)
         fin >> RZ[i].first >> RZ[i].second;
      fin.close();

      fin.open("elems.txt");
      fin >> Ne;
      elems.resize(Ne);
      for (int i = 0; i < Ne; i++)
      {
         elems[i].resize(5);
         fin >> elems[i][0] >> elems[i][1] >> elems[i][2] >> elems[i][3] >> elems[i][4];
      }
      fin.close();
      
      int N;
      fin.open("mat.txt");
      fin >> N;
      mat.resize(N);
      for (int i = 0; i < N; i++)
         fin >> mat[i].first >> mat[i].second;
      fin.close();

      fin.open("b1.txt");
      fin >> N;
      b1.resize(N);
      for (int i = 0; i < N; i++)
         fin >> b1[i].first;
      fin.close();

      fin.open("time.txt");
      fin >> N;
      t.resize(N);
      for (size_t i(0); i < N; i++)
         fin >> t[i];
      fin.close();

      q_h.resize(t.size());
      for (vector<double> &vct : q_h)
         vct.resize(RZ.size());
   }

   void Difur::Portrait(Matrix &matrix)
   {
      vector<vector<int>> adjacency;
      adjacency.resize(Nn);
      for (int i = 0; i < elems.size(); i++)
      {
         int n = adjacency[elems[i][1] - 1].size();
         adjacency[elems[i][1] - 1].resize(n + 1);
         adjacency[elems[i][1] - 1][n] = elems[i][0];

         n = adjacency[elems[i][2] - 1].size();
         adjacency[elems[i][2] - 1].resize(n + 2);
         adjacency[elems[i][2] - 1][n] = elems[i][0];
         adjacency[elems[i][2] - 1][n + 1] = elems[i][1];

         n = adjacency[elems[i][3] - 1].size();
         adjacency[elems[i][3] - 1].resize(n + 3);
         adjacency[elems[i][3] - 1][n] = elems[i][0];
         adjacency[elems[i][3] - 1][n + 1] = elems[i][1];
         adjacency[elems[i][3] - 1][n + 2] = elems[i][2];
      }

      for (int i = 0; i < Nn; i++)
      {
         sort(adjacency[i].begin(), adjacency[i].end());
         auto last = unique(adjacency[i].begin(), adjacency[i].end());
         adjacency[i].erase(last, adjacency[i].end());
      }

      matrix.di.resize(Nn);
      B.resize(Nn);
      matrix.ig.resize(Nn + 1);
      matrix.ig[0] = matrix.ig[1] = 1;
      for (int i = 2; i < Nn + 1; i++)
         matrix.ig[i] = matrix.ig[i - 1] + adjacency[i - 1].size();
      for (int i = 1; i < Nn; i++)
         matrix.jg.insert(matrix.jg.end(), adjacency[i].begin(), adjacency[i].end());
      matrix.ggu.resize(matrix.ig.back() - 1);
      matrix.ggl.resize(matrix.ig.back() - 1);

      //for (vector<int> row : adjacency)
      //{
      //   for (int elem : row)
      //      cout << elem << " ";
      //   cout << endl;
      //}
   }

   void Difur::localG(Matrix& matrix)
   {
      vector<vector<double>> G;
      G.resize(4);
      for (auto& row : G)
         row.resize(4);
      for (int i = 0; i < Ne; i++)
      {
         double rk = RZ[elems[i][0] - 1].first,
            _r0 = RZ[elems[i][0] - 1].first,
            _r1 = RZ[elems[i][3] - 1].first,
            _z0 = RZ[elems[i][0] - 1].second,
            _z1 = RZ[elems[i][3] - 1].second,
            l = mat[elems[i][4] - 1].first;

         G[0][0] = l * (Integ(dR1dR1, _r0, _r1) * Integ(Z1Z1, _z0, _z1) + Integ(R1R1, _r0, _r1) * Integ(dZ1dZ1, _z0, _z1));
         G[2][2] = l * (Integ(dR1dR1, _r0, _r1) * Integ(Z2Z2, _z0, _z1) + Integ(R1R1, _r0, _r1) * Integ(dZ2dZ2, _z0, _z1));
         G[0][1] = G[1][0] = l * (Integ(dR1dR2, _r0, _r1) * Integ(Z1Z1, _z0, _z1) + Integ(R1R2, _r0, _r1) * Integ(dZ1dZ1, _z0, _z1));
         G[0][2] = G[2][0] = l * (Integ(dR1dR1, _r0, _r1) * Integ(Z1Z2, _z0, _z1) + Integ(R1R1, _r0, _r1) * Integ(dZ1dZ2, _z0, _z1));
         G[0][3] = G[1][2] = G[2][1] = G[3][0] = l * (Integ(dR1dR2, _r0, _r1) * Integ(Z1Z2, _z0, _z1) + Integ(R1R2, _r0, _r1) * Integ(dZ1dZ2, _z0, _z1));
         G[1][1] = l * (Integ(dR2dR2, _r0, _r1) * Integ(Z1Z1, _z0, _z1) + Integ(R2R2, _r0, _r1) * Integ(dZ1dZ1, _z0, _z1));
         G[3][3] = l * (Integ(dR2dR2, _r0, _r1) * Integ(Z2Z2, _z0, _z1) + Integ(R2R2, _r0, _r1) * Integ(dZ2dZ2, _z0, _z1));
         G[1][3] = G[3][1] = l * (Integ(dR2dR2, _r0, _r1) * Integ(Z1Z2, _z0, _z1) + Integ(R2R2, _r0, _r1) * Integ(dZ1dZ2, _z0, _z1));
         G[2][3] = G[3][2] = l * (Integ(dR1dR2, _r0, _r1) * Integ(Z2Z2, _z0, _z1) + Integ(R1R2, _r0, _r1) * Integ(dZ2dZ2, _z0, _z1));

         for (int j = 0; j < 4; j++)
         {
            matrix.di[elems[i][j] - 1] += G[j][j];
         }

         int column, line;
         int n = 0;
         for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
            {
               n = 0;
               if (elems[i][j] > elems[i][k])
               {
                  line = elems[i][j] - 1;
                  column = elems[i][k];
                  while (matrix.jg[matrix.ig[line] - 1 + n] != column)
                     n++;
                  matrix.ggl[matrix.ig[line] - 1 + n] += G[j][k];
               }
               else if (elems[i][j] < elems[i][k])
               {
                  line = elems[i][k] - 1;
                  column = elems[i][j];
                  while (matrix.jg[matrix.ig[line] - 1 + n] != column)
                     n++;
                  matrix.ggu[matrix.ig[line] - 1 + n] += G[k][j];
               }
            }
      }

   }

   void Difur::localM(Matrix& matrix)
   {
      vector<vector<double>> C;

      C.resize(4);

      for (auto& row : C)
         row.resize(4);

      for (int i = 0; i < Ne; i++)
      {

         double rk = RZ[elems[i][0] - 1].first,
            _r0 = RZ[elems[i][0] - 1].first,
            _r1 = RZ[elems[i][3] - 1].first,
            _z0 = RZ[elems[i][0] - 1].second,
            _z1 = RZ[elems[i][3] - 1].second,
            g = mat[elems[i][4] - 1].second;

         C[0][0] = Integ(R1R1, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[2][2] = Integ(R1R1, _r0, _r1) * Integ(Z2Z2, _z0, _z1);
         C[0][1] = C[1][0] = Integ(R1R2, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[0][2] = C[2][0] = Integ(R1R1, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[0][3] = C[1][2] = C[2][1] = C[3][0] = Integ(R1R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[1][1] = Integ(R2R2, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[3][3] = Integ(R2R2, _r0, _r1) * Integ(Z2Z2, _z0, _z1);
         C[1][3] = C[3][1] = Integ(R2R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[2][3] = C[3][2] = Integ(R1R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);

         for (int j = 0; j < 4; j++)
         {
            matrix.di[elems[i][j] - 1] += g * C[j][j];
         }

         int column, line;
         int n = 0;
         for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
            {
               n = 0;
               if (elems[i][j] > elems[i][k])
               {
                  line = elems[i][j] - 1;
                  column = elems[i][k];
                  while (matrix.jg[matrix.ig[line] - 1 + n] != column)
                     n++;
                  matrix.ggl[matrix.ig[line] - 1 + n] += g * C[j][k];
               }
               else if (elems[i][j] < elems[i][k])
               {
                  line = elems[i][k] - 1;
                  column = elems[i][j];
                  while (matrix.jg[matrix.ig[line] - 1 + n] != column)
                     n++;
                  matrix.ggu[matrix.ig[line] - 1 + n] += g * C[k][j];
               }
            }
      }

   }

   void Difur::localVR(double time)
   {
      vector<vector<double>> C;

      C.resize(4);

      for (auto& row : C)
         row.resize(4);

      for (int i = 0; i < Ne; i++)
      {

         double rk = RZ[elems[i][0] - 1].first,
            _r0 = RZ[elems[i][0] - 1].first,
            _r1 = RZ[elems[i][3] - 1].first,
            _z0 = RZ[elems[i][0] - 1].second,
            _z1 = RZ[elems[i][3] - 1].second,
            g = mat[elems[i][4] - 1].second;

         C[0][0] = Integ(R1R1, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[2][2] = Integ(R1R1, _r0, _r1) * Integ(Z2Z2, _z0, _z1);
         C[0][1] = C[1][0] = Integ(R1R2, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[0][2] = C[2][0] = Integ(R1R1, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[0][3] = C[1][2] = C[2][1] = C[3][0] = Integ(R1R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[1][1] = Integ(R2R2, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[3][3] = Integ(R2R2, _r0, _r1) * Integ(Z2Z2, _z0, _z1);
         C[1][3] = C[3][1] = Integ(R2R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[2][3] = C[3][2] = Integ(R1R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);


         vector<double> b;// ��������� ������ ������ �����
         b.resize(4);

         for (int k = 0; k < 4; k++) //��������� ��������� ������
            for (int j = 0; j < 4; j++)
               b[k] += C[k][j] * Func(RZ[elems[i][k] - 1].first, RZ[elems[i][k] - 1].second, time);


         for (int j = 0; j < 4; j++)//��������� ���������� ������ � ������� ��������� ���������� �������
            B[elems[i][j] - 1] += b[j];
      }
   }

   void Difur::ImplicitScheme(int layer)
   {

      globalMatrix.clear_matrix();
      for (int i = 0; i < B.size(); i++)
         B[i] = 0;

      double t_ = t[layer];
      double t_1 = t[layer - 1];
      double dt = 0.0;

      if (layer == 1)
      {
         dt = t_ - t_1;
         globalMatrix.copy_matrix(globalMatrixM);
         globalMatrix.multiply_by_coef(1.0 / dt);
         globalMatrix.sum_matrix(globalMatrixG);

         localVR(t_);

         vector<double> temp;
         temp = globalMatrixM.multiply_matrix_by_vector(q_h[layer - 1], RZ.size());
         for (int i = 0; i < B.size(); i++)
            B[i] += 1.0 / dt * temp[i];
      }
      else
      {
         double t_2 = t[layer - 2];

         dt = t_ - t_2;
         double dt0 = t_ - t_1,
            dt1 = t_1 - t_2;

         globalMatrix.copy_matrix(globalMatrixM);
         globalMatrix.multiply_by_coef((dt + dt0) / (dt * dt0));
         globalMatrix.sum_matrix(globalMatrixG);

         localVR(t_);

         vector<double> temp;
         temp = globalMatrixM.multiply_matrix_by_vector(q_h[layer - 2], RZ.size());
         for (int i = 0; i < B.size(); i++)
            B[i] -= dt0 / (dt1 * dt) * temp[i];

         temp = globalMatrixM.multiply_matrix_by_vector(q_h[layer - 1], RZ.size());
         for (int i = 0; i < B.size(); i++)
            B[i] += dt / (dt1 * dt0) * temp[i];
      }
   }

   void Difur::output()
   {
      //for (int i = 0; i < q_h.size(); i++)
      //{
      //   for (int j = 0; j < q_h[i].size(); j++)
      //      cout << scientific << setprecision(15) << q_h[i][j] << endl;
      //   cout << endl;
      //}

      //for (int i = 0; i < q_h.size(); i++)
      //{
      //      cout << scientific << setprecision(7) << abs(U(XY[4].first, XY[4].second, t[i]) - q_h[i][4]) / U(XY[4].first, XY[4].second, t[i]) << endl;
      //}

      // Абсолютная погрешность
      //for (int i = 0; i < q_h.size(); i++)
      //{
      //   for (int j = 0; j < q_h[i].size(); j++)
      //      switch (j)
      //      {
      //      case 5:
      //      case 6:
      //      case 9:
      //      case 10:
      //      {
      //         cout << scientific << setprecision(7) << abs(FuncReal(RZ[j].first, RZ[j].second, t[i]) - q_h[i][j]) << endl;
      //         break;
      //      }
      //      default:
      //         break;
      //      }
      //   cout << endl;
      //}
      
      // Относительная погрешность
      for (int i = 0; i < q_h.size(); i++)
      {
         for (int j = 0; j < q_h[i].size(); j++)
            switch (j)
            {
            case 5:
            case 6:
            case 9:
            case 10:
            {
               cout << scientific << setprecision(7) << abs(FuncReal(RZ[j].first, RZ[j].second, t[i]) - q_h[i][j]) /
                  FuncReal(RZ[j].first, RZ[j].second, t[i]) << endl;
               break;
            }
            default:
               break;
            }
         cout << endl;
      }
      double sum = 0;

      for (int i = 0; i < q_h.size(); i++)
      {
            for (int j = 0; j < q_h[i].size(); j++)
               switch (j)
               {
               case 5:
               case 6:
               case 9:
               case 10:
               {
                  sum += abs(FuncReal(RZ[j].first, RZ[j].second, t[i]) - q_h[i][j]);
                  break;
               }
               default:
                  break;
               }
      }
      sum /= q_h.size();
      cout << scientific << setprecision(16) << sum << endl;
   }
   /*
   void Difur::GaMbo()
   {
      vector<vector<double>> G;
      vector<vector<double>> C;
      vector<double> b;
      G.resize(4); 
      C.resize(4); 

      for (auto& row : G)
         row.resize(4);

      for (auto& row : C)
         row.resize(4);

      b.resize(4);
      

      for (int i = 0; i < Ne; i++)
      {
         for (auto& elem : b)
            elem = 0;
         
         double rk = RZ[elems[i][0] - 1].first,
               _r0 = RZ[elems[i][0] - 1].first,
               _r1 = RZ[elems[i][3] - 1].first,
               _z0 = RZ[elems[i][0] - 1].second,
               _z1 = RZ[elems[i][3] - 1].second,
                 l = mat[elems[i][4] - 1].first,
                 g = mat[elems[i][4] - 1].second;
         
         G[0][0] =           l * (Integ(dR1dR1, _r0, _r1) * Integ(Z1Z1, _z0, _z1) + Integ(R1R1, _r0, _r1) * Integ(dZ1dZ1, _z0, _z1));
         G[2][2] =           l * (Integ(dR1dR1, _r0, _r1) * Integ(Z2Z2, _z0, _z1) + Integ(R1R1, _r0, _r1) * Integ(dZ2dZ2, _z0, _z1));
         G[0][1] = G[1][0] = l * (Integ(dR1dR2, _r0, _r1) * Integ(Z1Z1, _z0, _z1) + Integ(R1R2, _r0, _r1) * Integ(dZ1dZ1, _z0, _z1));
         G[0][2] = G[2][0] = l * (Integ(dR1dR1, _r0, _r1) * Integ(Z1Z2, _z0, _z1) + Integ(R1R1, _r0, _r1) * Integ(dZ1dZ2, _z0, _z1));
         G[0][3] = G[1][2] = G[2][1] = G[3][0] = l * (Integ(dR1dR2, _r0, _r1) * Integ(Z1Z2, _z0, _z1) + Integ(R1R2, _r0, _r1) * Integ(dZ1dZ2, _z0, _z1));
         G[1][1] =           l * (Integ(dR2dR2, _r0, _r1) * Integ(Z1Z1, _z0, _z1) + Integ(R2R2, _r0, _r1) * Integ(dZ1dZ1, _z0, _z1));
         G[3][3] =           l * (Integ(dR2dR2, _r0, _r1) * Integ(Z2Z2, _z0, _z1) + Integ(R2R2, _r0, _r1) * Integ(dZ2dZ2, _z0, _z1));
         G[1][3] = G[3][1] = l * (Integ(dR2dR2, _r0, _r1) * Integ(Z1Z2, _z0, _z1) + Integ(R2R2, _r0, _r1) * Integ(dZ1dZ2, _z0, _z1));
         G[2][3] = G[3][2] = l * (Integ(dR1dR2, _r0, _r1) * Integ(Z2Z2, _z0, _z1) + Integ(R1R2, _r0, _r1) * Integ(dZ2dZ2, _z0, _z1));

         C[0][0] =           Integ(R1R1, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[2][2] =           Integ(R1R1, _r0, _r1) * Integ(Z2Z2, _z0, _z1);
         C[0][1] = C[1][0] = Integ(R1R2, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[0][2] = C[2][0] = Integ(R1R1, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[0][3] = C[1][2] = C[2][1] = C[3][0] = Integ(R1R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[1][1] =           Integ(R2R2, _r0, _r1) * Integ(Z1Z1, _z0, _z1);
         C[3][3] =           Integ(R2R2, _r0, _r1) * Integ(Z2Z2, _z0, _z1);
         C[1][3] = C[3][1] = Integ(R2R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         C[2][3] = C[3][2] = Integ(R1R2, _r0, _r1) * Integ(Z1Z2, _z0, _z1);
         
         for (int k = 0; k < 4; k++)   
            for (int j = 0; j < 4; j++)
               b[k] += C[k][j] * Func(RZ[elems[i][j] - 1].first, RZ[elems[i][j] - 1].second, 0);

         for (int j = 0; j < 4; j++)
         {
            B[elems[i][j] - 1] += b[j];
            di[elems[i][j] - 1] += G[j][j] + g * C[j][j];
         }

         int column, line;
         int n = 0;
         for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
            {
               n = 0;
               if (elems[i][j] > elems[i][k])
               {
                  line = elems[i][j] - 1;
                  column = elems[i][k];
                  while (jg[ig[line] - 1 + n] != column)
                     n++;
                  ggl[ig[line] - 1 + n] += G[j][k] + g *  C[j][k];
               }
               else if (elems[i][j] < elems[i][k])
               {
                  line = elems[i][k] - 1;
                  column = elems[i][j];
                  while (jg[ig[line] - 1 + n] != column)
                     n++;
                  ggu[ig[line] - 1 + n] += G[k][j] + g * C[k][j];
               }
            }
      }
   }
   */
   
   void Difur::Boundary_conditions(double time)
   {
      for (int i = 0; i < b1.size(); i++)
      {
         double n = b1[i].first;
         globalMatrix.di[n - 1] = 1;
         for (int j = 0; j < globalMatrix.ig[n] - globalMatrix.ig[n - 1]; j++)
            globalMatrix.ggl[globalMatrix.ig[n - 1] + j - 1] = 0;

         for (int i = n; i < Nn; i++)
         {
            int k = -1;
            int j = 0;
            do
            {
               if (globalMatrix.jg[globalMatrix.ig[i] - 1 + j] == n)
                  k = j;
               j++;
            } while (j < globalMatrix.ig[i + 1] - globalMatrix.ig[i]);
            if (k != -1)
               globalMatrix.ggu[globalMatrix.ig[i] - 1 + k] = 0;
         }
         B[n - 1] = FuncReal(RZ[n - 1].first, RZ[n - 1].second, time);
      }
   }

   void Difur::LU()
   {
      d.resize(globalMatrix.di.size());
      gl.resize(globalMatrix.ggl.size());
      gu.resize(globalMatrix.ggu.size());
      int i0, i1, j, kk, ll, ll_1, lll, kkk;
      double sd, sl, su;
      for (int i = 0; i < Nn; i++)
      {
         i0 = globalMatrix.ig[i] - 1;
         i1 = globalMatrix.ig[i + 1] - 1;
         sd = 0;
         for (int jj = i0; jj < i1; jj++)
         {
            j = globalMatrix.jg[jj] - 1;
            kk = i0;
            ll = globalMatrix.ig[j] - 1;
            ll_1 = globalMatrix.ig[j + 1] - 1;
            sl = 0;
            su = 0;
            while (ll < ll_1 && kk < jj)
            {
               lll = globalMatrix.jg[ll] - 1;
               kkk = globalMatrix.jg[kk] - 1;

               if (lll == kkk)
               {
                  su += gl[ll] * gu[kk];
                  sl += gu[ll] * gl[kk];
                  kk++;
                  ll++;
               }
               else
               {
                  if (kkk < lll)
                     kk++;
                  else
                     ll++;
               }
            }
            gl[jj] = (globalMatrix.ggl[jj] - sl);
            gu[jj] = (globalMatrix.ggu[jj] - su) / d[j];
            sd += gu[jj] * gl[jj];
         }
         d[i] = globalMatrix.di[i] - sd;
      }

      for (auto &elem : P)
         elem = 0;
   }

   void Difur::InitiateConditions(double time)
   {
      for (int i = 0; i < RZ.size(); i++)
         q_h[time][i] = FuncReal(RZ[i].first, RZ[i].second, t[time]);
   }

   void Difur::mult_A(vector<double>& vect, vector<double>& res)
   {
      int i0, i1, j;
      for (int i = 0; i < Nn; i++)
      {
         i0 = globalMatrix.ig[i] - 1;
         i1 = globalMatrix.ig[i + 1] - 1;
         res[i] = globalMatrix.di[i] * vect[i];
         for (int k = i0; k < i1; k++)
         {
            j = globalMatrix.jg[k];
            res[i] += globalMatrix.ggl[k] * vect[j - 1];
            res[j - 1] += globalMatrix.ggu[k] * vect[i];
         }
      }
   }

   void Difur::Direct(vector<double>& vec, vector<double>& res)
   {
      for (int i = 0; i < Nn; i++)
      {
         res[i] = vec[i];
         double s = 0;
         int j;
         for (int k = globalMatrix.ig[i] - 1; k < globalMatrix.ig[i + 1] - 1; k++)
         {
            j = globalMatrix.jg[k] - 1;
            s += res[j] * gl[k];
         }
         res[i] -= s;
         res[i] /= d[i];
      }
   }

   void Difur::Reverse(vector<double>& vec, vector<double>& res)
   {
      for (int i = 0; i < Nn; i++)
         res[i] = vec[i];
      for (int i = Nn - 1; i >= 0; i--)
      {
         for (int k = globalMatrix.ig[i + 1] - 2; k >= globalMatrix.ig[i] - 1; k--)
         {
            int j = globalMatrix.jg[k] - 1;
            res[j] -= res[i] * gu[k];
         }
      }
   }

   double norm_vector(int n, vector<double>& vec)
   {
      double sum = 0;
      for (auto& elem : vec)
         sum += elem * elem;
      return sqrt(sum);
   }

   double scal_prod(vector<double>& a, vector<double>& b)
   {
      double sum = 0;
      for (int i = 0; i < a.size(); i++)
         sum += a[i] * b[i];
      return sum;
   }

   void mult_coef(vector<double>& vec, double k, vector<double>& res, int n)
   {
      for (int i = 0; i < n; i++)
         res[i] = k * vec[i];
   }

   void sum_vector(vector<double>& a, vector<double>& b, vector<double>& res, int n)
   {
      for (int i = 0; i < n; i++)
         res[i] = a[i] + b[i];
   }

   void Difur::LOS()
   {
      vector<double> r, p, z, LAUr, L, res;
      r.resize(Nn);
      p.resize(Nn);
      z.resize(Nn);
      LAUr.resize(Nn);
      L.resize(Nn);
      res.resize(Nn);
      double a, b, norm;
      int k;
      mult_A(P, res);
      for (int i = 0; i < Nn; i++)
         r[i] = B[i] - res[i];
      Direct(r, r);
      Reverse(r, z);
      mult_A(z, res);
      Direct(res, p);
      double normV = norm_vector(Nn, B);
      for (k = 0; (norm_vector(Nn, r) / normV) > eps && k < max_iter; k++)
      {
         norm = scal_prod(p, p);
         a = scal_prod(p, r) / norm;
         mult_coef(z, a, res, Nn);
         sum_vector(P, res, P, Nn);
         mult_coef(p, -a, res, Nn);
         sum_vector(r, res, r, Nn);
         Reverse(r, res);
         mult_A(res, LAUr);
         Direct(LAUr, LAUr);
         b = -scal_prod(p, LAUr) / norm;
         mult_coef(z, b, z, Nn);
         sum_vector(res, z, z, Nn);
         mult_coef(p, b, p, Nn);
         sum_vector(p, LAUr, p, Nn);
      }
   }
}