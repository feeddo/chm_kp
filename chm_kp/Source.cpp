#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

void computeLocalMatrixAndRHS( // ���������� ��������� ������� � ������ �����
   vector<vector<double>>& local_a, // ��������� ������� ��������
   vector<double>& local_b, // ��������� ������ ������ �����
   vector<double>& x, // ���� �����
   vector<double>& lambda,
   vector<double>& f, // ������ ����� ��
   uint32_t k, // ����� �������� ��������
   double gamma_k)
{
   double coeff_g, coeff_m;

   coeff_g = (lambda[k] + lambda[k + 1]) / (2 * (x[k + 1] - x[k]));
   coeff_m = gamma_k * (x[k + 1] - x[k]) / 6;

   local_a[0][0] = coeff_g + 2 * coeff_m;
   local_a[0][1] = coeff_m - coeff_g;
   local_a[1][0] = coeff_m - coeff_g;
   local_a[1][1] = coeff_g + 2 * coeff_m;

   local_b[0] = (x[k + 1] - x[k]) / 6 * (2 * f[k] + f[k + 1]);
   local_b[1] = (x[k + 1] - x[k]) / 6 * (f[k] + 2 * f[k + 1]);
}

void computeGlobalMatrixAndRHS( // ������ ���������� ������� � ������ �����
   vector<double>& di,
   vector<double>& ig,
   vector<double>& jg,
   vector<double>& ggl,
   vector<double>& ggu,
   vector<double>& b,
   vector<double>& x,
   vector<double>& lambda,
   vector<double>& gamma,
   vector<double>& f,
   uint32_t n)
{
   vector<vector<double> > local_a = vector<vector<double> >(2, vector<double>(2));
   vector<double> local_b(2, 0);

   jg.resize(n - 1);
   ggl.resize(n - 1);
   ggu.resize(n - 1);

   di[0] = 0;
   b[0] = 0;
   ig[0] = 1;

   for (uint32_t k = 0; k < n - 1; k++)
   {
      computeLocalMatrixAndRHS(local_a, local_b, x, lambda, f, k, gamma[k]);

      // ������������ ������ ���������� �������
      di[k] += local_a[0][0];
      di[k + 1] = local_a[1][1];
      ig[k + 1] = k + 1;
      jg[k] = k;
      ggl[k] = local_a[1][0];
      ggu[k] = local_a[0][1];
      b[k] += local_b[0];
      b[k + 1] = local_b[1];
   }
   ig[n] = n;
}

void includeBoundary(
   vector<double>& di,
   vector<double>& b,
   uint32_t n,
   double left_bound,
   double right_bound)
{
   double c = 10e11; // ������� ����� ��� ���������� ������������� ����. ����
   di[0] = c;
   di[n - 1] = c;
   b[0] = c * left_bound;
   b[n - 1] = c * right_bound;
}

void solveLLt( // ������� ���� � �������������� ���������� ���������
   vector<double>& di,
   vector<double>& ig,
   vector<double>& jg,
   vector<double>& ggl,
   vector<double>& ggu,
   vector<double>& b,
   vector<double>& u,
   uint32_t n)
{
   vector<double> ld(n, 0); // ������� ��������� ������ ����������� �������
   vector<double> ll(n - 1, 0); // ��������� ��� �������
   vector<double> y(n, 0);

   // ���������� LLt
   ld[0] = sqrt(di[0]);
   ll[0] = ggl[0] / ld[0];
   for (uint32_t i = 1; i < n - 1; i++)
   {
      ld[i] = sqrt(di[i] - ll[i - 1] * ll[i - 1]);
      ll[i] = ggl[i] / ld[i];
   }
   ld[n - 1] = sqrt(di[n - 1] - ll[n - 2] * ll[n - 2]);

   // ������ ��� (L y = b)
   y[0] = b[0] / ld[0];
   for (uint32_t i = 1; i < n; i++)
      y[i] = (b[i] - ll[i - 1] * y[i - 1]) / ld[i];

   // �������� ��� (Lt u = y)
   u[n - 1] = y[n - 1] / ld[n - 1];
   for (int i = n - 2; i >= 0; i--)
      u[i] = (y[i] - ll[i] * u[i + 1]) / ld[i];
}

double arbitrarySolution( // ���������� ������� � ������������ �����
   vector<double>& x,
   vector<double>& u,
   uint32_t n,
   double x_arbtr)
{
   uint32_t i;
   for (i = 0; i < n - 1; i++)
      if (x_arbtr < x[i + 1])
         break;
   return (u[i] * (x[i + 1] - x_arbtr) + u[i + 1] * (x_arbtr - x[i])) / (x[i + 1] - x[i]);
}

int main()
{
   vector<double> x; // ���� �����
   vector<double> lambda; // �����. �������� (�������� � �����)
   vector<double> gamma;
   vector<double> f; // ������ ����� ��
   double left_bound, right_bound; // ����. ������� �� �������� ����. �������
   uint32_t n; // ����� ����� �����, ����� �� = n - 1

   vector<double> di; // ��������� ���������� ������� ����
   vector<double> ig;
   vector<double> jg;
   vector<double> ggl; // ������ �����������
   vector<double> ggu; // ������� �����������
   vector<double> b; // ������ ������ ����� ����

   vector<double> u; // �������
   double x_arbtr; // ������������ ����� ��� ������ �������

   ifstream file("mesh.txt"); // ����, ���������� ����� ����� � �������� ������� � ��-���
   file >> n;

   f.resize(n);
   x.resize(n);
   lambda.resize(n);
   gamma.resize(n - 1);
   u.resize(n);
   b.resize(n);
   di.resize(n);
   ig.resize(n + 1);

   for (uint32_t i = 0; i < n; i++)
      file >> x[i];
   file.close();
   file.clear();

   ifstream file2("f.txt");
   for (uint32_t i = 0; i < n; i++)
      file2 >> f[i];
   file2.close();
   file2.clear();

   ifstream file3("lambda.txt");
   for (uint32_t i = 0; i < n; i++)
      file3 >> lambda[i];
   file3.close();
   file3.clear();

   ifstream file4("gamma.txt");
   for (uint32_t i = 0; i < n - 1; i++)
      file4 >> gamma[i];
   file4.close();
   file4.clear();

   ifstream file5("boundary.txt");
   file5 >> left_bound >> right_bound;
   file5.close();
   file5.clear();

   computeGlobalMatrixAndRHS(di, ig, jg, ggl, ggu, b, x, lambda, gamma, f, n);

   includeBoundary(di, b, n, left_bound, right_bound);

   solveLLt(di, ig, jg, ggl, ggu, b, u, n);

   ofstream fileout("solution.txt");
   for (uint32_t i = 0; i < n; i++)
      fileout << u[i] << endl;
   fileout.close();
   fileout.clear();

   cin >> x_arbtr;
   cout << arbitrarySolution(x, u, n, x_arbtr);

   return 0;
}
