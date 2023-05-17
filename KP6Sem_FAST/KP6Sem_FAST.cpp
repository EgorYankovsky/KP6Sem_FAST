#include"Header.h"
using namespace solution;

int main()
{
   setlocale(LC_ALL, "Russian");
   Difur object;

   object.TimeGenerator(1.0, 5, 0.0, 1.0);
   object.Grid_Generator();
   object.Read();
   object.Portrait(object.globalMatrixG);
   object.localG(object.globalMatrixG);
   object.Portrait(object.globalMatrixM); // список смежности элементов (портрет)?
   object.localM(object.globalMatrixM);
   object.Portrait(object.globalMatrix);
   
   //object.GaMbo();

   object.InitiateConditions(0);
   object.InitiateConditions(1);
   for (size_t i(2); i < object.t.size(); i++)
   {
      object.ImplicitScheme(i);
      object.Boundary_conditions(object.t[i]);
      object.LU();
      object.LOS();
      object.q_h[i] = object.P;
   }
   object.output();
   return 0;
}