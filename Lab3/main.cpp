#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0];
    double dy = x[1] - 1;
    values[0] = 500 * exp(-(dy*dy) * 1000);
    values[1] = 2000 * exp(-(dy*dy) * 1000);
  }
};


class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double y = x[1] - 1;
    values[0] = 20 * exp(-1*(y*y))*sin(y*y)*sin(y*y);
  }
};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    // return on_boundary;
    return (x[0]*x[0] + (x[1]-1)*(x[1]-1) < 0.25) || ((x[0]+2)*(x[0]+2) + (x[1]-1)*(x[1]-1) < 0.25) || ((x[0]-2)*(x[0]-2) + (x[1]-1)*(x[1]-1) < 0.25);
  }
};


int main()
{
  auto mesh = std::make_shared<Mesh>("./mesh/crown.xml");
  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  Function u(V);
  solve(a == L, u, bc);

  File file("res/poisson.pvd");
  file << u;

  return 0;
}
