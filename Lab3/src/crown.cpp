#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("crown");
  gmsh::option::setNumber("Mesh.MeshSizeMax", 0.1);

  double lc = 0.02;
  
  unsigned points[7];
  points[0] = gmsh::model::geo::addPoint(-3, 0, 0);
  points[1] = gmsh::model::geo::addPoint(-3, 4, 0);
  points[2] = gmsh::model::geo::addPoint(-1, 2, 0);
  points[3] = gmsh::model::geo::addPoint(0, 5, 0);
  points[4] = gmsh::model::geo::addPoint(1, 2, 0);
  points[5] = gmsh::model::geo::addPoint(3, 4, 0);
  points[6] = gmsh::model::geo::addPoint(3, 0, 0);

  for (unsigned i = 0; i < 6; ++i) {
    gmsh::model::geo::addLine(points[i], points[i + 1], i);
  }
  gmsh::model::geo::addLine(points[6], points[0], 6);

  gmsh::model::geo::addCurveLoop({0, 1, 2, 3, 4, 5, 6}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);
  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(2);
  gmsh::write("./mesh/crown.msh");
  gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}