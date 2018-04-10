#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3.1): Handle collisions with spheres.
  Vector3D dir = pm.position - this->origin;
  if (dir.norm2() < this->radius2) {
    dir.normalize();
    Vector3D correction = dir * this->radius + this->origin;
    correction -= pm.last_position;
    pm.position = (1.0 - this->friction) * correction + pm.last_position;
  }
  
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  Misc::draw_sphere(shader, origin, radius * 0.92);
}
