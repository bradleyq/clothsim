#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1.1): Build a grid of masses.
  double dw = this->width / this->num_width_points;
  double dh = this->height / this->num_height_points;
  if (this->orientation == HORIZONTAL) {
    for (int h = 0; h < this->num_height_points; h++) {
      for (int w = 0; w < this->num_width_points; w++) {
        bool pin = false;
        for (vector<int>& v: this->pinned) {
          if (v[0] == h && v[1]== w) {
            pin = true;
          }
        }
        this->point_masses.push_back(PointMass(Vector3D(w*dw, 1.0f, h*dh), pin));
        //cout << Vector3D(w*dw, 1.0f, h*dh) << "\n";
      }
    }
  } else {
    for (int h = 0; h < this->num_height_points; h++) {
      for (int w = 0; w < this->num_width_points; w++) {
        bool pin = false;
        for (vector<int>& v: this->pinned) {
          if (v[0] == h && v[1]== w) {
            pin = true;
          }
        }
        this->point_masses.push_back(PointMass(Vector3D(w*dw, h*dh, rand() / double(RAND_MAX) * 0.002 - 0.001), pin));
      }
    }
  }

  // TODO (Part 1.2): Add springs 
  for (int h = 0; h < this->num_height_points; h++) {
    for (int w = 0; w < this->num_width_points; w++) {
      PointMass* curr = &this->point_masses[w + this->num_width_points * h];
      if (h > 0) {
        this->springs.push_back(Spring(curr, curr - this->num_width_points, STRUCTURAL));
        if (w < this->num_width_points - 1) {
          this->springs.push_back(Spring(curr, curr - this->num_width_points + 1, SHEARING));
        }
        if (h > 1) {
          this->springs.push_back(Spring(curr, curr - 2 * this->num_width_points, BENDING));
        }
      }
      if (w > 0) {
        this->springs.push_back(Spring(curr, curr - 1, STRUCTURAL));
        if (h > 0) {
          this->springs.push_back(Spring(curr, curr - this->num_width_points - 1, SHEARING));
        }
        if (w > 1) {
          this->springs.push_back(Spring(curr, curr - 2, BENDING));
        }
      }
    }
  }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2.1): Compute total force acting on each point mass.
  for (PointMass& p : this->point_masses) {
    p.forces = Vector3D(0.0f, 0.0f, 0.0f);
    for (Vector3D& a : external_accelerations) {
      p.forces += mass * a;
    }
  }
  for (Spring& s : this->springs) {
    if (s.spring_type == STRUCTURAL && cp->enable_structural_constraints ||
        s.spring_type == SHEARING && cp->enable_shearing_constraints ||
        s.spring_type == BENDING && cp->enable_bending_constraints) {
      double factor = s.spring_type == BENDING ? 0.2 : 1.0; 
      PointMass* m1 = s.pm_a;
      PointMass* m2 = s.pm_b;
      Vector3D dir = m2->position - m1->position;
      double F = factor * cp->ks * (abs(dir.norm()) - s.rest_length);
      //cout << F << "\n";
      dir.normalize();
      m1->forces += F * dir;
      m2->forces += -F * dir;
    }  
  }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass positions
  for (PointMass& p : this->point_masses) {
    if (!p.pinned) {
      Vector3D newpos = p.position + (1.0 - cp->damping * 0.01) * (p.position - p.last_position) + p.forces * delta_t * delta_t / mass;
      p.last_position = p.position;
      p.position = newpos;
    }
  }


  // This won't do anything until you complete Part 4.
  build_spatial_map();
  for (PointMass &pm : point_masses) {
    self_collide(pm, simulation_steps);
  }

  // This won't do anything until you complete Part 3.
  for (PointMass &pm : point_masses) {
    for (CollisionObject *co : *collision_objects) {
      co->collide(pm);
    }
  }


  // TODO (Part 2.3): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  double constrain = 0.1;
  cout << "New iter:\n";
  for (Spring& s : this-> springs) {
    if (s.spring_type == STRUCTURAL && cp->enable_structural_constraints ||
        s.spring_type == SHEARING && cp->enable_shearing_constraints ||
        s.spring_type == BENDING && cp->enable_bending_constraints) {
      PointMass* m1 = s.pm_a;
      PointMass* m2 = s.pm_b;
      Vector3D nd = m2->position - m1->position;
      //double ratio = abs(nd.norm()) / s.rest_length;
      if (abs(nd.norm()) > s.rest_length * (1.0 + constrain)) {
        nd.normalize();
        Vector3D shift = (1.0 + constrain) * s.rest_length * nd;
        if (!m1->pinned && !m2->pinned) {
          //double corr = (ratio - 1.0 - constrain) / 2.0;
          //cout << nd << m1->position << m2->position << corr;
          Vector3D mp = (m1->position + m2->position) * 0.5;
          m1->position = -0.5 * shift + mp;
          m2->position = 0.5 * shift + mp;
          //cout << m1->position << m2->position << "\n";
        } else if (!m1->pinned && m2->pinned) {
          m1->position = -shift + m2->position;
        } else if (m1->pinned && !m2->pinned) {
          m2->position = shift + m1->position;
        }
      }
      else if (abs(nd.norm()) < s.rest_length * (1.0 - constrain)) {
        nd.normalize();
        Vector3D shift = (1.0 - constrain) * s.rest_length * nd;
        if (!m1->pinned && !m2->pinned) {
          Vector3D mp = (m1->position + m2->position) * 0.5;
          m1->position = -0.5 * shift + mp;
          m2->position = 0.5 * shift + mp;
        } else if (!m1->pinned && m2->pinned) {
          m1->position = -shift + m2->position;
        } else if (m1->pinned && !m2->pinned) {
          m2->position = shift + m1->position;
        }
      }
    }  
  }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4.2): Build a spatial map out of all of the point masses.
  for (PointMass& p : this->point_masses) {
    float hval = hash_position(p.position);
    if (this->map.find(hval) == this->map.end()) {
      vector<PointMass *>* vec = new vector<PointMass *>();
      this->map[hval] = vec;
    }
    vector<PointMass *>* v = this->map[hval];
    (*v).push_back(&p);
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4.3): Handle self-collision for a given point mass.
  vector<PointMass *>* v = this->map[hash_position(pm.position)];
  //cout << (*v).size() << "\n";
  Vector3D correction = Vector3D();
  int num = 0;
  for (PointMass* pp : *v) {
    Vector3D dir = pm.position - pp->position;
    if (pp != &pm && dir.norm2() < 4.0 * this->thickness * this->thickness) {
      dir.normalize();
      num += 1;
      Vector3D tmp = pp->position + dir * 2.0 * this->thickness;
      correction += tmp - pm.position;
    }
  }
  if (num > 0) {
    pm.position = pm.position + correction / (num * simulation_steps);
  }
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.
  float dx,dy,dz;
  if (this->orientation == HORIZONTAL) {
    dx = 3 * this->width / this->num_width_points;
    dz = 3 * this->height / this->num_height_points;
    dy = max(dx, dz);
  }
  else {
    dx = 3 * this->width / this->num_width_points;
    dy = 3 * this->height / this->num_height_points;
    dz = max(dx, dy);
  }
  return (floor(pos.x / dx) * 31 + floor(pos.y / dy)) * 31 + floor(pos.z / dz);
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm, pm + num_width_points, pm + 1));
      triangles.push_back(new Triangle(pm + 1, pm + num_width_points,
                                       pm + num_width_points + 1));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
