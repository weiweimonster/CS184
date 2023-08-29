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
    unsigned int pos_buffer;
    glGenBuffers(1, &pos_buffer);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, pos_buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, num_width_points * num_height_points * sizeof(PointMass), NULL, GL_DYNAMIC_DRAW);
    unsigned int bufMask = GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT;
    PointMass *points = (PointMass *)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, num_width_points * num_height_points * sizeof(PointMass), bufMask);
    int count = 0;
    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            if (this->orientation == HORIZONTAL) {
                // xz plane and set y to 1 for all point masses
                // top left (first) points mass be at (0, 0)
                // bottom right (last) point mass be at coordinate (width, hegiht)
                auto position = Vector3D(i * width / (num_width_points - 1.0), 1.0, j * height / (num_height_points - 1.0));
                points[count] = PointMass(position, false);
                /*auto position = Vector3D(i * width / (num_width_points - 1.0), 1.0, j * height / (num_height_points - 1.0));
                point_masses.emplace_back(PointMass(position, false));*/
            }
            else if (this->orientation == VERTICAL) {
                // xy plane and set z to a random number between -0.001 and 0.001 (rand())
                auto offset = -0.001 + rand() / (RAND_MAX / 0.002);
                auto position = Vector3D(i * width / (num_width_points - 1.0), j * height / (num_height_points - 1.0), offset);
                points[count] = PointMass(position, false);
                /*auto position = Vector3D(i * width / (num_width_points - 1.0), j * height / (num_height_points - 1.0), offset);
                point_masses.emplace_back(PointMass(position, false));*/
            }
            count++; 
        }
    }
    memcpy(&point_masses, &points, num_width_points * num_height_points * sizeof(PointMass));
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER); 


  // TODO (Part 1): Build a grid of masses and springs.
    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            if (this->orientation == HORIZONTAL) {
                // xz plane and set y to 1 for all point masses
                // top left (first) points mass be at (0, 0)
                // bottom right (last) point mass be at coordinate (width, hegiht)
                auto position = Vector3D(i * width / (num_width_points - 1.0), 1.0, j * height / (num_height_points - 1.0));
                point_masses.emplace_back(PointMass(position, false));
            }
            else if (this -> orientation == VERTICAL) {
                // xy plane and set z to a random number between -0.001 and 0.001 (rand())
                auto offset = -0.001 + rand() / (RAND_MAX / 0.002);
                auto position = Vector3D(i * width / (num_width_points - 1.0), j * height / (num_height_points - 1.0), offset);
                point_masses.emplace_back(PointMass(position, false));
            }
        }
    }
    for (auto idx : this->pinned) {
        point_masses[idx[1] * num_width_points + idx[0]].pinned = true;
    }
    // populate springs
    // structural constraints between self,left and self,top
    // shearing constraints between self,upperleft and self,upperright
    // bending constraints between self, two left and self, two top
    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            //structural constraints
            if (i >= 1) {
                auto *point_a = &point_masses[j * num_width_points + i], *point_b = &point_masses[j * num_width_points + i - 1];
                auto s = Spring(point_a, point_b,STRUCTURAL);
                this->springs.emplace_back(s);
            }
            if (j >= 1) {
                auto *point_a = &point_masses[j * num_width_points + i], *point_b = &point_masses[(j - 1) * num_width_points + i];
                auto s = Spring(point_a, point_b, STRUCTURAL);
                this->springs.emplace_back(s);
            }
            //shearing constraints
            if (j >= 1) {
                if (i >= 1) {
                    auto *point_a = &point_masses[j * num_width_points + i], *point_b = &point_masses[(j - 1) * num_width_points + i - 1];
                    auto s = Spring(point_a, point_b, SHEARING);
                    this->springs.emplace_back(s);
                }
                if (i < num_width_points - 1) {
                    auto *point_a = &point_masses[j * num_width_points + i], *point_b = &point_masses[(j - 1) * num_width_points + i + 1];
                    auto s = Spring(point_a, point_b, SHEARING);
                    this->springs.emplace_back(s);
                }
            }
            //bending constraints
            if (i >= 2) {
                auto *point_a = &point_masses[j * num_width_points + i], *point_b = &point_masses[j * num_width_points + i - 2];
                auto s = Spring(point_a, point_b, BENDING);
                this->springs.emplace_back(s);
            }
            if (j >= 2) {
                auto *point_a = &point_masses[j * num_width_points + i], *point_b = &point_masses[(j - 2) * num_width_points + i];
                auto s = Spring(point_a, point_b, BENDING);
                this->springs.emplace_back(s);
            }
        }
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
  auto total_a = Vector3D();
  for (Vector3D &a : external_accelerations) {
      total_a += a;
  }
  for (PointMass &point : this->point_masses) {
      point.forces = total_a * mass;
  }
  for (Spring &spring : this->springs) {
      if ((cp->enable_bending_constraints && spring.spring_type == BENDING) || (cp->enable_shearing_constraints && spring.spring_type == SHEARING) || (cp->enable_structural_constraints && spring.spring_type == STRUCTURAL)) {
          double ks;
          if (spring.spring_type == BENDING) { ks = cp->ks * 0.2; }
          else { ks = cp->ks; }
          auto direction = (spring.pm_b->position - spring.pm_a->position).unit();
          auto f_s = direction * ks * ((spring.pm_b->position - spring.pm_a->position).norm() - spring.rest_length);
          spring.pm_a->forces += f_s;
          spring.pm_b->forces -= f_s;
      }
  }
  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  auto d = cp->damping / 100.0;
  for (PointMass &point : this->point_masses) {
      if (!point.pinned) {
          auto temp = point.position + (1.0 - d) * (point.position - point.last_position) + point.forces / mass * delta_t * delta_t;
          point.last_position = point.position;
          point.position = temp;
      }
      
  }

  // TODO (Part 4): Handle self-collisions.
  build_spatial_map();
  // TODO (Part 3): Handle collisions with other primitives.
  for (auto& point : this->point_masses) {
      this->self_collide(point, simulation_steps);
      for (auto o : *collision_objects) {
          o->collide(point);
      }
  }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for (Spring &spring : this->springs) {
      auto elong = (spring.pm_a->position - spring.pm_b->position).norm() - spring.rest_length;
      auto dir = (spring.pm_a->position - spring.pm_b->position).unit();
      if (elong > spring.rest_length * 0.1) {
          if (!spring.pm_a->pinned && !spring.pm_b->pinned) {
              spring.pm_b->position += (elong - spring.rest_length * 0.1) * 0.5 * dir;
              spring.pm_a->position -= (elong - spring.rest_length * 0.1) * 0.5 * dir;
          }
          else if (!spring.pm_a->pinned && spring.pm_b->pinned){ spring.pm_a->position -= (elong - spring.rest_length * 0.1) * dir; }
          else if (spring.pm_a->pinned && !spring.pm_b->pinned){ spring.pm_b->position += (elong - spring.rest_length * 0.1) * dir; }
      }
  }
}

void Cloth::build_spatial_map() {
    for (const auto &entry : map) {
        delete(entry.second);
    }
    map.clear();
  // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (auto& point : this->point_masses) {
        auto key = hash_position(point.position);
        if (this->map[key]) { this->map[key]->push_back(&point); }
        else { this->map[key] = new vector <PointMass *>; }
    }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.
    auto key = hash_position(pm.position);
    auto *value = this->map[key];
    auto total_correction = Vector3D();
    auto count = 0;
    for (auto *point : *value) {
        auto dist = (point->position - pm.position).norm();
        if (point == &pm) { continue; }
        if (dist < 2 * this->thickness) {
            count++;
            auto correction = (2 * this->thickness - dist) * (pm.position - point->position).unit();
        }
    }
    if (count > 0) { pm.position += (total_correction / float(count) / simulation_steps); }
    
}

float Cloth::hash_position(Vector3D pos) {
    // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    auto w = 3 * width / num_width_points;
    auto h = 3 * height / num_height_points;
    auto t = max(w, h);
    auto x_index = floor(pos.x / w), y_index = floor(pos.y / h), z_index = floor(pos.z / t);
    //return 19 * x_index * x_index + 7 * y_index + 3 * z_index;
    return 0;
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
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
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
