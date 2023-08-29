#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.

	//Moller Trumbore Algorithm
	auto e1 = p2 - p1, e2 = p3 - p1, s = r.o - p1;
	auto s1 = cross(r.d, e2), s2 = cross(s, e1);
	auto result = (1 / dot(s1, e1)) * Vector3D(dot(s2, e2), dot(s1, s), dot(s2, r.d));
	auto t = result[0], b1 = result[1], b2 = result[2];
	if (r.min_t > t || t > r.max_t || b1 < 0 || b2 < 0 || b1 + b2>1|| t< 0 ) {
		return false;
	}
	r.max_t = t;
	return true;
   
}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
	auto e1 = p2 - p1, e2 = p3 - p1, s = r.o - p1;
	auto s1 = cross(r.d, e2), s2 = cross(s, e1);
	auto result = (1 / dot(s1, e1)) * Vector3D(dot(s2, e2), dot(s1, s), dot(s2, r.d));
	auto t = result[0], b1 = result[1], b2 = result[2];
	if (t < r.min_t|| b1 < 0 || b2 < 0 || b1 + b2 > 1 || t < 0 || t > r.max_t) {
		return false;
	}
	
	r.max_t = t;
	isect->bsdf = get_bsdf();
	isect->t =t;
	isect->primitive = this;
	isect->n = b1 * n2 + b2 * n3 + (1.0 - b1 - b2) * n1;
	
	return true;

}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
