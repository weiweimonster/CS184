#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {
bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
	auto txmin = (min.x - r.o.x) / r.d.x, txmax = (max.x - r.o.x) / r.d.x;
	auto tymin = (min.y - r.o.y) / r.d.y, tymax = (max.y - r.o.y) / r.d.y;
	auto tzmin = (min.z - r.o.z) / r.d.z, tzmax = (max.z - r.o.z) / r.d.z;
	auto temp = txmin;
	if (txmin > txmax) { txmin = txmax, txmax = temp;}
	if (tymin > tymax) { temp = tymin, tymin = tymax, tymax = temp; }
	if (tzmin > tzmax) { temp = tzmin, tzmin = tzmax, tzmax = temp; }
	if (txmin > tymax || txmin > tzmax || tymin > tzmax || tymin > txmax || tzmin > txmax || tzmin > tymax) {
		return false;
	}
	t0 = std::max(txmin, std::max(tymin, tzmin));
	t1 = std::min(txmax, std::min(tymax, tzmax));
	return true;

}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
