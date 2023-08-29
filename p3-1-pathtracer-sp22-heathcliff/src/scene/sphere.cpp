#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
	auto a = dot(r.d, r.d), b = 2 * dot((r.o - o), r.d), c = dot((r.o - o), (r.o - o)) - r2;
	auto D = b * b - 4 * a * c;
	if (D < 0) {return false;}
	else { 
		auto s_max = (-b + sqrt(D)) / (2 * a), s_min = (-b - sqrt(D)) / (2 * a); 
		if (s_min >= r.min_t && s_max <= r.max_t) {
			t1 = s_min;
			t2 = s_max;
			r.max_t = s_min;
			return true;
		}
		/*else if (r.max_t >= s_min >= r.min_t) {
			t1 = s_min;
			r.max_t = s_min;
			return true;
		}
		else if (r.max_t >= s_min >= r.min_t) {
			t1 = t2 = s_max;
			r.max_t = s_max;
			return true;
		}*/
	}
	return false;
}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
	double t1, t2;
	return test(r, t1, t2);

}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
	double t1, t2;
	auto flag =  test(r, t1, t2);
	if (!flag) { return false; }
	i->t = r.max_t;
	i->primitive = this;
	i->bsdf = get_bsdf();
	Vector3D normal = r.d * r.max_t +r.o;
	i->n = normal - this->o;
	i->n.normalize();
	/*i->n = normal;*/
	return flag;


}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
