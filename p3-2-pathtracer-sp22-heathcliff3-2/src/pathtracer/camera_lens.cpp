#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Project 3-2: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
	auto top_right = Vector3D(tan(radians(hFov) * 0.5), tan(radians(vFov) * 0.5), -1), bottom_left = Vector3D(-tan(radians(hFov) * 0.5), -tan(radians(vFov) * 0.5), -1);
	auto curr_pixel = Vector3D(x * top_right.x + (1.0 - x) * bottom_left.x, y * top_right.y + (1.0 - y) * bottom_left.y, -1);
	auto focus = curr_pixel * this->focalDistance;
	auto sample_lens = Vector3D(this->lensRadius * sqrt(rndR) * cos(rndTheta), this->lensRadius * sqrt(rndR) * sin(rndTheta), 0);
	auto r = this->c2w * (focus - sample_lens);
	r.normalize();
	auto ray = Ray(this->c2w * sample_lens + this->pos, r);
	ray.max_t = fClip, ray.min_t = nClip;
	return ray;
}


} // namespace CGL
