#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray& r,
    const Intersection& isect) {
    // Estimate the lighting from this intersection coming directly from a light.
    // For this function, sample uniformly in a hemisphere.

    // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
    // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

    // make a coordinate system for a hit point
    // with N aligned with the Z direction.
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();

    // w_out points towards the source of the ray (e.g.,
    // toward the camera if this is a primary ray)
    const Vector3D hit_p = r.o + r.d * isect.t;
    const Vector3D w_out = w2o * (-r.d);

    // This is the same number of total samples as
    // estimate_direct_lighting_importance (outside of delta lights). We keep the
    // same number of samples for clarity of comparison.
    int num_samples = scene->lights.size() * ns_area_light;
    Vector3D L_out;

    // TODO (Part 3): Write your sampling loop here
    // TODO BEFORE YOU BEGIN
    // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 
    for (int i = 0; i < num_samples; i++) {
        auto wj = UniformHemisphereSampler3D().get_sample();
        auto d = o2w * wj;
        Ray ray = Ray(hit_p, d);
        ray.min_t = EPS_F;
        Intersection isect_light;
        bool hit = this->bvh->intersect(ray, &isect_light);
        if (hit) {
            wj.normalize();
            L_out += isect.bsdf->f(w_out, wj) * isect_light.bsdf->get_emission() * cos_theta(wj);
        }
    }
    L_out = L_out * 2.0 * PI / double(num_samples);
    return L_out;
}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;
  
  
  for (auto *light : scene->lights) {
      int num_sample;
      if (light->is_delta_light()) { num_sample = 1; }
      else { num_sample = ns_area_light; }
      Vector3D L;
      double dist_to_light, pdf;
      Vector3D wi;
      for (int i = 0; i < num_sample; i++) {             
          auto radiance = (*light).sample_L(hit_p, &wi, &dist_to_light, &pdf);
          auto d = w2o * wi;
          if (d.z >= 0) {
              Ray ray = Ray(hit_p, wi);
              ray.min_t = EPS_F;
              ray.max_t = dist_to_light-EPS_F;
              Intersection isect_light;
              bool  hit = this->bvh->intersect(ray, &isect_light);
              if (!hit) {
                 /* d.normalize();*/
                  L += (isect.bsdf->f(w_out, d) * radiance * cos_theta(d)) / pdf;
              }
          }
      }
      L /= num_sample;
      L_out += L;
      
  }
  return L_out;

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light


  return isect.bsdf->get_emission();


}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

    if (direct_hemisphere_sample) { return estimate_direct_lighting_hemisphere(r, isect); }
    return estimate_direct_lighting_importance(r, isect);
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
    Matrix3x3 o2w;
    make_coord_space(o2w, isect.n);
    Matrix3x3 w2o = o2w.T();

    Vector3D hit_p = r.o + r.d * isect.t;
    Vector3D w_out = w2o * (-r.d);

    Vector3D L_out = one_bounce_radiance(r, isect);
    double CPDF = 0.35;
    // TODO: Part 4, Task 2
    // Returns the one bounce radiance + radiance from extra bounces at this point.
    // Should be called recursively to simulate extra bounces.
    Vector3D wi;
    double pdf;
    auto fraction = isect.bsdf->sample_f(w_out, &wi, &pdf);
    bool stop = coin_flip(CPDF);
    bool on = max_ray_depth > 1;
    bool max_bounce = r.depth <= max_ray_depth;
    bool max_recursion = r.depth > 1;
    if ((!stop) && on&& max_bounce&& max_recursion) {
        Ray out_ray = Ray(hit_p, o2w * wi);
        out_ray.min_t = EPS_F;
        out_ray.depth = r.depth - 1;
        Intersection isect1;
        bool hit = this->bvh->intersect(out_ray, &isect1);
        if (hit) {
            auto radiance = at_least_one_bounce_radiance(out_ray, isect1);
            if (r.depth == this->max_ray_depth)
                L_out += (wi.z * radiance * fraction) / pdf;
            else
                L_out += isect.bsdf->sample_f(w_out, &wi, &pdf)* radiance* cos_theta(wi) / pdf / (1.0 - CPDF);
        }

    }
    return L_out;

}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;


  L_out = normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.
  /*return at_least_one_bounce_radiance(r, isect);*/
  return zero_bounce_radiance(r, isect) +  at_least_one_bounce_radiance(r, isect);
  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct

  /*return L_out;*/
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

    int num_samples = ns_aa;          // total samples to evaluate
    Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
    auto w = sampleBuffer.w;
    auto h = sampleBuffer.h;
    auto total = 0.0, total_square = 0.0;
    Vector3D result;
    if (num_samples == 1) {
        Ray ray = this->camera->generate_ray((x + 0.5) / w, (y + 0.5) / h);
        ray.depth = max_ray_depth;
        result = this->est_radiance_global_illumination(ray);
        sampleBuffer.update_pixel(result, x, y);
        sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;
    }
    else {
        auto actual_sample = 0;
        for (int i = 0.0; i < num_samples; i++, actual_sample++) {
            if (i % samplesPerBatch == 0) {
                auto mean = total / double(i);
                auto std = (1 / double(i - 1)) * (total_square - total * total / double(i));
                auto illum = 1.96 * sqrt(std / double(i));
                if (illum <= maxTolerance * mean) { break; }
            }
            Vector2D sample = this->gridSampler->get_sample();
            Ray ray = this->camera->generate_ray((x + sample.x) / w, (y + sample.y) / h);
            ray.depth = max_ray_depth;
            auto radiance = this->est_radiance_global_illumination(ray);
            result += radiance;
            total += radiance.illum();
            total_square += (radiance.illum() * radiance.illum());
        }
        sampleBuffer.update_pixel(result/double(actual_sample), x, y);
        sampleCountBuffer[x + y * sampleBuffer.w] =actual_sample;
    }
  /*Vector3D result;
  for (int i = 0; i < num_samples; i++) {
      Vector2D sample = this->gridSampler->get_sample();
      Ray ray = this->camera->generate_ray((x+sample.x)/w , (y+ sample.y)/h);
      ray.depth = max_ray_depth;
      result += this->est_radiance_global_illumination(ray);
  }
  result /= num_samples;
  sampleBuffer.update_pixel(result, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = num_samples;*/

}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
