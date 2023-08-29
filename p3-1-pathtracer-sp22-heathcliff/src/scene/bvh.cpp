#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

 
      BBox bbox, centroid;
      int count = 0;
      for (auto p = start; p != end; p++) {
        BBox bb = (*p)->get_bbox();
        bbox.expand(bb);
        centroid.expand(bb.centroid());
        count++;
      }

      auto *node = new BVHNode(bbox);
      node->start = start;
      node->end = end;
      // return if there are no more than max leaf
      // split along the biggest axis if there are more than max leaf
      if (count > max_leaf_size) {
          auto width = centroid.extent.x, height = centroid.extent.y, depth = centroid.extent.z;
          std::vector<Primitive*>* l = new std::vector<Primitive*>;
          std::vector<Primitive*>* r = new std::vector<Primitive*>;
          if (width >= height && width >= depth) {
              auto split = centroid.min.x + width * 0.5;
              for (auto p = start; p != end; p++) {
                  if ((*p)->get_bbox().centroid().x <= split) { (*l).push_back(*p); }
                  else { (*r).push_back(*p); }
              }
          }
          else if (height >= width && height >= depth) {
              auto split = centroid.min.y + height * 0.5;
              for (auto p = start; p != end; p++) {
                  if ((*p)->get_bbox().centroid().y <= split) { (*l).push_back(*p); }
                  else { (*r).push_back(*p); }
              }
          }
          else{
              auto split = centroid.min.z + depth * 0.5;
              for (auto p = start; p != end; p++) {
                  if ((*p)->get_bbox().centroid().z <= split) { (*l).push_back(*p); }
                  else { (*r).push_back(*p); }
              }
          }
          node->l = construct_bvh((*l).begin(), (*l).end(), max_leaf_size);
          node->r = construct_bvh((*r).begin(), (*r).end(), max_leaf_size);
      }
      

      return node;


}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
    double t1, t2;
    if (!node->bb.intersect(ray, t1, t2)) { return false; }
    /*if (ray.min_t > t1 || ray.max_t < t2) { return false; }*/
    if (node->isLeaf()) {
        for (auto p = node->start; p != node->end; p++) {
            total_isects++;
            if ((*p)->has_intersection(ray))
                return true;
        }
    }
    else {
        total_isects++;
        return has_intersection(ray, node->l) || has_intersection(ray, node->r);
    }
    return false;


}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
    // safetly return false if the t interval has an empty intersection with the ray's valid min and max.
    // have to return the closest intersection so must check every primitive in the box
    // if a ray starts or ends within a bbox, that counts as a valid intersection
    double t1, t2;
    if (!node->bb.intersect(ray, t1, t2)) { return false; }
    if (node->isLeaf()) {
        bool hit = false;
        for (auto p = node->start; p != node->end; p++) {
            total_isects++;
            hit = (*p)->intersect(ray, i) || hit;
        }
        return hit;
    }
    else {
        total_isects++;
        auto hit1 = intersect(ray, i, node->r);
        auto hit2 = intersect(ray, i, node->l);
        return hit1 || hit2;
    }
}
  /*bool hit = false;
  for (auto p : primitives) {
    total_isects++;
    hit = p->intersect(ray, i) || hit;
  }
  return hit;*/

} // namespace SceneObjects
} // namespace CGL
