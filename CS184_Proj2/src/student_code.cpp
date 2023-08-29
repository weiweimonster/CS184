#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{   
    
  template <class T>
  inline T lerp(const T &b1,const T &b2, float t) {
      return (1 - t) * b1 + t * b2;
  }
  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
      vector<Vector2D> next_level;
      int num = points.size();
      for (int i = 0; i < num-1; i++) {
          Vector2D new_b = lerp<Vector2D>(points[i], points[i + 1], t);
          next_level.push_back(new_b);
      }
    return next_level;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      vector<Vector3D> next_level;
      int num = points.size();
      for (int i = 0; i < num - 1; i++) {
          Vector3D new_b = lerp<Vector3D>(points[i], points[i + 1], t);
          next_level.push_back(new_b);
      }
      return next_level;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      vector<Vector3D>next_level = points;
      next_level = evaluateStep(next_level, t);
      while (next_level.size() != 1) {
          next_level = evaluateStep(next_level, t);
      }
      return next_level[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
        // TODO Part 2.
       vector<Vector3D> moving_curve;
       for (int i = 0; i < controlPoints.size(); i++) {
           Vector3D new_b = evaluate1D(controlPoints[i], u);
           moving_curve.push_back(new_b);
       }
       return evaluate1D(moving_curve, v);
   }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      HalfedgeCIter h = this->halfedge()->twin();
      Vector3D normal;
      do {
          normal += h->face()->normal();
          h = h->next()->twin();
      } while (h != this->halfedge()->twin() && !h->isBoundary());
      normal.normalize();
      return normal;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
      if (!e0->isBoundary()) {
          HalfedgeIter h0 = e0->halfedge();
          HalfedgeIter h1 = h0->next();
          HalfedgeIter h2 = h1->next();
          HalfedgeIter h3 = h0->twin();
          HalfedgeIter h4 = h3->next();
          HalfedgeIter h5 = h4->next();
          FaceIter f0 = h0->face();
          FaceIter f1 = h3->face();
          VertexIter v0 = h0->vertex();
          VertexIter v1 = h1->vertex();
          VertexIter v2 = h2->vertex();
          VertexIter v3 = h5->vertex();


          h0->setNeighbors(h5, h3, v2, e0, f0);
          h1->setNeighbors(h0, h1->twin(), v1, h1->edge(), f0);
          h2->setNeighbors(h4, h2->twin(), v2, h2->edge(), f1);
          h3->setNeighbors(h2, h0, v3, e0, f1);
          h4->setNeighbors(h3, h4->twin(), v0, h4->edge(), f1);
          h5->setNeighbors(h1, h5->twin(), v3, h5->edge(), f0);

          v0->halfedge() = h4;
          v1->halfedge() = h1;
          v2->halfedge() = h0;
          v3->halfedge() = h3;

          f0->halfedge() = h0;
          f1->halfedge() = h4;

          
      }
      return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
      if (!e0->isBoundary()) {
          HalfedgeIter h0 = e0->halfedge();
          HalfedgeIter h1 = h0->next();
          HalfedgeIter h2 = h1->next();
          HalfedgeIter h3 = h0->twin();
          HalfedgeIter h4 = h3->next();
          HalfedgeIter h5 = h4->next();
          FaceIter f0 = h0->face();
          FaceIter f1 = h3->face();
          VertexIter v0 = h0->vertex();
          VertexIter v1 = h1->vertex();
          VertexIter v2 = h2->vertex();
          VertexIter v3 = h5->vertex();

          HalfedgeIter h6 = newHalfedge();
          HalfedgeIter h7 = newHalfedge();
          HalfedgeIter h8 = newHalfedge();
          HalfedgeIter h9 = newHalfedge();
          HalfedgeIter h10 = newHalfedge();
          HalfedgeIter h11 = newHalfedge();

          VertexIter midpoint = newVertex();
          midpoint->position = (v0->position + v1->position)/2.0;
          midpoint->halfedge() = h0;
         

          FaceIter f2 = newFace();
          FaceIter f3 = newFace();

          EdgeIter e1 = newEdge();
          EdgeIter e2 = newEdge();
          EdgeIter e3 = newEdge();
          

          h0->setNeighbors(h1, h3, midpoint, e0, f0);
          h1->setNeighbors(h6, h1->twin(), v1, h1->edge(), f0);
          h2->setNeighbors(h8, h2->twin(), v2, h2->edge(), f1);
          h3->setNeighbors(h11, h0, v1, e0, f3);
          h4->setNeighbors(h10, h4->twin(), v0, h4->edge(), f2);
          h5->setNeighbors(h3, h5->twin(), v3, h5->edge(), f3);
          h6->setNeighbors(h0, h7, v2, e1, f0);
          h7->setNeighbors(h2, h6, midpoint, e1, f1);
          h8->setNeighbors(h7, h9, v0, e2, f1);
          h9->setNeighbors(h4, h8, midpoint, e3, f2);
          h10->setNeighbors(h9, h11, v3, e3, f2);
          h11->setNeighbors(h5, h10, midpoint, e3, f3);

          v0->halfedge() = h8;
          v1->halfedge() = h1;
          v2->halfedge() = h2;
          v3->halfedge() = h5;
          

          f0->halfedge() = h0;
          f1->halfedge() = h2;
          f2->halfedge() = h4;
          f3->halfedge() = h3;

          e1->halfedge() = h6;
          e2->halfedge() = h8;
          e3->halfedge() = h10;
          midpoint->isNew = true;
          e0->isNew = false;
          e1->isNew = true;
          e2->isNew = false;
          e3->isNew = true;
          

      }
      return e0->halfedge()->vertex();
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          HalfedgeIter h;
          VertexIter v0, v1, v2, v3;
          h = e->halfedge();
          v0 = h->vertex();
          v1 = h->next()->vertex();
          v2 = h->next()->next()->vertex();
          v3 = h->twin()->next()->next()->vertex();
          e->newPosition = (v0->position + v1->position) * 0.375 + (v2->position + v3->position) * 0.125;
      }
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->isNew = false;
          HalfedgeIter h = v->halfedge();
          int degree = 0;
          Vector3D position_sum;
          do {          
              degree += 1;
              HalfedgeIter twin = h->twin();
              VertexIter v_next = twin->vertex();
              position_sum += v_next->position;
              h = twin->next();
          } while (h != v->halfedge());
          double u;
          if (degree == 3) {
              u == 3.0 / 16.0;
          }
          else {
              u = 3.0 / (8.0 * double(degree));
          }
          v->newPosition = (1.0 - double(degree) * u)*v->position + u * position_sum; 
      }
      EdgeIter e = mesh.edgesBegin();
      int edge_count = 0;
      int edge_num = mesh.nEdges();
      while (e != mesh.edgesEnd() && edge_count < edge_num) {
          edge_count++;
          EdgeIter next_e = e;
          next_e++;
          VertexIter v = mesh.splitEdge(e);
          v->newPosition = e->newPosition;
          e = next_e;
      }
      //for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      //    VertexIter v;
      //    v = mesh.splitEdge(e);
      //    v->newPosition = e->newPosition;
      //}
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          if (e->isNew) { 
              HalfedgeIter  h = e->halfedge();
              if (h->vertex()->isNew != h->twin()->vertex()->isNew) {
                  mesh.flipEdge(e);
              }
          }
      }
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->position = v->newPosition;
      }
  }
}
