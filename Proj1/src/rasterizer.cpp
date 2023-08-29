#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    double rate = sqrt(sample_rate);
    int num = 0;
    for (int z = 0; z < rate; z++) {
        for (int w = 0; w < rate; w++) {
            sample_buffer[sample_rate * (y * width + x) + num] = c;
            num++;
        }
    }
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }
  float line_test(float x, float y, float x1, float y1, float x2, float y2) {
      return -(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1); 
  }
  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
      float x_max =min(max(max(x0, x1), x2) , float(width)-1);
      float x_min = max(min(min(x0, x1), x2) , float(0.0));
      float y_max = min(max(max(y0, y1), y2) , float(height)-1);
      float y_min = max(min(min(y0, y1), y2) , float(0.0));
      float rate = sqrt(sample_rate);
      float step_size = 1 / rate;
      float offset = 1 / (2 * rate);
      for (int i = x_min; i <= x_max; i++) {
          for (int j = y_min; j <= y_max; j++) {
              /*int line1 = line_test(floor(i) + 0.5, floor(j) + 0.5, x0, y0, x1, y1);
              int line2 = line_test(floor(i) + 0.5, floor(j) + 0.5, x1, y1, x2, y2);
              int line3 = line_test(floor(i) + 0.5, floor(j) + 0.5, x2, y2, x0, y0);
              if ((line1 >= 0 && line2 >= 0 && line3 >= 0) || (line1 <= 0 && line2 <= 0 && line3 < 0)) {
                  rasterize_point(i, j, color);
              }*/
              int num = 0;
              for (int z = 0; z < rate; z++) {
                  for (int w = 0; w < rate; w++) {
                      float line1 = line_test(float(i) + float(z)*step_size + offset, float(j) + float(w)*step_size + offset, x0, y0, x1, y1);
                      float line2 = line_test(float(i) + float(z)*step_size + offset, float(j) + float(w)*step_size + offset, x1, y1, x2, y2);
                      float line3 = line_test(float(i) + float(z)*step_size + offset, float(j) + float(w)*step_size + offset, x2, y2, x0, y0);
                      if ((line1 >= 0.0 && line2 >= 0.0 && line3 >= 0.0) || (line1 <= 0.0 && line2 <= 0.0 && line3 <= 0.0)) {
                          sample_buffer[sample_rate * (j* width +i) + num] = color;
                      }
                      num++;
                  }
                  
              }

          }
      }
    // TODO: Task 2: Update to implement super-sampled rasterization



  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
      float x1, float y1, Color c1,
      float x2, float y2, Color c2)
  {
      // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
      // Hint: You can reuse code from rasterize_triangle
      float x_max = min(max(max(x0, x1), x2), float(width) - 1);
      float x_min = max(min(min(x0, x1), x2), float(0.0));
      float y_max = min(max(max(y0, y1), y2), float(height) - 1);
      float y_min = max(min(min(y0, y1), y2), float(0.0));
      float rate = sqrt(sample_rate);
      float step_size = 1 / rate;
      float offset = 1 / (2 * rate);
      for (int i = x_min; i <= x_max; i++) {
          for (int j = y_min; j <= y_max; j++) {
              int num = 0;
              for (int z = 0; z < rate; z++) {
                  for (int w = 0; w < rate; w++) {
                      float line1 = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset, x0, y0, x1, y1);
                      float line2 = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset, x1, y1, x2, y2);
                      float line3 = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset, x2, y2, x0, y0);
                      if ((line1 >= 0.0 && line2 >= 0.0 && line3 >= 0.0) || (line1 <= 0.0 && line2 <= 0.0 && line3 <= 0.0)) {
                          float alpha, beta, gamma;
                          alpha = line1 / line_test(x2, y2, x0, y0, x1, y1);
                          beta = line2 / line_test(x0, y0, x1, y1, x2, y2);
                          gamma = 1.0 - alpha - beta;
                          Color color = alpha * c2 + beta * c0 + gamma * c1;
                          sample_buffer[sample_rate * (j * width + i) + num] = color;
                      }
                      num++;
                  }

              }

          }


      }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
      SampleParams sp;
      sp.psm = psm;
      sp.lsm = lsm;
      float x_max = min(max(max(x0, x1), x2), float(width) - 1);
      float x_min = max(min(min(x0, x1), x2), float(0.0));
      float y_max = min(max(max(y0, y1), y2), float(height) - 1);
      float y_min = max(min(min(y0, y1), y2), float(0.0));
      float rate = sqrt(sample_rate);
      float step_size = 1 / rate;
      float offset = 1 / (2 * rate);
      for (int i = x_min; i <= x_max; i++) {
          for (int j = y_min; j <= y_max; j++) {
              int num = 0;
              for (int z = 0; z < rate; z++) {
                  for (int w = 0; w < rate; w++) {
                      float line1 = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset, x0, y0, x1, y1);
                      float line2 = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset, x1, y1, x2, y2);
                      float line3 = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset, x2, y2, x0, y0);
                      if ((line1 >= 0.0 && line2 >= 0.0 && line3 >= 0.0) || (line1 <= 0.0 && line2 <= 0.0 && line3 <= 0.0)) {
                          float alpha, beta, gamma, alpha_dx, beta_dx, gamma_dx, alpha_dy, beta_dy, gamma_dy, alpha_line1, beta_line2;
                          alpha_line1 = line_test(x2, y2, x0, y0, x1, y1);
                          beta_line2 = line_test(x0, y0, x1, y1, x2, y2);
                          alpha = line1 / alpha_line1;
                          beta = line2 / beta_line2;
                          gamma = 1.0 - alpha - beta;
                          float line1_dx = line_test(float(i) + float(z) * step_size + offset + 1.0, float(j) + float(w) * step_size + offset, x0, y0, x1, y1);
                          float line2_dx = line_test(float(i) + float(z) * step_size + offset + 1.0, float(j) + float(w) * step_size + offset, x1, y1, x2, y2);
                          alpha_dx = line1_dx / alpha_line1;
                          beta_dx = line2_dx / beta_line2;
                          gamma_dx = 1.0 - alpha_dx - beta_dx;
                          float line1_dy = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset + 1.0, x0, y0, x1, y1);
                          float line2_dy = line_test(float(i) + float(z) * step_size + offset, float(j) + float(w) * step_size + offset + 1.0, x1, y1, x2, y2);
                          alpha_dy = line1_dy / alpha_line1;
                          beta_dy = line2_dy / beta_line2;
                          gamma_dy = 1.0 - alpha_dy - beta_dy;
                          sp.p_uv = Vector2D(alpha * u2 + beta * u0 + gamma * u1, alpha * v2 + beta * v0 + gamma * v1);
                          sp.p_dx_uv = Vector2D(alpha_dx * u2 + beta_dx * u0 + gamma_dx * u1, alpha_dx * v2 + beta_dx * v0 + gamma_dx * v1);
                          sp.p_dy_uv = Vector2D(alpha_dy * u2 + beta_dy * u0 + gamma_dy * u1, alpha_dy * v2 + beta_dy * v0 + gamma_dy * v1);
                          sample_buffer[sample_rate * (j * width + i) + num] = tex.sample(sp);
                      }
                      num++;
                  }

              }

          }


      }



  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
    float rate = sqrt(sample_rate);
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col;
        int num = 0;
        for (int z = 0; z < rate; z++) {
            for (int w = 0; w < rate; w++) {
                col += (sample_buffer[sample_rate * (y * width + x) + num]) * (1 / double(sample_rate));
                num++;
            }
        }
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }
    
  }

  Rasterizer::~Rasterizer() { }


}// CGL
