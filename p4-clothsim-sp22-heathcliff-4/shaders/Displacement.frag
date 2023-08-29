#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
 return texture(u_texture_2, uv).r;
}

void main() {
  // YOUR CODE HERE
  mat3 TBN = mat3(normalize(v_tangent.xyz), normalize(cross(v_normal.xyz, v_tangent.xyz)), normalize(v_normal.xyz));
  float dU = (h(vec2(v_uv.x + 1.0/u_texture_2_size.x, v_uv.y)) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  float dV = (h(vec2(v_uv.x, v_uv.y + 1.0/u_texture_2_size.y)) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  vec3 n_0 = vec3(-dU, -dV, 1) ;
  vec3 n_d = TBN * n_0;


  vec3 l = u_light_pos - v_position.xyz;
  vec3 l_normalized = normalize(l);
  vec3 v = u_cam_pos - v_position.xyz;
  vec3 v_normalized = normalize(v);
  vec3 h = normalize(v_normalized + l_normalized);
  vec3 v_normal_normalized =  n_d;
  float specular_coeff = 0.6;
  float ambient_coeff = 0.05;
  float diffuse_coeff = 0.3;
  float p = 64;
  float light_dist_2 = dot(l, l);
  vec4 diffuse_shading = vec4(diffuse_coeff * u_light_intensity / light_dist_2 * clamp(dot(v_normal_normalized, l_normalized), 0.0, 1.0), 1) ;
  vec4 ambient_shading = vec4(ambient_coeff * u_light_intensity, 1);
  vec4 specular_shading = vec4(specular_coeff * u_light_intensity / light_dist_2 * pow(clamp(dot(v_normal_normalized, h), 0.0, 1.0), p), 1);
  // (Placeholder code. You will want to replace it.)
  out_color = diffuse_shading + ambient_shading + specular_shading;
  out_color.a = 1;
}

