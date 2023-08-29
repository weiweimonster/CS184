#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  vec3 l = u_light_pos - v_position.xyz;
  vec3 l_normalized = normalize(l);
  vec3 v = u_cam_pos - v_position.xyz;
  vec3 v_normalized = normalize(v);
  vec3 h = normalize(v_normalized + l_normalized);
  vec3 v_normal_normalized = normalize(v_normal.xyz);
  float specular_coeff = 0.6; //0.6
  float ambient_coeff = 0.05; //0.05
  float diffuse_coeff = 0.3;  //0.3
  float p = 8;
  float light_dist_2 = dot(l, l);
  vec4 diffuse_shading = vec4(diffuse_coeff * u_light_intensity / light_dist_2 * clamp(dot(v_normal_normalized, l_normalized), 0.0, 1.0), 1) ;
  vec4 ambient_shading = vec4(ambient_coeff * u_light_intensity, 1);
  vec4 specular_shading = vec4(specular_coeff * u_light_intensity / light_dist_2 * pow(clamp(dot(v_normal_normalized, h), 0.0, 1.0), p), 1);
  
  // (Placeholder code. You will want to replace it.)
  out_color = diffuse_shading + ambient_shading + specular_shading;
  out_color.a = 1;
}

