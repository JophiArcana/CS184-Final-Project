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
  
  // (Placeholder code. You will want to replace it.)
  // out_color = (vec4(1, 1, 1, 0) + v_normal) / 2;
  // out_color.a = 1;

  vec3 vv = u_light_pos - vec3(v_position);
  vec3 ll = u_cam_pos - vec3(v_position);
  vec3 h = (vv / length(vv)) + (ll / length(ll));
  h = h / length(h);

  float multiplier = max(0.0, dot(vv, vec3(v_normal)));
  float multiplier2 = pow(max(0.0, dot(vec3(v_normal), h)), 100);

  out_color = vec4(0.1, 0.1, 0.1, 0.1) + 0.8 * u_color * multiplier / length(vv) / length(vv) + 4 * u_color * multiplier2 / length(vv) / length(vv);
  // out_color = 4 * u_color * multiplier2 / length(vv) / length(vv);
  out_color[3] = 1.0;
}

