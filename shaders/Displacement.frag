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
  return texture(u_texture_2, uv)[0];
}

void main() {
  // YOUR CODE HERE
  
      vec3 b = cross(vec3(v_normal), vec3(v_tangent));
      mat3 tbn = mat3(vec3(v_tangent), b, vec3(v_tangent));

      float dU = (h(vec2(v_uv[0] + 1.0 / float(u_texture_2_size[0]), v_uv[1])) - h(v_uv)) * u_height_scaling * u_normal_scaling;
      float dV = (h(vec2(v_uv[0], v_uv[1] + 1.0 / float(u_texture_2_size[1]))) - h(v_uv)) * u_height_scaling * u_normal_scaling;

      vec3 no = vec3(-dU, -dV, 1);
      vec3 nd = tbn * no;

      nd = nd / length(nd);

      vec3 vv = u_light_pos - vec3(v_position);
      vec3 ll = u_cam_pos - vec3(v_position);
      vec3 h = (vv / length(vv)) + (ll / length(ll));
      h = h / length(h);

      float multiplier = max(0.0, dot(vv, vec3(nd)));
      float multiplier2 = pow(max(0.0, dot(vec3(nd), h)), 100);
      out_color = vec4(0.1, 0.1, 0.1, 0.1) + 0.8 * u_color * multiplier / length(vv) / length(vv) + 4 * u_color * multiplier2 / length(vv) / length(vv);
      out_color[3] = 1.0;
}

