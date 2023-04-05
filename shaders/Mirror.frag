#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  // YOUR CODE HERE

  vec3 wo = u_cam_pos - vec3(v_position);
  wo /= length(wo); // normalize
  vec3 wi = 2 * vec3(v_normal) - wo;
  // wi should already be normalized

  out_color = texture(u_texture_cubemap, wi);
  out_color.a = 1;
}
